!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE atmosphere_mod
  USE planet_mod
  USE lagrange_mod
  USE units_mod
  USE scaling_mod
  USE index_mod
  USE kinds_mod
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    integer:: id
    character(len=20) :: initial_condition
    real:: b0(3)
    real:: gamma
    real(8):: size(3)
    logical:: mhd = .false.
    character(len=64):: solver
    procedure(planet_atmosphere), pointer:: condition => null()
  contains
    procedure init
  end type
  logical, save:: test=.false.
  public planet_atmosphere
CONTAINS

!===============================================================================
!> Setup an initial state, with head wind speed and accretion speed.  Note
!> that the value of initial_condition is (also) set by the 'type' parameter
!> in the experiment_mod.  This is kept only in order not to break existing
!> setups -- the proper place to set the type of IC is here, via the namelist
!> below
!===============================================================================
SUBROUTINE init (self, solver, gamma)
   class(initial_t):: self
   character(len=64):: solver
   real, optional:: gamma
   !............................................................................
   logical, save:: first_time=.true.
   character(len=20), save:: initial_condition = 'planet_atmosphere'
   namelist /IC_params/ initial_condition
   integer :: iostat
   !----------------------------------------------------------------------------
   call trace%begin('initial_t%init')
   call scaling%init
   call io%lock%set
   if (first_time) then
     first_time = .false.
     initial_condition = self%initial_condition
     rewind (io%input); read (io%input, IC_params, iostat=iostat)
     if (io%master) write (*, IC_params)
   end if
   call io%lock%unset
   self%initial_condition = initial_condition
   self%solver = trim(solver)
   select case(trim(self%initial_condition))
     case('star_and_planet')
       self%condition => planet_atmosphere
     case('planet_atmosphere')
       self%condition => planet_atmosphere
     case('disk')
       self%condition => disk
     case default
       error stop "initial condition not understood"
     end select
   call trace%end()
END SUBROUTINE init

!===============================================================================
!> Interpolate in a hydrostatic planet atmosphere, computed in atmosphere_mod
!===============================================================================
SUBROUTINE planet_atmosphere (self, m, mem, idx)
  class(initial_t) :: self
  type(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: mem
  class(mesh_t), pointer, dimension(:) :: m
  !.............................................................................
  integer :: ix, iy, iz, nr, l(3), u(3)
  real(8) :: x, y, z, rp, d, ds, d0, T0, p, ux, uy, uz, ds2(3)
  real(8) :: xs, ys, zs, xp, yp, zp, R_Hill, f_Hill
  !-----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  real, save:: omega(3)=0.0
  real, save:: b0(3) = 0.0, gamma=1.4, v_head_cgs=0.0, v_accrete_cgs=0.0
  namelist /planet_atm_params/ omega, gamma, b0, v_head_cgs, v_accrete_cgs, test
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('initial_t%planet_atmosphere', itimer=itimer)
  call io%lock%set
  if (first_time) then
    first_time = .false.
    omega = planet%omega*cgs%yr/(2.*cgs%pi*scaling%t)
    rewind (io%input)
    read (io%input, planet_atm_params)
    if (io%master) write (*, planet_atm_params)
  end if
  call io%lock%unset
  planet%omega   = omega*2.*cgs%pi*scaling%t/cgs%yr
  planet%v_head    =    v_head_cgs/scaling%u
  planet%v_accrete = v_accrete_cgs/scaling%u
  T0 = atmosphere%external%temperature
  nr = atmosphere%nr
  if (self%solver(1:7)=='stagger') then
    ds2 = 0.5*m%d
  else if (self%solver(1:6)=='zeus3d') then
    ds2 = 0.5*m%d
  else
    ds2 = 0.0
  endif
  if (test) then
   l = m%li
   u = m%ui
  else
   l = m%lb
   u = m%ub
  end if
  R_Hill = cgs%au/scaling%l*(cgs%m_earth/(3.*cgs%m_sun))**(1./3.)
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
  do iz=l(3),u(3)
    z = m(3)%p + r3(iz)
    zp = z  - planet%position(3)
    zs = zp - ds2(3)
    do iy=l(2),u(2)
      y = m(2)%p + r2(iy)
      yp = y  - planet%position(2)
      ys = yp - ds2(2)
      do ix=l(1),u(1)
        x = m(1)%p + r1(ix)
        xp = x  - planet%position(1)
        xs = xp - ds2(1)
        rp = sqrt(xp**2+yp**2+zp**2)
        d = lagrange%interpolate1d (rp, nr ,atmosphere%rr, atmosphere%rho)
        p = lagrange%interpolate1d (rp, nr, atmosphere%rr, atmosphere%pg)
        !-----------------------------------------------------------------------
        ! Vertical stratification is added by adding Gaussian profile and
        ! subtracting off the otherwise constant external density
        !-----------------------------------------------------------------------
        d0 = d
        d = d*exp(-(planet%omega(3)**2*0.5*(zp**2+2.*planet%trap*xp**2)/T0))
        p = p*(d/d0)
        mem(ix,iy,iz,idx%d) = d
        !-----------------------------------------------------------------------
        ! Non-zero initial velocity components, to be used in tests only
        !-----------------------------------------------------------------------
        if (test) then
          ux = -planet%v_accrete + planet%omega(3)*(yp-xs+zp)*0.01
          uy = -(1.5+planet%trap)*planet%omega(3)*(xp-(ys+zp))*0.01 - planet%v_head
          uz = planet%omega(3)*(zs+xp-yp)*0.01
        !-----------------------------------------------------------------------
        ! A rotation profile, where the rotation center is at x=-1 AU, and
        ! the direction of rotation is +y, so the head wind is in the -y
        ! direction.  The slower rotation at larger radii means that Uy is
        ! proportional to -x.  This gives a speed difference of only 130 cm/s
        ! across the Earth, and plus/minus 0.15 km/s across two Hill radii.
        ! The effect is not enough to significantly flatten the outer part of
        ! the Hill sphere, so we will probably only be seeing effects of the
        ! Coriolis force beyond the Hill radius.
        !-----------------------------------------------------------------------
        else
          ux = -planet%v_accrete
          uy = -(1.5+planet%trap)*planet%omega(3)*xp - planet%v_head
          uz = 0.0
        end if
        f_Hill = (1.-exp(-rp**2/(2.*R_Hill**2)))
        uy = uy*f_Hill
        !-----------------------------------------------------------------------
        ! Momenta from (possibly) staggered densities
        !-----------------------------------------------------------------------
        rp = sqrt(xs**2+yp**2+zp**2)
        ds = lagrange%interpolate1d (rp, nr ,atmosphere%rr, atmosphere%rho)
        ds = ds*exp(-(planet%omega(3)**2*0.5*(zp**2+2.*planet%trap*xs**2)/T0))
        mem(ix,iy,iz,idx%px) = ds*ux
        rp = sqrt(xp**2+ys**2+zp**2)
        ds = lagrange%interpolate1d (rp, nr ,atmosphere%rr, atmosphere%rho)
        ds = ds*exp(-(planet%omega(3)**2*0.5*(zp**2+2.*planet%trap*xp**2)/T0))
        mem(ix,iy,iz,idx%py) = ds*uy
        rp = sqrt(xp**2+yp**2+zs**2)
        ds = lagrange%interpolate1d (rp, nr ,atmosphere%rr, atmosphere%rho)
        ds = ds*exp(-(planet%omega(3)**2*0.5*(zs**2+2.*planet%trap*xp**2)/T0))
        mem(ix,iy,iz,idx%pz) = ds*uz
        if (self%solver(1:7)=='ramses_') then
          mem(ix,iy,iz,idx%e) = p/(gamma-1.0) + 0.5*d*(ux**2+uy**2+uz**2)
        else if (self%solver(1:9)=='stagger2e') then
          mem(ix,iy,iz,idx%e) = p/(gamma-1.0)
        else if (self%solver(1:9)=='stagger2_') then
          mem(ix,iy,iz,idx%s) = d*(log(p)-gamma*log(d))/(gamma-1.0)
        else
          print *,trim(self%solver)
          error stop "initial_mod (1): solver type not known"
        end if
        if (self%mhd) then
          mem(ix,iy,iz,idx%bx) = b0(1)
          mem(ix,iy,iz,idx%by) = b0(2)
          mem(ix,iy,iz,idx%bz) = b0(3)
        end if
      end do
    end do
  end do
  end associate
  call trace%end (itimer)
END SUBROUTINE planet_atmosphere

!===============================================================================
!> Interpolate in a hydrostatic planet atmosphere, computed in atmosphere_mod
!===============================================================================
SUBROUTINE disk (self, m, mem, idx)
  class(initial_t) :: self
  type(index_t):: idx
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: mem
  class(mesh_t), pointer, dimension(:) :: m
  !.............................................................................
  integer :: ix, iy, iz, nr
  real(8) :: x, y, z, r, d, ds, d0, T0, p, ux, uy, uz, ds2(3)
  real(8) :: xs, ys, zs, xp, yp, zp, rp
  !-----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  real, save:: omega(3)=0.0
  real, save:: b0(3) = 0.0, gamma=1.4, v_head_cgs=0.0, v_accrete_cgs=0.0
  namelist /planet_atm_params/ omega, gamma, b0, v_head_cgs, v_accrete_cgs
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('initial_t%disk', itimer=itimer)
  call io%lock%set
  if (first_time) then
    first_time = .false.
    omega = planet%omega*cgs%yr/(2.*cgs%pi*scaling%t)
    rewind (io%input)
    read (io%input, planet_atm_params)
    if (io%master) write (*, planet_atm_params)
  end if
  call io%lock%unset
  planet%omega   = omega*2.*cgs%pi*scaling%t/cgs%yr
  planet%v_head    =    v_head_cgs/scaling%u
  planet%v_accrete = v_accrete_cgs/scaling%u
  nr = atmosphere%nr
  T0 = atmosphere%external%temperature
  if (self%solver(1:7)=='stagger') then
    ds2 = 0.5*m%d
  else
    ds2 = 0.0
  endif
  associate (r1=>m(1)%r, r2=>m(2)%r, r3=>m(3)%r)
    do iz=m(3)%lb,m(3)%ub
    z = m(3)%p + r3(iz)
    zp = z  - planet%position(3)
    zs = zp - ds2(3)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p + r2(iy)
      yp = y  - planet%position(2)
      ys = yp - ds2(2)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p + r1(ix)
        xp = x  - planet%position(1)
        xs = xp - ds2(1)
        d = atmosphere%external%density
        p = atmosphere%external%pressure
        !-----------------------------------------------------------------------
        ! Vertical stratification is added by adding Gaussian profile and
        ! subtracting off the otherwise constant external density
        !-----------------------------------------------------------------------
        d0 = d
        d = d*exp(-(planet%omega(3)**2*0.5*(z**2+2.*planet%trap*xp**2)/T0))
        p = p*(d/d0)
        mem(ix,iy,iz,idx%d) = d
        !-----------------------------------------------------------------------
        ! A rotation profile, where the rotation center is at x=-1 AU, and
        ! the direction of rotation is +y, so the head wind is in the -y
        ! direction.  The slower rotation at larger radii means that Uy is
        ! proportional to -x.  This gives a speed difference of only 130 cm/s
        ! across the Earth, and plus/minus 0.15 km/s across two Hill radii.
        ! The effect is not enough to significantly flatten the outer part of
        ! the Hill sphere, so we will probably only be seeing effects of the
        ! Coriolis force beyond the Hill radius.
        !-----------------------------------------------------------------------
        ux = -planet%v_accrete
        uy = -(1.5+planet%trap)*planet%omega(3)*xp - planet%v_head
        uz = 0.0
        !-----------------------------------------------------------------------
        ! Non-zero initial velocity components, to be used in tests only
        !-----------------------------------------------------------------------
        !ux = -planet%v_accrete + planet%omega(3)*(yp-xs)
        !uy = -(1.5+planet%trap)*planet%omega(3)*(xp-ys) - planet%v_head
        !uz = planet%omega(3)*zs*0.01
        !-----------------------------------------------------------------------
        ! Momenta from (possibly) staggered densities
        !-----------------------------------------------------------------------
        !rp = sqrt(xs**2+yp**2+zp**2)
        !ds = lagrange%interpolate (rp, nr ,atmosphere%rr, atmosphere%rho)
        !ds = ds*exp(-(planet%omega(3)**2*0.5*(zp**2+2.*planet%trap*xs**2)/T0))
        mem(ix,iy,iz,idx%px) = d*ux
        !rp = sqrt(xp**2+ys**2+zp**2)
        !ds = lagrange%interpolate (rp, nr ,atmosphere%rr, atmosphere%rho)
        !ds = ds*exp(-(planet%omega(3)**2*0.5*(zp**2+2.*planet%trap*xp**2)/T0))
        mem(ix,iy,iz,idx%py) = d*uy
        !rp = sqrt(xp**2+yp**2+zs**2)
        !ds = lagrange%interpolate (rp, nr ,atmosphere%rr, atmosphere%rho)
        !ds = ds*exp(-(planet%omega(3)**2*0.5*(zs**2+2.*planet%trap*xp**2)/T0))
        mem(ix,iy,iz,idx%pz) = d*uz
        if (self%solver(1:7)=='ramses_') then
          mem(ix,iy,iz,idx%e) = p/(gamma-1.0) + 0.5*d*(ux**2+uy**2+uz**2)
        else if (self%solver(1:9)=='stagger2e') then
          mem(ix,iy,iz,idx%e) = p/(gamma-1.0)
        else if (self%solver(1:9)=='stagger2_') then
          mem(ix,iy,iz,idx%s) = d*(log(p)-gamma*log(d))/(gamma-1.0)
        else
          print *,trim(self%solver)
          error stop "initial_mod (1): solver type not known"
        end if
        if (self%mhd) then
          mem(ix,iy,iz,idx%bx) = b0(1)
          mem(ix,iy,iz,idx%by) = b0(2)
          mem(ix,iy,iz,idx%bz) = b0(3)
        end if
      end do
    end do
  end do
  end associate
  call trace%end (itimer)
END SUBROUTINE disk

END MODULE
