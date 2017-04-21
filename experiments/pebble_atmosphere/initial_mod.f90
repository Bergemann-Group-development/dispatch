!> =============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!> =============================================================================
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
  implicit none
  private
  real(8), parameter:: pi2=8d0*atan(1d0)
  type, public:: initial_t
    character(len=20) :: initial_condition
    real:: b0(3)
    real:: gamma
    real(8):: size(3)
    logical:: mhd = .false.
    character(len=64):: solver
    procedure(planet_atmosphere), pointer:: selected => null()
  contains
    procedure init
  end type
  public planet_atmosphere
  type(lagrange_t):: interp
CONTAINS

!> =============================================================================
!> Setup an initial state, with head wind speed and accretion speed
!> =============================================================================
SUBROUTINE init (self, solver, type)
   class(initial_t):: self
   character(len=64):: solver
   character(len=*), optional:: type
   logical, save:: first_time=.true.
   character(len=20), save:: initial_condition = 'planet_atmosphere'
   namelist /initial_params/ initial_condition
   integer :: iostat
   !----------------------------------------------------------------------------
   call trace%begin('initial_t%init')
   if (present(type)) then
     initial_condition = type
   end if
   call scaling%init
   if (first_time) then
     first_time = .false.
     rewind (io%input); read (io%input, initial_params, iostat=iostat)
     if (io%master) write (*, initial_params)
   end if
   self%initial_condition = initial_condition
   self%solver = trim(solver)
   select case(trim(self%initial_condition))
     case('star_and_planet')
       self%selected => planet_atmosphere
     case('planet_atmosphere')
       self%selected => planet_atmosphere
     case('disk')
       self%selected => disk
     case default
       error stop "initial condition not understood"
     end select
   call trace%end
END SUBROUTINE init

!> =============================================================================
!> Interpolate in a hydrostatic planet atmosphere, computed in atmosphere_mod
!> =============================================================================
SUBROUTINE planet_atmosphere (self, m, mem)
  class(initial_t) :: self
  real, dimension(:,:,:,:), pointer  :: mem
  class(mesh_t), pointer, dimension(:) :: m
  !.............................................................................
  integer :: ix, iy, iz, nr
  real(8) :: x, y, z, rp, d, d0, T0, p, ux, uy, uz, ds2(3)
  real(8) :: xs, ys, zs, xp, yp, zp
  type(index_t):: i
  !-----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  real, save:: omega(3)=0.0
  real, save:: b0(3) = 0.0, gamma=1.4, v_head_cgs=0.0, v_accrete_cgs=0.0
  namelist /planet_atm_params/ omega, gamma, b0, v_head_cgs, v_accrete_cgs
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('initial_t%planet_atmosphere', itimer=itimer)
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    omega = planet%omega*cgs%yr/(2.*cgs%pi*scaling%t)
    rewind (io%input)
    read (io%input, planet_atm_params)
    if (io%master) write (*, planet_atm_params)
  end if
  !$omp end critical (input_cr)
  planet%omega   = omega*2.*cgs%pi*scaling%t/cgs%yr
  planet%v_head    =    v_head_cgs/scaling%u
  planet%v_accrete = v_accrete_cgs/scaling%u
  T0 = atmosphere%external%temperature
  nr = atmosphere%nr
  select case(trim(self%solver))
  case('ramses_hd_patch','ramses_mhd_patch')
    ds2 = 0
  case('stagger2e_hd_patch','stagger2e_mhd_patch')
    ds2 = 0.5*m%d
  end select
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
        rp = sqrt(xp**2+yp**2+zp**2)
        d = interpolate (rp, nr ,atmosphere%rr, atmosphere%rho)
        p = interpolate (rp, nr, atmosphere%rr, atmosphere%pg)
        !-----------------------------------------------------------------------
        ! Vertical stratification is added by adding Gaussian profile and
        ! subtracting off the otherwise constant external density
        !-----------------------------------------------------------------------
        d0 = d
        d = d*exp(-(planet%omega(3)**2*0.5*(z**2+2.*planet%trap*xp**2)/T0))
        p = p*(d/d0)
        mem(ix,iy,iz,i%d) = d
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
        uy = -(1.5+planet%trap)*planet%omega(3)*xs - planet%v_head
        uz = 0.0    
        select case (trim(self%solver))
          case('ramses_hd_patch','ramses_mhd_patch')
            mem(ix,iy,iz,i%e) = p/(gamma-1.0) + 0.5*d*(ux**2+uy**2+uz**2)
          case('stagger2e_hd_patch','stagger2e_mhd_patch')
            mem(ix,iy,iz,i%e) = p/(gamma-1.0)
          case('stagger2_hd_patch','stagger2_mhd_patch')
            mem(ix,iy,iz,i%s) = d*(log(p)-gamma*log(d))/(gamma-1.0)
          case default
            print *,trim(self%solver)
            error stop "initial_mod (1): solver type not known"
        end select
        if (self%mhd) then
          mem(ix,iy,iz,i%bx) = b0(1)
          mem(ix,iy,iz,i%by) = b0(2)
          mem(ix,iy,iz,i%bz) = b0(3)
        end if
        mem(ix,iy,iz,i%px) = d*ux
        mem(ix,iy,iz,i%py) = d*uy
        mem(ix,iy,iz,i%pz) = d*uz
      end do
    end do
  end do
  
  end associate
  call trace%end (itimer)
END SUBROUTINE planet_atmosphere

!> =============================================================================
!> Interpolate in a hydrostatic planet atmosphere, computed in atmosphere_mod
!> =============================================================================
SUBROUTINE disk (self, m, mem)
  class(initial_t) :: self
  real, dimension(:,:,:,:), pointer  :: mem
  class(mesh_t), pointer, dimension(:) :: m
  !.............................................................................
  integer :: ix, iy, iz, nr
  real(8) :: x, y, z, r, d, d0, T0, p, ux, uy, uz
  real(8) :: xs, ys, zs, xp, yp, zp, ds2(3)
  type(index_t):: i
  !-----------------------------------------------------------------------------
  logical, save:: first_time=.true.
  real, save:: omega(3)=0.0
  real, save:: b0(3) = 0.0, gamma=1.4, v_head_cgs=0.0, v_accrete_cgs=0.0
  namelist /planet_atm_params/ omega, gamma, b0, v_head_cgs, v_accrete_cgs
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('initial_t%planet_atmosphere', itimer=itimer)
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    omega = planet%omega*cgs%yr/(2.*cgs%pi*scaling%t)
    rewind (io%input)
    read (io%input, planet_atm_params)
    if (io%master) write (*, planet_atm_params)
  end if
  !$omp end critical (input_cr)
  planet%omega   = omega*2.*cgs%pi*scaling%t/cgs%yr
  planet%v_head    =    v_head_cgs/scaling%u
  planet%v_accrete = v_accrete_cgs/scaling%u
  nr = atmosphere%nr
  T0 = atmosphere%external%temperature
  select case(trim(self%solver))
  case('ramses_hd_patch','ramses_mhd_patch')
    ds2 = 0
  case('stagger2e_hd_patch','stagger2e_mhd_patch')
    ds2 = 0.5*m%d
  end select
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
        mem(ix,iy,iz,i%d) = d
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
        uy = -(1.5+planet%trap)*planet%omega(3)*xs - planet%v_head
        uz = 0.0    
        select case (trim(self%solver))
          case('ramses_hd_patch','ramses_mhd_patch')
            mem(ix,iy,iz,i%e) = p/(gamma-1.0) + 0.5*d*(ux**2+uy**2+uz**2)
          case('stagger2e_hd_patch','stagger2e_mhd_patch')
            mem(ix,iy,iz,i%e) = p/(gamma-1.0)
          case('stagger2_hd_patch','stagger2_mhd_patch')
            mem(ix,iy,iz,i%s) = d*(log(p)-gamma*log(d))/(gamma-1.0)
          case default
            print *,trim(self%solver)
            error stop "initial_mod (1): solver type not known"
        end select
        if (self%mhd) then
          mem(ix,iy,iz,i%bx) = b0(1)
          mem(ix,iy,iz,i%by) = b0(2)
          mem(ix,iy,iz,i%bz) = b0(3)
        end if
        mem(ix,iy,iz,i%px) = d*ux
        mem(ix,iy,iz,i%py) = d*uy
        mem(ix,iy,iz,i%pz) = d*uz
      end do
    end do
  end do
  end associate
  call trace%end (itimer)
END SUBROUTINE disk

!> =============================================================================
!> Use default Lagrange interpolation
!> =============================================================================
FUNCTION interpolate (x0, n, x, y) RESULT(f)
  implicit none
  integer:: n
  real(kind=8):: x0, x(n), y(n), f
  !.............................................................................
  interp%verbose = max(io%verbose-4,0)
  if (x0 < x(1)) then
    f = y(1)
  else if (x0 > x(n)) then
    f = y(n)
  else
    f = interp%sequence (x0, x, y)
  end if
END FUNCTION

END MODULE
