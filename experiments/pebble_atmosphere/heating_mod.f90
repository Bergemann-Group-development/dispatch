MODULE heating_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE kinds_mod
  USE mesh_mod
  USE atmosphere_mod
  USE planet_mod
  USE eos_mod
  USE rt_mod
  USE scaling_mod
  USE units_mod
  USE bits_mod
  USE lagrange_mod
  implicit none
  private
  !.............................................................................
  type, public, extends(gpatch_t):: heating_t
    real:: cdtd_l = 0.5
    real:: mdot   = 0.0
    type(rt_t), pointer:: rt => null()
  contains
  end type
  real, save:: epsilon=1.0
  real, save:: diffusion_radius=10.0
  real, save:: diffuse_radius=10.0
  real, save:: disk_heat_cgs=0.0
  integer, save:: verbose=0
  logical, save:: use_attenuation=.true.
  ! real, save:: temp0=0.0
  ! real, save:: epsilon=1.0
  ! real, save:: diffusion_radius=10.0
  ! real, save:: diffuse_radius=10.0
  ! real, save:: disk_heat_cgs=0.0
  ! logical, save:: spherical_bc_corners=.false.
  ! logical, save:: use_point_source=.false.
  ! logical, save:: use_attenuation=.true.
  ! logical :: sph_bnd(8)
  ! real(8):: me_per_year=0.0
  ! real(8):: r2, x1, x2, y1, y2, z1, z2
  ! self%mdot = me_per_year*cgs%m_earth/cgs%yr*scaling%t/scaling%m
CONTAINS

SUBROUTINE init (self)
  !rt_boundaries%diffuse_radius = diffuse_radius
END SUBROUTINE init

!===============================================================================
!> Compute radius squared
!===============================================================================
SUBROUTINE radius2 (self, r2)
  class(heating_t)           :: self
  real(var_k), dimension(:,:,:) :: r2
  class(mesh_t), pointer        :: m(:)
  real(8)                       :: x, y, z
  integer                       :: ix, iy, iz
  !-----------------------------------------------------------------------------
  call trace%begin('heating_t%radius2')
  m => self%mesh
  associate (rx=>m(1)%r, ry=>m(2)%r, rz=>m(3)%r)
  do iz=m(3)%lb,m(3)%ub
    z = m(3)%p - planet%position(3) + rz(iz)
    do iy=m(2)%lb,m(2)%ub
      y = m(2)%p - planet%position(2) + ry(iy)
      do ix=m(1)%lb,m(1)%ub
        x = m(1)%p - planet%position(1) + rx(ix)
        r2(ix,iy,iz) = max(x**2+y**2+z**2,m(1)%d)
      end do
    end do
  end do
  end associate
  call trace%end()
END SUBROUTINE radius2

!===============================================================================
!> If RT is on, add the radiative cooling / heating to the entropy time derivative
!===============================================================================
SUBROUTINE q_diffuse (self, q, d, pg)
  class (heating_t) :: self
  real(var_k), dimension(:,:,:), pointer :: q, d, pg
  real(var_k), dimension(:,:,:), allocatable :: fmax
  real(var_k), dimension(:,:,:), pointer :: dedt
  real:: u_max
  integer:: ix
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (eos%do_rt.and.associated(self%rt)) then
    call trace%begin ('heating_t%q_diffuse', itimer=itimer)
    self%rt%epsilon = epsilon
    call self%rt%update
    associate(l => self%li, u => self%ui)
    dedt => q(l(1):u(1),l(2):u(2),l(3):u(3))
    call self%rt%get_q (self%time, dedt)
    !if (self%id==286) print *, 'dedt min,max:', minval(dedt), maxval(dedt)
    !if (self%id==286) print *, 'dedt min,max:', minval(q   ), maxval(q   )
    end associate
    !---------------------------------------------------------------------------
    ! Time step limit
    !---------------------------------------------------------------------------
    allocate(fmax(self%gn(1), self%gn(2), self%gn(3)))
    fmax = (pg/d)**4/pg*(eos%kappa*d*minval(self%ds))
    fmax = fmax/(1.+(eos%kappa*d*minval(self%ds))**2)
    u_max = scaling%stefan*16./3.*self%fmaxval(fmax) &
          * self%courant/self%cdtd_l*epsilon
    deallocate(fmax)
    if (self%debug(2)) print 1, self%id, 'u_max(RT)', u_max, self%u_max
    1 format(i6,3x,a,1p,2g12.3)
    self%u_max = max(self%u_max,u_max)
    call trace%end (itimer)
  end if
END SUBROUTINE q_diffuse

!===============================================================================
!> Use diffusion approximation Q = d2S/dtau2 to compute the heating exchange
!===============================================================================
SUBROUTINE q_diffusion (self, rr2, q)
  class (heating_t) :: self
  real(var_k), dimension(:,:,:) :: rr2
  real(var_k), dimension(:,:,:), pointer :: q
  !.............................................................................
  integer:: ix, iy, iz, l(3), u(3), jx, jy, jz
  real:: dtau1, dtau2, dtau, rad2, c, c2, r(3), f, qd
  integer, save:: nprint=50
  !-----------------------------------------------------------------------------
  call eos%lookup (e=self%rt%e, d=self%rt%d, rk=self%rt%rk, src=self%rt%src)
  jx = self%mesh(1)%li - self%rt%mesh(1)%li
  jy = self%mesh(2)%li - self%rt%mesh(2)%li
  jz = self%mesh(3)%li - self%rt%mesh(3)%li
  if ( (sqrt(sum(self%position**2)))*0.5 < diffusion_radius) then
    if (nprint>0) then
      print *, self%id, 'DIFFUSION', self%position
      nprint = nprint-1
    end if
    l = self%rt%mesh%li
    u = self%rt%mesh%ui
    c = 1./3.
    c2 = log(1e4)
    rad2 = diffusion_radius**2
    associate (r1=>self%rt%mesh(1)%r, r2=>self%rt%mesh(2)%r, r3=>self%rt%mesh(3)%r)
    do iz=l(3),u(3)
      r(3) = self%rt%mesh(3)%p - planet%position(3) + r3(iz)
      do iy=l(2),u(2)
        r(2) = self%rt%mesh(2)%p - planet%position(2) + r2(iy)
        do ix=l(1),u(1)
          r(1) = self%rt%mesh(1)%p - planet%position(1) + r1(ix)
          f = exp(-(sum(r**2)/rad2)**3)
          ! --- x ---
          dtau1 = 0.5*(self%rt%rk(ix-1,iy,iz,1)+self%rt%rk(ix,iy,iz,1))*self%mesh(1)%d
          dtau2 = 0.5*(self%rt%rk(ix+1,iy,iz,1)+self%rt%rk(ix,iy,iz,1))*self%mesh(1)%d
          dtau  = 0.5*(dtau1+dtau2)
          qd = &
             + ((self%rt%src(ix+1,iy,iz,1)-self%rt%src(ix  ,iy,iz,1))/dtau2 &
              - (self%rt%src(ix  ,iy,iz,1)-self%rt%src(ix-1,iy,iz,1))/dtau1)/dtau*c
          ! --- x ---
          dtau1 = 0.5*(self%rt%rk(ix,iy-1,iz,1)+self%rt%rk(ix,iy,iz,1))*self%mesh(2)%d
          dtau2 = 0.5*(self%rt%rk(ix,iy+1,iz,1)+self%rt%rk(ix,iy,iz,1))*self%mesh(2)%d
          dtau  = 0.5*(dtau1+dtau2)
          qd = qd &
             + ((self%rt%src(ix,iy+1,iz,1)-self%rt%src(ix,iy  ,iz,1))/dtau2 &
              - (self%rt%src(ix,iy  ,iz,1)-self%rt%src(ix,iy-1,iz,1))/dtau1)/dtau*c
          ! --- x ---
          dtau1 = 0.5*(self%rt%rk(ix,iy,iz-1,1)+self%rt%rk(ix,iy,iz,1))*self%mesh(3)%d
          dtau2 = 0.5*(self%rt%rk(ix,iy,iz+1,1)+self%rt%rk(ix,iy,iz,1))*self%mesh(3)%d
          dtau  = 0.5*(dtau1+dtau2)
          qd = qd &
             + ((self%rt%src(ix,iy,iz+1,1)-self%rt%src(ix,iy,iz  ,1))/dtau2 &
              - (self%rt%src(ix,iy,iz  ,1)-self%rt%src(ix,iy,iz-1,1))/dtau1)/dtau*c
          ! --- combine ---
          q(ix+jx,iy+jy,iz+jz) = (1.-f)*q(ix+jx,iy+jy,iz+jz) + f*qd
        end do
      end do
    end do
    end associate
  end if
END SUBROUTINE q_diffusion

!===============================================================================
!> Point source heating from the central embryo, approximated as the flux from
!> a block body with the given surface temperature, attenuated by the surface
!> area factor and the optical depth attenuation. That intensity is then absorbed
!> in a fraction \Delta\tau (which is manifestly smaller than one in practice).
!===============================================================================
SUBROUTINE q_planet (self, r2, q, d, pg)
  class (heating_t) :: self
  real(var_k), dimension(:,:,:) :: r2, q, d, pg
  !.............................................................................
  real(var_k), dimension(:,:,:), pointer :: s
  real :: i_norm, tau(2), g1, rk, rt_int, temp, temp_ext, rt_src, u_max, w_omega
  real(8) :: rp, rad2, w_z, stefan_pi
  integer :: ix, iy, iz
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('heating_t%q_planet', itimer=itimer)
  !---------------------------------------------------------------------------
  ! The intensity at radius r is the surface intensity, \sigma*T**4/\pi,
  ! times planet%radius**2/r**2, times exp(-tau), where tau is the difference
  ! in optical depth from the surface to r.
  !---------------------------------------------------------------------------
  stefan_pi = scaling%stefan / cgs%pi
  i_norm = stefan_pi*atmosphere%surface%temperature**4
  s => self%mem(:,:,:,self%idx%s,self%it,1)
  g1 = eos%gamma-1.0
  rad2 = (5.*planet%radius)**2
  u_max = 0.0
  do iz=self%mesh(3)%li,self%mesh(3)%ui
  do iy=self%mesh(2)%li,self%mesh(2)%ui
    do ix=self%mesh(1)%li,self%mesh(1)%ui
      rp = sqrt(r2(ix,iy,iz))
      !-------------------------------------------------------------------------
      ! Reduced intensity, from optical depth distance from the surface and
      ! radial attenuation
      !-------------------------------------------------------------------------
      rt_int = i_norm
      if (use_attenuation) then
        tau(1) = lagrange%interpolate1d (rp, atmosphere%nr ,atmosphere%rr, atmosphere%tau)
        tau(2) = max(0.0,atmosphere%surface%tau-tau(1))
        rt_int = rt_int*exp(-tau(2))
      end if
      !-------------------------------------------------------------------------
      ! Let the source term cool only towards the external temperature, and
      ! taper that effect off toward the embryo, to avoid edge heating
      !-------------------------------------------------------------------------
      temp = d(ix,iy,iz)**g1*exp(s(ix,iy,iz)/d(ix,iy,iz)*g1)
      !temp_ext = atmosphere%external%temperature*min(1.0,max(0.0,0.1*(rp-planet%radius)/planet%radius))
      temp_ext = 0.0
      rt_src = stefan_pi*(temp-temp_ext)**4
      !-------------------------------------------------------------------------
      ! Weight of point source relative to diffuse light, normalized to 1
      !-------------------------------------------------------------------------
      if (eos%do_rt) then
        w_omega = 0.25*planet%radius**2/(r2(ix,iy,iz)-0.5*planet%radius**2)
      else
        w_omega = 1.0
      end if
      !-------------------------------------------------------------------------
      ! Add the contributions.  The incoming q(ix,iy,iz) is from q_diffuse,
      ! and is already normalized with 4.*pi*rk.  It already contains a term
      ! -4.*pi*rk*source, but since we need to normalize the contribution 
      ! with (1.-w_omega), we need to add -w_omega*4.*pi*rk*src
      !-------------------------------------------------------------------------
      rk = eos%kappa*d(ix,iy,iz)*4.*cgs%pi*epsilon
      w_z = exp(-rad2/max(rad2/log(1e9),r2(ix,iy,iz)))
      !if (iy==17.and.iz==17.and.self%id==285) print *, ix, i_norm, rt_int, rt_src, w_omega, w_z, &
      !  q(ix,iy,iz)*(1.-w_omega)*w_z + rk*(rt_int*w_omega-rt_src*(1.-w_omega)*(1.-w_z))
      q(ix,iy,iz) = q(ix,iy,iz)*(1.-w_omega)*w_z + rk*(rt_int*w_omega-rt_src*(1.-w_omega))*exp(-tau(1))*(1.-w_z)
      !-------------------------------------------------------------------------
      ! Estimate the max rate of change, from assuming no cooling
      !-------------------------------------------------------------------------
      if (rp > planet%radius-self%mesh(1)%d) then
        u_max = max(u_max,rk*rt_src/pg(ix,iy,iz)) !!!*w_omega
      end if
    end do
  end do
  end do
  u_max = u_max*self%mesh(1)%d*self%courant/self%cdtd_l
  if (self%boundaries%is_set(bits%spherical) .and. verbose>1) then
    print *, self%id, 'u_max(PS)', u_max, self%u_max
  end if
  self%u_max = max(self%u_max, u_max)
  call trace%end (itimer)
contains
  function trad (q)
    real:: q, trad
    trad = (max(q,0.0)/stefan_pi)**0.25
  end function
END SUBROUTINE q_planet

!===============================================================================
!> Add accretion heating
!===============================================================================
SUBROUTINE q_accretion (self, r2, q)
  class(heating_t)                    :: self
  real(var_k), dimension(:,:,:)          :: r2, q
  real(var_k), dimension(:,:,:), pointer :: d
  real(var_k)                            :: heat(2), pl_r2, d0
  real(8), dimension(:), allocatable     :: lnrr
  real(8)                                :: lnr, cx
  integer                                :: ix, iy, iz, nr
  integer, save                          :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%mdot>0.0 .or. disk_heat_cgs>0.0) then
    call trace%begin('heating_t%q_accretion', itimer=itimer)
    !---------------------------------------------------------------------------
    ! q is the heating rate from accretion, Q = Mdot G M / (4\pi r^4), in energy
    ! units per volume unit and time unit.
    !---------------------------------------------------------------------------
    pl_r2 = planet%radius**2
    heat(1) = self%mdot*scaling%grav*planet%mass/(4.*cgs%pi)
    heat(2) = disk_heat_cgs/(scaling%e/scaling%m/scaling%t)
    d => self%mem(:,:,:,self%idx%d,self%it,1)
    nr = atmosphere%nr
    allocate (lnrr(nr))
    lnrr = log(atmosphere%rr)
    cx = (nr-1)/(lnrr(nr)-lnrr(1))
    do iz=self%mesh(3)%li,self%mesh(3)%uo
    do iy=self%mesh(2)%li,self%mesh(2)%uo
    do ix=self%mesh(1)%li,self%mesh(1)%uo
      if (r2(ix,iy,iz) > pl_r2) then
        lnr = log(r2(ix,iy,iz))*0.5
        d0 = lagrange%interpolate1d (lnr, nr, lnrr, atmosphere%rho)
        !if (iy==17.and.iz==17.and.self%id==285) &
        !  print *,'xxxx', ix, d(ix,iy,iz)/d0, q(ix,iy,iz), heat(1)/d0/r2(ix,iy,iz)**2*d(ix,iy,iz), heat(2)*d(ix,iy,iz)
        q(ix,iy,iz) = q(ix,iy,iz) + (heat(1)/d0/r2(ix,iy,iz)**2 + heat(2))*d(ix,iy,iz)
      end if
    end do
    end do
    end do
    deallocate (lnrr)
    call trace%end (itimer)
  end if
contains
  !-----------------------------------------------------------------------------
  ! Interpolate in equidistant table
  !-----------------------------------------------------------------------------
  function interpolate (x0, nx, x, y) result (out)
    integer:: nx, ix
    real(8):: x0, x(nx), y(nx), px, out
    px = (x0-x(1))*cx + 1d0
    ix = px
    px = ix-1
    out = y(ix)*(1d0-px) + y(ix+1)*px
  end function interpolate
END SUBROUTINE q_accretion

END MODULE heating_mod
