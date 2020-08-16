!===============================================================================
!> Setup the initial conditions and experiment update
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE kinds_mod
  USE pebble_patch_mod
  USE initial_mod
  USE atmosphere_mod
  USE units_mod
  USE scaling_mod
  USE planet_mod
  USE mesh_mod
  USE boundaries_mod
  USE scalar_mod
  USE eos_mod
  USE lagrange_mod
  USE timestep_mod
  USE rt_boundaries_mod
  implicit none
  private
  type, public, extends(pebble_patch_t):: experiment_t
    real:: mdot           ! accretion rate in code units
    real(8):: x_sponge    ! distance from the physical boundary where
                          ! "wave-killing" boundaries are turned on.
  contains
    procedure :: init
    procedure :: update
    procedure :: rotate
    procedure :: radius2
    procedure :: q_diffuse
    procedure :: q_diffusion
    procedure :: q_planet
    procedure :: q_accretion
    procedure :: output
    procedure :: boundary_conditions
    procedure :: edge_boundary_conditions
    procedure :: density_boundary_condition
    procedure :: is_relevant
  end type
  real, save:: temp0=0.0
  real, save:: epsilon=1.0
  real, save:: diffusion_radius=10.0
  real, save:: diffuse_radius=10.0
  real, save:: disk_heat_cgs=0.0
  integer, save:: track_level=99
  logical, save:: first_time=.true.
  logical, save:: spherical_bc_corners=.false.
  logical, save:: use_point_source=.false.
  logical, save:: use_attenuation=.true.
  character(len=32), save:: type='planet_atmosphere'
CONTAINS

!===============================================================================
!> The initialization must be done in the right order.  The rotation speed
!> omega must be known before the initial%init runs, for example.
!>
!> dispatch
!>   setup::init
!>     scaling::init                            ! code units
!>     planet::init                             ! planet size, orbit
!>     atmosphere::init                         ! construct hydrostatic state
!>     rubiks_mod::init                         ! makes all patches
!>       experiment_mod::init                   ! calls other inits
!>         mhd_mod::init                        ! reads gamma
!>           force_mod::init                    ! initialize forcing
!>         pebbles_mod::init                    ! reads omega
!>         initial_mod::init                    ! selects ICs
!>         initial_mod::condition               ! computes ICs
!>         boundaries_mod::condition            ! apply BCs to ICs
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  ! ....................................................................
  integer:: iv, iostat
  logical, save:: first_time=.true.
  logical :: sph_bnd(8)
  real(8):: x_sponge=-1.0e30
  real(8):: me_per_year=0.0
  real(8):: r2, x1, x2, y1, y2, z1, z2
  namelist /experiment_params/ type, x_sponge, me_per_year, spherical_bc_corners, &
    use_point_source, use_attenuation, epsilon, diffusion_radius, diffuse_radius, &
    disk_heat_cgs, track_level
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace%begin('experiment_t%init')
  !$omp critical (io_cr)
  if (first_time) then
    first_time = .false.
    if (eos%do_rt) then
      use_attenuation = .true.
    else
      use_attenuation = .false.
    end if
    rewind (io%input)
    read (io%input,experiment_params) !, iostat=iostat)
    if (io%master) write (*,experiment_params)
  end if
  !$omp end critical (io_cr)
  self%x_sponge = x_sponge
  !----------------------------------------------------------------------------
  ! Allocate memory and set IC values.  The initial%init (type) specifies the
  ! type of IC to be applied when the call to gpatch_t%init is made.
  !----------------------------------------------------------------------------
  self%nv = 5
  self%mhd = .false.
  self%mdot = me_per_year*cgs%m_earth/cgs%yr*scaling%t/scaling%m
  self%msplit = 3
  call scaling%init
  self%nv = 5
  self%mhd = .false.
  !----------------------------------------------------------------------------
  ! Through the hierarchy of class(task_t)%init calls, the call to %init below
  ! sets the IC values, either from reading in snapshot values, or by calling
  ! the IC condition procedure self%initial%condition() from gpatch_t%init.
  !----------------------------------------------------------------------------
  self%periodic = .false.
  self%periodic(2) = .true.
  !----------------------------------------------------------------------------
  ! To avoid changing a lot of input files (including many not under GIT control),
  ! we set the initial%initial_condition parameter  (which should really be set
  ! via a namelist in initial_mod) from here.  Then, when mhd_t%init calls
  ! initial%init, the type of IC is known, and the procedure pointer is set
  ! correctly.
  !----------------------------------------------------------------------------
  self%initial%initial_condition = type
  self%initial%solver = 'planet_atmosphere'
  !----------------------------------------------------------------------------
  ! Initialize the variable indices, and copy the self%idx instance to a public
  ! instance, for backwards compatibility
  !----------------------------------------------------------------------------
  call self%idx%init (2, .false.)
  call self%idx%copy
  call self%pebble_patch_t%init
  call self%set (bits%static)
  call self%boundaries%init (self)
  if (self%boundaries%is_set (bits%spherical)) then
    print *, 'experiment_t%init: setting spherical BC bits on patch', self%id
    call self%set (bits%has_BC)
    do iv=self%idx%px,self%idx%pz
      call self%boundaries%spherical (self%mem(:,:,:,iv,self%it,1), 0.0, self%mesh, iv)
    end do
  end if
  call self%boundary_conditions (self%it, self%it)
  !----------------------------------------------------------------------------
  ! Choose interpolation types. The smallest initial glitch (positive) is with
  !
  ! self%unsigned(self%idx%d) = .true.
  ! self%pervolume(self%idx%s) = .false.
  !
  ! All other combinations give a stronger negative initial glitch in p1 at
  ! mesh%li. Both types of glitches come from the different precision in the
  ! initial hydrostatic state, and go away with time. The glitch is visible
  ! in u1 at index li when the pressure on one side comes from lower resolution.
  !
  ! Note that the current default in stagger2/mhd_mod.f90 is
  !
  ! self%unsigned(self%idx%d) = .true.
  ! self%pervolume(self%idx%s) = .true.
  !
  ! Tests show that these should either both be true, or both be false.
  !----------------------------------------------------------------------------
  self%unsigned(self%idx%d) = .true.
  self%pervolume(self%idx%s) = .false.
  !----------------------------------------------------------------------------
  ! Export the patch ids for debugging purposes
  !----------------------------------------------------------------------------
  if (self%level >= track_level) self%track = .true.
  self%initial%id = self%id
  rt_boundaries%diffuse_radius = diffuse_radius
  call trace%end()
END SUBROUTINE init


!===============================================================================
!> Choice of update steps for this experiment
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d, e, dsdt
  real(kind=KindScalarVar), dimension(:,:,:), allocatable :: r2
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: q
  real(8), dimension(3):: h
  integer:: it, iv, nb, i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%update', itimer=itimer)
  call self%check_density ('before update')
  call self%output
  !-----------------------------------------------------------------------------
  ! After MHD, if we are doing pebbles we have the gas pressure.  If not, get it
  !-----------------------------------------------------------------------------
  call self%mhd_t%pre_update
  call self%mhd_t%pde
  if (self%do_pebbles) then
    self%pg => self%pressure(:,:,:,self%it)
  else
    call allocate_scalars (self%gn, self%pg)
    self%pg = self%gas_pressure()
  end if
  self%d => self%mem(:,:,:,self%idx%d,self%it,1)
  !-----------------------------------------------------------------------------
  ! Radiative transfer and accretion heating, both need radius**2
  !-----------------------------------------------------------------------------
  if (trim(type) /= 'disk') then
    allocate (r2(self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn))
    allocate (q (self%mesh(1)%gn, self%mesh(2)%gn, self%mesh(3)%gn))
    q = 0.0
    if (eos%do_rt) call self%q_diffuse (q)
    call self%radius2 (r2)
    if (use_point_source) then
      call self%q_planet (r2, q)
    else if (eos%do_rt) then
      call self%q_diffusion (r2, q)
    end if
    if (eos%do_rt) call self%rt%finish ()
    call self%q_accretion (r2, q)
    !---------------------------------------------------------------------------
    ! Update dsdt and apply spherical boundary condition
    !---------------------------------------------------------------------------
    dsdt => self%mem(:,:,:,self%idx%s,self%it,2)
    dsdt = dsdt + q*self%d/self%pg
    deallocate (r2, q)
    call self%boundaries%spherical_jagged_s0 (dsdt, self%mesh, self%idx%s, 0.0)
  end if
  !-----------------------------------------------------------------------------
  ! Do the actual time step update and apply BCs to the result
  !-----------------------------------------------------------------------------
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh, self%lock)
  if (eos%do_rt) then
    self%rt%time  = self%time
    self%rt%dtime = self%dtime
    self%rt%t(self%it) = self%t(self%it)
  end if
  call self%boundary_conditions (self%it, self%new)
  !-----------------------------------------------------------------------------
  ! Pebble update, done after the time step update, to have access to both states
  !-----------------------------------------------------------------------------
  if (self%do_pebbles) then
    call self%pebble_patch_t%update_pebbles
  else
    call deallocate_scalars (self%pg)
  end if
  call self%mhd_t%post_update
  call self%check_nan ('after timestep', always=.true.)
  call self%counter_update
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Make sure the RT task rotates slots correctly, in perfect sync with this task
!> =============================================================================
SUBROUTINE rotate (self)
  class(experiment_t)                                     :: self
  !-----------------------------------------------------------------------------
  ! Check (already here -- to cover for rt%rotate) that we don't do it twice,
  ! and that we actually have the same information in the rt%task_t data type
  !-----------------------------------------------------------------------------
  if (self%rotated) return
  if (eos%do_rt) then
    self%rt%rotated = .false.
    self%rt%time  = self%time
    self%rt%dtime = self%dtime
    call self%rt%rotate()
  end if
  call self%pebble_patch_t%rotate()
END SUBROUTINE rotate

!> =============================================================================
!> Compute radius squared
!> =============================================================================
SUBROUTINE radius2 (self, r2)
  class(experiment_t)                                     :: self
  real(kind=KindScalarVar), dimension(:,:,:)              :: r2
  class(mesh_t), pointer                                  :: m(:)
  real(8)                                                 :: x, y, z
  integer                                                 :: ix, iy, iz
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%radius2')
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

!> =============================================================================
!> If RT is on, add the radiative cooling / heating to the entropy time derivative
!> =============================================================================
SUBROUTINE q_diffuse (self, q)
  class (experiment_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: q
  real(kind=KindScalarVar), dimension(:,:,:), allocatable :: fmax
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: dedt
  real:: u_max
  integer:: ix
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (eos%do_rt.and.associated(self%rt)) then
    call trace%begin ('experiment_t%q_diffuse', itimer=itimer)
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
    fmax = (self%pg/self%d)**4/self%pg*(eos%kappa*self%d*minval(self%ds))
    fmax = fmax/(1.+(eos%kappa*self%d*minval(self%ds))**2)
    u_max = scaling%stefan*16./3.*self%fmaxval(fmax)*self%courant/self%cdtd*epsilon
    deallocate(fmax)
    if (self%debug(2)) print 1, self%id, 'u_max(RT)', u_max, self%u_max
    1 format(i6,3x,a,1p,2g12.3)
    self%u_max = max(self%u_max,u_max)
    call trace%end (itimer)
  end if
END SUBROUTINE q_diffuse

!> =============================================================================
!> Use diffusion approximation Q = d2S/dtau2 to compute the heating exchange
!> =============================================================================
SUBROUTINE q_diffusion (self, rr2, q)
  class (experiment_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:) :: rr2
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: q
  !.............................................................................
  integer:: ix, iy, iz, l(3), u(3), jx, jy, jz
  real:: dtau1, dtau2, dtau, rad2, c, c2, r(3), f, qd
  integer, save:: nprint=50
  !-----------------------------------------------------------------------------
  call eos%lookup (shape(self%rt%d), x=self%rt%e, y=self%rt%d, rk=self%rt%rk, src=self%rt%src)
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

!> =============================================================================
!> Point source heating from the central embryo, approximated as the flux from
!> a block body with the given surface temperature, attenuated by the surface
!> area factor and the optical depth attenuation. That intensity is then absorbed
!> in a fraction \Delta\tau (which is manifestly smaller than one in practice).
!> =============================================================================
SUBROUTINE q_planet (self, r2, q)
  class (experiment_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:) :: r2, q
  !.............................................................................
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: s
  real :: i_norm, tau(2), g1, rk, rt_int, temp, temp_ext, rt_src, u_max, w_omega
  real(8) :: rp, rad2, w_z, stefan_pi
  integer :: ix, iy, iz
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%q_planet', itimer=itimer)
  !---------------------------------------------------------------------------
  ! The intensity at radius r is the surface intensity, \sigma*T**4/\pi,
  ! times planet%radius**2/r**2, times exp(-tau), where tau is the difference
  ! in optical depth from the surface to r.
  !---------------------------------------------------------------------------
  stefan_pi = scaling%stefan / cgs%pi
  i_norm = stefan_pi*atmosphere%surface%temperature**4
  s => self%mem(:,:,:,self%idx%s,self%it,1)
  g1 = self%gamma-1.0
  rad2 = (5.*planet%radius)**2
  u_max = 0.0
  do iz=self%mesh(3)%li,self%mesh(3)%ui
  do iy=self%mesh(2)%li,self%mesh(2)%ui
    do ix=self%mesh(1)%li,self%mesh(1)%ui
      rp = sqrt(r2(ix,iy,iz))
      !---------------------------------------------------------------------
      ! Reduced intensity, from optical depth distance from the surface and
      ! radial attenuation
      !---------------------------------------------------------------------
      rt_int = i_norm
      if (use_attenuation) then
        tau(1) = lagrange%interpolate1d (rp, atmosphere%nr ,atmosphere%rr, atmosphere%tau)
        tau(2) = max(0.0,atmosphere%surface%tau-tau(1))
        rt_int = rt_int*exp(-tau(2))
      end if
      !---------------------------------------------------------------------
      ! Let the source term cool only towards the external temperature, and
      ! taper that effect off toward the embryo, to avoid edge heating
      !---------------------------------------------------------------------
      temp = self%d(ix,iy,iz)**g1*exp(s(ix,iy,iz)/self%d(ix,iy,iz)*g1)
      !temp_ext = atmosphere%external%temperature*min(1.0,max(0.0,0.1*(rp-planet%radius)/planet%radius))
      temp_ext = 0.0
      rt_src = stefan_pi*(temp-temp_ext)**4
      !---------------------------------------------------------------------
      ! Weight of point source relative to diffuse light, normalized to 1
      !---------------------------------------------------------------------
      if (eos%do_rt) then
        w_omega = 0.25*planet%radius**2/(r2(ix,iy,iz)-0.5*planet%radius**2)
      else
        w_omega = 1.0
      end if
      !---------------------------------------------------------------------
      ! Add the contributions.  The incoming q(ix,iy,iz) is from q_diffuse,
      ! and is already normalized with 4.*pi*rk.  It already contains a term
      ! -4.*pi*rk*source, but since we need to normalize the contribution 
      ! with (1.-w_omega), we need to add -w_omega*4.*pi*rk*src
      !---------------------------------------------------------------------
      rk = eos%kappa*self%d(ix,iy,iz)*4.*cgs%pi*epsilon
      w_z = exp(-rad2/max(rad2/log(1e9),r2(ix,iy,iz)))
      !if (iy==17.and.iz==17.and.self%id==285) print *, ix, i_norm, rt_int, rt_src, w_omega, w_z, &
      !  q(ix,iy,iz)*(1.-w_omega)*w_z + rk*(rt_int*w_omega-rt_src*(1.-w_omega)*(1.-w_z))
      q(ix,iy,iz) = q(ix,iy,iz)*(1.-w_omega)*w_z + rk*(rt_int*w_omega-rt_src*(1.-w_omega))*exp(-tau(1))*(1.-w_z)
      !---------------------------------------------------------------------
      ! Estimate the max rate of change, from assuming no cooling
      !---------------------------------------------------------------------
      if (rp > planet%radius-self%mesh(1)%d) then
        u_max = max(u_max,rk*rt_src/self%pg(ix,iy,iz)) !!!*w_omega
      end if
    end do
  end do
  end do
  u_max = u_max*self%mesh(1)%d*self%courant/self%cdtd
  if (self%boundaries%is_set(bits%spherical) .and. io%verbose>1) then
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

!> =============================================================================
!> Add accretion heating
!> =============================================================================
SUBROUTINE q_accretion (self, r2, q)
  class(experiment_t)                                 :: self
  real(kind=KindScalarVar), dimension(:,:,:)          :: r2, q
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: d
  real(kind=KindScalarVar)                            :: heat(2), pl_r2, d0
  real(8), dimension(:), allocatable                  :: lnrr
  real(8)                                             :: lnr, cx
  integer                                             :: ix, iy, iz, nr
  integer, save                                       :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%mdot>0.0 .or. disk_heat_cgs>0.0) then
    call trace%begin('experiment_t%q_accretion', itimer=itimer)
    !-----------------------------------------------------------------------------
    ! q is the heating rate from accretion, Q = Mdot G M / (4\pi r^4), in energy
    ! units per volume unit and time unit.
    !-----------------------------------------------------------------------------
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

!===============================================================================
!> Interface to patch_t%output, for debug printout
!===============================================================================
SUBROUTINE output (self)
  class(experiment_t):: self
  integer:: iv, it, m(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%output', itimer=itimer)
  if (self%id==io%id_debug) then
    m = self%mesh%gn/2
    it = self%iit(self%nt-2)
    do iv=2,4
      print *,'OUT: t, iv, u(z) =', real(self%t(it)), iv, &
        self%mem(m(1),m(2),m(3)-2:m(3)+2,iv,it,1)/ &
        self%mem(m(1),m(2),m(3)-2:m(3)+2, 1,it,1)
    end do
  end if
  !if (self%time > self%out_next) &
  !  call self%boundary_conditions (self%it)
  call self%pebble_patch_t%output
  call trace%end (itimer)
END SUBROUTINE

!===============================================================================
!> Boundary conditions for this experiment
!===============================================================================
SUBROUTINE boundary_conditions (self, it, new)
  class(experiment_t) :: self
  integer             :: it, new
  !.............................................................................
  real(8)             :: h(3)
  integer             :: iv, nb,i
  real, pointer       :: d(:,:,:), e(:,:,:), ss(:,:,:), p(:,:,:,:), ss2(:,:,:)
  real                :: s0
  integer, save       :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%boundaries%bits /= 0) then
    call trace%begin ('experiment_t%boundary_conditions', itimer=itimer)
    call self%check_density ('before boundary_conditions')
    if (self%boundaries%is_set (bits%spherical)) then
      if (io%verbose>1) print *, self%id, 'BOUNDARY CONDITION: spherical'
      if (io%verbose>0) print *, self%id, 'BOUNDARY CONDITION: spherical, t =', self%t(self%it)
      if (spherical_bc_corners) then
       !------------------------------------------------------------------------
       ! Jagged boundary conditions
       !------------------------------------------------------------------------
        call self%boundaries%spherical_jagged_d (self%mem(:,:,:,self%idx%d,:,1), self%mesh, self%idx%d, it, new)
        d => self%mem(:,:,:,self%idx%d,new,1)
        allocate (ss(self%gn(1),self%gn(2),self%gn(3)))
        allocate (ss2(self%gn(1),self%gn(2),self%gn(3)))
        if (self%kind(1:9)=='stagger2_') then
          ss = self%mem(:,:,:,self%idx%s,new,1)/d
          ss2= self%mem(:,:,:,self%idx%s,it,1)/self%mem(:,:,:,self%idx%d,it,1)
          call self%boundaries%spherical_jagged_s (ss, ss2, self%mesh)
          self%mem(:,:,:,self%idx%s,new,1) = ss*d
        else
          call self%e2ss (ss, new)
          call self%e2ss (ss2, it)
          call self%boundaries%spherical_jagged_s (ss, ss2,self%mesh)
          call self%ss2e (ss, new)
        end if
        deallocate (ss,ss2)
        do iv=self%idx%px,self%idx%pz
          call self%boundaries%spherical_jagged_p (self%mem(:,:,:,iv,new,1), self%mesh, iv)
        end do
      else
       !------------------------------------------------------------------------
       ! Spherical inner BCs on mass flux, which automatically makes density and
       ! entropy evolve consistently
       !------------------------------------------------------------------------
        p => self%mem(:,:,:,self%idx%px:self%idx%pz,new,1)
        call self%boundaries%spherical_p_ns (p, self%mesh, planet%position, planet%radius)
      end if
    end if
    !---------------------------------------------------------------------------
    ! Edge boundary conditions
    !---------------------------------------------------------------------------
    if (self%boundaries%bits /= 0) then
      if (io%verbose>2) then
        print *, self%id, 'MK x-BC:', self%boundaries%is_set (bits%xl), self%boundaries%is_set (bits%xu)
        print *, self%id, 'MK y-BC:', self%boundaries%is_set (bits%yl), self%boundaries%is_set (bits%yu)
        print *, self%id, 'MK z-BC:', self%boundaries%is_set (bits%zl), self%boundaries%is_set (bits%zu)
      end if
      call self%edge_boundary_conditions(new)
    end if
    call self%check_density ('after boundary_conditions')
    call trace%end (itimer)
  end if
END SUBROUTINE boundary_conditions

!===============================================================================
!> Boundary conditions at the edges of the outermost rubiks cubes; extrapolate
!> the density and pressure, and enforce vanishing derivative on the perpendicular
!> velocities. To what extend the boundary conditions are applied in the x,y,z
!> directions depends on the self%boundaries%bits bits. Note that the same
!> condition is applied in all directions that have boundary bits set!
!===============================================================================
SUBROUTINE edge_boundary_conditions (self, it)
  class(experiment_t):: self
  integer:: it
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e_th, uu, ss
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: u
  integer:: iv, select
  real:: ee
  integer, save       :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%edge_boundary_conditions', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allocate and load values into velocity, consistent with current density
  !-----------------------------------------------------------------------------
  d => self%mem(:,:,:,self%idx%d,it,1)
  allocate (u   (size(d,1),size(d,2),size(d,3),3))
  call self%p2U(u,it)
  !-----------------------------------------------------------------------------
  ! Density boundary conditions
  !-----------------------------------------------------------------------------
  call self%check_density ('before density_BC')
  call self%density_boundary_condition (it)
  !call self%boundaries%log_extrapolate (d, self%idx%d)
  call self%check_density ('after density_BC')
  !-----------------------------------------------------------------------------
  ! Velocity boundary conditions
  !-----------------------------------------------------------------------------
  do iv=1,3
    uu => u(:,:,:,iv)
    !---------------------------------------------------------------------------
    ! The azimuthal velocity should be safe to extrapolate in directions
    ! where it it parallel to the boundary
    !---------------------------------------------------------------------------
    if (iv==2) then
      call self%boundaries%extrapolate (uu, self%idx%px+iv-1)
    !---------------------------------------------------------------------------
    ! Use vanishing perpendicular derivative in other cases
    !---------------------------------------------------------------------------
    else
      call self%boundaries%symmetric (uu, self%idx%px+iv-1)
    end if
  end do
  call self%U2p(u,it)
  deallocate (u)
  !-----------------------------------------------------------------------------
  ! Energy BCs: Set thermal_energy per unit mass to initial value
  !-----------------------------------------------------------------------------
  allocate (e_th(size(d,1),size(d,2),size(d,3)  ))
  call self%e2E_th(e_th, it)
  call self%boundaries%constant (e_th, self%idx%e, atmosphere%external%energy)
  call self%E_th2e(e_th, it)
  deallocate (e_th)

  ! apply the "wave killing" boundaries to all variables near the x-boundaries;
  if (self%x_sponge > -1.0e-32) then
    call apply_sponge_boundaries (self)
  end if

  call trace%end (itimer)
END SUBROUTINE edge_boundary_conditions

!===============================================================================
!> Constrain the density to follow an anchored Gaussian profile
!===============================================================================
SUBROUTINE density_boundary_condition (self, it)
  class(experiment_t):: self
  integer:: it
  class(mesh_t), pointer, dimension(:) :: m
  integer:: iv, iz, iy, ix, lz, uz, ly, uy, lx, ux, ix0, iy0, iz0
  real(8) :: x, x0, z, z0, r, d1, d0, tt, c2
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, dd
  !-----------------------------------------------------------------------------
  ! Vertical edge boundary conditions
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%density_boundary_conditions')
  m => self%mesh
  d => self%mem(:,:,:,1,it,1)
  if (self%nw>1) dd => self%mem(:,:,:,self%idx%d,it,2)
  associate (r1=>m(1)%r, r3=>m(3)%r)
  tt = atmosphere%external%temperature
  if(self%boundaries%is_set(bits%xl)) then
    lx = m(1)%lb; ux = m(1)%lo
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%ub
    ix0 = m(1)%li
    x0 = m(1)%p + r1(ix0) - planet%position(1)
    call apply_x
  end if
  if (self%boundaries%is_set(bits%xu)) then
    lx = m(1)%uo; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%ub
    ix0 = m(1)%ui
    x0 = m(1)%p + r1(ix0) - planet%position(1)
    call apply_x
  end if
  if(self%boundaries%is_set(bits%yl)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%lo
    lz = m(3)%lb; uz = m(3)%ub
    iy0 = m(2)%li
    call apply_y
  end if
  if (self%boundaries%is_set(bits%yu)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%uo; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%ub
    iy0 = m(2)%ui
    call apply_y
  end if
  if(self%boundaries%is_set(bits%zl)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%lb; uz = m(3)%lo
    iz0 = m(3)%li
    z0 = m(3)%p + r3(iz0) - planet%position(3)
    call apply_z
  end if
  if (self%boundaries%is_set(bits%zu)) then
    lx = m(1)%lb; ux = m(1)%ub
    ly = m(2)%lb; uy = m(2)%ub
    lz = m(3)%uo; uz = m(3)%ub
    iz0 = m(3)%ui
    z0 = m(3)%p + r3(iz0) - planet%position(3)
    call apply_z
  end if
  !if (self%id==io%id_debug) &
  !  print '(a,2i5,2l4,3(2x,2i5))', 'MK', self%id, it, &
  !    self%boundaries%is_set(bits%zl), &
  !    self%boundaries%is_set(bits%zu), lx,ux, ly,uy, lz,uz
  call self%check_density ('after density BC')
  end associate
  call trace%end()
  return
contains
  subroutine apply_x
    if (self%boundaries%check_filled) call trace%begin ('apply_x')
    if (io%verbose>2) &
    print *, self%id, 'MK: x-density BC', lx, ux, ly, uy, lz, uz
    c2 = planet%omega(3)**2*planet%trap/tt
    associate (r1=>m(1)%r)
    do iz=lz,uz
      do iy=ly,uy
        do ix=lx,ux
          x = m(1)%p + r1(ix) - planet%position(1)
          d0 = self%mem(ix0,iy,iz,self%idx%d,it,1)
          self%mem(ix,iy,iz,self%idx%d,it,1) = d0*exp(-c2*(x**2-x0**2))
        end do
        if (self%boundaries%check_filled) then
          do ix=lx,ux
            self%boundaries%filled(ix,iy,iz,self%idx%d) = .true.
          end do
        end if
      end do
    end do
    end associate
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%idx%d,it,2) = 0.0
          end do
        end do
      end do
    end if
    if (self%boundaries%check_filled) call trace%end()
  end subroutine apply_x
  subroutine apply_y
    if (self%boundaries%check_filled) call trace%begin ('apply_y')
    if (io%verbose>2) &
      print *, self%id, 'MK: y-density BC', lx, ux, ly, uy, lz, uz
    do iz=lz,uz
      do iy=ly,uy
        do ix=lx,ux
          self%mem(ix,iy,iz,self%idx%d,it,1) = self%mem(ix,iy0,iz,self%idx%d,it,1)
        end do
        if (self%boundaries%check_filled) then
          do ix=lx,ux
            self%boundaries%filled(ix,iy,iz,self%idx%d) = .true.
          end do
        end if
      end do
    end do
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%idx%d,it,2) = 0.0
          end do
        end do
      end do
    end if
    if (self%boundaries%check_filled) call trace%end()
  end subroutine apply_y
  subroutine apply_z
    if (self%boundaries%check_filled) call trace%begin ('apply_z')
    if (io%verbose>2) &
    print *, self%id, 'MK: z-density BC', lx, ux, ly, uy, lz, uz
    c2 = 0.5*planet%omega(3)**2/tt
    associate (r3=>m(3)%r)
    do iz=lz,uz
      z = m(3)%p + r3(iz) - planet%position(3)
      do iy=ly,uy
        do ix=lx,ux
          d0 = self%mem(ix,iy,iz0,self%idx%d,it,1)
          d1 = d0*exp(-c2*(z**2-z0**2))
          if (self%id==io%id_debug.and.iy==uy/2.and.iz==uz) &
            print '(a,i4,1p,8e14.6)', &
            'BCd: ix,d0,p0/do,exp,do,d =',ix, d0, tt, (z-z0), &
            exp(-c2*(z**2-z0**2)), &
            self%mem(ix,iy,iz,self%idx%d,it,1), d1
          self%mem(ix,iy,iz,self%idx%d,it,1) = d1
        end do
        if (self%boundaries%check_filled) then
          do ix=lx,ux
            self%boundaries%filled(ix,iy,iz,self%idx%d) = .true.
          end do
        end if
      end do
    end do
    end associate
    if (self%nw>1) then
      do iz=lz,uz
        do iy=ly,uy
          do ix=lx,ux
            self%mem(ix,iy,iz,self%idx%d,it,2) = 0.0
          end do
        end do
     end do
    end if
    if (self%boundaries%check_filled) call trace%end()
  end subroutine apply_z
END SUBROUTINE density_boundary_condition

!===============================================================================
!> Apply the "wave killing" boundaries of de Val-Borro et al. (2006) to all
!> variables.
!> Currently implemented for x-boundaries only.
!===============================================================================
SUBROUTINE apply_sponge_boundaries (self)
  USE units_mod, only: cgs
  class(experiment_t):: self
  integer:: iv, i, j, k, eye
  real(kind=KindScalarVar):: x, z, xkill_lo, xkill_hi, xbndy_lo, xbndy_hi, tau_kill
  real(kind=KindScalarVar):: d0, d, p0, p, T0, entropy0, px0, py0, pz0, &
                             asymptote
  class(mesh_t), pointer:: m1, m2, m3
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: mem
  !-----------------------------------------------------------------------------
  if (.not. (self%boundaries%is_set(bits%xl) .or. self%boundaries%is_set(bits%xu))) &
    return
  call trace%begin ("experiment_t::apply_sponge_boundaries")

  if (io%debug(2) .or. self%track .and. io%verbose>0) print*, self%id, ": Wave killing ON"
  mem => self%mem(:,:,:,:,self%it,1)
  m1 => self%mesh(1)
  m2 => self%mesh(2)
  m3 => self%mesh(3)
  tau_kill = 2.0 * cgs%pi / planet%omega(3)
  d0 = atmosphere%external%density
  p0 = atmosphere%external%pressure
  T0 = atmosphere%external%temperature
  px0 = 0.0
  pz0 = 0.0
  !
  xbndy_lo = m1%p + m1%r(m1%li) - 0.5 * m1%d
  xkill_lo = xbndy_lo + self%x_sponge
  xbndy_hi = m1%p + m1%r(m1%ui) + 0.5 * m1%d
  xkill_hi = xbndy_hi - self%x_sponge
  !
  do iv=1,self%nv
    do k=m3%lb,m3%ub
      z = m3%p + m3%r(k)
      do j=m2%lb,m2%ub
        do i=m1%lb,m1%ub
          x = m1%p + m1%r(i)
          d = d0 * exp(-(planet%omega(3)**2 * 0.5 * (z**2 + 2.0 * planet%trap * (x - planet%position(1))**2) / T0))
          p = p0*(d/d0)
          if (iv == self%idx%d) then
            asymptote = d
          else if (iv == self%idx%py) then
            asymptote = -d * ((1.5+planet%trap)*planet%omega(3)*x - planet%v_head)
          else if (iv == self%idx%e) then
            if (self%kind(1:7)=='ramses_') then
              py0 = -d * ((1.5+planet%trap)*planet%omega(3)*x - planet%v_head)
              asymptote = p / (self%gamma - 1.0) + 0.5 * d * py0**2
            else if (self%kind(1:9)=='stagger2e') then
              asymptote = p / (self%gamma - 1.0)
            else if (self%kind(1:9)=='stagger2_') then
              asymptote = d * (log(p) - self%gamma * log(d)) / (self%gamma - 1.0)
            else if (self%kind(1:6)=='zeus3d') then
              asymptote = p / (self%gamma - 1.0)
            end if
          else if (iv == self%idx%px) then
            asymptote = px0
          else if (iv == self%idx%pz) then
            asymptote = pz0
          end if
          !
          ! apply damping procedure to each variable.
          if (self%boundaries%is_set(bits%xl) .and. (x <= xkill_lo)) then
            mem(i,j,k,iv) = mem(i,j,k,iv) + self%dtime &
                          * self%boundaries%sponge(mem(i,j,k,iv), asymptote, tau_kill, xkill_lo, xbndy_lo, x)
          else if (self%boundaries%is_set(bits%xu) .and. (x >= xkill_hi)) then
            mem(i,j,k,iv) = mem(i,j,k,iv) + self%dtime &
                          * self%boundaries%sponge(mem(i,j,k,iv), asymptote, tau_kill, xkill_hi, xbndy_hi, x)
          else if (self%boundaries%is_set(bits%xl) .and. (i <= m1%lo)) then
            ! reset ghost cell values; in the case of the x-velocity, make it anti-symmetric
            ! lower x-boundary
            if (iv == self%idx%px) then
              eye = 2*m1%lo - i + 1
              mem(i,j,k,iv) = -mem(eye,j,k,iv)
            else
              mem(i,j,k,iv) = asymptote
            end if
          else if (self%boundaries%is_set(bits%xu) .and. (i >= m1%uo)) then
            ! upper x-boundary
            if (iv == self%idx%px) then
              eye = 2*m1%uo - i - 1
              mem(i,j,k,iv) = -mem(eye,j,k,iv)
            else
              mem(i,j,k,iv) = asymptote
            end if
          end if
        end do
      end do
    end do
  end do
  call trace%end()
END SUBROUTINE apply_sponge_boundaries

!===============================================================================
!===============================================================================
FUNCTION is_relevant (self) RESULT (ok)
  class(experiment_t):: self
  logical:: ok
  real(8):: r(3), rmax
  integer:: ix, iy, iz
  !-----------------------------------------------------------------------------
  ok = .false.
  rmax = 0.0
  do iz=-1,1,2
    r(3) = self%position(3) -planet%position(3) + iz*self%size(3)*0.5_8
    do iy=-1,1,2
      r(2) = self%position(2) -planet%position(2) + iy*self%size(2)*0.5_8
      do ix=-1,1,2
        r(1) = self%position(1) -planet%position(1) + ix*self%size(1)*0.5_8
        rmax = max(rmax,sqrt(sum(r**2)))
        if (rmax > planet%radius) ok=.true.
      end do
    end do
  end do
  if (io%verbose>1) &
    print '("is_relevant:",i6,1p,2e11.3,l3)', self%level, rmax, planet%radius, ok
END FUNCTION is_relevant

END MODULE experiment_mod
