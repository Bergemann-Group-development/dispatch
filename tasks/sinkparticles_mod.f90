!> =============================================================================
!> A sink particle is represented as primarily a particle task, whose position
!> is updated with a particle solver.  In addition, it carries along a patch,
!> which is usually much smaller than an MHD patch, and whose values are
!> downloaded from contribution MHD patches.
!>
!> The test procedure is called from refine_t%check_sinkparticle, which
!> handles creating a task list link, adding an nbor list, and adding the
!> task to task list.
!> =============================================================================
MODULE sinkparticles_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE index_mod
  USE particle_history_mod
  USE download_mod
  USE data_hub_mod
  USE link_mod
  USE list_mod
  USE bits_mod
  USE scaling_mod
  implicit none
  private
  type, public, extends (patch_t):: sinkparticles_t
    real(8):: mass=0d0, mom(3)
    real(8):: acceleration(3)=0d0
    real(8), dimension(:,:)  , allocatable:: p, v, a
    real(4), dimension(:,:,:), allocatable:: r
    logical, dimension(:,:,:), allocatable:: mask
    logical:: do_sink
  contains
    procedure:: init
    procedure:: check
    procedure:: create
    procedure:: dnload
    procedure:: update
    procedure:: accrete
    procedure:: move
  end type
  !-----------------------------------------------------------------------------
  logical:: do_sink=.false.           ! turn feature on/off
  logical:: do_translate=.true.       ! make a transformation to the velocity of the sink
  logical:: do_refine=.true.          ! insist on refining down to levelmax at star
  logical:: only_high_levels=.true.   ! only accrete at high levels
  integer:: id=0                      ! sink id
  integer:: center_star=0             ! index of star to zoom in to
  integer:: verbose=2                 ! verbosity level
  integer:: max_sinks=10000           ! for testing
  real:: accretion_rate=0.01          ! The rate of accretion, measured in Omega_Kepler
  real:: accretion_efficiency=0.5     ! The accretion radius, measured in cells
  real:: accretion_fraction=1e-6      ! The accretion radius, measured in cells
  real:: accretion_radius=4.          ! The accretion radius, measured in cells
  real:: exclusion_radius=32.         ! The radius, measured in cells, inside which no other sink are born
  real:: jeans_resolution=2.          ! Minimum, umber of cells per Jeans' length
  real:: rho_limit=1e6                ! The limit above which sink_particles are created
  real:: rho_limit_factor=8.          ! Express rho_limit relative to Jeans' ladder top
  real:: rho_fractionR=1e-6           ! The fration of rho_limit above which accretion happens
  real:: rho_cut=1.0                  ! ??
  real:: softening_length=1.0         ! Softening of gravitational potential
  real(8):: out_time=0.1              ! cadence of output into sink.dat
  character(len=64):: file ='sink.dat'! name of sink data file
  !-----------------------------------------------------------------------------
  integer, save:: n_sinks=0           ! counter
  type(sinkparticles_t), public:: sinkparticles
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class (sinkparticles_t):: self
  !.............................................................................
  logical, save:: first_time=.true.
  integer:: iostat
  namelist /sinkparticles_params/ &
    do_sink, do_translate, do_refine, only_high_levels, center_star, verbose, &
    accretion_rate, accretion_efficiency, accretion_fraction, accretion_radius, &
    exclusion_radius, jeans_resolution, rho_limit, rho_limit_factor, rho_cut, &
    max_sinks, softening_length, out_time, file
  !-----------------------------------------------------------------------------
  call trace%begin ('sinkparticles_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, sinkparticles_params, iostat=iostat)
    if (io_unit%master) write (io_unit%output, sinkparticles_params)
  end if
  !$omp end critical (input_cr)
  self%do_sink = do_sink
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Test if a sink particle should be created in a given patch
!> =============================================================================
FUNCTION check (self, patch)
  class(sinkparticles_t):: self
  class(patch_t), target:: patch
  logical:: check
  !-----------------------------------------------------------------------------
  check = .false.
  if (n_sinks >= max_sinks) return
  call trace%begin ('sinkparticles_t%check')
  check = .true.
  check = check .and. patch%fmaxval(patch%idx%d) > rho_limit
  if (check) then
    print *, 'sinkparticles_t%check: dmax in', &
      maxloc (patch%mem(:,:,:,patch%idx%d,patch%it,1))
    !$omp atomic
    n_sinks = n_sinks+1
    print *,'n_sinks =', n_sinks
  end if
  call trace%end
END FUNCTION check

!> =============================================================================
!> Create a sink partícle 
!> =============================================================================
SUBROUTINE create (self, patch)
  class(sinkparticles_t):: self
  class(patch_t), target:: patch
  !.............................................................................
  class(link_t), pointer:: nbor, prev
  real, dimension(:,:,:), pointer:: d
  real(8), dimension(:), pointer:: r
  integer:: i, l(3), u(3), m(3), ix, iy, iz
  real:: r2, dmin
  !-----------------------------------------------------------------------------
  call trace%begin ('sinkparticles_t%create')
  allocate (self%p(3,self%nt))
  allocate (self%v(3,self%nt))
  allocate (self%a(3,self%nt))
  l = patch%mesh%li
  u = patch%mesh%ui
  d => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)
  m = maxloc(d)+l-1
  do i=1,3
    r => patch%mesh(i)%r
    self%position(i) = patch%position(i) + r(m(i))
    self%ds(i) = patch%mesh(i)%d
  end do
  self%n  = 2*ceiling(accretion_radius)
  self%ng = 2
  self%nv = patch%nv
  self%nt = patch%nt
  self%nw = 1
  !-----------------------------------------------------------------------------
  ! We should use no_mans_land=.false., and an odd number for self%n
  !-----------------------------------------------------------------------------
  self%no_mans_land = .false.
  self%n  = (self%n/2)*2 + 1
  self%size = self%ds*(self%n-1)
  !-----------------------------------------------------------------------------
  ! Set a high level on the sink particle patch, to force downloading of interior
  ! data
  !-----------------------------------------------------------------------------
  call self%patch_t%setup (self%position, self%size, self%n, self%ng, self%nt, self%nv, self%nw)
  if (self%id==io%id_debug) then
    print *, 'sinkparticles_t%create: id, patch%id, m, dmax', self%id, patch%id, m, maxval(d)
    print *, 'sinkparticles_t%create: (pos-patch%pos)/size', (self%position-patch%position)/patch%size
    print *, 'sinkparticles_t%create: patch%pos, d', patch%myposition(m), patch%mem(m(1),m(2),m(3),patch%idx%d,patch%it,1)
    print *, 'sinkparticles_t%create:  self%n   , patch%n   ', self%n   , patch%n
  end if
  self%quality = 1024
  self%level   = patch%level
  self%time    = patch%time
  do i=1,3
    self%mesh(i)%h  = patch%mesh(i)%h
  end do
  self%mass = 0d0
  self%velocity = 0d0
  self%unsigned(self%idx%d) = .true.
  call self%set (bits%internal)
  !-----------------------------------------------------------------------------
  ! Pre-load radius and mask array for use in accretion calculations
  !-----------------------------------------------------------------------------
  allocate (self%mask(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  allocate (self%r   (self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  do iz=1,self%mesh(3)%gn
  do iy=1,self%mesh(2)%gn
  do ix=1,self%mesh(1)%gn
    r2 = (ix-self%mesh(1)%o)**2 &
       + (iy-self%mesh(2)%o)**2 &
       + (iz-self%mesh(3)%o)**2
    self%r(ix,iy,iz) = max(sqrt(r2),1e-3)
    self%mask(ix,iy,iz) = merge(.true., .false., r2<accretion_radius**2)
  end do
  end do
  end do
  call trace%end
END SUBROUTINE create

!> =============================================================================
!> Download values to all the sinkparticle%patch cells, but only to the 
!> density and momentum slots.  This relies on that the nbor list is 
!> comprehensive, including also the parent patch
!> =============================================================================
SUBROUTINE dnload (self, only)
  class (sinkparticles_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('sinkparticles_t%dnload')
  call download%download_link (self%link, all_cells=.true., only=self%idx%d)
  call download%download_link (self%link, all_cells=.true., only=self%idx%px)
  call download%download_link (self%link, all_cells=.true., only=self%idx%py)
  call download%download_link (self%link, all_cells=.true., only=self%idx%pz)
  call trace%end
END SUBROUTINE dnload

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class (sinkparticles_t):: self
  call trace%begin ('sinkparticles_t%update')
  call self%accrete
  call self%move
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  call trace%end
END SUBROUTINE update

!> =============================================================================
!> =============================================================================
SUBROUTINE accrete (self)
  class (sinkparticles_t):: self
  !.............................................................................
  real, dimension(:,:,:), pointer:: d, px, py, pz
  real, dimension(:,:,:), allocatable:: prf
  real, dimension(:,:,:), allocatable:: u2
  real:: v2_kepler
  integer:: ix, iy, iz
  real:: omega, dv
  real(8):: p, q, dm
  !-----------------------------------------------------------------------------
  call trace%begin ('sinkparticles_t%accrete')
  allocate (prf(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  allocate (u2 (self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  !-----------------------------------------------------------------------------
  ! Short hand
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  u2 = (px/d)**2
  u2 = (py/d)**2 + u2
  u2 = (px/d)**2 + u2
  self%u_max = sqrt(self%fmaxval(u2))
  call self%courant_condition
  dm = 0d0
  dv = product(self%ds)
  !print *, 'XX0', minval(d), maxval(d), self%fmaxval(d)
  !print *, 'XX1', dv, self%mesh%n, self%mesh%gn
  if (self%mass==0d0) then
    do iz=self%mesh(3)%li,self%mesh(3)%ui
    do iy=self%mesh(2)%li,self%mesh(2)%ui
    do ix=self%mesh(1)%li,self%mesh(1)%ui
      q = min(rho_limit/d(ix,iy,iz),1d0)
      q = merge (q, 1d0, self%mask(ix,iy,iz))
      p = 1d0-q
      if (self%mask(ix,iy,iz)) then
        dm          = dm          + d (ix,iy,iz)*dv*p
        if (verbose>1) &
          print '(3i4,3f12.6,l4)', ix, iy, iz, d(ix,iy,iz), self%r(ix,iy,iz), p, self%mask(ix,iy,iz)
        self%mom(1) = self%mom(1) + px(ix,iy,iz)*dv*p
        self%mom(2) = self%mom(2) + py(ix,iy,iz)*dv*p
        self%mom(3) = self%mom(3) + pz(ix,iy,iz)*dv*p
        d (ix,iy,iz) = d (ix,iy,iz)*q
        px(ix,iy,iz) = px(ix,iy,iz)*q
        py(ix,iy,iz) = py(ix,iy,iz)*q
        pz(ix,iy,iz) = pz(ix,iy,iz)*q
      end if
    end do
    end do
    end do
    !print *, '1 id,dm:', self%id, dm
  !-----------------------------------------------------------------------------
  ! Remove a fraction p of per-volume quantities and keep q=(1-p)
  !-----------------------------------------------------------------------------
  else
    !print *, 'XX', scaling%grav, self%mass, minval(self%r), maxval(self%r)
    do iz=1,self%mesh(3)%gn
    do iy=1,self%mesh(2)%gn
    do ix=1,self%mesh(1)%gn
      omega = sqrt(scaling%grav*self%mass/self%r(ix,iy,iz))/self%r(ix,iy,iz)
      p = merge (accretion_rate*self%dtime*omega, 0d0, self%mask(ix,iy,iz))
      q = 1d0-p
      dm          = dm          + d (ix,iy,iz)*dv*p
      self%mom(1) = self%mom(1) + px(ix,iy,iz)*dv*p
      self%mom(2) = self%mom(2) + py(ix,iy,iz)*dv*p
      self%mom(3) = self%mom(3) + pz(ix,iy,iz)*dv*p
      d (ix,iy,iz) = d (ix,iy,iz)*q
      px(ix,iy,iz) = px(ix,iy,iz)*q
      py(ix,iy,iz) = py(ix,iy,iz)*q
      pz(ix,iy,iz) = pz(ix,iy,iz)*q
    end do
    end do
    end do
    !print *, '2 id,dm:', self%id, dm
  end if
  self%mass = self%mass + dm
  self%velocity = self%mom/self%mass
  if (verbose > 0) &
    print '(4(a,1p,e12.4,2x),2(a,3e12.4,2x))', &
      'sink: time=', self%time, 'mass=', self%mass, &
      'dmin=',self%fminval(self%idx%d), 'dmax=',self%fmaxval(self%idx%d), &
      'vel=',self%velocity, 'pos=',self%position
  call trace%end
END SUBROUTINE accrete

!> =============================================================================
!> Move the sink particle one time step self%dt, using kidk-drift-kick
!> =============================================================================
SUBROUTINE move (self)
  class (sinkparticles_t):: self
  real(8):: acceleration(3)
  real, dimension(:,:,:), pointer:: phi
  integer:: ix, iy, iz, ii(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sinkparticles_t%move')
  !-----------------------------------------------------------------------------
  ! Kick. The sink is assumed to be in the exact grid position floor(mesh%o)
  !-----------------------------------------------------------------------------
  !phi = self%mem(:,:,:,self%idx%phi,self%it,1)
  ii = floor(self%mesh%o)
  ix = ii(1); iy=ii(2); iz=ii(3)
  self%acceleration = 0.0
  !self%acceleration = [(phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  ))/(2.*self%ds(1)), &
  !                     (phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  ))/(2.*self%ds(2)), &
  !                     (phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1))/(2.*self%ds(3))]
  self%velocity = self%velocity + 0.5d0*self%acceleration*self%dtime
  !-----------------------------------------------------------------------------
  ! Drift + update phi to the new space-time position
  !-----------------------------------------------------------------------------
  self%time     = self%time     + self%dtime
  self%position = self%position + self%dtime*self%velocity    
  !call download%download_link (self%link, all_cells=.true., only=self%idx%phi)
  !-----------------------------------------------------------------------------
  ! Kick
  !-----------------------------------------------------------------------------
  self%acceleration = 0.0
  !self%acceleration = [(phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  ))/(2.*self%ds(1)), &
  !                     (phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  ))/(2.*self%ds(2)), &
  !                     (phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1))/(2.*self%ds(3))]
  self%velocity = self%velocity + 0.5d0*self%acceleration*self%dtime
  call trace%end
END SUBROUTINE move

END MODULE sinkparticles_mod
