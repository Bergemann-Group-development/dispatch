!==============================================================================
!> Particles data types for dust-gas interaction.  This module expects to have
!> the permanent particle proporties listed below defined from outside, as well
!> as the particle accelerations, gas density and velocity (and sound speed in
!> the supersonic Epstein regime case).  These prepatations should happen in the
!> (possibly extended) particle_solver_t data type.
!===============================================================================
!define DBG
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE timer_mod
  USE omp_lock_mod
  USE gpatch_mod
  USE link_mod
  USE eos_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! nbor particles
  !-----------------------------------------------------------------------------
  type, public:: pnbor_t
    type(pnbor_t), pointer:: next => null()
    type(particles_t), pointer:: imports
  end type
  !-----------------------------------------------------------------------------
  ! Arrays of particle attributes
  !-----------------------------------------------------------------------------
  type, public:: particles_t                               ! patch particles
    integer                         :: id                  ! patch ID
    integer                         :: lb(3), ub(3)        ! index bounds
    integer                         :: m = 0               ! max particles
    integer                         :: n = 0               ! number of particles
    integer                         :: l = 0               ! length of used part
    integer                         :: seek = 1            ! next empty slot
    integer                         :: offset(3)           ! index offsets
    real                            :: ds                  ! cell size
    real(8)                         :: t = 0d0             ! particle time
    real(8)                         :: momentum(3)         ! total momentum
    type(particles_t), pointer      :: exports => null()   ! particles to export
    type(particles_t), pointer      :: imports => null()   ! particles to import
    class(pnbor_t), pointer         :: pnbors  => null()   ! nbor list
    class(lock_t), pointer          :: lock    => null()   ! OMP lock
    class(gpatch_t), pointer        :: patch   => null()   ! patch pointer
    !---------------------------------------------------------------------------
    ! Permanent particle arrays
    !---------------------------------------------------------------------------
    real, pointer, contiguous       :: p(:,:)  => null()   ! float position 12b
    real, pointer, contiguous       :: v(:,:)  => null()   ! velocity       12b
    integer(2), pointer, contiguous :: q(:,:)  => null()   ! index position  6b
    integer(4), pointer, contiguous :: i(:)    => null()   ! index position  4b
    real, pointer, contiguous       :: w(:)    => null()   ! weight          4b
    real, pointer, contiguous       :: s(:)    => null()   ! size            4b
    !                                                 total bytes/particle: 42b
    !---------------------------------------------------------------------------
    ! Temporary particle arrays
    !---------------------------------------------------------------------------
    real, pointer, contiguous       :: a(:,:)      => null()  ! acceleration
    real, pointer, contiguous       :: c(:)        => null()  ! sound speed
    real, pointer, contiguous       :: dgas(:)     => null()  ! gas dens
    real, pointer, contiguous       :: vgas(:,:)   => null()  ! gas vel
    !---------------------------------------------------------------------------
    ! Temporary field arrays
    !---------------------------------------------------------------------------
    real, pointer, contiguous       :: dg(:,:,:)   => null()  ! gas density
    real, pointer, contiguous       :: eg(:,:,:)   => null()  ! gas internal en
    real, pointer, contiguous       :: cs(:,:,:)   => null()  ! sound speed
    real, pointer                   :: pg(:,:,:,:) => null()  ! gas momentum
    real, pointer, contiguous       :: vg(:,:,:,:) => null()  ! gas velocity
    real, pointer, contiguous       :: gf(:,:,:,:) => null()  ! generic force
    !---------------------------------------------------------------------------
    ! Cell-averaged particle values
    !---------------------------------------------------------------------------
    real,       pointer, contiguous :: dd(:,:,:)   => null()  ! dust density
    real,       pointer, contiguous :: pd(:,:,:,:) => null()  ! dust momentum
    integer(2), pointer, contiguous :: pc(:,:,:)   => null()  ! particle count
    integer   , pointer, contiguous :: ic(:,:,:)   => null()  ! cell indices
    !---------------------------------------------------------------------------
    ! Particle attributes
    !---------------------------------------------------------------------------
    integer:: verbose = 0                       ! noise level
    integer:: per_cell = 50                     ! particles_per_cell
    integer:: seed = 11                         ! random seed
    integer:: iout = 0                          ! output snapshot
    integer:: ip_dbg = 0                        ! particle ID to debug
    integer:: id_dbg = 0                        ! patch ID to debug
    integer:: i_dbg(5) = 0                      ! cell index triplet to debug
    logical:: on = .false.                      ! on or not
    logical:: backreaction = .true.             ! onto gas?
    logical:: supersonic = .false.              ! supersonic Epstein?
    logical:: stokes = .false.                  ! Stokes number instead of size
    integer:: ppc = 50, ppc_min = 25            ! particles-per-cell
    real:: solid_density                        ! dust solid density
    real:: dust_to_gas = 0.01                   ! dust-to-gas ratio
    real:: alpha_exponent = -1.0                ! dust power law distribution
    real:: min_stokes = 0.1                     ! minimum Stokes number
    real:: max_stokes = 0.1                     ! maximum Stokes number
    real:: tolerance = 1e-4                     ! in conserving momentum
    real(8):: out_next = 0.0                    ! output cadence
    real(8):: out_time = -1.0                   ! output cadence
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_light
    procedure:: initial_condition
    procedure:: init_indices
    procedure:: pre_update
    procedure:: post_update
    procedure:: update
    procedure:: export
    procedure:: import
    procedure:: reset
    procedure:: move
    procedure:: cell_averages
    procedure:: renormalize
  end type
  real:: density_yr_cgs = 3e-12               ! drag scaling constant
  real:: solid_density_cgs = 2                ! dust solid density
  real:: dust_to_gas = 0.01                   ! dust-to-gas ratio
  real:: max_size_cgs = 1e-1                  ! max particle size
  real:: min_size_cgs = 1e-4                  ! min particle size
  real:: max_init_size_cgs = 0.0              ! single initial size
  real:: alpha_exponent = -3.5                ! dust power law distribution
  real:: min_fraction = 0.7                   ! lower limit
  real:: max_fraction = 1.4                   ! upper limit
  real:: drag_coeff                           ! scaled value
  real:: v_stream = 10                        ! Two-stream speed
  real:: min_stokes = 0.1                     ! minimum Stokes number
  real:: max_stokes = 0.1                     ! maximum Stokes number
  real(8):: out_time = -1.0                   ! output cadence
  integer:: seed = 11                         ! random seed
  integer:: test_nr = 0                       ! test nr (0=off)
  integer:: output_every=1000                 ! output only that fraction
  integer:: ip_dbg = 0                        ! particle ID to debug
  integer:: id_dbg = 0                        ! patch ID to debug
  integer:: i_dbg(5) = 0                      ! cell index triplet to debug
  logical:: particle_io = .false.             ! output particle positions?
  logical:: backreaction = .true.             ! backreaction on gas
  logical:: detailed_timer = .false.          ! precision timing
  logical:: do_supersonic = .true.            ! supersonic dust friction
  logical:: debug=.false.                     ! check things
  logical:: on = .false.                      ! on or not
  logical:: stokes = .false.                  ! Stokes number instead of size
  integer:: per_cell = 50                     ! particles per cell
  real:: tolerance = 1e-4
  integer:: verbose=1
CONTAINS

!===============================================================================
!> Initialize particles, reading parameters and adding task link.  Having
!> patch as a parameter allows the module to initialize things independently
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(particles_t), target  :: self
  class(gpatch_t), pointer    :: patch
  !.............................................................................
  integer                     :: iostat
  real                        :: solid_density, dln_s
  namelist /particles_params/ &
    on, verbose, detailed_timer, particle_io, out_time, output_every, &
    per_cell, seed, ip_dbg, id_dbg, i_dbg, min_fraction, max_fraction, &
    max_size_cgs, min_size_cgs, density_yr_cgs, solid_density_cgs, &
    dust_to_gas, debug, alpha_exponent, stokes, min_stokes, max_stokes
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%pre_init')
  self%patch => patch                           ! make the patch known
  self%id = patch%id                            ! for now, inherit
  !-----------------------------------------------------------------------------
  allocate (self%lock)
  call self%lock%init ('ppart')
  call io%lock%set ('particles_t%pre_init')
  if (first_time) then
    rewind (io%input)
    read (io%input, particles_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, particles_params)
    end if
    call io%namelist_warning ('particles_params', iostat)
    write (io%namelist, particles_params)
    if (max_init_size_cgs == 0.0) &
      max_init_size_cgs = max_size_cgs
    first_time = .false.
  end if
  call io%lock%unset ('particles_t%pre_init')
  !-----------------------------------------------------------------------------
  ! The might be useful in the particle_solver
  !-----------------------------------------------------------------------------
  self%dust_to_gas = dust_to_gas
  self%per_cell = per_cell
  self%out_next = 0d0
  self%out_time = out_time
  self%stokes = stokes
  self%min_stokes = min_stokes
  self%max_stokes = max_stokes
  self%iout = 0
  self%on = on
  !---------------------------------------------------------------------------
  ! Pre-factor for drag force.  The input parameter density_yr_cgs is the
  ! density (in CGS) for which the stopping time in 1 yr for a 1 cm particle.
  ! Since we keep size in cm, we just need to multiply with the density
  ! scaling constant, and divide with the value of one year in code units.
  !---------------------------------------------------------------------------
  drag_coeff = scaling%d/(density_yr_cgs*cgs%yr/scaling%t)
  solid_density = solid_density_cgs/(scaling%m/scaling%l**3)
#define BACKREACTION
#ifdef BACKREACTION
  self%backreaction = .true.
#else
  self%backreaction = .false.
#endif
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize, allocating arrays for particles and their cell attributes
!===============================================================================
SUBROUTINE init (self, patch, m)
  class(particles_t):: self
  class(gpatch_t), pointer:: patch
  integer:: m, l(3), u(3)
  !.............................................................................
  if (m==0) return
  call trace%begin ('particles_t%init')
  call self%init_light (patch, m)
  !-----------------------------------------------------------------------------
  ! Allocate cell arrays
  !-----------------------------------------------------------------------------
  l = self%lb
  u = self%ub
  allocate (self%dd(  l(1):u(1),l(2):u(2),l(3):u(3)))     ! dust density
  allocate (self%pd(3,l(1):u(1),l(2):u(2),l(3):u(3)))     ! dust momentum
  allocate (self%vg(3,l(1):u(1),l(2):u(2),l(3):u(3)))     ! gas velocity
  allocate (self%pc(  l(1):u(1),l(2):u(2),l(3):u(3)))     ! particles/cell
  allocate (self%ic(  l(1):u(1),l(2):u(2),l(3):u(3)))     ! start indices
  call io%bits_mem (storage_size(self%dd), product(shape(self%dd)), 'dd')
  call io%bits_mem (storage_size(self%pd), product(shape(self%pd)), 'pd')
  call io%bits_mem (storage_size(self%vg), product(shape(self%vg)), 'vg')
  call io%bits_mem (storage_size(self%pc), product(shape(self%pc)), 'pc')
  call io%bits_mem (storage_size(self%ic), product(shape(self%ic)), 'ic')
  self%dd = 1.0
  self%pd = 0.0
  if (associated(self%gf)) then
    self%gf = 0.0
  end if
  call self%init_indices ()
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize, allocating arrays for particles and their cell attributes
!===============================================================================
SUBROUTINE init_light (self, patch, n)
  class(particles_t):: self
  class(gpatch_t), pointer:: patch
  integer:: n
  integer:: l(3), u(3), iostat, m
  !.............................................................................
  if (n==0) return
  call trace%begin ('particles_t%init_light')
  self%m = n
  self%n = 0
  self%seek = 1
  self%patch => patch
  self%id = patch%id
  self%t  = patch%time
  if (.not. associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init
  endif
  self%offset = patch%n * patch%key(2:4)
  !-----------------------------------------------------------------------------
  ! Allocate particle arrays, plus a margin = size of exports
  !-----------------------------------------------------------------------------
  m = self%m
  allocate (self%w(m), self%s(m), self%i(m))
  allocate (self%p(m,3), self%v(m,3), self%q(m,3))
  self%v = 0.0
  self%w = 0.0
  self%s = 1.0
  self%q = self%lb(1)
  call trace%end ()
END SUBROUTINE init_light

!==============================================================================
!> Deallocate particles arrays and reset counters
!===============================================================================
SUBROUTINE reset (self)
  class(particles_t):: self
  self%seek = 1
  self%n = 0
  self%m = 0
  if (associated(self%a)) deallocate (self%a)
  if (associated(self%p)) deallocate (self%p)
  if (associated(self%q)) deallocate (self%q)
  if (associated(self%w)) deallocate (self%w)
  if (associated(self%dgas)) deallocate (self%dgas)
  if (associated(self%vgas)) deallocate (self%vgas)
  if (associated(self%dd)) deallocate (self%dd)
  if (associated(self%pd)) deallocate (self%pd)
  if (associated(self%vg)) deallocate (self%vg)
  if (associated(self%pc)) deallocate (self%pc)
  if (associated(self%ic)) deallocate (self%ic)
END SUBROUTINE reset

!===============================================================================
!> Initialize particles
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particles_t)        :: self
  class(gpatch_t),  pointer :: patch
  integer                   :: m(3), i1, i2, i3, i, j, n, nv, l, nb, id(5)
  real                      :: d, w, w1, v(3), chance, lns_min, lns_max
  real                      :: ran_u, wtot
  logical                   :: draw
  real, contiguous, pointer :: s(:)                 ! particle size
  integer, save             :: lid = 0
  integer, save             :: itimer=0, extra_bits
  !$omp threadprivate (lid)
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%initial_condition', itimer=itimer)
  !-----------------------------------------------------------------------------
  patch => self%patch
  !-----------------------------------------------------------------------------
  ! Initialize dust.  At boundaries, add particles also in guard zones
  !-----------------------------------------------------------------------------
#ifdef __GFORTRAN__
  !$omp critical (alloc_cr)
#endif
  lns_min = log(min_size_cgs)
  lns_max = log(max_size_cgs)
  self%t = 0
  i = 0
  self%dg => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  self%pg => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  allocate (s(per_cell))
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
    n = 0
    do i1=self%lb(1),self%ub(1)
      !-------------------------------------------------------------------------
      ! Compute normalizatoin factor in sum over particles in cell
      !-------------------------------------------------------------------------
      wtot=0
      do j=1,per_cell
        ran_u = patch%random%ran1()
        s(j) = exp(lns_min + ran_u*(lns_max-lns_min))
        wtot = wtot+s(j)**(alpha_exponent+4.0)
      end do
      d = self%dg(i1,i2,i3)
      v = self%pg(i1,i2,i3,:)/(patch%ds*d)
      d = d*dust_to_gas
      !-------------------------------------------------------------------------
      ! Apply factors to get sum(self%w) = d in cell
      !-------------------------------------------------------------------------
      do j=1,per_cell
        i = i+1
        self%i(i)   = patch%id*self%m + i
        self%q(i,:) = [i1,i2,i3]
        self%p(i,:) = patch%random%uniform (3)
        self%v(i,:) = v
        self%w(i)   = s(j)**(alpha_exponent+4.0)/wtot * d
        self%s(i)   = s(j)/scaling%l
        n = n+1
      end do
    end do
  end do
  end do
  self%n = i
  self%l = i
  deallocate(s)
#ifdef __GFORTRAN__
  !$omp end critical (alloc_cr)
#endif
  call trace%end (itimer)
END SUBROUTINE initial_condition

!===============================================================================
!> Initialize a counter of new particles per cell, and the first address before
!> the paricle array for each cell
!===============================================================================
SUBROUTINE init_indices (self)
  class(particles_t):: self
  integer:: i1, i2, i3, n(3), l(3), u(3)
  !-----------------------------------------------------------------------------
  l = self%lb
  u = self%ub
  n = u-l+1
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    self%pc(i1,i2,i3) = 0
    self%ic(i1,i2,i3) = i1-l(1) + (i2-l(2) + (i3-l(3))*n(2))*n(1)
    self%ic(i1,i2,i3) = self%ic(i1,i2,i3)*per_cell
  end do
  end do
  end do
END SUBROUTINE init_indices

!===============================================================================
!> Prepare for update
!===============================================================================
SUBROUTINE pre_update (self)
  class(particles_t)       :: self
  class(gpatch_t), pointer :: patch
  real(8)                  :: ds
  integer                  :: i, l(3), u(3)
  integer,            save :: itimer=0
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%pre_update', itimer=itimer)
  patch => self%patch
  self%dg => patch%mem(:,:,:,patch%idx%d,patch%it,1)
  self%eg => patch%mem(:,:,:,patch%idx%e,patch%it,1)
  self%pg => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  ds = patch%ds(1)
  if (first_time) call timer%tic
  !.............................................................................
  if (first_time) then
    call self%lock%set
    if (first_time) then
      first_time = .false.
      call io%print_hl
      call timer%toc ('computing vgas', self%n)
      call io%print_hl
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Allocate temporary particle arrays
  !-----------------------------------------------------------------------------
  l = self%lb
  u = self%ub
  allocate (self%c(self%m), self%dgas(self%m))
  allocate (self%a(self%m,3), self%vgas(self%m,3))
  call io%bits_mem (storage_size(self%c)   , product(shape(self%c)), 'c')
  call io%bits_mem (storage_size(self%a)   , product(shape(self%a)), 'c')
  call io%bits_mem (storage_size(self%dgas), product(shape(self%dgas)), 'c')
  call io%bits_mem (storage_size(self%vgas), product(shape(self%vgas)), 'c')
  self%dgas = 1.0
  self%vgas = 0.0
  !-----------------------------------------------------------------------------
  if (self%supersonic) then
    allocate (self%cs(l(1):u(1),l(2):u(2),l(3):u(3)))     ! gas sound speed
  end if
  call trace%end (itimer)
END SUBROUTINE pre_update

!===============================================================================
!> Prepare for update
!===============================================================================
SUBROUTINE post_update (self)
  class(particles_t)       :: self
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%post_update')
  call io%bits_mem (-storage_size(self%c)   , product(shape(self%c)), 'c')
  call io%bits_mem (-storage_size(self%a)   , product(shape(self%a)), 'c')
  call io%bits_mem (-storage_size(self%dgas), product(shape(self%dgas)), 'c')
  call io%bits_mem (-storage_size(self%vgas), product(shape(self%vgas)), 'c')
  deallocate (self%c, self%dgas)
  deallocate (self%a, self%vgas)
  if (self%supersonic) then
    deallocate (self%cs)
  end if
  !-----------------------------------------------------------------------------
  ! Reset the exports buffer to the empty state
  !-----------------------------------------------------------------------------
  call self%exports%reset ()
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Particle mover, unrolled
!===============================================================================
SUBROUTINE update (self, dtime)
  class(particles_t)    :: self
  real                  :: dtime, dt
  integer               :: i, n, ne
  real                  :: d, e, cs, pc0, pcm, pwd, pc, w
  integer, dimension(3) :: q, ii, o, l, u
  real,    dimension(3) :: vg, vp, pq, pa, dv, gp, fd, pd, vd
  type(particles_t)     :: pt
  logical               :: outside
  integer,         save :: itimer=0
  logical,         save :: first_time=.true.
!-------------------------------------------------------------------------------
  call trace%begin ('particles_t%update', itimer=itimer)
  call self%import ()
  n = self%n
  !-----------------------------------------------------------------------------
  if (first_time) call timer%tic
  ne = 0
  l = self%lb
  u = self%ub
  o = self%patch%mesh%o
  self%momentum(:) = 0d0
  self%dg => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  do i = 1,self%l
    q  = self%q(i,:)                                    ! initical cell indices
    vg = self%vgas(i,:)                                 ! gas velocity
    d  = self%dgas(i)                                   ! gas density
    w  = self%w(i)                                      ! gas density
    vd = self%v(i,:)                                    ! dust velocity
    pd = self%p(i,:)                                    ! dust position
    dt = merge (0.0, dtime, w<0.0)                      ! freeze holes
    pc0 = drag_coeff*d/self%s(i)                        ! drag coeff size s
    pq = (pd + q - o)*self%ds                           ! physical position
    !---------------------------------------------------------------------------
    ! FIXME: the particle forces should be preloaded, in a separate loop
    !---------------------------------------------------------------------------
#ifdef SUPERSONIC
    cs = self%c(i)                                      ! sound speed
    dv = vg-vd                                          ! speed difference
    pc = pc0*sqrt(1.0 + sum(dv**2)/cs**2)               ! supersonic Epstein
    fd = pc*dv                                          ! drag force
#else
    pc = pc0                                            ! constant Epstein
    fd = pc*(vg-vd)                                     ! drag force
#endif
    pwd = self%w(i)/d                                   ! density ratio
    !---------------------------------------------------------------------------
    ! To allow extreme dust:gas ratios, compute the largest action first
    !---------------------------------------------------------------------------
    pa = self%a(i,:)                                    ! acceleration
#ifdef BACKREACTION
    dv = dt*(pa+fd)/(1.+pc*dt*(1.+pwd))                 ! velocity delta
    vd = vd + dv                                        ! particle update
    vg = vg - dv*pwd                                    ! gas update
    pd = pd + 2.*dt*vd                                  ! forward two dr
    fd = pc*(vg-vd)                                     ! drag force
    dv = dt*(pa+fd)/(1.+pc*dt*(1.+pwd))                 ! velocity delta
    vd = vd + dv                                        ! particle update
    vg = vg - dv*pwd                                    ! gas ufdate
#else
    vd = vd + dt*(pa+fd)/(1.+pc*dt)                     ! particle update
    pd = pd + 2.*dt*vd                                  ! forward two dr
    self%t = self%patch%time
    fd = pc*(vg-vd)                                     ! drag force
    vd = vd + dt*(pa+fd)/(1.+pc*dt)                     ! particle update
#endif
    !---------------------------------------------------------------------------
    ! Update gas velocity, q+r representation and patch particle momentum
    !---------------------------------------------------------------------------
    self%vg(:,q(1),q(2),q(3)) = vg                      ! save new gas v
    q = floor(pd)                                       ! integer shift
    pd = pd - q                                         ! subtract
    q = q + self%q(i,:)                                 ! new particle 3D index
    outside = any(q < l .or. q > u)                     ! detect exports
    w  = merge(-w,w,outside)                            ! mark weights
    ne = merge(ne+1,ne,outside)                         ! count number
    self%p(i,:) = pd                                    ! update position
    self%v(i,:) = vd                                    ! update velocity
    self%q(i,:) = q                                     ! update indices
    self%w(i)   = w                                     ! udpate weight
if (verbose>0 .and. outside) PRINT *,ne,i,self%w(i),self%q(i,:)
  end do
  self%exports%n = ne
  !.............................................................................
  if (first_time) then
    call io%lock%set
    if (first_time) then
      call io%print_hl
      call timer%toc ('update', self%n)
      call io%print_hl
      first_time = .false.
    end if
    call io%lock%unset
  end if
  call self%export ()
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Copy all particles that have left the patch to the import buffers of nbors
!===============================================================================
SUBROUTINE export (self)
  class(particles_t):: self
  class(gpatch_t), pointer:: patch, nbpatch
  class(pnbor_t), pointer:: pnbor
  class(particles_t), pointer:: to, fr, exports
  class(nbor_t), pointer:: nbor
  integer:: i, offset(3), nbox(3), l(3), u(3), q(3), n, j
  !-----------------------------------------------------------------------------
  if (self%exports%n < 1) &
    return
  call trace%begin ('particles_t%export')
  patch => self%patch
  exports => self%exports
if (verbose>0) PRINT *,exports%n,' particles in exports',exports%id
  call exports%init_light (patch, exports%n)
  exports%id = self%id
  self%imports%id = self%id
  !-----------------------------------------------------------------------------
  ! Move particles that fall outside our region of authority to exports
  !-----------------------------------------------------------------------------
  l = self%lb
  u = self%ub
  n = 0
  do i = 1,self%m
    if (self%w(i) < 0.0) then
      n = n+1
if (verbose>0) PRINT *,i,self%w(i),self%q(i,:)
    end if
  end do
if (verbose>0) PRINT *,n,' particles with negative w',exports%id
  n = 0
  j = 0
  do i = 1,self%m
    if (self%w(i) < 0.0) then
      self%w(i) = -self%w(i)
      n = n+1
      j = i
      call self%move (exports, i)
    end if
  end do
if (verbose>0) PRINT *,exports%n,' particles in exports',exports%id, n, j
  !-----------------------------------------------------------------------------
  ! Loop over nbor patches and copy exports particles that fall inside an nbor
  ! patch to its imports buffer
  !-----------------------------------------------------------------------------
  nbox = nint(patch%box/patch%ds)
  pnbor => self%pnbors
  do while(associated(pnbor))    
    to => pnbor%imports
    offset = exports%offset - to%offset
                                                     call to%lock%set ('export')
    do i=1,exports%n
      q = exports%q(i,:)
      q = modulo (q + offset, nbox)
      if (all(q >= l .and. q <= u)) then
        call exports%move (to, i)
      end if
    end do
if (verbose>0) PRINT *,to%n,' particles in imports to', to%id
                                                   call to%lock%unset ('export')
    pnbor => pnbor%next
  end do
if (verbose>0) PRINT *,exports%n,' particles left in exports before BCs'
  call trace%end ()
END SUBROUTINE export

!===============================================================================
!> Copy all particles that have been left in the import buffer
!===============================================================================
SUBROUTINE boundary_conditions (self, ib)
  class(particles_t):: self
  integer:: ib(3)
END SUBROUTINE boundary_conditions

!===============================================================================
!> Copy all particles that have been left in the import buffer
!===============================================================================
SUBROUTINE import (self)
  class(particles_t):: self
  class(particles_t), pointer:: fr
  integer:: i, j, n
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%import')
  fr => self%imports
                                                     call fr%lock%set ('export')
  i = self%seek
  do j = 1,fr%n
    call fr%move (self, i)
  end do
  self%seek = i
                                                   call fr%lock%unset ('export')
  call trace%end ()
END SUBROUTINE import

!==============================================================================
!> Copy particle i from self to recv and increment n, if the slot is empty
!===============================================================================
SUBROUTINE move (self, recv, j)
  class(particles_t):: self, recv
  integer:: j, i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%n < 1) return
  call trace%begin ('particles_t%move', itimer=itimer)
  i = recv%seek
  do while (recv%w(i) > 0.0)
    i = i + 1    
  end do
  if (i <= recv%m) then
    recv%p(i,:) = self%p(j,:)
    recv%q(i,:) = self%q(j,:)
    recv%v(i,:) = self%v(j,:)
    recv%w(i)   = self%w(j)
    recv%s(i)   = self%s(j)
    recv%i(i)   = self%i(j)
    self%w(j)   = 0
    self%q(j,:) = self%lb
    recv%seek = i
    recv%l = max(i,recv%l)
    recv%n = recv%n + 1
    self%n = self%n - 1
  else
    write (stderr,*) 'particles_t%move: out of space'
  end if
  call trace%end (itimer)
END SUBROUTINE move

!==============================================================================
!> Pack the self pool, by removing slots with zero weight
!===============================================================================
SUBROUTINE cell_averages (self)
  class(particles_t):: self
  integer:: i, q(3), l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%cell_averages', itimer=itimer)
  l = self%patch%li
  u = self%patch%ui
  self%dd(  l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  self%pd(:,l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
  self%pc(  l(1):u(1),l(2):u(2),l(3):u(3)) = 0
  do i=1,self%n
    q = self%q(i,:)
    !if (all(q >= l .and. q <=u)) then
      self%dd(  q(1),q(2),q(3)) = &
      self%dd(  q(1),q(2),q(3)) + self%w(i)
      self%pd(:,q(1),q(2),q(3)) = &
      self%pd(:,q(1),q(2),q(3)) + self%w(i)*self%v(i,:)
      self%pc  (q(1),q(2),q(3)) = &
      self%pc  (q(1),q(2),q(3)) + 1
    !else
    !  write (stderr,*) 'out-of-bounds particle: id,q=', self%patch%id, q
    !end if
  end do
  call trace%end (itimer)
END SUBROUTINE cell_averages

!===============================================================================
!> Strategy: Copy particles from the existing arrays to new arrays, then
!>           deallocate the old arrays and stich in the new ones.
!>
!>   Method: Use two sweeps. The first one measures density, momentum, and
!>           number of particles per cell in the current state. The second
!>           one adds new particles in cells that have to few, in a way that
!>           maintains the distribution properties.  A third sweep copies
!>           particles directly in place, with weights adjusted to accomodate
!>           the added particles or removed particles, as needed.
!>
!> When removing particles, their mass and momentum are added to the remaining
!> ones, with the choice of which particle to give the mass and momentum to
!> dictated by their existing mass. To maintain the mass distribution, the
!> remaining particle that is most similar to the removed one is given a new
!> mass that is the weighted average of the two it now represents.
!>
!> When adding particles, the new ones are given the average density and
!> momentum, and their size (or mass) is drawn from the local sample of
!> existing ones, picked by random choice.
!===============================================================================
SUBROUTINE renormalize (self)
  class(particles_t):: self
  type(particles_t):: new
  real:: dd, pd
  class(gpatch_t), pointer:: patch
  integer:: i, n, q(3), l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particles_t%renormalize', itimer=itimer)
  patch => self%patch
  l = patch%li
  u = patch%ui
  !-----------------------------------------------------------------------------
  ! Allocate a new particles_t instance
  !-----------------------------------------------------------------------------
  call new%init (self%patch, self%n)
  !-----------------------------------------------------------------------------
  ! First sweep, computing mass, momentum, and particle-per-cell
  !-----------------------------------------------------------------------------
  call self%cell_averages ()  
  !-----------------------------------------------------------------------------
  ! Second sweep, copying over the particles that can be there
  !-----------------------------------------------------------------------------
  call self%init_indices ()
  n = self%n
  do i=1,n
    q = self%q(i,:)
    if (new%pc(q(1),q(2),q(3)) < per_cell*min_fraction) then
      call copy1
    end if
  end do
  !-----------------------------------------------------------------------------
  ! Third sweep
  !-----------------------------------------------------------------------------
  do i=1,n
    q = self%q(i,:)
    if (new%pc(q(1),q(2),q(3)) < per_cell*min_fraction) then
      call add (i)
    else if (new%pc(q(1),q(2),q(3)) == per_cell*min_fraction) then
      call remove (i)
    end if
  end do
  call trace%end (itimer)
contains
  !=============================================================================
  ! Copy a particle over, and increment the npc counter
  !=============================================================================
  subroutine copy1
    integer:: j
    j = self%ic(q(1),q(2),q(3)) + 1
    new%p(:,j) = self%p(i,:)
    new%v(:,j) = self%v(i,:)
    new%q(:,j) = self%q(i,:)
    new%w(j)   = self%w(i)
    new%s(j)   = self%s(i)
    new%i(j)   = self%i(i)
    self%ic(q(1),q(2),q(3)) = j
  end subroutine copy1
  !=============================================================================
  ! Remove a particle, adding its mass and momentum to existing ones
  !=============================================================================
  subroutine remove (i)
    integer:: i
    integer,    pointer:: ic_p
    integer(2), pointer:: nc_p
    integer:: j, j0, j1
    real:: f, dd0, dd1, pd1(3), dv(3)
    !---------------------------------------------------------------------------
    ic_p => new%ic(q(1),q(2),q(3))      ! pointer to particle index
    nc_p => new%pc(q(1),q(2),q(3))      ! pointer to particle count
    j1      = ic_p                      ! last particle index
    j0      = j1 - nc_p + 1             ! first particle index
    !---------------------------------------------------------------------------
    ! Loop over exiting particles and sum up their mass and momentum
    !---------------------------------------------------------------------------
    dd1 = 0.0
    pd1 = 0.0
    do j=j0,j1
      dd1 = dd1 + new%w(j)
      pd1 = pd1 + new%w(j)*new%v(:,j)
    end do
    !---------------------------------------------------------------------------
    ! Add the mass and momentum of the removed particle 
    !---------------------------------------------------------------------------
    dd1 = dd1 + new%w(i)
    pd1 = pd1 + new%w(i)*new%v(:,j)
    !---------------------------------------------------------------------------
    ! Compensate by reducing the weight and ajusting the velocity
    !---------------------------------------------------------------------------
    f = dd0/dd1
    pd1 = pd1 - self%pd(:,q(1),q(2),q(3))
    dv = pd1/dd0
    do j=j0,j1
      new%w(j)   = new%w(j)*f
      new%v(:,j) = new%v(:,j) + dv
    end do
    !---------------------------------------------------------------------------
    ! After removing the final one, renormalize
    !---------------------------------------------------------------------------
    self%pc(q(1),q(2),q(3)) = &
    self%pc(q(1),q(2),q(3)) - 1
    if (self%pc(q(1),q(2),q(3)) == per_cell) then
      call normalize
    end if
  end subroutine
  !=============================================================================
  ! Add a particle, representative of the existing one -- drawing randomly
  ! from the existing ones, adopting its mass and velocity.  When adding the
  ! last particle, renormalize to the given total mass and momentum
  !=============================================================================
  subroutine add (i)
    integer:: i
    integer:: j0, j1, j2, jr
    integer, pointer:: ic_p
    integer(2), pointer:: nc_p
    !---------------------------------------------------------------------------
    ic_p => new%ic(q(1),q(2),q(3))          ! pointer to particle index
    nc_p => new%pc(q(1),q(2),q(3))          ! pointer to particle count
    j1      = ic_p                          ! last particle index
    j0      = j1 - nc_p + 1                 ! first particle index
    jr = (j1-j0)*patch%random%ran1()        ! random integet, 0...(j1+j0-1)
    j2      = j0 + jr                       ! random particle choice
    j1      = j1 + 1                        ! new particle index
    new%w(j1)   = new%w(j2)                 ! adopt weight
    new%s(j1)   = new%s(j2)                 ! adopt size
    new%v(:,j1) = new%v(:,j2)               ! adopt velocity
    new%p(:,j1) = patch%random%uniform(3)   ! random position in cell
    ic_p    = j1 + 1                        ! increment cell particle index
    nc_p    = nc_p + 1                      ! increment cell particle count
    !---------------------------------------------------------------------------
    ! If this was the last particle to be added, renormalize
    !---------------------------------------------------------------------------
    if (nc_p == per_cell) then
      call normalize
    end if
  end subroutine add
  !=============================================================================
  ! Adjust mass and momentum to what it should be
  !=============================================================================
  subroutine normalize
    integer:: j, j0, j1
    real:: dd1, dd0, fd, pd1(3), dv(3), wtot
    integer, pointer:: ic_p
    integer(2), pointer:: nc_p
    !-------------------------------------------------------------------------
    ic_p => new%ic(q(1),q(2),q(3))      ! pointer to particle index
    nc_p => new%pc(q(1),q(2),q(3))      ! pointer to particle count
    j1      = ic_p                      ! last particle index
    j0      = j1 - nc_p + 1             ! first particle index
    !-------------------------------------------------------------------------
    ! Sum up existing weights and compute normalization factor
    !-------------------------------------------------------------------------
    dd1 = 0.0
    do j=j0,j1
      dd1 = dd1 + new%w(j)
    end do
    dd0 = self%dd(q(1),q(2),q(3))
    fd = dd0/dd1
    !-------------------------------------------------------------------------
    ! Sum up existing momentum and compute excess
    !-------------------------------------------------------------------------
    pd1 = 0.0
    dd1 = 0.0
    do j=j0,j1
      new%w(j) = new%w(j)*fd
      pd1 = pd1 + new%w(j)*new%v(:,j)
      dd1 = dd1 + new%w(j)
    end do
    pd1 = pd1 - self%pd(:,q(1),q(2),q(3))
    !-------------------------------------------------------------------------
    ! Subtract the excess, by making a change of velocity of each particle
    !-------------------------------------------------------------------------
    dv = pd1/dd1
    do j=j0,j1
      new%v(:,j) = new%v(:,j) - dv
    end do
    call check
  end subroutine
  !=============================================================================
  ! Check conservation of mass and momentum
  !=============================================================================
  subroutine check
    integer:: j, j0, j1
    real:: dd1, dd0, pd1(3), pd0(3)
    !---------------------------------------------------------------------------
    if (verbose > 0) then
      dd1 = 0.0
      pd1 = 0.0
      do j=j0,j1
        dd1 = dd1 + new%w(j)
        pd1 = pd1 + new%w(j)*new%v(:,j)
      end do
      dd0 = self%dd(q(1),q(2),q(3))
      pd0 = self%pd(:,q(1),q(2),q(3))
      if (abs(dd1/dd0-1.) > tolerance) then
        call io%abort ('not consevering momentum')
      end if
      if (sqrt(sum(pd-pd0)**2)/dd0 > tolerance) then
        call io%abort ('not consevering momentum')
      end if
    end if
  end subroutine
END SUBROUTINE renormalize

END MODULE particles_mod
