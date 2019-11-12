!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for patches, which adds pointers to memory and mesh, and
!> number of dimensions and variables.  The number of dimensions is hardwired
!> to 3, since we can just choose the actual dimensions to get 2-D and 1-D.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE iso_fortran_env, only: int8
  USE io_mod
  USE mpi_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mesg_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mod
  USE kinds_mod
  USE task_mod
  USE link_mod
  USE mesh_mod
  USE boundaries_mod
  USE index_mod
  USE shared_mod
  USE connect_mod
  implicit none
  private
  integer, parameter:: ndim = 3
  type, extends(task_t), public:: patch_t
    real(kind=KindScalarVar), dimension(:,:,:,:,:,:), pointer:: mem      => null()
    real(kind=KindScalarVar), dimension(:,:,:,:,:),   pointer:: vgas     => null()
    real(kind=KindScalarVar), dimension(:,:,:,:),     pointer:: density  => null()
    real(kind=KindScalarVar), dimension(:,:,:,:),     pointer:: pressure => null()
    !dir$ attributes align: 64 :: mem, vgas, density, pressure
    !---------------------------------------------------------------------------
    ! These allocatables are intended to deposite external forces and heating in,
    ! typically used by methods called from the extras_t data type procedures.
    !---------------------------------------------------------------------------
    real(kind=KindScalarVar), dimension(:,:,:,:),     allocatable:: force_per_unit_mass
    real(kind=KindScalarVar), dimension(:,:,:,:),     allocatable:: force_per_unit_volume
    real(kind=KindScalarVar), dimension(:,:,:),       allocatable:: heating_per_unit_mass
    real(kind=KindScalarVar), dimension(:,:,:),       allocatable:: heating_per_unit_volume
    integer(kind=int8)      , dimension(:,:,:),       allocatable:: filled
    !dir$ attributes align: 64 :: force_per_unit_mass, force_per_unit_volume
    !dir$ attributes align: 64 :: heating_per_unit_mass, heating_per_unit_volume
    !dir$ attributes align: 64 :: filled
    !--- anonymous pointer, which may be used to connect extra info to the solvers
    type(connect_t):: connect
    class(mesh_t), dimension(:), pointer:: mesh => null()
    real:: courant                      ! courant number
    real(8):: pdf_next=0.0              ! next PDF output time
    real(8):: gamma=1.4                 ! gas gamma
    real(8):: dt_fixed                  ! for fixed time step Courant condition
    real:: u_fixed, u_max               ! for fixed time step Courant condition
    integer:: ndim = ndim               ! number of dimensions
    integer:: nv=0                      ! number of variables
    integer:: nw=2                      ! number of work indices
    integer:: msplit=0                  ! how to split a patch when refining
    integer:: max_files                 ! prevent disk overrun when NaN happens
    integer:: check_refine_next=0
    integer:: time_derivs               ! for restart changes
    integer:: not_filled=0              ! guard zones not filled
    integer:: refine_ratio=2            ! may differ btw patches
    logical:: guard_zones               ! for restart changes
    logical:: mem_allocated=.false.
    logical:: do_pebbles=.false.
    logical:: staggered=.true.
    logical:: no_mans_land=.false.
    logical:: use_data_hub=.false.
    logical:: manual_refine=.false.
    logical:: get_umax_location=.false.
    logical, dimension(:), pointer:: unsigned  => null()
    logical, dimension(:), pointer:: pervolume => null()
    class(link_t), pointer:: link       ! this could be useful
    type(boundaries_t):: boundaries
    character(len=32):: eos='ideal'
    character(len=32):: opacity='none'
    integer:: file_offset(4)
    type(index_t):: idx
    real:: safety_factor=1.0
    !--- used in extras ---
    logical, allocatable:: mask(:,:,:)
    !--- used in refine ---
    integer(int8), allocatable:: irefine(:,:,:)
    type(lock_t), pointer:: mem_lock(:) => null()
    !---------------------------------------------------------------------------
    ! HACKs below:  FIXME
    !---------------------------------------------------------------------------
    real(8):: theta=0d0                 ! latitude
    real(8):: phi=0d0                   ! longitude at t=0.0
    real:: csound                       ! isothermal sound speed, for tests
    real:: grav                         ! gravitational consant (code units)
    integer:: n(3)=0, ng(3)=0           ! defaults, to detect unchanged values
    integer:: li(3), ui(3), gn(3)       ! convenience -- avoid patch%mesg%li etc
    integer:: ncell(3)                  ! number of cells per patch
    integer:: ipos(3)                   ! integer position, 3D
    real(8):: offset(3)                 ! this is used to compute indices, locations
    real(8):: llc_cart(3)=0d0           ! the lower left corner of the patch in *Cartesian* coords.
    real(8):: llc_nat(3)=0d0            ! the lower left corner of the patch in its *own* coords.
    real(8):: centre_cart(3)=0d0        ! the geometric centre of the patch in *Cartesian* coords.
    real(8):: centre_nat(3)=0d0         ! the geometric centre of the patch in its *own* coords.
    integer :: mesh_type = 1            ! default: Cartesian mesh; FIXME (there must be a better way)
    real(8):: grid_vel(3) = 0d0         ! velocity of the patch in its native coords. and
                                        ! units (e.g. rad/s for the phi coord.); for
                                        ! Cartesian coords., equal to `velocity`.
  contains
    procedure, nopass:: cast2patch
    procedure:: init
    procedure:: init_default
    procedure:: init_bdries
    procedure:: setup
    procedure:: dealloc
    procedure:: clone_mem_accounting
    procedure:: rotate
    procedure:: allocate_mesg           ! allocate a mesg buffer
    procedure:: pack                    ! pack patch into buffer
    procedure, nopass:: unpack_id       ! unpack patch id
    procedure:: unpack                  ! unpack buffer into patch
    procedure:: info                    ! print info to stdout
    procedure, private:: contains1      ! patch contains a point
    procedure, private:: contains2      ! patch contains a patch
    generic:: contains => contains1, contains2
    procedure:: overlaps                ! patch overlaps test
    procedure:: index                   ! compute the index for a single point
    procedure:: nearest                 ! the nearest index that maps to patch point
    procedure:: lowest                  ! the lowest index that maps to patch point
    procedure:: highest                 ! the highest index that maps to patch point
    procedure:: myposition
    procedure:: myposition_staggered
    procedure:: periodic_grid
    procedure:: extrapolate_guards
    procedure:: is_periodic
    procedure, private :: make_periodic4
    procedure, private :: make_periodic8
    generic, public :: make_periodic  => make_periodic4 , make_periodic8
    procedure, private :: fsum4
    procedure, private :: fsum8
    generic, public :: fsum    => fsum4, fsum8        ! interior sum
    procedure, private :: faver4
    procedure, private :: faver8
    generic, public :: faver    => faver4, faver8     ! interior average
    procedure, private :: fminval4
    procedure, private :: fminval8
    procedure, private :: fminvali
    generic, public :: fminval => fminval4, fminval8, fminvali
    procedure, private :: fmaxval4
    procedure, private :: fmaxval8
    procedure, private :: fmaxvali
    generic, public :: fmaxval => fmaxval4, fmaxval8, fmaxvali
    generic, public :: fminvalvar => fminval4, fminval8  ! these two overloaded functions are
    generic, public :: fmaxvalvar => fmaxval4, fmaxval8  ! used by `fminvali`/`fmaxvali` only.
    procedure:: courant_condition
    procedure:: check_density
    procedure:: check_nan
    procedure:: custom_refine
    procedure:: patch_to_header
    procedure:: header_to_patch
    procedure, private:: stats_patch
    procedure, private:: stats_scalar
    procedure, private:: stats_vector
    generic,  public:: stats => stats_patch, stats_scalar, stats_vector
    procedure:: MapVectors
    procedure:: beyond_patch_edge
    procedure:: myfloatindex_scalar
    procedure:: log_interpolate
    procedure:: time_interval
    procedure:: distance => distance_curvilinear
    procedure:: index_stagger
    procedure:: interpolate
    procedure:: interpolate_specific
    procedure:: index_space             ! full mapping to index space
    procedure:: index_space_of
    procedure:: index_space2
    procedure:: index_only
    procedure:: index_only2
    procedure:: count_edges
    procedure:: update_position
    procedure:: get_overlap
    procedure:: init_level
  end type
  !-----------------------------------------------------------------------------
  ! Sequential header, for MPI and I/O
  !-----------------------------------------------------------------------------
  integer, public, save:: n_header = (14*4 + 4*8 + 5*8*3 + 1*8 + 2*8*8 + 1*4*8 + 2*4*3)/4
  type, public:: header_t
    sequence
    integer:: version=1
    integer:: n_header
    integer:: n_data
    integer:: id
    integer:: status
    integer:: rank
    integer:: istep
    integer:: iout
    integer:: level
    integer:: nt
    integer:: nv                        ! number of variables
    integer:: it
    integer:: nq=0
    integer:: npad=0
    real(8):: time
    real(8):: out_next
    real(8):: print_next
    real(8):: dtime
    real(8):: position(3)
    real(8):: velocity(3)
    real(8):: size(3)
    real(8):: box(3)                    ! the periodic wrapping size
    real(8):: ds(3)
    real(8):: send_time
    real(8), dimension(8):: t
    real(8), dimension(8):: dt
    integer, dimension(8):: iit
    integer:: n(3)=16                   ! convenience
    integer:: ng(3)=2                   ! convenience
  end type
  !-----------------------------------------------------------------------------
  ! Template task, with procedure to place it it MPI-space
  !-----------------------------------------------------------------------------
  type box_t
    real:: size(3)=1d0
    real:: num(3)
    integer:: dim(3)
  end type
  integer, save:: print_every=0, iprint=0           ! print cadence
  logical, save:: extrapolate_gz=.false.            ! extrapolate guard zones?
  logical, save:: do_check_nan=.false.              ! check for NaN?
  logical, save:: zero_order_extrap=.true.          ! extrap to no-mans-land
  public task2patch
CONTAINS

!===============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!===============================================================================
SUBROUTINE patch_to_header (self, head)
  class(patch_t):: self
  type(header_t):: head
  integer:: nt
  !-----------------------------------------------------------------------------
  head%n_header   = n_header
  head%n_data     = product(self%gn)*self%nv
  head%id         = self%id
  head%status     = self%status
  head%rank       = self%rank
  head%istep      = self%istep
  head%iout       = self%iout
  head%level      = self%level
  head%time       = self%time
  head%out_next   = self%out_next
  head%print_next = self%print_next
  head%dtime      = self%dtime
  head%position   = self%position
  head%velocity   = self%velocity
  head%size       = self%size
  head%it         = self%it
  head%nt         = self%nt
  nt              = self%nt
  head%t(1:nt)    = self%t
  head%dt(1:nt)   = self%dt
  head%iit(1:nt)  = self%iit
  head%nq         = self%nq
  head%send_time  = wallclock()
  if (io%verbose > 1) &
    write (io_unit%log,'(f12.6,3x,a,i7,i4,l3,2i4,f8.4)') &
      wallclock(), 'patch_to_header: id, nq =', self%id, &
      self%nq, self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n
END SUBROUTINE patch_to_header

!===============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!===============================================================================
SUBROUTINE header_to_patch (self, head)
  class(patch_t):: self
  type(header_t):: head
  integer:: nt, new
  !-----------------------------------------------------------------------------
  ! Swap the roles of boundary and virtual on incoming patches.  Also clear the
  ! ready bit, which is an indicator of a patch being in the ready queue, which
  ! an incoming patch is not.  By doing this here, we avoid having to do it in
  ! a critical region, otherwise needed to protect an operation on the patch
  !-----------------------------------------------------------------------------
  if (iand(head%status, bits%boundary) /= 0) then
    head%status = ior(head%status, bits%virtual)
    head%status = iand(head%status, not(bits%boundary))
    if (io%verbose>2) &
      write (io_unit%mpi, *) head%id, 'swapping to virtual'
  else if (iand(head%status, bits%virtual) /= 0) then
    head%status = ior(head%status, bits%boundary)
    head%status = iand(head%status, not(bits%virtual))
    if (io%verbose>2) &
      write (io_unit%mpi, *) head%id, 'swapping to boundary'
  end if
  head%status = iand(head%status, not(bits%ready+bits%busy))
  if (head%id /= self%id) then
    call io%abort('wrong message ID unpacked')
  end if
  !-----------------------------------------------------------------------------
  ! Except for the first and last appearance, istep should increase be one
  !-----------------------------------------------------------------------------
  if ((self%istep > 0) .and. (iand(head%status, bits%remove) == 0)) then
    if (head%istep > self%istep+1) then
      write(io_unit%log,1) &
        'WARNING: early arrival in time reversal id, rank =', &
        self%id, mpi%rank, head%istep, self%istep, head%time, self%time
    1 format(a,4i6,2f12.6)
      flush (io_unit%log)
    else if (head%istep < self%istep+1) then
      write(io_unit%log,1) &
        'WARNING: late arrival in time reversal id, rank =', &
        self%id, mpi%rank, head%istep, self%istep, head%time, self%time
      flush (io_unit%log)
    end if
  end if
  self%id         = head%id
  self%status     = head%status
  self%rank       = head%rank
  self%istep      = head%istep
  self%iout       = head%iout
  self%level      = head%level
  self%time       = head%time
  self%out_next   = head%out_next
  self%print_next = head%print_next
  self%dtime      = head%dtime
  self%position   = head%position
  self%velocity   = head%velocity
  self%size       = head%size
  !$omp atomic write
  self%it         = head%it
  new             = mod(head%it,head%nt) + 1
  !$omp atomic write
  self%new        = new
  self%nt         = head%nt
  nt              = head%nt
  if (nt > 0) then
    self%t          = head%t(1:nt)
    self%dt         = head%dt(1:nt)
    self%iit        = head%iit(1:nt)
    !self%new        = head%iit(nt)
  else
    write (io_unit%log,*) omp%thread,'header_to_patch WARNING: nt non-positive', nt
    flush (io_unit%log)
  end if
  self%nq         = head%nq
  self%latency    = wallclock() - head%send_time
  !$omp atomic
  timer%latency%max  = max(timer%latency%max,self%latency)
  !$omp atomic
  timer%latency%aver = timer%latency%aver+self%latency
  !$omp atomic
  timer%latency%n    = timer%latency%n+1
  !.............................................................................
END SUBROUTINE header_to_patch

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2patch (task) RESULT(patch)
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  !.............................................................................
  select type (task)
  class is (patch_t)
  patch => task
  class default
  nullify(patch)
  call io%abort ('patch_t%cast: failed to cast a task to patch_t')
  end select
END FUNCTION cast2patch

!===============================================================================
!> Public function (deprecated)
!===============================================================================
FUNCTION task2patch (task) RESULT (patch)
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  !.............................................................................
  patch => cast2patch (task)
END FUNCTION task2patch

!===============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be told the
!> patch dimensions (which include the number of variable and the number of time
!> slices), the patch position, and the patch size
!===============================================================================
SUBROUTINE init (self)
  class(patch_t):: self
  character(len=120):: id = &
    '$Id$ tasks/patch_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%init')
  call trace%print_id (id)
  call self%setup
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize default values for the patch memory and patch dimensions. The 
!> order of setting values is that the calling code, typically experiment_mod,
!> may set preferred defaults that differ from the ones set here.  The input
!> namelist then allows these to be changed.  This is done (only) once, on first
!> call, after which the saved properties are reused on subsequent calls.
!>
!> However, in cases where not all patches should have identical properties it
!> is necessary be able to override the saved default.  In the case of n and ng
!> -- the most likely properties to vary -- this can be done by changing the patch
!> values before the call, since only properties with default (zero) values
!> are set to the saved values.  It is possible to change other parameters, by
!> instead of calling patch_t%init, first call patch_t%init_default, and then
!> call patch_t%setup (nt=..., nv=..., ...).
!>
!> Values that are expected to vary from patch to patch, such as size and
!> position, must be set before calling patch_t%init, or else must be set by
!> instead calling patch_t%init_defult + patch_t%setup.
!>
!> Any experiment that should start up with non-default values of parameters
!> given default values here (such as 'periodic'), should use the method:
!> call init_default + set values + call setup, since otherwise it would have
!> to rely on the values being set correctly in the input namelist.
!===============================================================================
SUBROUTINE init_default (self)
  class(patch_t):: self
  !.............................................................................
  integer, save:: n(3)=16, ng(3)=2, nv=5, nt=5, nw, max_files=1000
  logical, save:: no_mans_land=.false., use_data_hub=.false., limit_dtime=.false.
  real, save:: u_fixed=-1.0, grace=0.0
  real(8), save:: end_time=1.0, out_time=0.1, print_time=0.1, dt_fixed=-1.0
  real(8):: position(3), size(3)
  integer:: iostat, i, n_alloc=0, n_ext=0
  namelist /patch_params/ n, nv, ng, nt, u_fixed, dt_fixed, &
    grace, extrapolate_gz, no_mans_land, use_data_hub, do_check_nan, &
    zero_order_extrap, limit_dtime
  namelist /out_params/ out_time, end_time, print_time, print_every, &
    max_files
  logical, save:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%init_default')
  call self%task_t%init                         ! sets unique id
  !-----------------------------------------------------------------------------
  ! Read and set patch parameters; the preferred defalts may have been set by
  ! the calling code -- typically in experiment_mod.f90
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    ! -- if values have been imposed from above, FI show them as default
    if (any(self%n  /= 0)) n  = self%n          ! patch dimensions
    if (any(self%ng /= 0)) ng = self%ng         ! patch guard zones
    if (    self%nv /= 0 ) nv = self%nv         ! number of variables
    if (    self%nt /= 0 ) nt = self%nt         ! number of time slices
    first_time = .false.
    rewind (io%input)
    read (io%input, patch_params)
    if (dt_fixed > 0.0) grace=max(grace,1e-6)
    if (io%master) write (*, patch_params)
    rewind (io%input)
    read (io%input, out_params)
    if (io%master) write (*, out_params)
    if (mpi%rank > 0) &
      self%track = .true.
    iprint = print_every                        ! print 1st task
  end if
  !-----------------------------------------------------------------------------
  ! To allow patches to have n, ng, nv, and nt values that differ, values that
  ! differ from the data type zero defaults take precedence over saved values
  !-----------------------------------------------------------------------------
  if (any(self%n  == 0)) self%n  = n            ! patch dimensions
  if (any(self%ng == 0)) self%ng = ng           ! patch guard zones
  if (    self%nv == 0 ) self%nv = nv           ! number of variables
  if (    self%nt == 0 ) self%nt = nt           ! number of time slices
  !-----------------------------------------------------------------------------
  ! Output related parameters
  !-----------------------------------------------------------------------------
  io%end_time      = end_time
  io%out_time      = out_time
  io%print_time    = print_time
  if (io%nv==0) &
    io%nv = self%nv
  !$omp end critical (input_cr)
  self%restart     = io%restart                 ! restart snapshot
  self%grace       = grace                      ! grace fraction of dtime
  self%dt_fixed    = dt_fixed
  self%u_fixed     = u_fixed
  self%max_files   = max_files                  ! prevent run-away NaN
  self%use_data_hub = use_data_hub              ! guard zone handling
  self%no_mans_land = no_mans_land              ! mesh centering
  self%print_next  = print_time
  self%limit_dtime = limit_dtime
  call trace%end()
END SUBROUTINE init_default

SUBROUTINE init_bdries (self)
  class(patch_t):: self
  real(8):: position(3), limit(3)
  !-----------------------------------------------------------------------------
  ! If at edge, set boundary bits
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%init_bdries')
  position = self%position - self%origin - 0.5d0 * self%box
  limit = 0.5d0 * (self%box - self%size - self%ds)
  if (.not.self%periodic(1)) then
    if (position(1)  <= -limit(1)) then
      call self%boundaries%set (bits%xl)
      self%mesh(1)%lower_boundary = .true.
    end if
    if (position(1) >= +limit(1)) then
      call self%boundaries%set (bits%xu)
      self%mesh(1)%upper_boundary = .true.
    end if
  end if
  if (.not.self%periodic(2)) then
    if (position(2)  <= -limit(2)) then
      call self%boundaries%set (bits%yl)
      self%mesh(2)%lower_boundary = .true.
    end if
    if (position(2) >= +limit(2)) then
      call self%boundaries%set (bits%yu)
      self%mesh(2)%upper_boundary = .true.
    end if
  end if
  if (.not.self%periodic(3)) then
    if (position(3)  <= -limit(3)) then
      call self%boundaries%set (bits%zl)
      self%mesh(3)%lower_boundary = .true.
    end if
    if (position(3) >= +limit(3)) then
      call self%boundaries%set (bits%zu)
      self%mesh(3)%upper_boundary = .true.
    end if
  end if
  call trace%end()
END SUBROUTINE init_bdries

!===============================================================================
!> Setup a patch, possibly with non-standard parameters. First setup with
!> defaults, then possibly override, then do the actual setup
!===============================================================================
SUBROUTINE setup (self, position, size, n, ng, nt, nv, nw)
  class(patch_t)    :: self
  real(8), optional :: position(3), size(3)
  integer, optional :: n(3), ng(3)
  integer, optional :: nt, nv, nw
  !.............................................................................
  type(header_t):: header
  character(len=4) kind
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%setup')
  call self%init_default                        ! read default params
  if (present(position)) self%position = position
  if (present(size    )) self%size     = size
  if (present(n       )) self%n        = n
  if (present(ng      )) self%ng       = ng
  if (present(nt      )) self%nt       = nt
  if (present(nv      )) self%nv       = nv
  if (present(nw      )) self%nw       = nw
!print '(a,3(2x,3i4),2(2x,1p,3e11.2))', &
!'patch_t%setup: n,ng,nt,nv,nw,size,pos=', &
!self%n, self%ng, self%nv, self%nt, self%nw, self%size, self%position
  !-----------------------------------------------------------------------------
  ! Pass the desired patch dimension via the shared_mod to refine_mod
  !-----------------------------------------------------------------------------
  shared%patch_n = self%n
  !-----------------------------------------------------------------------------
  self%llc_cart = self%position - 0.5 * self%size
  self%llc_nat = self%llc_cart
  !-----------------------------------------------------------------------------
  ! Detect root grid
  !-----------------------------------------------------------------------------
  if (all(self%size==self%box)) then
    !$omp critical (log_cr)
    write (io_unit%log,*) mpi%rank, 'task', self%id, ' is root grid'
    call self%set(bits%root_grid)
    !$omp end critical (log_cr)
  end if
  !-----------------------------------------------------------------------------
  ! Special handling of collapsed dimensions
  !-----------------------------------------------------------------------------
  where (self%n==1)                             ! for collapsed dimensions
    self%ng = 0                                 ! no guard zones
    self%box = 1                                ! make sure it is non-zero
  end where
  self%llc_cart = self%position - 0.5_8*self%size
  if (self%mesh_type == mesh_types%Cartesian) then
    self%llc_nat = self%llc_cart
    self%centre_nat = self%position
  end if
  !-----------------------------------------------------------------------------
  ! Initialize and allocate the meshes, afterwards imposing the external values
  ! box size and position
  !-----------------------------------------------------------------------------
  call omp%set_stacksize (self%nv)
  call MeshFactory(self%mesh_type, self%mesh, n=self%n, ng=self%ng, s=self%size, &
                   p=self%position, llc_native=self%llc_nat, nml=self%no_mans_land)
  ! intialise staggering indices; it is the responsibility of the *solver* to
  ! assign non-zero values to this variable (if any).
  do i=1,3
    allocate(self%mesh(i)%h(self%nv))
    self%mesh(i)%h = 0.0
    if (self%kind == 'zeus_mhd_patch' .and. &
        self%staggered .and. self%n(i) /= 1) self%mesh(i)%gn = self%mesh(i)%gn + 1
    if (self%kind == 'zeus_tw_mhd_patch' .and. &
        self%staggered .and. self%n(i) /= 1) self%mesh(i)%gn = self%mesh(i)%gn + 1
  end do
  self%ds = self%mesh%d
  self%ncell = self%mesh%nc
  self%mesh%id = self%id
  self%mesh%b = self%box
  self%mesh%llc_cart = self%llc_cart
  self%centre_nat = [self%mesh(1)%centre_nat,self%mesh(2)%centre_nat,self%mesh(3)%centre_nat]
  self%centre_cart = CurrentToCartesian(self%mesh_type, &
                       self%centre_nat(1), self%centre_nat(2), self%centre_nat(3))
  !-----------------------------------------------------------------------------
  ! Hack for backwards compatibility
  !-----------------------------------------------------------------------------
  self%gn = self%mesh%gn
  self%li = self%mesh%li
  self%ui = self%mesh%ui
  self%gsize = (self%size*self%gn/self%n)
  self%offset = (self%mesh%lb+self%mesh%ub)/2.0 ! mid point index offset
  if (self%id == io%id_debug) then
    print 1,'patch_t%init: position    =', self%position
    print 1,'patch_t%init: m%p         =', self%mesh%p
    print 1,'patch_t%init: offset      =', self%offset
    print 1,'patch_t%init: m%o         =', self%mesh%o
    print 1,'patch_t%init: centre_nat  =', self%centre_nat
    print 1,'patch_t%init: centre_cart =', self%centre_cart
    print 1,'patch_t%init: llc_nat     =', self%llc_nat
    print 1,'patch_t%init: m%llc_nat   =', self%mesh%llc_nat
    print 1,'patch_t%init: llc_cart    =', self%llc_cart
    print 1,'patch_t%init: m%llc_cart  =', self%mesh%llc_cart
    1 format(a,1p,3e15.6)
  end if
  !-----------------------------------------------------------------------------
  ! Measure levels in terms of the refinement factor -- make sure this doesn't 
  ! come in conflict with levels set by e.g. rubiks_mod
  !-----------------------------------------------------------------------------
  call self%init_level
  if (io%verbose>2) then
    !$omp critical (log_cr)
    write (io_unit%log,'(1x,a,i7,2(1p,3e12.4,2x))') 'patch_mod::init: id, size, pos =', &
      self%id, self%size, self%position
    if (io%verbose>3) &
      write (io_unit%log,*) 'id, LEVEL: ', self%id, self%level, maxval(self%box), minval(self%ds)
    !$omp end critical (log_cr)
  end if
  !-----------------------------------------------------------------------------
  ! Allocate memory, if not already allocated
  !-----------------------------------------------------------------------------
  if (.not.self%mem_allocated) then
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
!print '(a,3(2x,3i4))', &
!'patch_t: gn, nv, nt, nw =', self%gn, self%mesh%gn, self%nv, self%nt, self%nw
    allocate (self%mem(self%mesh(1)%gn, &
                       self%mesh(2)%gn, &
                       self%mesh(3)%gn, &
                       self%nv,self%nt,self%nw))
    allocate (self%mem_lock(self%nt))
    do i=1,self%nt
      write (kind,'("mem",i1)') i
      call self%mem_lock(i)%init (kind)
    end do
    self%mem_allocated = .true.               ! mark
    call io%bits_mem (storage_size(self%mem), product(shape(self%mem)), 'mem')
    allocate (self%unsigned (self%nv))
    allocate (self%pervolume(self%nv))
    self%unsigned (:) = .false.
    self%pervolume(:) = .false.
  end if
  self%dt = 0.0
  self%t = -1.0                                 ! mark as "no guard zones"
  self%time = -1.0                              ! mark as "no guard zones"
  self%iit = 1                                  ! initial memory slots
  self%it = 1                                   ! short-hand (=1 initially)
  self%new = min(2,self%nt)                     ! slot for updates
  self%iit(self%nt) = self%new                  ! first output slot
  self%new = self%iit(self%nt)                  ! short-hand
  self%t(self%it) = 0.0                         ! the only valid slot (no. 1)
  self%wallclock = wallclock()                  ! starting time
  n_header = storage_size(header)/32             ! 32 bits per word
  if (n_header*4 /= storage_size(header)/8) then
    write(io%output,*) n_header*4, storage_size(header)/8
    call mpi%abort ('n_header is incorrect')
  end if
  if (all(self%size == self%box)) then
    call self%set (bits%root_grid)
    if (io%verbose > 1) write(io%output,*) 'set bits%root_grid id =', self%id
  end if
  call trace%end()
END SUBROUTINE setup

!===============================================================================
!> Deallocate a patch and its allocated parts
!===============================================================================
SUBROUTINE dealloc (self)
  class (patch_t):: self
  integer                  :: i
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%dealloc')
  if (associated(self%mesh)) &
    call MeshRecycler (self%mesh)
  call io%bits_mem (-storage_size(self%mem), product(shape(self%mem)),'-mem')
  if (associated(self%mem)) deallocate (self%mem)
  nullify (self%mem)
  self%mem_allocated = .false.
  if (allocated(self%force_per_unit_mass)) then
    call io%bits_mem (-storage_size(self%force_per_unit_mass), &
                      product(shape(self%force_per_unit_mass)),'-fpm')
    deallocate (self%force_per_unit_mass)
  end if
  if (allocated(self%force_per_unit_volume)) then
    call io%bits_mem (-storage_size(self%force_per_unit_volume), &
                      product(shape(self%force_per_unit_volume)),'-fpv')
    deallocate (self%force_per_unit_volume)
  end if
  if (allocated(self%heating_per_unit_mass)) then
    call io%bits_mem (-storage_size(self%heating_per_unit_mass), &
                      product(shape(self%heating_per_unit_mass)),'-qpm')
    deallocate (self%heating_per_unit_mass)
  end if
  if (allocated(self%heating_per_unit_volume)) then
    call io%bits_mem (-storage_size(self%heating_per_unit_volume), &
                      product(shape(self%heating_per_unit_volume)),'-qpv')
    deallocate (self%heating_per_unit_volume)
  end if
  if (associated(self%vgas)) then
    call io%bits_mem (-storage_size(self%vgas), &
                      product(shape(self%vgas)),'-vgas')
    deallocate (self%vgas)
  end if
  if (associated(self%density)) then
    call io%bits_mem (-storage_size(self%density), &
                      product(shape(self%density)),'-density')
    deallocate (self%density)
  end if
  if (associated(self%pressure)) then
    call io%bits_mem (-storage_size(self%pressure), &
                      product(shape(self%pressure)),'-pressure')
    deallocate (self%pressure)
  end if
  if (allocated(self%mask)) then
    call io%bits_mem (-storage_size(self%mask), &
                      product(shape(self%mask)), '-mask')
    deallocate (self%mask)
  end if
  if (allocated(self%irefine)) then
    call io%bits_mem (-storage_size(self%irefine), &
                      product(shape(self%irefine)), '-irefine')
    deallocate (self%irefine)
  end if
  call self%task_t%dealloc
  call trace%end()
END SUBROUTINE dealloc

SUBROUTINE clone_mem_accounting (self)
  class(patch_t):: self
  if (allocated(self%mask)) &
    call io%bits_mem (storage_size(self%mask), &
                     product(shape(self%mask)), 'mask')
  if (allocated(self%irefine)) &
    call io%bits_mem (storage_size(self%irefine), &
                     product(shape(self%irefine)), 'irefine')
  if (allocated(self%force_per_unit_mass)) &
    call io%bits_mem (storage_size(self%force_per_unit_mass), &
                     product(shape(self%force_per_unit_mass)),'fpm')
  if (allocated(self%force_per_unit_volume)) &
    call io%bits_mem (storage_size(self%force_per_unit_volume), &
                     product(shape(self%force_per_unit_volume)),'fpv')
  if (allocated(self%heating_per_unit_mass)) &
    call io%bits_mem (storage_size(self%heating_per_unit_mass), &
                     product(shape(self%heating_per_unit_mass)),'qpm')
  if (allocated(self%heating_per_unit_volume)) &
    call io%bits_mem (storage_size(self%heating_per_unit_volume), &
                     product(shape(self%heating_per_unit_volume)),'qpv')
END SUBROUTINE clone_mem_accounting

!===============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!===============================================================================
SUBROUTINE rotate (self)
  class(patch_t):: self
  integer:: iv, nv, new
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Periodic wrapping, time extrapolation
  !-----------------------------------------------------------------------------
  if (self%rotated) return
  call trace%begin ('patch_t%rotate', itimer=itimer)
  if (self%is_periodic()) then
    call self%periodic_grid                     ! periodic root grid
  else if (extrapolate_gz) then
    call self%extrapolate_guards                ! time-extrapolate guard cells
  end if
  !-----------------------------------------------------------------------------
  ! One line summary
  !-----------------------------------------------------------------------------
  if (io%verbose>1) then
    write (io_unit%log,'(a,i6,2i3,2e13.6,2x,a,1p,6e10.2)') 'update: id,j,k,l,n,dt,t', &
      self%id, self%it, self%new, self%dtime, self%time, 'minmax:', &
      minval(self%mem(:,:,:,1,self%it ,self%nw)), &
      maxval(self%mem(:,:,:,1,self%it ,self%nw)), &
      minval(self%mem(:,:,:,1,self%new,1)), &
      maxval(self%mem(:,:,:,1,self%new,1))
  end if
  !-----------------------------------------------------------------------------
  ! Update patch position for non-zero patch velocities;
  ! This happens *before* `rotate` because `dtime` stores the time step just taken.
  !-----------------------------------------------------------------------------
  call self%update_position
  !-----------------------------------------------------------------------------
  ! Task time slot rotation
  !-----------------------------------------------------------------------------
  call self%task_t%rotate
  call trace%end (itimer)
END SUBROUTINE rotate

!===============================================================================
!> This function gives the physical position of an index triple (patch point)
!> in absolute Cartesian coords.
!===============================================================================
FUNCTION myposition (self, ii, iv, pp) RESULT(out)
  class(patch_t)      :: self
  integer, intent(in) :: ii(3)
  integer, optional   :: iv
  real, optional      :: pp(3)
  !.............................................................................
  real(8)             :: out(3), cellpos(3), h(3)
  !.............................................................................
  if (present(pp)) then
    cellpos = self%position + (real(ii)-self%offset + pp)*self%ds
  else
    cellpos = self%position + (real(ii)-self%offset)*self%ds
  end if
  if (present(iv)) then
    h = self%index_stagger(iv)
    cellpos = cellpos + h*self%mesh%d
  end if
  if (self%id==io%id_debug .and. io%verbose>1) &
    write(io%output,'(a,3i4,4(3x,3f10.2))') &
      'mypos:',ii,self%position,self%offset,self%ds,cellpos
  ! convert to Cartesian coords.
  if (self%mesh_type == mesh_types%cylindrical) then
    cellpos = CylindricalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  else if (self%mesh_type == mesh_types%spherical) then
    cellpos = SphericalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  end if
  out = cellpos ! absolute position
END FUNCTION myposition

!===============================================================================
!> Scalar version of the previous function, but **including staggering factor**.
!> For the conversion to Cartesian coords., all three indices are required!
!===============================================================================
FUNCTION myposition_staggered (self, ii, idir, h) RESULT(out)
  class(patch_t):: self
  integer, intent(in):: ii(3), idir
  real(8):: out, cellpos(3)
  real(8), optional :: h
  real(4) :: hlocal
  !.............................................................................
  if (present(h)) then
    hlocal = h
  else
    hlocal = 0.0
  end if

  cellpos = (real(ii-self%li,kind=8) + merge(0.5 + hlocal,0.0,self%n(idir)>1)) &
          * self%ds + self%llc_nat ! in patch coords.
  ! convert to Cartesian coords.
  if (self%mesh_type == mesh_types%cylindrical) then
    cellpos = CylindricalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  else if (self%mesh_type == mesh_types%spherical) then
    cellpos = SphericalToCartesian(cellpos(1), cellpos(2), cellpos(3))
  end if
  out = cellpos(idir)  ! absolute position

END FUNCTION myposition_staggered

!===============================================================================
FUNCTION myfloatindex_scalar (self, cpos, idir, h) RESULT(out)
  class(patch_t) :: self
  real(8), intent(in) :: cpos(3)
  integer, intent(in) :: idir
  real(8), intent(in) :: h
  real(8) :: pp(3), ppCart(3), out, intermezzo
  !.............................................................................
  ! convert to native patch coords.
  if (self%mesh_type == mesh_types%Cartesian) then
    ! do nothing
    pp = cpos
  else if (self%mesh_type == mesh_types%spherical) then
    ppCart = cpos - self%position ! Cart. coords to centre of patch's coord. system
    pp = CartesianToSpherical(ppCart(1), ppCart(2), ppCart(3))
  else if (self%mesh_type == mesh_types%cylindrical) then
    ppCart = cpos - self%position ! Cart. coords to centre of patch's coord. system
    pp = CartesianToCylindrical(ppCart(1), ppCart(2), ppCart(3))
  end if

  intermezzo = (pp(idir) - self%llc_nat(idir)) / self%ds(idir) + self%li(idir) &
             - merge(0.5d0 + h,0.0d0,self%n(idir)>1)
  out = intermezzo

END FUNCTION myfloatindex_scalar

!===============================================================================
!> Periodic mapping of a grid -- typically the root grid. This is typically 
!> only used in tests, so does not need to be efficient.
!===============================================================================
SUBROUTINE periodic_grid (self, only)
  class(patch_t):: self
  integer, optional:: only
  !.............................................................................
  integer:: ix, iy, iz, i(3), j(3), lb(3), ub(3), li(3), ui(3), n(3)
  logical:: periodic(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%periodic_grid')
  periodic = self%periodic .and. self%size+self%ds > self%box
  if (.not.any(periodic)) call mpi%abort ('patch_t%periodic_grid')
  lb = self%mesh%lb
  ub = self%mesh%ub
  li = self%mesh%li
  ui = self%mesh%ui
  where (self%mesh%lower_boundary)
    li = self%mesh%lb
  end where
  where (self%mesh%upper_boundary)
    ui = self%mesh%ub
  end where
  n = self%ncell
  do iz=lb(3),ub(3)
  do iy=lb(2),ub(2)
  do ix=lb(1),ub(1)
    i = [ix,iy,iz]
    if (all(i >= li .and. i <= ui)) cycle
    where (n > 0)
      j = modulo(i-li,n)+li
    else where
      j = 1
    end where
    !---------------------------------------------------------------------------
    ! Prevent wrapping in non-periodic directions
    !---------------------------------------------------------------------------
    where (.not.periodic)
      j = i
    end where
    !---------------------------------------------------------------------------
    ! Make sure horizontal periodic wrapping continues over the whole vertical
    ! range when physical boundary conditions are present.  Note the use of
    ! mesh%ui and mesh%li below; this is because the physical boundaries are
    ! inside the normal patch boundaries.
    !---------------------------------------------------------------------------
    where (self%mesh%lower_boundary .and. i < self%mesh%li+4)
      j = i
    end where
    where (self%mesh%upper_boundary .and. i > self%mesh%ui-4)
      j = i
    end where
    if (present(only)) then
      self%mem(ix,iy,iz,only,self%it ,1) = self%mem(j(1),j(2),j(3),only,self%it ,1)
      self%mem(ix,iy,iz,only,self%new,1) = self%mem(j(1),j(2),j(3),only,self%new,1)
    else
      self%mem(ix,iy,iz,:,self%it ,1) = self%mem(j(1),j(2),j(3),:,self%it ,1)
      self%mem(ix,iy,iz,:,self%new,1) = self%mem(j(1),j(2),j(3),:,self%new,1)
    end if
  end do
  end do
  end do
  if (io%verbose > 1) write(io%output,*) &
    'periodic_grid: id, it, new =', self%id, self%it, self%new
  call trace%end()
END SUBROUTINE periodic_grid

!===============================================================================
FUNCTION is_periodic(self) RESULT(out)
  class(patch_t):: self
  logical:: out
  !.............................................................................
  out = any(self%periodic .and. self%size+self%ds > self%box)
END FUNCTION is_periodic

!=======================================================================
!> Optimized periodic boundary enforcer
!=======================================================================
SUBROUTINE make_periodic4 (patch, f)
  class(patch_t):: patch
  class(mesh_t), pointer:: m(:)
  real(4):: f(:,:,:)
  integer:: ix, iy, iz, n(3)
  logical:: periodic(3)
  !-----------------------------------------------------------------------------
  periodic = patch%periodic .and. patch%size+patch%ds > patch%box
  m => patch%mesh
  n = patch%ncell
  !-----------------------------------------------------------------------------
  ! x periodicity; make periodic inside the yz patch limits
  !-----------------------------------------------------------------------------
  if (n(1) > 1 .and. periodic(1)) then
    do iz=m(3)%li,m(3)%ui
      do ix=m(1)%lb,m(1)%lo
        do iy=m(2)%li,m(2)%ui
          f(ix,iy,iz) = f(ix+n(1),iy,iz)
        end do
      end do
      do ix=m(1)%uo,m(1)%ub
        do iy=m(2)%li,m(2)%ui
          f(ix,iy,iz) = f(ix-n(1),iy,iz)
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! y periodicity, full range in x now
  !-----------------------------------------------------------------------------
  if (n(2) > 1 .and. periodic(2)) then
    do iz=m(3)%li,m(3)%ui
      do iy=m(2)%lb,m(2)%lo
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy+n(2),iz)
        end do
      end do
      do iy=m(2)%uo,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy-n(2),iz)
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! z periodicity, full range in xy
  !-----------------------------------------------------------------------------
  if (n(3) > 1 .and. periodic(3)) then
    do iy=m(2)%lb,m(2)%ub
      do iz=m(3)%lb,m(3)%lo
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy,iz+n(3))
        end do
      end do
      do iz=m(3)%uo,m(3)%ub
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy,iz-n(3))
        end do
      end do
    end do
  end if
END SUBROUTINE make_periodic4

!=======================================================================
!> Optimized periodic boundary enforcer
!=======================================================================
SUBROUTINE make_periodic8 (patch, f)
  class(patch_t):: patch
  class(mesh_t), pointer:: m(:)
  real(8):: f(:,:,:)
  integer:: ix, iy, iz, n(3)
  logical:: periodic(3)
  !-----------------------------------------------------------------------------
  periodic = patch%periodic .and. patch%size+patch%ds > patch%box
  m => patch%mesh
  n = patch%ncell
  !-----------------------------------------------------------------------------
  ! x periodicity; make periodic inside the yz patch limits
  !-----------------------------------------------------------------------------
  if (n(1) > 1 .and. periodic(1)) then
    do iz=m(3)%li,m(3)%ui
      do ix=m(1)%lb,m(1)%lo
        do iy=m(2)%li,m(2)%ui
          f(ix,iy,iz) = f(ix+n(1),iy,iz)
        end do
      end do
      do ix=m(1)%uo,m(1)%ub
        do iy=m(2)%li,m(2)%ui
          f(ix,iy,iz) = f(ix-n(1),iy,iz)
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! y periodicity, full range in x now
  !-----------------------------------------------------------------------------
  if (n(2) > 1 .and. periodic(2)) then
    do iz=m(3)%li,m(3)%ui
      do iy=m(2)%lb,m(2)%lo
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy+n(2),iz)
        end do
      end do
      do iy=m(2)%uo,m(2)%ub
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy-n(2),iz)
        end do
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! z periodicity, full range in xy
  !-----------------------------------------------------------------------------
  if (n(3) > 1 .and. periodic(3)) then
    do iy=m(2)%lb,m(2)%ub
      do iz=m(3)%lb,m(3)%lo
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy,iz+n(3))
        end do
      end do
      do iz=m(3)%uo,m(3)%ub
        do ix=m(1)%lb,m(1)%ub
          f(ix,iy,iz) = f(ix,iy,iz-n(3))
        end do
      end do
    end do
  end if
END SUBROUTINE make_periodic8

!===============================================================================
!> Time extrapolate guard cell values, to allow using also guard zone values to
!> get boundary conditions for higher level patches.  This is done before time
!> slot rotation, so the slot that has just been updated is iit(nt).
!===============================================================================
SUBROUTINE extrapolate_guards (self)
  class(patch_t)        :: self
  !.............................................................................
  integer               :: it0, it1, it2, nt, iit(self%nt)
  integer, dimension(3) :: lb, lo, li, ub, uo, ui, i
  integer               :: iv
  integer, save         :: itimer=0
  real(8)               :: t(self%nt)
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%extrapolate_guards', 3, itimer=itimer)
  nt = self%nt
  call self%timeslots (iit, t)
  it0 = iit(nt)                                         ! current time slot
  it1 = iit(nt-1)                                       ! previous time slot
  it2 = iit(nt-2)                                       ! the one before that
  !-----------------------------------------------------------------------------
  ! Before someone asks: Associate on the non-changing 6 doesn't please gfort
  !-----------------------------------------------------------------------------
  lb = self%mesh%lb
  li = self%mesh%li
  lo = self%mesh%lo
  ub = self%mesh%ub
  ui = self%mesh%ui
  uo = self%mesh%uo
  where (self%mesh%lower_boundary)
    li = self%mesh%lb
  end where
  where (self%mesh%upper_boundary)
    ui = self%mesh%ub
  end where
  if (io%id_debug>0.and.(self%id==io%id_debug.or.abs(self%id-io%id_debug)==27)) then
    write(io%output,'(i5,2x,a,6i7)') self%id, &
      'MK extrapolating guard zone values', it2, it1, it0, self%it, self%new
  end if
  do iv=1,self%nv
    if (self%unsigned(iv)) then
      !-------------------------------------------------------------------------
      ! XY-layers (Z varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(3) > 1) then
            self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it1,1))*2 - &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it2,1)) )
            self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it1,1))*2 - &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it2,1)) )
      end if
      !-------------------------------------------------------------------------
      ! ZX-layers (Y varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(2) > 1) then
            self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it1,1))*2 - &
        log(self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it2,1)) )
            self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it1,1))*2 - &
        log(self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it2,1)) )
      end if
      !-------------------------------------------------------------------------
      ! YZ-layers (X varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(1) > 1) then
            self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it1,1))*2 - &
        log(self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it2,1)) )
            self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it1,1))*2 - &
        log(self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it2,1)) )
      end if
    else
      !-------------------------------------------------------------------------
      ! XY-layers (Z varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(3) > 1) then
        self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it1,1)*2 - &
        self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it2,1)
        self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it1,1)*2 - &
        self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it2,1)
      end if
      !-------------------------------------------------------------------------
      ! ZX-layers (Y varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(2) > 1) then
        self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it1,1)*2 - &
        self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it2,1)
        self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it1,1)*2 - &
        self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it2,1)
      end if
      !-------------------------------------------------------------------------
      ! YZ-layers (X varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(1) > 1) then
        self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it1,1)*2 - &
        self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it2,1)
        self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it1,1)*2 - &
        self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it2,1)
      end if
    end if
  end do
  call trace%end (itimer)
END SUBROUTINE extrapolate_guards

!===============================================================================
!> Pack patch info into a buffer, for MPI communication
!===============================================================================
SUBROUTINE allocate_mesg (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  allocate (self%mesg)
  self%mesg%nbuf = product(self%gn)      ! 4-bytes per element; `anonymous_copy` also uses 4-bytes.
  if (kind(self%mem)==8) &
    self%mesg%nbuf = self%mesg%nbuf * 2  ! 8-bytes per element, but `anonymous_copy` uses 4-bytes!
  self%mesg%nbuf = n_header + self%nv*self%mesg%nbuf
  allocate (self%mesg%buffer(self%mesg%nbuf))
  call io%bits_mem (storage_size(self%mesg%buffer), product(shape(self%mesg%buffer)))
END SUBROUTINE allocate_mesg

!===============================================================================
!> Pack patch info into a buffer, for MPI communication
!===============================================================================
SUBROUTINE pack (self, mesg)
  class(patch_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  type(header_t):: header
  integer:: n_buf, ibuf, iv, it
  real:: pt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer and allocate -- it will be freed by the writer
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%pack', itimer=itimer)
  allocate (mesg)
  n_buf = product(self%gn)      ! 4-bytes per element; `anonymous_copy` also uses 4-bytes.
  if (kind(self%mem)==8) &
    n_buf = n_buf * 2           ! 8-bytes per element, but `anonymous_copy` uses 4-bytes!
  mesg%nbuf = n_header + self%nv*n_buf
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%boundary)) then
    mesg%nbuf = mesg%nbuf+ self%nv*n_buf
  end if
  allocate (mesg%buffer(mesg%nbuf))
  if (mpi_mesg%nbuf==0) then
    !$omp critical (nbuf_cr)
    if (mpi_mesg%nbuf==0) then
      mpi_mesg%nbuf = mesg%nbuf
      write (io_unit%log,*) &
        'patch_t%pack: setting mpi_mesg%nbuf =', mpi_mesg%nbuf
    end if
    !$omp end critical (nbuf_cr)
  end if
  call io%bits_mem (storage_size(mesg%buffer), product(shape(mesg%buffer)))
  allocate (mesg%reqs(self%n_nbors))
  mesg%id = self%id                                             ! MPI message tag
  mesg%nreq = 0                                                 ! reset counter
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info to sequenced header, and copy that to the buffer
  !-----------------------------------------------------------------------------
  if (self%debug(2)) then
    write (io_unit%log,1) wallclock(), &
    'DBG patch_t%pack, id, nq, time, status:', &
    self%id, self%nq, self%time, &
    self%is_set (bits%internal), &
    self%is_set (bits%boundary), &
    self%is_set (bits%virtual), &
    self%is_set (bits%external), &
    self%is_set (bits%swap_request)
    1 format(f12.6,2x,a,2i9,f12.6,2x,5l1)
    flush (io_unit%log)
  end if
  call patch_to_header (self, header)
  ibuf = 1
  call anonymous_copy (n_header, header, mesg%buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer.  If the bdry+swap bit are set,
  ! copy as a minimum also the previous time slot
  !-----------------------------------------------------------------------------
  do iv=1,self%nv
    call anonymous_copy (n_buf, self%mem(:,:,:,iv,self%it,1), mesg%buffer(ibuf))
    ibuf = ibuf + n_buf
  end do
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%boundary)) then
    it = mod(self%it-1+self%nt-1,self%nt)+1
    if (io%verbose>1) then
      write (io_unit%log,*) self%id, 'extra time slot packed', it
      flush (io_unit%log)
    end if
    do iv=1,self%nv
      call anonymous_copy (n_buf, self%mem(:,:,:,iv,it,1), mesg%buffer(ibuf))
      ibuf = ibuf + n_buf
    end do
  end if
  if (io%verbose > 2) &
    write (io_unit%log,*) self%id, 'after pack', minval(self%mem(:,:,:,1,self%it,1))
  call trace%end (itimer)
END SUBROUTINE pack

!===============================================================================
!> Unpack patch id from an MPI patch_t buffer
!===============================================================================
INTEGER FUNCTION unpack_id (mesg)
  class(mesg_t), pointer:: mesg
  !.............................................................................
  type(header_t):: header
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%unpack_id')
  call anonymous_copy (n_header, mesg%buffer, header)
  unpack_id = header%id
  call trace%end ()
END FUNCTION unpack_id

!===============================================================================
!> Unpack patch info from a buffer, for MPI communication
!===============================================================================
SUBROUTINE unpack (self, mesg)
  class(patch_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  type(header_t):: header
  integer:: n_buf, ibuf, iv, it
  integer, dimension(:), pointer:: buffer
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%unpack', itimer=itimer)
  call self%set (bits%busy)
  self%unpack_time = wallclock()
  buffer => mesg%buffer
  n_buf = product(self%gn)      ! 4-bytes per element; `anonymous_copy` also uses 4-bytes.
  if (kind(self%mem)==8) &
    n_buf = n_buf * 2           ! 8-bytes per element, but `anonymous_copy` uses 4-bytes!
  !-----------------------------------------------------------------------------
  ! Copy relevant patch info from sequenced header
  !-----------------------------------------------------------------------------
  ibuf = 1
  call anonymous_copy (n_header, buffer(ibuf), header)
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Lock the %new slot, which will become the %it slot after header_to_patch
  !-----------------------------------------------------------------------------
  call self%header_to_patch (header)
  if (io%verbose > 1) &
    write (io_unit%log,'(f12.6,3x,a,i7,2i4,2x,3l1,2i4,f8.4,f11.6)') &
      wallclock(), 'header_to_patch: id, nq, BVS, nsnd, nrcv, latency, time =', &
      self%id, mesg%sender, self%nq, &
      self%is_set(bits%boundary), &
      self%is_set(bits%virtual), &
      self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n, self%latency, self%time
  self%wc_last = wallclock()
  if (io%verbose > 1) write (io_unit%log,1) 'patch_mod::unpk id, nq, time  after:', &
    self%id, self%nq, self%time, &
    self%is_set (bits%internal), &
    self%is_set (bits%boundary), &
    self%is_set (bits%virtual), &
    self%is_set (bits%external), &
    self%is_set (bits%swap_request), &
    self%wc_last
    1 format(a,2i9,f12.6,2x,5l1,f12.1)
  !-----------------------------------------------------------------------------
  ! Copy the variables from the output buffer.  If this is swap of patch roles
  ! expect at least one more time slot.
  !-----------------------------------------------------------------------------
  call self%mem_lock(self%it)%set ('unpack')
  do iv=1,self%nv
    call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,:,iv,self%it,1))
!write(io_unit%log,'(a,i6,2i4,1p,e14.5)') &
!'patch_t%unpack: id, iv, it, it, minval =', &
!self%id, iv, self%it, minval(self%mem(:,:,:,iv,it,1))
    ibuf = ibuf + n_buf
  end do
  call self%mem_lock(self%it)%unset ('unpack')
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%virtual)) then
    it = mod(self%it-1+self%nt-1,self%nt)+1
    call self%mem_lock(it)%set ('unpack')
    if (io%verbose>0) then
      write (io_unit%log,*) self%id, 'extra time slot unpacked', it
      flush (io_unit%log)
    end if
    do iv=1,self%nv
      call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,:,iv,it,1))
      ibuf = ibuf + n_buf
    end do
    call self%mem_lock(it)%unset ('unpack')
  end if
  if (io%verbose > 2) &
    write (io_unit%log,*) self%id, 'after unpack', minval(self%mem(:,:,:,1,self%it,1))
  call self%clear (bits%busy)
  call trace%end (itimer)
END SUBROUTINE unpack

!===============================================================================
!> Print info to stdout.  If
!>   verbose==0:        print only id=1 to stdout and io_unit%log
!>   verbose==1:        print only id=1 to stdout but all to io_unit%log
!>   verbose==2:        print iv=1 to stdout and io_unit%log for all id
!>   verbose==3:        print iv=1,nv to stdout and io_unit%log for all id
!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(patch_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !.............................................................................
  integer:: iv, nv, iprint_l
  integer:: it0, it1, it
  real:: fmin, fmax
  real, dimension(:,:,:), pointer:: df
  integer, parameter:: max_lines=50
  integer, save:: counts(6)=0, id_prv=0
  integer:: counts_l(6), i
  character(len=120):: fmt
  real(8):: time, print_next
  logical:: print_it
  type(lock_t), save:: lock
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%info', itimer=itimer)
  if (io_unit%do_validate) then
    it0 = merge (self%iit(self%nt-2), 1, self%nt>2)
    it1 = merge (self%iit(self%nt-1), 1, self%nt>1)
    if (self%nw==1) then
      allocate (df(size(self%mem,1),size(self%mem,2),size(self%mem,3)))
    end if
    do iv=1,self%nv
      if (self%nw==1) then
        df = self%mem(:,:,:,iv,it1,1) &
           - self%mem(:,:,:,iv,it0,1)
        fmin = self%fminval(df)
        fmax = self%fmaxval(df)
      else
        fmin = self%fminval(self%mem(:,:,:,iv,it0,self%nw))
        fmax = self%fmaxval(self%mem(:,:,:,iv,it0,self%nw))
      end if
      write(io%output,2) &
      'upd',self%id, self%level, iv, self%istep, self%time-self%dtime, &
      self%dtime, self%u_max, fmin, fmax, &
      self%fminval(self%mem(:,:,:,iv,it0,1)), &
      self%fmaxval(self%mem(:,:,:,iv,it0,1)), &
      real(timer%n_update), ntask, nq, self%is_set(bits%boundary), self%is_set(bits%virtual)
    2 format(a,i8,2i4,i3,1p,e14.5,e12.4,0p,f7.2,2x,1p,2e12.4,1x,2e12.4,2x,e10.2,2i5,l3,l1,o12)
    end do
    if (self%nw==1) then
      deallocate (df)
    end if
    call trace%end (itimer)
    return
  end if
  1 format(a,i8,2i4,i3,f12.6,1p,e12.4,0p,f7.2,2x,1p,2e10.2,2x,2e11.3,2x,e10.2,9i5,2i4,l3,l1,o12)
  if (self%time < 1d4) then
    fmt = '(a,i8,2i4,i3,f12.6,1p,e12.4,0p,f7.2,2x,1p,2e10.2,2x,2e11.3,2x,e10.2,i6,3i5,i6,4i5,2i4,l3,l1,o12)'
  else
    fmt = '(a,i8,2i4,i3,1p,2e12.5,     0p,f7.2,2x,1p,2e10.2,2x,2e11.3,2x,e10.2,i6,3i5,i6,4i5,2i4,l3,l1,o12)'
  end if
  !-----------------------------------------------------------------------------
  if ((io%verbose>=0).and.(self%is_clear(bits%no_density))) then
    !$omp atomic
    io%dmin = min(io%dmin, real(self%fminval(self%mem(:,:,:,self%idx%d,self%it,1)),kind=4))
    !$omp atomic
    io%dmax = max(io%dmax, real(self%fmaxval(self%mem(:,:,:,self%idx%d,self%it,1)),kind=4))
    !---------------------------------------------------------------------------
    ! After coming back to the saved ID, reset it, so we stop searching
    !---------------------------------------------------------------------------
    if (self%id == id_prv) then
      !$omp atomic write
      id_prv = 0
    end if
    !---------------------------------------------------------------------------
    ! Allow print every print_every update, if print_every > 0
    !---------------------------------------------------------------------------
    if (print_every > 0) then
      !$omp atomic read
      iprint_l = iprint
      if (iprint_l >= print_every .and. print_every > 0) then
        !$omp atomic write
        iprint = 1
        print_it = .true.
      else
        !$omp atomic update
        iprint = iprint+1
        print_it = .false.
      end if
    !---------------------------------------------------------------------------
    ! If print_every==0, then check if print_time controls cadence, using
    ! an if-lock-if sequence to ensure only one thread prints and updates
    !---------------------------------------------------------------------------
    else if (io%print_time > 0.0) then
      if (print_every == 0 .and. self%time >= io%print_next) then
        call lock%set ('print')
        if (self%time >= io%print_next) then
          print_it = .true.
          print_next = (nint(self%time/io%print_time)+1)*io%print_time
          io%print_next = print_next
        end if
        call lock%unset ('print')
      end if
    !---------------------------------------------------------------------------
    ! Otherwise, fall back on printing the task update with the smallest dt
    !---------------------------------------------------------------------------
    else if (self%id==io%id) then
      print_it = .true.
    end if
    if (print_it .or. io%verbose>0) then
      nv = 1
      if (io%verbose>2) nv = self%nv
      !-------------------------------------------------------------------------
      counts_l = &
        [mpi_mesg%n_check, mpi_mesg%n_ready, mpi_mesg%n_update, &
         mpi_mesg%n_send , mpi_mesg%n_recv , mpi_mesg%n_unpk  ]
      do i=1,size(counts)
        !$omp atomic update
        counts(i) = counts(i) + counts_l(i)
      end do
      !$omp atomic update
      timer%n_lines = timer%n_lines-1
      if (timer%n_lines==0) then
        !$omp critical (info_cr2)
        !$omp atomic write
        timer%n_lines = max_lines
        !$omp end atomic
        write(io%output,'(a)') &
        '0....+....1....+....2....+....3....+....4....+....5....+....' &
        //'6....+....7....+....8....+....9....+....0....+....1  ntsk   nq ' &
        //' chk  rdy   upd  snt  prb  unp  snq rcq unq  BV'
        !$omp end critical (info_cr2)
      end if
      !-----------------------------------------------------------------------
      ! time-dtime is the task time when update started; it is guaranteed to
      ! be increasing monotomically, since the ready queue is time sorted.
      ! To allow using the column for plotting, we show this time when all
      ! times are printed.
      !-----------------------------------------------------------------------
      if (io%print_time == 0.0 .and. io%verbose > 0) then
        time = self%time - self%dtime
      else
        time = self%time
      end if
      do iv=1,nv
        write(io%output,fmt) &
          'upd',self%id, self%level, iv, omp_mythread, time, &
          self%dtime, log10(max(self%u_max,1e-30)), &
          io%dmin, io%dmax, &
          self%fminval(self%mem(:,:,:,iv,self%it ,1)), &
          self%fmaxval(self%mem(:,:,:,iv,self%it ,1)), &
          real(timer%n_update), ntask, nq, counts, &
          mpi_mesg%nq_send, mpi_mesg%nq_recv, mpi_mesg%unpk_list%n, &
          self%is_set(bits%boundary), self%is_set(bits%virtual)
      end do
      flush (io%output)
      !$omp atomic write
      io%dmin = huge(1.)
      !$omp atomic write
      io%dmax = tiny(1.)
      !$omp atomic write
      io%dtime = huge(1d0)
      !$omp atomic write
      id_prv = io%id
      counts(:) = 0
      !-----------------------------------------------------------------------
      ! Reset the message counters after every timestep, if relevant
      !-----------------------------------------------------------------------
      mpi_mesg%n_check=0; mpi_mesg%n_ready=0; mpi_mesg%n_update=0
      mpi_mesg%n_send=0 ; mpi_mesg%n_recv=0 ; mpi_mesg%n_unpk=0
      mpi_mesg%n_delay=0
    end if
    !---------------------------------------------------------------------------
    ! Keep track of the shortest times step between updates of the io%id local
    ! task.  Use this task in patch%info, rather than the one with ip==1.  As
    ! long as shorter dtime are being found, io%id will not match the next self%id.
    ! But if the same io%id survices one update cycle, it will match.
    !---------------------------------------------------------------------------
    if (self%dtime < io%dtime .and. id_prv /= 0) then
      !$omp atomic write
      io%id = self%id
      !$omp atomic write
      io%dtime = self%dtime
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE info

!===============================================================================
!> Return true if a point is inside the self interior
!===============================================================================
LOGICAL FUNCTION contains1 (self, point)
  class(patch_t):: self
  real(8):: point(3)
  contains1 = all(abs(point-self%position) < 0.5_8*self%size)
END FUNCTION contains1

!===============================================================================
!> Return true if a patch position is inside the self interior
!===============================================================================
LOGICAL FUNCTION contains2 (self, patch)
  class(patch_t):: self, patch
  contains2 = contains1 (self, patch%position)
END FUNCTION contains2

!===============================================================================
!> Check if there is overlap, using the size plus an extra cell margin. The
!> definition of overlap MUST be symmetric, since it is being used to define
!> neighbor relations.  This function overloads task_t%overlaps, in order to
!> force use of the patch_t%distance function.
!===============================================================================
FUNCTION overlaps (self, task)
  class(patch_t):: self
  class(task_t):: task
  logical:: overlaps
  real(8):: distance, box, limit
  integer:: i
  !.............................................................................
  !call trace%begin ('patch_t%overlaps', 2)
  select type (task)
  class is (patch_t)
    overlaps = .true.
    do i=1,3
      distance = self%centre_nat(i) - task%centre_nat(i)
      box = self%box(i)
      if (self%periodic(i)) then
        distance = modulo (distance+0.5_8*box, box) - 0.5_8*box
      end if
      limit = 0.5_8*(self%size(i)+task%size(i)+self%ds(i)+task%ds(i))
      if (abs(distance) > limit) then
        overlaps = .false.
        exit
      end if
    end do
  class default
    overlaps = self%task_t%overlaps(task)
  end select
  !call trace%end ()
END FUNCTION overlaps

!===============================================================================
!> Signed distance between the centers of two tasks in the periodic box
!> This version is adjusted to account for curvilinear coords. (the position of
!> a curvilinear patch is set to the coordinate *origin) and requires components
!> defined in `patch_t`.
!===============================================================================
FUNCTION distance_curvilinear (self, task) RESULT (out)
  class(patch_t):: self
  class(task_t):: task
  real(8):: out(3)
  !.............................................................................
  call trace%begin ('patch_t%distance_curvilinear',2)
  !write (io_unit%log,*) 'distance'; flush(io_unit%log)
  select type (task)
  class is (patch_t)
    out = self%centre_nat-task%centre_nat
    !write (io_unit%log,*) out; flush(io_unit%log)
  class default
    out = self%position-task%position
  end select
  ! account for periodic wrapping
  !write (io_unit%log,*) self%box; flush(io_unit%log)
  where (self%periodic .and. self%box > 0d0)
    out = modulo (out+0.5_8*self%box, self%box) - 0.5_8*self%box
  end where
  if (io%verbose>1 .and. (task%id==io%id_debug .or. self%id==io%id_debug)) then
    write(io%output,1) &
      'distance: self', self%id, self%position, self%size, self%box
    write(io%output,1) &
      'distance: task', task%id, task%position, task%size, task%box, out
    1 format(a,i5,4(3x,3f10.2))
  end if
  call trace%end()
END FUNCTION distance_curvilinear

!===============================================================================
!> The nearest smaller index position in the self patch coordinate system.  The
!> point with index self%mesh%lb is half a cell inside the patch guard size, so the
!> floating point position is counted relative to that point, and is split into
!> an integer part and a remainder.  For patches smaller than the box, the index
!> may become both smaller than self%mesh%lb and larger than self%mesh%ub.  For the root
!> grid patch, the index is in the range (li-1,ui) inside the interior.
!===============================================================================
FUNCTION index (self, other, oindex, f)
  class(patch_t):: self, other
  real, dimension(3):: position, f, pos, pos_lb
  integer, dimension(3):: index, oindex
  !.............................................................................
  pos = other%position + (oindex-other%offset)*other%ds         ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  index = f                                                     ! smaller integer
  f = f - index                                                 ! remainder
END FUNCTION index

!===============================================================================
!> The nearest index position in the self patch of the point in another patch.
!> The cost of this operation is 4 plus and 2 multiply per direction = 18 flops
!===============================================================================
FUNCTION nearest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset)*other%ds          ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! nearest integer
END FUNCTION nearest

!===============================================================================
!> The lowest index position in the self patch of the point in another patch
!===============================================================================
FUNCTION lowest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset-0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! lowest integer
END FUNCTION lowest
!===============================================================================
!> The highest index position in the self patch of the point in another patch
!===============================================================================
FUNCTION highest (self, other, index) RESULT (out)
  class(patch_t):: self, other
  integer, dimension(3):: out, index
  real, dimension(3):: position, pos_lb, f, pos
  !.............................................................................
  pos = other%position + (index-other%offset+0.49)*other%ds     ! position in other
  pos = pos-self%position                                       ! relative
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! highest integer
END FUNCTION highest

!===============================================================================
!> Courant condition, one way or another
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(patch_t):: self
  logical, optional:: detailed_timer
  real(8):: dtime, ds(3), courant
  integer:: i2
  integer, save:: itimer=0
  logical, parameter:: save_dbg=.false.
  !.............................................................................
  call trace%begin('patch_t%courant_condition', itimer=itimer)
  if (io%verbose>1) &
    write (io_unit%log,*) self%id, ' courant: u_max, ds', self%u_max, minval(self%ds)
  ds(1) = self%ds(1)
  ds(2) = self%ds(2)*self%mesh(1)%h2c(self%mesh(1)%li)
  ds(3) = self%ds(3)*self%mesh(1)%h31c(self%mesh(1)%li)*abs(self%mesh(2)%h32c(self%mesh(2)%li))
  if (trim(self%mesh(3)%type) == 'spherical_mesh') then
    do i2=self%mesh(2)%li,self%mesh(2)%ui
      ds(3) = min(ds(3),self%ds(3)*self%mesh(1)%h31c(i2)*abs(self%mesh(2)%h32c(i2)))
    end do
  end if
  if (self%dt_fixed > 0.0) then
    dtime = self%dt_fixed
    courant = dtime*self%u_max/minval(ds)
  else if (self%u_fixed > 0.0) then
    dtime = self%courant*minval(ds)/self%u_fixed
    courant = self%courant
  else
    dtime = self%courant*minval(ds)/self%u_max*self%safety_factor
    courant = self%courant
  endif
  !$omp atomic write
  self%dtime = dtime
  if (self%nt > 1 .and. self%istep > 5) then
    if (dtime < 0.5*self%dt(self%iit(self%nt-1)) .and. self%dt_fixed < 0.0) then
      self%safety_factor = 0.5
    end if
    if (self%safety_factor < 0.6 .and. save_dbg) then
      !$omp critical (fort88_cr)
      if (self%safety_factor == 0.5) then
        write (io_unit%dbg) self%id, self%gn, self%ipos
        write (io_unit%dbg) self%t(self%nt-3), self%mem(:,:,:,:,self%iit(self%nt-3),1)
        write (io_unit%dbg) self%id, self%gn, self%ipos
        write (io_unit%dbg) self%t(self%nt-2), self%mem(:,:,:,:,self%iit(self%nt-2),1)
      end if
      write (io_unit%dbg) self%id, self%gn, self%ipos
      write (io_unit%dbg) self%time, self%mem(:,:,:,:,self%it,1)
      flush (io_unit%dbg)
      if (io%verbose>2 .or. (self%istep>5)) &
        write(io%output,'(a,i6,i3,f12.6,1p,3g11.2,0p,1x,3f12.6,1x,1p,7g12.3)') &
          ' courant: ', self%id, self%it, &
          self%time, dtime, courant, self%u_max, self%position, self%dt(self%iit)
      !$omp end critical (fort88_cr)
    end if
    self%safety_factor = self%safety_factor**0.9
    self%get_umax_location = (self%safety_factor < 0.9)
  end if
  !-----------------------------------------------------------------------------
  ! Force the time step to give exactly sync_time if within reach
  !-----------------------------------------------------------------------------
  if (self%time < self%sync_time .and. self%time+self%dtime > self%sync_time) then
    self%dtime = self%sync_time - self%time
    self%syncing = .true.
    if (io%verbose>2) &
      write (io_unit%mpi,*) 'patch_t%courant_condition: syncing turned on for', &
        self%id, self%time
  end if
  !-----------------------------------------------------------------------------
  ! Prevent large (cheap) patches from taking time steps longer than the out_time
  !-----------------------------------------------------------------------------
  if (io%do_output .and. io%out_time > 0.0) &
    self%dtime = min (self%dtime, io%out_time)
  self%dt(self%it) = self%dtime
  if (io%verbose>3) &
    write(io%output,'(a,i7,f12.6,f8.4,1p,e10.2)') 'courant:', &
      self%id, self%dtime, self%courant, self%u_max
  call trace%end (itimer)
END SUBROUTINE courant_condition

!===============================================================================
FUNCTION fsum4 (self, f) RESULT (out)
  class(patch_t):: self
  real:: out
  real, dimension(:,:,:), intent(in):: f
  integer:: n(3)
  !-----------------------------------------------------------------------------
  if (self%is_periodic()) then
    n = self%ncell-1
    out = sum(real(f(self%mesh(1)%li:self%mesh(1)%li+n(1), &
                     self%mesh(2)%li:self%mesh(2)%li+n(2), &
                     self%mesh(3)%li:self%mesh(3)%li+n(3)),kind=8))
  else
    out = sum(real(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui),kind=8))
  end if
END FUNCTION fsum4
!===============================================================================
FUNCTION fsum8 (self, f) RESULT (out)
  class(patch_t):: self
  real(8):: out
  real(8), dimension(:,:,:), intent(in):: f
  integer:: n(3)
  !-----------------------------------------------------------------------------
  if (self%is_periodic()) then
    n = self%ncell-1
    out = sum(real(f(self%mesh(1)%li:self%mesh(1)%li+n(1), &
                     self%mesh(2)%li:self%mesh(2)%li+n(2), &
                     self%mesh(3)%li:self%mesh(3)%li+n(3)),kind=8))
  else
    out = sum(real(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui),kind=8))
  end if
END FUNCTION fsum8

!===============================================================================
FUNCTION faver4 (self, f) RESULT (out)
  class(patch_t):: self
  real:: out
  real, dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  out = self%fsum4(f)/product(self%ncell)
END FUNCTION faver4
!===============================================================================
FUNCTION faver8 (self, f) RESULT (out)
  class(patch_t):: self
  real(8):: out
  real(8), dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  out = self%fsum8(f)/product(self%ncell)
END FUNCTION faver8

!===============================================================================
FUNCTION fminvali (self, iv) RESULT (fminval)
  class(patch_t):: self
  integer:: iv
  real:: fminval
  !-----------------------------------------------------------------------------
  fminval = self%fminvalvar (self%mem(:,:,:,iv,self%it,1))
END FUNCTION fminvali
!===============================================================================
FUNCTION fmaxvali (self, iv) RESULT (fmaxval)
  class(patch_t):: self
  integer:: iv
  real:: fmaxval
  !-----------------------------------------------------------------------------
  fmaxval = self%fmaxvalvar (self%mem(:,:,:,iv,self%it,1))
END FUNCTION fmaxvali
!===============================================================================
FUNCTION fmaxval4 (self, f, outer) RESULT (fmaxval)
  real:: fmaxval
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  logical, optional:: outer
  logical:: outer_l
  !-----------------------------------------------------------------------------
  if (present(outer)) then
    outer_l = outer
  else
    outer_l = .false.
  end if
  if (self%boundaries%is_set (bits%spherical)) then
    fmaxval = maxval(f, mask=self%boundaries%mask(self%mesh))
  else if (outer_l) then
    fmaxval = maxval(f(self%mesh(1)%lo:self%mesh(1)%uo, &
                       self%mesh(2)%lo:self%mesh(2)%uo, &
                       self%mesh(3)%lo:self%mesh(3)%uo))
  else
    fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fmaxval4
!===============================================================================
FUNCTION fminval4 (self, f) RESULT(fminval)
  real:: fminval
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%spherical)) then
    fminval = minval(f, mask=self%boundaries%mask(self%mesh))
  else
    fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fminval4
!===============================================================================
FUNCTION fmaxval8 (self, f, outer) RESULT (fmaxval)
  real(8):: fmaxval
  class(patch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  logical, optional:: outer
  logical:: outer_l
  !-----------------------------------------------------------------------------
  if (present(outer)) then
    outer_l = outer
  else
    outer_l = .false.
  end if
  if (self%boundaries%is_set (bits%spherical)) then
    fmaxval = maxval(f, mask=self%boundaries%mask(self%mesh))
  else if (outer_l) then
    fmaxval = maxval(f(self%mesh(1)%lo:self%mesh(1)%uo, &
                       self%mesh(2)%lo:self%mesh(2)%uo, &
                       self%mesh(3)%lo:self%mesh(3)%uo))
  else
    fmaxval = maxval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fmaxval8
!===============================================================================
FUNCTION fminval8 (self, f) RESULT(fminval)
  real(8):: fminval
  class(patch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  if (self%boundaries%is_set (bits%spherical)) then
    fminval = minval(f, mask=self%boundaries%mask(self%mesh))
  else
    fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                       self%mesh(2)%li:self%mesh(2)%ui, &
                       self%mesh(3)%li:self%mesh(3)%ui))
  end if
END FUNCTION fminval8

!===============================================================================
!> Check for NaN
!===============================================================================
SUBROUTINE check_nan (self, panic, label, always)
  class(patch_t):: self
  logical:: panic
  character(len=*), optional:: label
  logical, optional:: always
  !.............................................................................
  integer:: ix, iy, iz, iv, ii(3)
  logical:: nan
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. do_check_nan) then
    if (io%verbose < 2 .and. .not. present(always)) return
  end if
  if (.not. associated(self%mem)) return
  call timer%begin ('patch_t%check_nan', itimer=itimer)
  if (io%verbose>1) then
    if (present(label)) then
      write(io%output,*) self%id, trim(label)//', min(d), max(d):', &
        minval(self%mem(:,:,:,self%idx%d,self%it,1)), &
        maxval(self%mem(:,:,:,self%idx%d,self%it,1))
    else
      write(io%output,*) self%id, 'min(d), max(d):', &
        minval(self%mem(:,:,:,self%idx%d,self%it,1)), &
        maxval(self%mem(:,:,:,self%idx%d,self%it,1))
    end if
  end if
  nan = .false.
  do iv=1,self%nv
  do iz=1,self%gn(3)
  do iy=1,self%gn(2)
  do ix=1,self%gn(1)
    if (isinf(self%mem(ix,iy,iz,iv,self%it,1))) then
      if (io%verbose>0) then
        if (present(label)) then
          write(io%output,1) trim(label)//' ERROR: NaN id,ix,iy,iz,iv:', &
          self%id,ix,iy,iz,iv
          1 format(a,i6,3i4,i5)
        else
          write(io%output,1) 'ERROR: NaN ix,iy,iz,iv:', ix,iy,iz,iv
        end if
        flush (io%output)
      end if
      ii = [ix,iy,iz]
      nan = .true.
    end if
  end do
  end do
  end do
  end do
  if (nan) then
    flush (io%output)
    panic = .true.
    if (present(label)) then
      write(io%output,*) self%id, trim(label)//' ERROR NaN, min(d), max(d):', &
        minval(self%mem(:,:,:,self%idx%d,self%it,1)), &
        maxval(self%mem(:,:,:,self%idx%d,self%it,1)), ii
    else
      write(io%output,*) self%id, ' ERROR NaN, min(d), max(d):', &
        minval(self%mem(:,:,:,self%idx%d,self%it,1)), &
        maxval(self%mem(:,:,:,self%idx%d,self%it,1)), ii
    end if
    write(io%output,*) self%unsigned
    flush (io%output)
    call io%abort ('found NaN -- check rank & thread logs!')
  end if
  call timer%end (itimer)
END SUBROUTINE check_nan

LOGICAL FUNCTION isinf(a)
  use ieee_arithmetic, only : ieee_is_nan
  real(kind=KindScalarVar):: a
  isinf = ieee_is_nan(a) .or. (abs(a) > huge(1.0_KindScalarVar))
END FUNCTION isinf

!===============================================================================
!> Check if the density is negative somewhere, counting points in the interior
!> and guard zones separately.  Attempt to patch the guard zones, if applicable.
!===============================================================================
SUBROUTINE check_density (self, label)
  class(patch_t):: self
  character(len=*), optional:: label
  real(kind=KindScalarVar):: mind, posd
  real:: tiny = 1.0e-32
  integer:: ix, iy, iz, iv, ii(3), ntot, nins, l(3), u(3), it, jt, iw
  logical:: inside, nan, panic
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if ((io%verbose < 2 .or. io%print_time > 0.0 .or. print_every > 1) .and. &
      .not.do_check_nan) return
  if (self%is_set(bits%no_density+bits%frozen)) return
  if (.not. associated(self%mem)) return
  call timer%begin ('patch_t%check_density', itimer=itimer)
  l = self%mesh%lb
  u = self%mesh%ub
  mind = minval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,1))
  if (present(label) .and. io%verbose>1) write (io_unit%log,*) label, self%id, mind
  panic = .false.
  if (mind <= 0.0) then
    panic = .true.
    !---------------------------------------------------------------------------
    ! Negative densities appear to develop at a few points in the guard zones,
    ! so this tries to patch up, sampling densitites in the guard zone.
    !---------------------------------------------------------------------------
    posd = tiny
    ntot = 0
    nins = 0
    do iz=1,self%gn(3)
    do iy=1,self%gn(2)
    do ix=1,self%gn(1)
      ii = [ix,iy,iz]
      inside = all(ii>self%mesh%lo .and. ii<self%mesh%uo)
      posd = merge (posd, max(posd,self%mem(ix,iy,iz,self%idx%d,self%it,1)), inside)
      ntot = ntot + merge(1,0,self%mem(ix,iy,iz,self%idx%d,self%it,1)<=0.0)
      nins = nins + merge(1,0,self%mem(ix,iy,iz,self%idx%d,self%it,1)<=0.0 .and. inside)
      !-------------------------------------------------------------------------
      ! Pick the sampled positive value if a negative one is found in the
      ! guard zone
      !-------------------------------------------------------------------------
      !self%mem(ix,iy,iz,1,self%it,1) = merge (self%mem(ix,iy,iz,1,self%it,1), &
      !  posd, self%mem(ix,iy,iz,1,self%it,1) > 0.0 .or. inside)
    end do
    end do
    end do
    1 format(a,i8,3x,a,i8,a,i8,a)
    if (present(label)) then
      if (nins > 0) then
        write (io_unit%output,1) 'ID ',self%id, &
         'ERROR: density non-positive in',ntot,' points,',nins,' inside, '// &
         trim(label)
      else if (io%verbose>0) then
          write (io_unit%output,1) 'ID ',self%id, &
           'WARNING: density non-positive in',ntot,' points,',nins,' inside, '// &
           trim(label)
      end if
    else
      if (nins > 0) then
        write (io_unit%output,1) 'ID ',self%id, &
         'WARNING: density non-positive in',ntot,' points,',nins,' inside'
      else if (io%verbose>0) then
        write (io_unit%output,1) 'ID ',self%id, &
         'WARNING: density non-positive in',ntot,' points,',nins,' inside'
      end if
    end if
    if (nins > 0) &
      panic = .true.
  end if
  call self%check_nan (panic, label)
  !-----------------------------------------------------------------------------
  ! If error indicated by check_density or check_nan, dump the patch and abort
  !-----------------------------------------------------------------------------
  if (panic) then
    !$omp critical (panic_cr)
    write (io_unit%dump) self%id, self%gn, self%nv, self%nt, self%nw, self%istep
    write (io_unit%dump) self%time, self%dtime
    do iw=1,self%nw
    do it=1,self%nt
      !$omp atomic read
      jt = self%iit(it)
      !$omp end atomic
      write (io_unit%dump) self%mem(:,:,:,:,jt,iw)
    end do
    end do
    !$omp end critical (panic_cr)
    if (nins > 0 .or. io%verbose > 0) &
      print *, mpi%rank, self%id, &
      'WARNING: found non-positive density points -- check rank and thread logs'
    if (nins > 0) then
      write (stderr,*) 'check dump, rank and thread logs for rank', mpi%rank
      flush (stdout)
      flush (stderr)
      flush (io_unit%log)
      flush (io_unit%dump)
      call io%abort ('found non-positive internal density points')
    end if
  end if
  call timer%end (itimer)
END SUBROUTINE check_density

!===============================================================================
!> Statistics for all vars in a patch
!===============================================================================
SUBROUTINE stats_patch (self, label)
  class (patch_t):: self
  character(len=*):: label
  integer:: iv
  real(8), dimension(4):: si, sb
  do iv=1,self%nv
    call scalar_stats (self%mem(:,:,:,iv,self%it,1), self%mesh%li, self%mesh%ui, si)
    call scalar_stats (self%mem(:,:,:,iv,self%it,1), self%mesh%lb, self%mesh%ub, sb)
    write (io_unit%log,'(i7,2x,a,i3,1p,2(3x,4g12.4))') &
      self%id, trim(label)//' stats: iv, min, aver, rms, max =', iv, si, sb
  end do
END SUBROUTINE stats_patch

!===============================================================================
!> Statistics for a scalar field
!===============================================================================
SUBROUTINE stats_scalar (self, var, label)
  class (patch_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: var
  character(len=*):: label
  !.............................................................................
  real(8), dimension(4):: si, sb
  !-----------------------------------------------------------------------------
  call scalar_stats (var, self%mesh%li, self%mesh%ui, si)
  call scalar_stats (var, self%mesh%lb, self%mesh%ub, sb)
  write (io_unit%output,'(i7,2x,a,1p,2(3x,4g12.4))') &
    self%id, trim(label)//' stats: min, aver, rms, max =', si, sb
END SUBROUTINE stats_scalar

!===============================================================================
!> Statistics for a vector field
!===============================================================================
SUBROUTINE stats_vector (self, var, label)
  class (patch_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:):: var
  character(len=*):: label
  !.............................................................................
  integer:: iv
  real(8), dimension(4):: si, sb
  !-----------------------------------------------------------------------------
  do iv=1,size(var,4)
    call scalar_stats (var(:,:,:,iv), self%mesh%li, self%mesh%ui, si)
    call scalar_stats (var(:,:,:,iv), self%mesh%lb, self%mesh%ub, sb)
    write (io_unit%output,'(i7,2x,a,i3,1p,2(3x,4g12.4))') &
      self%id, trim(label)//' stats: iv, min, aver, rms, max =', iv, si, sb
  end do
END SUBROUTINE stats_vector

!===============================================================================
!===============================================================================
SUBROUTINE scalar_stats (f, l, u, s)
  real(kind=KindScalarVar), dimension(:,:,:):: f
  integer:: l(3), u(3)
  real(8), dimension(4):: s
  integer:: ix, iy, iz
  s(2) = 0d0
  s(3) = 0d0
  s(1) = f(1,1,1)
  s(4) = f(1,1,1)
  do iz=l(3),u(3)
  do iy=l(2),u(2)
  do ix=l(1),u(1)
    s(2) = s(2) + f(ix,iy,iz)
    s(3) = s(3) + f(ix,iy,iz)**2
    s(1) = min(s(1),f(ix,iy,iz))
    s(4) = max(s(4),f(ix,iy,iz))
  end do
  end do
  end do
  s(2) = s(2)/product(u-l+1)
  s(3) = sqrt(max(0d0,s(3)/product(u-l+1)-s(2)**2))
END SUBROUTINE scalar_stats

!===============================================================================
!> default version; only supports Cartesian coords.
!> see the `solvers/zeus3d` version for curvilinear coord. suppport.
!===============================================================================
SUBROUTINE MapVectors (self, in, target, ii, jj, out)
  class(patch_t) :: self
  class(patch_t), pointer :: target
  real(kind=KindScalarVar), dimension(:), intent(in):: in
  real(kind=KindScalarVar), intent(out):: out(size(in))
  integer, intent(in) :: ii(3), jj(3)
  !.............................................................................
  call trace%begin('patch_t%MapVectors', 2)

  select case (self%mesh_type)
  case (1) ! Cartesian
    if (target%mesh_type == mesh_types%Cartesian) then
      ! no conversion necessary
      out = in
    else
      call mpi%abort('Not implemented error! patch_t%MapVectors')
    end if
  case default ! other; not implemented
    call mpi%abort('Not implemented error! patch_t%MapVectors')
  end select

  call trace%end()
END SUBROUTINE MapVectors

!===============================================================================
FUNCTION beyond_patch_edge (self, idir, position) result (outside)
  class(patch_t) :: self
  integer, intent(in) :: idir
  real(8), intent(in) :: position(3)
  logical :: outside
  real(8) :: edge, curvpos(3)
  !.............................................................................

  if (self%n(idir) == 1) then
    outside = .false.
    return
  endif

  select case (self%mesh_type)
  case (1) ! Cartesian
    edge = 0.5 * self%gsize(idir)
    outside = abs(position(idir)-self%position(idir)) > edge
  case (2) ! spherical
    curvpos = position - self%centre_cart
    curvpos = CartesianToSpherical(curvpos(1), curvpos(2), curvpos(3))
    select case (idir)
    case (1)
      edge = self%gsize(idir) + self%llc_nat(idir)
      outside = abs(curvpos(idir)) > edge .or. abs(curvpos(idir)) < self%llc_nat(idir)
    case (2)
      if (associated(self%mesh)) then
        edge = self%gsize(idir) + self%mesh(2)%xi
      else
        edge = self%gsize(idir)
      end if
      outside = abs(curvpos(idir)) > edge
    case (3)
      if (associated(self%mesh)) then
        edge = self%gsize(idir) + self%mesh(3)%zeta
      else
        edge = self%gsize(idir)
      end if
      outside = abs(curvpos(idir)) > edge
    end select
    outside = abs(curvpos(idir)) > edge
  case (3) ! cylindrical
    curvpos = position - self%centre_cart
    curvpos = CartesianToCylindrical(curvpos(1), curvpos(2), curvpos(3))
    select case (idir)
    case (1)
      edge = 0.5 * self%gsize(idir)
      outside = abs(curvpos(idir)) > edge
    case (2)
      edge = self%gsize(idir) + self%llc_nat(idir)
      outside = abs(curvpos(idir)) > edge .or. abs(curvpos(idir)) < self%llc_nat(idir)
    case (3)
      if (associated(self%mesh)) then
        edge = self%gsize(idir) + self%mesh(3)%zeta
      else
        edge = self%gsize(idir)
      end if
      outside = abs(curvpos(idir)) > edge
    end select
  end select
END FUNCTION beyond_patch_edge

!=======================================================================
!> Interpolate -- only in time for now -- from the MHD time slots save in
!> this module to the time requested from the solve procedure
!=======================================================================
SUBROUTINE log_interpolate (self, time, i, a)
  class (patch_t)                 :: self
  real(8)                         :: time
  integer                         :: i
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: a, b
  !.....................................................................
  integer                         :: jt(2)
  real                            :: pt(2)
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace%begin ('patch_t%log_interpolate', itimer=itimer)
  call self%time_interval (time, jt, pt)
  if (io%verbose > 1) then
    b => self%mem(:,:,:,i,jt(1),1)
    write(io%output,*) self%id, 'log_interp(1)', jt(1), minval(b), maxval(b)
    b => self%mem(:,:,:,i,jt(2),1)
    write(io%output,*) self%id, 'log_interp(2)', jt(2), minval(b), maxval(b)
  end if

  a = log(self%mem(:,:,:,i,jt(1),1))*pt(1) &
    + log(self%mem(:,:,:,i,jt(2),1))*pt(2)
  call trace%end (itimer)
END SUBROUTINE log_interpolate

!=======================================================================
!> Search for a time interval to interpolate in, and return the time
!> slot indices and weight factors
!=======================================================================
SUBROUTINE time_interval (self, time, jt, pt, all)
  class (patch_t)                 :: self
  real(8)                         :: time
  integer                         :: jt(2)
  real                            :: pt(2)
  logical, optional               :: all
  !.....................................................................
  integer                         :: i, j, n, iit(self%nt-1)
  real(8)                         :: t(self%nt-1), dt
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  !call trace%begin ('patch_t%time_interval', 2, itimer=itimer)
  if (self%nt==1) then
    jt = 1
    pt(2) = 0.0
  else
    call self%timeslots (iit, t)
    j = 1
    n = self%nt-2
    !-------------------------------------------------------------------
    ! Seach and set the j index to the last point that has t(j)<time.
    ! Note that this does NOT guarantee that t(j+1) >= time.  In cases
    ! where both t(nt-2) < time and t(n-1) < time, extrapolation in time
    ! will result.  What should NOT happen is that a memory slot that
    ! has not yet received an update is chosen.
    !-------------------------------------------------------------------
    do i=1,n
      if (t(i) < time) j=i
    end do
    jt = iit(j:j+1)
    dt = t(j+1)-t(j)
    if (dt==0.0) then
      pt(2) = 0.0
    else
      pt(2) = (time-t(j))/dt
    end if
  end if
  pt(1) = 1.-pt(2)
  if ((io%verbose>3 .or. io_unit%do_validate) .and. pt(2)>0.0) then
    write(io%output,'(a,i6,i6,1p,e14.6,i3,0p,2f8.4,2x,1p,7e14.6)') &
      'time_interval: times  =', self%id, self%istep, time, j, pt, t
  end if
  !call trace%end (itimer)
END SUBROUTINE time_interval

!===============================================================================
!===============================================================================
FUNCTION index_stagger (self, iv) RESULT (out)
  class(patch_t) :: self
  integer        :: iv
  integer        :: i
  real(8)        :: out(3)
  !-----------------------------------------------------------------------------
  do i=1,3
    out(i) = self%mesh(i)%h(iv)
  end do
END FUNCTION index_stagger

!===============================================================================
!> Return indices and fractions for 3D interpolation.  Only called from
!> patch_mod::interpolation, where a possible issue is to ensure the correct
!> behavior when the source and target points are essentially the same, except
!> for roundoff.  We then want to make sure that the index becomes the same,
!> with fraction vanishing, so there is no possible issue with extrapolations.
!> Hence we should add a small fraction before taking the floor value, and then
!> subtract the same.
!===============================================================================
SUBROUTINE index_space (self, pos, iv, i, p)
  class(patch_t):: self
  real(8)           :: pos(3)
  integer           :: iv
  integer           :: i(3)
  real              :: p(3)
  !.............................................................................
  real              :: h(3)
  real, parameter   :: eps=0.001
  integer, save     :: itimer=0
  !-----------------------------------------------------------------------------
  !call trace%begin ('patch_t%index_space', itimer=itimer)
  h = self%index_stagger(iv)
  p = pos - self%position                               ! relative position
  p = modulo(p+0.5_8*self%mesh%b,self%mesh%b) &         ! wrap
              -0.5_8*self%mesh%b
  p = p/self%mesh%d + self%mesh%o - h                   ! float index
  p = p + eps                                           ! ensure same points ok
  i = p                                                 ! integer part
  i = max(self%mesh%li, &
      min(self%mesh%ui-1,i))                            ! valid range
  p = p - eps - i                                       ! remaining fraction
  !-----------------------------------------------------------------------------
  ! The slight allowance for extrapolation below is to detect it in tests
  !-----------------------------------------------------------------------------
  if (zero_order_extrap .and. .not.io_unit%do_validate) &
    p = min(p,max(p,-0.01),1.01)                        ! 0th order extrapolation
  !call trace%end (itimer)
END SUBROUTINE index_space

!===============================================================================
!> Optimized case, where it is known that target and source (self) overlap, so
!> no periodic mapping is needed, and iv is always present
!===============================================================================
SUBROUTINE index_space_of (self, target, ii, iv, jj, pp)
  class(patch_t)          :: self
  class(patch_t), pointer :: target
  integer                 :: ii(3)
  integer                 :: iv
  real                    :: pp(3)
  integer                 :: jj(3)
  !.............................................................................
  class(mesh_t), pointer  :: mt, ms
  integer                 :: i, j
  real                    :: p
  real, parameter         :: eps=0.001
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  do i=1,3
    mt => target%mesh(i)                                  ! target mesh
    ms => self%mesh(i)                                    ! source mesh
    p = mt%p - ms%p + (ii(i) - mt%o + mt%h(iv))*mt%d      ! real position
    p = p/ms%d + ms%o - ms%h(iv) + eps                    ! float index
    j = p                                                 ! integer part
    j = max(ms%li, &
        min(ms%ui-1,j))                                   ! valid range
    p = p - eps - j                                       ! remaining fraction
    if (zero_order_extrap .and. .not.io_unit%do_validate) &
      p = min(max(p,-0.01),1.01)                          ! 0th order extrapolation
    pp(i) = p
    jj(i) = j
  end do
END SUBROUTINE index_space_of

!===============================================================================
!> Return only the index in self corresponding to a point in native space.
!> "pos" is a relative position in target coordinate space, and since we are
!> not informed about the target grid spacing it must be adjusted for staggering
!> before the call.
!===============================================================================
FUNCTION index_only (self, pos, iv, roundup, nowrap) RESULT (ii)
  class(patch_t)      :: self
  real(8), intent(in) :: pos(3)
  integer, optional   :: iv
  logical, optional   :: roundup, nowrap
  integer             :: ii(3)
  !.............................................................................
  integer             :: i
  real(8)             :: p(3), h(3)
  class(mesh_t), pointer :: mesh
  !-----------------------------------------------------------------------------
  if (present(iv)) &
    h = self%index_stagger(iv)
  do i=1,3
    mesh => self%mesh(i)
    p(i) = pos(i) - self%position(i)                            ! relative position
    if (self%periodic(i)) then
      if (.not.present(nowrap)) &
        p(i) = modulo(p(i)+0.5_8*mesh%b, mesh%b)-0.5_8*mesh%b   ! periodic wrap
    end if
    p(i) = p(i)/mesh%d + self%offset(i)                         ! float index
    if (present(iv)) then
      p(i) = p(i) + mesh%h(iv)                                  ! staggering
    end if
    if (present(roundup)) then
      ii(i) = ceiling(p(i))                                     ! round up
    else
      ii(i) = p(i) + 0.0001                                     ! floor w margin
    end if
  end do
END FUNCTION index_only

!===============================================================================
!> Return indices and fractions for 3D interpolation.  Only called from
!> patch_mod::interpolation, where a possible issue is to ensure the correct
!> behavior when the source and target points are essentially the same, except
!> for roundoff.  We then want to make sure that the index becomes the same,
!> with fraction vanishing, so there is no possible issue with extrapolations.
!> Hence we should add a small fraction before taking the floor value, and then
!> subtract the same.
!===============================================================================
SUBROUTINE index_space2 (self, nbor, pos, iv, ii, p)
  class(patch_t):: self
  class(patch_t), pointer :: nbor
  real(8)           :: pos(3)
  integer           :: iv
  integer           :: ii(3)
  real              :: p(3)
  !.....................................................................
  class(mesh_t), pointer :: m,mn
  integer           :: i
  real(8)           :: pp
  real(8), parameter:: eps=0.001_8
  integer, save     :: itimer=0
  real, parameter:: sml = -1.0e-5
  !-----------------------------------------------------------------------------
  !call trace%begin ('patch_t%index_space', itimer=itimer)
  do i=1,3
    m => self%mesh(i)
    mn => nbor%mesh(i)
    pp = mn%p - m%p                                                     ! relative position
    pp = pp + pos(i) + m%h(iv)*(mn%d-m%d)                               ! position in the local frame
    if (self%periodic(i)) &
      pp = modulo(pp+0.5_8*m%b,m%b)-0.5_8*m%b                           ! wrap
    pp = pp/m%d + m%o                                                   ! float index
    !pp = round(pp,rnd)
    !pp = pp + eps                                                      ! ensure same points ok
    ii(i)= floor(pp+eps)                                                ! integer part
    ii(i) = max(m%lb, &
        min(m%ub,ii(i)))                                                ! valid range
    p(i) = pp - real(ii(i),kind=4)
    if ((p(i)<0.0).and.(p(i)>sml)) p(i) = 0.0
    !p(i) = round(p(i),rnd)                                              ! remaining fraction
  end do
  !call trace%end (itimer)
END SUBROUTINE index_space2

!===============================================================================
!> Return only the index in self corresponding to a point in native space
!===============================================================================
FUNCTION index_only2 (self, nbor, pos, iv, roundup, closest) RESULT (ii)
  class(patch_t)    :: self
  class(patch_t), pointer :: nbor
  real(8)           :: pos(3)
  integer           :: iv
  integer           :: ii(3)
  logical, optional :: roundup
  logical, optional :: closest
  !.............................................................................
  class(mesh_t), pointer :: m, mn
  integer           :: i
  real(8), parameter:: eps=0.001_8
  real(8)           :: p(3)
  !-----------------------------------------------------------------------------
  do i=1,3
    m => self%mesh(i)
    mn => nbor%mesh(i)
    p(i) = mn%p - m%p                                                   ! relative position
    p(i) = p(i) + pos(i) + m%h(iv)*(mn%d-m%d)                           ! position in the local frame
    if (self%periodic(i)) &
      p(i) = modulo(p(i)+0.5_8*m%b,m%b)-0.5_8*m%b                       ! wrap
    p(i) = p(i)/m%d+m%o                                                 ! float index
    !p(i) = round(p(i),m%rnd)
    if (present(roundup)) then
      ii(i) = ceiling(p(i)-eps)                                         ! round up, with epsilon substracted
    else if (present(closest)) then
      ii(i) = nint(p(i))
    else
      ii(i) = floor(p(i)+eps)                                           ! integer part with added epsilon
    end if
    ii(i) = max(m%lb,min(m%ub,ii(i)))                                   ! valid range
  end do
END FUNCTION index_only2

!===============================================================================
!> How many edges do `self` and its neighbor `source` share?
!===============================================================================
FUNCTION count_edges (self, source) RESULT (out)
  class(patch_t)          :: self
  class(patch_t), pointer :: source
  integer:: out, i
  real(8):: dist(3)
  !-----------------------------------------------------------------------------
  dist = self%distance (source)
  out = 0
  do i=1,3
    if  (abs(dist(i)) > self%ds(i)) out = out+1
  end do
END FUNCTION count_edges

!===============================================================================
!> Interpolate in space and time, for a case where it is safe = the source point
!> is inside the valid domain (inside a distance of half a cell from the source
!> patch boundary for non-staggered points)
!===============================================================================
SUBROUTINE interpolate (self, target, ii, iv, jt, pt)
  class(patch_t)           :: self
  class(patch_t), pointer  :: target
  integer                  :: ii(3), iv, jt(2)
  real                     :: pt(2)
  !.............................................................................
  real(8)                  :: pos(3)
  integer                  :: j(3)
  real                     :: p(3), q(3)
  integer, save            :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%pervolume(iv)) then
    call self%interpolate_specific (target, ii, iv, jt, pt)
  else
    pos = target%myposition (ii, iv)
    call self%index_space (pos, iv, j, p)
    q = 1.0-p
    !if (io%verbose > 1) &
    !  write(io_unit%log,'(a,2i6,2x,3i4,2x,3f7.3)') &
    !    'interpolate: target, source, j, p =', target%id, self%id, j, p
    !call trace%begin ('patch_t%interpolate', itimer=itimer)
    if (self%unsigned(iv)) then
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = exp( &
       pt(1)*(q(3)*(q(2)*(q(1)*log(self%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(1),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(1),1)))   + &
                    p(2)*(q(1)*log(self%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(1),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(1),1))))  + &
              p(3)*(q(2)*(q(1)*log(self%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(1),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(1),1)))   + &
                    p(2)*(q(1)*log(self%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(1),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(1),1))))) + &
       pt(2)*(q(3)*(q(2)*(q(1)*log(self%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(2),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(2),1)))   + &
                    p(2)*(q(1)*log(self%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(2),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(2),1))))  + &
              p(3)*(q(2)*(q(1)*log(self%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(2),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(2),1)))   + &
                    p(2)*(q(1)*log(self%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(2),1))    + &
                          p(1)*log(self%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(2),1))))))
    else
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
       pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(1),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(1),1))   + &
                         (   p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(1),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(1),1)))  + &
              (   p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(1),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(1),1))   + &
                         (   p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(1),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(1),1)))) + &
       pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)  ,j(3)  ,iv,jt(2),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)  ,j(3)  ,iv,jt(2),1))   + &
                         (   p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)+1,j(3)  ,iv,jt(2),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)+1,j(3)  ,iv,jt(2),1)))  + &
              (   p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)  ,j(3)+1,iv,jt(2),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)  ,j(3)+1,iv,jt(2),1))   + &
                         (   p(2))*((1.0-p(1))*self%mem(j(1)  ,j(2)+1,j(3)+1,iv,jt(2),1)    + &
                                    (    p(1))*self%mem(j(1)+1,j(2)+1,j(3)+1,iv,jt(2),1))))
    end if
  end if
  if (target%id==io%id_debug .and. io%verbose>2) &
    write(io%output,1) target%id, ii, iv, self%id, j, jt, p, pt, &
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1)
    1 format(i6,1x,3i3,i4,2x,"DBG patch_t%interpolate src:",i6, &
      "   j,p;",3i4,1x,2i4,3f6.2,1x,2f6.2,"   out:",1p,e15.6)
  !call trace%end (itimer)
END SUBROUTINE interpolate

!===============================================================================
!> Interpolate in space and time, for a case where it is safe = the source point
!> is inside the valid domain (inside a distance of half a cell from the source
!> patch boundary for non-staggered points)
!===============================================================================
SUBROUTINE interpolate_specific (self, target, ii, iv, jt, pt)
  class(patch_t)           :: self
  class(patch_t), pointer  :: target
  integer                  :: ii(3), iv, jt(2)
  real                     :: pt(2)
  real                     :: tmp(2,2,2,2)
  !.............................................................................
  real(8)                  :: pos(3)
  integer                  :: j(3)
  real                     :: p(3)
  integer, save            :: itimer=0
  !-----------------------------------------------------------------------------
  pos = target%myposition (ii, iv)
  call self%index_space (pos, iv, j, p)
  !call trace%begin ('patch_t%interpolate', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Pick up 16 values for tri-linear interpolation, divided by density
  !-----------------------------------------------------------------------------
  associate (id=>self%idx%d)
  tmp(:,:,:,:) = self%mem(j(1):j(1)+1,j(2):j(2)+1,j(3):j(3)+1,iv,jt(1:2),1) &
               / self%mem(j(1):j(1)+1,j(2):j(2)+1,j(3):j(3)+1,id,jt(1:2),1)
  target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
   pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,1,1)    + &
                                (    p(1))*tmp(2,1,1,1))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,1,1)    + &
                                (    p(1))*tmp(2,2,1,1)))  + &
          (   p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,2,1)    + &
                                (    p(1))*tmp(2,1,2,1))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,2,1)    + &
                                (    p(1))*tmp(2,2,2,1)))) + &
   pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,1,2)    + &
                                (    p(1))*tmp(2,1,1,2))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,1,2)    + &
                                (    p(1))*tmp(2,2,1,2)))  + &
          (   p(3))*((1.-p(2))*((1.0-p(1))*tmp(1,1,2,2)    + &
                                (    p(1))*tmp(2,1,2,2))   + &
                     (   p(2))*((1.0-p(1))*tmp(1,2,2,2)    + &
                                (    p(1))*tmp(2,2,2,2))))
  !-----------------------------------------------------------------------------
  ! Multiply by target denssity
  !-----------------------------------------------------------------------------
  target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
    target%mem(ii(1),ii(2),ii(3),iv,target%it,1) * &
    target%mem(ii(1),ii(2),ii(3),id,target%it,1)
  end associate
  if (target%id==io%id_debug .and. io%verbose>2) &
    write(io%output,1) target%id, ii, iv, self%id, j, jt, p, pt, &
      target%mem(ii(1),ii(2),ii(3),iv,target%it,1)
    1 format(i6,1x,3i3,i4,2x,"DBG patch_t%interpolate_specific src:", &
      i6,"   j,p;",3i4,1x,2i4,3f6.2,1x,2f6.2,"   out:",1p,e15.6)
  !call trace%end (itimer)
END SUBROUTINE interpolate_specific

SUBROUTINE update_position (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  if (all(self%velocity == 0.0)) return
  call trace%begin("patch_t%update_position")
  !
  self%position = self%position + self%velocity * self%dtime
  self%llc_cart = self%llc_cart + self%velocity * self%dtime
  self%centre_cart = self%centre_cart + self%velocity * self%dtime
  !
  self%mesh(:)%p = self%position(:)
  self%mesh(:)%llc_cart = self%llc_cart(:)
  !
  ! adjust the patch position for periodic wrapping
  where(self%position > 0.5*self%box)
    self%position = self%position - self%box
    self%mesh(:)%p = self%mesh(:)%p - self%box(:)
    self%llc_cart = self%llc_cart - self%box
    self%mesh(:)%llc_cart = self%mesh(:)%llc_cart - self%box
    self%centre_cart = self%centre_cart - self%box
  end where
  call trace%end()
END SUBROUTINE update_position

!===============================================================================
SUBROUTINE custom_refine (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  call io%abort ("patch_t:: custom refine is triggered, but none provided. aborting.")
END SUBROUTINE custom_refine

!===============================================================================
!> Check overlap and return overlap range in a patch, in terms self indices
!===============================================================================
FUNCTION get_overlap (self, patch, guards, l, u) RESULT(overlap)
  class (patch_t):: self
  class (patch_t), pointer:: patch
  integer:: l(3), u(3)
  logical:: guards, overlap
  !.............................................................................
  real(8):: dist(3), la(3), ua(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('download_t%overlap_range')
  overlap = self%overlaps (patch)
  if (overlap) then
    !---------------------------------------------------------------------------
    ! First take care of mapping the distance between the patch and self
    ! centers, periodically, so this can be ignored in the next steps.
    !---------------------------------------------------------------------------
    dist = patch%position - self%position
    where (patch%periodic)
      dist = modulo(dist+0.5_8*patch%mesh%b,patch%mesh%b) - 0.5_8*patch%mesh%b
    end where
    !---------------------------------------------------------------------------
    ! Add the distances from the patch center to the start and end of guard
    ! zones, and convert to indices in self index space, making sure the range
    ! is legal
    !---------------------------------------------------------------------------
    if (guards) then
      la = dist - patch%size/2d0 - (patch%mesh%ng+1)*patch%mesh%d
      ua = dist + patch%size/2d0 + (patch%mesh%ng+1)*patch%mesh%d
    else
      la = dist - patch%size/2d0
      ua = dist + patch%size/2d0
    end if
    l =   floor(la/self%mesh%d + self%mesh%o)
    u = ceiling(ua/self%mesh%d + self%mesh%o)
    l = max(l,self%mesh%lb)
    u = min(u,self%mesh%ub)
    !---------------------------------------------------------------------------
    ! Check that there really is overlap with the authoritative region
    !---------------------------------------------------------------------------
    if (any(l > self%mesh%ui .or. u < self%mesh%li .or. l > u)) then
      write (stderr,'(a,2(3x,a,i6,":",i2),3x,a,3f8.3,3x,a,2(3i4,2x))') &
      'patch_t%overlap_range false overlap, with ', &
      'self:', self%id, self%level, &
      'patch:', patch%id, patch%level, &
      'dist:', dist/self%ds, &
      'l, u:', l, u
      overlap = .false.
    end if
    !---------------------------------------------------------------------------
    ! Account for collapsed dimensions
    !---------------------------------------------------------------------------
    where (self%mesh%n == 1)
      l = min(l,self%mesh%ub)
      u = max(u,self%mesh%lb)
    end where
  end if
  call trace%end()
END FUNCTION get_overlap

!===============================================================================
SUBROUTINE init_level (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  self%level = maxval(nint(log(self%box/self%ds) / &
    log(real(self%refine_ratio))), mask=self%n > 1)
END SUBROUTINE init_level

END MODULE patch_mod
