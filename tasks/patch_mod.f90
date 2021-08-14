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
  USE index_mod
  USE shared_mod
  USE region_mod
  USE lagrange_mod
  implicit none
  private
  integer, parameter:: ndim = 3
  type, extends(task_t), public:: patch_t
    real(kind=KindScalarVar), dimension(:,:,:,:,:,:), pointer, contiguous:: mem      => null()
    real(kind=KindScalarVar), dimension(:,:,:,:,:),   pointer, contiguous:: vgas     => null()
    real(kind=KindScalarVar), dimension(:,:,:,:),     pointer            :: density  => null()
    real(kind=KindScalarVar), dimension(:,:,:,:),     pointer, contiguous:: pressure => null()
    !dir$ attributes align: 64 :: mem, vgas, density, pressure
    !---------------------------------------------------------------------------
    ! These allocatables are intended to deposite external forces and heating in,
    ! typically used by methods called from the extras_t data type procedures.
    !---------------------------------------------------------------------------
    real              , pointer, contiguous:: frefine(:,:,:) => null()
    integer(kind=int8), pointer, contiguous:: filled (:,:,:) => null()
    integer(kind=int8), pointer, contiguous:: request(:,:,:) => null()
    logical           , pointer, contiguous:: outside(:,:,:) => null()
    real(var_k)       , pointer, contiguous:: u_phase(:,:,:) => null()
    !dir$ attributes align: 64 :: filled
    class(mesh_t), dimension(:), pointer:: mesh => null()
    real:: courant                      ! courant number
    real(8):: pdf_next=0.0              ! next PDF output time
    real(8):: dt_fixed                  ! for fixed time step Courant condition
    !real(8):: omega(3)=0d0              ! generic angular velocity
    real(8):: radius=0d0                ! generic radius
    real(8), pointer:: e(:,:) => null() ! local unit vectors
    real(8), pointer:: einv(:,:) => null() ! local unit vectors
    real:: u_max=0.0                    ! for Courant condition
    real:: u_fixed                      ! for fixed time step Courant condition
    real:: d_min, d_max                 ! for info
    integer:: ndim = ndim               ! number of dimensions
    integer:: nv=1                      ! number of variables
    integer:: mv=0                      ! number of memory slots
    integer:: nw=2                      ! number of work indices
    integer:: msplit=0                  ! how to split a patch when refining
    integer:: max_files                 ! prevent disk overrun when NaN happens
    integer:: check_refine_next=0       ! next step to check refinement
    integer:: time_derivs               ! for restart changes
    integer:: not_filled=0              ! guard zones not filled
    integer:: record=0                  ! record number in snapshot file
    logical:: guard_zones               ! for restart changes
    logical:: mem_allocated=.false.
    logical:: do_pebbles=.false.
    logical:: staggered=.true.
    logical:: no_mans_land=.true.
    logical:: use_data_hub=.false.
    logical:: manual_refine=.false.
    logical:: get_umax_location=.false.
    logical, dimension(:), pointer:: unsigned  => null()
    logical, dimension(:), pointer:: pervolume => null()
    class(link_t),         pointer:: link      => null()
    character(len=32):: eos='ideal'
    character(len=32):: opacity='none'
    integer:: file_offset(4)
    type(index_t):: idx
    real:: safety_factor=1.0
    real:: ee_min=tiny(1.0)
    !--- used in extras ---
    logical, contiguous, pointer :: mask(:,:,:) => null()
    !--- used in refine ---
    type(lock_t),        pointer :: request_lock   => null()
    type(lock_t),        pointer :: mem_lock(:)    => null()
    !---------------------------------------------------------------------------
    ! HACKs below:  FIXME
    !---------------------------------------------------------------------------
    real(8):: theta=0d0                 ! latitude
    real(8):: phi=0d0                   ! longitude at t=0.0
    real:: csound                       ! isothermal sound speed, for tests
    real:: grav                         ! gravitational consant (code units)
    integer:: n(3)=0, ng(3)=0           ! defaults, to detect unchanged values
    integer:: li(3), ui(3), gn(3)=0     ! convenience -- avoid patch%mesg%li etc
    integer:: ncell(3)=0                ! number of cells per patch
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
    procedure:: pre_init
    procedure:: setup
    procedure:: reset
    procedure:: dealloc
    procedure:: clone_mem_accounting
    procedure:: allocate_mesg           ! allocate a mesg buffer
    procedure:: pack                    ! pack patch into buffer
    procedure, nopass:: unpack_id       ! unpack patch id
    procedure:: unpack                  ! unpack buffer into patch
    procedure:: index                   ! compute the index for a single point
    procedure:: nearest                 ! the nearest index that maps to patch point
    procedure:: lowest                  ! the lowest index that maps to patch point
    procedure:: highest                 ! the highest index that maps to patch point
    procedure:: myposition
    procedure:: myposition_staggered
    procedure:: local_position
    procedure:: extrapolate_guards
    procedure:: time_extrapolate
    procedure:: is_periodic
    procedure, private :: make_periodic4
    procedure, private :: make_periodic8
    generic, public :: make_periodic  => make_periodic4 , make_periodic8
    procedure, private :: fsum4
    procedure, private :: fsum8
    generic, public :: fsum    => fsum4, fsum8        ! interior sum
    procedure, private :: faver4
    procedure, private :: faver8
    procedure, private :: faveri
    generic, public :: faver    => faver4, faver8, faveri   ! interior average
    procedure, private :: fminloc4
    procedure, private :: fminloc8
    generic, public :: fminloc => fminloc4, fminloc8
    procedure, private :: fmaxloc4
    procedure, private :: fmaxloc8
    generic, public :: fmaxloc => fmaxloc4, fmaxloc8
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
    procedure:: nbor_courant_conditions
    procedure:: check_density
    procedure:: check_nan
    procedure:: custom_refine
    procedure:: patch_to_header
    procedure:: header_to_patch
    procedure, private:: stats_patch
    procedure, private:: stats_scalar
    procedure, private:: stats_scalar8
    procedure, private:: stats_vector
    generic,  public:: stats => stats_patch, stats_scalar, stats_scalar8, stats_vector
    procedure:: MapVectors
    procedure:: beyond_patch_edge
    procedure:: myfloatindex_scalar
    procedure:: periodic_grid
    procedure:: log_interpolate
    procedure:: time_interval
    procedure:: distance
    procedure:: distance_curvilinear
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
    procedure:: dumpm
    procedure, private:: dump1r
    procedure, private:: dump1d
    procedure, private:: dump2r
    procedure, private:: dump2d
    procedure, private:: dump3r
    procedure, private:: dump3d
    generic, public:: dump => dumpm, dump1r, dump1d, dump2r, dump2d, dump3r, dump3d
    procedure:: cost_f
    procedure:: overlap
    procedure:: cover
    procedure:: support_arrays
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
    integer:: rank_steps=0
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
  logical, save:: do_check_nan=.false.              ! check for NaN?
  logical, save:: zero_order_extrap=.true.          ! extrap to no-mans-land
  logical, save:: courant_ramp=.false.              ! ramp up Courant param
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
  head%nq         = shared%n_queue
  head%rank_steps = self%rank_steps
  head%send_time  = wallclock()
  if (io%verbose > 1) &
    write (io_unit%log,'(f12.6,3x,a,i7,i4,l3,2i4,f8.4)') &
      wallclock(), 'patch_to_header: id, nq =', self%id, &
      self%nq, self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n
END SUBROUTINE patch_to_header

!===============================================================================
!> Copy he most important fields from a sequenced header, for communication and I/O
!===============================================================================
SUBROUTINE header_to_patch (self, head)
  class(patch_t):: self
  type(header_t):: head
  integer:: nt, new, i, save_bits
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Swap the roles of boundary and virtual on incoming patches.  Also clear the
  ! ready bit, which is an indicator of a patch being in the ready queue, which
  ! an incoming patch is not.  By doing this here, we avoid having to do it in
  ! a critical region, otherwise needed to protect an operation on the patch
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%header_to_patch', itimer=itimer)
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
  !-----------------------------------------------------------------------------
  ! Lock task while modifying metadata
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! Except for the first and last appearance, istep should increase be one
  !-----------------------------------------------------------------------------
  if (head%id /= self%id) then
    call io%abort('wrong message ID unpacked')
  end if
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
  !-----------------------------------------------------------------------------
  ! If the task is virtual, preserve bits%init_nbors and bits%update_status
  ! across the task update
  !-----------------------------------------------------------------------------
                                        call self%lock%set ('header_to_patch')
  if (self%is_set(bits%virtual)) then
    save_bits = iand(self%status,bits%init_nbors+bits%update_status)
    self%status = ior(self%status,save_bits)
  else
    self%status     = head%status
  end if
  self%id         = head%id
  self%rank       = head%rank
  self%istep      = head%istep
  self%iout       = head%iout
  self%level      = head%level
  self%out_next   = head%out_next
  self%print_next = head%print_next
  self%dtime      = head%dtime
  self%position   = head%position
  self%velocity   = head%velocity
  self%size       = head%size
  self%rank_steps = head%rank_steps
  self%it         = head%it
  self%nt         = head%nt
  self%nq         = head%nq
  self%t          = head%t  (1:self%nt)
  self%dt         = head%dt (1:self%nt)
  self%iit        = head%iit(1:self%nt)
  self%time       = head%time
  !
  if (head%istep == 0) then
    write (io_unit%log,'(1x,a,i6,f12.6,2x,9f12.6)') &
      'patch_t%header_to_patch: id, time, t =', &
      self%id, self%time, self%t
    write (io_unit%log,'(1x,a,i6,2x,3i4,2x,9i4)') &
      'patch_t%header_to_patch: id, istep, it, new, iit =', &
      self%id, self%istep, self%it, self%new, self%iit
  end if
  if (self%it > self%nt .or. self%it < 1) then
    write (io_unit%log,*) 'header_to_patch ERROR: it =', self%it, self%iit(nt-1)
    self%it = self%iit(nt-1)
  end if
  self%nq         = head%nq
  self%latency    = wallclock() - head%send_time
                                        call self%lock%unset ('header_to_patch')
  !-----------------------------------------------------------------------------
  ! Latency stats
  !-----------------------------------------------------------------------------
  !$omp atomic
  timer%latency%max  = max(timer%latency%max,self%latency)
  !$omp atomic
  timer%latency%aver = timer%latency%aver+self%latency
  !$omp atomic
  timer%latency%n    = timer%latency%n+1
  !.............................................................................
  call trace%end (itimer)
END SUBROUTINE header_to_patch

!===============================================================================
!> Cast a generic task_t to patch_t -- the calling code should test for null
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
!> Initialize default values for the patch memory and patch dimensions. This
!> is the 1st stop in the cascading %pre_init() calls, where initial values
!> are assinged to basic attributes.  These may be adjuste by higher level
!> pre_init() calls, after the point where they made the call that ended up
!> here.  The final verdict on all attributes are acted upon by the %init()
!> chain of calls.
!>
!> The order of calls is (should be):
!>  solver_t%pre_init
!>    mhd_t%pre_init
!>      self%nv = ..              ! solver preference = minimum values
!>      extras_t%pre_init
!>        gpatch_t%pre_init
!>          patch_t%pre_init
!>            self%mv = self%nv   ! first chance to copy solver choice
!>        minions_t%pre_init      ! might increase self%nv and/or self%mv
!>  solver_t%init
!>    mhd_t%init
!>      extras_t%init             ! controls minion initialization
!>        gpatch_t%init           ! pass through
!>          patch_t%init          ! allocates memory with %mv slots
!>        minions_t%init          ! initializes requested slots
!===============================================================================
SUBROUTINE pre_init (self)
  class(patch_t), target:: self
  !.............................................................................
  integer       :: iostat, i, n_alloc=0, n_ext=0
  !.............................................................................
  integer, save :: n(3)=16, ng(3)=0, nv=1, nt=5, nw=1, max_files=1000, &
                   print_every=0
  logical, save :: no_mans_land=.true., use_data_hub=.false., limit_dtime=.false.
  real,    save :: u_fixed=-1.0, grace=0.0
  real(8), save :: start_time=0d0, end_time=1.0, out_time=0.1, print_time=0.0
  real(8), save :: dt_fixed=-1.0
  logical, save :: first_time = .true.
  namelist /patch_params/ n, nv, nt, u_fixed, dt_fixed, grace, &
                          no_mans_land, use_data_hub, do_check_nan, &
                          zero_order_extrap, limit_dtime, courant_ramp
  namelist /out_params/ start_time, out_time, end_time, print_time, &
                        print_every, max_files
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%pre_init')
  !-----------------------------------------------------------------------------
  ! Read and set patch parameters; the preferred defalts may have been set by
  ! the calling code -- typically in experiment_mod.f90
  !-----------------------------------------------------------------------------
  call io%lock%set ('patch_t%pre_init')
  if (first_time) then
    nw = self%nw
    nv = self%nv
    ng = self%ng
    rewind (io%input)
    read (io%input, patch_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, patch_params)
    end if
    call io%namelist_warning ('patch_params', iostat)
    if (dt_fixed > 0.0) grace=max(grace,1e-6)
    write (io%namelist, patch_params)
    !...........................................................................
    rewind (io%input)
    read (io%input, out_params,iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, out_params)
    end if
    call io%namelist_warning ('out_params', iostat)
    write (io%namelist, out_params)
    !...........................................................................
    if (mpi%rank > 0) &
      self%track = .true.
    io%print_every = print_every                ! print 1st task
    io%iprint = print_every                     ! print 1st task
    first_time = .false.
  end if
  !-----------------------------------------------------------------------------
  ! Output related parameters
  !-----------------------------------------------------------------------------
  io%start_time    = start_time
  io%end_time      = end_time
  io%out_time      = out_time
  io%print_time    = print_time
  !-----------------------------------------------------------------------------
  ! Multi-threading below -- only self% references allowed
  !-----------------------------------------------------------------------------
  self%n            = n
  self%nt           = nt
  self%ng           = ng                         ! cf warning
  self%nv           = nv                         ! ditto
  self%mv           = self%nv                    ! make sure %mv is set
  self%nw           = nw                         ! ditto
  self%restart      = io%restart                 ! restart snapshot
  self%grace        = grace                      ! grace fraction of dtime
  self%dt_fixed     = dt_fixed
  self%u_fixed      = u_fixed
  self%max_files    = max_files                  ! prevent run-away NaN
  self%use_data_hub = use_data_hub               ! guard zone handling
  self%no_mans_land = no_mans_land               ! mesh centering
  self%print_next   = print_time
  self%limit_dtime  = limit_dtime
  !-----------------------------------------------------------------------------
  ! Preliminary region information, relative to the task position.  This may
  ! change in patch_t%init, but the settings below allow overlap() calls to work
  !-----------------------------------------------------------------------------
  self%authority%l = - self%size*0.5_8
  self%authority%u = + self%size*0.5_8
  self%interest%l  = self%authority%l - self%ng*self%size/self%n
  self%interest%u  = self%authority%u + self%ng*self%size/self%n
  !-----------------------------------------------------------------------------
  allocate (self%request_lock)                   ! allocate lock for frefine
  call self%request_lock%init                    ! initialize lock
  !-----------------------------------------------------------------------------
  ! Preliminary level and hash key assignment
  !-----------------------------------------------------------------------------
  self%ds           = self%size/self%n           ! for init_level
  call self%init_level                           ! redone in %init()
  call self%task_t%pre_init ()                   ! creates hash key
  call io%lock%unset ('patch_t%pre_init')
  !-----------------------------------------------------------------------------
  ! The solver sets the initial values of self%n{w,v,g}, and these are also set
  ! as default values in the namelist.  It is possible to change them here,
  ! but normally this should only be done in connection with test, and one or
  ! more extras modules may also modify any one of them afterwards.
  !-----------------------------------------------------------------------------
  if (any([nw,nv,ng] /= [self%nw,self%nv,self%ng])) then
    write (stdout,'(a,2(2x,2i4,1x,3i4))') &
      'patch_t%pre_init WARNING: nw,nv,ng old/new =', &
      nw,nv,ng, self%nw,self%nv,self%ng
  end if
  call omp%set_stacksize (self%nv)
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Mesh-related patch initialization.  Previously part of the %pre_init-chain
!> of calls, but now part of the %init-chain
!===============================================================================
SUBROUTINE mesh_init (self)
  class(patch_t):: self
  logical, save:: warn=.true.
  integer:: i
  !-----------------------------------------------------------------------------
  ! Initialize the meshes, in order to be able to do the staggering, which is
  ! needed in the mhd_t%init(), before the call to extras_t%init()
  !-----------------------------------------------------------------------------
  call MeshFactory(self%mesh_type, self%mesh, n=self%n, ng=self%ng, &
    s=self%size, p=self%position, llc_native=self%llc_nat, nml=self%no_mans_land)
  !-----------------------------------------------------------------------------
  ! Pass on a few things to the mesh
  !-----------------------------------------------------------------------------
  self%mesh%id = self%id
  self%mesh%b = self%box
  self%mesh%llc_cart = self%llc_cart
  !-----------------------------------------------------------------------------
  ! Cache a few things from the mesh
  !-----------------------------------------------------------------------------
  self%gn = self%mesh%gn
  self%li = self%mesh%li
  self%ui = self%mesh%ui
  self%ds = self%mesh%d
  self%ncell = self%mesh%nc
  self%centre_nat = self%mesh%centre_nat
  self%centre_cart = CurrentToCartesian(self%mesh_type, &
                     self%centre_nat(1), self%centre_nat(2), self%centre_nat(3))
  !-----------------------------------------------------------------------------
  ! Make dead sure these crucial values are correct for Cartesian meshes
  !-----------------------------------------------------------------------------
  self%llc_cart = self%position - 0.5_8*self%size
  if (self%mesh_type == mesh_types%Cartesian) then
    if (any(self%centre_nat /= self%position)) then
      call io%lock%set ('warn')
      if (warn) then
        call io%warning ('patch_t%centre_nat not equal to task_t%position')
        write (stderr,'(a,1p,3e16.6)') 'centre_nat =', self%centre_nat
        write (stderr,'(a,1p,3e16.6)') '  position =', self%position
        warn = .false.
      end if
      call io%lock%unset ('warn')
    end if
    self%llc_nat         = self%llc_cart
    self%mesh%llc_cart   = self%llc_cart
    self%mesh%llc_nat    = self%llc_cart
    self%centre_nat      = self%position
    self%centre_cart     = self%position
    self%mesh%p          = self%position
    self%mesh%centre_nat = self%position
  end if
  !-----------------------------------------------------------------------------
  ! Hacks for backwards compatibility
  !-----------------------------------------------------------------------------
  self%gsize = self%ds*self%gn
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
  ! Intialise staggering indices; it is the responsibility of the *solver* to
  ! assign non-zero values to this variable (if any).
  !-----------------------------------------------------------------------------
  do i=1,3
    allocate(self%mesh(i)%h(self%nv))
    self%mesh(i)%h = 0.0
    if (self%kind == 'zeus_mhd_patch' .and. &
        self%staggered .and. self%n(i) /= 1) self%mesh(i)%gn = self%mesh(i)%gn + 1
    if (self%kind == 'zeus_tw_mhd_patch' .and. &
        self%staggered .and. self%n(i) /= 1) self%mesh(i)%gn = self%mesh(i)%gn + 1
  end do
END SUBROUTINE mesh_init

!===============================================================================
!> Deprecated way to impose some key parameters from e.g. a test, coded before
!> the introduction of %pre_init(), and kept here only for the purpose of not
!> breaking those tests completely (FIXME: they should be rewritten in any case).
!===============================================================================
SUBROUTINE setup (self, position, size, n, ng, nt, nv, nw)
  class(patch_t)    :: self
  real(8), optional :: position(3), size(3)
  integer, optional :: n(3), ng(3)
  integer, optional :: nt, nv, nw
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%setup')
  if (present(position)) self%position = position
  if (present(size    )) self%size     = size
  if (present(n       )) self%n        = n
  if (present(ng      )) self%ng       = ng
  if (present(nt      )) self%nt       = nt
  if (present(nv      )) self%nv       = nv
  if (present(nw      )) self%nw       = nw
  call io%warning('WARNING: calling deprecated patch_t%setup!')
END SUBROUTINE setup

!===============================================================================
!> Make sure pointers are nullified (FIXME: make this redundant)
!===============================================================================
SUBROUTINE reset (self)
  class (patch_t):: self
  !.............................................................................
  self%mem_allocated = .false.                      ! make sure to get ..
  !--- task pointers ---
  nullify (self%lock)                               ! .. new lock
  nullify (self%t)                                  ! .. new t
  nullify (self%dt)                                 ! .. new dt
  nullify (self%iit)                                ! .. new iit
  nullify (self%mesg)                               ! .. new message
  nullify (self%sent_list)                          ! .. new message list
  !--- patch pointers ---
  nullify (self%mem)                                ! .. new mem
  nullify (self%mem_lock)                           ! .. new mem_lock
  nullify (self%mesh)                               ! .. new mesh
  nullify (self%density)                            ! gas density
  nullify (self%pressure)                           ! gas pressure
  nullify (self%vgas)                               ! gas velocity
  nullify (self%u_phase)                            ! phase speed
  nullify (self%frefine)                            ! refine condition
  nullify (self%request)                            ! refine support
  nullify (self%mask)                               ! used in extras
  nullify (self%filled)                             ! guard zone check
  nullify (self%outside)                            ! FIXME: needed??
  self%id = 0                                       ! ensure a new ID
  self%seq = -1                                     ! triggers tag=0 in send
  self%n_needed = 1                                 ! reset to initial
  self%istep = 0                                    ! for diagnostics
  self%rotated = .false.                            ! force rotate
END SUBROUTINE reset

!===============================================================================
!> Initialize the patch memory and patch dimensions.  For this we need to be
!> told the patch dimensions (which include the number of variable and the
!> number of time slices), the patch position, and the patch size
!===============================================================================
SUBROUTINE init (self)
  class(patch_t):: self
  !.............................................................................
  integer            :: i, gn(3)
  type(header_t)     :: header
  character(len=120) :: id = &
    '$Id$'
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%init')
  call io_unit%print_id (id)
  call self%task_t%init
  !-----------------------------------------------------------------------------
  ! Everything mesh-related -- call to be moved to %init()-chain later
  !-----------------------------------------------------------------------------
  call mesh_init (self)
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
  where (self%n==1.and.self%box==0)             ! for collapsed dimensions
    self%box = 1                                ! make sure it is non-zero
  end where
  !-----------------------------------------------------------------------------
  ! Region information, relative to the task position
  !-----------------------------------------------------------------------------
  self%authority%l = - self%size*0.5_8
  self%authority%u = + self%size*0.5_8
  self%interest%l  = self%authority%l - self%ng*self%ds
  self%interest%u  = self%authority%u + self%ng*self%ds
  !.............................................................................
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
    print 1,'patch_t%init: ds          =', self%ds
    print 1,'patch_t%init: authority%u =', self%authority%u
    print 1,'patch_t%init: interest%u  =', self%interest%u
    1 format(a,1p,3e15.6)
  end if
  !-----------------------------------------------------------------------------
  ! Recalculate level and hash key, in case size or position was changed
  !-----------------------------------------------------------------------------
  call self%init_level                           ! redo, in case of changes
  call self%task_t%pre_init ()                   ! redefine hash key
  !-----------------------------------------------------------------------------
  if (io%verbose>2) then
    !$omp critical (log_cr)
    write (io_unit%log,'(1x,a,i7,2(1p,3e12.4,2x))') 'patch_mod::init: id, size, pos =', &
      self%id, self%size, self%position
    if (io%verbose>3) &
      write (io_unit%log,*) 'id, LEVEL: ', self%id, self%level, maxval(self%box), minval(self%ds)
    !$omp end critical (log_cr)
  end if
  !-----------------------------------------------------------------------------
  ! Allocate memory, if not already associated
  !-----------------------------------------------------------------------------
  if (.not.self%mem_allocated) then
    self%mem_allocated = .true.
    !---------------------------------------------------------------------------
    ! extras%pre_init could increases self%mv = the actual number of memory slots
    !---------------------------------------------------------------------------
    if (self%id == 1 .and. self%mv > self%nv) then
      write (stderr,'(a)')    'NOTE:'
      write (stderr,'(a,i3)') '     Number of variables nv =', self%nv
      write (stderr,'(a,i3)') '     Number of mem slots mv =', self%mv
    end if
    self%mv = max (self%mv, self%nv)
    gn = self%mesh%gn
    allocate (self%mem    (gn(1), gn(2), gn(3), self%mv, self%nt, self%nw))
    allocate (self%u_phase(gn(1), gn(2), gn(3)))
    self%u_phase = 0.0
    allocate (self%mem_lock(self%nt))
    do i=1,self%nt
      call self%mem_lock(i)%init ('meml')
    end do
    call io%bits_mem (storage_size(self%mem), product(shape(self%mem)), 'mem')
    !---------------------------------------------------------------------------
    ! The unsigned and pervolume flags may also be modified by the solver. They
    ! are allocated here, but need to be reallocated by solvers modifying nv.
    ! Tracking even the small memory use, in case deallocation is forgotten.
    !---------------------------------------------------------------------------
    allocate (self%e(3,3), self%einv(3,3))
    self%e = 0.0
    do i=1,3
      self%e(i,i) = 1.0
      self%einv(i,i) = 1.0
    end do
    !---------------------------------------------------------------------------
    allocate (self%unsigned (self%mv))
    allocate (self%pervolume(self%mv))
    self%unsigned (:) = .false.
    self%pervolume(:) = .false.
    call io%bits_mem (storage_size(self%unsigned)  , &
                    product(shape(self%unsigned)) , 'uns')
    call io%bits_mem (storage_size(self%pervolume) , &
                    product(shape(self%pervolume)), 'perv')
    if (shared%refine_on) then
      allocate (self%request(gn(1),gn(2),gn(3)))
      call io%bits_mem (storage_size(self%request) , &
                       product(shape(self%request)), 'req')
      self%request = 0
    end if
  end if
  self%wallclock = wallclock()                  ! starting time
  n_header = storage_size(header)/32            ! 32 bits per word
  if (n_header*4 /= storage_size(header)/8) then
    write(io%output,*) n_header*4, storage_size(header)/8
    call mpi%abort ('n_header is incorrect')
  end if
  if (all(self%size == self%box)) then
    call self%set (bits%root_grid)
    if (io%verbose > 1) write(io%output,*) 'set bits%root_grid id =', self%id
  end if
  call self%set (bits%has_output)
  !-----------------------------------------------------------------------------
  ! Write out the variable indices to metadata
  !-----------------------------------------------------------------------------
  call io%lock%set ('patch_t%init')
  if (first_time) then
    write (stdout,'(a,6i4)') &
      ' memory shape:', shape(self%mem)
    call self%idx%print (stdout)
    !---------------------------------------------------------------------------
    ! If no other context has already set io%nv -- the number of memory slots
    ! to save in a snapshot, set it here to all memory slots
    !---------------------------------------------------------------------------
    if (io%nv==0) then
      io%nv = self%mv
    end if
    first_time = .false.
  end if
  call io%lock%unset ('patch_t%init')
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE init

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
  if (associated(self%mask)) then
    call io%bits_mem (-storage_size(self%mask), &
                      product(shape(self%mask)), '-mask')
    deallocate (self%mask)
  end if
  if (associated(self%frefine)) then
    call io%bits_mem (-storage_size(self%frefine), &
                      product(shape(self%frefine)), '-%frefine')
    deallocate (self%frefine)
  end if
  if (associated(self%unsigned)) then
    call io%bits_mem (-storage_size(self%unsigned), &
                      product(shape(self%unsigned)), '-unsigned')
    deallocate (self%unsigned)
  end if
  if (associated(self%pervolume)) then
    call io%bits_mem (-storage_size(self%pervolume), &
                      product(shape(self%pervolume)), '-pervolume')
    deallocate (self%pervolume)
  end if
  if (associated(self%mesg)) then
    call self%mesg%dealloc
    deallocate (self%mesg)
  end if
  call self%task_t%dealloc
  self%pre_initialized = .false.
  self%mem_allocated = .false.
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> Assiting refine_t%make_child_patch()
!===============================================================================
SUBROUTINE clone_mem_accounting (self)
  class(patch_t):: self
  if (associated(self%mask)) &
    call io%bits_mem (storage_size(self%mask), &
                     product(shape(self%mask)), 'mask')
  if (associated(self%frefine)) &
    call io%bits_mem (storage_size(self%frefine), &
                     product(shape(self%frefine)), '%frefine')
END SUBROUTINE clone_mem_accounting

!===============================================================================
!> This function gives the physical position of an index triple (patch point)
!> in absolute Cartesian coords.   Note that in some cases the periodic
!> wrapping must, as here, be done on a point-by-point basis.
!===============================================================================
FUNCTION myposition (self, ii, iv, pp) RESULT(out)
  class(patch_t)      :: self
  integer, intent(in) :: ii(3)
  integer, optional   :: iv
  real, optional      :: pp(3)
  real(8)             :: out(3)
  !.............................................................................
  real(8)             :: h(3)
  !-----------------------------------------------------------------------------
  if (present(pp)) then
    out = self%position + (real(ii)-self%offset + pp)*self%ds
  else
    out = self%position + (real(ii)-self%offset)*self%ds
  end if
  if (present(iv)) then
    h = self%index_stagger(iv)
    out = out + h*self%mesh%d
  end if
  !.............................................................................
  if (self%id==io%id_debug .and. io%verbose>1) &
    write(io%output,'(a,3i4,4(3x,3f10.2))') &
      'mypos:',ii,self%position,self%offset,self%ds,out
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
!> Return position of cell index, relative to the local patch
!===============================================================================
FUNCTION local_position (self, ii, iv) RESULT(out)
  class(patch_t)      :: self
  integer, intent(in) :: ii(3)
  integer, optional   :: iv
  real(8)             :: out(3)
  !.............................................................................
  real(8)             :: h(3)
  !-----------------------------------------------------------------------------
  out = (real(ii)-self%offset)*self%ds
  if (present(iv)) then
    h = self%index_stagger(iv)
    out = out + h*self%mesh%d
  end if
END FUNCTION local_position

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
  periodic = self%periodic .and. self%size+self%ds > self%box
  if (.not.any(periodic)) return
  call trace%begin ('patch_t%periodic_grid')
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
    !----------------------------------------------------------------------------
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
!> FIXME: Consider removal
!===============================================================================
SUBROUTINE extrapolate_guards (self)
  class(patch_t)        :: self
  !.............................................................................
  integer               :: it0, it1, it2, nt, iit(self%nt-1)
  integer, dimension(3) :: lb, lo, li, ub, uo, ui, i
  integer               :: iv
  integer, save         :: itimer=0
  real(8)               :: t(self%nt-1)
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%extrapolate_guards', 3, itimer=itimer)
  nt = self%nt
  call self%timeslots (iit, t)
  it0 = iit(nt-1)                                       ! current time slot
  it1 = iit(nt-2)                                       ! previous time slot
  it2 = iit(nt-3)                                       ! the one before that
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
  if (self%id==io%id_debug) then
    write(io%output,'(i5,2x,a,6i7)') self%id, &
      'MK extrapolating guard zone values', it2, it1, it0, self%it, self%new
  end if
  do iv=1,self%mv
    if (self%unsigned(iv)) then
      !-------------------------------------------------------------------------
      ! XY-layers (Z varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(3) > 1) then
            self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it1,1))*2. - &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it2,1)) )
            self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it1,1))*2. - &
        log(self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it2,1)) )
      end if
      !-------------------------------------------------------------------------
      ! ZX-layers (Y varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(2) > 1) then
            self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it1,1))*2. - &
        log(self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it2,1)) )
            self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it1,1))*2. - &
        log(self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it2,1)) )
      end if
      !-------------------------------------------------------------------------
      ! YZ-layers (X varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(1) > 1) then
            self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it1,1))*2. - &
        log(self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it2,1)) )
            self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it0,1) = exp(  &
        log(self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it1,1))*2. - &
        log(self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it2,1)) )
      end if
    else
      !-------------------------------------------------------------------------
      ! XY-layers (Z varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(3) > 1) then
        self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it1,1)*2. - &
        self%mem(lb(1):ub(1),lb(2):ub(2),lb(3):lo(3),iv,it2,1)
        self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it1,1)*2. - &
        self%mem(lb(1):ub(1),lb(2):ub(2),uo(3):ub(3),iv,it2,1)
      end if
      !-------------------------------------------------------------------------
      ! ZX-layers (Y varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(2) > 1) then
        self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it1,1)*2. - &
        self%mem(lb(1):ub(1),lb(2):lo(2),li(3):ui(3),iv,it2,1)
        self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it1,1)*2. - &
        self%mem(lb(1):ub(1),uo(2):ub(2),li(3):ui(3),iv,it2,1)
      end if
      !-------------------------------------------------------------------------
      ! YZ-layers (X varies the least)
      !-------------------------------------------------------------------------
      if (self%gn(1) > 1) then
        self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it1,1)*2. - &
        self%mem(lb(1):lo(1),li(2):ui(2),li(3):ui(3),iv,it2,1)
        self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it0,1) =   &
        self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it1,1)*2. - &
        self%mem(uo(1):ub(1),li(2):ui(2),li(3):ui(3),iv,it2,1)
      end if
    end if
  end do
  call trace%end (itimer)
END SUBROUTINE extrapolate_guards

!===============================================================================
!> Extrapolate (or interpolate) the given variable slot to the given time,
!> using Lagrange weigts, and storing the result in the %new time slot
!===============================================================================
SUBROUTINE time_extrapolate (self, iv, order, time)
  class(patch_t):: self
  integer:: iv, order
  real(8):: time
  !.............................................................................
  real(8):: w(order+1), t(self%nt-1)
  integer:: it(self%nt-1), i, i0, i2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%time_extrapolate', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Get patch times in increasing order, and their time slots
  !-----------------------------------------------------------------------------
  call self%timeslots (it, t)
  !-----------------------------------------------------------------------------
  ! Compute Lagrange interpolation weights for the sequence, with given order
  ! (needs to less than self%nt-1)
  !-----------------------------------------------------------------------------
  order = min (order, self%nt-2, self%istep)
  call lagrange%sequence_weights (time, t, i0, i2, w, order)
  !-----------------------------------------------------------------------------
  ! Sum up into slow %new, under OMP lock
  !-----------------------------------------------------------------------------
  call self%lock%set ('time_extrapolate')
  self%mem(:,:,:,iv,self%new,1) = &
  self%mem(:,:,:,iv,it(i0),1)*w(1)
  do i=i0+1,i2
    self%mem(:,:,:,iv,self%new,1) = &
    self%mem(:,:,:,iv,self%new,1) + &
    self%mem(:,:,:,iv,it(i),1)*w(i-i0+1)
  end do
  call self%lock%unset ('time_extrapolate')
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE time_extrapolate

!===============================================================================
!> Allocate a message structure, for communication
!===============================================================================
SUBROUTINE allocate_mesg (self)
  class(patch_t):: self
  class(mesg_t), pointer:: mesg
  integer:: n_buf
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%allocate_mesg')
  allocate (self%mesg)
  mesg => self%mesg
  mesg%id = self%id
  call mesg%init
  !-----------------------------------------------------------------------------
  ! Store the standard buffer size (not valid for swap requests!)
  !-----------------------------------------------------------------------------
  n_buf = product(self%gn)      ! 4-bytes; `anonymous_copy` also uses 4-bytes.
  if (kind(self%mem)==8) &
    n_buf = n_buf * 2           ! 8-bytes, but `anonymous_copy` uses 4-bytes!
  mesg%nbuf = n_header + self%nv*n_buf
  if (mpi_mesg%nbuf==0) then
    !$omp critical (nbuf_cr)
    if (mpi_mesg%nbuf==0) then
      mpi_mesg%nbuf = mesg%nbuf
      write (io_unit%log,*) &
        'patch_t%pack: setting mpi_mesg%nbuf =', mpi_mesg%nbuf
    end if
    !$omp end critical (nbuf_cr)
  end if
  !-----------------------------------------------------------------------------
  ! Add extra space in the buffer for swap requests, and allocate
  !-----------------------------------------------------------------------------
  if (self%is_set(bits%swap_request) .and. self%is_set(bits%boundary)) then
    mesg%nbuf = mesg%nbuf+ self%nv*n_buf
  end if
  allocate (mesg%buffer(mesg%nbuf))
  call io%bits_mem (storage_size(self%mesg%buffer), &
                    product(shape(self%mesg%buffer)))
  call trace%end ()
END SUBROUTINE allocate_mesg

!===============================================================================
!> Pack patch info into a buffer, for MPI communication
!> FIXME: Consider replacing with mpi_buffer() calls
!===============================================================================
SUBROUTINE pack (self)
  class(patch_t):: self
  !.............................................................................
  class(mesg_t), pointer:: mesg
  type(header_t):: header
  integer:: n_buf, ibuf, iv, it
  real:: pt
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%pack', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Allocate a task%mesg instance, mesg%reqs() and similar array
  !-----------------------------------------------------------------------------
  call self%allocate_mesg
  if (omp%nthreads > 1) then
    allocate (self%mesg%lock)
    call self%mesg%lock%init ('smsg')
  end if
  mesg => self%mesg
  allocate (mesg%reqs(mpi%size))
  allocate (mesg%fail(mpi%size))
  allocate (mesg%ranks(mpi%size))
  mesg%fail = 0
  mesg%reqs = 0
  mesg%ranks = 0
  call io%bits_mem (storage_size(mesg%reqs), &
                    product(shape(mesg%reqs)))
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
  !-----------------------------------------------------------------------------
  ! Packing meta-data is a task-owner read operation, and does not need locking
  !-----------------------------------------------------------------------------
  call patch_to_header (self, header)
  ibuf = 1
  call anonymous_copy (n_header, header, mesg%buffer(ibuf))
  ibuf = ibuf + n_header
  !-----------------------------------------------------------------------------
  ! Copy the variables to the output buffer.  If the bdry+swap bit are set,
  ! copy as a minimum also the previous time slot
  !-----------------------------------------------------------------------------
  n_buf = product(self%gn)
  if (kind(self%mem)==8) n_buf = n_buf*2
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
  !-----------------------------------------------------------------------------
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
!> Unpack patch info from an MPI buffer.  To ensure consistency, the task meta-
!> data should not be updated until after the time slot data are in place.
!> Any task that has this task as virtual nbor and is waiting for the time slot
!> will thus not be scheduled until both the meta-data and the time slot data
!> are up-to-date.  This also reduces the chance of lock contention with respect
!> to the new time slot.
!> FIXME: Consider replacing with mpi_buffer() calls
!===============================================================================
SUBROUTINE unpack (self, mesg)
  class(patch_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  type(header_t):: header
  type(lock_t), pointer:: lock
  integer:: n_buf, ibuf, iv, it, i
  integer, dimension(:), pointer:: buffer
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Compute size of buffer
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%unpack', itimer=itimer)
  call self%set (bits%busy)
  self%unpack_time = wallclock()
  buffer => mesg%buffer
  !--- 4-bytes per element; `anonymous_copy` also uses 4-bytes ---
  n_buf = product(self%gn)
  !--- 8-bytes per element, but `anonymous_copy` uses 4-bytes! ---
  if (kind(self%mem)==8) n_buf = n_buf * 2
  !-----------------------------------------------------------------------------
  ! Copy buffer into sequenced header, to get the new %it slot
  !-----------------------------------------------------------------------------
  ibuf = 1
  call anonymous_copy (n_header, buffer(ibuf), header)
  ibuf = ibuf + n_header
  call self%header_to_patch (header)
  !-----------------------------------------------------------------------------
  ! Copy array data from the buffer to the time slot, under time slot lock
  !-----------------------------------------------------------------------------
  lock => self%mem_lock(header%it)
                                                        call lock%set ('unpack')
  do iv=1,self%nv
    call anonymous_copy (n_buf, buffer(ibuf), self%mem(:,:,:,iv,header%it,1))
    ibuf = ibuf + n_buf
  end do
                                                      call lock%unset ('unpack')
  !-----------------------------------------------------------------------------
  ! Update the patch meta-data under lock, now that the time slices are in place
  !-----------------------------------------------------------------------------
                                                   call self%lock%set ('unpack')
  self%it         = header%it
  self%time       = header%time
  self%nt         = header%nt
  if (header%nt > 0) then
    self%new      = mod(header%it,header%nt) + 1
    do i=1,self%nt
      self%t(i)   = header%t(i)
      self%dt(i)  = header%dt(i)
      self%iit(i) = header%iit(i)
    end do
  else
    write (io_unit%log,*) omp%thread, &
      'header_to_patch WARNING: nt non-positive', header%nt
    flush (io_unit%log)
  end if
!.............................................................................
                                               call self%lock%unset ('unpack')
  call self%clear (bits%busy)
  self%wc_last = wallclock()
  !-----------------------------------------------------------------------------
  ! Diagnostics
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (io_unit%log,'(f12.6,3x,a,i7,2i4,2x,3l1,2i4,f8.4,f11.6)') &
      wallclock(), 'patch_t%unpack: id, nq, BVS, nsnd, nrcv, latency, time =', &
      self%id, mesg%sender, self%nq, &
      self%is_set(bits%boundary), &
      self%is_set(bits%virtual), &
      self%is_set(bits%swap_request), &
      mpi_mesg%sent_list%n, mpi_mesg%recv_list%n, self%latency, self%time
    write (io_unit%log,1) &
      'patch_t%unpack id, it, nq, time  after:', &
      self%id, self%it, self%nq, self%time, &
      self%is_set (bits%internal), &
      self%is_set (bits%boundary), &
      self%is_set (bits%virtual), &
      self%is_set (bits%external), &
      self%is_set (bits%swap_request), &
      self%wc_last
    1 format(a,3i6,f12.6,2x,5l1,f12.1)
    if (io%verbose > 2) then
      write (io_unit%log,*) &
        self%id, 'after unpack', minval(self%mem(:,:,:,1,self%it,1))
    end if
  end if
  !-----------------------------------------------------------------------------
    call trace%end (itimer)
END SUBROUTINE unpack

!===============================================================================
!> Signed distance from task%position to self%position, wrapped in the case of a
!> periodic box.
!===============================================================================
FUNCTION distance (self, task) RESULT (out)
  class(patch_t):: self
  class(task_t):: task
  real(8):: out(3)
  integer, save:: itimer=0
  !.............................................................................
  if (timer%detailed) &
    call trace%begin ('patch_t%distance')
  if (self%mesh_type == mesh_types%Cartesian) then
    out = self%position-task%position
  else
    out = distance_curvilinear (self, task)
  end if
  !--- account for periodic wrapping ---
  where (self%periodic .and. self%box > 0d0)
    out = modulo (out+0.5_8*self%box, self%box) - 0.5_8*self%box
  end where
  if (timer%detailed) &
    call trace%end (itimer)
END FUNCTION distance

!===============================================================================
!> Signed distance between the centers of two task, adjusted to account for
!> curvilinear coordnates (the position of a curvilinear patch is set to the
!> coordinate *origin) and requires components defined in `patch_t`.
!===============================================================================
FUNCTION distance_curvilinear (self, task) RESULT (out)
  class(patch_t):: self
  class(task_t):: task
  real(8):: out(3)
  integer, save:: itimer=0
  !.............................................................................
  if (timer%detailed) &
    call trace%begin ('patch_t%distance_curvilinear', itimer=itimer)
  select type (task)
  class is (patch_t)
    out = self%centre_nat-task%centre_nat
  class default
    out = self%position-task%position
  end select
  ! account for periodic wrapping
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
  if (timer%detailed) &
    call trace%end (itimer)
END FUNCTION distance_curvilinear

!===============================================================================
!> The nearest smaller index position in the self patch coordinate system.  The
!> point with index self%mesh%lb is half a cell inside the patch guard size, so the
!> floating point position is counted relative to that point, and is split into
!> an integer part and a remainder.  For patches smaller than the box, the index
!> may become both smaller than self%mesh%lb and larger than self%mesh%ub.  For
!> the root grid patch, the index is in the range (li-1,ui) inside the interior.
!===============================================================================
FUNCTION index (self, other, oindex, f)
  class(patch_t):: self, other
  real, dimension(3):: position, f, pos, pos_lb
  integer, dimension(3):: index, oindex
  !.............................................................................
  pos = other%position - self%position                          ! relative pos
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  pos = pos + (oindex-other%offset)*other%ds                    ! position in other
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
  pos = other%position - self%position                          ! relative pos
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  pos = pos + (index-other%offset)*other%ds                     ! position in other
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
  pos = other%position - self%position                          ! relative pos
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  pos = pos + (index-other%offset-0.49)*other%ds     ! position in other
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
  pos = other%position - self%position                          ! relative pos
  pos = modulo (pos+0.5*self%box,self%box) - 0.5*self%box       ! wrap
  pos = pos + (index-other%offset+0.49)*other%ds                ! position in other
  pos = pos-self%position                                       ! relative
  f = self%offset + pos/self%ds                                 ! floating integer
  out = nint(f)                                                 ! highest integer
END FUNCTION highest

!===============================================================================
!> Courant condition, one way or another
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(patch_t):: self
  logical, optional:: detailed_timer
  real(8):: dtime, ds(3), courant, dsmin
  integer:: i2
  logical, parameter:: save_dbg=.false.
  !.............................................................................
  call trace%begin('patch_t%courant_condition')
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
  if (all(self%n==1)) then
    dsmin = minval(ds)                  ! single cell tests
  else
    dsmin = minval(ds,self%n>1)         ! ignore n==1 dimensions
  end if
  courant = self%courant
  if (self%dt_fixed > 0.0) then
    dtime = self%dt_fixed
    self%u_max = self%courant*dsmin/dtime
  else if (self%u_fixed > 0.0) then
    dtime = self%courant*dsmin/self%u_fixed
  else
    if (courant_ramp) then
      courant = self%courant*min(self%istep+1,self%nt)/real(self%nt)
    end if
    dtime = courant*dsmin/self%u_max*self%safety_factor
  endif
  !-----------------------------------------------------------------------------
  !> Prevent too rapid increase of time steps
  !-----------------------------------------------------------------------------
  if (self%nt > 1 .and. self%istep > 5) then
    if (dtime < 0.5*self%dt(self%iit(self%nt-1)) .and. self%dt_fixed < 0.0) then
      self%safety_factor = 0.5
    end if
    self%safety_factor = self%safety_factor**0.9
    dtime = dtime*self%safety_factor
  end if
  !$omp atomic write
  self%dtime = dtime
  !-----------------------------------------------------------------------------
  ! Force the time step to give exactly sync_time if within reach (active only
  ! if/when self%sync_time is being updated).  Syncing to a common time makes
  ! creates a point in time where all tasks are "ready to update".  This may
  ! improve performance eg when the number of tasks per core is relatively low.
  !-----------------------------------------------------------------------------
  if (self%time < self%sync_time .and. self%time+self%dtime > self%sync_time) then
    self%dtime = self%sync_time - self%time
    self%syncing = .true.
    if (io%verbose>2) &
      write (io_unit%mpi,*) 'patch_t%courant_condition: syncing turned on for', &
        self%id, self%time
  end if
  !-----------------------------------------------------------------------------
  ! Apply constraints that should hold for tasks of any kind, and set %dt(it)
  !-----------------------------------------------------------------------------
  call self%task_t%courant_condition (detailed_timer)
  self%dt(self%it) = self%dtime
  !-----------------------------------------------------------------------------
  if (io%verbose> 3 .or. self%id == io%id_debug) &
    write(io%output,'(a,i7,f12.6,f8.4,1p,e10.2)') 'courant:', &
      self%id, self%dtime, courant, self%u_max
  call trace%end ()
END SUBROUTINE courant_condition

!===============================================================================
!> Compute nbor%u_max -- the maximum phase speed in the guard zones of the nbor
!> that reside inside the current patch.
!===============================================================================
SUBROUTINE nbor_courant_conditions (self)
  class(patch_t):: self
  !.............................................................................
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: nbtask
  integer:: i1, i2, i3, l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%nbor_courant_conditions', itimer=itimer)
  !-----------------------------------------------------------------------------
  nbor => self%link%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (patch_t)
      call overlap_range (self, nbtask, l, u)
      nbor%u_max = maxval(self%u_phase(l(1):u(1),l(2):u(2),l(3):u(3)))
      if (self%id==io%id_debug) then
        write(io_unit%log,'(1x,a,2(i6,f12.6),2(2x,3i4),1p,2e12.2)') &
         'overlap_range:', self%id, self%time, nbtask%id, nbtask%time, &
          l, u, nbor%u_max, maxval(self%u_phase)
      end if
    end select
    nbor => nbor%next
  end do
  call trace%end (itimer)
END SUBROUTINE nbor_courant_conditions

!===============================================================================
!> Find the self index ranges that correspond to the overlap region, for
!> self with respect to source.
!> FIXME: Consider replacing by region_t operations
!===============================================================================
SUBROUTINE overlap_range (self, source, l, u)
  class(patch_t)          :: self, source
  integer                 :: l(3), u(3)
  !.............................................................................
  real(8), dimension(3)   :: dist, la, ua
  integer                 :: i
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%overlap_range', 2)
  dist = source%position - self%position
  where (source%periodic)
    dist = modulo(dist+0.5_8*source%mesh%b,source%mesh%b) - 0.5_8*source%mesh%b
  end where
  la = dist - source%gsize/2d0
  ua = dist + source%gsize/2d0
  l = ceiling(la/self%mesh%d + self%mesh%o)
  u =   floor(ua/self%mesh%d + self%mesh%o)
  l = max(l,self%mesh%li)
  u = min(u,self%mesh%ui)
  where (self%n == 1)
    l = min(l,self%mesh%ui)
    u = max(u,self%mesh%li)
  end where
  call trace%end()
END SUBROUTINE overlap_range

!===============================================================================
!> Overloaded sum over internal cells only
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
!> Overloaded average over internal cells only
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
FUNCTION faveri (self, iv) RESULT (out)
  class(patch_t):: self
  real(8):: out
  integer:: iv, l(3), u(3)
  !-----------------------------------------------------------------------------
  l = self%li
  u = self%ui
  out = sum(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,self%it,1)) &
      / product(u-l+1)
END FUNCTION faveri

!===============================================================================
!> Overloaded min/max over internal cells only
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
  integer:: l(3), u(3)
  !-----------------------------------------------------------------------------
  if (present(outer)) then
    outer_l = outer
  else
    outer_l = .false.
  end if
  l = self%mesh%li
  u = self%mesh%ui
  if (outer_l) then
    where(self%n==1)
      l = self%mesh%li
      u = self%mesh%ui
    end where
  end if
  fmaxval = maxval(f(l(1):u(1),l(2):u(2),l(3):u(3)))
END FUNCTION fmaxval4
!===============================================================================
FUNCTION fminval4 (self, f) RESULT(fminval)
  real:: fminval
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
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
  if (outer_l) then
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
  fminval = minval(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui))
END FUNCTION fminval8

!===============================================================================
FUNCTION fminloc4 (self, f) RESULT(fminloc)
  integer, dimension(3):: fminloc
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  fminloc = minloc(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui)) + self%ng
END FUNCTION fminloc4
!===============================================================================
FUNCTION fmaxloc4 (self, f) RESULT(fmaxloc)
  integer, dimension(3):: fmaxloc
  class(patch_t):: self
  real, dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  fmaxloc = maxloc(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui)) + self%ng
END FUNCTION fmaxloc4
!===============================================================================
FUNCTION fminloc8 (self, f) RESULT(fminloc)
  integer, dimension(3):: fminloc
  class(patch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  fminloc = minloc(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui)) + self%ng
END FUNCTION fminloc8
!===============================================================================
FUNCTION fmaxloc8 (self, f) RESULT(fmaxloc)
  integer, dimension(3):: fmaxloc
  class(patch_t):: self
  real(8), dimension(:,:,:), intent(in):: f
  !-----------------------------------------------------------------------------
  fmaxloc = maxloc(f(self%mesh(1)%li:self%mesh(1)%ui, &
                     self%mesh(2)%li:self%mesh(2)%ui, &
                     self%mesh(3)%li:self%mesh(3)%ui)) + self%ng
END FUNCTION fmaxloc8

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
    if (isinf(self%mem(ix,iy,iz,iv,self%it,1)) .or. &
        isnan(self%mem(ix,iy,iz,iv,self%it,1))) then
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
  ! If I/O reducing measures are in effect, or self%mem is not associated, or
  ! the task is frozen, skip the test
  !-----------------------------------------------------------------------------
  if ((io%verbose < 2 .or. io%print_time > 0.0 .or. io%print_every > 1) .and. &
      .not.do_check_nan) return
  if (self%is_set(bits%no_density+bits%frozen)) return
  if (.not. associated(self%mem)) return
  !-----------------------------------------------------------------------------
  call timer%begin ('patch_t%check_density', itimer=itimer)
  l = self%mesh%lb
  u = self%mesh%ub
  mind = minval(self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d,self%it,1))
  if (present(label) .and. io%verbose>1) write (io_unit%log,*) label, self%id, mind
  panic = .false.
  if (mind <= 0.0) then
    panic = .true.
    call self%dump (2, .true.)
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
SUBROUTINE stats_patch (self, label, unit, verbose)
  class (patch_t):: self
  character(len=*):: label
  integer, optional:: unit, verbose
  !.............................................................................
  integer:: iv, unit_l
  real(8), dimension(4):: si, sb
  class(mesh_t), pointer:: mesh(:)
  integer, save:: tag=0
  !-----------------------------------------------------------------------------
  if (present(unit)) then
    unit_l = unit
  else
    unit_l = io%output
  end if
  mesh => self%mesh
  do iv=1,self%mv
    call scalar_stats (self%mem(:,:,:,iv,self%it,1), mesh%li, mesh%ui, si)
    call scalar_stats (self%mem(:,:,:,iv,self%it,1), mesh%lb, mesh%ub, sb)
    write (unit_l, '(i7,2x,a,i3,1p,2(3x,4g12.4))') &
      self%id, trim(label)//' stats: iv, min, aver, rms, max =', iv, si, sb
  end do
  !.............................................................................
  if (self%id == io%id_debug) then
    tag = tag+1
    if (present(verbose)) then
      if (verbose > 0) &
        print *, 'dump: id,iv,tag=', self%id,   tag 
    end if
    call self%dump (tag, it=self%it, ad_hoc=.true.)
  end if
  flush (unit_l)
END SUBROUTINE stats_patch

!===============================================================================
!> Statistics for a scalar field
!===============================================================================
SUBROUTINE stats_scalar (self, var, label, unit)
  class (patch_t):: self
  real(kind=KindScalarVar), dimension(:,:,:):: var
  character(len=*):: label
  integer, optional:: unit
  !.............................................................................
  integer:: unit_l
  real(8), dimension(4):: si, sb
  !-----------------------------------------------------------------------------
  if (present(unit)) then
    unit_l = unit
  else
    unit_l = io%output
  end if
  call scalar_stats (var, self%mesh%li, self%mesh%ui, si)
  call scalar_stats (var, self%mesh%lb, self%mesh%ub, sb)
  write (unit_l, '(i7,2x,a,1p,2(3x,4g12.4),2x,a)') &
    self%id, ' stats: min, aver, rms, max =', si, sb,label
  flush (unit_l)
END SUBROUTINE stats_scalar

!===============================================================================
!> Statistics for a scalar field
!===============================================================================
SUBROUTINE stats_scalar8 (self, var, label, unit)
  class (patch_t):: self
  real(8), dimension(:,:,:):: var
  character(len=*):: label
  integer, optional:: unit
  !.............................................................................
  integer:: unit_l
  real(8), dimension(4):: si, sb
  !-----------------------------------------------------------------------------
  if (present(unit)) then
    unit_l = unit
  else
    unit_l = io%output
  end if
  call scalar_stats (real(var), self%mesh%li, self%mesh%ui, si)
  call scalar_stats (real(var), self%mesh%lb, self%mesh%ub, sb)
  write (unit_l, '(i7,2x,a,1p,2(3x,4g12.4),2x,a)') &
    self%id, ' stats: min, aver, rms, max =', si, sb,label
  flush (unit_l)
END SUBROUTINE stats_scalar8

!===============================================================================
!> Statistics for a vector field
!===============================================================================
SUBROUTINE stats_vector (self, var, label, unit)
  class (patch_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:):: var
  character(len=*):: label
  integer, optional:: unit
  integer:: unit_l
  !.............................................................................
  integer:: iv
  real(8), dimension(4):: si, sb
  !-----------------------------------------------------------------------------
  if (present(unit)) then
    unit_l = unit
  else
    unit_l = io%output
  end if
  do iv=1,size(var,4)
    call scalar_stats (var(:,:,:,iv), self%mesh%li, self%mesh%ui, si)
    call scalar_stats (var(:,:,:,iv), self%mesh%lb, self%mesh%ub, sb)
    write (io_unit%output,'(i7,2x,a,i3,1p,2(3x,4g12.4))') &
      self%id, trim(label)//' stats: iv, min, aver, rms, max =', iv, si, sb
  end do
  flush (unit_l)
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

!===============================================================================
!> Interpolate -- only in time for now -- from the MHD time slots save in
!> this module to the time requested from the solve procedure
!===============================================================================
SUBROUTINE log_interpolate (self, time, i, a)
  class (patch_t)                 :: self
  real(8)                         :: time
  integer                         :: i
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: a, b
  !.............................................................................
  integer                         :: jt(2)
  real                            :: pt(2)
  integer, save                   :: itimer=0
  !-----------------------------------------------------------------------------
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

!===============================================================================
!> Search for a time interval to interpolate in, and return the time
!> slot indices and weight factors
!===============================================================================
SUBROUTINE time_interval (self, time, jt, pt, all, id, extrap)
  class (patch_t)                 :: self
  real(8)                         :: time
  integer                         :: jt(2)
  real                            :: pt(2)
  logical, optional               :: all
  integer, optional               :: id
  logical, optional               :: extrap
  !.............................................................................
  integer                         :: i, j, n, iit(self%nt-1)
  real(8)                         :: t(self%nt-1), dt
  real                            :: grace
  integer, save                   :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('patch_t%time_interval', 2, itimer=itimer)
  grace = self%grace + shared%extra_grace
  if (self%nt==1) then
    jt = 1
    pt(2) = 0.0
  else
    call self%timeslots (iit, t)
    if (io%verbose > 0) then
      !-------------------------------------------------------------------------
      ! Self test
      !-------------------------------------------------------------------------
      call omp_lock%check_pointer (self%lock, 'time_interval', self%id)
      call self%timeslots (iit, t)
      if (self%verbose > 2) then
        do i=1,self%nt-1
          if (iit(i) /= self%iit(i)) then
            write (stderr,*) 'patch_t%time_interval ERROR: i, iit, self%iit =', &
            i, iit(i), self%iit(i)
          end if
          if (t(i) /= self%t(self%iit(i))) then
            write (stderr,*) 'patch_t%time_interval ERROR: i, times(i), iit(i), t', &
              i, t(i), self%iit(i), self%t(self%iit(i))
          end if
        end do
      end if
    end if
    j = 1
    n = self%nt-2
    !---------------------------------------------------------------------------
    ! Seach and set the j index to the last point that has t(j)<time.
    ! Note that this does NOT guarantee that t(j+1) >= time.  In cases
    ! where both t(nt-2) < time and t(n-1) < time, extrapolation in time
    ! will result.  What should NOT happen is that a memory slot that
    ! has not yet received an update is chosen.
    !---------------------------------------------------------------------------
    do i=1,n
      if (t(i) < time) j=i
    end do
    jt = iit(j:j+1)
    dt = t(j+1)-t(j)
    if (dt==0.0) then
      pt(2) = 0.0
    else
      pt(2) = (time-t(j))/dt
      if (pt(2) > 1.01d0 + grace .and. .not.present(extrap)) then
        pt(1) = 1.-pt(2)
        if (present(id)) then
          write (stderr,'(2a,2i6,":",i2.2,f12.6,f7.2,9f12.6)') &
            'patch_t%time_interval: excessive extrapolation prevented, ', &
            'id,src,time,pt,t=', id, self%id, self%level, time, pt(2), t
          write (io_unit%mpi,'(f8.3,i4,2x,a,2i6,1p,4e11.2)') &
            wallclock(), omp%thread, &
            'patch_t%time_interval: excessive id,id,t,tm,t,f =', &
            self%id, id, time, time-t(j), t(j+1)-t(j), (time-t(j))/dt
          write (io_unit%mpi,'(f8.3,i4,2x,a,2i2,2f5.2,1p,11e11.3)') &
            wallclock(), omp%thread, &
            'patch_t%time_interval: jt,pt,tm,t(:) =', &
            jt, pt, time, self%t-time
        else
          write (stderr,'(2a,i6,":",i2.2,f12.6,f7.2,9f12.6)') &
            'patch_t%time_interval: excessive extrapolation prevented, ', &
            'id,time,pt,t=', self%id, self%level, time, pt(2), t
          write (io_unit%mpi,'(i4,2x,a,i6,1p,6e11.2)') omp%thread, &
            'patch_t%time_interval: excessive id,t,tm,t,f =', &
            self%id, time, t(j+1), time-t(j), t(j+1)-t(j), (time-t(j))/dt
        end if
        if (io%verbose > 1) then
          write (stdlog,'(a,2i6,2f6.2      )') 'time_interval:', jt, pt
          write (stdlog,'(a,2i6            )') '   nbor, targ:', self%id, id
          write (stdlog,'(a,9i4            )') '     self%iit:', self%iit
          write (stdlog,'(a,9i4            )') '          iit:', iit
          write (stdlog,'(a,f12.6,2x,9f12.6)') '  time, times:', time, self%time
          write (stdlog,'(a,f12.6,2x,9f12.6)') '     grace, t:', grace, t
          flush (stdlog)
        end if
        if (io%task_logging > 2) then
          block
          character(len=24):: label
          write (label,'(a,i6,f10.5)') 'excess ', id, time
          call self%log (label)
          end block
        end if
      end if
      pt(2) = max(0d0-grace, min(1d0+grace, pt(2)))
    end if
  end if
  pt(1) = 1.-pt(2)
  if (io%verbose > 3) then
    if (io%verbose > 4 .or. (io_unit%do_validate .and. pt(2)>0.0)) then
      write(io%output,'(a,i6,i6,1p,e14.6,i3,0p,2f8.4,2x,1p,7e14.6)') &
        'time_interval: times  =', self%id, self%istep, time, j, pt, t
      write(io%output,'(a,9i10  )') 'time_interval: iit =', iit
      write(io%output,'(a,9f10.4)') 'time_interval:   t =', t
    end if
  end if
  call trace%end (itimer)
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
FUNCTION index_only (self, pos, iv, roundup, rounddn, nearest, nowrap, local) &
    RESULT (ii)
  class(patch_t)      :: self
  real(8), intent(in) :: pos(3)
  integer, optional   :: iv
  logical, optional   :: roundup, rounddn, nearest, nowrap, local
  integer             :: ii(3)
  !.............................................................................
  integer             :: i
  real(8)             :: p(3), h(3)
  class(mesh_t), pointer :: mesh
  !-----------------------------------------------------------------------------
  if (present(iv)) &
    h = self%index_stagger(iv)
  do i=1,3
    if (self%n(i) == 1) then
      ii(i) = self%li(i)
      cycle
    end if
    mesh => self%mesh(i)
    if (.not.present(local)) then
      p(i) = pos(i) - self%position(i)                          ! relative position
    else
      p(i) = pos(i)
    end if
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
    else if (present(rounddn)) then
      ii(i) = floor(p(i))                                       ! round dn
    else if (present(nearest)) then
      ii(i) = nint(p (i))                                       ! round nearest
    else
      ii(i) = floor(p(i) + 0.0001)                              ! floor w margin
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
  class(patch_t)           :: target
  integer                  :: ii(3), iv, jt(2)
  real                     :: pt(2)
  !.............................................................................
  real(8)                  :: pos(3)
  integer                  :: j(3), i1, j1, k1
  real                     :: p(3), q(3)
  integer, save            :: itimer=0
  !-----------------------------------------------------------------------------
  if (self%pervolume(iv)) then
    call self%interpolate_specific (target, ii, iv, jt, pt)
  else
    i1 = merge(1,0,size(self%mem,1) > 1)
    j1 = merge(1,0,size(self%mem,2) > 1)
    k1 = merge(1,0,size(self%mem,3) > 1)
    pos = target%myposition (ii, iv)
    call self%index_space (pos, iv, j, p)
    if (any(j < self%li .or. j > self%ui-[i1,j1,k1])) then
      write (io_unit%mpi,'(1x,a,3i4,2x,3f6.3)') &
        'patch_t%interpolate ERROR: j, p =', j, p
      where (j < self%li)
        p = p + (j-self%li)
        j = self%li
      end where
      where (j > self%ui-[i1,j1,k1])
        p = p + (j-(self%ui-[i1,j1,k1]))
        j = (self%ui-[i1,j1,k1])
      end where
      p = max(0.0,min(1.0,p))
    end if
    q = 1.0-p
    if (self%unsigned(iv)) then
      if (pt(2) == 0.0) then
        target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = exp( &
          q(3)*(q(2)*(q(1)*log(self%mem(j(1)   ,j(2)   ,j(3)   ,iv,jt(1),1))    + &
                      p(1)*log(self%mem(j(1)+i1,j(2)   ,j(3)   ,iv,jt(1),1)))   + &
                p(2)*(q(1)*log(self%mem(j(1)   ,j(2)+j1,j(3)   ,iv,jt(1),1))    + &
                      p(1)*log(self%mem(j(1)+i1,j(2)+j1,j(3)   ,iv,jt(1),1))))  + &
          p(3)*(q(2)*(q(1)*log(self%mem(j(1)   ,j(2)   ,j(3)+k1,iv,jt(1),1))    + &
                      p(1)*log(self%mem(j(1)+i1,j(2)   ,j(3)+k1,iv,jt(1),1)))   + &
                p(2)*(q(1)*log(self%mem(j(1)   ,j(2)+j1,j(3)+k1,iv,jt(1),1))    + &
                      p(1)*log(self%mem(j(1)+i1,j(2)+j1,j(3)+k1,iv,jt(1),1)))))
      else
        target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = exp( &
         pt(1)*(q(3)*(q(2)*(q(1)*log(self%mem(j(1)   ,j(2)   ,j(3)   ,iv,jt(1),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)   ,j(3)   ,iv,jt(1),1)))   + &
                      p(2)*(q(1)*log(self%mem(j(1)   ,j(2)+j1,j(3)   ,iv,jt(1),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)+j1,j(3)   ,iv,jt(1),1))))  + &
                p(3)*(q(2)*(q(1)*log(self%mem(j(1)   ,j(2)   ,j(3)+k1,iv,jt(1),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)   ,j(3)+k1,iv,jt(1),1)))   + &
                      p(2)*(q(1)*log(self%mem(j(1)   ,j(2)+j1,j(3)+k1,iv,jt(1),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)+j1,j(3)+k1,iv,jt(1),1))))) + &
         pt(2)*(q(3)*(q(2)*(q(1)*log(self%mem(j(1)   ,j(2)   ,j(3)   ,iv,jt(2),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)   ,j(3)   ,iv,jt(2),1)))   + &
                      p(2)*(q(1)*log(self%mem(j(1)   ,j(2)+j1,j(3)   ,iv,jt(2),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)+j1,j(3)   ,iv,jt(2),1))))  + &
                p(3)*(q(2)*(q(1)*log(self%mem(j(1)   ,j(2)   ,j(3)+k1,iv,jt(2),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)   ,j(3)+k1,iv,jt(2),1)))   + &
                      p(2)*(q(1)*log(self%mem(j(1)   ,j(2)+j1,j(3)+k1,iv,jt(2),1))    + &
                            p(1)*log(self%mem(j(1)+i1,j(2)+j1,j(3)+k1,iv,jt(2),1))))))
      end if
    else
      if (pt(2) == 0.0) then
        target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
          (1.-p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)   ,j(3)   ,iv,jt(1),1)    + &
                                (    p(1))*self%mem(j(1)+i1,j(2)   ,j(3)   ,iv,jt(1),1))   + &
                     (   p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)+j1,j(3)   ,iv,jt(1),1)    + &
                                (    p(1))*self%mem(j(1)+i1,j(2)+j1,j(3)   ,iv,jt(1),1)))  + &
          (   p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)   ,j(3)+k1,iv,jt(1),1)    + &
                                (    p(1))*self%mem(j(1)+i1,j(2)   ,j(3)+k1,iv,jt(1),1))   + &
                     (   p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)+j1,j(3)+k1,iv,jt(1),1)    + &
                                (    p(1))*self%mem(j(1)+i1,j(2)+j1,j(3)+k1,iv,jt(1),1)))
      else
        target%mem(ii(1),ii(2),ii(3),iv,target%it,1) = &
         pt(1)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)   ,j(3)   ,iv,jt(1),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)   ,j(3)   ,iv,jt(1),1))   + &
                           (   p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)+j1,j(3)   ,iv,jt(1),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)+j1,j(3)   ,iv,jt(1),1)))  + &
                (   p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)   ,j(3)+k1,iv,jt(1),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)   ,j(3)+k1,iv,jt(1),1))   + &
                           (   p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)+j1,j(3)+k1,iv,jt(1),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)+j1,j(3)+k1,iv,jt(1),1)))) + &
         pt(2)*((1.-p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)   ,j(3)   ,iv,jt(2),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)   ,j(3)   ,iv,jt(2),1))   + &
                           (   p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)+j1,j(3)   ,iv,jt(2),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)+j1,j(3)   ,iv,jt(2),1)))  + &
                (   p(3))*((1.-p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)   ,j(3)+k1,iv,jt(2),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)   ,j(3)+k1,iv,jt(2),1))   + &
                           (   p(2))*((1.0-p(1))*self%mem(j(1)   ,j(2)+j1,j(3)+k1,iv,jt(2),1)    + &
                                      (    p(1))*self%mem(j(1)+i1,j(2)+j1,j(3)+k1,iv,jt(2),1))))
      end if
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
  class(patch_t)           :: target
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
  ! Multiply by target density
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
!> Check overlap and return overlap range in a patch, in terms interior region
!> self indices.  In other words:  Get the self (source) index range that the
!> interior of patch (target) covers.  If guards==.true., then include also the
!> guard zones of the target patch.
!===============================================================================
FUNCTION get_overlap (self, patch, guards, l, u) RESULT(overlap)
  class (patch_t):: self
  class (patch_t):: patch
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
      la = dist - patch%size/2d0 - (patch%mesh%ng-0.5d0)*patch%mesh%d
      ua = dist + patch%size/2d0 + (patch%mesh%ng-0.5d0)*patch%mesh%d
    else
      la = dist - patch%size/2d0
      ua = dist + patch%size/2d0
    end if
    !---------------------------------------------------------------------------
    ! This works for both NML and NNML meshes:
    !---------------------------------------------------------------------------
    l =   floor(la/self%mesh%d + self%mesh%o)
    u = ceiling(ua/self%mesh%d + self%mesh%o)
    l = min(max(l,self%mesh%li),self%mesh%ub)
    u = min(max(u,self%mesh%lb),self%mesh%ui)
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
  if (any(self%n > 1)) then
    self%level = maxval(nint(log(self%box/self%ds) / &
      log(real(shared%refine_ratio))), mask=self%n > 1)
  end if
END SUBROUTINE init_level

!===============================================================================
!> Dump the patch values into data/run/rank_rrrrr.dump
!===============================================================================
SUBROUTINE dumpm (self, tag, ad_hoc, it)
  class(patch_t):: self
  integer:: tag
  logical, optional:: ad_hoc
  integer, optional:: it
  !.............................................................................
  integer:: nt, jt
  !-----------------------------------------------------------------------------
  if (self%dump_it .or. present(ad_hoc)) then
    call io%lock%set ('dump')
    nt = merge (self%nt, 1, present(ad_hoc) .and. .not.present(it))
    if (self%verbose > 0) then
      write (io_unit%log,*) 'patch_t%dump: id, time =', self%id, self%time
    end if
    write (io_unit%dump) 0, kind(self%mem), tag, self%id, self%level, self%mv, &
                         self%nt, nt, self%nw, self%istep, self%it, self%gn, &
                         self%n
    write (io_unit%dump) self%time, self%dtime, self%position, self%size, &
                         self%t, self%dt, self%iit
    write (io_unit%dump) kind(self%mem), size(shape(self%mem(:,:,:,:,1,1)))
    write (io_unit%dump) shape(self%mem(:,:,:,:,1,1))
    !---------------------------------------------------------------------------
    ! Write nt time slots, in reverse order, starting from the most recent one
    !---------------------------------------------------------------------------
    if (present(it)) then
      write (io_unit%dump) self%mem(:,:,:,:,it,1)
    else
      do jt=1,nt
        write (io_unit%dump) self%mem(:,:,:,:,self%iit(self%nt-jt+1),1)
      end do
    end if
    call io%lock%unset ('dump')
  end if
END SUBROUTINE dumpm

!===============================================================================
!> Dump an array into data/run/rank_rrrrr.dump
!===============================================================================
SUBROUTINE dump1r (self, tag, a)
  class(patch_t):: self
  integer:: tag
  real:: a(:)
  !.............................................................................
  integer:: nt, nv, gn(3), n(3)
  !-----------------------------------------------------------------------------
  call io%lock%set ('dump')
  nt = 1
  nv = 1
  n  = 1
  gn = 1
  n(1) = size(a)
  gn(1) = size(a)
  if (self%verbose > 0) then
    write (io_unit%log,*) 'patch_t%dumpa: id, time =', self%id, self%time
  end if
  write (io_unit%dump) 1, kind(a), tag, self%id, self%level, nv, nt, nt, &
                       self%nw, self%istep, self%it, gn, n
  write (io_unit%dump) self%time, self%dtime, self%position, self%size, &
                       self%t(self%it), self%dt(self%it), self%iit(self%it)
  write (io_unit%dump) a
  call io%lock%unset ('dump')
END SUBROUTINE dump1r

!===============================================================================
!> Dump an array into data/run/rank_rrrrr.dump
!===============================================================================
SUBROUTINE dump2r (self, tag, a)
  class(patch_t):: self
  integer:: tag
  real:: a(:,:)
  !.............................................................................
  integer:: nt, nv, n(3), gn(3)
  !-----------------------------------------------------------------------------
  call io%lock%set ('dump')
  nt = 1
  nv = 1
  n  = 1
  gn = 1
  n(1:2) = shape(a)
  gn(1:2) = shape(a)
  if (self%verbose > 0) then
    write (io_unit%log,*) 'patch_t%dumpa: id, time =', self%id, self%time
  end if
  write (io_unit%dump) 2, kind(a), tag, self%id, self%level, nv, nt, nt, &
                       self%nw, self%istep, self%it, gn, n
  write (io_unit%dump) self%time, self%dtime, self%position, self%size, &
                       self%t(self%it), self%dt(self%it), self%iit(self%it)
  write (io_unit%dump) a
  call io%lock%unset ('dump')
END SUBROUTINE dump2r

!===============================================================================
!> Dump an array into data/run/rank_rrrrr.dump
!===============================================================================
SUBROUTINE dump3r (self, tag, a)
  class(patch_t):: self
  integer:: tag
  real:: a(:,:,:)
  !.............................................................................
  integer:: nt, nv
  !-----------------------------------------------------------------------------
  call io%lock%set ('dump')
  nt = 1
  nv = 1
  if (self%verbose > 0) then
    write (io_unit%log,*) 'patch_t%dumpa: id, time =', self%id, self%time
  end if
  write (io_unit%dump) 3, kind(a), tag, self%id, self%level, nv, nt, nt, &
                       self%nw, self%istep, self%it, self%gn, self%n
  write (io_unit%dump) self%time, self%dtime, self%position, self%size, &
                       self%t(self%it), self%dt(self%it), self%iit(self%it)
  write (io_unit%dump) a
  call io%lock%unset ('dump')
END SUBROUTINE dump3r

!===============================================================================
!> Dump an array into data/run/rank_rrrrr.dump
!===============================================================================
SUBROUTINE dump1d (self, tag, a)
  class(patch_t):: self
  integer:: tag
  real(8):: a(:)
  !.............................................................................
  integer:: nt, nv, n(3), gn(3)
  !-----------------------------------------------------------------------------
  call io%lock%set ('dump')
  nt = 1
  nv = 1
  n  = 1
  gn = 1
  n(1) = size(a)
  gn(1) = size(a)
  if (self%verbose > 0) then
    write (io_unit%log,*) 'patch_t%dumpa: id, time =', self%id, self%time
  end if
  write (io_unit%dump) 1, kind(a), tag, self%id, self%level, nv, nt, nt, &
                       self%nw, self%istep, self%it, gn, n
  write (io_unit%dump) self%time, self%dtime, self%position, self%size, &
                       self%t(self%it), self%dt(self%it), self%iit(self%it)
  write (io_unit%dump) a
  call io%lock%unset ('dump')
END SUBROUTINE dump1d

!===============================================================================
!> Dump an array into data/run/rank_rrrrr.dump
!===============================================================================
SUBROUTINE dump2d (self, tag, a)
  class(patch_t):: self
  integer:: tag
  real(8):: a(:,:)
  !.............................................................................
  integer:: nt, nv, n(3), gn(3)
  !-----------------------------------------------------------------------------
  call io%lock%set ('dump')
  nt = 1
  nv = 1
  n  = 1
  gn = 1
  n(1:2) = shape(a)
  gn(1:2) = shape(a)
  if (self%verbose > 0) then
    write (io_unit%log,*) 'patch_t%dumpa: id, time =', self%id, self%time
  end if
  write (io_unit%dump) 2, kind(a), tag, self%id, self%level, nv, nt, nt, &
                       self%nw, self%istep, self%it, gn, n
  write (io_unit%dump) self%time, self%dtime, self%position, self%size, &
                       self%t(self%it), self%dt(self%it), self%iit(self%it)
  write (io_unit%dump) a
  call io%lock%unset ('dump')
END SUBROUTINE dump2d

!===============================================================================
!> Dump an array into data/run/rank_rrrrr.dump
!===============================================================================
SUBROUTINE dump3d (self, tag, a)
  class(patch_t):: self
  integer:: tag
  real(8):: a(:,:,:)
  !.............................................................................
  integer:: nt, nv
  !-----------------------------------------------------------------------------
  call io%lock%set ('dump')
  nt = 1
  nv = 1
  if (self%verbose > 0) then
    write (io_unit%log,*) 'patch_t%dumpa: id, time =', self%id, self%time
  end if
  write (io_unit%dump) 3, kind(a) ,tag, self%id, self%level, nv, nt, nt, &
                       self%nw, self%istep, self%it, self%gn, self%n
  write (io_unit%dump) self%time, self%dtime, self%position, self%size, &
                       self%t(self%it), self%dt(self%it), self%iit(self%it)
  write (io_unit%dump) a
  call io%lock%unset ('dump')
END SUBROUTINE dump3d

!===============================================================================
!> Return a measure of the update cost
!===============================================================================
REAL FUNCTION cost_f (self)
  class(patch_t):: self
  !.............................................................................
  if (self%dtime > 0.0) then
    cost_f = 1d-6*product(self%n)/self%dtime
  else
    cost_f = 0.0
  end if
END FUNCTION cost_f

!===============================================================================
!> Get target and source index ranges necessary for linear interpolation from
!> source to target meshes.  If guards is true, use the full target region.
!> If all_src is true, use the full source region.
!===============================================================================
SUBROUTINE overlap (self, lt, ut, srce, ls, us, all_tgt, all_src, verbose)
  class (patch_t)        :: self, srce
  integer, dimension(3)  :: lt, ut, ls, us
  logical, optional      :: all_tgt, all_src
  integer, optional      :: verbose
  !.............................................................................
  class(mesh_t), pointer :: tm(:), sm(:)
  real(8), dimension(3)  :: tl, tu, sl, su, fl, fu, sm_p
  integer                :: verbose_l
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%begin ('patch_t%get_overlap', itimer=itimer)
  if (present(verbose)) then
    verbose_l = verbose
  else
    verbose_l = io%verbose
  end if
  !-----------------------------------------------------------------------------
  ! Self is the assumed target, patch is the assumed source
  !-----------------------------------------------------------------------------
  tm => self%mesh
  sm => srce%mesh
  !-----------------------------------------------------------------------------
  ! Compute target and source corners
  !-----------------------------------------------------------------------------
  if (present(all_tgt) .and. all_tgt) then
    tl = tm%p - tm%s*0.5d0 - tm%ng*tm%d      ! tgt lower corner
    tu = tm%p + tm%s*0.5d0 + tm%ng*tm%d      ! tgt upper corner
  else
    tl = tm%p - tm%s*0.5d0                   ! tgt lower corner
    tu = tm%p + tm%s*0.5d0                   ! tgt upper corner
  end if
  sm_p = tm%p + self%distance_to (srce)      ! source position relative target
  if (present(all_src) .and. all_src) then
    sl = sm_p - sm%s*0.5d0 - sm%ng*sm%d      ! src lower corner
    su = sm_p + sm%s*0.5d0 + sm%ng*sm%d      ! src upper corner
  else
    sl = sm_p - sm%s*0.5d0                   ! src lower corner
    su = sm_p + sm%s*0.5d0                   ! src upper corner
  end if
  !.............................................................................
  if (verbose_l > 2) then
    write (stdlog,'(a,i6,2x,3i5,5(2x,3f8.3))') &
      'patch_t%get_overlap: src mesh id,s,p,b=', &
      sm(1)%id, sm%n, sm%s, sm%p, sm%b
    write (stdlog,'(a,i6,2x,3i5,5(2x,3f8.3))') &
      'patch_t%get_overlap; tgt mesh id,s,p,b=', &
      tm(1)%id, tm%n, tm%s, tm%p, tm%b
  end if
  !-----------------------------------------------------------------------------
  ! Floating point corners of the overlap region
  !-----------------------------------------------------------------------------
  fl = max(tl,sl)                               ! overlap lower corner
  fu = min(tu,su)                               ! overlap upper corner
  if (self%level < srce%level) then
    !-------------------------------------------------------------------------
    ! Restrict operation:  Use indices that are strictly inside intersection
    !-------------------------------------------------------------------------
    ls = index_only (sm, fl, roundup=.true.)   ! src lower index
    us = index_only (sm, fu, rounddn=.true.)   ! src upper index
    tl = myposition (sm, ls)                   ! position of ls
    tu = myposition (sm, us)                   ! position of us
    lt = index_only (tm, tl, roundup=.true.)   ! target lower corner
    ut = index_only (tm, tu, rounddn=.true.)   ! target upper corner
  else
    !-------------------------------------------------------------------------
    ! Prolongation operation: Use tgt-indices that are strictly inside, and
    ! the correspondingly necessary source indices (which may fall outside
    ! the common intersection).
    !-------------------------------------------------------------------------
    lt = index_only (tm, fl, roundup=.true.)   ! target lower corner
    ut = index_only (tm, fu, rounddn=.true.)   ! target upper corner
    tl = myposition (tm, lt)                   ! position of lt
    tu = myposition (tm, ut)                   ! position of ut
    ls = index_only (sm, fl, rounddn=.true.)   ! src lower index
    us = index_only (sm, fu, roundup=.true.)   ! src upper index
  end if
  if (present(all_src) .and. all_src) then
    ls = max(ls,sm%lb)                          ! inside source halo
    us = min(us,sm%ub)                          ! inside source halo
  else
    ls = max(ls,sm%li)                          ! inside source interior
    us = min(us,sm%ui)                          ! inside source interior
  end if
  !-----------------------------------------------------------------------------
  ! This should be redundant
  !-----------------------------------------------------------------------------
  if (any(lt < tm%lb .or. ut > tm%ub)) then
    call io%warning ('target limits outside range')
  end if
  lt = max(lt,tm%lb)                            ! inside target halo
  ut = min(ut,tm%ub)                            ! inside target halo
  !-----------------------------------------------------------------------------
  ! Take care of collapsed dimensions
  !-----------------------------------------------------------------------------
  where (self%n == 1)
    lt = min(lt,tm%ui)
    ut = max(ut,tm%li)
    ls = min(ls,sm%ui)
    us = max(us,sm%li)
  end where
  !.............................................................................
  if (timer%detailed) &
    call trace%end (itimer)
!===============================================================================
contains
  function index_only (m, p, roundup, rounddn) result(i)
    class(mesh_t)    :: m(:)
    real(8)          :: p(3)
    logical, optional:: roundup, rounddn
    integer          :: i(3)
    !---------------------------------------------------------------------------
    if      (present(roundup)) then
      i = ceiling((p-m%p)/m%d + m%o)
    else if (present(rounddn)) then
      i =   floor((p-m%p)/m%d + m%o)
    else
      i =    nint((p-m%p)/m%d + m%o)
    end if
  end function
  !-----------------------------------------------------------------------------
  function myposition (m, i) result(p)
    class(mesh_t) :: m(:)
    real(8)       :: p(3)
    integer       :: i(3)
    !---------------------------------------------------------------------------
    p = m%p + m%d*(i-m%o)
  end function
END SUBROUTINE overlap

!===============================================================================
!> Mark bits where a source ROA ovelaps with a self ROI
!===============================================================================
SUBROUTINE cover (self, source, distance, bits)
  class(patch_t):: self, source
  real(8), dimension(3):: distance, fl, fu
  integer(int8):: bits(:,:,:)
  integer:: l(3), u(3)
  !-----------------------------------------------------------------------------
  fl = max (self%interest%l, source%authority%l + distance)
  fu = min (self%interest%u, source%authority%u + distance)
  l = ceiling (fl/self%mesh%d + self%mesh%o)
  u =   floor (fu/self%mesh%d + self%mesh%o)
  l = max (l, self%mesh%lb)
  u = min (u, self%mesh%ub)
  bits(l(1):u(1),l(2):u(2),l(3):u(3)) = 0
END SUBROUTINE cover

!===============================================================================
!> Check that %frefine is allocated, and allocate under lock if not
!===============================================================================
SUBROUTINE support_arrays (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  select type (self)
  class is (patch_t)
  if (.not.associated(self%frefine)) then
    call self%lock%set ('patch_t%support_arrays')
    if (.not.associated(self%frefine)) then
      allocate (self%frefine(self%gn(1),self%gn(2),self%gn(3)))
    end if
    call self%lock%unset ('patch_t%support_arrays')
  end if
  if (.not.associated(self%filled)) then
    call self%lock%set ('patch_t%support_arrays 2')
    if (.not.associated(self%filled)) then
      allocate (self%filled(self%gn(1),self%gn(2),self%gn(3)))
    end if
    call self%lock%unset ('patch_t%support_arrays 2')
  end if
  if (.not.associated(self%request)) then
    call self%lock%set ('patch_t%support_arrays 2')
    if (.not.associated(self%request)) then
      allocate (self%request(self%gn(1),self%gn(2),self%gn(3)))
    end if
    call self%lock%unset ('patch_t%support_arrays 2')
  end if
  if (.not.associated(self%outside)) then
    call self%lock%set ('patch_t%support_arrays 3')
    if (.not.associated(self%outside)) then
      allocate (self%outside(self%gn(1),self%gn(2),self%gn(3)))
    end if
    call self%lock%unset ('patch_t%support_arrays 3')
  end if
end select
END SUBROUTINE support_arrays

END MODULE patch_mod
