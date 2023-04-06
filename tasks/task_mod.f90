!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE io_mod
  USE os_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mod
  USE bits_mod
  USE trace_mod
  USE mpi_mesg_mod
  USE random_mod
  USE dll_mod
  USE aux_mod
  USE region_mod
  USE running_mean_mod
  USE shared_mod
  USE kinds_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type:: in_t
    integer:: rank=-1
  end type
  type, public:: task_t
    integer:: id = 0
    integer:: status = 0                ! 32-bit status word
    integer:: previous_status = 0       ! to detect changes
    integer:: rank = 0
    integer:: n_check = 0
    integer:: n_nbors = 0
    integer:: istep = 0
    integer:: iout = 0
    integer:: level = 0
    integer:: nt=0
    integer:: it=1
    integer:: new=1
    integer:: verbose=1
    integer:: nq=999
    integer:: n_behind=0
    integer:: restart=-9
    integer:: state=0
    integer:: parentid=0
    integer:: thread=-1                 ! OpenMP thread
    integer:: place=-1                  ! NUMA place
    integer:: n_dump=0                  ! number of dump files
    integer:: ip                        ! rank local task id
    integer:: n_needed=1                ! number of tasks that needs it
    integer:: seq=0                     ! incremental sequence number
    integer:: n_failed                  ! number of times check_ready() failed
    integer:: logging=0                 ! level of logging
    integer:: next_lb                   ! next LB step
    integer:: rank_steps=0              ! steps on rank
    integer:: check_step=0              ! step when cleared check_ready
    integer:: tag=0                     ! MPI message tag
    real(8):: cost = 0d0                ! task update cost
    integer:: key(4)=0                  ! hash key
    real(8):: time = 0d0
    real(8):: atime = 0d0               ! for active queue
    real(8):: out_next = 0d0
    real(8):: print_next = 0d0
    real(8):: dtime = 0d0
    real(4):: courant_info = 0.0        ! info for gpatch_t%info
    real(8):: min_dtime = 0d0
    real(4):: grace = 0.0
    real(8):: position(3) = 0d0
    real(8):: origin(3)=0d0             ! the cartesian box origin (for BCs)
    real(8):: radius=1d0                ! generic radius (do NOT set to 0d0)
    real(8):: velocity(3) = 0d0         ! velocity of task in Cartesian coords.
    real(8):: size(3) = 0.9999d0
    real(8):: gsize(3)                  ! convenience
    real(8):: ds(3)=0d0                 ! convenience
    real(8):: box(3)=1d0                ! the periodic wrapping size
    real(8):: wallclock = 0d0           ! wallclock time
    real(8):: update_time = 0d0         ! update time
    real(8):: unpack_time = 0d0         ! unpack time
    real(8):: dnload_time = -1d0        ! last dnload time
    real(8):: mulligan_time = -1d0      ! time to rewind to
    real(8):: dt_small=1d-7             ! smallest allowed dtime
    real(8):: dt_fixed=-1d0             ! for fixed time step Courant condition
    real(8):: dt_extra=0d0              ! extra step (cf. courant_condition)
    real(8):: mb=0d0                    ! MB memory
    real   :: quality = 0.0             ! floating point level
    real   :: safety_factor=1.0         ! time step change
    character(len=64):: kind='task'     ! solver name
    character(len=12):: type='task_t'   ! task type
    real:: latency=0.0                  ! MPI message latency
    real(8):: wc_last=0.0_8
    real(8):: sync_time=0d0
    real(8):: update_last=0.0_8
    real:: update_cadence=0.0
    logical:: syncing=.false.
    logical:: track=.false.
    logical:: periodic(3)=.true.
    logical:: rotated=.false.
    logical:: limit_dtime=.false.
    logical:: dump_it=.false.
    logical:: finished=.false.
    logical:: allow_extrapolate=.false.
    logical:: pre_initialized=.false.
    type(random_t):: random             ! pseudo-random number generator
    type(aux_t):: aux                   ! auxiliary output
    integer(8):: amr_offset=0_8         ! byte offset in amr_io
    type(region_t):: authority          ! region of authority
    type(region_t):: interest           ! region of interest
    type(running_mean_t):: running_cost ! running averge update cost
    type(in_t):: in
    !---------------------------------------------------------------------------
    ! Pointer attributes -- make sure to add to task_t%reset when adding, to
    ! avoid that cloned instances inadvertenly refer to parent pointers
    !---------------------------------------------------------------------------
    character(len=64),     pointer:: name      => null()
    class(task_t),         pointer:: parent    => null()
    class(mesg_t),         pointer:: mesg      => null()
    class(mesg_list_t),    pointer:: sent_list => null()
    integer, dimension(:), pointer:: iit       => null()
    real(8), dimension(:), pointer:: t         => null()
    real(8), dimension(:), pointer:: dt        => null()
    type(lock_t),          pointer:: lock      => null()
    type(lock_t),          pointer:: mulligan_lock => null()
  contains
    procedure:: pre_init
    procedure:: pre_update
    procedure:: post_update
    procedure:: init
    procedure:: init_unique
    procedure:: init_level
    procedure:: dealloc
    procedure:: dnload
    procedure:: update => void
    procedure:: header_to_task
    procedure:: rotate
    procedure:: derotate
    procedure:: info
    procedure:: timeslots
    procedure:: is_ahead_of
    procedure:: set
    procedure:: clear
    procedure:: print
    procedure:: is_set
    procedure:: is_clear
    procedure:: allocate => void
    procedure:: load_balance
    procedure:: has_finished
    procedure:: overlaps
    procedure:: overlaps_point
    procedure:: nbor_relations
    procedure:: distance_to
    procedure:: distance_from
    procedure:: distance
    procedure:: point_distance
    procedure:: pack
    procedure:: unpack
    procedure:: is_relevant
    procedure:: test_update
    procedure:: solver_is
    procedure:: debug
    procedure:: contains1      ! patch contains a point
    procedure:: contains2      ! patch contains a patch
    generic:: contains => contains1, contains2
    procedure:: task_info
    procedure:: log
    procedure:: run => void
    procedure:: statistics
    procedure:: bndries
    procedure:: periodic_grid
    procedure:: courant_condition
    procedure:: output
    procedure:: output_metadata
    procedure:: input_metadata
    procedure:: allocate_mesg
    procedure:: roi
    procedure:: roa
    procedure:: update_position
    procedure:: get_position
    procedure:: velocity_shift
    procedure:: hash_key
    procedure:: hash_position
    procedure:: set_properties
    procedure:: check_refine
    procedure:: write_nml
    procedure:: reset
    procedure:: periodic_wrap4
    procedure:: periodic_wrap8
    generic:: periodic_wrap => periodic_wrap4, periodic_wrap8
    procedure:: set_rank
    ! ----- Void overloadable subroutines
    procedure:: allocate_generic => void
    procedure:: deallocate_generic => void
    procedure:: to_generic => void
    procedure:: from_generic => void
    procedure:: init_bdries => void
  end type
  !-----------------------------------------------------------------------------
  ! Sequential header, for MPI and I/O
  !-----------------------------------------------------------------------------
  integer, save:: n_header = (14*4 + 4*8 + 5*8*3 + 1*8 + 2*8*8 + 1*4*8 + 2*4*3)/4
  type:: header_t
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
    integer:: key(4)                    ! task hash key
  end type
  logical:: first_time=.true.
  integer:: track_id=0
  integer:: verbose=0
  type(dll_t), pointer:: deleted_list
CONTAINS

!===============================================================================
!> For timing
!===============================================================================
SUBROUTINE pre_update (self)
  class(task_t), target:: self
  self%wallclock = wallclock()
END SUBROUTINE pre_update

!===============================================================================
SUBROUTINE post_update (self)
  class (task_t):: self
  call trace%begin ('task_t%post_update')
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
SUBROUTINE void (self)
  class (task_t):: self
  !.............................................................................
  call trace%begin ('task_t%void')
  call io%abort('task_t%void called')
  call trace%end ()
END SUBROUTINE void

!===============================================================================
SUBROUTINE dnload (self, only)
  class (task_t):: self
  integer, optional:: only
  call trace%begin ('task_t%dnload')
  call io%warning('task_t%dnload called')
  call trace%end ()
END SUBROUTINE dnload

!===============================================================================
!> When restarting the global shared%rank_id must be set to a value larger
!> than or equal to the previously largest value
!===============================================================================
SUBROUTINE pre_init (self)
  class (task_t), target:: self
  !.............................................................................
  self%out_next = io%start_time
  self%aux%aux_next = io%start_time
  call self%set (bits%internal)                   ! to be revised later
END SUBROUTINE pre_init

!===============================================================================
!> Define unique task id, unless already set before entry.
!===============================================================================
SUBROUTINE init (self)
  class (task_t):: self
  integer:: id, id_master, id_tmp
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%init', itimer=itimer)
  call init_params (self)
  if (self%id == 0) then
    call init_unique (self)
    if (verbose > 0) &
      write (stdlog,*) 'task_t%init: id =', self%id
  else
    !---------------------------------------------------------------------------
    ! If a task already has an id, update the MPI counter if its value is less.
    ! The id of a virtual copy of a new AMR task should NOT trigger this.
    !---------------------------------------------------------------------------
    id_master = mpi%id%update (0)
    id_tmp = id_master
    if (self%id > id_master) then
      id_master = mpi%id%update (self%id-id_master+1)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Locks, random numbers, etc
  !-----------------------------------------------------------------------------
  call io%lock%set ('task_t%init')
  if (.not.associated(self%lock)) then
    allocate (self%lock, self%mulligan_lock)
    call self%lock%init ('task')
    call self%mulligan_lock%init ('mull')
  end if
  if (.not.self%lock%initialized) then
    call self%lock%init ('task')
  end if
  call io%lock%unset ('task_t%init')
  call self%random%init
  call self%aux%init
  !-----------------------------------------------------------------------------
  ! Only initialize time slot values when first allocated
  !-----------------------------------------------------------------------------
  if (.not.associated(self%t)) then
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
    self%dt = self%dtime                          ! in case already set
    self%t  = self%time                           ! in case already set
    self%iit = 1                                  ! initial slots
    self%it = 1                                   ! slot for current values
    self%new = 2                                  ! slot for new values
    self%iit(self%nt) = self%new                  ! slot for new values
    self%time = self%t(self%it)                   ! initial time
  end if
  if (self%id == track_id) then
    self%track = .true.
  end if
  self%thread = omp%thread                        ! thread affinity
  self%place = omp%place                          ! memory affinity
  self%key = self%hash_key ()                     ! cache the hash key
  call self%log ('init', 3)
  call trace%end (itimer)
END SUBROUTINE init

!===============================================================================
!> Read task_params, allowing to tag task IDs for dumping
!===============================================================================
SUBROUTINE init_params (self)
  class (task_t):: self
  character(len=120):: ids= &
    '$Id$'
  integer:: iostat, id=0, id_master
  logical, save:: do_dumps=.false., dump_all=.false.
  integer, save:: id_max=0
  namelist /task_params/ track_id, do_dumps, dump_all, verbose
  namelist /dump_params/ id
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%init_params')
  call trace%print_id (ids)
  call io%lock%set ('task_t%init_params')
  !-----------------------------------------------------------------------------
  ! First thread to enter reads task_params
  !-----------------------------------------------------------------------------
  if (first_time) then
    rewind (io%input)
    read (io%input, task_params, iostat=iostat)
    call io%namelist_warning ('task_params', iostat)
    write (io%namelist, task_params)
    first_time = .false.
  end if
  self%out_next = io%start_time
  self%aux%aux_next = io%start_time
  !-----------------------------------------------------------------------------
  ! If dumps are requested, read a list of task ids for which dumps should be on
  !-----------------------------------------------------------------------------
  if (dump_all) then
    self%dump_it = .true.
    self%verbose = verbose
  else if (do_dumps) then
    rewind (io%input)
    do
      iostat = 0
      verbose = 0
      read (io%input, dump_params, iostat=iostat)
      if (iostat /= 0) then
        exit
      else if (id == self%id) then
        self%verbose = verbose
        self%dump_it = .true.
        write (stdlog,*) 'task_t%init: dump_id =', id
        write (io%namelist, dump_params)
        exit
      end if
    end do
  end if
  call io%lock%unset ('task_t%init_params')
  call trace%end ()
END SUBROUTINE init_params

!===============================================================================
!> Provide a unique task ID, and call init_params.  The optinal "same" parameter
!> results in the same IDs being assigned on all ranks.  This may only be used
!> initially, since it would collide with the normal, unique set of IDs if used
!> after those have started to be used.
!===============================================================================
SUBROUTINE init_unique (self, same)
  class (task_t):: self
  logical, optional:: same
  logical:: same_l
  integer:: id, id_master
  integer, save:: id_same=0
  !-----------------------------------------------------------------------------
  ! Ranks that call this routine in sync  with same=.true. get a matching set
  ! of IDs to prune from, while the master rank updates the global counter
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%init_unique')
  call init_params (self)
  if (present(same)) then
    same_l = same
  else
    same_l = .false.
  end if
  if (same_l) then
    if (mpi%master) &
      id_master = mpi%id%update (1)
    !$omp atomic capture
    id_same = id_same+1
    id = id_same
    !$omp end atomic
    if (verbose > 0) &
      write (stdlog,*) 'task_t%init_unique(same): id =', &
        id, id_same, id_master
  !-----------------------------------------------------------------------------
  ! Calls without the same option return unique IDs, which do not overlap with
  ! matching set of IDs
  !-----------------------------------------------------------------------------
  else
    id_master = mpi%id%update (1)
    do while (id_master == id_same)
      write (stderr,*) 'task_t%init_unique ERROR: same value returned', id_same
      write (io_unit%mpi,*) 'task_t%init_unique ERROR: same value returned', id_same
      id_master = mpi%id%update (1)
    end do
    !$omp atomic write
    id_same = id_master-1
    id = id_master
    if (verbose > 0) &
      write (stdlog,*) 'task_t%init_unique(else): id =', &
        id, id_same, id_master
  end if
  self%id = id
  call trace%end ()
END SUBROUTINE init_unique

!===============================================================================
SUBROUTINE init_level (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  call io%abort ('task_t%init_level called')
END SUBROUTINE init_level

!===============================================================================
!> Deallocate allocated attributes
!===============================================================================
SUBROUTINE dealloc (self)
  class(task_t):: self
  class(mesg_t), pointer:: mesg, next
  character(2):: label
  !.............................................................................
  call trace%begin ('task_t%dealloc')
  if (associated(self%t  )) call os%deallocate (self%t  )
  if (associated(self%dt )) call os%deallocate (self%dt )
  if (associated(self%iit)) call os%deallocate (self%iit)
  if (associated(self%sent_list)) then
    write (label,'(i2)') self%sent_list%n
    call self%log ('dealloc mesg_list '//label, 3)
    call self%sent_list%dealloc (self%mesg)
    call os%mem_register (-storage_size(self%sent_list), 1, 'sent_list')
    deallocate (self%sent_list)
  end if
  if (associated(self%mesg)) then
    call self%log ('dealloc mesg', 3)
    call self%mesg%dealloc ()
    deallocate (self%mesg)
  end if
  call trace%end ()
END SUBROUTINE dealloc

!===============================================================================
!> Save the id of a deleted task, for resuse
!===============================================================================
SUBROUTINE save_id (self)
  class (task_t):: self
  !.............................................................................
  class(dll_node_t), pointer:: node
  integer, pointer:: id
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%save_id')
  if (.not.associated(deleted_list%head)) then
    call deleted_list%init
  end if
  allocate (node, id)
  id = self%id
  node%car => id
  call deleted_list%append (node)
  call trace%end ()
END SUBROUTINE save_id

!===============================================================================
!> Copy he most important fields to a sequenced header, for communication and I/O
!===============================================================================
SUBROUTINE header_to_task (self, head)
  class(task_t):: self
  type(header_t):: head
  integer:: nt, new, i, save_bits
  !-----------------------------------------------------------------------------
  ! Swap the roles of boundary and virtual on incoming tasks.  Also clear the
  ! ready bit, which is an indicator of a task being in the ready queue, which
  ! an incoming task is not.  By doing this here, we avoid having to do it in
  ! a critical region, otherwise needed to protect an operation on the task
  !-----------------------------------------------------------------------------
  call trace%begin ('task%header_to_task')
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
  ! Except for the first and last appearance, istep should increase be one
  !-----------------------------------------------------------------------------
  if (head%id /= self%id) then
    call io%abort('wrong message ID unpacked')
  end if
  if ((self%istep > 0) .and. (iand(head%status, bits%remove) == 0)) then
    if (head%istep > self%istep+1) then
      write(stdlog,1) &
        'WARNING: early arrival in time reversal id, rank =', &
        self%id, mpi%rank, head%istep, self%istep, head%time, self%time
    1 format(a,4i6,2f12.6)
      flush (stdlog)
    else if (head%istep < self%istep+1) then
      write(stdlog,1) &
        'WARNING: late arrival in time reversal id, rank =', &
        self%id, mpi%rank, head%istep, self%istep, head%time, self%time
      flush (stdlog)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! If the task is virtual, preserve bits%init_nbors and bits%update_status
  ! across the task update
  !-----------------------------------------------------------------------------
                                         call self%lock%set ('header_to_task')
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
  self%time       = head%time
  self%out_next   = head%out_next
  self%print_next = head%print_next
  self%dtime      = head%dtime
  self%position   = head%position
  self%velocity   = head%velocity
  self%size       = head%size
  self%rank_steps = head%rank_steps
  !$omp atomic write
  self%it         = head%it
  new             = mod(head%it,head%nt) + 1
  !$omp atomic write
  self%new        = new
  self%nt         = head%nt
  nt              = head%nt
  if (nt > 0) then
    do i=1,nt
      !$omp atomic write
      self%t(i)   = head%t(i)
      !$omp atomic write
      self%dt(i)  = head%dt(i)
      !$omp atomic write
      self%iit(i) = head%iit(i)
    end do
    !self%new      = head%iit(nt)
  else
    write (stdlog,*) omp%thread,'header_to_task WARNING: nt non-positive', nt
    flush (stdlog)
  end if
  if (head%istep == 0) then
    write (stdlog,'(1x,a,i6,f12.6,2x,9f12.6)') &
      'task%header_to_task: id, time, t =', &
      self%id, self%time, self%t
    write (stdlog,'(1x,a,i6,2x,3i4,2x,9i4)') &
      'task%header_to_task: id, istep, it, new, iit =', &
      self%id, self%istep, self%it, self%new, self%iit
  end if
  if (self%it > self%nt .or. self%it < 1) then
    write (stdlog,*) 'header_to_task ERROR: it =', self%it, self%iit(nt-1)
    self%it = self%iit(nt-1)
  end if
  self%nq         = head%nq
  self%latency    = wallclock() - head%send_time
                                         call self%lock%unset ('header_to_task')
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
  call trace%end ()
END SUBROUTINE header_to_task

!===============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!===============================================================================
SUBROUTINE rotate (self)
  class(task_t):: self
  integer:: iv, nv, new, i
  real(8):: time
  integer, save:: itimer=0
  real(8), pointer:: rptr
  integer, pointer:: iptr
  !-----------------------------------------------------------------------------
  ! Implement setting of per-patch parameters here
  !-----------------------------------------------------------------------------
  if (shared%benchmark > 0) return
  if (io%processing > 0d0) then
    if (io%out_next > 0.0) self%out_next = io%out_next
    if (io%grace    > 0.0) self%grace    = io%grace
  end if
  if (self%rotated) return
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%rotate',itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Call timestep update procedure to fill the next time slot with values
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    write (stdlog,'(a,7i10)')   'task_t%rotate iit =', self%iit
    flush (stdlog)
  end if
  !-----------------------------------------------------------------------------
  ! If this task has self%syncing set, then the time step has been set so as to
  ! arrive exactly ast sync_time.  To precent problems with round-off we set the
  ! time exactly.
  !-----------------------------------------------------------------------------
  if (self%syncing) then
    time = self%sync_time
    self%syncing = .false.
    if (io%verbose>2) &
      write (io_unit%mpi,*) 'task_t%rotate:: turning off syncing for now', self%id
  else
    time = self%time + self%dtime
  end if
  !-----------------------------------------------------------------------------
  ! Time slot rotation.
  !-----------------------------------------------------------------------------
  ! call                                                  self%lock%set ('rotate')
  !-----------------------------------------------------------------------------
  ! First update the t and dt of the new time slot, and only then the index.
  ! This ensures that other tasks always see a consistent state.  The new
  ! state is NOT to be referenced by any other task, and hence those updates
  ! do not have to be atomic!
  !-----------------------------------------------------------------------------
  self%dt(self%new) = self%dtime                      ! save dt
  self%t (self%new) = time                            ! initial time
  self%it = self%new                                  ! update time slot, new
  new = mod(self%new,self%nt) + 1                     ! this thread owns it
  self%new = new                                      ! circular increment
  !-----------------------------------------------------------------------------
  ! Now rotate the time slot array, which needs to be consistent with the new
  ! %it and %new values as quikly as possible, if/when not under lock. These
  ! are used when finding the interpolation slots and weights, and jence that
  ! operation (guard zone downloading) must really be done under lock (using
  ! an nbor list copy to avoid deadlock)
  !-----------------------------------------------------------------------------
  self%iit(1:self%nt-1) = self%iit(2:self%nt)         ! rotate left
  self%iit(self%nt) = self%new                        ! top slot
  self%dt(self%new)= self%dtime                       ! estimate
  !-----------------------------------------------------------------------------
  ! By waiting with update task%time until after the time slots are rotated
  ! one prevents the next guard zone downloading from starting with until
  ! the slot values are correct
  !-----------------------------------------------------------------------------
  !$omp atomic write
  self%time = time                                    ! new time
  !-----------------------------------------------------------------------------
  ! These updates also must be done atomically, but are not operation critical
  !-----------------------------------------------------------------------------
  self%istep = self%istep + 1                         ! increment task step
  self%rank_steps = self%rank_steps + 1               ! increment rank step
  !-----------------------------------------------------------------------------
  ! These are not time critical
  !-----------------------------------------------------------------------------
  self%rotated = .true.
  ! call                                                self%lock%unset ('rotate')
  !-----------------------------------------------------------------------------
  ! Compute update cost
  !-----------------------------------------------------------------------------
  time = wallclock()
  self%update_time = time - self%wallclock
  !if (self%dtime > 0d0) then
  !  self%cost = self%running_cost%value (real(self%update_time/self%dtime))
  !else
  !  self%cost = self%running_cost%average
  !end if
  !.............................................................................
  if (io%verbose>2) then
   write (stdlog,*) 'task id:',self%id,'  it:',self%it,'  new:',self%new
   if (io%verbose>3) then
     write (stdlog,*) 'task iit:', self%iit
     write (stdlog,*) 'task t(iit):', self%t(self%iit)
    end if
  end if
  call self%log ('rotate', 3)
  call trace%end (itimer)
END SUBROUTINE rotate

!===============================================================================
!> Derotate task time slots, effectively backing up one time step.  Task locking
!> should be done by the calling procedure.
!===============================================================================
SUBROUTINE derotate (self)
  class(task_t):: self
  integer:: i
  real(8):: out_next
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%deroate')
  if (self%time > 0d0 .and. self%istep > 0) then
    self%iit = modulo(self%iit-2,self%nt) + 1
    self%istep = self%istep - 1
    self%it = self%iit(self%nt-1)
    self%new = self%iit(self%nt)
    self%time = self%t(self%it)
    self%dtime = self%dt(self%it)
    out_next = self%out_next
    self%out_next = (floor(self%time/io%out_time)+1)*io%out_time
    if (self%out_next < out_next) then
      self%iout = self%iout - nint((out_next-self%out_next)/io%out_time)
      write (stderr,*) 'derotate: id, iout =', self%id, self%iout
    end if
  end if
  call trace%end ()
END SUBROUTINE derotate

!===============================================================================
!> Print info to stdout
!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(task_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !.............................................................................
END SUBROUTINE info

!===============================================================================
!> Get timeslot info atomically.
!===============================================================================
SUBROUTINE timeslots (self, slots, times)
  class(task_t):: self
  real(8):: times(self%nt-1)
  integer:: slots(self%nt-1)
  integer:: i, j
  !-----------------------------------------------------------------------------
  ! Cache the current slot index, and the times.  If the times happen to
  ! advance in the middle of the loop, the time may jump over one slot
  !-----------------------------------------------------------------------------
                                                call self%lock%set ('timeslots')
  do i=1,self%nt-1
    j = self%iit(i)
    times(i)  = self%t(j)
    slots(i) = j
  end do
                                              call self%lock%unset ('timeslots')
END SUBROUTINE timeslots

!===============================================================================
!> Get timeslot info atomically.  Note that, while the thread updating the task
!> set the time and dtime of slot (nt) which, in terms of the indices retrieved
!> here is outside the range of concern, in rare cases it may be updated again,
!> in which case it will access the current slot (1) the next time, then (2)
!> and so on; therefore, it is prudent to use a task lock shared with rotate
!===============================================================================
SUBROUTINE timeslots2 (self, slots, times)
  class(task_t):: self
  real(8):: times(self%nt-1), t(self%nt)
  integer:: slots(self%nt-1)
  integer:: i, j, o, it, nt
  !-----------------------------------------------------------------------------
  ! Cache the current slot index, and the times.  If the times happen to
  ! advance in the middle of the loop, the time may jump over one slot
  !-----------------------------------------------------------------------------
  !$omp atomic read
  it = self%it
  !-----------------------------------------------------------------------------
  ! Below, index i corresponds to the index into iit(:), and j corresponds to
  ! iit(i); the index that makes t(iit(i)) monotonous.
  !-----------------------------------------------------------------------------
  do i=1,self%nt
    !$omp atomic read
    t(i)  = self%t(i)
  end do
  nt = self%nt
  i = it                                ! current slot
  do j=nt-1,1,-1
    slots(j) = i                        ! slot value
    times(j) = t(i)                     ! time value
    i = 1 + modulo(i-2,nt)              ! decrease, wrapping
  end do
  !-----------------------------------------------------------------------------
  ! Check if the times are actually equal; if so, the slots should be as well
  !-----------------------------------------------------------------------------
  do j=nt-2,1,-1
    if (times(j)==times(j+1)) then
      slots(j) = slots(j+1)
    end if
  end do
END SUBROUTINE timeslots2

!===============================================================================
!> Check if self (which is an nbor task) is ahead of target (which is the one to
!> possibly move to the ready queue), using self%dtime*target%grace as the grace
!> period, since we want to limit the extrapolation in the nbor task to at most
!> target%grace*self%dtime.
!===============================================================================
LOGICAL FUNCTION is_ahead_of (self, target)
  class (task_t):: self
  class (task_t), pointer:: target
  !.............................................................................
  real:: grace
  real(8):: nbtime, nbdtime, tgtime
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! dt(self%dt) is always consistent with t(:), while dtime changes in
  ! courant_condition(), before other attributes change in rotate()
  !-----------------------------------------------------------------------------
  !$omp atomic read
  tgtime = target%time
                                         call self%lock%set_read ('is_ahead_of')
  nbtime  = self%t (self%it)
  nbdtime = max (self%dt(self%it), 1d-30)   ! ensure non-zero
  if (self%is_set(bits%frozen)) then
    is_ahead_of = .true.
  else if (self%level < target%level) then  ! lower levels go first
    is_ahead_of = nbtime >= tgtime
  else if (self%istep < 3) then
    is_ahead_of = nbtime >= tgtime
  else
    grace = self%grace
    is_ahead_of = nbtime + nbdtime*grace >= tgtime
  end if
                                       call self%lock%unset_read ('is_ahead_of')
END FUNCTION is_ahead_of

SUBROUTINE print (self)
  class (task_t):: self
END SUBROUTINE

!===============================================================================
!> Set status bits
!===============================================================================
SUBROUTINE set (self, bits)
  class (task_t):: self
  integer:: bits
  !.............................................................................
!  call trace%begin ('task_t%set')
  !$omp atomic
  self%status = ior(self%status, bits)
!  call trace%end
END SUBROUTINE

!===============================================================================
!> Clear status bits
!===============================================================================
SUBROUTINE clear (self, bits)
  class (task_t):: self
  integer:: bits
  !.............................................................................
!  call trace%begin ('task_t%clear')
  !$omp atomic
  self%status = iand(self%status, not(bits))
!  call trace%end
END SUBROUTINE

!===============================================================================
!> Check if ANY of the bits is set (argument can be sum of bits)
!===============================================================================
LOGICAL FUNCTION is_set (self, bits)
  class (task_t):: self
  integer:: bits, status
  !.............................................................................
!  call trace%begin ('task_t%is_set')
  !$omp atomic read
  status = self%status
  is_set = (iand(status, bits) /= 0)
!  call trace%end
END FUNCTION

!===============================================================================
!> Check that ALL bits are clear (argument can be sum of bits)
!===============================================================================
LOGICAL FUNCTION is_clear (self, bits)
  class (task_t):: self
  integer:: bits, status
  !.............................................................................
!  call trace%begin ('task_t%is_clear')
  !$omp atomic read
  status = self%status
  is_clear = (iand(status, bits) == 0)
!  call trace%end
END FUNCTION

!===============================================================================
LOGICAL FUNCTION has_finished (self)
  class (task_t):: self
  !.............................................................................
  has_finished = (self%time > io%end_time)
END FUNCTION

!===============================================================================
!===============================================================================
SUBROUTINE load_balance (self)
  class(task_t):: self
END SUBROUTINE load_balance

!===============================================================================
!> Check if there is overlap. The definition of overlap MUST be symmetric, since
!> it is being used to define neighbor relations.
!===============================================================================
FUNCTION overlaps (self, task)
  class(task_t):: self
  class(task_t):: task
  logical:: overlaps
  !-----------------------------------------------------------------------------
  overlaps = self%interest%overlaps (task%interest, self%distance_to (task)) &
        .or. task%interest%overlaps (self%interest, task%distance_to (self))
END FUNCTION overlaps

!===============================================================================
!> Check if point p overlaps with the internal region of self
!===============================================================================
FUNCTION overlaps_point (self, p)
  class(task_t):: self
  real(8):: p(3), d(3)
  logical:: overlaps_point
  integer, save:: itimer=0
  !.............................................................................
  if (timer%detailed) &
    call trace%begin ('task_t%overlaps_point', itimer=itimer)
  d = self%point_distance (self%position, p)
  overlaps_point = all(abs(d) < self%size*0.5d0)
  if (timer%detailed) &
    call trace%end (itimer)
END FUNCTION overlaps_point

!===============================================================================
!> Decide if another task is needed, if it needs_me, and if it should be
!> downloaded (= used in one way or another to prepare the task for update)
!> This is the lowest, most generic relations call. The checks are based solely
!> on ranks.  A proxy task is used to determine the status bits from scratch,
!> which allows determination of transition from boundary to internal, or from
!> virtual to external.
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(task_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.associated(another)) then
    write (stderr,*) mpi%rank, 'another not associated'
    call io%abort ('task_t another not associated')
  end if
  if (timer%detailed) &
    call trace%begin ('task_t%nbor_relations', 2, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If this thread owns the task, start by assuming either internal or external
  !-----------------------------------------------------------------------------
  if (self%rank == mpi%rank) then
    !---------------------------------------------------------------------------
    ! nbor is frozen
    !---------------------------------------------------------------------------
    if (another%is_set (bits%frozen)) then
      !$omp atomic write
      needs_me = .false.
      !$omp atomic write
      needed = .false.
    !---------------------------------------------------------------------------
    ! nbor is virtual
    !---------------------------------------------------------------------------
    else if (another%rank /= mpi%rank) then
      ! -- use in check_ready(), but not in check_nbors()
      !$omp atomic write
      needed = .true.
      !$omp atomic write
      needs_me = .false.
    !---------------------------------------------------------------------------
    ! nbor is real
    !---------------------------------------------------------------------------
    else
      ! -- use in check_ready(), and in check_nbors()
      !$omp atomic write
      needed = .true.
      !$omp atomic write
      needs_me = .true.
    end if
    !---------------------------------------------------------------------------
    ! real tasks need guard zones
    !---------------------------------------------------------------------------
    !$omp atomic write
    download = .true.
  !-----------------------------------------------------------------------------
  ! task is virtual
  !-----------------------------------------------------------------------------
  else
    ! write (stdout,*) self%id, another%id, another%is_set (bits%frozen)
    !---------------------------------------------------------------------------
    ! nbor is frozen
    !---------------------------------------------------------------------------
    if (.not.associated(another)) then
      write (stderr,*) mpi%rank, 'another not associated'
      call io%abort ('task_t another not associated')
    end if
    if (another%is_set (bits%frozen)) then
      !$omp atomic write
      needs_me = .false.
      !$omp atomic write
      needed = .false.
    end if
    !---------------------------------------------------------------------------
    ! task is virtual, nbor is real
    !---------------------------------------------------------------------------
    if (another%rank==mpi%rank) then
      ! -- don't use in check_ready(), use in check_nbors()
      !$omp atomic write
      needed = .false.
      !$omp atomic write
      needs_me = .true.
    !---------------------------------------------------------------------------
    ! task is virtual, nbor is virtual
    !---------------------------------------------------------------------------
    else
      !$omp atomic write
      needed = .false.
      !$omp atomic write
      needs_me = .false.
    end if
    !---------------------------------------------------------------------------
    ! virtual tasks don't need guard zones
    !---------------------------------------------------------------------------
    download = .false.
  end if
  !-----------------------------------------------------------------------------
  ! nbor is more than one level different; these are normally not needed,
  ! and are only used in support considerations, so should not take part
  ! in check_nbors() and check_ready(), but may be needed exceptionally
  !-----------------------------------------------------------------------------
  if (abs(another%level-self%level) > 1) then
    !$omp atomic write
    needed   = .false.    ! Do NOT comment this out -- check_ready fails!
    !$omp atomic write
    needs_me = .false.    ! Do NOT comment this out -- check_ready fails!
    !$omp atomic write
    download = .false.
  end if
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%end (itimer)
END SUBROUTINE nbor_relations

!===============================================================================
!> Disambiguation: distance from self to task -- applies (only) to Cartesian
!> meshes
!===============================================================================
FUNCTION distance_to (self, task) RESULT (d)
  class(task_t):: self
  class(task_t):: task
  real(8):: d(3), p(3)
  !-----------------------------------------------------------------------------
  d = task%position - self%position
  p = self%box
  where (self%periodic)
    d = modulo (d + 0.5_8*p, p) - 0.5_8*p
  end where
END FUNCTION distance_to

!===============================================================================
!> Disambiguation
!===============================================================================
FUNCTION distance_from (self, task)
  class(task_t):: self
  class(task_t):: task
  real(8):: distance_from(3)
  distance_from = distance (self, task)
END FUNCTION distance_from

!===============================================================================
!> Signed distance between the centers of two tasks in a possibly periodic box
!>
!> NOTE: This function is normally NOT used (but remains so it can be overloaded
!> or explicitly called), since it is overloaded by the patch_t%distance, which
!> is patch geometry aware.
!===============================================================================
FUNCTION distance (self, task) RESULT (out)
  class(task_t):: self
  class(task_t):: task
  real(8):: out(3)
  integer, save:: itimer=0
  !.............................................................................
  if (timer%detailed) then
    call trace%begin ('task_t%distance', itimer=itimer)
    out = self%point_distance (self%position, task%position)
    call trace%end (itimer)
  else
    out = self%point_distance (self%position, task%position)
  end if
END FUNCTION distance

!===============================================================================
!> Signed distance between two points in a possibly periodic box
!===============================================================================
FUNCTION point_distance (self, point1, point2) RESULT (out)
  class(task_t):: self
  real(8), intent(in):: point1(3), point2(3)
  real(8):: out(3)
  integer, save:: itimer=0
  !.............................................................................
  if (timer%detailed) then
    call trace%begin ('task_t%point_distance', itimer=itimer)
    out = point1-point2
    where (self%periodic .and. self%box /= 0d0)
      out = modulo (out+0.5_8*self%box, self%box) - 0.5_8*self%box
    end where
    call trace%end (itimer)
  else
    out = point1-point2
    where (self%periodic .and. self%box /= 0d0)
      out = modulo (out+0.5_8*self%box, self%box) - 0.5_8*self%box
    end where
  end if
END FUNCTION point_distance

SUBROUTINE pack (self)
  class(task_t):: self
  call mpi%abort ('task_t%pack called')
END SUBROUTINE

SUBROUTINE unpack (self, buffer)
  class(task_t):: self
  integer, pointer:: buffer(:)
  call mpi%abort ('task_t%unpack called')
END SUBROUTINE

!===============================================================================
!> Stub, to be overloaded in experiment_mod to drop patches
!===============================================================================
FUNCTION is_relevant (self)
  class(task_t):: self
  logical:: is_relevant
  !-----------------------------------------------------------------------------
  is_relevant = .true.
END FUNCTION is_relevant

!===============================================================================
!> Test update procedure for dispatcher tests
!===============================================================================
SUBROUTINE test_update (self)
  class(task_t):: self
  self%dtime = self%random%ran1()
  self%dt(self%it) = self%dtime
END SUBROUTINE test_update

!===============================================================================
!> Test if solver leading characters match a given string
!===============================================================================
FUNCTION solver_is (self, kind)
  class(task_t):: self
  logical solver_is
  character(len=*):: kind
  solver_is = self%kind(1:len(kind)) == kind
END FUNCTION solver_is

!===============================================================================
!> Return true if a point is inside the self interior
!===============================================================================
LOGICAL FUNCTION contains1 (self, point)
  class(task_t):: self
  real(8):: point(3), delta(3)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) then
    call trace%begin ('task_t%contains1', itimer=itimer)
    delta = self%point_distance (self%position, point)
    contains1 = all(abs(delta) <= (0.5_8+tiny(1d0))*self%size)
    call trace%end (itimer)
  else
    contains1 = all(abs(self%point_distance (self%position, point)) &
                    <= (0.5_8+tiny(1d0))*self%size)
  end if
END FUNCTION contains1

!===============================================================================
!> Return true if a patch position is inside the self interior
!===============================================================================
LOGICAL FUNCTION contains2 (self, patch)
  class(task_t):: self, patch
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) then
    call trace%begin ('task_t%contains2', itimer=itimer)
    contains2 = contains1 (self, patch%position)
    call trace%end (itimer)
  else
    contains2 = contains1 (self, patch%position)
  end if
END FUNCTION contains2

!===============================================================================
!> Turn on debug output at and beyond a given level, and for id_debug task
!===============================================================================
FUNCTION debug (self, level)
  class(task_t):: self
  integer:: level
  logical debug
  !-----------------------------------------------------------------------------
  debug = (io%verbose >= level) .or. (self%id == io%id_debug)
END FUNCTION debug

!===============================================================================
SUBROUTINE task_info (self, unit)
  class(task_t):: self
  integer, optional:: unit
  integer:: parent, unit_L
  !-----------------------------------------------------------------------------
  unit_L = shared%option (unit, io_unit%mpi)
  if (associated(self%parent)) then
    parent = self%parent%id
  else
    parent = 0
  end if
  write(unit_L,'(1x,a,4x,i6,":",i2.2,3i4,3i6,5x,2l1,2x,3f9.3)') &
    'task_t%task_info: id,key,nnb,pid,rank,BV,pos=', self%id, self%key, &
    self%n_nbors, parent, self%rank, self%is_set(bits%boundary), &
    self%is_set(bits%virtual), self%position
END SUBROUTINE task_info

!===============================================================================
SUBROUTINE log (self, label, level)
  class(task_t):: self
  character(len=*):: label
  integer, optional:: level
  real(8):: now
  character(128):: s
  real(8), save:: last=0d0
  integer, save:: ilog=1
  !-----------------------------------------------------------------------------
  if (io%task_logging > 0) then
    if (present(level)) then
      if (level > io%task_logging) &
        return
    end if
    !$omp critical (tasklog_cr)
    now = wallclock()
    if (now-last > io%log_interval) then
      write (io_unit%task,'(f12.6,4x,a)') now, 'swapping task log file'
      close (io_unit%task)
      if (ilog==0) then
        call io%open_rank_file (io_unit%task, ".task1", 'formatted')
        ilog = 1
      else
        call io%open_rank_file (io_unit%task, ".task2", 'formatted')
        ilog = 0
      end if
      last = now
    end if
    !---------------------------------------------------------------------------
    ! Write a log line to io_uni%task, and echo to stdlog if verbose > 0
    !---------------------------------------------------------------------------
    if (os%verbose > 0) then
      write (s,'(f12.6,i5,i4,i7,":",i2.2,i6,2f12.6,1p,e11.2,4x,a)') &
        wallclock(), self%rank, omp%thread, self%id, self%level, &
        self%istep, os%gb, self%time, self%dtime, trim(label)
    else
      write (s,'(f12.6,i5,i4,i7,":",i2.2,i6,f12.6,1p,e11.2,4x,a)') &
        wallclock(), self%rank, omp%thread, self%id, self%level, &
        self%istep, self%time, self%dtime, trim(label)
    end if
    write (io_unit%task,'(a)') trim(s)
    if (verbose > 0) then
      write (stdlog,'(a)') 'task_t%log:'//trim(s)
      flush (stdlog)
    end if
    if (io_unit%flush_logs .or. os%isnan(self%time)) &
      flush (io_unit%task)
    !$omp end critical (tasklog_cr)
  end if
END SUBROUTINE log

!===============================================================================
!> Statistics collection stub; this call is made to only one task (the first
!> on the task list), allowing it to report solver statistics commont to, and
!> available to, all tasks
!===============================================================================
SUBROUTINE statistics (self)
  class(task_t):: self
  write (stdlog,*) 'task_t%statistics: stub'
END SUBROUTINE statistics

!===============================================================================
!> Stub for boundary condition implementations
!===============================================================================
SUBROUTINE bndries (self)
  class(task_t):: self
END SUBROUTINE bndries

SUBROUTINE periodic_grid (self, only)
  class(task_t):: self
  integer, optional:: only
  ! ----------------------------------------------------------------------------
  call io%abort('task_t%periodic grid called')
END SUBROUTINE periodic_grid

!===============================================================================
!> Make sure the time step is not longer than io%out_time, and set exactly
!> syncronized output time from all tasks, if possible
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(task_t):: self
  logical, optional:: detailed_timer
  real(8):: start_time, dtime
  integer, save:: nprint=10
  logical:: snapshot_time
  logical, parameter:: ramp_time=.false.
  !-----------------------------------------------------------------------------
  start_time = io%start_time - io%out_time
  call trace%begin ('task_t%courant_condition')
  snapshot_time = .false.
  dtime = self%dtime
  !if (io%do_output .and. io%out_time > 0.0 .and. self%time > start_time) then
  if (io%do_output .and. io%out_time > 0.0 .and. &
      self%is_set (bits%has_output)) then
    if (self%time < self%out_next) then
      if (self%time + dtime > self%out_next) then
        dtime = self%out_next - self%time
        if (dtime == 0d0) then
          dtime = 1d-10*self%time
          write (stdlog,'(a,2i6,1p,e115.6)') &
            'task_t%courant_condition: fix zero dtime =', &
            self%id, self%istep, dtime
        end if
        do while (self%time + dtime < self%out_next)
          dtime = dtime*1.00001d0
          write (stdlog,'(a,2i6,1p,e115.6)') &
            'task_t%courant_condition: iterating, dtime =', &
            self%id, self%istep, dtime
        end do
        !-----------------------------------------------------------------------
        ! If dtime is reduced with more than a factor of 2, set dt_extra, which
        ! may trigger the solver to take that extra step.
        !-----------------------------------------------------------------------
        if (dtime < self%dtime/2d0) then
          self%dt_extra = self%dtime
        end if
        !.......................................................................
        if (io%verbose > 0) &
          write (stdlog,'(a,2i6,1p,4e15.6)') &
              'task_t%courant_condition: snapshot dtime =', &
              self%id, self%istep, self%dt(self%it), dtime, &
              self%time, self%out_next
        !$omp atomic write
        !.......................................................................
        self%dtime = dtime
        snapshot_time = .true.
      end if
    end if
    if (dtime > io%out_time) then
      dtime = io%out_time
      !$omp atomic write
      self%dtime = dtime
    end if
  end if
  !---------------------------------------------------------------------------
  ! After a few sime steps, and with time slices:
  !---------------------------------------------------------------------------
  if (self%nt > 1 .and. self%istep >= self%nt .and. self%dt_fixed < 0.0) then
    !-------------------------------------------------------------------------
    ! Check for small timesteps
    !-------------------------------------------------------------------------
    if (dtime < self%dt_small .and. .not.snapshot_time) then
      if (nprint > 0) then
        print '(a,i6,1p,2g15.6)', &
          'task_t%courant_condition WARNING: id,dtime=', &
          self%id, dtime, self%time
        nprint = nprint-1
      else
        call io%abort ('task_t%courant_condition: dt < dt_small')
      end if
    !-------------------------------------------------------------------------
    ! Allow at most a factor two larger dt than in the previous step
    !-------------------------------------------------------------------------
    else if (ramp_time .and. dtime > 2d0*self%dt(self%it)) then
      if (io%verbose > 0) &
        write (stdlog,'(a,2i6,1p,2e15.6)') &
          'task_t%courant_condition: ramp dtime =', &
          self%id, self%istep, self%dt(self%it), dtime
      dtime = 2d0*self%dt(self%it)
      !$omp atomic write
      self%dtime = dtime
    !-------------------------------------------------------------------------
    !> Prevent too rapid increase of time steps after a sharp decrease
    !-------------------------------------------------------------------------
    else if (dtime < 0.5*self%dt(self%it) .and. .not.snapshot_time) then
      self%safety_factor = 0.5
      write (stdlog,'(a,2i6,1p,2e15.6)') &
            'task_t%courant_condition: setting safety, id,dt(it),dt=', &
            self%id, self%istep, self%dt(self%it), dtime
    end if
    !-------------------------------------------------------------------------
    ! Event memory
    !-------------------------------------------------------------------------
    if (self%safety_factor < 0.95) then
      self%safety_factor = self%safety_factor**0.9
      dtime = dtime*self%safety_factor
      !$omp atomic write
      self%dtime = dtime
    else
      self%safety_factor = 1.0
    end if
  end if
  call trace%end()
END SUBROUTINE courant_condition

!===============================================================================
!> Template for output
!===============================================================================
SUBROUTINE output (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%output')
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Template for metadata output
!===============================================================================
SUBROUTINE output_metadata (self, unit, first_time)
  class(task_t):: self
  integer:: unit
  logical, optional:: first_time
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%output_metadata')
  call io%abort ('task_t%output_metadata called')
  call trace%end ()
END SUBROUTINE output_metadata

!===============================================================================
!> Template for metadata output
!===============================================================================
SUBROUTINE input_metadata (self, unit)
  class(task_t):: self
  integer:: unit
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%input_metadata')
  call io%abort ('task_t%input_metadata called')
  call trace%end ()
END SUBROUTINE input_metadata

!===============================================================================
SUBROUTINE allocate_mesg (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  call io%abort ('task_t%allocate_mesg called')
END SUBROUTINE allocate_mesg

!===============================================================================
!> Return region of interest, adjusted for distance
!===============================================================================
FUNCTION roi (self, task) RESULT (region)
  class(task_t)           :: self
  class(task_t), optional :: task
  type(region_t)          :: region
  real(8)                 :: distance(3)
  !-----------------------------------------------------------------------------
  if (present(task)) then
    distance = task%distance_to (self)
    region%l = distance + self%interest%l
    region%u = distance + self%interest%u
  else
    region%l = self%interest%l
    region%u = self%interest%u
  end if
END FUNCTION roi

!===============================================================================
!> Return region of authority, adjusted for distance
!===============================================================================
FUNCTION roa (self, task) RESULT (region)
  class(task_t)           :: self
  class(task_t), optional :: task
  type(region_t)          :: region
  real(8)                 :: distance(3)
  !-----------------------------------------------------------------------------
  if (present(task)) then
    distance = task%distance_to (self)
    region%l = distance + self%authority%l
    region%u = distance + self%authority%u
  else
    region%l = self%authority%l
    region%u = self%authority%u
  end if
END FUNCTION roa

SUBROUTINE update_position (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  call io%abort ("task_t%update_position")
END SUBROUTINE update_position

SUBROUTINE get_position (self, source, p, v, e)
  class(task_t):: self, source
  real(8):: time, p(3), v(3), e(3,3)
  optional:: p, v, e
  !-----------------------------------------------------------------------------
  call io%abort ("task_t%get_position_at")
END SUBROUTINE get_position

!===============================================================================
!> Return the velocity shift from source to target
!===============================================================================
FUNCTION velocity_shift (self, source, r_tgt) RESULT (v)
  class(task_t):: self, source
  real(8), dimension(3):: r_tgt, v
  !-----------------------------------------------------------------------------
  v = 0d0
  call io%abort ('task_t%velocity_shift')
END FUNCTION velocity_shift

!===============================================================================
!> Return a unique hash key for the patch, encoding position and level
!===============================================================================
FUNCTION hash_key (self, position, size, level, box, origin) RESULT (key)
  class(task_t):: self
  real(8), dimension(3), optional:: position, origin, size, box
  integer, optional:: level
  real(8), dimension(3):: p, o, s, b, c, rkey
  integer:: key(4), n(3), l
  !-----------------------------------------------------------------------------
  p = shared%option (position, self%position)
  s = shared%option (size    , self%size    )
  o = shared%option (origin  , self%origin  )
  b = shared%option (box     , self%box     )
  l = shared%option (level   , self%level   )
  n = nint(b/s)
  key(1) = l
  !-----------------------------------------------------------------------------
  ! A position half the box size from the origin (left corner) should have
  ! position index 0,0,0, with 1,1,1 if half a patch size displaced.
  ! The position is always encoded as the key times half the patch size.
  !-----------------------------------------------------------------------------
  c = o + b/2d0                         ! box center
  rkey = (p-c)/s                        ! distance from center in size units
  where (mod(n,2)==1)                   ! for odd number of patches
    key(2:4) = 2*nint(rkey)             ! round to zero, with increment 2
  else where
    key(2:4) = 2*floor(rkey) + 1        ! round to one, with increment 2
  end where
  where (self%periodic)
    key(2:4) = modulo (key(2:4)+n, 2*n) - n
  end where
  if (verbose > 0) &
    write (stdlog,'(a,i6,":",i2.2,3i4,4(2x,3f10.4),i4)') &
      'hash_key: id,key,pos,size,ori,box,level=', &
      self%id, self%key, self%position, self%size, self%origin, self%box, &
      self%level
END FUNCTION hash_key

!===============================================================================
!> Return the exact position corresponding to a hash key
!===============================================================================
FUNCTION hash_position (self, key) RESULT (out)
  class(task_t):: self
  integer:: key(4)
  real(8):: out(3)
  !.............................................................................
  out = self%origin + self%box/2d0 + key(2:4)*self%size/2d0 
END FUNCTION hash_position

!===============================================================================
!> Set the properties from the argument list
!===============================================================================
SUBROUTINE set_properties (self, key, position, size, level)
  class(task_t):: self
  integer, optional:: key(4), level
  real(8), optional:: size(3), position(3)
  !.............................................................................
  if (present(level)) then
    self%level = level
  end if
  if (present(size)) then
    self%size = size
  end if
  if (present(key)) then
    self%level = key(1)
    self%position = self%origin + self%size/2d0 + key(2:4)*self%size/2d0 
  end if
  if (present(position)) then
    self%position = position
  end if
END SUBROUTINE set_properties

!===============================================================================
!> Interface to refinement
!===============================================================================
INTEGER function check_refine (self, test)
  class(task_t), target:: self
  real, dimension(:,:,:), pointer, contiguous:: test
  !.............................................................................
  check_refine = -1                                                    ! keep
END FUNCTION check_refine

!===============================================================================
!> Write metadata namelist
!===============================================================================
SUBROUTINE write_nml (self, unit)
  class(task_t):: self
  integer:: unit
  integer:: id, id1, level, istep, rank
  real(8):: position(3), size(3), dtime, ds(3), velocity(3), cost
  real:: quality
  character(len=16):: kind
  namelist /task_nml/ istep, position, size, ds, velocity, dtime, &
    level, quality, rank, cost, kind
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%write_nml')
  id        = self%id
  istep     = self%istep
  position  = self%position
  size      = self%size
  ds        = self%ds
  velocity  = self%velocity
  level     = self%level
  quality   = self%quality
  dtime     = self%dtime
  rank      = self%rank
  cost      = self%cost
  kind      = self%kind
  !-----------------------------------------------------------------------------
  ! Work-aound, to fix a problem with gfortran namelist output giving ID=0
  !-----------------------------------------------------------------------------
  write (unit,'(a)') '&TASK_NML'
  write (unit,'(" ID=",i6,",")') self%id
  write (unit,'(a)') '/'
  write (unit, task_nml)
  flush (unit)
  call trace%end ()
END SUBROUTINE write_nml

!===============================================================================
!> Reset all pointers (except %name) in cloned task_t instances
!===============================================================================
SUBROUTINE reset (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%reset')
  nullify (self%parent)
  nullify (self%mesg)
  nullify (self%sent_list)
  nullify (self%iit)
  nullify (self%t)
  nullify (self%dt)
  nullify (self%lock)
  nullify (self%mulligan_lock)
  call self%clear (bits%ready + bits%busy)
  call trace%end ()
END SUBROUTINE reset

!===============================================================================
!> Protected periodic wrapping
!===============================================================================
SUBROUTINE periodic_wrap4 (self, position)
  class(task_t):: self
  real(4):: position(3)
  where (self%periodic .and. self%box > 0d0)
    position = modulo (position + self%box*0.5d0, self%box) - self%box*0.5d0
  end where
END SUBROUTINE periodic_wrap4
!===============================================================================
SUBROUTINE periodic_wrap8 (self, position)
  class(task_t):: self
  real(8):: position(3)
  where (self%periodic .and. self%box > 0d0)
    position = modulo (position + self%box*0.5d0, self%box) - self%box*0.5d0
  end where
END SUBROUTINE periodic_wrap8

!===============================================================================
!> Stub
!===============================================================================
SUBROUTINE set_rank (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%set_rank')
  call trace%end ()
END SUBROUTINE set_rank

END MODULE task_mod
