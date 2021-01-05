!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE io_mod
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
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type, public:: task_t
    integer:: id = 0
    integer:: status = 0
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
    integer:: thread=-1
    integer:: mem_thread=-1
    integer:: n_dump=0                  ! number of dump files
    integer:: ip                        ! rank local task id
    integer:: n_needed=1                ! number of tasks that needs it
    integer:: seq=0                     ! incremental sequence number
    integer:: n_failed                  ! number of times check_ready() failed
    integer:: logging=0                 ! level of logging
    integer:: next_lb                   ! next LB step
    integer:: rank_steps=0              ! steps on rank
    real(8):: time = 0d0
    real(8):: atime = 0d0               ! for active queue
    real(8):: out_next = 0d0
    real(8):: print_next = 0d0
    real(8):: dtime = 0d0
    real(8):: min_dtime = 0d0
    real(4):: grace = 0.0
    real(8):: position(3) = 0d0
    real(8):: origin(3)=0d0             ! the cartesian box origin (for BCs)
    real(8):: velocity(3) = 0d0         ! velocity of task in Cartesian coords.
    real(8):: size(3) = 0.9999d0
    real(8):: gsize(3)                  ! convenience
    real(8):: ds(3)=0d0                 ! convenience
    real(8):: box(3)=1d0                ! the periodic wrapping size
    real(8):: wallclock
    real(8):: unpack_time=0d0
    real(8):: dnload_time=-1d0          ! last dnload time
    real(8):: mulligan_time=-1d0        ! time to rewind to
    real(4):: quality = 0.0
    character(len=64):: kind='task'     ! solver name
    character(len=12):: type='task_t'   ! task type
    character(len=64),     pointer:: name      => null()
    class(task_t),         pointer:: parent    => null()
    class(mesg_t),         pointer:: mesg      => null()
    class(mesg_list_t),    pointer:: sent_list => null()
    integer, dimension(:), pointer:: iit       => null()
    real(8), dimension(:), pointer:: t         => null()
    real(8), dimension(:), pointer:: dt        => null()
    type(lock_t),          pointer:: lock      => null()
    type(lock_t),          pointer:: mulligan_lock => null()
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
    type(random_t):: random             ! pseudo-random number generator
    type(aux_t):: aux
    integer(8):: amr_offset=0_8         ! byte offset in amr_io
  contains
    procedure:: init
    procedure:: init_unique
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
    procedure:: allocate_mesg => void
    procedure:: unique_id
    procedure:: load_balance
    procedure:: has_finished
    procedure:: overlaps
    procedure:: overlaps_point
    procedure:: nbor_relations
    procedure:: distance_from
    procedure:: distance
    procedure:: point_distance
    procedure:: pack
    procedure:: unpack
    procedure:: is_relevant
    procedure:: test_update
    procedure:: solver_is
    procedure:: debug
    procedure, private:: contains1      ! patch contains a point
    procedure, private:: contains2      ! patch contains a patch
    generic:: contains => contains1, contains2
    procedure:: task_info
    procedure:: log
    procedure:: run => void
    procedure:: statistics
    procedure:: bndries
    procedure:: courant_condition
    procedure:: output_metadata
    procedure:: pack_output
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
  end type
  logical:: first_time=.true.
  integer:: track_id=0
  integer:: verbose=0
  integer, dimension(:), pointer, save:: id => null()
  type(dll_t), pointer:: deleted_list
CONTAINS

!===============================================================================
SUBROUTINE void (self)
  class (task_t):: self
  call trace%begin ('task_t%void')
  print*, 'ERROR: task_t%void called'
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
!> Define unique task id, unless already set before entry.
!===============================================================================
SUBROUTINE init (self)
  class (task_t):: self
  integer:: id, id_master
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%init')
  call init_params (self)
  if (self%id == 0) then
    call init_unique (self)
    if (verbose > 0) &
      write (io_unit%log,*) 'task_t%init: id =', id
  else
    !---------------------------------------------------------------------------
    ! If a task already has an id, update the MPI counter if its value is less.
    ! The id of a virtual copy of a new AMR task should NOT trigger this.
    !---------------------------------------------------------------------------
    id_master = mpi%id%update(0)
    if (verbose > 0) &
      write (io_unit%log,*) 'task_t%init: id, id_master =', self%id, id_master
    if (id_master < self%id) then
      write (io_unit%log,*) 'task_t%init: updating master id', self%id, id_master
      id_master = mpi%id%update (self%id-id_master)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Locks, random numbers, etc
  !-----------------------------------------------------------------------------
  if (.not.associated(self%lock)) then
    allocate (self%lock, self%mulligan_lock)
    call self%lock%init (self%kind(1:4))
    call self%mulligan_lock%init ('mull')
  end if
  if (.not.self%lock%initialized) then
    call self%lock%init (self%kind(1:4))
  end if
  call self%random%init
  call self%aux%init
  !-----------------------------------------------------------------------------
  ! Only initialize time slot values when first allocated
  !-----------------------------------------------------------------------------
  if (.not.associated(self%t)) then
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))
    self%dt = -1d0                                ! mark as "not known"
    self%t  = 0d0                                 ! initial time
    self%iit = 1                                  ! initial slots
    self%it = 1                                   ! slot for current values
    self%new = 2                                  ! slot for new values
    self%iit(self%nt) = self%new                  ! slot for new values
    self%time = self%t(self%it)                   ! initial time
  end if
  if (self%id == track_id) then
    self%track = .true.
  end if
  call self%log ('init', 2)
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Define unique task id, unless already set before entry.
!===============================================================================
SUBROUTINE init_params (self)
  class (task_t):: self
  character(len=120):: ids= &
    '$Id$'
  integer:: iostat, id, id_master
  logical, save:: do_dumps
  integer, save:: id_max=0
  namelist /task_params/ track_id, do_dumps, verbose
  namelist /dump_params/ id, iostat
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
  if (do_dumps) then
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
        write (io_unit%log,*) 'task_t%init: dump_id =', id
        write (io%namelist, dump_params)
        exit
      end if
    end do
  end if
  call io%lock%unset ('task_t%init_params')
  call trace%end ()
END SUBROUTINE init_params

!===============================================================================
!> Define unique patch ids.  A critical region is OK; this is only done once.
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
      write (io_unit%log,*) 'task_t%init_unique(same): id =', &
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
      write (io_unit%log,*) 'task_t%init_unique(else): id =', &
        id, id_same, id_master
  end if
  self%id = id
  call trace%end ()
END SUBROUTINE init_unique

!===============================================================================
!> Deallocate permanent arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class (task_t):: self
  !.............................................................................
  call trace%begin ('task_t%dealloc')
  if (associated(self%t  )) deallocate (self%t  )
  if (associated(self%dt )) deallocate (self%dt )
  if (associated(self%iit)) deallocate (self%iit)
  if (associated(self%mesg)) then
    call self%mesg%dealloc
    deallocate (self%mesg)
  end if
  if (associated(self%sent_list)) then
    call self%sent_list%dealloc
    deallocate (self%sent_list)
  end if
  call trace%end
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
  call trace%end
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
    write (io_unit%log,*) omp%thread,'header_to_task WARNING: nt non-positive', nt
    flush (io_unit%log)
  end if
  if (head%istep == 0) then
    write (io_unit%log,'(1x,a,i6,f12.6,2x,9f12.6)') &
      'task%header_to_task: id, time, t =', &
      self%id, self%time, self%t
    write (io_unit%log,'(1x,a,i6,2x,3i4,2x,9i4)') &
      'task%header_to_task: id, istep, it, new, iit =', &
      self%id, self%istep, self%it, self%new, self%iit
  end if
  if (self%it > self%nt .or. self%it < 1) then
    write (io_unit%log,*) 'header_to_task ERROR: it =', self%it, self%iit(nt-1)
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
    write (io_unit%log,'(a,7i10)')   'task_t%rotate iit =', self%iit
    flush (io_unit%log)
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
                                                 call self%lock%set ('rotate')
  self%dt(self%it)= self%dtime                          ! save dt
  self%dt(self%new)= self%dtime                         ! estimate
  self%t(self%new) = time                               ! initial time
  self%it = self%new                                    ! update time slot, new
  self%new = mod(self%new,self%nt) + 1                  ! circular increment
  self%iit(1:self%nt-1) = self%iit(2:self%nt)           ! rotate left
  self%iit(self%nt) = self%new                          ! top slot
  self%istep = self%istep + 1                           ! increment task step
  self%rank_steps = self%rank_steps + 1                 ! increment rank step
  self%time = time                                      ! reveal new time
  self%rotated = .true.
                                                 call self%lock%unset ('rotate')
  !.............................................................................
  if (io%verbose>2) then
   write (io_unit%log,*) 'task id:',self%id,'  it:',self%it,'  new:',self%new
   if (io%verbose>3) then
     write (io_unit%log,*) 'task iit:', self%iit
     write (io_unit%log,*) 'task t(iit):', self%t(self%iit)
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
  class (task_t):: self, target
  real(8):: nbtime, nbdtime, tgtime
  !.............................................................................
  if (omp_lock%tasks) then
                                            call self%lock%set ('is_ahead_of')
    nbtime = self%time
    nbdtime = self%dtime
    tgtime = target%time
    if (self%is_set(bits%frozen)) then
      is_ahead_of = .true.
    else if (self%istep < 3) then
      is_ahead_of = nbtime >= tgtime
    else
      is_ahead_of = nbtime + nbdtime*self%grace > tgtime
    end if
                                            call self%lock%unset ('is_ahead_of')
  else
    !$omp atomic read
    nbtime = self%time
    !$omp atomic read
    nbdtime = self%dtime
    !$omp atomic read
    tgtime = target%time
    ! -- if a task is frozen (in time) assume it is forever ahead of other tasks
    if (self%is_set(bits%frozen)) then
      is_ahead_of = .true.
    ! -- a negative nbor time signals an invalid (e.g. initial RT) state
  !  else if (nbtime < 0.0) then
  !    is_ahead_of = .false.
    ! -- do not use a grace interval for the first few steps
    else if (self%istep < 3) then
      is_ahead_of = nbtime >= tgtime
    ! -- use a grace interval that is fraction of the nbor time step
    else
      is_ahead_of = nbtime + nbdtime*self%grace > tgtime
    end if
  end if
  if (io_unit%verbose>1) then
    if (target%id==io%id_debug.or.io_unit%verbose>4) &
    print'(i6,i4,2x,a,i6,3f9.5,l3)', self%id, omp_mythread, 'mk is_ahead_of: ', &
    target%id, self%time, self%dtime*target%grace, target%time, is_ahead_of
  end if
END FUNCTION

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
FUNCTION unique_id (self)
  class (task_t):: self
  integer:: unique_id
  !$omp atomic
  id(self%rank) = id(self%rank) + mpi%size
  unique_id = id(self%rank)
END FUNCTION unique_id

!===============================================================================
!> The load balancing procedure is part of the negotiation between a patch and
!> its virtual neighbors.  We use the pro-active approach, where the decision
!> to migrate to a different rank is taken by the task itself.  It can then be
!> implemented with a minimum of complication, in that one of the last things
!> the task does before broadcasting itself to its neighbors is to change its
!> rank.  The new rank should check if an incoming virtual patch in fact has
!> opted to become a boundart patch instead.  It comes with a complete list of
!> neighbors already.
!>
!> To decide whether to migrate to another rank, the following parameters should
!> be considered, for both tasks:
!>
!> task%load    : the load of the task, in milliseconds per update
!> task%dt      : the time step per update
!> task%time    : the current time of the task
!> rank%time    : the current time of the neighbor rank
!> rank%surplus : the extra capacity available on the nbor rank
!>
!> If the time difference (rank%time-task%time)/task%dt is positive it means
!> that even the most delayed task on the other rank is ahead of the current
!> task, which speaks in favor of migrating.
!>
!> If the nbor rank has extra capacity, this also speaks in favor of migration
!>
!> A good measure of the load of a rank is the number of wall clock seconds per
!> simulation time unit.  This can be computed by running through the list of
!> tasks on the rank, adding for each task its cost in ms/time_unit, and finally
!> dividing by the number of (real) cores available to work on the tasks.
!===============================================================================
SUBROUTINE load_balance (self)
class(task_t):: self

END SUBROUTINE load_balance

!===============================================================================
!> Check if there is overlap, using the size plus an extra cell margin. The
!> definition of overlap MUST be symmetric, since it is being used to define
!> neighbor relations.
!>
!> NOTE 1: This function was in the past NOT used by e.g. list_t%init_nbors, which
!> was instead using list_t%overlaps, based on the same criterion as below.
!>
!> NOTE 2: Since overlap and nbors can have a generalized definition it makes
!> sense that is valid for the whole task_t class.  If needed, the definition
!> can be overloaded with alternative definitions for particular tasks.
!===============================================================================
FUNCTION overlaps (self, task)
  class(task_t):: self
  class(task_t):: task
  logical:: overlaps
  real(8):: dr(3)
  integer:: itimer=0
  !.............................................................................
  if (timer%detailed) then
    call trace%begin ('task_t%overlaps', itimer=itimer)
    dr = distance (self, task)
    overlaps = all(abs(dr) <= 0.5_8*(self%size+task%size+self%ds+task%ds))
    call trace%end (itimer)
  else
    overlaps = all(abs(distance (self, task)) <= &
                   0.5_8*(self%size+task%size+self%ds+task%ds))
  end if
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
!> on ranks and some bits.
!> Initially all 3 logicals are set to true and if any criterion is not met,
!> the particular flag is set to .false.
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(task_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%begin ('task_t%nbor_relations', 2, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! task is real
  !-----------------------------------------------------------------------------
  if (self%rank == mpi%rank) then
    call self%clear(bits%external)
    !---------------------------------------------------------------------------
    ! nbor is frozen
    !---------------------------------------------------------------------------
    if (another%is_set (bits%frozen)) then
      needs_me = .false.
      needed = .false.
    !---------------------------------------------------------------------------
    ! nbor is virtual
    !---------------------------------------------------------------------------
    else if (another%rank /= mpi%rank) then
      call self%set(bits%boundary)
      call self%clear(bits%internal)
      ! -- use in check_ready(), but not in check_nbors()
      needed = .true.
      needs_me = .false.
    !---------------------------------------------------------------------------
    ! nbor is real
    !---------------------------------------------------------------------------
    else
      ! -- use in check_ready(), and in check_nbors()
      needed = .true.
      needs_me = .true.
    end if
    !---------------------------------------------------------------------------
    ! real tasks need guard zones
    !---------------------------------------------------------------------------
    download = .true.
  !-----------------------------------------------------------------------------
  ! task is virtual
  !-----------------------------------------------------------------------------
  else
    call self%clear(bits%internal)
    !---------------------------------------------------------------------------
    ! nbor is frozen
    !---------------------------------------------------------------------------
    if (another%is_set (bits%frozen)) then
      needs_me = .false.
      needed = .false.
    !---------------------------------------------------------------------------
    ! task is virtual, nbor is real
    !---------------------------------------------------------------------------
    else if (another%rank==mpi%rank) then
      call self%set(bits%virtual)
      call self%clear(bits%external)
      ! -- don't use in check_ready(), use in check_nbors()
      needed = .false.
      needs_me = .true.
    !---------------------------------------------------------------------------
    ! task is virtual, nbor is virtual
    !---------------------------------------------------------------------------
    else
      needed = .false.
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
  !  needed = .false.
  !  needs_me = .false.
    download = .false.
  end if
  if (self%verbose > 2) &
    write(io_unit%log,'(a,i6,i4,3x,3l1)') &
      'set_nbor_relations: nbor, rank, IBV =', another%id, another%rank, &
      another%is_set(bits%internal), &
      another%is_set(bits%boundary), &
      another%is_set(bits%virtual)
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%end (itimer)
END SUBROUTINE nbor_relations

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

SUBROUTINE unpack (self, mesg)
  class(task_t):: self
  class(mesg_t), pointer:: mesg
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
SUBROUTINE task_info (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(1x,a,40x,i6,2i4,2x,5x,2l1,2x,3f9.3)') &
    'patch_t%task_info: id, BV, pos =', self%id, self%rank, self%level, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
END SUBROUTINE task_info

!===============================================================================
SUBROUTINE log (self, label, level)
  class(task_t):: self
  character(len=*):: label
  integer, optional:: level
  real(8):: now
  real(8), save:: last=0d0
  integer, save:: ilog=0
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
    write (io_unit%task,'(f12.6,i5,i4,i6,i4,f12.6,1p,e11.2,4x,a)') &
      wallclock(), mpi%rank, omp%thread, self%id, self%level, &
      self%time, self%dtime, trim(label)
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

!===============================================================================
!> Make sure the time step is not longer than io%out_time, and set exactly
!> syncronized output time from all tasks, if possible
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(task_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%courant_condition')
  if (io%do_output .and. io%out_time > 0.0) then
    self%dtime = min (self%dtime, io%out_time)
    if (self%time < self%out_next) then
      if (self%time + self%dtime*(1.+self%grace) > self%out_next) then
        self%dtime = self%out_next - self%time
        do while (self%time + self%dtime*(1.+self%grace) < self%out_next)
          self%dtime = self%dtime*1.00001d0
          write (stdlog,'(a,i6,1p,e115.6)') &
            'task_t%courant_condition: iterating, dtime =', self%id, self%dtime
        end do
      end if
    end if
  end if
  call trace%end()
END SUBROUTINE courant_condition

!===============================================================================
!> Template for metadata output
!===============================================================================
SUBROUTINE output_metadata (self)
  class(task_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%write_metadata')
  call trace%end ()
END SUBROUTINE output_metadata

!===============================================================================
!> Template for self-packing routine
!===============================================================================
SUBROUTINE pack_output (self, p)
  class(task_t):: self
  class(*), pointer:: p
  !-----------------------------------------------------------------------------
  call trace%begin ('task_t%pack_output')
  nullify(p)
  call trace%end ()
END SUBROUTINE pack_output

END MODULE task_mod
