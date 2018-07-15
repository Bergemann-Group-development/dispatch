!> =============================================================================
!> Execute a task list.  Since only the master thread is calling check_mpi(), we
!> need to make sure that it handles all send/recv, hence send_priv=.false. and
!> recv_priv=.false., and we can use the simplest receive mechanism (recv_active
!> =.true.), with no buffering needed (queue_unpack=.false.).  These values are
!> imposed after the input namelists are read, to make sure the choices are
!> correct.
!>
!> Strategy: Keep track of, separately, the number of ready tasks, the number
!> of currently spawned tasks, and the number of currently busy tasks.  To this
!> end, define the following states:
!>
!> state=0: task is inactive or new
!> state=1: the task is ready to update, but not yet spawned
!> state=2: the task has been spawned
!> state=3: the task is busy updating
!> state=4: the task has finished updating
!>
!> The simplest strategy for counting is to keep track of the number of tasks
!> in each state, by incrementing the count for the new state and decrementing
!> the count for the old state, whenever the state changes. These counts should
!> be part of the dispatcher data type.
!>
!> On the balance between spawning tasks for other treads and actually working
!> on an updete:  If a process has N threads, and has N tasks active, evenly
!> spread in their progress, then essentially all of these tasks will have
!> finished updating while the master thread works on its update.  Hence there
!> must be approximately N additional tasks spawned, but not yet busy updating.
!> Hence, the criterion for the master thread to take on work itself is that
!> there are at least N tasks in state 2 (not counting those in state 3 or 4).
!>
!> If the master thread is unable to find enough tasks in state 1 for it to be
!> able to reach N in state 2, then it is certainly not a good idea to take on
!> one of the tasks, since there is alreadt a surplus of threads available.
!> =============================================================================
MODULE dispatcher3_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE timer_mod
  USE link_mod
  USE mpi_mesg_mod
  USE task_list_mod
  USE bits_mod
  USE global_mod
  USE task_mod
  implicit none
  private
  type, public:: dispatcher3_t
    integer:: verbose=0
    integer:: n_spawn=0
    logical:: use_locks=.false.
    logical:: use_critical=.true.
    type(task_list_t), pointer:: task_list => null()
  contains
    procedure:: execute
    procedure:: update
    procedure:: spawn_update
    procedure:: set_state
    procedure:: is_ready
    procedure:: print_state
  end type
  type(global_t):: average_density
  integer:: n_fail=0
  integer:: n_state(0:4)=0
  real(8):: wc_failed=0d0
  real(8):: start
  character(len=32):: fmt0='(    f12.6,i4,2x,a,2x,i7,2i4)'
  character(len=32):: fmt1='(40x,f12.6,i4,2x,a,2x,i7,2i4)'
  type(dispatcher3_t), public:: dispatcher3
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self, task_list)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  !.............................................................................
END SUBROUTINE init

!> =============================================================================
!> Simple dispatcher, which in each loop through the task list takes on the
!> first task that is ready, and then spawns other threads that handle the
!> remaining tasks that are ready.
!> =============================================================================
SUBROUTINE execute (self, task_list, use_locks, use_critical, test, verbose)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  integer:: verbose
  logical:: use_locks, use_critical, test, fell_through
  !.............................................................................
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  logical:: mytask
  real(8):: start_fail=0d0
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  self%task_list    => task_list
  self%use_locks    = use_locks
  self%use_critical = use_critical
  self%verbose      = verbose
  call task_list%init_all_nbors
  n_state = 0; n_state(0) = io%ntask
  call tic (time=start)
  call timer%print
  !$omp parallel private(link,nbor,task,mytask) shared(n_state) default(shared)
  !$omp master
  print *, 'DISPATCHER3', n_state
  call trace%begin ('dispatcher3_t%execute', itimer=itimer)
  do while (task_list%na>0 .and. wallclock()-start < io%job_seconds)
    !---------------------------------------------------------------------------
    ! Check if there is MPI work
    !---------------------------------------------------------------------------
    call task_list%check_mpi()
    !---------------------------------------------------------------------------
    ! Run through the task list, looking for and counting tasks that are ready
    !---------------------------------------------------------------------------
    link => task_list%head
    do while (associated(link))
      if (.not.link%task%is_set (bits%virtual)) then
        !-----------------------------------------------------------------------
        ! If task is updated, set it to idle; this should ONLY be done by the
        ! master thread
        !-----------------------------------------------------------------------
        if (link%task%state==4) then
          call self%set_state (link%task, 0)
        end if
        !-----------------------------------------------------------------------
        ! If task is idle, check if it is ready to update
        !-----------------------------------------------------------------------
        if (link%task%state==0) then
          if (self%is_ready (link)) then
            call self%set_state (link%task, 1)
          end if
        end if
      end if
      link => link%next
    end do
    if (self%verbose>1) call self%print_state ('1')
    !---------------------------------------------------------------------------
    ! Run through the task list and start updates of tasks that are ready
    !---------------------------------------------------------------------------
    fell_through = .true.
    link => task_list%head
    do while (associated(link))
      task => link%task
      if (.not.task%is_set (bits%virtual)) then
        !-----------------------------------------------------------------------
        ! Spawn a background task, unless enough has been spawned, or there is
        ! only one left
        !-----------------------------------------------------------------------
        if (task%state==1) then
          !$omp critical (set_state_cr)
          mytask = omp%nthreads==1 .or. n_state(1)==1 .or. n_state(2)>omp%nthreads
          !$omp end critical (set_state_cr)
          if (mytask) then
            call self%set_state (task, 2)
            call self%update (task_list, link, test)
          else
            call self%spawn_update (task_list, link, test)
          end if
          fell_through = .false.
        end if
      end if
      link => link%next
    end do
    if (fell_through) then
      if (self%verbose>3) print *,'task yield'
      n_fail = n_fail+1
      if (n_fail==1) start_fail=wallclock()
      !!omp taskyield
    else
      if (n_fail>0) then
        wc_failed = wc_failed + (wallclock()-start_fail)
        if (self%verbose>0) call self%print_state ('2')
      end if
      n_fail = 0
    end if
  end do
  call trace%end (itimer)
  !$omp end master
  !$omp end parallel
  call timer%print
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', n_update, time=start)
END SUBROUTINE execute

!> =============================================================================
!> =============================================================================
SUBROUTINE print_state (self, label)
  class(dispatcher3_t):: self
  character(len=*):: label
  !$omp critical (set_state_cr)
  print '("DISPATCHER4: ",a,f12.6,f9.3,4x,a,i8,6(4x,a,i4))', &
    label, wallclock(), &
    wc_failed/max(1d-10,wallclock()-start), &
    'failed:'  , n_fail, &
    'idle:'    , n_state(0), &
    'ready:'   , n_state(1), &
    'spawned:' , n_state(2), &
    'updating:', n_state(3), &
    'updated:' , n_state(4)
  !$omp end critical (set_state_cr)
END SUBROUTINE

!> =============================================================================
!> A task that is being acticely updated is set in state 3, and when finished
!> it switched to state 4.   We then use the same OMP task to check if any if
!> the nbors have beome ready to update, switching the state if it has, and
!> leaving it to the master thread to actually activate the task.
!> =============================================================================
SUBROUTINE update (self, task_list, link, test)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link, nbor
  logical:: test
  !...........................................................................
  call self%set_state (link%task, 3)
  link%task%thread = omp%thread
  if (self%verbose>3) print fmt1,wallclock(),omp%thread,'updating task',link%task%id
  call task_list%update (link, test)
  if (self%verbose>3) print fmt1,wallclock(),omp%thread,'updated  task',link%task%id
  nbor => link%nbor
  do while (associated(nbor))
    if (self%is_ready(nbor%link)) then
      call self%set_state (nbor%link%task, 1)
      if (self%verbose>3) print fmt1,wallclock(),omp%thread,'found(2) task',nbor%task%id
    end if
    nbor => nbor%next
  end do
  call self%set_state (link%task, 4)
  link%task%thread = -1
END SUBROUTINE

!> =============================================================================
!> Spawn an OpenMP task, for updating a DISPATCH task
!> =============================================================================
SUBROUTINE spawn_update (self, task_list, link, test)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link
  logical:: test
  !...........................................................................
  call self%set_state (link%task, 2)
  if (self%verbose>2) print fmt0,wallclock(),omp%thread,'spawning task',link%task%id, link%task%state, link%task%thread
  !$omp task default(none) firstprivate(self,link,test,task_list)
    if (self%use_locks) then
      if (link%task%lock%test()) then
        call self%update (task_list, link, test)
        call link%task%lock%unset
      else
        print *,'WARNING, lock conflict, thread, id =', omp%thread, link%task%id, link%task%thread
      end if
    else
      call self%update (task_list, link, test)
    end if
  !$omp end task
END SUBROUTINE

!> =============================================================================
!> Set a new task state, while keeping track of the number of tasks in each state
!> =============================================================================
SUBROUTINE set_state (self, task, in)
  class(dispatcher3_t):: self
  class(task_t), pointer:: task
  integer:: in, io
  !...........................................................................
  if (self%use_critical) then
    !$omp critical (set_state_cr)
    io = task%state
    if (self%verbose>1) &
      print '(a,i4,f12.6,2x,5i4,2x,2i3)', ' dispatcher3:', omp%thread, wallclock(), n_state, io, in
    n_state(io) = n_state(io) - 1
    n_state(in) = n_state(in) + 1
    task%state = in
    !$omp end critical (set_state_cr)
  else
   !$omp atomic read
   io = task%state
   !$omp atomic
   n_state(io) = n_state(io) - 1
   !$omp atomic
   n_state(in) = n_state(in) + 1
   !$omp atomic write
   task%state = in
   !print *,'dispatcher3:',omp%thread,wallclock(),n_state
  end if
END SUBROUTINE

!> =============================================================================
!> Check if a task is ready to update, by changing if all nbor task are "ahead
!> of" the task in question
!> =============================================================================
LOGICAL FUNCTION is_ready (self, link)
  class(dispatcher3_t):: self
  class(link_t), pointer:: link, nbor
  integer:: state
  !...........................................................................
  !$omp atomic read
  state = link%task%state
  if (link%task%is_set (bits%virtual)) then
    is_ready = .false.
  else if (state>=1) then
    is_ready = .false.
  else
    is_ready = .true.
    nbor => link%nbor
    do while (associated(nbor))
      !$omp atomic read
      state = nbor%task%state
      if (state==3) then
        is_ready = .false.
      else
        is_ready = is_ready .and. nbor%task%is_ahead_of (link%task)
      end if
      if (.not.is_ready) then
        exit
      end if
      nbor => nbor%next
    end do
  end if
END FUNCTION

END MODULE dispatcher3_mod
