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
!> one of the tasks, since there is already a surplus of threads available.
!>
!> If, on the other hand, there is only single task available in state 2, then
!> it IS a good idea to let the master thread update it, since it then has some
!> useful work to do, while other tasks are finishing, so when the master thread
!> looks at the situation again, there are likely to be a number of new tasks
!> ready to be updated.
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
    type(task_list_t), pointer:: task_list => null()
  contains
    procedure:: init
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
  integer:: verbose=0
  integer(8):: n_check=0, n_ready=0
  real:: f_over=2.0
  real(8):: wc_failed=0d0
  real(8):: start
  logical:: use_locks=.false.
  logical:: use_critical=.true.
  logical:: use_taskyield=.false.
  character(len=32):: fmt0='(    f12.6,i4,2x,a,2x,i7,2i4)'
  character(len=32):: fmt1='(40x,f12.6,i4,2x,a,2x,i7,2i4)'
  character(len=48):: fmt2='(40x,f12.6,i4,2x,a,2x,i7,2x,a,1p,g14.5)'
  type(dispatcher3_t), public:: dispatcher3
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self, task_list)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  integer:: iostat
  namelist /dispatcher3_params/ use_locks, use_critical, use_taskyield, f_over, &
    verbose
  !.............................................................................
  rewind (io%input); read (io%input, dispatcher3_params, iostat=iostat)
  if (io%master) write (io%output, dispatcher3_params)
  self%verbose = verbose
END SUBROUTINE init

!> =============================================================================
!> Simple dispatcher, which in each loop through the task list takes on the
!> first task that is ready, and then spawns other threads that handle the
!> remaining tasks that are ready.
!> =============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  logical:: test, fell_through
  !.............................................................................
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  logical:: mytask
  integer:: n1, n2
  real(8):: start_fail, start_iter
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher3_t%execute', itimer=itimer)
  call self%init (task_list)
  self%task_list    => task_list
  ! -- this spoils a list that has been doctored by e.g. rt_setup
  !call task_list%init_all_nbors
  n_state = 0; n_state(0) = io%ntask
  call tic (time=start)
  call timer%print
  !$omp parallel private(link,nbor,task,mytask) shared(n_state) default(shared)
  !$omp master
  do while (task_list%na>0 .and. wallclock()-start < io%job_seconds)
    start_iter = wallclock()
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
          task_list%nq = n_state(2)
          if (use_critical) then
            !$omp critical (set_state_cr)
            mytask = omp%nthreads==1 .or. n_state(1)==1 .or. n_state(2)>f_over*omp%nthreads
            !$omp end critical (set_state_cr)
          else
            !$omp atomic read
            n1 = n_state(1)
            !$omp atomic read
            n2 = n_state(2)
            mytask = omp%nthreads==1 .or. n1==1 .or. n2>f_over*omp%nthreads
          end if
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
      if (self%verbose>3) print *,'dispatcher3: fell through'
      wc_failed = wc_failed + (wallclock()-start_fail - start_iter)
      if (use_taskyield) then
        n_fail = n_fail+1
        !$omp taskyield
      else
        n_fail = n_fail+1
      end if
    else
      if (n_fail>0 .and. self%verbose>0) &
        call self%print_state ('2')
      n_fail = 0
    end if
  end do
  !$omp end master
  !$omp end parallel
  call timer%print
  call finalize
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=start)
  call trace%end (itimer)
END SUBROUTINE execute

!> =============================================================================
!> =============================================================================
SUBROUTINE print_state (self, label)
  class(dispatcher3_t):: self
  character(len=*):: label
  !$omp critical (set_state_cr)
  print '("DISPATCHER3: ",a,f12.6,f9.3,4x,a,i8,6(4x,a,i4))', &
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
!> A task that is being actively updated is set in state 3, and when finished
!> it switched to state 4.   We then use the same OMP task to check if any if
!> the nbors have beome ready to update, switching the state if it has, and
!> leaving it to the master thread to actually activate the task.
!> =============================================================================
SUBROUTINE update (self, task_list, link, test)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  logical:: test
  integer:: n_ready
  !.............................................................................
  task => link%task
  !-----------------------------------------------------------------------------
  ! If task mem is OMP placed, ignore tasks allocated on a different thread,
  ! passing them on to a matching thead -- unless there are no more tasks to
  ! choose from, or the current thread is the master thread
  !-----------------------------------------------------------------------------
  !$omp atomic read
  n_ready = n_state(1)
  !$omp atomic
  n_ready = n_ready + n_state(2)
  if (omp%thread==task%mem_thread .or. n_ready<omp%nthreads .or.&
      omp%thread==0 .or. task%mem_thread==-1) then
    call self%set_state (task, 3)
    task%thread = omp%thread
    if (self%verbose>3) print fmt2, &
      wallclock(), omp%thread, 'updating task', task%id, 'from', task%time
    call task_list%update (link, test)
    if (self%verbose>3) print fmt2, &
      wallclock(), omp%thread, 'updated  task', task%id, '  to', task%time
    nbor => link%nbor
    do while (associated(nbor))
      if (self%is_ready(nbor%link)) then
        call self%set_state (nbor%link%task, 1)
        if (self%verbose>3) print fmt1,wallclock(),omp%thread,'found(2) task',nbor%task%id
      end if
      nbor => nbor%next
    end do
    if (verbose > 2) &
      print *, omp%thread, 'match  ', task%id
  else
    if (verbose > 2) &
      print *, omp%thread, 'skipped', task%id, task%mem_thread
    !$omp taskyield
  end if
  link%task%thread = -1
  call self%set_state (link%task, 4)
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
  !$omp task firstprivate(self,link,test,task_list) default(shared)
    if (use_locks) then
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
  if (use_critical) then
    !$omp critical (set_state_cr)
    io = task%state
    task%state = in
    n_state(io) = n_state(io) - 1
    n_state(in) = n_state(in) + 1
    if (self%verbose>1) &
      print '(a,i4,f12.6,2x,5i4,2x,2i3)', ' dispatcher3:', omp%thread, wallclock(), n_state, io, in
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
    if (self%verbose>1) &
      print '(a,i4,f12.6,2x,5i4,2x,2i3)', ' dispatcher3:', omp%thread, wallclock(), n_state, io, in
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
  !$omp atomic
  n_check = n_check+1
  call io%assert (associated(link), 'is_ready: link missing')
  call io%assert (associated(link%task), 'is_ready: link%task missing')
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
        if (verbose > 2) &
          print *, link%task%id, 'is_ready: failed  on', nbor%task%id, &
            nbor%task%time, link%task%time
        exit
      end if
      nbor => nbor%next
    end do
  end if
  if (is_ready) then
    !$omp atomic
    n_ready = n_ready+1
  end if
END FUNCTION

!> =============================================================================
!> Print statistics
!> =============================================================================
SUBROUTINE finalize
  if (omp%master) &
    print '(a,i6,4(a,f7.4))', ' dispatcher3_t%finalize: rank =',mpi%rank, &
      ', fraction of ready tasks =', real(n_ready)/real(n_check), &
      ', fraction of dispatcher idle time =', wc_failed/max(1d-10,wallclock()-start)
END SUBROUTINE finalize


END MODULE dispatcher3_mod
