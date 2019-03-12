!===============================================================================
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
!> on an update:  If a process has N threads, and has N tasks active, evenly
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
!===============================================================================
MODULE dispatcher4_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_io_mod
  USE timer_mod
  USE link_mod
  USE mpi_mesg_mod
  USE task_list_mod
  USE bits_mod
  USE global_mod
  USE task_mod
  implicit none
  private
  type, public:: dispatcher4_t
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
    procedure:: diagnostics
  end type
  type(global_t):: average_density
  integer:: n_fail=0
  integer:: n_state(0:4)=0
  integer:: d_state(0:4)=0
  integer:: verbose=0
  integer:: n_write=2000
  integer(8):: n_check=0, n_ready=0
  real:: f_over=2.0
  real(8):: wc_failed=0d0
  real(8):: start
  logical:: master_works=.true.
  logical:: use_locks=.false.
  logical:: use_critical=.true.
  logical:: use_taskyield=.false.
  character(len=32):: fmt0='(    f12.6,i4,2x,a,2x,i7,2i4)'
  character(len=32):: fmt1='(40x,f12.6,i4,2x,a,2x,i7,2i4)'
  character(len=48):: fmt2='(40x,f12.6,i4,2x,a,2x,i7,2x,a,1p,g14.5)'
  type(dispatcher4_t), public:: dispatcher4
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, task_list)
  class(dispatcher4_t):: self
  type(task_list_t), pointer:: task_list
  integer:: iostat
  namelist /dispatcher4_params/ use_locks, use_critical, use_taskyield, f_over, &
    master_works, n_write, verbose
  !.............................................................................
  rewind (io%input); read (io%input, dispatcher4_params, iostat=iostat)
  if (io%master) write (io%output, dispatcher4_params)
  self%verbose = verbose
END SUBROUTINE init

!===============================================================================
!> Simple dispatcher, which in each loop through the task list takes on the
!> first task that is ready, and then spawns other threads that handle the
!> remaining tasks that are ready.
!===============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher4_t):: self
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
  call trace%begin ('dispatcher4_t%execute', itimer=itimer)
  call io%header('begin dispatcher4_t%execute: Any opportunity dispatcher')
  call self%init (task_list)
  self%task_list    => task_list
  call task_list%info
  n_state = 0; n_state(0) = io%ntask
  call tic (time=start)
  call timer%print()
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
    fell_through = .true.
    link => task_list%head
    do while (associated(link))
      task => link%task
      if (.not.task%is_set (bits%virtual)) then
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
        if (task%state==0) then
          if (self%is_ready (link)) then
            call self%set_state (link%task, 1)
          else if (link%task%time == 0.0_8) then
            write(io_unit%queue,*) link%task%id, 'ERROR: at time=0, but not ready!'
            nbor => link%nbor
            do while (associated(nbor))
              write(io_unit%queue,'(i7,i3,2f10.6,2x,2l3)') &
                nbor%task%id, nbor%task%state, (nbor%task%time-link%task%time), nbor%task%dtime, &
                nbor%needed, nbor%task%is_ahead_of (link%task)
              nbor => nbor%next
            end do
          end if
        end if
        !-----------------------------------------------------------------------
        ! If task is ready to update, pawn a background task
        !-----------------------------------------------------------------------
        if (task%state==1) then
          fell_through = .false.
          !$omp atomic read
          task_list%nq = n_state(1)
          !$omp atomic
          task_list%nq = task_list%nq  + n_state(2)
          if (omp%nthreads > 1) then
            call self%spawn_update (task_list, link, test)
          else
            call self%update (task_list, link, test)
          end if
        end if
      end if
      link => link%next
    end do
    call write_state (self, io_unit%queue, '1')
    if (self%verbose>1) call self%print_state ('1')
  end do
  !$omp end master
  !$omp end parallel
  call timer%print()
  call finalize
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=start)
  call trace%end (itimer)
END SUBROUTINE execute

!===============================================================================
!===============================================================================
SUBROUTINE print_state (self, label)
  class(dispatcher4_t):: self
  character(len=*):: label
  !-----------------------------------------------------------------------------
  call write_state (self, io_unit%output, label)
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE write_state (self, unit, label)
  class(dispatcher4_t):: self
  integer:: unit
  character(len=*):: label
  !-----------------------------------------------------------------------------
  if (n_write <= 0) return
  !$omp critical (set_state_cr)
  n_write = n_write-1
  write(unit,'("DISPATCHER3: ",a,f12.6,f9.3,2x,3(2x,a,i6),9(2x,a,i5))') &
    label, wallclock(), &
    wc_failed/max(1d-10,wallclock()-start), &
    'fail:',n_fail    , &
    'idl:', n_state(0), &
    'rdy:', n_state(1), &
    'spw:', n_state(2), &
    'act:', n_state(3), &
    'upd:', n_state(4), &
    'idl:', d_state(0), &
    'rdy:', d_state(1), &
    'spw:', d_state(2), &
    'act:', d_state(3), &
    'upd:', d_state(4)
  d_state = 0
  !$omp end critical (set_state_cr)
END SUBROUTINE write_state

!===============================================================================
!> A task that is being actively updated is set in state 3, and when finished
!> it switched to state 4.   We then use the same OMP task to check if any if
!> the nbors have beome ready to update, switching the state if it has, and
!> leaving it to the master thread to actually activate the task.
!===============================================================================
SUBROUTINE update (self, task_list, link, test)
  class(dispatcher4_t):: self
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  logical:: test
  integer:: n_ready
  !.............................................................................
  call self%set_state (link%task, 3)
  call task_list%update (link, test)
  call self%set_state (link%task, 4)
END SUBROUTINE update

!===============================================================================
!> Spawn an OpenMP task, for updating a DISPATCH task
!===============================================================================
SUBROUTINE spawn_update (self, task_list, link, test)
  class(dispatcher4_t):: self
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link
  logical:: test
  !.............................................................................
  call self%set_state (link%task, 2)
  !$omp task firstprivate(link) default(shared)
  call self%set_state (link%task, 3)
  call task_list%update (link, test)
  call self%set_state (link%task, 4)
  !$omp end task
END SUBROUTINE spawn_update

!===============================================================================
!> Set a new task state, while keeping track of the number of tasks in each state
!===============================================================================
SUBROUTINE set_state (self, task, in)
  class(dispatcher4_t):: self
  class(task_t), pointer:: task
  integer:: in, ip
  !...........................................................................
  if (use_critical) then
    !$omp critical (set_state_cr)
    ip = task%state
    task%state = in
    n_state(ip) = n_state(ip) - 1
    n_state(in) = n_state(in) + 1
    d_state(in) = d_state(in) + 1
    if (self%verbose>1) &
      print '(a,i4,f12.6,2x,5i4,2x,2i3)', ' dispatcher4:', omp%thread, wallclock(), n_state, ip, in
    !$omp end critical (set_state_cr)
  else
   !$omp atomic read
   ip = task%state
   !$omp atomic
   n_state(ip) = n_state(ip) - 1
   !$omp atomic
   n_state(in) = n_state(in) + 1
   !$omp atomic
   d_state(in) = d_state(in) + 1
   !$omp atomic write
   task%state = in
   if (self%verbose>1) &
     write (io%output,'(a,i4,f12.6,2x,5i4,2x,2i3)') &
       ' dispatcher4:', omp%thread, wallclock(), n_state, ip, in
  end if
END SUBROUTINE

!===============================================================================
!> Check if a task is ready to update, by changing if all nbor task are "ahead
!> of" the task in question
!===============================================================================
LOGICAL FUNCTION is_ready (self, link)
  class(dispatcher4_t):: self
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
  else
    is_ready = .true.
    nbor => link%nbor
    do while (associated(nbor))
      if (nbor%needed) then
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
      end if
      nbor => nbor%next
    end do
  end if
  if (is_ready) then
    !$omp atomic
    n_ready = n_ready+1
  end if
END FUNCTION is_ready

!===============================================================================
!> Print statistics
!===============================================================================
SUBROUTINE finalize
  if (omp%master) &
    write(io%output,'(a,i6,4(a,f7.4))') &
      ' dispatcher4_t%finalize: rank =',mpi%rank, &
      ', fraction of ready tasks =', real(n_ready)/real(n_check), &
      ', fraction of dispatcher idle time =', wc_failed/max(1d-10,wallclock()-start)
END SUBROUTINE finalize

!===============================================================================
!===============================================================================
SUBROUTINE diagnostics (self, task_list)
  class(dispatcher4_t):: self
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link, nbor, culprit
  !-----------------------------------------------------------------------------
  !write (io_unit%queue,*) 'begin starved analysis: nstate =', nstate
  link => task_list%head
  do while (associated(link))
    nbor => link%nbor
    link => link%next
  end do
  !write (io_unit%queue,*) 'end starved analysis: nstate =', nstate
END SUBROUTINE diagnostics

END MODULE dispatcher4_mod
