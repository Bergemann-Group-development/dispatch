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
!> =============================================================================
MODULE dispatcher3_mod
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
  USE patch_mod
  USE experiment_mod
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
  type(dispatcher3_t), public:: dispatcher3
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self, task_list)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  integer:: iostat
  namelist /dispatcher3_params/ use_locks, use_critical, use_taskyield, f_over, &
    master_works, n_write, verbose
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
  call io%header('begin dispatcher3_t%execute: Any opportunity dispatcher')
  call io%gb_print
  call self%init (task_list)
  self%task_list => task_list
  n_state = 0; n_state(0) = io%ntask
  call info (task_list) 
  call check (task_list)
  call tic (time=start)
  call timer%print
  !$omp parallel private(link,nbor,task,mytask) shared(n_state) default(shared)
  !$omp master
  do while (task_list%na>0 .and. wallclock()-start < io%job_seconds)
    start_iter = wallclock()
    !---------------------------------------------------------------------------
    ! Run through the task list, looking for and counting tasks that are ready
    !---------------------------------------------------------------------------
    fell_through = .true.
    link => task_list%head
    do while (associated(link))
      task => link%task
      !-------------------------------------------------------------------------
      ! For virtual tasks, check for incoming messages
      !-------------------------------------------------------------------------
      if ((mpi%size > 1) .and. task%is_set (bits%virtual)) then
        call check_recv (task, task_list)
      else
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
          if (experiment%is_ready (link)) then
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
        ! If task is ready to update, spawn a background task
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
    !---------------------------------------------------------------------------
    ! Check if there are MPI package buffers to reclaim
    !---------------------------------------------------------------------------
    if (mpi%size > 1) &
      call check_sent (mpi_mesg%sent_list)
    if (self%verbose>1) call self%print_state ('1')
    !---------------------------------------------------------------------------
    ! Diagnostics
    !---------------------------------------------------------------------------
    if (fell_through) then
      if (self%verbose>3) print *,'dispatcher3: fell through'
      call self%diagnostics (task_list)
      !wc_failed = wc_failed + (wallclock()-start_fail - start_iter)
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
  !-----------------------------------------------------------------------------
  call write_state (self, io_unit%output, label)
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE write_state (self, unit, label)
  class(dispatcher3_t):: self
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
  class(experiment_t), pointer:: task
  logical:: test
  integer:: n_ready
  !.............................................................................
  task => task2experiment (link%task)
  !-----------------------------------------------------------------------------
  ! If task mem is OMP placed, ignore tasks allocated on a different thread,
  ! passing them on to a matching thread -- unless there are no more tasks to
  ! choose from, or the current thread is the master thread
  !-----------------------------------------------------------------------------
  !$omp atomic read
  n_ready = n_state(1)
  !$omp atomic
  n_ready = n_ready + n_state(2)
  call self%set_state (link%task, 3)
!  !$omp atomic write
!  task%thread = omp%thread
  if (self%verbose>3) write(io%output,fmt2) &
    wallclock(), omp%thread, 'updating task', task%id, 'from', task%time
  !-----------------------------------------------------------------------------
  ! The link%task is updated, which includes calling link%send_to_vnbors(),
  ! which call mpi_mesg%sent() that adds an item to the list mpi_mesg%sent_list,
  ! in a critical region (sent_cr)
  !-----------------------------------------------------------------------------
  call task_list%update (link, test)
  !-----------------------------------------------------------------------------
  if (self%verbose>3) write(io%output,fmt2) &
    wallclock(), omp%thread, 'updated  task', task%id, '  to', task%time
  ! nbor => link%nbor
  ! do while (associated(nbor))
    ! if (nbor%needs_me) then
      ! if (experiment%is_ready(nbor%link)) then
        ! call self%set_state (nbor%link%task, 1)
        ! if (self%verbose>3) &
          ! write(io%output,fmt1) &
            ! wallclock(), omp%thread, 'found(2) task', nbor%task%id
      ! end if
    ! end if
    ! nbor => nbor%next
  ! end do
  if (verbose > 2) &
    write(io%output,*) omp%thread, 'match  ', task%id
!  link%task%thread = -1
  call self%set_state (link%task, 4)
END SUBROUTINE update

!> =============================================================================
!> Upgrade task pointer to experiment level
!> =============================================================================
FUNCTION task2experiment (task) RESULT (exper)
  class(task_t), pointer:: task
  class(experiment_t), pointer:: exper
  !...........................................................................
  select type (task)
  class is (experiment_t)
  exper => task
  end select
END FUNCTION task2experiment

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
  if (self%verbose>2) &
    write(io%output,fmt0) &
      wallclock(), omp%thread, 'spawning task', link%task%id, link%task%state, &
      link%task%thread
  !$omp task firstprivate(self,link,test,task_list) default(shared)
    if (use_locks) then
      if (link%task%lock%test()) then
        call self%update (task_list, link, test)
        call link%task%lock%unset
      else
        write(io%output,*) &
          'WARNING, lock conflict, thread, id =', &
          omp%thread, link%task%id, link%task%thread
      end if
    else
      call self%update (task_list, link, test)
    end if
  !$omp end task
END SUBROUTINE spawn_update

!> =============================================================================
!> Set a new task state, while keeping track of the number of tasks in each state
!> =============================================================================
SUBROUTINE set_state (self, task, in)
  class(dispatcher3_t):: self
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
    if (self%verbose > 1) &
      write (io%output, '(a,i4,f12.6,2x,5i4,2x,2i3)') &
        ' dispatcher3:', omp%thread, wallclock(), n_state, ip, in
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
        ' dispatcher3:', omp%thread, wallclock(), n_state, ip, in
  end if
END SUBROUTINE

!> =============================================================================
!> Print statistics
!> =============================================================================
SUBROUTINE finalize
  if (omp%master) &
    write(io%output,'(a,i6,4(a,f7.4))') &
      ' dispatcher3_t%finalize: rank =',mpi%rank, &
      ', fraction of ready tasks =', real(n_ready)/real(n_check), &
      ', fraction of dispatcher idle time =', wc_failed/max(1d-10,wallclock()-start)
END SUBROUTINE finalize

!> =============================================================================
!> =============================================================================
SUBROUTINE diagnostics (self, task_list)
  class(dispatcher3_t):: self
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

!> =============================================================================
!> Check messages on the mpi_mesg%sent_list for completeness.  Since threads
!> add items to the sent_list in a critical region (sent_cr), we use the same
!> =============================================================================
SUBROUTINE check_sent (list)
  type(mesg_list_t):: list
  !.............................................................................
  class(mesg_t), pointer:: mesg, next
  logical:: flag
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher3_t%check_sent', itimer=itimer)
  if (associated(list%head)) then
    !write(io_unit%mpi,*) wallclock(), 'check_sent: n =', mpi_mesg%sent_list%n
    !$omp critical (sent_cr)
    mesg => list%head                                           ! top of list
    do while (associated(mesg))
      next => mesg%next                                         ! in case removed
      call mesg%test_all (flag)                                 ! done with?
      if (flag) then
        write(io_unit%mpi,*) 'check_sent:     done', mesg%id
        call list%remove (mesg)                                 ! yes, remove
        call list%delete (mesg, send=.true.)                    ! deallocate
      else
        write(io_unit%mpi,*) wallclock(), 'check_sent: not done', mesg%id
      end if
      mesg => next                                              ! true next
    end do
    !$omp end critical (sent_cr)
  end if
  call trace%end (itimer)
END SUBROUTINE check_sent

!> =============================================================================
!> Check that task has a message data type, make sure it has an active recv req,
!> check if there is a message, unpack if so, and issue a new request.
!> =============================================================================
SUBROUTINE check_recv (task, task_list)
  class(task_t), pointer:: task
  type(task_list_t), pointer:: task_list
  class(patch_t), pointer:: patch
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('dispatcher3_t%check_recv', itimer=itimer)
  !write(io_unit%mpi,*) wallclock(), 'check_recv: incoming', task%id 
  if (.not.associated(task%mesg)) then
    select type (task)
    class is (patch_t)
    call task%allocate_mesg
    end select
    write(io_unit%mpi,*) wallclock(), 'check_recv: allocated', task%mesg%id 
    call task%mesg%irecv  (task%rank, task%id)
    write(io_unit%mpi,*) wallclock(), 'check_recv: irecv', task%mesg%id 
  end if
  if (task%mesg%test_recv()) then
    write(io_unit%mpi,*) wallclock(), 'check_recv: unpacking', task%mesg%id 
    call task_list%unpack (task%mesg)
    write(io_unit%mpi,*) wallclock(), 'check_recv: irecv', task%mesg%id 
    call task%mesg%irecv (task%rank, task%id)
  end if
  call trace%end (itimer)
END SUBROUTINE check_recv

!> =============================================================================
!> Dump nbor relations, flags and status bits around io%id_debug
!> =============================================================================
SUBROUTINE nbors_info (task_list, id_debug)
  type(task_list_t), pointer:: task_list
  integer, optional:: id_debug
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (experiment_t)
    if (present(id_debug)) then
      if (task%id == id_debug) then
          call task%nbors_info()
        exit
      end if
    else
      if (task%i_omega==0) then
        call task%nbors_info()
      end if
    end if
    end select
    link => link%next
  end do
END SUBROUTINE nbors_info

!> =============================================================================
!> Dump nbor relations, flags and status bits around io%id_debug
!> =============================================================================
SUBROUTINE info (task_list)
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link
  class(experiment_t), pointer:: task
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => task2experiment(link%task)
    call task%nbor_info(link)
    link => link%next
  end do
END SUBROUTINE info

!> =============================================================================
!> Check that all nbor task needed are in the task list
!> =============================================================================
SUBROUTINE check (task_list)
  class(task_list_t):: task_list
  class(link_t), pointer:: link, nbor
  logical:: ok
  !-----------------------------------------------------------------------------
  ok = .true.
  link => task_list%head
  do while (associated(link))
    nbor => link%nbor
    do while (associated(nbor))
      if (nbor%needed .and. .not. found(nbor%link)) then
        write (io_unit%mpi,*) &
          'link task, missing nbor task:', link%task%id, nbor%task%id
        ok = .false.
      end if
      nbor => nbor%next
    end do
    link => link%next
  end do
  if (ok) then
    write (io_unit%mpi,*) 'INFO: task_list checked: no missing tasks'
  else
    call mpi%abort ('ERROR: missing tasks')
  end if
contains
  !-----------------------------------------------------------------------------
  function found (test)
    logical:: found
    class(link_t), pointer:: link, test
    !...........................................................................
    found = .false.
    link => task_list%head
    do while (associated(link))
      if (associated(link,test)) then
        found = .true.
        exit
      end if
      link => link%next
    end do
  end function found
END SUBROUTINE check

END MODULE dispatcher3_mod
