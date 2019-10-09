!===============================================================================
!> Dispatcher method that relies on all threads maintaining a "ready queue", with
!> tasks ready for updating.  Each thread picks the task at the head of the queue,
!> and checks all nbors for tasks that are ready, putting them back in the queue.
!> The method is similar to method5, except the state of each task is given by
!> the value of task_t%state, rather than by status bits in task_t%status, and
!> locks are used to protect on the one hand the ready queue, and on the other
!> hand each individual task.
!>
!> READY QUEUE: The ready queue consists of links starting with task_list%queue,
!> and continuing with link%next_time.  Any change to one of these quantities
!> requires that the thread acquires the task_list%lock.  All threads will be
!> competing to do that, so the locked time must be as short as possible.
!>
!> TASK LOCKS:  Each task can change state from stat 0 (dormant) to state 1 
!> ready to become updated, state 2 (busy = being updated), state 3 (just updated),
!> and state 9 (finished).  Each such change of state must be done by first 
!> acquiring the task%lock, changing the state variable, and releasing the lock.
!> The state of task needs to be checked in various contexts, and to get a unique
!> answer on such a check the task%lock must be acquired during the check.  Some
!> of these checks are anyway done in connection with changing the state, so the
!> lock is anyway already acquired.  Other checks (e.g. the check on whether the
!> task is ahead or behind another task) do not occur in connection with a change
!> of state, and then the task%lock must be explicitly acquired before the check,
!> and then released.  Each lock set/unset may take of the order 0.1-0.3 micro-
!> seconds, so up to several dozen of these can be done without significantly
!> increasing the computing time, since each task update typically uses several
!> tens of milliseconds.  A task is typically acquired as many times as there are
!> nbors in the nbor list (or at most that many times), so there should be no
!> chance that the set/unset of locks could take significant time, and since a
!> specific task%lock has no impact on other tasks, this should be scalable to
!> any number of threads per process.
!>
!> When should a task be locked?  The most conservative choice is that a task is
!> locked the entire time it is being updated, so while being in state 2, while
!> the least conservative choice is that a task is locked only when it changes
!> state, and when the state is being tested.  The task cycle consists of a long
!> dormant time, a short (e.g. few percent) update time, and the very brief times
!> when it changes state.  The main outside reference to task data is when the
!> task time is compared to that of another task. The task time should only be
!> updated in locked state, going from state 2 to state 3.  While the task is in
!> state 0, 1, or 2, it may still have a sufficiently advanced time to be able
!> to serve guard zones, so there is no reason to exclude any state, or require
!> a particular state, while it is being tested.  It just needs to be locked by
!> the thread performing the comparison.
!>
!> When a task is used as a source for guard zone loading it may be wise and 
!> conservative to lock it, so the time slots cannot be rotated during the guard
!> zone date acquisition.  In principle this may not be necessary, at least if
!> one does NOT make use of the indirect addressing array iit((), which is being
!> changed during a rotate.  However, it may still be risky, since a time slot
!> that is being used during guard zone interpolations might be overwritten with
!> new data during the computations, if the task is not locked.  Since the time
!> a task is needed as a source for guard zone values is only 5-10% of the update
!> time, locking it briefly during this time may be ok.  There are, however, about
!> 26 such accesses per time step, so that could start to become a problem, if
!> the source task is locked during the whole interpolation time.  A better 
!> approach is to lock briefly, compute the memory slot needed, check if there is
!> a risk that it will become overwritten (this would be the case only if the
!> slot is nr 1 in the iit() array, and keep the task locked only if this is the
!> case.
!===============================================================================
MODULE dispatcher6_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_lib
  USE timer_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_io_mod
  USE mpi_mesg_mod
  USE link_mod
  USE list_mod
  USE task_mod
  USE experiment_mod
  USE bits_mod
  USE task_list_mod
  USE refine_mod
  USE patch_mod
  USE load_balance_mod
  USE task_mesg_mod
  !USE rt_mod
  implicit none
  private
  type, public:: dispatcher6_t
    integer:: verbose=0
    integer:: n_spawn=0
  contains
    procedure:: init
    procedure:: execute
  end type
  integer, save:: n_state(0:4)=0
  integer, save:: d_state(0:4)=0
  integer, save:: verbose=0
  integer, save:: stalled=0, max_stalled=10000, retry_stalled=100
  integer, save:: min_nq=2**30
  logical, save:: do_delay=.false.
  logical, save:: track_active=.false.
  type(dispatcher6_t), save:: virtual_list
  !$omp threadprivate (virtual_list)
  type(dispatcher6_t), public:: dispatcher6
CONTAINS

!===============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!===============================================================================
SUBROUTINE init (self, name)
  class(dispatcher6_t):: self
  character(len=*), optional:: name
  !.............................................................................
  class(link_t), pointer:: link
  logical:: queue_unpack, send_priv, recv_active, recv_priv
  integer:: iostat
  namelist /dispatcher6_params/ verbose, max_stalled, retry_stalled, do_delay
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  call trace%begin('dispatcher6_t%init')
  call mpi_mesg%init
  rewind (io%input)
  read(io%input, dispatcher6_params, iostat=iostat)
  write (io%output, dispatcher6_params)
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Execute the task list, updating it until it is empty.
!===============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher6_t):: self
  type(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  real(8):: sec
  integer:: dims(4)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher6_t%execute', itimer=itimer)
  call startup (task_list)
  call tic (time=sec)
  call timer%print
  !-----------------------------------------------------------------------------
  ! Loop over tasks
  !-----------------------------------------------------------------------------
  !$omp parallel
  do while (task_list%na > 0 .and. wallclock() < io%job_seconds)
    call update_list (task_list, test)
    if (io%do_stop) call mpi%abort ('stop flag')
    if (task_list%na == task_list%n_tasks) then
      !$omp atomic
      min_nq = min(min_nq,task_list%nq)
    end if
  end do
  write (io_unit%log,*) 'thread',omp%thread,' arrived'
  flush (io_unit%log)
  !$omp barrier
  !$omp end parallel
  call timer%print
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=sec)
  call mpi%barrier ('end')
  write (io%output,*) "task list finished, min_nq =", min_nq
  call trace%end (itimer)
END SUBROUTINE execute

!===============================================================================
!> Update the state of the task list, taking the steps necessary to update the
!> head task on the ready_queue, and check for consequences.  Two strategies:
!>
!> 1) The threads pick up tasks from the queue themselves, and put itself or
!>    other tasks back onto the queue, based on the results of list_t%check_ready
!>
!> In this case the threads are prevented from messing up for each other by
!> using critical regions when manipulating the queue, and by using status bits
!> to indicate the state each task is in.  These can be 1) in the queue, and
!> ready to be updated, 2) not in the queue, and busy, 3) not in the queue, and
!> not busy. The 'ready' bit is set while the task is in the queue, and the
!> 'busy' but is set while the task is busy being updated.  A task should be
!> checked for being ready only when in state 3, meaning only when none of the
!> two bits are set.
!>
!> 2) The master thread takes care of picking tasks from the queue, starting an
!>    OMP thread to update it, and checking for nbors that are ready to be
!>    updated as a consquence
!===============================================================================
SUBROUTINE update_list (task_list, test)
  type(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  class(link_t), pointer:: head, prev
  class(task_t), pointer:: task, otask
  class(mesg_t), pointer:: mesg
  logical:: already_busy, was_refined
  real(8):: wc
  real(8), save:: time, otime=0d0
  integer, save:: oid=0
  integer:: id, nq
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin('dispatcher6::update_list', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Check the sent_list for completed sends
  !-----------------------------------------------------------------------------
  if (mpi%size > 1) then
    !$omp atomic read
    nq = task_list%nq
    call mpi_mesg%sent_list%check_sent (nq)
  end if
  call mpi_io%iwrite_list%check                         ! I/O check
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in a very brief critical region
  !-----------------------------------------------------------------------------
  call task_list%lock%set()
  head => task_list%queue                                    ! queue start
  if (associated(head)) then
    task_list%queue => head%next_time                        ! chop head off
  end if
  call task_list%lock%unset()
  if (associated(head)) then
    task => head%task
    task%atime = task%time
    call task%lock%set()
    id = task%id
    if (task%state==1) then                                  ! verify ready state
      call set_state (task, 2)                               ! set busy state
      !$omp atomic
      task_list%nq = task_list%nq-1                          ! decrement queue count
      task%nq = task_list%nq                                 ! for info print
      call task%lock%unset()
      if (task%is_set(bits%virtual)) then
        call unpack (task_list, task%mesg, head)
      else
        call update_task (task_list, head, test, was_refined)
      end if
      if (.not. was_refined) then
        call check_nbors (task_list, head)                   ! any nbors ready?
      end if
      call task%lock%set()
      call set_state (task, 0)                               ! mark done
    else
      write (io%output,*) 'WARNING: queue head task was in state', task%state
    end if
    call task%lock%unset
  !-----------------------------------------------------------------------------
  ! If the queue is empty, check all tasks; this may be due to late packages
  !-----------------------------------------------------------------------------
  else
    call check_all (task_list)
  end if
  call trace_end (itimer)
END SUBROUTINE update_list

!===============================================================================
!> Update the task
!===============================================================================
SUBROUTINE update_task (task_list, head, test, was_refined)
  class(task_list_t):: task_list
  class(link_t), pointer:: head
  logical:: test
  logical, optional:: was_refined
  !.............................................................................
  class(task_t), pointer:: task
  real(8):: wc
  logical:: refined, derefined
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin('dispatcher6::update_task', itimer=itimer)
  task => head%task
  !-----------------------------------------------------------------------------
  ! Check if refinement is needed on the task; if so this will push new tasks
  ! onto the queue, with the same task time; i.e., to the head of the queue.
  ! If the task is virtual, it will not be checked by this rank.
  !-----------------------------------------------------------------------------
  call refine%check_current(task_list, head, refined, derefined)
  if (derefined) then
    call trace%end (itimer)
    return
  end if
  if (present(was_refined)) then
    was_refined = refined
  end if
  !-----------------------------------------------------------------------------
  ! Download nbor info -- this may be an experiment_t procedure, but if not,
  ! the call should be answered by a solver-specific procedure (which may or
  ! may not choose to call the generic patch guard zone handler in tasks/)
  !-----------------------------------------------------------------------------
  call task%dnload                                    ! download nbor info
  !-----------------------------------------------------------------------------
  ! Update the task, whatever that means (may include call to task%output)
  !-----------------------------------------------------------------------------
  task%rotated = .false.
  call task%update                                    ! update the task
  if (.not.task%rotated) then
    call rotate (task)
  end if
  call task%info (task_list%nq, task_list%na)                   ! print info on stdout
  !-----------------------------------------------------------------------------
  ! If the task is a boundary patch, first check if it should be given to
  ! another rank.  A patch that has been given to another rank becomes
  ! a virtual patch, and the load balance procedure has already sent those
  ! of its nbors that became new boundary patches over to relevant ranks,
  ! while the patch ittask_list is sent here, with bits%virtual+swap_reqest set.
  !-----------------------------------------------------------------------------
  if (task%is_set(bits%boundary)) then                ! boundary patch?
    head%task%nq = task_list%nq                            ! make sure to pass on
    if (load_balance%check_load (head)) then          ! sell?
      !$omp atomic
       task_list%na = task_list%na-1
      !$omp atomic
       task_list%nb = task_list%nb-1
      !$omp atomic
       task_list%nv = task_list%nv-1
       call task_list%count_status                    ! redundant? (FIXME)
    end if
    if (task%id == io%id_debug) &
      write(io_unit%mpi,*) &
        'DBG task_list_t%update: calling send_to_vnbors', task%id
    call task_list%send_to_vnbors (head)              ! send to virtual nbors
  end if
  !-----------------------------------------------------------------------------
  ! Tasks that have finished are subtracted from the task list count but are
  ! not removed, since their data may be needed by other tasks (including on
  ! other ranks).  They are, however, not added back onto the ready queue.
  !-----------------------------------------------------------------------------
  call task%lock%set()
  if (task%has_finished()) then                       ! finished:
    call set_state (task, 4)
    call task%lock%unset()
    call load_balance%active (.false.)                ! turn off load balancing
    !$omp atomic
    task_list%na = task_list%na-1
  else
    call task%lock%unset()
  end if
  call trace%end (itimer)
END SUBROUTINE update_task

!===============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready. 
!>
!> The link pointer and everything it points to are private to this task, and
!> are not at this point in time accessible from the ready queue.
!===============================================================================
SUBROUTINE check_nbors (task_list, link)
  class(list_t):: task_list
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('dispatcher6::check_nbors', itimer=itimer)
  task => link%task                                   ! main task
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    !print *, nbor%task%id, nbor%task%lock%id, nbor%task%lock%thread, omp_get_thred_num()
    call check_ready (task_list, nbor%link)         ! pointer back
    nbor => nbor%next                                 ! next nbor
  end do
  if (link%task%state==0) &
    call check_ready (task_list, link)                     ! finally check link task
  call trace_end (itimer)
END SUBROUTINE check_nbors

!===============================================================================
!> Check if the link task is ready, and if so, add it to the ready_queue
!===============================================================================
SUBROUTINE check_ready (task_list, link)
  class(list_t):: task_list
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  logical:: ok
  integer:: state
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher6::check_ready', itimer=itimer)
  task => link%task
  call task%lock%set()
  state = task%state
  call task%lock%unset()
  !-----------------------------------------------------------------------------
  ! Only check tasks that are not in the queue, busy, frozen, or external
  !-----------------------------------------------------------------------------
  if (task%state==0) then
    if (task%is_set(bits%virtual)) then
      !$omp critical (irecv_cr)
      ok = task%mesg%is_complete()
      if (ok) then
        if (task%mesg%req==0) then
          call task%mesg%irecv (task%rank, task%id)
        end if
      end if
      !$omp end critical (irecv_cr)
    else 
      ok = .true.
      nbor => link%nbor                                     ! start on nbor list
      do while (associated (nbor))                          ! keep going until end
        if (nbor%needed) then
          if (.not. is_ahead_of (nbor%task, task)) then     ! fail?
            ok = .false.
            exit
          end if
        end if
        nbor => nbor%next                                   ! next nbor
      end do
    end if
    if (ok) then
      call task%lock%set()
      call set_state (task, 1)                              ! go to state 1
      call queue_by_time (task_list, link)                  ! add task to queue
      call task%lock%unset()
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE check_ready

!===============================================================================
!===============================================================================
SUBROUTINE queue_by_time (task_list, this)
  class(list_t):: task_list
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  class(task_t), pointer:: task
  integer, save:: itimer=0
  integer:: nit
  !-----------------------------------------------------------------------------
  call trace_begin ('dispatcher6::queue_by_time ', itimer=itimer)
  call task_list%lock%set
  task => this%task
  mpi_mesg%n_ready = mpi_mesg%n_ready+1
  nullify (prev)
  next => task_list%queue
  do while (associated(next))
    if (associated(next%task, task)) then
      write (io_unit%log,*) omp_mythread, ' WARNING: task', task%id, ' is already in ready queue'
      go to 9
    else if (next%task%time > task%time) then
      this%next_time => next
      if (associated(prev)) then
        prev%next_time => this
      else
        task_list%queue => this
      end if
      !$omp atomic
      task_list%nq = task_list%nq+1
      go to 9
    end if
    prev => next
    next => next%next_time
  end do
  task_list%nq = task_list%nq+1
  if (associated(prev)) then
    prev%next_time => this
  else
    task_list%queue => this
  end if
  nullify (this%next_time)
9 continue
  call task_list%lock%unset
  call trace_end (itimer)
END SUBROUTINE queue_by_time

!===============================================================================
!> Check if source (which is a nbor task) is ahead of target (which is the one to
!> possibly move to the ready queue), using source%dtime*target%grace as the grace
!> period, since we want to limit the extrapolation in the nbor task to at most
!> target%grace*source%dtime.
!===============================================================================
LOGICAL FUNCTION is_ahead_of (source, target)
  class (task_t):: source, target
  !.............................................................................
  real(8):: nbtime, nbdtime, tgtime
  integer:: state, istep
  !-----------------------------------------------------------------------------
  call source%lock%set()
  nbtime = source%time
  nbdtime = source%dtime
  state = source%state
  istep = source%istep
  call source%lock%unset()
  !-----------------------------------------------------------------------------
  call target%lock%set()
  tgtime = target%time
  call target%lock%unset()
  !-----------------------------------------------------------------------------
  ! -- if a task is frozen (in time) assume it is forever ahead of other tasks
  if (state==4) then
    is_ahead_of = .true.
  ! -- do not use a grace interval for different levels or for the first steps
  else if (source%level /= target%level .or. istep < 3) then
    is_ahead_of = nbtime >= tgtime
  ! -- use a grace interval that is fraction of the nbor time step
  else
    is_ahead_of = nbtime + nbdtime*target%grace > tgtime
  end if
  if (io_unit%verbose>1) then
    if (target%id==io%id_debug.or.io_unit%verbose>4) &
    print'(i6,i4,2x,a,i6,3f9.5,l3)', source%id, omp_mythread, 'mk is_ahead_of: ', &
    target%id, nbtime, nbdtime*target%grace, tgtime, is_ahead_of
  end if
END FUNCTION is_ahead_of

!===============================================================================
!> Unpack a message, where the MPI tag is the task id.  Use that to search
!> for the task, apply its unpack method, and check if any nbors become ready.
!> This entire operation should be threadsafe, since no other thread should be
!> working on the same message and the same patch.
!===============================================================================
SUBROUTINE unpack (self, mesg, link)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer:: link
  class(link_t), pointer:: link2
  class(task_t), pointer:: task
  logical:: failed
  integer:: id
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher6::unpack', itimer=itimer)
  if (mesg%nbuf < 40) then
    call load_balance%unpack (mesg%buffer)
    return
  end if
  !$omp critical (unpack_cr)
  task => link%task
  failed = .false.
  !-----------------------------------------------------------------------------
  ! Guard against lingering extra messages to a swapped patch
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%boundary)) then
    write (io_unit%log,'(f12.6,2x,a,i9,3x,5l1)') wallclock(), &
      'task_mesg_mod::unpack ERROR, received mpi_mesg for boundary task:', task%id, &
      task%is_set(bits%internal), &
      task%is_set(bits%boundary), &
      task%is_set(bits%virtual), &
      task%is_set(bits%external), &
      task%is_set(bits%swap_request)
    failed = .true.
  end if
  if (task%is_set (bits%ready)) then
    write (io_unit%log,'(f12.6,2x,a,i9,3x,5l1)') wallclock(), &
      'task_mesg_mod::unpack ERROR, received mpi_mesg for task with ready bit:', &
      task%id, &
      task%is_set(bits%internal), &
      task%is_set(bits%boundary), &
      task%is_set(bits%virtual), &
      task%is_set(bits%external), &
      task%is_set(bits%swap_request)
    failed = .true.
  end if
  !---------------------------------------------------------------------------
  ! Unpack a patch message (which includes swapping the roles of boundary bits).
  ! Since an already existing patch may, at any one time, be under investigation
  ! by check_nbors, it must be protected by a critical region (or an OMP
  ! lock) while it is being updated here
  !---------------------------------------------------------------------------
  id = task%id
  call task%unpack (mesg)
  if (mpi_mesg%debug) &
    write (io_unit%log,'(f12.6,2x,a,i9,1p,e18.6)') wallclock(), &
      'unpk: id, time =', task%id, task%time
  if (id /= mesg%id) then
    write(io%output,'(i6,i4,2x,a,3i6)') mpi%rank, omp%thread, &
      'unpack ERROR: wrong mesg%id', id, task%id, mesg%id
    write(io_unit%log,'(f12.6,i6,i4,2x,a,4i6)') wallclock(), mpi%rank, &
      omp%thread, 'unpack ERROR: wrong mesg%id', mesg%sender, id, task%id, mesg%id
  end if
  if (.not. failed) then
    !$omp atomic
    mpi_mesg%n_unpk = mpi_mesg%n_unpk+1
    !$omp end atomic
    !---------------------------------------------------------------------------
    ! If the boundary+swap bits are set, this is a task that has just changed
    ! rank, and it needs to have its nbor relations re-initialized. This includes
    ! resorting (removing + re-adding) the nbor's nbor lists in rank order.
    ! FIXME: The load balancing steps should be checked for threadsafe operation
    !---------------------------------------------------------------------------
    if (task%is_set(bits%swap_request) .and. task%is_set(bits%boundary)) then
      self%na = self%na+1; self%nb = self%nb+1; self%nv = self%nv-1
      call self%init_nbors (link)
      call task%clear (bits%swap_request+bits%ready)
      call self%update_nbor_status (link)
      call self%count_status
      if (verbose>1) &
        write(io%output,'(f12.6,2x,a,i6,a,i9,a,i6)') &
        wallclock(),'LB: rank',mpi%rank,' given patch',task%id,' by',mesg%sender
      if (verbose>0) &
        write (io_unit%log,*) 'task_mesg_t%unpack: swapped virtual to boundary:', task%id
    !---------------------------------------------------------------------------
    ! If the link has no nbors it is a newly created virtual task. Does it need
    ! an nbor list?  At least we can use the nbor list to check that the link is
    ! in its nbors nbor lists.  A new virtual task (where no task existed) means
    ! that some nbor of it has changed from internal to boundary, which will be
    ! checked by the test_nbor_status call below, but only if an nbor list exists.
    !---------------------------------------------------------------------------
    !else if (.not.associated(link%nbor)) then
    else if (task%is_set(bits%swap_request) .and. task%is_set(bits%virtual)) then
      if (verbose>0) &
        write (io_unit%log,*) 'task_mesg_t%unpack: new virtual patch:', task%id
      self%nv = self%nv+1
      call self%init_nbors (link)
      call task%clear (bits%swap_request+bits%ready)
      call self%update_nbor_status (link)
      call self%count_status
    end if
  end if
  !$omp end critical (unpack_cr)
  call trace%end (itimer)
END SUBROUTINE unpack

!===============================================================================
!> Start-up preparation; initialize task message, look for updateable tasks
!===============================================================================
SUBROUTINE startup (task_list)
  class(task_list_t):: task_list
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher6::startup')
  call task_list%info
  !-----------------------------------------------------------------------------
  ! Initialize the task message, and request the first package
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%virtual)) then
      call task%allocate_mesg
      task%wc_last = wallclock()
      call task%mesg%irecv  (task%rank, task%id)
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Look for tasks ready to update
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    call set_state (task, 0)
    call check_ready (task_list, link)
    link => link%next
  end do
  call timer%print
  task_list%n_tasks = task_list%na
  !-----------------------------------------------------------------------------
  ! This may not be needed, but might avoid initial load balance excursions
  !-----------------------------------------------------------------------------
  call mpi%barrier ('task_list%execute')
  call trace%end()
END SUBROUTINE startup

!===============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!===============================================================================
SUBROUTINE check_all (list)
  class(list_t):: list
  class(link_t), pointer:: link
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('dispatcher6::check_all', itimer=itimer)
  call list%lock%set
  link => list%head
  do while (associated (link))                          ! keep going until end
    call check_ready (list, link)                       ! link ready?
    link => link%next
  end do
  call list%lock%unset
  call trace_end (itimer)
END SUBROUTINE check_all

!===============================================================================
!> Rotate time slots.  The initial conditions are in slot 1, and the first time
!> step puts new values in the 'new' slot 2, while saving the time step used in
!> dt(1). Then the current slot (it) becomes 2, and the new one becomes 3, etc.
!> This way, there is no need to copy memory btw time steps.
!===============================================================================
SUBROUTINE rotate (self)
  class(task_t):: self
  integer:: i
  integer, save:: itimer=0
  !.............................................................................
  if (self%rotated) return
  call trace_begin ('dispatcher6::rotate',itimer=itimer)
  call self%lock%set()
  !-----------------------------------------------------------------------------
  self%dt(self%it) = self%dtime                         ! just updated
  ! ZEUS determines a new time step at the end of its update; this is to
  ! prevent the clobbering of the updated time step; FIXME!
  if (trim(self%kind) /= 'zeus_mhd_patch') then
    self%dt(self%new)= self%dtime                       ! next estimate
  end if
  self%time = self%time + self%dtime                    ! time update
  self%t(self%new) = self%time                          ! initial time
  self%it = self%new                                    ! update time slot, new
  self%new = mod(self%new,self%nt)+1                    ! increment / rotate
  do i=1,self%nt-1
    self%iit(i) = self%iit(i+1)
  end do
  self%iit(self%nt) = self%new                          ! new right-most slot
  self%istep = self%istep + 1
  self%rotated = .true.
  !$omp flush
  call self%lock%unset()
  call trace_end (itimer)
END SUBROUTINE rotate

!===============================================================================
!> Set a new task state, while keeping track of the number of tasks in each state
!===============================================================================
SUBROUTINE set_state (task, in)
  class(task_t), pointer:: task
  integer:: in, ip
  !.............................................................................
  ip = task%state
  task%state = in
  n_state(ip) = n_state(ip) - 1
  n_state(in) = n_state(in) + 1
  d_state(in) = d_state(in) + 1
  if (verbose>1) &
    print '(a,i4,f12.6,2x,5i4,2x,2i3)', ' dispatcher6:', omp%thread, wallclock(), n_state, ip, in
END SUBROUTINE set_state

END MODULE
