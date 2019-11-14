!===============================================================================
!> Dispatcher method that relies on all threads maintaining a "ready queue", with
!> tasks ready for updating.  Each thread picks the task at the head of the queue,
!> and checks all nbors for tasks that are ready, putting them back in the queue.
!===============================================================================
MODULE dispatcher0_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
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
  USE data_io_mod
  USE omp_lock_mod
  USE validate_mod
  USE refine_mod
  implicit none
  private
  type, public:: dispatcher0_t
    integer:: verbose=0
    integer:: n_spawn=0
    type(lock_t):: lock
  contains
    procedure:: init
    procedure:: execute
    procedure:: update
  end type
  real(8):: stall_start=0d0, max_stalled=600d0, retry_stalled=30d0
  integer, save:: mpi_only_master=100
  integer, save:: verbose=0
  integer, save:: stalled=0
  integer, save:: n_spin=0
  integer, save:: n_update=0
  integer, save:: min_nq=2**30
  logical, save:: debug=.false.
  logical, save:: do_delay=.false.
  logical, save:: track_active=.false.
  logical, save:: omp_pick=.false.
  logical, save:: detailed_timer=.false.
  type(dispatcher0_t), save:: virtual_list
  !$omp threadprivate (virtual_list)
  type(dispatcher0_t), public:: dispatcher0
CONTAINS

!===============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!===============================================================================
SUBROUTINE init (self, name)
  class(dispatcher0_t):: self
  character(len=*), optional:: name
  !.............................................................................
  class(link_t), pointer:: link
  integer:: iostat
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  namelist /dispatcher0_params/ verbose, max_stalled, retry_stalled, do_delay, &
    mpi_only_master, debug, detailed_timer
  character(len=120):: ids = &
  '$Id$ dispatchers/dispatcher0_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%print_id (ids)
  call trace%begin('dispatcher0_t%init')
  rewind (io%input)
  read(io%input, dispatcher0_params, iostat=iostat)
  write (io%output, dispatcher0_params)
  call mpi_mesg%init
  call self%lock%init ('disp')
  call self%lock%append
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!===============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher0_t):: self
  type(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  real(8):: sec
  integer:: dims(4)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher0_t%execute', itimer=itimer)
  call io%header('begin dispatcher0_t%execute:')
  call task_list%startup
  call tic (time=sec)
  call timer%print()
  !-----------------------------------------------------------------------------
  ! Loop over tasks
  !-----------------------------------------------------------------------------
  !$omp parallel
  do while (task_list%na > 0 .and. wallclock() < io%job_seconds)
    call self%update (task_list, test)
    if (task_list%na == task_list%n_tasks) then
      !$omp atomic
      min_nq = min(min_nq,task_list%nq)
    end if
  end do
  write (io_unit%log,*) 'thread',omp%thread,' arrived'
  flush (io_unit%log)
  !$omp barrier
  !$omp end parallel
  call timer%print()
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=sec)
  write (io_unit%mpi,*) 'at mpi%barrier'
  flush (io_unit%mpi)
  call mpi%barrier ('end')
  write (io%output,*) "task list finished, min_nq =", min_nq
  if (validate%mode == "write") then
    call io%print_hl()
    write (io%output,'(a)') &
      ' validate file '//trim(io%outputname)//'/rank_00000.val written'
    call io%print_hl()
  else if (validate%mode == "compare") then
    call io%print_hl()
    write (io%output,*) "validate%ok =", validate%ok
    call io%print_hl()
  end if
  call trace%end (itimer)
  write (io_unit%mpi,*) 'end dispatcher0_t%execute'
  flush (io_unit%mpi)
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
SUBROUTINE update (self, task_list, test)
  class(dispatcher0_t):: self
  type(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  class(link_t), pointer:: head, prev
  class(task_t), pointer:: task, otask
  class(mesg_t), pointer:: mesg
  logical:: already_busy, was_refined, was_derefined
  real(8):: wc, start
  real(8), save:: time, otime=0d0
  integer, save:: oid=0
  integer:: i, id, nq, n_unpk
  integer:: stalled_l
  integer, save:: itimer(5)=0
  !.............................................................................
  i = 1
  call trace%begin('dispatcher0_t%update(1)', itimer=itimer(i))
  start = wallclock()
  !$omp atomic update
  n_update = n_update+1
  !-----------------------------------------------------------------------------
  ! Check incoming MPI, which may free up tasks for execution
  !-----------------------------------------------------------------------------
  call task_list%check_mpi (n_unpk)                     ! check incoming MPI
  call mpi_io%iwrite_list%check                         ! I/O check
  if (omp%nthreads >= mpi_only_master .and. omp%master) then
    call trace%end (itimer(i))
    return
  end if
  !-----------------------------------------------------------------------------
  ! Now is a convenient time to check if a file should be opened or closed.  No
  ! harm is done if that makes the thread hangs for a while; other threads work
  !-----------------------------------------------------------------------------
  call data_io%open_and_close ()
  if (task_list%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' waiting for tasklist(1)'
  if (detailed_timer) then
    call trace%end (itimer(i))
    i = i+1
    call trace%begin('dispatcher0_t%update(2)', itimer=itimer(i))
  end if
  !-----------------------------------------------------------------------------
  ! Check atomically if the ready queue is empty and return w/o locking if so,
  ! to avoid hammering the task list lock, delaying other threads
  !-----------------------------------------------------------------------------
  !$omp atomic read
  nq = task_list%nq
  if (nq == 0) then
    call stall_handler
    call trace%end(itimer(i))
    return
  end if
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in what is normally a very brief OMP lock region
  !-----------------------------------------------------------------------------
  call task_list%lock%set ('dispatch0_t%update 1')
  head => task_list%queue                                    ! queue start
  if (detailed_timer) then
    call trace%end (itimer(i))
    i = i+1
    call trace%begin('dispatcher0_t%update(3)', itimer=itimer(i))
  end if
  nullify (prev)
  if (associated(head) .and. .not.task_list%syncing) then
    if (verbose > 0) &
      write(io_unit%log,'(f12.6,2x,a,i6,i7)') wallclock(), &
        'dispather0_t%update: na, id, time =', &
        task_list%na, head%task%id, head%task%time
    !---------------------------------------------------------------------------
    ! If the patch mem was OMP placed, search ready queue for matching thread
    !---------------------------------------------------------------------------
    if (omp_pick) then
      do while (associated(head))
        task => head%task
        if (task%mem_thread==-1) exit
        if (omp%thread == task%mem_thread) then
          !$omp atomic update
          timer%mem_hit = timer%mem_hit+1
          exit
        end if
        prev => head
        head => head%next_time
      end do
      !$omp atomic update
      timer%mem_test = timer%mem_test+1
      !-------------------------------------------------------------------------
      ! If not found, fall back on first task in queue
      !-------------------------------------------------------------------------
      if (.not. associated(head)) then
        nullify (prev)
        head => task_list%queue
        task => head%task                                   ! head task
      end if
    else
      task => head%task
    end if
    call task%log ('dispatcher')
    if (debug) then
      !$omp critical (wrt_cr)
      write (io_unit%queue,'(f12.6,2i7,3i5,f12.6,i6,2i9,7i7)') &
        wallclock(), task%id, task%istep, task%it, task%new, omp%thread, &
        task%time, task_list%nq, n_spin, timer%n_master
      timer%n_master(:) = 0
      n_spin = 0
      !$omp end critical (wrt_cr)
    end if
    !---------------------------------------------------------------------------
    if (debug) then
      !$omp atomic read
      stalled_l = stalled
      if (stalled_l > 0) then
        write (io_unit%log,*) wallclock(), '  stall ended', stalled_l
      end if
    end if
    !$omp atomic write
    stalled = 0
    !---------------------------------------------------------------------------
    task%atime = task%time
    if (verbose > 1) then
      write (io_unit%log,'(f12.6,2x,a,i4,2x,a,i6,2x,a,1p,g16.6,2x,a,i5)') &
        wallclock(), 'thread', omp%thread,'takes', task%id, &
        'time', task%time, 'nq', task%nq
      flush (io_unit%log)
    end if
    !---------------------------------------------------------------------------
    ! Set the head of the queue to the next task in time.  The ready bit needs
    ! to remain on for now, to prevent other threads from checking this task
    ! before it is updated.  After the update, the clearing of the ready bit
    ! allows the state of the current task to be evaluated by other threads, as
    ! well as by the current thread.  Locking should be used to allow only a
    ! single thread and task to check and change the ready bit status, at any
    ! one time.
    !
    ! The busy bit is set on tasks that a thread is updating, to check that a
    ! a task that has already been taken off the queue by another thread is not
    ! taken again.  This should in principle not be possible, since removing
    ! the task from the queue is done in a critical region, and removing a
    ! task can only be done once.
    !
    ! A task generates a critical region three times: 1) when taken off the queue,
    ! 2) when its time info is updated, 3) when it is added back to the queue.
    !
    ! To ensure 100% consistency, the task time information should not be
    ! allowed to change during the testing going on in list_t%check_ready().
    ! This could possibly be critical, if the fact that a task is / becomes
    ! ready is missed, because of sychronization issue.  Several of the tasks
    ! that form the nbor list of a task are typically undergoing updates at the
    ! same time.  In rare cases, but only if the update time of task time info
    ! is not uniquely defined, a task could thus be deemed not ready by the
    ! last of its nbors that has just been updated, because another task in the
    ! nbor list had not yet propagated its updated time info.
    !---------------------------------------------------------------------------
    already_busy = task%is_set (bits%busy)
    id = task%id
    if (.not.already_busy) then
      call task%set (bits%busy)                              ! mark task busy
      if (associated(prev)) then
        !print *, omp%thread, 'match'
        !call prev%qlock%set ('dipatcher0')
        prev%next_time => head%next_time                     ! skip over
        !call prev%qlock%unset ('dipatcher0')
      else
        !print *, omp%thread, 'fall back'
        !call task_list%queue%qlock%set ('dipatcher0')
        task_list%queue => head%next_time                    ! chop head off
        !call task_list%queue%qlock%unset ('dipatcher0')
      end if
      call task%clear (bits%ready)                           ! not in queue
      if (track_active) &
        call task_list%queue_active (head)
      !$omp atomic
      task_list%nq = task_list%nq-1                          ! decrement queue count
      task%nq = task_list%nq                                 ! for info print
      if (verbose > 0) then
        if (track_active) then
          otask => task_list%active%task
          time = otask%atime
          if (time < otime) &
            write (*,'(a,2(f12.6,i6))') 'TIME ERROR: otime, oid, time, id =', otime, oid, time, otask%id
          !$omp atomic write
          otime = time
          !$omp atomic write
          oid = otask%id
          if (omp%master) then
            write (io_unit%queue,'(f12.6,2i7,2f12.6,2i5)') &
              wallclock(), task%id, task%istep, task%time, &
              otask%atime, task_list%nq, task_list%nac
            flush (io_unit%queue)
          end if
        else
          if (omp%master) then
            write (io_unit%queue,'(f12.6,2i7,f12.6,2i5)') &
              wallclock(), task%id, task%istep, task%time, &
              task_list%nq, task_list%nac
            flush (io_unit%queue)
          end if
        end if
      end if
    end if
    !---------------------------------------------------------------------------
    ! The first time a task is at the head of the ready queue with time =
    ! sync_time we can be sure that all other tasks are also being upated
    ! (by other threads) to arrive at this time.  We can then temporarily
    ! halt updating, and wait for this task to finish its update, at a barrier
    ! where all other ranks are doing the same.
    !---------------------------------------------------------------------------
    if (task%time == task_list%sync_next) then
      task_list%syncing = .true.
      write (io_unit%mpi,*) task%id,omp%thread, &
        ' is triggering a sync at t =', task_list%sync_next
    end if
  else
    if (verbose > 0) &
      write(io_unit%log,'(f12.6,i6,2x,a)') wallclock(), task_list%na, 'no queue'
  end if
  call task_list%lock%unset ('dispatch0_t%update 1')
  if (task_list%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(1)'
  !-----------------------------------------------------------------------------
  ! As long as we are syncing, skip the update and come back.  The on thread
  ! that hit the sync time last of all threads has set this flag, and will clear
  ! it as soon as all other ranks have arrived at the same time.
  !-----------------------------------------------------------------------------
  if (associated(head)) then
    if (detailed_timer) then
      call trace%end (itimer(i))
      i = i+1
      call trace%begin('dispatcher0_t%update(4)', itimer=itimer(i))
    end if
    if (already_busy) then
      !$omp critical (stderr_cr)
      write (io_unit%log,*) mpi%rank,' WARNING: thread', omp%thread, &
        ' tried to update busy task', task%id
      write (stderr,*) mpi%rank,' WARNING: thread', omp%thread, &
        ' tried to update busy task', task%id
      !$omp end critical (stderr_cr)
      call trace%end (itimer(i))
      return
    end if
    call task_list%update (head, test, was_refined, was_derefined)
    !---------------------------------------------------------------------------
    ! The task should now again be checked, before it can be returned to the
    ! queue, and its neighbors should also be checked for return to the queue.
    ! If the task was turned into a virtual task by the load balance, its virtual
    ! bit will prevent it being checked, but its nbors still need to be checked
    !---------------------------------------------------------------------------
    if (.not. was_derefined) then
      call task%clear (bits%busy)
      call task_list%check_nbors (head)   ! any nbors ready?
    end if
    !$omp atomic update
    timer%busy_time = timer%busy_time + (wallclock()-start)
  else
    call stall_handler
  end if
  call trace%end (itimer(i))
contains
!-------------------------------------------------------------------------------
!> Stall handling, when the queue is empty
!-------------------------------------------------------------------------------
subroutine stall_handler
    !$omp atomic update
    n_spin = n_spin+1
    if (detailed_timer) then
      call trace%end (itimer(i))
      i = i+1
      call trace%begin('dispatcher0_t%update(5)', itimer=itimer(i))
    end if
    !---------------------------------------------------------------------------
    ! As long as there are tasks, the stalled value remains zero.  The test
    ! below is triggered the first time a thread enters here, and the start
    ! time of stalling is registered, before incrementing the counter, which
    ! prevents the stall_start time from being set again.
    !---------------------------------------------------------------------------
    !$omp atomic read
    stalled_l = stalled
    !if (n_unpk > 0 .and. stalled_l > 0) then              ! unpack reset
    !  write (io_unit%log,'(a,i7,i4,f12.6)') &
    !    'dispatcher0_t%update: unpack reset, n, time =', &
    !    stalled_l, n_unpk, wallclock()-stall_start
    !  !$omp atomic write
    !  stalled = 0                                         ! reset stalled
    !  stalled_l = 0                                       ! new stall interval
    !end if
    if (stalled_l == 0) then
      !$omp atomic write
      stall_start = wallclock()
      if (debug) &
        write (io_unit%log,*) wallclock(), 'queue stalled', stalled_l
    end if
    !$omp atomic update
    stalled = stalled+1
    if (wallclock()-stall_start > max_stalled) then
      print *, mpi%rank, 'STALLED diagnostics'
      call mpi_mesg%diagnostics (1)
      print *, mpi%rank, 'STALLED bailing out'
      call mpi%abort ('exceeded max_stalled')
    else if (wallclock()-stall_start > retry_stalled) then
      !$omp critical (stall_cr)
      if (wallclock()-stall_start > retry_stalled) then
        !print *, mpi%rank, 'STALLED:'
        call task_list%check_all
        if (associated(task_list%queue)) then
          write(stderr,1) mpi%rank, omp%thread, 'check_all', wallclock()
        1 format("rank:",i5,2x,"thread:",i4,3x,"STALL revived by ",a," at",f12.3)
          if (verbose > 1) &
            call io%abort ('STALLED revided by check_all -- check thread log')
          stall_start = wallclock()
        else
          call task_list%check_oldest
          if (associated(task_list%queue)) then
            write(stderr,1) mpi%rank, omp%thread, 'check_oldest', wallclock()
            if (verbose > 1) &
              call io%abort ('STALLED revided by check_oldest -- check thread log')
            stall_start = wallclock()
          else
            write (stderr,*) mpi%rank, omp%thread, 'STALLED, revived FAILED'
            call io%abort ('STALLED, revived FAILED')
          end if
        end if
      end if
      !$omp end critical (stall_cr)
    else if (do_delay) then
      !$omp atomic update
      timer%spin_time = timer%spin_time + (wallclock()-start)
      call mpi_mesg%delay (stalled)
      return
    end if
    !$omp atomic update
    timer%spin_time = timer%spin_time + (wallclock()-start)
end subroutine stall_handler
END SUBROUTINE update

END MODULE
