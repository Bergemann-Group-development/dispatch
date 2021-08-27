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
  USE shared_mod
  implicit none
  private
  type, public:: dispatcher0_t
    integer:: verbose=0
    integer:: n_spawn=0
    type(lock_t), pointer:: lock
  contains
    procedure:: init
    procedure:: execute
    procedure:: update
    procedure:: no_queue
  end type
  real(8):: stall_start=0d0, max_stalled=600d0, retry_stalled=30d0
  real:: delay_ms=1.0
  real:: extra_grace=0.0            ! extra grace max, at nstall2
  integer, save:: mpi_only_master=100
  integer, save:: verbose=0
  integer, save:: stalled=0
  integer, save:: n_spin=0
  integer, save:: n_update=0
  integer, save:: nstall1=10000     ! after nstall1*delay_ms, increase grace
  integer, save:: nstall2=30000     ! after nstall2*delay_ms, call stall_handler
  integer, save:: min_nq=2**30
  logical, save:: debug=.false.
  logical, save:: do_delay=.true.
  logical, save:: track_active=.false.
  logical, save:: omp_pick=.true.
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
    delay_ms, mpi_only_master, debug, detailed_timer, nstall1, nstall2, &
    extra_grace, omp_pick
  character(len=120):: ids = &
  '$Id$'
  !-----------------------------------------------------------------------------
  call trace%print_id (ids)
  call trace%begin('dispatcher0_t%init')
  rewind (io%input)
  read(io%input, dispatcher0_params, iostat=iostat)
  call io%namelist_warning ('dispatcher0_params', iostat)
  write (io%namelist, dispatcher0_params)
  allocate (self%lock)
  call self%lock%init ('dispather0_t%init')
  call trace_end ()
END SUBROUTINE init

!===============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!===============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher0_t):: self
  class(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  real(8):: sec
  integer:: dims(4), nq
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher0_t%execute', itimer=itimer)
  call io%header ('dispatcher0: begin execution:')
  call task_list%startup
  call task_list%list_t%info
  call tic (time=sec)
!  call timer%print()
  !-----------------------------------------------------------------------------
  ! Loop over tasks
  !-----------------------------------------------------------------------------
  !$omp parallel
  call task_list%head%task%statistics
  do while (.not.timer%done .and. wallclock() < io%job_seconds)
    call self%update (task_list, test)
    if (task_list%na == task_list%n_tasks) then
      !$omp atomic read
      nq = task_list%nq
      !$omp atomic update
      min_nq = min(min_nq,nq)
    end if
  end do
  call task_list%head%task%statistics
  write (io_unit%log,*) 'thread',omp%thread,' arrived'
  flush (io_unit%log)
  !$omp barrier
  !$omp end parallel
!  call timer%print()
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
  class(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  class(link_t), pointer:: head, prev
  class(task_t), pointer:: task, otask
  class(mesg_t), pointer:: mesg
  real(8):: wc, start, time_incr
  real(8), save:: time, otime=0d0, spin_start=0d0
  logical:: already_busy
  integer, save:: oid=0
  integer:: i, id, nq, n_unpk
  integer:: stalled_l
  integer, save:: n_stalled=0
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
  if (task_list%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' waiting for tasklist(1)'
  if (detailed_timer) then
    call trace%end (itimer(i))
    i = i+1
    call trace%begin('dispatcher0_t%update(2)', itimer=itimer(i))
  end if
  !-----------------------------------------------------------------------------
  ! Check atomically if the ready queue is empty and return w/o locking if so,
  ! to avoid hammering the task list lock, delaying other threads.  Most of the
  ! time, the call to stall_handler() will come from here, but occasionally nq
  ! may be non-zero at the time of the test, but another thread picks up the
  ! last task before the test on head below, possibly resulting in a call to
  ! stall_handler, if n_stalled is already larger than 10,000.
  !-----------------------------------------------------------------------------
  !$omp atomic read
  nq = task_list%nq
  if (nq == 0) then
    if (n_stalled == 0) then
      !$omp atomic write
      spin_start = wallclock()
    end if
    !$omp atomic update
    n_stalled = n_stalled + 1
    call stall_handler
    call trace%end(itimer(i))
    return
  end if
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in what is normally a very brief OMP lock region
  !-----------------------------------------------------------------------------
                               call task_list%qlock%set ('dispatch0_t%update 1')
  head => task_list%queue                                    ! queue start
  if (detailed_timer) then
    call trace%end (itimer(i))
    i = i+1
    call trace%begin('dispatcher0_t%update(3)', itimer=itimer(i))
  end if
  nullify (prev)
  if (associated(head) .and. .not.task_list%syncing) then
    if (verbose > 0) &
      write(io_unit%log,'(f12.6,2x,a,i6,i7,1p,e12.4)') wallclock(), &
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
          exit
        end if
        prev => head
        head => head%next_time
      end do
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
    if (omp%thread == task%mem_thread) then
      !$omp atomic update
      timer%mem_hit = timer%mem_hit+1
    end if
    !$omp atomic update
    timer%mem_test = timer%mem_test+1
    call task%log ('dispatcher')
    !...........................................................................
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
        prev%next_time => head%next_time                     ! skip over
      else
        !print *, omp%thread, 'fall back'
        task_list%queue => head%next_time                    ! chop head off
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
          if (omp%master .and. verbose > 0) then
            write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,2i6)') &
              wallclock(), omp%thread, 'take', task%id, task%time, &
              task%istep, task_list%nq
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
                             call task_list%qlock%unset ('dispatch0_t%update 1')
  if (task_list%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(1)'
  !-----------------------------------------------------------------------------
  ! As long as we are syncing, skip the update and come back.  The on thread
  ! that hit the sync time last of all threads has set this flag, and will clear
  ! it as soon as all other ranks have arrived at the same time.
  !-----------------------------------------------------------------------------
  if (associated(head)) then
    nullify (head%next_time)
    if (detailed_timer) then
      call trace%end (itimer(i))
      i = i+1
      call trace%begin('dispatcher0_t%update(4)', itimer=itimer(i))
    end if
    if (already_busy) then
      !$omp critical (stderr_cr)
      write (stdlog,*) mpi%rank,' WARNING: thread', omp%thread, &
        ' tried to update busy task', task%id
      write (stderr,'(f12.6,i5,2x,a,i4,1x,a,i6)') &
        wallclock(), mpi%rank,' WARNING: thread', omp%thread, &
        'tried to update busy task', task%id
      !$omp end critical (stderr_cr)
      call trace%end (itimer(i))
      return
    end if
    call task_list%update (head, test)
    !---------------------------------------------------------------------------
    ! This taskyield is to encourage I/O tasks created in the mean time to run
    !---------------------------------------------------------------------------
    !$omp taskwait
    !---------------------------------------------------------------------------
    ! The first update after a stall period finds stalled non-zero, and sets it
    ! to zero, to prevent the if-block from being executed repeatedly
    !---------------------------------------------------------------------------
    !$omp atomic update
    timer%busy_time = timer%busy_time + (wallclock()-start)
    !$omp atomic read
    stalled_l = stalled
    !...........................................................................
    if (debug) then
      wc = wallclock()
      write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,2i6,f12.6)') &
        wc, omp%thread, 'updated', task%id, task%time, stalled_l, &
        task%istep, task%dtime
      !$omp atomic write
      n_spin = 0
    end if
    !...........................................................................
    if (stalled_l > 0) then
      !$omp atomic write
      stalled = 0
      !-------------------------------------------------------------------------
      ! We come here the first time the queue is non-empty, after being empty
      !-------------------------------------------------------------------------
      time_incr = wallclock()-spin_start
      !$omp atomic update
      timer%spin_time = timer%spin_time + time_incr
      !$omp atomic write
      n_stalled = 0
      !$omp atomic write
      spin_start = 0d0
    end if
    !---------------------------------------------------------------------------
    ! As the last thing in a successful update, reset the stall_start time, for
    ! use in diagnostic printout.  The actual start of a stall period is set
    ! the first time stall_handler() is called.
    !---------------------------------------------------------------------------
    !$omp atomic write
    stall_start = wallclock()
  else if (n_stalled == 0) then
    !---------------------------------------------------------------------------
    ! We come here the first time the queue is empty, after being non-empty,
    ! since n_stalled remains empty until we come here
    !---------------------------------------------------------------------------
    !$omp atomic write
    spin_start = wallclock()
    !$omp atomic update
    n_stalled = n_stalled + 1
  else if (n_stalled >= nstall1) then
    !$omp atomic write
    shared%extra_grace = (n_stalled-nstall1)*extra_grace/(nstall2-nstall1)
    write (io_unit%mpi,'(f12.6,2x,a,f10.3)') wallclock(), &
      'dispatcher0_t%update: extra_grace=', shared%extra_grace
    call task_list%check_all
  else if (n_stalled >= nstall2) then
    call stall_handler
    time_incr = wallclock()-spin_start
    !$omp atomic update
    timer%spin_time = timer%spin_time + time_incr
    !$omp atomic write
    n_stalled = 0
  else
    !---------------------------------------------------------------------------
    ! Continuing empty streak
    !---------------------------------------------------------------------------
    !$omp atomic update
    n_stalled = n_stalled+1
  end if
  call trace%end (itimer(i))
  return
contains
!-------------------------------------------------------------------------------
!> Stall handling, when the queue is empty
!-------------------------------------------------------------------------------
subroutine stall_handler
    real(8):: now
    integer, save:: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin ('dispatcher0_t%stall_handler', itimer=itimer)
    !$omp atomic update
    n_spin = n_spin+1
    !---------------------------------------------------------------------------
    ! As long as there are tasks, the stalled value remains zero.  The test
    ! below is triggered the first time a thread enters here, and the start
    ! time of stalling is registered, before incrementing the counter, which
    ! prevents the stall_start time from being set again.  The value of the
    ! 'stalled' counter is at least 1 when this procedure ends, and it should
    ! only be possible to find a zero value after one or more updates have been
    ! done, which means that logically, any 'stalled = 0' statement must come
    ! after the task_list%update, not before.
    !---------------------------------------------------------------------------
    !$omp atomic capture
    stalled_l = stalled
    stalled = stalled+1
    !$omp end atomic
    if (stalled_l == 0) then
      !$omp atomic write
      stall_start = wallclock()
      if (debug) then
        write (io_unit%queue,'(f12.6,i4,1x,a12)') &
          stall_start, omp%thread, 'noqueue'
      end if
    !---------------------------------------------------------------------------
    ! An unpacked MPI message counts as an update, and hence resets the stall
    !---------------------------------------------------------------------------
    else if (n_unpk > 0 .and. stalled_l > 0) then              ! unpack reset
      wc = wallclock()
      if (debug) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'unpack', &
          stalled_l, wc-stall_start, n_unpk
      end if
      !$omp atomic write
      stall_start = wc
      !$omp atomic write
      stalled = 0                                         ! reset stalled
    end if
    now = wallclock()
    !---------------------------------------------------------------------------
    ! To avoid that the task log is overwritten during a stall, stop task
    ! logging towards the of a log_interval
    !---------------------------------------------------------------------------
    if (io%task_logging > 0) then
      if (stall_start > 0d0 .and. now-stall_start > io%log_interval*0.9) then
        !$omp critical (tasklog_cr)
        if (io%task_logging > 0) then
          write (io_unit%task,'(f12.6,i4,4x,a)') &
            now, omp%thread, 'shutting down logging'
          flush (io_unit%task)
          write (stderr,'(1x,a,i5,3f12.5)') 'task logging stopped on rank', &
            mpi%rank, now, stall_start, io%log_interval*0.9
          write (stderr,*) 'task logging stopped on rank', mpi%rank
          io%task_logging = .false.
          !---------------------------------------------------------------------
          ! Check what's holding up tasks on this rank
          !---------------------------------------------------------------------
          call task_list%check_all
          call task_list%check_oldest
        end if
        !$omp end critical (tasklog_cr)
      end if
    end if
    if (stall_start > 0d0 .and. now-stall_start > max_stalled) then
      !$omp single
      write (io_unit%task,'(f12.6,i4,4x,a)') &
        now, omp%thread, 'max_stalled reached'
      print *, mpi%rank, 'STALLED diagnostics'
      call task_list%check_all
      call task_list%check_oldest
      call mpi_mesg%diagnostics (1)
      print *, mpi%rank, 'STALLED bailing out'
      call mpi%abort ('exceeded max_stalled')
      !$omp end single
    else if (stall_start > 0d0 .and. wallclock()-stall_start > retry_stalled) then
      write (io_unit%task,'(f12.6,i4,4x,a)') now, omp%thread, 'retry_stalled reached'
                                      call task_list%qlock%set ('stall_handler')
      if (wallclock()-stall_start > retry_stalled) then
        !-----------------------------------------------------------------------
        ! Try reviving the queue by using first check_all, which checks the whole
        ! task list for ready tasks -- finding such tasks shows that check_nbors()
        ! is not working correctly.
        !-----------------------------------------------------------------------
        print '(f12.6,i6,2x,a,2i6,f12.6)', now, mpi%rank, &
          'STALLED:', stalled, n_stalled, now-stall_start
        call task_list%check_all
        if (associated(task_list%queue)) then
          write(stderr,1) mpi%rank, omp%thread, 'check_all', wallclock()
        1 format("rank:",i5,2x,"thread:",i4,3x,"STALL revived by ",a," at",f12.3)
          if (verbose > 1) &
            call io%abort ('STALLED revided by check_all -- check thread log')
          stall_start = 0d0
        !-----------------------------------------------------------------------
        ! If that did not result in queued tasks, find the oldest tasks and the
        ! oldest virtual task.  Discovering that the queue is held up by such
        ! tasks shows that MPI communications are not synced correctly between
        ! ranks.
        !-----------------------------------------------------------------------
        else
          call task_list%check_oldest (clear=.true.)
          if (associated(task_list%queue)) then
            write(stderr,1) mpi%rank, omp%thread, 'check_oldest', wallclock()
            if (verbose > 1) &
              call io%abort ('STALLED revided by check_oldest -- check thread log')
          else
            write (stderr,*) mpi%rank, omp%thread, 'STALLED, revived FAILED'
            call io%abort ('STALLED, revived FAILED')
          end if
          stall_start = wallclock()
        end if
      end if
                                    call task_list%qlock%unset ('stall_handler')
    else if (do_delay) then
      call self%no_queue (stalled, task_list)
    end if
    call trace%end (itimer)
end subroutine stall_handler
END SUBROUTINE update

!===============================================================================
!> When the queue is empty, delay for about a task update time when MPI is
!> active, to avoid pounding MPI with a very large number of MPI_TEST calls.
!> For single rank runs, just avoid overflowing the rank log.
!===============================================================================
SUBROUTINE no_queue (self, stalled, task_list)
  class(dispatcher0_t):: self
  integer:: stalled
  class(task_list_t), pointer:: task_list
  real:: ms
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher0%no_queue', itimer=itimer)
  if (delay_ms > 0.0) then
    ms = delay_ms
  else
    ms = task_list%update_time*1e3
    if (mpi%size == 1) then
      ms = ms*0.1
    end if
  end if
  if (verbose > 0) then
    write (io_unit%log,'(f12.6,2x,a,i8,3f7.3)') wallclock(), &
     'dispatcher0_t%no_queue: stalled, update_time_ms, delay_ms =', &
     stalled, task_list%update_time*1e3, delay_ms, ms
  end if
  call omp%delay (ms)
  call trace%end (itimer)
END SUBROUTINE no_queue

END MODULE
