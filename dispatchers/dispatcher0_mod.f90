!===============================================================================
!> Dispatcher method that relies on all threads maintaining a "ready queue", with
!> tasks ready for updating.  Each thread picks the task at the head of the queue,
!> and checks all nbors for tasks that are ready, putting them back in the queue.
!===============================================================================
MODULE dispatcher0_mod
  USE io_mod
  USE os_mod
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
  real(8):: stall_start=0d0, max_stalled=1e8, retry_stalled=1e8
  real(8):: soft_mem_limit=huge(1.d0)
  real:: delay_ms=0.05
  integer, save:: mpi_only_master=5
  integer, save:: verbose=0
  integer, save:: stalled=0
  integer, save:: n_update=0
  integer, save:: min_nq=2**30
  logical, save:: debug=.false.
  logical, save:: do_delay=.true.
  logical, save:: track_active=.false.
  logical, save:: omp_pick=.false.
  logical, save:: pick_thread=.false.
  logical, save:: detailed_timer=.false.
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
    delay_ms, mpi_only_master, debug, detailed_timer, omp_pick, pick_thread, soft_mem_limit
  character(len=120):: ids = &
  '$Id$'
  !-----------------------------------------------------------------------------
  call trace%print_id (ids)
  call trace%begin('dispatcher0_t%init')
  rewind (io%input)
  read(io%input, dispatcher0_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read(io%input, dispatcher0_params)
  end if
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
  integer:: dims(4), nq, n
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher0_t%execute', itimer=itimer)
  call io%header ('dispatcher0: begin execution:')
  ! call timer%tic (time=sec)
!  call timer%print()
  !-----------------------------------------------------------------------------
  ! Loop over tasks, keeping track of the shortest queue size until finish
  !-----------------------------------------------------------------------------
  !$omp parallel
  call task_list%head%task%statistics
  do while (.not.timer%done .and. wallclock() < io%job_seconds)
    call update (self, task_list, test)
    if (task_list%first_finished == 0d0) then
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
  n = timer%n_update
  ! call timer%toc ('wall time', n, time=sec, verbose=.true.)
  write (io_unit%mpi,*) 'at mpi%barrier'
  flush (io_unit%mpi)
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
  !call mpi%end ()
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
  class(link_t), pointer:: link, prev, first
  class(task_t), pointer:: task, otask
  class(mesg_t), pointer:: mesg
  real(8):: now, start, time_incr
  real(8), save:: time, otime=0d0, spin_start=0d0
  character(24):: label
  logical:: already_busy, reset
  integer, save:: oid=0
  integer:: i, id, nq, n_unpk, n_unpk_L
  integer, save:: itimer(5)=0
  !-----------------------------------------------------------------------------
  ! Check incoming MPI, which may free up tasks for execution
  !-----------------------------------------------------------------------------
  start = wallclock()
  n_unpk = 0
  if (mpi%size > 1 .and. omp%thread == 0) then
    call task_list%check_recv (n_unpk)                   ! check incoming MPI
    if (omp%nthreads > mpi_only_master) then            ! master exclusive?
      if (n_unpk > 1) then
        !$omp atomic write
        stalled = 0
        !$omp atomic write
        stall_start = start
        !$omp atomic write
        spin_start = 0d0
      end if
      return
    end if
  end if
  i = 1
  call trace%begin('dispatcher0_t%update(1)', itimer=itimer(i))
  if (stall_start==0d0) then
    !$omp atomic write
    stall_start = start
  end if
  !$omp atomic update
  n_update = n_update+1
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
  ! Pick a task off the queue, in what is normally a very brief OMP lock region
  !-----------------------------------------------------------------------------
  call                              task_list%qlock%set ('dispatch0_t%update 1')
  link => task_list%queue                                    ! queue start
  if (detailed_timer) then
    call trace%end (itimer(i))
    i = i+1
    call trace%begin('dispatcher0_t%update(3)', itimer=itimer(i))
  end if
  nullify (prev)
  nullify (first)
  !---------------------------------------------------------------------------
  ! If the patch mem was OMP placed, search ready queue for matching thread
  !---------------------------------------------------------------------------
  if (omp_pick) then
    do while (associated(link))
      task => link%task
      if (task%place == -1) then
        exit
      endif
      prev => link
      link => link%next_time
    end do
    !-------------------------------------------------------------------------
    ! If not found, fall back on first task in queue
    !-------------------------------------------------------------------------
    if (.not. associated(link)) then
      if (os%gb < soft_mem_limit .or. omp%thread == 0) then
        nullify (prev)
        link => task_list%queue
      endif
    end if
  end if
  if (associated(link) .and. .not.task_list%syncing) then
    if (verbose > 0) &
      write(io_unit%log,'(f12.6,2x,a,i6,i7,1p,e12.4)') wallclock(), &
        'dispather0_t%update: na, id, time =', &
        task_list%na, link%task%id, link%task%time
    task => link%task
    !---------------------------------------------------------------------------
    ! Bail out if the busy bit is already set
    !---------------------------------------------------------------------------
    already_busy = task%is_set (bits%busy)
    if (task%is_set(bits%busy)) then
      call task%log ('already busy')
      write (stdlog,'(a,i4,a,i9)') &
        ' WARNING: thread', omp%thread, &
        ' found task with busy-bit set:', task%id
      call trace%end (itimer(i))
      return
    end if
    !---------------------------------------------------------------------------
    ! Change the task state from ready to busy.  The task is owned by the
    ! current thread as long as it is under update or post-update, so there
    ! should be no ambiguity possible about the bit status
    !---------------------------------------------------------------------------
    task%atime = task%time
    id = task%id
    call task%set (bits%busy)                              ! mark task busy
    call task%clear (bits%ready)                           ! clear ready state
    call task%log ('set busy', 2)
    !...........................................................................
    write (label,'(a,i6,l2)') &
      trim(task%kind), task_list%nq, associated(task_list%queue)
    call task%log (trim(label), 1)
    call task%log (trim(label), 4)
    !---------------------------------------------------------------------------
    ! Keep track of memory affinity
    !---------------------------------------------------------------------------
    if (pick_thread) then
      if (omp%thread == task%thread) then
        !$omp atomic update
        timer%mem_hit = timer%mem_hit+1
      end if
    else
      if (omp%place == task%place) then
        !$omp atomic update
        timer%mem_hit = timer%mem_hit+1
      end if
    end if
    !$omp atomic update
    timer%mem_test = timer%mem_test+1
    if (verbose > 1) then
      write (io_unit%log,'(f12.6,2x,a,i4,2x,a,i6,2x,a,1p,g16.6,2x,a,i5)') &
        wallclock(), 'thread', omp%thread,'takes', task%id, &
        'time', task%time, 'nq', task%nq
      flush (io_unit%log)
    end if
    !---------------------------------------------------------------------------
    ! Remove the accepted task link from the queue
    !---------------------------------------------------------------------------
    if (associated(prev)) then
      prev%next_time => link%next_time                     ! skip over
    else
      task_list%queue => link%next_time                    ! chop head off
    end if
    if (track_active) &
      call task_list%queue_active (link)
    !$omp atomic capture
    task_list%nq = task_list%nq-1                          ! decrement queue count
    task%nq = task_list%nq                                 ! for info print
    !$omp end atomic
    !...........................................................................
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
    !$omp atomic write
    task_list%nq = 0
    if (verbose > 0) &
      write(io_unit%log,'(f12.6,i6,2x,a)') wallclock(), task_list%na, 'no queue'
  end if
  call                            task_list%qlock%unset ('dispatch0_t%update 1')
  if (task_list%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(1)'
  !-----------------------------------------------------------------------------
  ! Update task list or handle empty queue
  !-----------------------------------------------------------------------------
  reset = .false.
  now = wallclock()
  if (associated(link)) then
    nullify (link%next_time)
    call trace%end (itimer(i))
    i = i+1
    call trace%begin('dispatcher0_t%update(task)', itimer=itimer(i))
    !-------------------------------------------------------------------------
    ! Initially, the task%cost should be zero, so the first time the task
    ! list cost remains zero, but the next time a task updates, the previous
    ! cost of the task is subtracted and the new cost is added, which should
    ! keep the task_list%cost always equal to the sum of all task costs.
    !-------------------------------------------------------------------------
    !associate (task_list_cost=>task_list%cost)
    !!$omp atomic update
    !task_list_cost = task_list_cost - task%cost
    call task_list%update (link, test)
    !-------------------------------------------------------------------------
    ! The start time was set before the MPI and other service calls, to make
    ! sure that their contribution too the task list cost is counted.
    !-------------------------------------------------------------------------
    !cost = (wallclock()-start)/max(task%dtime,tiny(1d0))
    !task%cost = cost
    !!$omp atomic
    !task_list_cost = task_list_cost + task%cost
    !call task%log ('cost', 4, task_list%nq, r1=cost, d1=task_list_cost)
    !end associate
    !...........................................................................
    if (abs(io%task_logging) > 1) then
      write (label,'(a,i6,l2)') &
        'done', task_list%nq, associated(task_list%queue)
      call task%log (trim(label), 2)
    end if
    !---------------------------------------------------------------------------
    ! Force outstanding I/O tasks to complete (taskyield is not enough)
    !---------------------------------------------------------------------------
    !$omp taskwait
    !---------------------------------------------------------------------------
    ! The first update after a stall period finds stalled non-zero, and sets it
    ! to zero, to prevent the if-block from being executed repeatedly
    !---------------------------------------------------------------------------
    now = wallclock()
    !$omp atomic update
    timer%busy_time = timer%busy_time + (now-start)
    reset = .true.
  !-----------------------------------------------------------------------------
  ! An MPI unpack counts as a stall reset event
  !-----------------------------------------------------------------------------
  else if (n_unpk > 0) then
    reset = .true.
  !-----------------------------------------------------------------------------
  ! Call stall handler if time limit exceeded
  !-----------------------------------------------------------------------------
  else if (now-stall_start > retry_stalled) then
    call stall_handler (now, stall_start)
  !-----------------------------------------------------------------------------
  ! Continue a streak of no queue
  !-----------------------------------------------------------------------------
  else
    !$omp atomic update
    stalled = stalled + 1
    !---------------------------------------------------------------------------
    ! Accumulate spin time
    !---------------------------------------------------------------------------
    time_incr = now-spin_start
    !$omp atomic update
    timer%spin_time = timer%spin_time + time_incr
    !$omp atomic write
    spin_start = now
    !---------------------------------------------------------------------------
    ! Continuing empty streak, after registering a brief delay
    !---------------------------------------------------------------------------
    if (do_delay) then
      call self%no_queue (stalled, task_list)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Conditionally reset the local stall counter and increment spin time
  !-----------------------------------------------------------------------------
  if (reset) then
    !$omp atomic write
    stalled = 0
    !$omp atomic write
    stall_start = now
    !$omp atomic write
    spin_start = 0d0
  end if
  call trace%end (itimer(i))
contains
!-------------------------------------------------------------------------------
!> Stall handling, when the queue is empty
!-------------------------------------------------------------------------------
subroutine stall_handler (now, stall_start)
    real(8):: now, stall_start
    integer, save:: itimer=0
    !---------------------------------------------------------------------------
    call trace%begin ('dispatcher0_t%stall_handler', itimer=itimer)
    !---------------------------------------------------------------------------
    ! To avoid that the task log is overwritten during a stall, stop task
    ! logging towards the end of a log_interval
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
          io%task_logging = 0
        end if
        !$omp end critical (tasklog_cr)
      end if
    end if
    !---------------------------------------------------------------------------
    ! End of retry processing
    !---------------------------------------------------------------------------
    if (now-stall_start > max_stalled) then
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
    !---------------------------------------------------------------------------
    ! Retry processing
    !---------------------------------------------------------------------------
    else if (now-stall_start > retry_stalled) then
      !-------------------------------------------------------------------------
      ! Allow only one thread here, and also lock the queue
      !-------------------------------------------------------------------------
      call                                 task_list%qlock%set ('stall_handler')
      call                                       self%lock%set ('stall_handler')
      !-------------------------------------------------------------------------
      ! Try reviving the queue by using first check_all, which checks the whole
      ! task list for ready tasks -- finding such tasks shows that check_nbors()
      ! is not working correctly.
      !-------------------------------------------------------------------------
      if ((.not.associated(task_list%queue)) .and. &
          (now-stall_start > retry_stalled)) then
        write (io_unit%task,'(f12.6,i4,4x,a)') &
          now, omp%thread, 'retry_stalled reached'
        call task_list%check_all
        if (associated(task_list%queue)) then
          write(stderr,1) mpi%rank, omp%thread, 'check_all', wallclock()
        1 format("rank:",i5,2x,"thread:",i4,3x,"STALL revived by ",a," at",f12.3)
          if (verbose > 1) &
            call io%abort ('STALLED revided by check_all -- check thread log')
          stall_start = now
        end if
      end if
      !-------------------------------------------------------------------------
      ! If that did not result in queued tasks, find the oldest tasks and the
      ! oldest virtual task.  Finding that the queue is held up by such tasks
      ! shows that MPI communications are not synced correctly between ranks.
      !-------------------------------------------------------------------------
      if ((.not.associated(task_list%queue)) .and. &
          (now-stall_start > retry_stalled)) then
        call task_list%check_oldest (clear=.true.)
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
      call                                     self%lock%unset ('stall_handler')
      call                               task_list%qlock%unset ('stall_handler')
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
  if (mpi%size>1) timer%done = mpi%end_counter%update(0) >= mpi%size
  call trace%end (itimer)
END SUBROUTINE no_queue

END MODULE
