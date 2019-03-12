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
  integer, save:: verbose=0
  integer, save:: stalled=0
  integer, save:: min_nq=2**30
  logical, save:: debug=.false.
  logical, save:: do_delay=.false.
  logical, save:: track_active=.false.
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
  logical:: queue_unpack, send_priv, recv_active, recv_priv
  integer:: iostat
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  namelist /dispatcher0_params/ verbose, max_stalled, retry_stalled, do_delay, &
    queue_unpack, send_priv, recv_active, recv_priv, debug
  character(len=120):: ids = &
  '$Id$ dispatchers/dispatcher0_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%print_id (ids)
  call trace%begin('dispatcher0_t%init')
  call mpi_mesg%init
  recv_active  = .true.
  queue_unpack = .true.
  recv_priv    = .true.
  send_priv    = .true.
  rewind (io%input); read(io%input, dispatcher0_params, iostat=iostat)
  mpi_mesg%recv_active  = recv_active
  mpi_mesg%queue_unpack = queue_unpack
  mpi_mesg%recv_priv    = recv_priv
  mpi_mesg%send_priv    = send_priv
  if (io%master) then
    if (io%master) write (*, dispatcher0_params)
  end if
  call self%lock%init
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
  call task_list%startup
  call tic (time=sec)
  call timer%print()
  !-----------------------------------------------------------------------------
  ! Loop over tasks
  !-----------------------------------------------------------------------------
  !$omp parallel
  do while (task_list%na > 0 .and. wallclock() < io%job_seconds)
    call self%update (task_list, test)
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
  call timer%print()
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=sec)
  write (io_unit%mpi,*) 'at mpi%barrier'
  flush (io_unit%mpi)
  call mpi%barrier ('end')
  write (io%output,*) "task list finished, min_nq =", min_nq
  if (validate%mode == "write") then
    call io%print_hl
    write (io%output,'(a)') &
      ' validate file '//trim(io%outputname)//'/rank_00000.val written'
    call io%print_hl
  else if (validate%mode == "compare") then
    call io%print_hl
    write (io%output,*) "validate%ok =", validate%ok
    call io%print_hl
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
  real(8):: wc
  real(8), save:: time, otime=0d0
  integer, save:: oid=0
  integer:: id
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin('dispatcher0_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Check incoming MPI, which may free up tasks for execution
  !-----------------------------------------------------------------------------
  call task_list%check_mpi                              ! check incoming MPI
  call mpi_io%iwrite_list%check                         ! I/O check
  !-----------------------------------------------------------------------------
  ! Now is a convenient time to check if a file should be opened or closed.  No
  ! harm is done if that makes the thread hangs for a while; other threads work
  !-----------------------------------------------------------------------------
  call data_io%open_and_close ()
  if (task_list%verbose > 0) &
    write (io%output,*) wallclock(),' thread',omp%thread,' waiting for tasklist(1)'
  call task_list%lock%set
  if (task_list%verbose > 0) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(1)'
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in what is normally a very brief critical region
  !-----------------------------------------------------------------------------
  nullify (prev)
  head => task_list%queue                                    ! queue start
  if (associated(head) .and. .not.task_list%syncing) then
    if (verbose > 0) &
      write(io_unit%log,'(f12.6,2x,a,i6,i7)') wallclock(), &
        'dispather0::update: na, id =',task_list%na, head%task%id
    !---------------------------------------------------------------------------
    ! If the patch mem was OMP placed, search ready queue for matching thread
    !---------------------------------------------------------------------------
    do while (associated(head))
      task => head%task
      if (task%mem_thread==-1 .or. omp%thread == task%mem_thread) exit
      prev => head
      head => head%next_time
    end do
    !---------------------------------------------------------------------------
    ! If not found, fall back on first task in queue
    !---------------------------------------------------------------------------
    if (.not. associated(head)) then
      nullify (prev)
      head => task_list%queue
      task => head%task                                   ! head task
    end if
    !---------------------------------------------------------------------------
    if (debug) then
      if (stalled > 0) then
        write (io_unit%log,*) wallclock(), '  stall ended', stalled
      end if
    end if
    stalled = 0
    !---------------------------------------------------------------------------
    task%atime = task%time
    if (io%verbose > 1) then
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
    ! A task generats a critical region three times: 1) when taken off the queue,
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
      if (io%verbose >= 0) then
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
  call task_list%lock%unset
  if (task_list%verbose > 0) &
    write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(1)'
  !-----------------------------------------------------------------------------
  ! As long as we are syncing, skip the update and come back.  The on thread
  ! that hit the sync time last of all threads has set this flag, and will clear
  ! it as soon as all other ranks have arrived at the same time.
  !-----------------------------------------------------------------------------
  if (associated(head)) then
    if (already_busy) then
      !!$omp critical (io_cr)
      !write (io_unit%log,*) mpi%rank,' WARNING: thread',omp%thread, &
      !  ' tried to update busy task', task%id
      print *, mpi%rank,' WARNING: thread',omp%thread, &
        ' tried to update busy task', task%id
      !!$omp end critical (io_cr)
      call trace_end (itimer)
      return
    end if
    call task_list%update (head, test, was_refined, was_derefined)
    !-----------------------------------------------------------------------------
    ! The task should now again be checked, before it can be returned to the
    ! queue, and its neighbors should also be checked for return to the queue.
    ! If the task was turned into a virtual task by the load balance, its virtual
    ! bit will prevent it being checked, but its nbors still need to be checked
    !-----------------------------------------------------------------------------
    if (.not. was_derefined) then
      call task%clear (bits%busy)           ! clear busy bit (FIXME: check this!)
      if (.not. was_refined) then
        call task_list%check_nbors (head)   ! any nbors ready?
        !call mpi_mesg%diagnostics (0)
      end if
    end if
  else
    call task_list%lock%set
    if (stalled == 0) then
      stall_start = wallclock()
      if (debug) &
        write (io_unit%log,*) wallclock(), 'queue stalled', stalled
    end if
    stalled = stalled+1
    call task_list%lock%unset
    call self%lock%set
    if (wallclock()-stall_start > retry_stalled) then
       !print *, mpi%rank, 'STALLED:'
       call task_list%check_all
       if (associated(task_list%queue)) then
         print *, mpi%rank, 'STALLED revived by check_all'
         stall_start = wallclock()
       else
         call task_list%check_oldest
         if (associated(task_list%queue)) then
           print *, mpi%rank, 'STALLED revived by check_oldest'
           stall_start = wallclock()
         else
           print *, mpi%rank, 'STALLED, revived FAILED'
           call io%abort ('STALLED, revived FAILED')
         end if
       end if
    else if (do_delay) then
      call mpi_mesg%delay (stalled)
    end if
    if (wallclock()-stall_start > max_stalled) then
      print *, mpi%rank, 'STALLED diagnostics'
      call mpi_mesg%diagnostics (1)
      print *, mpi%rank, 'STALLED bailing out'
      call mpi%abort ('exceeded max_stalled')
    end if
    call self%lock%unset
  end if
  call trace_end (itimer)
END SUBROUTINE update

END MODULE
