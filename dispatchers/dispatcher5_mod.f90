!> =============================================================================
!> Dispatcher method that relies on all threads maintaining a "ready queue", with
!> tasks ready for updating.  Each thread picks the task at the head of the queue,
!> and checks all nbors for tasks that are ready, putting them back in the queue.
!> =============================================================================
MODULE dispatcher5_mod
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
  USE refine_mod
  USE patch_mod
  USE load_balance_mod
  USE task_mesg_mod
  implicit none
  private
  type, public:: dispatcher5_t
    integer:: verbose=0
    integer:: n_spawn=0
  contains
    procedure:: init
    procedure:: execute
    procedure:: update
  end type
  integer, save:: verbose=0
  integer, save:: stalled=0, max_stalled=10000, retry_stalled=100
  integer, save:: min_nq=2**30
  logical, save:: do_delay=.false.
  logical, save:: track_active=.false.
  type(dispatcher5_t), save:: virtual_list
  !$omp threadprivate (virtual_list)
  type(dispatcher5_t), public:: dispatcher5
CONTAINS

!> =============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!> =============================================================================
SUBROUTINE init (self, name)
  class(dispatcher5_t):: self
  character(len=*), optional:: name
  !.............................................................................
  class(link_t), pointer:: link
  logical:: queue_unpack, send_priv, recv_active, recv_priv
  integer:: iostat
  namelist /dispatcher5_params/ verbose, max_stalled, retry_stalled, do_delay, &
    queue_unpack, send_priv, recv_active, recv_priv
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  call trace%begin('dispatcher5_t%init')
  call mpi_mesg%init
  recv_active  = .true.
  queue_unpack = .true.
  recv_priv    = .true.
  send_priv    = .true.
  rewind (io%input); read(io%input, dispatcher5_params, iostat=iostat)
  mpi_mesg%recv_active  = recv_active
  mpi_mesg%queue_unpack = queue_unpack
  mpi_mesg%recv_priv    = recv_priv
  mpi_mesg%send_priv    = send_priv
  if (io%master) then
    if (io%master) write (*, dispatcher5_params)
  end if
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!> =============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher5_t):: self
  type(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  real(8):: sec
  integer:: dims(4)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher5_t%execute', itimer=itimer)
  call startup (task_list)
  call tic (time=sec)
  call timer%print
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
  call timer%print
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=sec)
  call mpi%barrier ('end')
  write (io%output,*) "task list finished, min_nq =", min_nq
  call trace%end (itimer)
END SUBROUTINE execute

!> =============================================================================
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
!> =============================================================================
SUBROUTINE update (self, task_list, test)
  class(dispatcher5_t):: self
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
  call trace%begin('dispatcher5_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Check incoming MPI, which may free up tasks for execution
  !-----------------------------------------------------------------------------
  !call task_list%check_mpi                              ! check incoming MPI
  if (mpi%size > 1) then
    !$omp atomic read
    nq = task_list%nq
    call mpi_mesg%sent_list%check_sent (nq)
  end if
  call mpi_io%iwrite_list%check                         ! I/O check
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in a very brief critical region
  !-----------------------------------------------------------------------------
  !$omp critical (queue_cr)
  nullify (prev)
  head => task_list%queue                                    ! queue start
  if (associated(head) .and. .not.task_list%syncing) then
    if (verbose > 0) &
      write(io_unit%log,'(f12.6,i6,i7)') wallclock(), task_list%na, head%task%id
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
            write (io_unit%queue,'(f12.6,i6,2f12.6,2i5)') &
              wallclock(), task%istep, task%time, &
              otask%atime, task_list%nq, task_list%nac
            flush (io_unit%queue)
            end if
        else
            if (omp%master) then
            write (io_unit%queue,'(f12.6,i6,f12.6,2i5)') &
              wallclock(), task%istep, task%time, &
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
  !$omp end critical (queue_cr)
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
    if (task%is_set(bits%virtual)) then
      call unpack (task_list, task%mesg, head)
    else
      call update2 (task_list, head, test, was_refined)
    end if
    !-----------------------------------------------------------------------------
    ! The task should now again be checked, before it can be returned to the
    ! queue, and its neighbors should also be checked for return to the queue.
    ! If the task was turned into a virtual task by the load balance, its virtual
    ! bit will prevent it being checked, but its nbors still need to be checked
    !-----------------------------------------------------------------------------
    call task%clear (bits%busy)                         ! clear busy bit (FIXME: check this!)
    if (.not. was_refined) then
      call check_nbors (task_list, head)                ! any nbors ready?
      !$omp atomic write
      stalled = 0
      !call mpi_mesg%diagnostics (0)
    end if
  else
    call check_all (task_list)
  end if
  call trace_end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Common part of task list update, used by dispatchers
!> =============================================================================
SUBROUTINE update2 (self, head, test, was_refined)
  class(task_list_t):: self
  class(link_t), pointer:: head
  logical:: test
  logical, optional:: was_refined
  !.............................................................................
  class(task_t), pointer:: task
  real(8):: wc
  logical:: refined
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  ! Check for flag files.  If io%out_time is set, set a new out_next.
  ! If task%time is slightly smaller than a multiple of out_time, then
  ! out_next will becoe that multiple.
  !----------------------------------------------------------------------------
  call trace%begin('dispatcher5::update2', itimer=itimer)
  task => head%task
  call io%check_flags
  if (task%id==io%id_track) then
    task%track =  .not.task%track
    io%id_track = 0
  end if
  !-----------------------------------------------------------------------------
  ! Check if refinement is needed on the task; if so this will push new tasks
  ! onto the queue, with the same task time; i.e., to the head of the queue.
  ! If the task is virtual, it will not be checked by this rank.
  !-----------------------------------------------------------------------------
  call refine%check_current(self, head, refined)
  if (present(was_refined)) then
    was_refined = refined
  end if
  !---------------------------------------------------------------------------
  ! If the task is frozen (e.g. because it has finished), return
  !---------------------------------------------------------------------------
  if (task%is_set (bits%frozen)) then
    !if (task%iout == 0) call task%output(self%name)
    ! if a frozen task is set to ready, make it not ready
    if (task%is_set(bits%ready)) then
      !$omp atomic
      self%na = self%na - 1
    end if
    call trace%end (itimer)
    return
  end if
  !----------------------------------------------------------------------------
  ! Download nbor info -- this may be an experiment_t procedure, but if not,
  ! the call should be answered by a solver-specific procedure (which may or
  ! may not choose to call the generic patch guard zone handler in tasks/)
  !----------------------------------------------------------------------------
  if (.not.test) &
    call task%dnload                                  ! download nbor info
  !----------------------------------------------------------------------------
  ! Update the task, whatever that means (may include call to task%output)
  !----------------------------------------------------------------------------
  if (mpi_mesg%delay_ms==0.0) then
    wc = wallclock()
  end if
  task%rotated = .false.
  if (test) then
    call task%test_update
  else
    call task%update                                  ! update the task
  end if
  !----------------------------------------------------------------------------
  ! Default delay time = first update time
  !----------------------------------------------------------------------------
  if (mpi_mesg%delay_ms==0.0) then
    !$omp atomic write
    mpi_mesg%delay_ms = 1d3*(wallclock()-wc)
  end if
  !$omp atomic
  mpi_mesg%n_update = mpi_mesg%n_update+1
  if (io%verbose>1) &
    write (io_unit%log,'(a,i4,2x,a,i7,2x,a,2x,a,i7,2x,a,1p,2g14.6,2x,a,2i5,l3)') &
    'thread', omp_mythread, 'task', task%id, trim(task%type), &
    'step', task%istep, 'dt, time:', task%dtime, task%time, &
    'n, nq', self%n, self%nq, associated(head%nbor)
  !-----------------------------------------------------------------------------
  ! The rotate procedure in the generic patch_t module is responsible for
  ! updating the time, and rotating the memory slots where information about
  ! previous time steps are saved.  Task may signal that they have already done
  ! the rotate internally, or else are not yet ready to do so, by setting the
  ! task%rotated flag
  !-----------------------------------------------------------------------------
  if (.not.(task%rotated)) then
    call task%rotate
  end if
  select type (task)
  class is (patch_t)
    timer%n_update = timer%n_update + product(task%n)
  end select
  call task%info (self%nq, self%na)                   ! print info on stdout
  !-----------------------------------------------------------------------------
  ! Tasks that have finished are subtracted from the task list count but are
  ! not removed, since their data may be needed by other tasks (including on
  ! other ranks).  They are, however, not added back onto the ready queue.
  !-----------------------------------------------------------------------------
  if (task%has_finished()) then                       ! finished:
    call task%set (bits%frozen)
    call load_balance%active (.false.)                ! turn off load balancing
    !$omp atomic
    self%na = self%na-1
  end if
  !-----------------------------------------------------------------------------
  ! If the task is a boundary patch, first check if it should be given to
  ! another rank.  A patch that has been given to another rank becomes
  ! a virtual patch, and the load balance procedure has already sent those
  ! of its nbors that became new boundary patches over to relevant ranks,
  ! while the patch itself is sent here, with bits%virtual+swap_reqest set.
  !-----------------------------------------------------------------------------
  if (task%is_set(bits%boundary)) then                ! boundary patch?
    head%task%nq = self%nq                            ! make sure to pass on
    if (load_balance%check_load (head)) then          ! sell?
      !$omp atomic
       self%na = self%na-1
      !$omp atomic
       self%nb = self%nb-1
      !$omp atomic
       self%nv = self%nv-1
       call self%count_status                         ! redundant? (FIXME)
    end if
    if (task%id == io%id_debug) &
      write(io_unit%mpi,*) &
        'DBG task_list_t%update: calling send_to_vnbors', task%id
    call head%send_to_vnbors                          ! send to virtual nbors
  end if
  !-----------------------------------------------------------------------------
  ! Periodic task syncronization:  The first time a task%time is exactly equal
  ! to load_balance%sync_next, we wait for all ranks to arrive at the same
  ! time. All tasks within a rank are then also at this time, since as long as
  ! one has not arrived there yet, the head%ask%time is smaller.
  !-----------------------------------------------------------------------------
  if (self%syncing) then
    write (io_unit%mpi,*) task%id,omp%thread, &
      ' is waiting on a sync at t =', self%sync_next
    !call trace%end; call trace%begin ('mpi%barrier', itimer=itimer)
    call mpi%barrier ('sync')
    !call trace%end (itimer); call trace%begin('dispatcher0_t%update')
    write (io_unit%mpi,*) task%id,omp%thread, &
      ' finished wating on a sync at t =', self%sync_next
    self%sync_next = self%sync_next + self%sync_time
    self%syncing = .false.
  end if
  task%sync_time = self%sync_next
  call trace%end (itimer)
END SUBROUTINE update2

!> =============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE check_ready (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  logical:: ok, debug, debug1
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher5::check_ready', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we first check if any bits are set indicating this to be the case
  !-----------------------------------------------------------------------------
  task => link%task
  debug = io%verbose > 1 .or. task%id==io%id_debug
  if (task%is_set(bits%virtual)) then
    ok = task%mesg%test_recv()
    if (ok) then
      call task%mesg%irecv (task%rank, task%id)
    end if
  else if (task%is_clear(bits%ready+bits%busy+bits%frozen+bits%external)) then
    ok = .true.
    nbor => link%nbor                                     ! start on nbor list
    do while (associated (nbor))                          ! keep going until end
      debug1 = nbor%task%id==io%id_debug
      if (debug1) write(io_unit%log,*) 'DBG check_ready: id, nbor, needed, ahead', &
        task%id, nbor%task%id, nbor%needed, nbor%task%is_ahead_of(task)
      if (nbor%needed) then
        if (nbor%task%is_ahead_of(task)) then             ! and is ahead in time
          if (debug1) &
            write (io_unit%log,'("DBG list_t%check_ready:",i5,f10.6,a,i6,2f10.6,f6.2,l4)') &
              task%id, task%time, '  is OK on', nbor%task%id, nbor%task%time, &
              nbor%task%dtime, nbor%task%grace, nbor%task%is_set(bits%virtual)
        else                                              ! not ahead in time
          if (debug1) &
            write (io_unit%log,'("DBG list_t%check_ready:",i5,f10.6,a,i6,2f10.6,f6.2,l4)') &
              task%id, task%time, ' failed on', nbor%task%id, nbor%task%time, &
              nbor%task%dtime, nbor%task%grace, nbor%task%is_set(bits%virtual)
          ok = .false.
          exit
        end if
      end if
      nbor => nbor%next                                   ! next nbor
    end do
  end if
  if (ok) then
    if (debug) &
      write (io_unit%output,'("list_t%check_ready:",i5,f10.6,a)') &
        task%id, task%time, ' succeded'
    if (task%is_clear(bits%ready+bits%busy+bits%frozen+bits%external)) &
      call self%queue_by_time (link)                      ! add task to queue
  end if
  call trace%end (itimer)
END SUBROUTINE check_ready

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready.  If the ready bit is already set it means the patch has already
!> been put in the ready queue, and should not be checked again.
!>
!> The link pointer and everything it points to are private to this task, and
!> are not at this point in time accessible from the ready queue.
!> =============================================================================
SUBROUTINE check_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('dispatcher5::check_nbors', itimer=itimer)
  if (io%verbose > 1) &
    write (io_unit%log,*) 'checking nbors of',link%task%id,associated(link%nbor)
  !.... critical (queue_cr)
  task => link%task                                   ! main task
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    if (nbor%task%id==io%id_debug) print *, &
      'task', task%id,' needs task ', nbor%task%id, nbor%needs_me
!    if (nbor%needs_me) then                          ! needs checking?
      if (io%verbose > 1) &
        write (io_unit%log,*) 'check nbor', nbor%task%id
      call check_ready (self, nbor%link)              ! pointer back
!    else
!      if (io%verbose > 1) &
!        write (io_unit%log,*) 'no need to', nbor%task%id
!    end if
    nbor => nbor%next                                 ! next nbor
  end do
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we check if any bits are set indicating this to be the case
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%external)) then
    call check_ready (self, link)                     ! finally check link task
  end if
  !.... end critical (queue_cr)
  call trace_end (itimer)
END SUBROUTINE check_nbors

!> =============================================================================
!> Unpack a message, where the MPI tag is the task id.  Use that to search
!> for the task, apply its unpack method, and check if any nbors become ready.
!> This entire operation should be threadsafe, since no other thread should be
!> working on the same message and the same patch.
!> =============================================================================
SUBROUTINE unpack (self, mesg, link)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer, optional:: link
  class(link_t), pointer:: link1, link2
  class(task_t), pointer:: task
  logical:: found, failed
  integer:: id
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher5::unpack', itimer=itimer)
  if (mesg%nbuf < 40) then
    call load_balance%unpack (mesg%buffer)
    return
  end if
  !$omp critical (unpack_cr)
  if (present(link)) then
    link1 => link
  else
    link1 => task_mesg%find_task (mesg%id)
  end if
  task => link1%task
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
    found = .false.
    link2 => self%queue
    do while (associated(link2))
      if (link2%task%id == task%id) then
        found = .true.
        exit
      end if
      link2 => link2%next_time
    end do
    write (io_unit%log,'(f12.6,2x,a,i9,l3,3x,5l1)') wallclock(), &
      'task_mesg_mod::unpack ERROR, received mpi_mesg for task with ready bit:', task%id, found, &
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
  if (mpi_mesg%debug .or. id==io%id_debug) then
    !$omp critical (debug_cr)
    write (io_unit%log,'(f12.6,2x,a,2i9,z12,2x,5l1)') wallclock(), &
      'DBG unpk: id, sender, req =', mesg%id, task%rank, mesg%req, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%external), &
      task%is_set (bits%swap_request)
    !$omp end critical (debug_cr)
  end if
  if (mpi_mesg%debug) &
    write (io_unit%log,'(f12.6,2x,a,i9,1p,e18.6)') wallclock(), 'unpk: id, time =', task%id, task%time
  if (id /= mesg%id) &
    write(io%output,*) 'unpack ERROR: wrong mesg%id', id, task%id, mesg%id
  if (.not. failed) then
    !$omp atomic
    mpi_mesg%n_unpk = mpi_mesg%n_unpk+1
    !---------------------------------------------------------------------------
    ! If the boundary+swap bits are set, this is a task that has just changed
    ! rank, and it needs to have its nbor relations re-initialized. This includes
    ! resorting (removing + re-adding) the nbor's nbor lists in rank order.
    ! FIXME: The load balancing steps should be checked for threadsafe operation
    !---------------------------------------------------------------------------
    if (task%is_set(bits%swap_request) .and. task%is_set(bits%boundary)) then
      self%na = self%na+1; self%nb = self%nb+1; self%nv = self%nv-1
      call self%init_nbors (link1)
      call task%clear (bits%swap_request+bits%ready)
      call link1%update_nbor_status (link1)
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
        write (io_unit%log,*) 'task_mesg_t%unpack: new virtal patch:', task%id
      self%nv = self%nv+1
      call self%init_nbors (link1)
      call task%clear (bits%swap_request+bits%ready)
      call link1%update_nbor_status (link1)
      call self%count_status
    end if
  end if
  !$omp end critical (unpack_cr)
  !-----------------------------------------------------------------------------
  ! As the task has now been updated, we need to check if any of the neighbors
  ! have become ready to update. If the task has just been swapped into being
  ! a boundary patch, it could possibly also be ready to update.
  !-----------------------------------------------------------------------------
  !if (verbose>1) &
  !  call link1%test_nbor_status
  if (self%method==0) &
    call check_nbors (self, link1)
  call trace%end (itimer)
END SUBROUTINE unpack

!> =============================================================================
!> =============================================================================
SUBROUTINE startup (self)
  class(task_list_t):: self
  !.............................................................................
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher5::startup')
  call init_virtual (self)
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call check_ready (self, link)
    link => link%next
  end do
  call mpi%barrier ('self%execute')
  call timer%print
  self%n_tasks = self%na
  call trace%end
END SUBROUTINE startup

!> =============================================================================
!> =============================================================================
SUBROUTINE init_virtual (self)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  !.............................................................................
  if (mpi%size==1) return
  call trace%begin ('dispatcher5::init_virtual')
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%virtual)) then
      call task%allocate_mesg
      task%wc_last = wallclock()
      call task%mesg%irecv  (task%rank, task%id)
    end if
    link => link%next
  end do
  call trace%end
END SUBROUTINE init_virtual

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_all (self)
  class(list_t):: self
  class(link_t), pointer:: link
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('dispatcer5::check_all', itimer=itimer)
  !$omp critical (list_cr)
  link => self%head
  do while (associated (link))                          ! keep going until end
    call link%task%clear (bits%ready)
    call check_ready (self, link)                       ! link ready?
    link => link%next
  end do
  !$omp end critical (list_cr)
  call trace_end (itimer)
END SUBROUTINE check_all

END MODULE
