!====================-===========================================================
!> $Id$
!> A task_list_t has the distinction that it must be able to call an initial
!> condition routine, which only exists in experiment_t tasks.  Procedures in
!> the module use these foreign procedures:
!>
!>    task_t%init         : to initialize a task
!>    task_t%update       : to update a task
!>    task_t%send_to      : send data to foreign task
!>    list_t%remove_task  : remove task from list
!>    list_t%remove_head  : remove head from list
!>
!===============================================================================
MODULE task_list_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE link_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE bits_mod
  USE trace_mod
  USE refine_mod
  USE load_balance_mod
  USE data_io_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
    logical:: syncing=.false.
    logical:: recv_active=.false.
    real(8):: sync_time=-1.0, sync_next=-1.0
    real(8):: job_seconds=1d30
  contains
    procedure:: init
    procedure:: updater
    procedure:: find_task
    procedure:: unpack
    procedure:: check_mpi
    procedure:: init_virtual
    procedure:: check_virtual
    procedure:: check_active
    procedure:: output
    procedure:: print => print_table
  end type
  real, save:: check_ms=10.
  integer, save:: stalled=0, max_stalled=10000, retry_stalled=100
  integer, save:: min_nq=2**30
  logical, save:: do_delay=.false.
  logical, save:: queue_unpack=.true.
  logical, save:: track_active=.false.
  real(8), save:: sync_time=0.0
  type(task_list_t), save:: virtual_list
  !$omp threadprivate (virtual_list)
CONTAINS

!===============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!===============================================================================
SUBROUTINE init (self, name)
  class(task_list_t):: self
  character(len=*), optional:: name
  class(link_t), pointer:: link
  integer:: iostat
  integer:: debug=0
  real(8):: job_seconds=1d30
  logical, save:: recv_active=.true.
  namelist /task_list_params/ debug, max_stalled, retry_stalled, check_ms, &
    do_delay, recv_active, queue_unpack, sync_time, job_seconds, track_active
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  call trace%begin('task_list_t%init')
  call self%list_t%init (name)
  call mpi_mesg%init
  rewind (io%input); read(io%input, task_list_params, iostat=iostat)
  if (io%master) then
    self%debug = debug
    if (io%master) write (*,task_list_params)
  end if
  if (sync_time > 0.0) &
    self%sync_next = sync_time
  self%sync_time = sync_time
  self%job_seconds = job_seconds
  self%recv_active = recv_active
  write (io_unit%mpi,*) 'task_list%init: next_sync =', self%sync_next, self%sync_time
  call load_balance%init
  call trace_end
END SUBROUTINE init

SUBROUTINE updater (self, head)
  implicit none
  class(task_list_t):: self
  class(link_t), pointer:: head
  class(task_t), pointer:: task
  real(8):: wc
  integer, save:: itimer=0
  !---------------------------------------------------------------------------
  ! Check for flag files.  If io%out_time is set, set a new out_next.
  ! If task%time is slightly smaller than a multiple of out_time, then
  ! out_next will becoe that multiple.
  !---------------------------------------------------------------------------
  call trace%begin('task_list_t%updater', itimer=itimer)
  task => head%task
  call io%check_flags
  if (task%id==io%id_track) then
    task%track =  .not.task%track
    io%id_track = 0
  end if
  !---------------------------------------------------------------------------
  ! If the task is frozen (e.g. because it has finished), return
  !---------------------------------------------------------------------------
  if (task%is_set (bits%frozen+bits%virtual)) then
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
  call task%dnload                                    ! download nbor info
  !----------------------------------------------------------------------------
  ! Update the task, whatever that means (may include call to task%output)
  !----------------------------------------------------------------------------
  if (mpi_mesg%delay_ms==0.0) then
    wc = wallclock()
  end if
  task%rotated = .false.
  call task%update                                    ! update the task
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
  !----------------------------------------------------------------------------
  ! The rotate procedure in the generic patch_t module is responsible for
  ! updating the time, and rotating the memory slots where information about
  ! previous time steps are saved.
  !----------------------------------------------------------------------------
  call task%rotate
  select type (task)
  class is (patch_t)
    !$omp atomic
    n_update = n_update + product(task%n)
  end select
  call task%info (self%nq, self%na)                   ! print info on stdout
  !-----------------------------------------------------------------------------
  ! Tasks that have finished are subtracted from the task list count but are
  ! not removed, since their data may be needed by other tasks (including on
  ! other ranks).  They are, however, not added back onto the ready queue.
  !-----------------------------------------------------------------------------
  if (task%has_finished()) then  ! finished:
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
    call head%send_to_vnbors                          ! send to virtual nbors
  end if
  call trace%end (itimer)
END SUBROUTINE updater

!===============================================================================
!> Find the link with a specified task id.  If there isn't any, create one!
!> If the patch already existed, it's nbor list should be OK, but the status
!> bit of it, and all its nbors, need to updated.
!===============================================================================
FUNCTION find_task (self, id) RESULT (link)
  class(task_list_t):: self
  integer:: id
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  class(experiment_t), pointer:: task
  !.............................................................................
  call trace%begin ('task_list_t%find_task')
  link => self%head
  do while (associated(link))
!write (io_unit%log,*) mpi%rank,' find task: id, link%id, virtual =', id, link%task%id, link%task%is_set(bits%virtual)
    patch => task2patch(link%task)
    if (patch%id == id) then
      if (io%verbose>1) &
        write (io_unit%log,*) mpi%rank,' find task: match', id, associated(patch%mem)
      flush(io_unit%log)
      call trace%end
      return
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! No task found, which means it's a new virtual task that needs to be created
  !-----------------------------------------------------------------------------
  allocate (link)
  allocate (task)
  link%task => task
  task%id = id
  task%box = self%size
  call task%patch_t%init
  if (io%verbose>1) &
    write (io_unit%log,*) mpi%rank,' find task: new task', id, associated(task%mem)
  flush(io_unit%log)
  !-----------------------------------------------------------------------------
  ! Append to the task_list.  The calling procedure must create an nbor list,
  ! and must also update the nbor list of the nbors
  !-----------------------------------------------------------------------------
  call task%set (bits%virtual)
  call self%append_link (link)
  call trace%end
END FUNCTION find_task

!===============================================================================
!> Unpack a message, where the MPI tag is the task id.  Use that to search
!> for the task, apply its unpack method, and check if any nbors become ready.
!> This entire operation should be threadsafe, since no other thread should be
!> working on the same message and the same patch.
!===============================================================================
SUBROUTINE unpack (self, mesg, link)
  class(task_list_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer, optional:: link
  class(link_t), pointer:: link1, link2
  class(task_t), pointer:: task
  logical:: found, failed
  integer:: id
  !-----------------------------------------------------------------------------
  if (mesg%nbuf < 40) then
    call load_balance%unpack (mesg%buffer)
    return
  end if
  !$omp critical (unpack_cr)
  if (present(link)) then
    link1 => link
  else
    link1 => self%find_task (mesg%id)
  end if
  task => link1%task
  failed = .false.
  !-----------------------------------------------------------------------------
  ! Guard against lingering extra messages to a swapped patch
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%boundary)) then
    write (io_unit%log,'(f12.6,2x,a,i9,3x,5l1)') wallclock(), &
      'task_list_mod::unpack ERROR, received mpi_mesg for boundary task:', task%id, &
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
      'task_list_mod::unpack ERROR, received mpi_mesg for task with ready bit:', task%id, found, &
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
  if (mpi_mesg%debug) then
    !$omp critical (debug_cr)
    write (io_unit%mpi,'(f12.6,2x,a,i9,1p,e18.6,2x,a,i5,2x,5l1)') &
      wallclock(), 'recv: id, time =', mesg%id, task%time, 'fr', task%rank, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%external), &
      task%is_set (bits%swap_request)
    !$omp end critical (debug_cr)
  end if
  if (mpi_mesg%debug) &
    write (io_unit%log,'(f12.6,2x,a,i9,1p,e18.6)') wallclock(), 'unpk: id, time =', task%id, task%time
  if (id /= mesg%id) print *,'unpack ERROR: wrong mesg%id', id, task%id, mesg%id
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
      if (io%verbose>1) &
        print '(f12.6,2x,a,i6,a,i9,a,i6)', &
        wallclock(),'LB: rank',mpi%rank,' given patch',task%id,' by',mesg%sender
      if (io%verbose>0) &
        write (io_unit%log,*) 'swapped virtual to boundary:', task%id
    !---------------------------------------------------------------------------
    ! If the link has no nbors it is a newly created virtual task. Does it need
    ! an nbor list?  At least we can use the nbor list to check that the link is
    ! in its nbors nbor lists.  A new virtual task (where no task existed) means
    ! that some nbor of it has changed from internal to boundary, which will be
    ! checked by the test_nbor_status call below, but only if an nbor list exists.
    !---------------------------------------------------------------------------
    !else if (.not.associated(link%nbor)) then
    else if (task%is_set(bits%swap_request) .and. task%is_set(bits%virtual)) then
      if (io%verbose>0) &
        write (io_unit%log,*) 'new virtal patch:', task%id
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
  !if (io%verbose>1) &
  !  call link1%test_nbor_status
  !call self%check_nbors (link1)
END SUBROUTINE unpack

!===============================================================================
!> Check for and unpack MPI messages.
!===============================================================================
SUBROUTINE check_mpi (self)
  class(task_list_t):: self
  class(mesg_t), pointer:: mesg, next
  type(mesg_list_t):: unpk_tmp, recv_tmp
  integer, save:: itimer=0
  integer:: i, n, max_recv, nq
  real(8):: wc(0:10)=0.0
  !.............................................................................
  if (mpi%size <= 1) return
  unpk_tmp%name = 'unpk_tmp'
  recv_tmp%name = 'recv_tmp'
  wc(0) = wallclock()
  !-----------------------------------------------------------------------------
  ! Step 1: Prune the send list -- critical regions are handled inside
  !-----------------------------------------------------------------------------
  !$omp atomic read
  nq = self%nq
  call mpi_mesg%sent_list%check_sent (nq)
  if (self%recv_active) then
    if (mpi_mesg%recv_priv) then
      call self%check_virtual
    else
      call self%check_active
    end if
    return
  end if
  if (mpi_mesg%recv_priv) then
    call check_priv (self)
    return
  end if
  call trace%begin ('task_list_t%check_mpi', 1, itimer=itimer)
  i=1; wc(i) = wallclock()-wc(0)
  if (wc(i) > check_ms*1e-3)  go to 9
  !-----------------------------------------------------------------------------
  ! Step 2: probe for new messages and check for completeness. Put messages on
  ! one of two thread-local lists, depending on if complete or not.  The probing
  ! can be done outside (recv_cr), but it may be necessary to serialize on systems
  ! that do not support fully multi-threaded MPI
  !-----------------------------------------------------------------------------
  call unpk_tmp%reset
  call recv_tmp%reset
  call timer%end (itimer=itimer)
  call mpi_mesg%get (mesg)
  call timer%begin (itimer=itimer)
  n = 0
  do while (associated(mesg))
    if (mesg%is_complete()) then
      call unpk_tmp%add (mesg)
    else
      call recv_tmp%add (mesg)
    end if
    n = n+1
    if (n > mpi_mesg%max_probe) exit
    call timer%end (itimer=itimer)
    call mpi_mesg%get (mesg)
    call timer%begin (itimer=itimer)
  end do
  !-----------------------------------------------------------------------------
  ! Move pending receive messages to the recv_list.  This part must be in a
  ! (recv_cr) critical region, but is very fast
  !-----------------------------------------------------------------------------
  !$omp critical (recv_cr)
  mesg => recv_tmp%head
  do while (associated(mesg))
    next => mesg%next
    call recv_tmp%remove (mesg)
    call mpi_mesg%recv_list%add (mesg)
    mesg => next
  end do
  !$omp end critical (recv_cr)
  !-----------------------------------------------------------------------------
  ! Now we can unpack messages from the thread-local list, outside of (recv_cr)
  !-----------------------------------------------------------------------------
  mesg => unpk_tmp%head
  do while (associated(mesg))
    next => mesg%next
    call timer%end (itimer=itimer)
    call self%unpack (mesg)
    call timer%begin (itimer=itimer)
    call unpk_tmp%remove (mesg)
    call unpk_tmp%delete (mesg, .false.)
    mesg => next
  end do
  i=i+1; wc(i) = wallclock()-wc(0)
  if (wc(i) > check_ms*1e-3)  go to 9
  !-----------------------------------------------------------------------------
  ! Step 3: if there are more than max_recv messages in recv_list, then wait
  ! for and unpack enough messages to reach the limit.  If task_list%nq is less
  ! then min_nq, reduce the max_recv to zero, forcing a wait on all receives.
  ! By splitting the operation into one part that moves the messages to a
  ! thead-local list, and a second part that does the wait and unpack, we
  ! minimize the blocking of other threads.
  !-----------------------------------------------------------------------------
  max_recv = merge(mpi_mesg%max_recv,0,self%nq >= mpi_mesg%min_nq)
  call unpk_tmp%reset
  !$omp atomic read
  n = mpi_mesg%recv_list%n
  if (n > max_recv) then
    !$omp critical (recv_cr)
    mesg => mpi_mesg%recv_list%head
    do while (mpi_mesg%recv_list%n > max_recv)
      next => mesg%next
      call mpi_mesg%recv_list%remove (mesg)
      call unpk_tmp%add (mesg)
      mesg => next
    end do
    !$omp end critical (recv_cr)
  end if
  mesg => unpk_tmp%head
  do while (associated(mesg))
    next => mesg%next
    call timer%end (itimer=itimer)
    call mesg%wait_for_completion ()
    call self%unpack (mesg)
    call timer%begin (itimer=itimer)
    call unpk_tmp%remove (mesg)
    call unpk_tmp%delete (mesg, .false.)
    mesg => next
  end do
  i=i+1; wc(i) = wallclock()-wc(0)
  if (wc(i) > check_ms*1e-3) go to 9
  !-----------------------------------------------------------------------------
  ! Step 4: check the rest of the recv_list for completed messages.  The (fast)
  ! traversal is protected inside (recv_cr), but unpacking from the local list
  ! can be done outside.
  !-----------------------------------------------------------------------------
  call unpk_tmp%reset
  !$omp critical (recv_cr)
  mesg => mpi_mesg%recv_list%head
  do while (associated(mesg))
    next => mesg%next
    if (mesg%is_complete()) then
      call mpi_mesg%recv_list%remove (mesg)
      call unpk_tmp%add (mesg)
    end if
    mesg => next
  end do
  !$omp end critical (recv_cr)
  mesg => unpk_tmp%head
  do while (associated(mesg))
    next => mesg%next
    call timer%end (itimer=itimer)
    call self%unpack (mesg)
    call timer%begin (itimer=itimer)
    call unpk_tmp%remove (mesg)
    call unpk_tmp%delete (mesg, .false.)
    mesg => next
  end do
  i=i+1; wc(i) = wallclock()-wc(0)
9 continue
  !write (io_unit%log, '(a,6f10.4)') 'check_mpi times:', wc(1:i)-wc(0:i-1)
  call trace%end (itimer)
END SUBROUTINE check_mpi

!===============================================================================
!> Check for and unpack MPI messages, using threadprivate recv_list and
!> unpk_list
!===============================================================================
SUBROUTINE check_priv (self)
  class(task_list_t):: self
  class(mesg_t), pointer:: mesg, next
  integer, save:: itimer=0
  integer:: n
  !.............................................................................
  if (mpi%size <= 1) return
  call trace%begin ('task_list_t%check_priv', 1, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Step 1: Add incoming messages to either recv_list or unpk_list
  !-----------------------------------------------------------------------------
  call mpi_mesg%get (mesg)
  n = 0
  do while (associated(mesg))
    if (mesg%is_complete()) then
      call unpk_list%add (mesg)
    else
      call recv_list%add (mesg)
    end if
    n = n+1
    if (n > mpi_mesg%max_probe) exit
    call mpi_mesg%get (mesg)
  end do
  !-----------------------------------------------------------------------------
  ! Step 2: check the recv_list for completed messages
  !-----------------------------------------------------------------------------
  mesg => recv_list%head
  do while (associated(mesg))
    next => mesg%next
    if (mesg%is_complete()) then
      call recv_list%remove (mesg)
      call unpk_list%add (mesg)
    end if
    mesg => next
  end do
  !-----------------------------------------------------------------------------
  ! Step 3: unpack messags in unpk_list
  !-----------------------------------------------------------------------------
  mesg => unpk_list%head
  do while (associated(mesg))
    next => mesg%next
    call self%unpack (mesg)
    call unpk_list%remove (mesg)
    call unpk_list%delete (mesg, .false.)
    mesg => next
  end do
  call trace%end (itimer)
END SUBROUTINE check_priv

!===============================================================================
!> Run through all virtual patches, count them, and create separate threadprivate
!> lists that together include all of them.
!===============================================================================
SUBROUTINE init_virtual (self)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: nvirt, i, i1, i2
  integer, save:: itimer=0
  !.............................................................................
  if (mpi%size==1) return
  call trace%begin ('task_list_t%init_virtual', itimer=itimer)
  !---------------------------------------------------------------------------
  ! Count the number of virtual tasks, and assign a number per thread
  !---------------------------------------------------------------------------
  nvirt = 0
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%virtual)) then
      nvirt = nvirt+1
    end if
    link => link%next
  end do
  nvirt = nvirt/omp%nthreads + 1
  !---------------------------------------------------------------------------
  ! Each thread copies its share of virtual patches to its private list, and
  ! initializes its private virtual task list.  Using this with succesive
  ! passage of each thread through a critical region worked OK with gfortran,
  ! but gives serious compiler related problems with ifort (17.0.1).  Below
  ! it has been rewriten to work with computed intervals for each thread.
  !---------------------------------------------------------------------------
  !$omp parallel shared(nvirt,mpi,bits,self,io) private(link,i,i1,i2,task) default(none)
  i1 = omp_get_thread_num()*nvirt
  i2 = i1 + nvirt
  i = 1
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%virtual)) then
      if (i >= i1 .and. i < i2) then
        call virtual_list%append (task, nbor=link%nbor, needed=link%needed, needs_me=link%needs_me)
        call task%allocate_mesg
        task%wc_last = wallclock()
        call task%mesg%irecv  (task%rank, task%id)
        if (io%verbose>0) &
          write (io_unit%log,*) 'init_virtual:', omp_get_thread_num(), task%id
        !if (mpi%rank==0) print *,mpi%rank,'MK',omp_get_thread_num(),i,task%id,'irec'
      end if
      i = i+1
    end if
    link => link%next
  end do
  !$omp end parallel
  call trace%end (itimer)
END SUBROUTINE init_virtual

!===============================================================================
!> Run through all virtual patches, unpack completed receives, and issue new
!> receive requests.  The virtual_list is threadprivate, so no critical region
!> is needed.
!===============================================================================
SUBROUTINE check_virtual (self)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: nunpk
  real(8):: wc
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('task_list_t%check_virtual', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Each thread runs through its list of virtual tasks, unpacks completed
  ! messages, and issues a new request
  !-----------------------------------------------------------------------------
  nunpk = 0
  link => virtual_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
    if (task%mesg%test_recv()) then
      if (io%verbose > 0) &
        write (io_unit%log,*) wallclock(),' check_virtual hit ', task%id
      call self%unpack (task%mesg, link=link)
      !$omp atomic
      mpi_mesg%n_recv = mpi_mesg%n_recv+1
      !$omp atomic
      timer%bytes_recv = timer%bytes_recv + 4.0_8*task%mesg%nbuf
      !$omp atomic
      timer%n_recv = timer%n_recv + 1_8
      wc = wallclock()
      if (wc-task%wc_last > 10.) print &
        '("WARNING: virtual patch not updated in",f5.1," sec")', wc
      task%wc_last = wc
      call task%mesg%irecv (task%rank, task%id)
      nunpk = nunpk+1
    else if (io%verbose > 2) then
      write (io_unit%log,*) wallclock(),' check_virtual miss', task%id
    end if
    end select
    link => link%next
  end do
  if (io%verbose>1) write (io_unit%log,*) 'check_virtual', nunpk
  call trace%end (itimer)
END SUBROUTINE check_virtual

!===============================================================================
!> Run through all virtual patches, and allocate message data types if not
!> allocated.  If an old request exists, test if it has completed.  If yes,
!> unpack and issue a new request.  If no request exists, issue a request.
!===============================================================================
SUBROUTINE check_active (self)
  class(task_list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: nvirt, nunpk
  real(8):: wc
  integer, save:: itimer=0
  type(mesg_list_t):: unpk_tmp
  class(mesg_t), pointer:: mesg
  !.............................................................................
  call trace%begin ('task_list_t%check_active', itimer=itimer)
  !$omp critical (active_cr)
  nvirt = 0
  nunpk = 0
  link => self%head
  if (queue_unpack) call unpk_tmp%reset
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%virtual)) then
      nvirt = nvirt+1
      select type (task)
      class is (patch_t)
      if (.not.associated(task%mesg)) then
        call task%allocate_mesg
        task%wc_last = wallclock()
      end if
      if (task%mesg%req == 0) then
        call task%mesg%irecv  (task%rank, task%id)
        if (mpi_mesg%debug) &
          write (io_unit%log,'(f12.6,2x,a,3i9)') wallclock(), 'recv: id, sender, req =', task%id, task%rank, task%mesg%req
      end if
      !-------------------------------------------------------------------------
      ! If a mesg was recvd, then either unpack directly and issue a new request,
      ! or put the mesg on a thread private unpack list, to unpack outside the
      ! critical region
      !-------------------------------------------------------------------------
      if (task%mesg%test_recv()) then
        if (queue_unpack) then
          call unpk_tmp%add (task%mesg)
        else
          call self%unpack (task%mesg)
          if (mpi_mesg%debug) &
            write (io_unit%log,'(f12.6,2x,a,3i9)') wallclock(), 'unpk: id, sender, req =', task%mesg%id, task%mesg%sender, task%mesg%req
          call task%mesg%irecv (task%rank, task%id)
          if (mpi_mesg%debug) &
            write (io_unit%log,'(f12.6,2x,a,3i9)') wallclock(), 'recv: id, sender, req =', task%id, task%rank, task%mesg%req
          if (task%rank /= task%mesg%sender) error stop 'rank ne sender'
          if (task%id   /= task%mesg%id)     error stop 'task%id ne mesg%id'
        end if
      end if
      end select
    end if
    link => link%next
  end do
  !$omp end critical (active_cr)
  if (queue_unpack) then
    mesg => unpk_tmp%head
    do while (associated(mesg))
      call self%unpack (mesg)
      if (mpi_mesg%debug) &
        write (io_unit%log,'(f12.6,2x,a,3i9)') wallclock(), 'unpk: id, time =', mesg%id, mesg%sender, mesg%req
      call unpk_tmp%remove (mesg)
      call mesg%irecv (mesg%sender, mesg%id)
      if (mpi_mesg%debug) &
        write (io_unit%log,'(f12.6,2x,a,3i9)') wallclock(), 'recv: id, time =', mesg%id, mesg%sender, mesg%req
      nunpk = nunpk+1
      mesg => mesg%next
    end do
  end if
  call trace%end (itimer)
END SUBROUTINE check_active

!===============================================================================
!> Print a table with patch id, positions, and all neighbors
!===============================================================================
SUBROUTINE output (self)
  class (task_list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call trace%begin('task_list%output')
  link => self%head
  do while (associated(link))
    !call link%task%output(self%name)
    link => link%next
  end do
  call trace%end
END SUBROUTINE output

!===============================================================================
!> Print a table with task id and all neighbors
!===============================================================================
SUBROUTINE print_table (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call trace%begin ('task_list_t%print_tasklist')
  if (io%master.and.present(label)) write (io_unit%log,*) label
  link => self%head
  do while (associated(link))
    call link%print_nbors ('')
    link => link%next
  end do
  !call self%print_queue
  call trace%end
END SUBROUTINE print_table

END MODULE task_list_mod
