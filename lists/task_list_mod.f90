!> ==================-===========================================================
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
!> =============================================================================
MODULE task_list_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE bits_mod
  USE trace_mod
  USE refine_mod
  USE download_mod
  USE buffered_output_mod
  USE load_balance_mod
  implicit none
  private
  type, public, extends(list_t):: task_list_t
  contains
    procedure:: init
    procedure:: execute
    procedure:: update
    procedure:: find_task
    procedure:: unpack
    procedure:: check_mpi
    procedure:: split
    procedure:: output
    procedure:: check_flags
    procedure:: print
    procedure:: print_tasks
  end type
  real, save:: check_ms=10.
  integer, save:: stalled=0, max_stalled=10000, retry_stalled=100
  integer, save:: min_nq=2**30
CONTAINS

!> =============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!> =============================================================================
SUBROUTINE init (self, name)
  class(task_list_t):: self
  character(len=*), optional:: name
  class(link_t), pointer:: link
  integer:: iostat
  integer:: debug=0
  namelist /task_list_params/ debug, max_stalled, retry_stalled, check_ms
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  call trace_begin('task_list_t%init ')
  call self%list_t%init (name)
  call mpi_mesg%init
  rewind (io%input); read(io%input, task_list_params, iostat=iostat)
  if (iostat==0 .and. io%master) then
    self%debug = debug
    if (io%master) write (*,task_list_params)
  end if
  call load_balance%init
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!> =============================================================================
SUBROUTINE execute (self)
  USE omp_timer_mod
  class(task_list_t):: self
  class(link_t), pointer:: link, rem
  real(8):: sec
  !-----------------------------------------------------------------------------
  ! Set virtual and boundary bits
  !-----------------------------------------------------------------------------
  call self%reset_status
  !-----------------------------------------------------------------------------
  ! Check the ready status of all tasks
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    call self%check_ready (link)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Start timer, and run until the task count is zero
  !-----------------------------------------------------------------------------
  call mpi%barrier ('task_list_mod::execute')
  call tic (time=sec)
  call timer%print
  self%n_tasks = self%na
  !$omp parallel
  do while (self%na > 0)
    call self%update
    if (self%na == self%n_tasks) then
      !$omp atomic
      min_nq = min(min_nq,self%nq)
    end if
  end do
  call timer%print
  call mpi_mesg%diagnostics(1)
  !$omp end parallel
  call toc ('wall time', n_update, time=sec)
  if (io%master) write (*,*) "task list finished, min_nq =", min_nq
  write (io_unit%log,*) 'thread',omp%thread,' arrived'
  flush (io_unit%log)

  !if (io%master) &
  !  write (io_unit%log,*) 'download_differ % =', 100.*download%n_differ/ &
  !    max(1,(download%n_differ+download%n_same))
  call buffered_output%close
END SUBROUTINE execute

!> =============================================================================
!> Update the state of the task list, taking the steps necessary to update the
!> head task on the ready_queue, and check for consequences
!> =============================================================================
SUBROUTINE update (self)
  class(task_list_t):: self
  !.............................................................................
  class(link_t), pointer:: head
  class(task_t), pointer:: task
  class(mesg_t), pointer:: mesg
  logical:: already_busy
  real(8):: wc
  !.............................................................................
  call trace_begin('task_list_t%update')
  !call self%check_flags
  !-----------------------------------------------------------------------------
  ! Check incoming MPI, which may free up tasks for execution
  !-----------------------------------------------------------------------------
  call self%check_mpi                                   ! check incoming MPI
  !-----------------------------------------------------------------------------
  ! Pick a task off the queue, in a very brief critical region
  !-----------------------------------------------------------------------------
  !$omp critical (queue_cr)
  head => self%queue                                    ! queue start
  if (associated(head)) then
    !---------------------------------------------------------------------------
    ! Set the head of the queue to the next task in time.  The ready bit needs
    ! to remain on for now, to prevent other threads from checking this task
    ! before it is updated.  After the update, the clearing of the ready bit
    ! allows the state of the current task to be evaluated by other threads, as
    ! well as by the current thread.  Locking should be used to allow only a
    ! single thread and task to check and change the ready bit status, at any
    ! one time.
    !---------------------------------------------------------------------------
    self%queue => head%next_time
    task => head%task                                   ! head task
    task%nq = self%nq                                   ! for info print
    !$omp atomic
    self%nq = self%nq-1                                 ! decrement queue count
    already_busy = task%is_set (bits%busy)
    if (.not.already_busy) call head%task%set (bits%busy)
  end if
  !$omp end critical (queue_cr)
  if (associated(head)) then
    if (io%verbose > 1) then
      write (io_unit%log,'(a,i4,a,i5,2x,a,2x,i4)') &
        'thread', omp%thread,' takes', task%id, trim(task%type), self%nq
      flush (io_unit%log)
    end if
    !---------------------------------------------------------------------------
    ! If the task is busy or frozen (e.g. because it has finished), return
    !---------------------------------------------------------------------------
    if (already_busy) then
      call trace_end
      return
    else if (task%is_set (bits%frozen+bits%virtual)) then
      if (task%iout == 0) call task%output(self%name)
      ! if a frozen task is set to ready, make it not ready
      if (task%is_set(bits%ready)) then
        !$omp atomic
        self%na = self%na - 1
      end if
      call trace_end
      return
    end if
    !----------------------------------------------------------------------------
    ! Download nbor info -- this may be an experiment_t procedure, but if not,
    ! the call should be answered by a solver-specific procedure (which may or
    ! may not choose to call the generic patch guard zone handler in tasks/)
    !----------------------------------------------------------------------------
    call task%dnload                                    ! download nbor info
    !----------------------------------------------------------------------------
    ! Update the task, whatever that means
    !----------------------------------------------------------------------------
    wc = wallclock()
    call task%update                                    ! update the task
    !$omp atomic write
    mpi_mesg%delay_ms = 1d3*(wallclock()-wc)
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
         self%na = self%na-1
         self%nb = self%nb-1
         self%nv = self%nv-1
         call self%count_status
      end if
      call head%send_to_vnbors                          ! send to virtual nbors
    end if
    !-----------------------------------------------------------------------------
    ! The task should now again be checked, before it can be returned to the
    ! queue, and its neighbors should also be checked for return to the queue.
    ! If the task was turned into a virtual task by the load balance, its virtual
    ! bit will prevent it being checked, but its nbors still need to be checked
    !-----------------------------------------------------------------------------
    call self%check_nbors (head)                        ! any nbors ready?
    !-----------------------------------------------------------------------------
    ! Since output of a task may be time consuming, wait with this until the task
    ! has been used for all other purposes
    !-----------------------------------------------------------------------------
    call task%output (self%name)                        ! check output
    call task%clear (bits%busy)                         ! clear busy bit
    !$omp atomic write
    stalled = 0
    !call mpi_mesg%diagnostics (0)
  else
    !$omp critical (stalled_cr)
    stalled = stalled+1
    if (stalled == retry_stalled) then
       call self%check_all
    else if (stalled == retry_stalled+1) then
       call self%check_oldest
    else
      call mpi_mesg%delay (stalled)
    end if
    if (stalled == max_stalled/2) call mpi_mesg%diagnostics (1)
    if (stalled > max_stalled) call mpi%abort ('exceeded max_stalled')
    !$omp end critical (stalled_cr)
  end if
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Find the link with a specified task id.  If there isn't any, create one!
!> If the patch already existed, it's nbor list should be OK, but the status
!> bit of it, and all its nbors, need to updated.
!> =============================================================================
FUNCTION find_task (self, id) RESULT (link)
  class(task_list_t):: self
  integer:: id
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  class(experiment_t), pointer:: task
  !.............................................................................
  call trace_begin ('task_list_mod::find_task')
  link => self%head
  do while (associated(link))
!write (io_unit%log,*) mpi%rank,' find task: id, link%id, virtual =', id, link%task%id, link%task%is_set(bits%virtual)
    patch => task2patch(link%task)
    if (patch%id == id) then
      if (io%verbose>1) &
        write (io_unit%log,*) mpi%rank,' find task: match', id, associated(patch%mem)
      flush(io_unit%log)
      call trace_end
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
  call trace_end
END FUNCTION find_task

!> =============================================================================
!> Unpack a message, where the MPI tag is the task id.  Use that to search
!> for the task, apply its unpack method, and check if any nbors become ready.
!> This entire operation should be threadsafe, since no other thread should be
!> working on the same message and the same patch.
!> =============================================================================
SUBROUTINE unpack (self, mesg)
  class(task_list_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  logical:: found
  !-----------------------------------------------------------------------------
  if (mesg%nbuf < 40) then
    call load_balance%unpack (mesg%buffer)
    return
  end if
  !$omp critical (list_cr)
  link => self%find_task (mesg%id)
  task => link%task
  !$omp end critical (list_cr)
  !-----------------------------------------------------------------------------
  ! Guard against lingering extra messages to a swapped patch
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%boundary)) then
    write (io_unit%log,'(f12.6,2x,a,i9)') wallclock(), &
      'task_list_mod::unpack ERROR, received mpi_mesg for boundary task:', task%id
    go to 9
  end if
  if (task%is_set (bits%ready)) then
    !$omp critical (list_cr)
    found = .false.
    link => self%queue
    do while (associated(link))
      if (link%task%id == task%id) then
        found = .true.
        exit
      end if
      link => link%next_time
    end do
    !$omp end critical (list_cr)
    write (io_unit%log,'(f12.6,2x,a,i9,l3)') wallclock(), &
      'task_list_mod::unpack ERROR, received mpi_mesg for task with ready bit:', task%id, found
    go to 9
  end if
  !---------------------------------------------------------------------------
  ! Unpack a patch message (which includes swapping the roles of boundary bits).
  ! Since an already existing patch may, at any one time, be under investigation
  ! by check_nbors, it must be protected by a critical region (or an OMP
  ! lock) while it is being updated here
  !---------------------------------------------------------------------------
  call task%unpack (mesg)
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
    call self%init_nbors (link)
    call task%clear (bits%swap_request+bits%ready)
    call link%update_nbor_status (link)
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
    call self%init_nbors (link)
    call task%clear (bits%swap_request+bits%ready)
    call link%update_nbor_status (link)
    call self%count_status
  end if
9 continue
  !-----------------------------------------------------------------------------
  ! As the task has now been updated, we need to check if any of the neighbors
  ! have become ready to update. If the task has just been swapped into being
  ! a boundary patch, it could possibly also be ready to update.
  !-----------------------------------------------------------------------------
  if (io%verbose>1) &
    call link%test_nbor_status
  call self%check_nbors (link)
END SUBROUTINE unpack

!> =============================================================================
!> Check for and unpack MPI messages.
!> =============================================================================
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
  call trace_begin ('task_list_mod::check_mpi', 1, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Step 1: Prune the send list -- critical regions are handled inside
  !-----------------------------------------------------------------------------
  !$omp atomic read
  nq = self%nq
  call timer%end (itimer=itimer)
  call mpi_mesg%sent_list%check_sent (nq)
  call timer%begin (itimer=itimer)
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
    call unpk_tmp%delete (mesg)
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
    call unpk_tmp%delete (mesg)
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
    call unpk_tmp%delete (mesg)
    mesg => next
  end do
  i=i+1; wc(i) = wallclock()-wc(0)
9 write (io_unit%log, '(a,6f10.4)') 'check_mpi times:', wc(1:i)-wc(0:i-1)
  call trace_end (itimer)
END SUBROUTINE check_mpi

!> =============================================================================
!> Split a patch into two (parameter) in each direction, merging the new patches
!> into the list.  Note that the link to the original patch is kept, with a "leaf"
!> flag being cleared, indicating it is not longer a "leaf" patch.
!> =============================================================================
SUBROUTINE split (self, this)
  class(task_list_t):: self
  class(link_t), pointer:: this, link
  integer:: nsplit(3)
  class(patch_t), pointer:: patch, this_patch
  class(task_t), pointer:: task
  class(link_t), pointer:: tail
  integer dir, i, ix, iy, iz
  real:: pos(3)
  !.............................................................................
  if (this%task%is_set(bits%not_leaf)) return
  call trace_begin('task_list_t%split', 1)
  nsplit = msplit
  task => this%task
  select type (task)
  class is (patch_t)
  this_patch => task
  end select
  where (this_patch%size==0 .or. this_patch%n==1)
    nsplit = 1
  end where
  do ix=1,nsplit(1)
    pos(1) = ix-(nsplit(1)+1)*0.5
    do iy=1,nsplit(2)
      pos(2) = iy-(nsplit(2)+1)*0.5
      do iz=1,nsplit(3)
        pos(3) = iz-(nsplit(3)+1)*0.5
        !------------------------------------------------------------------------
        ! Allocate a new task, 1/nsplit the size of original one, and adjust its
        ! position in each direction with "pos", in units of the new size.
        !------------------------------------------------------------------------
        allocate (patch, source=this_patch)
        patch%mem_allocated = .false.
        nullify (patch%mem)
        patch%size = this_patch%size/nsplit             ! new size
        patch%position = this_patch%position + pos*patch%size
        call patch%init                                 ! patch, for a new id
        patch%t   = this_patch%t                        ! newly allocated
        patch%dt  = this_patch%dt                       ! newly allocated
        patch%iit = this_patch%iit                      ! newly allocated
        !write (io_unit%log,*) 'mk1: id,it,it,new,new',patch%id,patch%it,this_patch%it,patch%new,this_patch%new
        !write (io_unit%log,*) 'mk2; id,it,iit',patch%id,this_patch%it, this_patch%iit
        patch%iit(1:patch%nt-1) = this_patch%it         ! fake history
        patch%out_next = patch%time                     ! force output now
        !------------------------------------------------------------------------
        ! Append links to the new patch in the patch list, sorted by quality
        !------------------------------------------------------------------------
        allocate (link)
        link%task => patch
        link%parent => this                     ! link is the link to the new patch
        call self%add_by_quality (link)         ! which is added to the patch list
        call link%make_new_nbors (link, this)   ! make a new nbors list for link
        !------------------------------------------------------------------------
        ! Download values from parent and neighbors, except at time==0.0, to
        ! allow recursive refinement there
        !------------------------------------------------------------------------
        if (this_patch%time > 0.0) call download%download_link (link, all_cells=.true.)
        !------------------------------------------------------------------------
        ! Insert the new patch at the head of the ready queue
        !------------------------------------------------------------------------
        call self%queue_by_time (link)          ! add patch to queue
        if (io%verbose>0) then
          print 1,'task_list_t%split:  patch, level, iout, time, size =', &
            patch%id, patch%level, patch%iout, patch%time, patch%size
          1 format(a,3i6,f10.6,2(2x,3f8.4))
          print 1,'task_list_t%split: parent, level, iout, time, size =', &
            this_patch%id, this_patch%level, this_patch%iout, this_patch%time, this_patch%size
        end if
      end do
    end do
  end do
  call self%check_queue
  call this_patch%set(bits%not_leaf)
  call trace_end
END SUBROUTINE split

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE output (self)
  class (task_list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call trace_begin('task_list%output')
  link => self%head
  do while (associated(link))
    call link%task%output(self%name)
    link => link%next
  end do
  call trace_end
END SUBROUTINE output

!> =============================================================================
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_tasks (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !if (io%verbose < 2) return
  if (present(label)) then
    write (io_unit%log,'(a,i4,":")') 'task_list: '//trim(label), self%n
  else
    write (io_unit%log,'(a,i4,":")') 'task_list', self%n
  end if
  link => self%head
  do while (associated(link))
    if (io%verbose < 3) then
      write (io_unit%log,'(i4,g15.6,l5)') link%task%id, link%task%time, link%task%is_set(bits%ready)
    else
      call link%print_nbors
    end if
    link => link%next
  end do
  write (io_unit%log,*) ''
END SUBROUTINE print_tasks

!> =============================================================================
!> Print a table with task id and all neighbors
!> =============================================================================
SUBROUTINE print (self, label)
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
  call self%print_queue
  call trace%end
END SUBROUTINE print

!> =============================================================================
!> Check debug flag files
!> =============================================================================
SUBROUTINE check_flags (self)
  class (task_list_t):: self
  logical:: changed
  !.............................................................................
  if (io%verbose > 0 .and. io%do_flags) then             ! avoid speed impact
    call trace_begin ('task_list_mod::check_flags')
    call icheck ('id_debug.flag', io%id_debug, changed)
    call icheck ('verbose.flag' , io%verbose , changed)
    call lcheck ('do_trace.flag', io%do_trace, changed)
    call trace_end
  end if
contains
  !=============================================================================
  logical function exists(file)
  character(len=*):: file
  inquire (file=file, exist=exists)
  end function
  !=============================================================================
  subroutine icheck (file, value, changed)
  character(len=*):: file
  integer:: value
  logical:: changed
  if (exists(trim(file))) then
    open (io_unit%flag, file=trim(file), form='formatted', status='old')
    read (io_unit%flag,*) value
    close (io_unit%flag,status='delete')
    write (io_unit%log,*) 'flag check:', file, value
    changed = .true.
  else
    changed = .false.
  end if
  end subroutine
  !=============================================================================
  subroutine lcheck (file, value, changed)
  character(len=*):: file
  logical:: value
  logical:: changed
  if (exists(trim(file))) then
    open (io_unit%flag, file=trim(file), form='formatted', status='old')
    read (io_unit%flag,*) value
    close (io_unit%flag,status='delete')
    write (io_unit%log,*) 'flag check:', file, value
    changed = .true.
  else
    changed = .false.
  end if
  end subroutine
END SUBROUTINE check_flags

END MODULE
