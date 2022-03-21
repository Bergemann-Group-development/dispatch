!===============================================================================
!> Simplified task mesg handling, relying on that threads working on task
!> updates should also handle sending and receiving task messages.  Each task
!> has a mesg instance, and a sent_list message list.  When a boundary task
!> updates, it sends a message with a unique tag to the rank nbors.  The tag
!> is made unique by the construction tag = task%id*100 + modulo(task%seq,100),
!> where task%seq starts with 1, and increases by 1 every time the task is sent.
!>
!> New boundary tasks, which do not yet have a virtual copy on the nbor ranks,
!> are given the unique tag 0, which causes the receiving task to get the new
!> task id from the message content.  Each thread maintains a thread-private
!> task template for new virtual task, with an active request waiting for
!> messages with tag==0.  When such a message is received, the message is
!> unpacked into the virtual task template and the thread then adds it to 
!> the task list, using the same procedure as was used on the sender side 
!> (list_t%add_new_link).
!>
!> Each virtual task starts updating from the owner rank by issuing a request
!> for a message with the predicted tag, which thus goes directly to the correct 
!> task.  Subsequent calls check if the message has been received, unpacks it
!> if it has, and then issues a request for the next package.
!>
!> Call hierarchy:          accesses                locks
!> --------------------------------------------------------------
!> check_mpi                entire task_list        task_list
!>   check_recv             task%mesg
!>     unpack               task%mem()              task%mem_lock 
!>   check_new              threadprivate (new)
!>     create_new_task
!> send_to_vnbors           nbor list               link
!>   check_sent             task%sent_list
!> recv_from_vnbords        nbor list               link
!>   check_recv             task%mesg               
!>     unpack               task%mem()
!>
!> Task locking:
!> -------------
!> Task execution on a single rank requires locking of the task list while it
!> is being manipulated, which for each task happens 1) when it is added to
!> the ready queue, and when it is taken off the ready queue.  Optionally, one
!> can lock access to the task time information, which is modified in task_t%
!> rotate(), and used in task_t%ahead_of(), but accessing this data atomically
!> should in principle be sufficient.  Patch data cubes need to be locked when
!> modified in the solver, and when being accessed in guard zone loading.  This
!> should normally not cause much contention, since modifications only happen
!> to a time slice that is about to be reused.  Immediately after update, the
!> time slots are rotated and the new slot is free to be used.  If significant
!> contention is measured, it may be reduced by increasing the number of time
!> slots.
!>
!> No locking is needed in connection with sending messages, since the
!> messages are added to a task private message list, which is checked later
!> by a single thread, which is the only one that can access the message.
!> 
!> The receive process needs locking, unless one decides to use thread-specific
!> sets of messages, using modulo(mpi%size).  One thread may be trying to
!> check a message for completeness, while another one is actually already
!> unpacking the message buffer.  Therefore, a message lock must be used while
!> checking and unpacking.
!>
!> While unpacking, the header_to_patch procedure modifies the task metadata,
!> which should be done under a task%lock, since the task metadat could also be 
!> used by e.g. check_ready() by other threads at the same time.   This could
!> possibly be done by instead using atomic constructs.
!>
!> Subsequently, the patch%unpack() procedure writes to a memory slot, which
!> must be protected by a task%mem_lock() for that slot, thus locking out access
!> by download_link, using the patch as a source, or by the solver (in case of
!> real tasks -- both real and virtual tasks could collide with download).
!===============================================================================
MODULE task_mesg_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE link_mod
  USE halo_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE bits_mod
  USE load_balance_mod
  USE shared_mod
  USE dll_mod
  implicit none
  private
  type, public, extends(list_t):: task_mesg_t
    integer:: method=0
  contains
    procedure:: init
    procedure:: check_mpi
    procedure:: send_to_vnbors
    procedure:: recv_from_vnbors
  end type
  integer, save:: verbose=0
  logical, save:: lock_recv=.false.
  logical, save:: await_new=.true.
  logical, save:: use_copy = .true.
  character(len=8), save:: method='irecv'
  public:: task2patch
  type(task_mesg_t), public:: task_mesg
  integer, parameter::tag_num=10
CONTAINS

!===============================================================================
!> Initialize parameters
!===============================================================================
SUBROUTINE init (self, name)
  class(task_mesg_t), target:: self
  character(len=*), optional:: name
  integer:: iostat
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  namelist /task_mesg_params/ verbose, method, lock_recv, use_copy, await_new
  character(len=120):: id = &
    '$Id$'
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) &
    return
  call trace%begin ('task_mesg_t%init')
  call trace%print_id (id)
  rewind (io%input)
  read (io%input, task_mesg_params, iostat=iostat)
  call io%namelist_warning ('task_mesg_params', iostat)
  write (io%namelist, task_mesg_params)
  !-----------------------------------------------------------------------------
  ! Choose a unique tag for each message and task ID, based on task%seq
  !-----------------------------------------------------------------------------
  mpi_mesg%uniq_mesg = .true.
  mpi_mesg%tag_type = 1
  !-----------------------------------------------------------------------------
  ! Initialize a task%mesg in all boundary and virtual tasks, and a
  ! task%sent_list in all boundary tasks
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%boundary+bits%virtual)) then
      if (.not.associated(task%mesg)) then
        select type (task)
        class is (patch_t)
          call task%allocate_mesg
        end select
      end if
      !-------------------------------------------------------------------------
      ! Make sure there's an outstanding initial recv request; this isn't
      ! strictly necessary, but helps avoid initial corner case stalls
      !-------------------------------------------------------------------------
      if (task%is_set (bits%virtual)) then
        call new_request (task)
      end if
    end if
    call link%init_verbose (verbose)
    link => link%next
  end do
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Check all virtual tasks.  Note that the accesses are grouped by thread in
!> check_recv(), so no mesg%lock is needed, and while only reading the task
!> list it does not need locking either -- a missed new task is included the
!> next time.
!===============================================================================
SUBROUTINE check_mpi (self, n)
  class(task_mesg_t):: self
  integer, optional:: n
  !.............................................................................
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  type(dll_t), pointer:: mpi_list
  class(dll_node_t), pointer:: item, next, new
  class(*), pointer:: ptr
  integer:: n_unpk
  logical:: lock
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  n_unpk = 0
  if (mpi%size > 1) then
    call trace%begin ('task_mesg_t%check_mpi', itimer=itimer)
    !---------------------------------------------------------------------------
    ! List locking is needed in AMR
    !---------------------------------------------------------------------------
    lock = shared%refine_on 
    !---------------------------------------------------------------------------
    if (use_copy) then
      if (lock)                                 call self%lock%set ('check_mpi')
      allocate (mpi_list)
      call mpi_list%init ()
      link => self%head
      do while (associated(link))
        if (associated(link%task)) then
          if (link%task%is_set(bits%virtual)) then
            allocate (item)
            item%car => link
            call mpi_list%append (item)
            ! write (stdlog,'(a,i6,1x,a,i4)') &
            !   'check_mpi:',link%task%id,'appended',mpi_list%n
          end if
        end if
        link => link%next
      end do
      item => mpi_list%head
      do while (associated(item))
        next => item%next
        ptr => item%car
        select type (ptr)
        class is (link_t)
          link => ptr
          call check_recv (self, link, n_unpk)
          ! write (stdlog,'(a,i6,1x,a,i4)') &
          !   'check_mpi:',link%task%id,'check_recv',mpi_list%n
        end select
        deallocate (item)
        item => next
      end do
      deallocate (mpi_list)
      if (lock)                               call self%lock%unset ('check_mpi')
    else
      if (lock)                                 call self%lock%set ('check_mpi')
      link => self%head
      if (lock)                               call self%lock%unset ('check_mpi')
      do while (associated(link))
        task => link%task
        if (associated(task)) then
          if (task%is_set(bits%virtual)) then
            call check_recv (self, link, n_unpk)
          end if
        end if
        if (lock)                               call self%lock%set ('check_mpi')
        link => link%next
        if (lock)                             call self%lock%unset ('check_mpi')
      end do
    end if
    call check_new (self)
    call trace%end (itimer)
  end if
  if (present(n)) then
    n = n_unpk
  end if
END SUBROUTINE check_mpi

!===============================================================================
!> Send copies of boundary tasks to all rank nbors, appending the message to
!> a message list with sent messages.   After sending, check the list for
!> completed messages, and remove and delete them.  Since only one thread at
!> a time is accessing the task here, this is thread-safe, and no locks are
!> required to handle the sending of messages.
!===============================================================================
SUBROUTINE send_to_vnbors (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  if (mpi%size == 1) return
  call trace%begin ('task_mesg_t%send_to_vnbors')
  !-----------------------------------------------------------------------------
  ! Sending with list_t%send_to_vnbors() creates a task%mesg instance
  !-----------------------------------------------------------------------------
  link%task%seq = link%task%seq + 1
  call self%list_t%send_to_vnbors (link)
  if (.not.associated(link%task%sent_list)) then
    allocate (link%task%sent_list)
    call link%task%sent_list%init
  end if
  call link%task%sent_list%add (link%task%mesg)
  !$omp atomic update
  mpi_mesg%nq_send = mpi_mesg%nq_send+1
  call check_sent (self, link%task)
  call trace%end ()
END SUBROUTINE send_to_vnbors

!===============================================================================
!> Check sent messages for this task
!===============================================================================
SUBROUTINE check_sent (self, task)
  class(task_mesg_t):: self
  class(task_t), pointer:: task
  !.............................................................................
  class(mesg_t), pointer:: mesg, next
  !-----------------------------------------------------------------------------
  ! The function mesg%test_sent() returns true when all ranks have received
  ! the message, which can then be removed and deleted, thus reclaiming memory.
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%check_sent')
  if (associated(task%sent_list)) then
    call self%lock%set ('task_mesg_t%check_sent')
    mesg => task%sent_list%head
    do while (associated(mesg))
      next => mesg%next
      if (mesg%test_sent()) then
        call task%sent_list%remove (mesg)
        call task%sent_list%delete (mesg, send=.true.)
      end if
      mesg => next
    end do
    call self%lock%unset ('task_mesg_t%check_sent')
  end if
  call trace%end()
END SUBROUTINE check_sent

!===============================================================================
!> Check the virtual nbors of a boundary task link for completed messages,
!> unpack the ones that are complete, and issue new requests for those task.
!===============================================================================
SUBROUTINE recv_from_vnbors (self, link, n_unpk)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  integer:: n_unpk
  !.............................................................................
  class(nbor_t), pointer:: nbor
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) return
  call trace%begin ('task_mesg_t%recv_from_vnbors', itimer=itimer)
  call link%lock%set ('recv_from_vnbors', itimer)
  nbor => link%nbor
  do while (associated(nbor))
    if (nbor%task%is_set (bits%virtual)) then
      call check_recv (self, nbor%link, n_unpk)
    end if
    nbor => nbor%next
  end do
  call check_new (self)
  call link%lock%unset ('recv_from_vnbors')
  call trace%end (itimer)
END SUBROUTINE recv_from_vnbors

!===============================================================================
!> Check a virtual task for completed messages.  Note that the procedure that
!> creates the task must make sure to also create a task message instance.
!===============================================================================
SUBROUTINE check_recv (self, link, n_unpk)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: n_unpk
  !.............................................................................
  integer:: seq, tag
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) return
  task => link%task
  if (lock_recv) then
    call task%mesg%lock%set ('check_recv')
  else
    if (mod(task%id,omp%nthreads) /= omp%thread) then
      return
    end if
  end if
  call trace%begin ('task_mesg_t%check_recv')
  if (verbose > 1) then
    write (io_unit%log,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
      wallclock(), omp%thread, 'check_recv', task%id, task%time, task%mesg%seq
  end if
!-----------------------------------------------------------------------------
  ! If the task has no pending request, issue one
  !-----------------------------------------------------------------------------
  if (task%mesg%req == 0) then
    call new_request (task)
    !$omp atomic update
    mpi_mesg%nq_recv = mpi_mesg%nq_recv+1
  end if
  !-----------------------------------------------------------------------------
  ! If a message is complete, unpack and issue a new request
  !-----------------------------------------------------------------------------
  if (task%mesg%is_complete()) then
    call unpack (self, task%mesg, link)
    n_unpk = n_unpk+1
    if (verbose > 1) then
      write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
        wallclock(), omp%thread, 'unpacked', task%id, task%time, n_unpk
    end if
    call new_request (task)
  end if
  if (lock_recv) then
    call task%mesg%lock%unset ('check_recv')
  end if
  call trace%end ()
END SUBROUTINE check_recv

!===============================================================================
!> Increment sequence value and issue a new request
!===============================================================================
SUBROUTINE new_request (task)
  class(task_t):: task
  !-----------------------------------------------------------------------------
  task%seq = task%seq + 1
  task%mesg%id = task%id
  task%mesg%tag = task%id*tag_num + mod(task%seq,tag_num)
  task%mesg%sender = task%rank
  call task%mesg%request
END SUBROUTINE new_request

!===============================================================================
!> Check for new virtual tasks
!===============================================================================
SUBROUTINE check_new (self)
  class(task_mesg_t):: self
  !.............................................................................
  class(mesg_t), pointer:: mesg
  class(experiment_t), save, pointer:: task => null()
  !$omp threadprivate(task)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Each thread has a thread-private new task allocated, ready to receive the
  ! next new virtual task, distinguished by having tag=MPI_ANY_TAG
  !-----------------------------------------------------------------------------
  if (.not. await_new) return
  call trace%begin ('task_mesg_t%check_new', itimer=itimer)
  if (.not.associated(task)) then
    call create_new (task)
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,1x,a,i6,i5,i8)') wallclock(), &
        'task_mesg_t%check_new: requested id, sender, tag =', &
        task%id, task%rank, task%mesg%tag
    end if
  end if
  !-----------------------------------------------------------------------------
  ! When a new task%mesg is complete, unpack it, issue a new request for the new
  ! task, and then create a new template task (also issues a new request)
  !-----------------------------------------------------------------------------
  if (task%mesg%is_complete()) then
    if (verbose > 1) then
      write (io_unit%log,'(f12.6,1x,a,2i6)') wallclock(), &
        'task_mesg_t%check_new: message id, rank =', &
        task%id, task%rank
    end if
    call unpack (self, task%mesg, task%link)
    call new_request (task)
    nullify(task)
    call create_new (task)
  else
    if (verbose > 2) then
      write (io_unit%log,'(f12.6,1x,a)') wallclock(), &
        'task_mesg_t%check_new: no message'
    end if
  end if
  call trace%end (itimer)
contains
  !=============================================================================
  ! Internal subroutine to create a new virtual task, including a task%mesg,
  ! and a task sequence starting at 1.
  !=============================================================================
  subroutine create_new (task)
    class(experiment_t), pointer:: task
    class(link_t), pointer:: task_link
    !---------------------------------------------------------------------------
    call trace%begin ('task_mesg_t%create_new_task')
    allocate (task)
    allocate (halo_t:: task_link)
    task%link => task_link
    task_link%task => task
    task%id = 1
    call task%init
    task%id = 0
    call task_link%init
    call task%allocate_mesg
    if (lock_recv) then
      allocate (task%mesg%lock)
      call task%mesg%lock%init ('mesg')
    end if
    call task%set (bits%virtual)
    task%id = 0
    task%rank = -1
    task%seq = 0
    task%mesg%tag = 0
    call task%mesg%request (any=.true.)
    call trace%end ()
  end subroutine create_new
END SUBROUTINE check_new

!===============================================================================
!> Unpack a message, where the MPI tag is the task id.  Use that to search
!> for the task, apply its unpack method, and check if any nbors become ready.
!>
!> The following cases are possible       bits: virt swap suic  new  fail
!> 1) a suicide task, to be removed              1    0    1     F    F
!> 2) a swap request on a virtual task           1    1    0     F    F
!> 3) a new virtual task                         1    0    0     T    F
!> 4) a normal update of a virtual task          1    0    0     F    F
!> 5) a belated, unwanted message                x    x    x     x    T
!>
!> Note that this procedure is sometimes handled exclusively by the master
!> thread, hence it must in principle make sure to properly lock items that it
!> modiefies.  However, it generally handles only virtual tasks, or task that
!> come in as virtual tasks and are transformed into boundary tasks, so it is
!> OK to modify task memory and generate new nbor lists, etc.
!>
!> This entire operation should be threadsafe, since no other thread should be
!> working on the same message and the same patch.
!===============================================================================
SUBROUTINE unpack (self, mesg, link)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer:: link
  !.............................................................................
  class(task_t), pointer:: task
  logical:: found, bonafide, new
  integer:: save_bits, it
  integer, save:: itimer=0
  character(len=24):: label
  type(patch_t):: patch
  real(8):: wc
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%unpack', itimer=itimer)
  task => link%task
  !-----------------------------------------------------------------------------
  ! If the task about to receive the message doesn't yet have a task id, get it
  ! from the message itself, and flag the task as new
  !-----------------------------------------------------------------------------
  if (task%id <= 0) then
    task%id = patch%unpack_id (mesg)
    task%rank = mesg%sender
    new = .true.
  else
    new = .false.
  end if
  !-----------------------------------------------------------------------------
  bonafide = .true.
  if (verbose > 1) &
    write (io_unit%log,'(a,i6,2x,6l1)') 'unpack: id, bits IBVRES =', mesg%id, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%ready), &
      task%is_set (bits%remove), &
      task%is_set (bits%swap_request)
  !-----------------------------------------------------------------------------
  ! Guard against lingering extra messages to a swapped patch
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%ready)) then
    write (io_unit%log,'(f12.6,2x,a,i9,l3,3x,5l1)') wallclock(), &
     'task_mesg_mod::unpack ERROR, received mpi_mesg for task with ready bit:',&
      task%id, found, &
      task%is_set(bits%internal), &
      task%is_set(bits%boundary), &
      task%is_set(bits%virtual), &
      task%is_set(bits%external), &
      task%is_set(bits%swap_request)
    bonafide = .false.
  end if
  !-----------------------------------------------------------------------------
  ! Measure the update cadence for this virtual task
  !-----------------------------------------------------------------------------
  wc = wallclock()
  task%update_cadence = wc - task%update_last
  task%update_last = wc
  !-----------------------------------------------------------------------------
  ! Unpack a patch message (which includes swapping the roles of boundary bits).
  !-----------------------------------------------------------------------------
  save_bits = iand (task%status, bits%init_nbors+bits%update_status)
  call task%unpack (mesg)
  task%status = ior (task%status, save_bits)
  !$omp atomic
  mpi_mesg%n_unpk = mpi_mesg%n_unpk+1
  !-----------------------------------------------------------------------------
  if (io%task_logging > 1) then
    write (label,'(a,i4,i8)') 'unpack ', mesg%sender, mesg%seq
    call task%log (label)
  end if
  !-----------------------------------------------------------------------------
  ! Since an already existing patch may, at any one time, be under investigation
  ! by check_nbors, it must be protected by a critical region (or an OMP
  ! lock) while it is being updated here
  !-----------------------------------------------------------------------------
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,3x,a,i6,f12.6,i3,2x,7l1,i8,10i3)') &
      wallclock(), 'task_mesg_t%unpack: task%id, time, IBVrRSN, tag =', &
      task%id, task%time, task%it, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%ready), &
      task%is_set (bits%remove), &
      task%is_set (bits%swap_request), &
      new, task%mesg%tag, task%iit
    flush (io_unit%mpi)
  end if
  !-----------------------------------------------------------------------------
  ! Logging and verbose output
  !-----------------------------------------------------------------------------
  if (io%log_sent > 0) then
    !$omp critical (log_sent_cr)
    call mpi_mesg%log_files()
    write (io_unit%sent,'(f16.6,i4,2x,a,i6,f16.6,8i5)') wallclock(), &
      omp%thread, 'unp', task%id, task%time, task%rank
    flush (io_unit%sent)
    !$omp end critical (log_sent_cr)
  end if
  if (mpi_mesg%debug .or. task%id==io%id_debug) then
    !$omp critical (debug_cr)
    write (io_unit%log,'(f12.6,2x,a,2i9,z12,2x,5l1)') wallclock(), &
      'DBG unpk: id, sender, req =', task%id, task%rank, mesg%req, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%external), &
      task%is_set (bits%swap_request)
    !$omp end critical (debug_cr)
  end if
  if (mpi_mesg%debug) &
    write (io_unit%log,'(f12.6,2x,a,i9,1p,e18.6)') wallclock(), &
      'unpk: id, time =', task%id, task%time
  !-----------------------------------------------------------------------------
  if (bonafide) then
    if (self%verbose > 1) write(io_unit%log,*) 'bonafide', new, &
      task%is_set(bits%init_nbors), task%is_set(bits%remove), &
      task%is_set(bits%swap_request)
    if (task%is_set(bits%init_nbors)) then
      !-------------------------------------------------------------------------
      ! If the task has the init_nbors bit set, regenerate the nbor list.
      ! Task sent has boundary tasks have had their bits changed to virtual
      ! in the unpack routine, so init_nbors() should clear the bits%init_nbors.
      !-------------------------------------------------------------------------
      if (self%verbose > 1) write(io_unit%log,*) task%id, 'init_nbors', &
        task%is_set(bits%boundary)
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vinitnbors', task%id, task%time
      end if
      call self%init_nbors (link)
      call task%clear (bits%init_nbors)
    end if
    if (task%is_set(bits%update_status)) then
      !-------------------------------------------------------------------------
      ! If the task has the update_status bit set, re-evaluate status bits
      !-------------------------------------------------------------------------
      if (self%verbose > 1) write(io_unit%log,*) task%id, 'update_nbors'
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vupnbors', task%id, task%time
      end if
      call self%update_status (link)
    end if
    if (task%is_set(bits%remove)) then
      if (self%verbose > 0) write(io_unit%log,*) task%id, 'remove'
      !-------------------------------------------------------------------------
      ! 1) Suicide note, remove the task and delete the message
      !-------------------------------------------------------------------------
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vremove', task%id, task%time
      end if
      call self%remove_and_reset (link)
    else if (task%is_set(bits%swap_request)) then
      if (self%verbose > 0) write(io_unit%log,*) task%id, 'swap'
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vswap', task%id, task%time
      end if
      !-------------------------------------------------------------------------
      ! 2) swap_request bit set; the task has just changed owner, and the task
      !    should become a boundary task on this rank.  Indeed, since it was
      !    sent as a virtual task, the flipping should have been done in the
      !    unpack procedure already.   Delete the message.
      !-------------------------------------------------------------------------
      call self%lock%set ('task_mesg_t%unpack')
      if (load_balance%is_on) &
        call load_balance%lock%set ('task_mesg_t%unpack')
      call self%swap (link, .false.)
      if (load_balance%is_on) &
        call load_balance%lock%unset ('task_mesg_t%unpack')
      call self%lock%unset ('task_mesg_t%unpack')
    else if (new) then
      if (self%verbose > 0) write(io_unit%log,*) task%id, 'new', task%it
      !-------------------------------------------------------------------------
      ! 3) the task is a new virtual task.  It has the init_nbors bit set, and
      !    since it is a virtual task that bit should be cleared the next time
      !    the task is updated, coming back not as new task but as a virtual
      !    task update.  But first, issue a request for the next message.
      !-------------------------------------------------------------------------
      it = task%it
      !-------------------------------------------------------------------------
      ! Let all time slots refer to the mem slot the new task data is in
      !-------------------------------------------------------------------------
      task%t(:) = task%time
      task%dt(:) = task%dtime
      task%iit(:) = it
      !-------------------------------------------------------------------------
      ! Add the new task to the task list, and give it an nbor list, etc
      !-------------------------------------------------------------------------
      task%key = task%hash_key()                    ! not transferred!
      call self%add_by_level (link, count=.true.)   ! add to task list
      call self%add_new_link (link)                 ! default actions
      call self%find_parent (link)                  ! find L-1 parent
      task%mesg%seq = task%mesg%seq-1               ! adjust sequence number
      if (verbose > 0) then
        if (associated(link%parent)) then
          write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,3i6)') &
            wallclock(), omp%thread, 'vnew', task%id, task%time, task%istep, &
            task%mesg%seq, link%parent%task%id
        else
          write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,3i6)') &
            wallclock(), omp%thread, 'vnew', task%id, task%time, task%istep, &
            task%mesg%seq
        end if
      end if
    else
      if (self%verbose > 1) write(io_unit%log,*) task%id, 'update'
      !-------------------------------------------------------------------------
      ! 4) task was updated, so check nbors for ready tasks, after issuing a
      !    request for the next message.
      !-------------------------------------------------------------------------
      call self%check_nbors (link)
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vupdate', task%id, task%time
      end if
    end if
  else
    if (self%verbose > 1) write(io_unit%log,*) task%id, 'crap'
    !---------------------------------------------------------------------------
    ! 5) not bonafide, so ignore
    !---------------------------------------------------------------------------
    if (verbose > 0) then
      write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
        wallclock(), omp%thread, 'vcrap', task%id, task%time
    end if
  end if
  call self%check_hash_table ('unpack')
  call trace%end (itimer)
END SUBROUTINE unpack

END MODULE task_mesg_mod
