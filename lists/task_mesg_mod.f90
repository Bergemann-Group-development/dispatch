!===============================================================================
!> Simplified task mesg handling, relying on that threads working on task
!> updates should also handle sending and receiving task messages.  Each task
!> has a mesg instance, and a sent_list message list.  When a boundary task
!> updates it sends a message with a unique tag to the rank nbors.  The tag
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
!> task.  Subsequent calls check if the message has been received, unpack it
!> if it has, and then issue a request for the next package.
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
!> No locking is needed in connection with sending messages, since the
!> messages are added to a task private message list, which is checked later
!> by a single thread, which is the only one that can access the message.
!> 
!> The receive process needs locking, unless one decides to use thread-specific
!> sets of messages (e.g. using modulo(mpi%size)).  One thread may be trying to
!> checking a message for completeness, while another one is actually already
!> unpacking the message buffer.  Therefore, a message lock must be used while
!> checking and unpacking.
!>
!> While unpacking, the header_to_patch procedure modifies the task metadata,
!> which should be done under a task%lock, which could also be used by e.g.
!> check_ready() by other threads at the same time.   This could possibly be
!> protected instead by atomic constructs.
!>
!> Subsequently, the patch%unpack() procedure writes to a memory slot, which
!> must be protected by a task%mem_lock() for that slot, thus locking out access
!> by download_link, using the patch as a source, or by the solver (in case of
!> real tasks -- real or virtual tasks could collide with download).
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
  USE list_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE bits_mod
  USE load_balance_mod
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
  logical, save:: use_hashtable=.true.
  character(len=8), save:: method='irecv'
  public:: task2patch
  type(task_mesg_t), public:: task_mesg
CONTAINS

!===============================================================================
!> Initialize parameters
!===============================================================================
SUBROUTINE init (self, name)
  class(task_mesg_t), target:: self
  character(len=*), optional:: name
  integer:: iostat
  class(link_t), pointer:: link
  namelist /task_mesg_params/ verbose, method
  character(len=120):: id = &
    '$Id$ lists/task_mesg_mod.f90'
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) &
    return
  call trace%begin ('task_mesg_t%init')
  call trace%print_id (id)
  rewind (io%input)
  read (io%input, task_mesg_params, iostat=iostat)
  call io%namelist_warning ('task_mesg_params', iostat)
  write (io%output, task_mesg_params)
  !-----------------------------------------------------------------------------
  ! Choose a unique tag for each message and task ID, based on task%seq
  !-----------------------------------------------------------------------------
  mpi_mesg%uniq_mesg = .true.
  mpi_mesg%tag_type = 1
  !-----------------------------------------------------------------------------
  ! Initialize task%mesg in all boundary and virtual tasks
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    if (link%task%is_set (bits%boundary+bits%virtual)) then
      associate (task => link%task)
      select type (task)
      class is (patch_t)
        if (.not.associated(task%mesg)) then
          call task%allocate_mesg
          task%mesg%sender = task%rank
          task%seq = 0
        end if
      end select
      end associate
    end if
    call link%init_verbose (verbose)
    link => link%next
  end do
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Check all virtual tasks.  Note that the accesses are grouped by thread, so
!> no mesg%lock is needed
!===============================================================================
SUBROUTINE check_mpi (self, n)
  class(task_mesg_t):: self
  integer, optional:: n
  !.............................................................................
  class(link_t), pointer:: link
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%check_mpi', itimer=itimer)
  !call self%lock%set ('check_mpi')
  link => self%head
  do while (associated(link))
    if (link%task%is_set(bits%virtual)) then
      call check_recv (self, link%task)
    end if
    link => link%next
  end do
  !call check_new(self)
  !call self%lock%unset ('check_mpi')
  call trace%end (itimer)
END SUBROUTINE check_mpi

!===============================================================================
!> Send copies of boundary tasks to all rank nbors that should have a copy,
!> appending the message to a message list with sent messages.   After sending,
!> check the list for completed messages, and remove and delete them.  Since
!> only one thread at a time is accessing the task here, this is thread-safe,
!> and no locks are required to handle the sending of messages.
!===============================================================================
SUBROUTINE send_to_vnbors (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(task_t), pointer:: task
  class(mesg_t), pointer:: mesg, next
  class(mesg_list_t), pointer:: sent_list
  integer:: req
  !-----------------------------------------------------------------------------
  ! Sending the task creates a task%mesg instance, which is appended to
  ! the task%sent_list -- created on the fly if not existing
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%send_to_vnbors')
  task => link%task
  task%seq = task%seq + 1
write (io_unit%log,*) link%task%seq, 'calling list_t%send_to_vnbors'  
  call self%list_t%send_to_vnbors (link)
write (io_unit%log,*) '      sent', task%id, task%mesg%id, task%mesg%tag 
  if (associated(task%sent_list)) then
    sent_list => task%sent_list
  else
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,1x,a,2i6)') wallclock(), &
        'task_mesg_t%send_to_vnbors: sent_list added to task', link%task%id
    end if
    allocate (sent_list)
    call sent_list%init
    task%sent_list => sent_list
  end if
  call sent_list%lock%set
  call sent_list%add (link%task%mesg)
  !$omp atomic update
  mpi_mesg%nq_send = mpi_mesg%nq_send+1
  call check_sent (self, task)
  call sent_list%lock%unset
  call trace%end ()
END SUBROUTINE send_to_vnbors

!===============================================================================
!> Check recv messages on all virtual task nbors
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
    mesg => task%sent_list%head
    do while (associated(mesg))
      next => mesg%next
      if (mesg%test_sent()) then
write (io_unit%log,*) '   deleted', task%id, mesg%id, mesg%tag, task%sent_list%n
flush (io_unit%log)
!write (io_unit%log,*) mesg%id, mesg%nreq, mesg%req, task%sent_list%n, mpi_mesg%nq_send, 'delete'
        call task%sent_list%remove (mesg)
        call task%sent_list%delete (mesg, send=.true.)
      else
write (io_unit%log,*) '      kept', task%id, mesg%id, mesg%tag
flush (io_unit%log)
!write (io_unit%log,*) mesg%id, mesg%nreq, mesg%req, task%sent_list%n, mpi_mesg%nq_send, 'keep'
      end if
      mesg => next
    end do
  end if
  call trace%end()
END SUBROUTINE check_sent

!===============================================================================
!> Check the virtual nbors of a boundary task link for completed messages,
!> unpack the ones that are complete, and issue new requests for those task.
!===============================================================================
SUBROUTINE recv_from_vnbors (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%recv_from_vnbors')
  call link%lock%set ('recv_from_vnbors')
  nbor => link%nbor
  do while (associated(nbor))
    if (nbor%task%is_set (bits%virtual)) then
      call check_recv (self, nbor%task)
    end if
    nbor => nbor%next
  end do
  call check_new (self)
  call link%lock%unset ('recv_from_vnbors')
  call trace%end ()
END SUBROUTINE recv_from_vnbors

!===============================================================================
!> Check a virtual task for completed messages.  Note that the procedure that
!> creates the task must make sure to also create a task message instance.
!===============================================================================
SUBROUTINE check_recv (self, task)
  class(task_mesg_t):: self
  class(task_t), pointer:: task
  !.............................................................................
  integer:: seq, tag
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%check_recv')
  !-----------------------------------------------------------------------------
  ! If the task has no pending request, issue one
  !-----------------------------------------------------------------------------
  call task%mesg%lock%set
  if (task%mesg%req == 0) then
    call new_request (task)
    !$omp atomic update
    mpi_mesg%nq_recv = mpi_mesg%nq_recv+1
  end if
  !-----------------------------------------------------------------------------
  ! If a message is complete, unpack and issue a new request
  !-----------------------------------------------------------------------------
  if (task%mesg%is_complete()) then
!if (mpi%rank==1) &
write (io_unit%log,*) '  received', task%id, task%mesg%id, task%mesg%tag
flush (io_unit%log)
    select type (task)
    class is (patch_t)
      call unpack (self, task%mesg, task%link)
    end select
    call new_request (task)
  else
write (io_unit%log,*) '  nothing', task%id, task%mesg%id, task%mesg%tag
flush (io_unit%log)
  end if
  call task%mesg%lock%unset
write (io_unit%log,*) '  lock%unset'
flush (io_unit%log)
  call trace%end ()
write (io_unit%log,*) '  trace%end'
flush (io_unit%log)
END SUBROUTINE check_recv



!===============================================================================
!> Increment sequence value and issue a new request
!===============================================================================
SUBROUTINE new_request (task)
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  select type (task)
  class is (patch_t)
    task%seq = task%seq + 1
    task%mesg%tag = task%id*100 + mod(task%seq,100)
    task%mesg%sender = task%rank
!if (mpi%rank==1) &
write (io_unit%log,*) ' requested', task%id, task%mesg%id, task%mesg%tag 
    call task%mesg%request
  end select
END SUBROUTINE new_request

!===============================================================================
!> Check for new virtual tasks
!===============================================================================
SUBROUTINE check_new (self)
  class(task_mesg_t):: self
  !.............................................................................
  class(mesg_t), pointer:: mesg
  class(experiment_t), pointer:: task => null()
  !$omp threadprivate(task)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Each thread has a thread-private new task allocated, ready to receive the
  ! next new virtual task, distinguished by having tag=MPI_ANY_TAG
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%check_new', itimer=itimer)
  if (.not.associated(task)) then
    call create_new (task)
    task%mesg%tag = 0
    task%id = 0
    call task%mesg%request (any=.true.)
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,1x,a,3i6)') wallclock(), &
        'task_mesg_t%check_new: requested id, sender, tag =', &
        task%id, task%rank, task%mesg%tag
    end if
  end if
  if (task%mesg%is_complete()) then
    task%id = task%unpack_id (task%mesg)
    task%rank = task%mesg%sender
    if (verbose > 1) then
      write (io_unit%log,'(f12.6,1x,a,2i6)') wallclock(), &
        'task_mesg_t%check_new: message id, rank =', &
        task%id, task%rank
    end if
    call unpack (self, task%mesg, task%link)
    call omp_lock%check_pointer (task%lock, 'task_mesg_t%check_new', task%id)
    call self%add_new_link (task%link)
    nullify(task)
    call create_new (task)
    task%mesg%tag = 0
    call task%mesg%request (any=.true.)
  else
    if (verbose > 1) then
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
    allocate (task_link)
    task%link => task_link
    task_link%task => task
    call task%init
    call omp_lock%check_pointer (task%lock, &
      'task_mesg_t%create_new task%lock', task%id)
    call task_link%init
    call task%allocate_mesg
write (io_unit%log,*) mpi%rank,' created new task', task%id
    task%rank = -1
    task%seq = 0
    call task%set (bits%virtual)
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
!===============================================================================
SUBROUTINE unpack (self, mesg, link)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer:: link
  !.............................................................................
  class(task_t), pointer:: task
  logical:: found, bonafide, new
  integer:: id
  integer, save:: itimer=0
  character(len=24):: label
  type(patch_t):: patch
  real(8):: wc
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%unpack', itimer=itimer)
  task => link%task
  select type (task)
  class is (experiment_t)
  !-----------------------------------------------------------------------------
  ! This entire operation should be threadsafe, since no other thread should be
  ! working on the same message and the same patch. And find_task uses only the
  ! static links in the task list -- not the dynamic next_time links.
  !-----------------------------------------------------------------------------
  !!omp critical (unpack_cr)
  new = .false.
  id = patch%unpack_id (mesg)
  !-----------------------------------------------------------------------------
  if (io%task_logging > 1) then
    write (label,'(a,i4,i8)') 'unpack ', mesg%sender, mesg%seq
    call task%log (label)
  end if
  !-----------------------------------------------------------------------------
  bonafide = .true.
  if (associated(task)) then
   if (verbose > 1) &
    write (io_unit%log,'(a,i6,2x,6l1)') 'unpack: id, bits IBVRES =', mesg%id, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%ready), &
      task%is_set (bits%remove), &
      task%is_set (bits%swap_request)
  end if
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
  ! Since an already existing patch may, at any one time, be under investigation
  ! by check_nbors, it must be protected by a critical region (or an OMP
  ! lock) while it is being updated here
  !-----------------------------------------------------------------------------
  call task%unpack (mesg)
  !$omp atomic
  mpi_mesg%n_unpk = mpi_mesg%n_unpk+1
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,3x,a,2i6,f12.6,2x,7l1)') &
      wallclock(), 'task_mesg_t%unpack: id, task%id, time, IBVrRSN =', &
      id, task%id, task%time, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%ready), &
      task%is_set (bits%remove), &
      task%is_set (bits%swap_request), &
      new
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
  if (mpi_mesg%debug .or. id==io%id_debug) then
    !$omp critical (debug_cr)
    write (io_unit%log,'(f12.6,2x,a,2i9,z12,2x,5l1)') wallclock(), &
      'DBG unpk: id, sender, req =', id, task%rank, mesg%req, &
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
      call self%init_nbors (link)
      call task%clear (bits%init_nbors)
    end if
    if (task%is_set(bits%update_status)) then
      !-------------------------------------------------------------------------
      ! If the task has the update_status bit set, re-evaluate status bits
      !-------------------------------------------------------------------------
      if (self%verbose > 1) write(io_unit%log,*) task%id, 'update_nbors'
      call self%update_status (task%link)
    end if
    if (task%is_set(bits%remove)) then
      if (self%verbose > 0) write(io_unit%log,*) task%id, 'remove'
      !-------------------------------------------------------------------------
      ! 1) Suicide note, remove the task and delete the message
      !-------------------------------------------------------------------------
      call self%remove_and_reset (link)
    else if (task%is_set(bits%swap_request)) then
      if (self%verbose > 0) write(io_unit%log,*) task%id, 'swap'
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
      if (self%verbose > 0) write(io_unit%log,*) task%id, 'new'
      !-------------------------------------------------------------------------
      ! 3) the task is a new virtual task.  It has the init_nbors bit set, and
      !    since it is a virtual task that bit should be cleared the next time
      !    the task is updated, coming back not as new task but as a virtual
      !    task update.  But first, issue a request for the next message.
      !-------------------------------------------------------------------------
      call task%setup
      call self%add_new_link (link)
    else
      if (self%verbose > 1) write(io_unit%log,*) task%id, 'update'
      !-------------------------------------------------------------------------
      ! 4) task was updated, so check nbors for ready tasks, after issuing a
      !    request for the next message.
      !-------------------------------------------------------------------------
      call self%check_nbors (link)
    end if
  else
    if (self%verbose > 1) write(io_unit%log,*) task%id, 'crap'
    !---------------------------------------------------------------------------
    ! 5) not bonafide, so ignore
    !---------------------------------------------------------------------------
  end if
  end select
  call trace%end (itimer)
END SUBROUTINE unpack

END MODULE task_mesg_mod
