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
!> task template for new virtual task, with an activ request waiting for
!> messags with tag==0.  When such a message is received, the message is
!> unpacked into the virtual task template and the thread then adds it to 
!> the task list, using the same procedure as was used on the sender side 
!> (list_t%add_new_link).
!>
!> Each virtual task starts updating from the owner rank by issuing a request
!> for a message with the predicted tag, which thus goes directly to the correct task
!> task.  Subsequent calls checks if the message has been received, unpacks it
!> if it has, and then issues a request for the next package.
!===============================================================================
MODULE task_mesg_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE omp_mod
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
  USE load_balance_mod
  USE refine_mod
  USE hash_table_mod
  implicit none
  private
  type, public, extends(list_t):: task_mesg_t
    integer:: method=0
    type(hash_table_t):: hash_table
  contains
    procedure:: init
    procedure:: check_mpi
    procedure, private:: unpack
    procedure:: send_to_vnbors
    procedure:: recv_from_vnbors
    procedure:: check_vnbors
  end type
  integer, save:: verbose=0
  logical, save:: use_hashtable=.true.
  type(mesg_list_t), save:: priv_list
  type(task_mesg_t), save:: private_list, master_list, remove_list
  !$omp threadprivate (private_list,priv_list)
  character(len=8), save:: method='irecv'
  public:: task2patch
  type(task_mesg_t), public:: task_mesg
CONTAINS

!===============================================================================
!> Check for and unpack MPI messages.
!===============================================================================
SUBROUTINE init (self, name)
  class(task_mesg_t), target:: self
  character(len=*), optional:: name
  integer:: iostat
  integer, save:: hash_table_size=10000
  type(link_t):: link
  namelist /task_mesg_params/ verbose, method, hash_table_size, use_hashtable
  character(len=120):: id = &
    '$Id$ lists/task_mesg_mod.f90'
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) &
    return
  call trace%begin ('task_mesg_t%init')
  call trace%print_id (id)
  rewind (io%input)
  read (io%input, task_mesg_params, iostat=iostat)
  call link%init_verbose (verbose)
  write (io%output, task_mesg_params)
  if (use_hashtable) &
    call self%hash_table%init (hash_table_size)
  !-----------------------------------------------------------------------------
  ! Choose a unique tag for each message and task ID
  !-----------------------------------------------------------------------------
  mpi_mesg%uniq_mesg = .true.
  mpi_mesg%tag_type = 1
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Stub
!===============================================================================
SUBROUTINE check_mpi (self, n)
  class(task_mesg_t):: self
  integer, optional:: n
END SUBROUTINE check_mpi

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
SUBROUTINE unpack (self, mesg, link1)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer:: link1, link2
  class(experiment_t), pointer:: task
  logical:: found, bonafide, new
  integer:: id
  integer, save:: itimer=0
  character(len=24):: label
  type(patch_t):: patch
  real(8):: wc
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%unpack', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! This entire operation should be threadsafe, since no other thread should be
  ! working on the same message and the same patch. And find_task uses only the
  ! static links in the task list -- not the dynamic next_time links.
  !-----------------------------------------------------------------------------
  !!omp critical (unpack_cr)
  new = .false.
  id = patch%unpack_id (mesg)
write (stdout,*) 'task_mesg_t%unpack: task id, seq =', id, mesg%seq
  !-----------------------------------------------------------------------------
  ! Make sure the task pointer has all experiment_t attributes
  !-----------------------------------------------------------------------------
  associate (ltask=>link1%task)
  select type (ltask)
  class is (experiment_t)
  task => ltask
  end select
  end associate
  !-----------------------------------------------------------------------------
  if (io%task_logging > 1) then
    write (label,'(a,i4,i8)') 'unpack ', mesg%sender, mesg%seq
    call task%log (label)
  end if
  !-----------------------------------------------------------------------------
  bonafide = .true.
  if (associated(link1%task)) then
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
    write (io_unit%log,'(f12.6,3x,a,2i6,f12.6,2x,3f10.4,2x,7l1)') &
      wallclock(), 'unpack: after task%unpack IBVrRSN =', &
      id, task%id, task%time, task%position, &
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
    write (io_unit%sent,'(f16.6,i4,2x,a,i6,f16.6,8i5)') wallclock(), omp%thread, 'unp', task%id, task%time, task%rank
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
    write (io_unit%log,'(f12.6,2x,a,i9,1p,e18.6)') wallclock(), 'unpk: id, time =', task%id, task%time
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
      call self%init_nbors (link1)
      call link1%task%clear (bits%init_nbors)
    end if
    if (task%is_set(bits%update_status)) then
      !-------------------------------------------------------------------------
      ! If the task has the update_status bit set, re-evaluate status bits
      !-------------------------------------------------------------------------
      if (self%verbose > 1) write(io_unit%log,*) task%id, 'update_nbors'
      call self%update_status (link1)
    end if
    if (task%is_set(bits%remove)) then
      if (self%verbose > 0) write(io_unit%log,*) task%id, 'remove'
      !-------------------------------------------------------------------------
      ! 1) Suicide note, remove the task and delete the message
      !-------------------------------------------------------------------------
      call self%remove_and_reset (link1)
      call mpi_mesg%unpk_list%delete (mesg, .false.)    !XX
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
      call self%swap (link1, .false.)
      if (load_balance%is_on) &
        call load_balance%lock%unset ('task_mesg_t%unpack')
      call self%lock%unset ('task_mesg_t%unpack')
      call mpi_mesg%unpk_list%delete (mesg, .false.)    !XX
    else if (new) then
      if (self%verbose > 0) write(io_unit%log,*) link1%task%id, 'new'
      !-------------------------------------------------------------------------
      ! 3) the task is a new virtual task.  It has the init_nbors bit set, and
      !    since it is a virtual task that bit should be cleared the next time
      !    the task is updated, coming back not as new task but as a virtual
      !    task update.  But first, issue a request for the next message.
      !-------------------------------------------------------------------------
      call mpi_mesg%unpk_list%delete (mesg, .false.)    !XX
      !XX call mpi_mesg%iget (mesg)
      !XX call mpi_mesg%recv_list%add (mesg)
      call task%setup
      call self%add_new_link (link1)
    else
      if (self%verbose > 1) write(io_unit%log,*) link1%task%id, 'update'
      !-------------------------------------------------------------------------
      ! Get a new message from the same source; this is the message with the
      ! next sequnce number from that source -- typically a message for a
      ! a different task ID.  To avoid that the message in recv_list starts to
      ! be used before it obtained, the %add() must follow the %iget()
      !-------------------------------------------------------------------------
      !call mpi_mesg%unpk_list%delete (mesg, .false.)    !XX
write (stdout,*) 'received seq:', mesg%seq, '  recv_list%n =', mpi_mesg%recv_list%n
      call mpi_mesg%iget (mesg)                         !XX
      call mpi_mesg%recv_list%add (mesg)                !XX
write (stdout,*) 'request  seq:', mesg%seq, '  recv_list%n =', mpi_mesg%recv_list%n
      !-------------------------------------------------------------------------
      ! 4) task was updated, so check nbors for ready tasks, after issuing a
      !    request for the next message.
      !-------------------------------------------------------------------------
      call self%check_nbors (link1)
    end if
  else
    if (self%verbose > 1) write(io_unit%log,*) link1%task%id, 'crap'
    !---------------------------------------------------------------------------
    ! 5) not bonafide, so ignore
    !---------------------------------------------------------------------------
    continue
  end if
  call trace%end (itimer)
END SUBROUTINE unpack

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
  !-----------------------------------------------------------------------------
  ! Sending the task creates a task%mesg instance, which is appended to
  ! the task%sent_list -- created on the fly if not existing
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%send_to_vnbors')
  call self%list_t%send_to_vnbors (link)
  if (.not.associated(link%task%sent_list)) then
    allocate (link%task%sent_list)
  end if
  call link%task%sent_list%add (link%task%mesg)
  call self%check_vnbors (link)
  call trace%end ()
END SUBROUTINE send_to_vnbors

!===============================================================================
!> Check recv messages on all virtual task nbors
!===============================================================================
SUBROUTINE check_vnbors (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(mesg_t), pointer:: mesg, next
  !-----------------------------------------------------------------------------
  ! The function mesg%test_sent() returns true when all ranks have received
  ! the message, which can then be removed and deleted, thus reclaiming memory.
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%check_vnbors')
  if (associated(link%task%sent_list)) then
write(io_unit%log,*) 1,associated(link%task); flush (io_unit%log)
write(io_unit%log,*) 2,associated(link%task%sent_list); flush (io_unit%log)
write(io_unit%log,*) 3,associated(link%task%sent_list%head); flush (io_unit%log)
    mesg => link%task%sent_list%head
    do while (associated(mesg))
      next => mesg%next
      if (mesg%test_sent()) then
        call link%task%sent_list%remove (mesg)
        call link%task%sent_list%delete (mesg, send=.true.)
      end if
      mesg => next
    end do
  end if
  call trace%end()
END SUBROUTINE check_vnbors

!===============================================================================
!> Check the virtual nbors of a boundary task link for completed messages,
!> unpack the ones that are complete, and issue new requests for those task.
!===============================================================================
SUBROUTINE recv_from_vnbors (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor
  class(mesg_t), pointer:: mesg
  class(experiment_t), pointer:: virtual => null()
  !$omp threadprivate(virtual)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%recv_from_vnbors', itimer=itimer)
  nbor => link%nbor
  do while (associated(nbor))
    write(io_unit%log,*) 1,associated(nbor%task); flush (io_unit%log)
    if (nbor%task%is_set (bits%virtual)) then
      !-------------------------------------------------------------------------
      ! The first time a task is checked there is no active request, try to 
      ! allocate one, punting if the buffer size (computed in patch%pack()) is 
      ! not yet available.
      !-------------------------------------------------------------------------
      if (.not.associated(nbor%task%mesg)) then
        allocate (nbor%task%mesg)
        mesg => nbor%task%mesg
        if (mpi_mesg%nbuf == 0) then
          write (stderr,*) mpi%rank, 'WARNING: mpi_mesg%nbuf is still zero'
          call trace%end (itimer)
          return
        end if
        mesg%nbuf = mpi_mesg%nbuf
        allocate (mesg%buffer(mesg%nbuf))
        mesg%sender = nbor%task%rank
        nbor%task%seq = 1
      end if
      write(io_unit%log,*) 2,associated(nbor%task%mesg); flush (io_unit%log)
      if (nbor%task%mesg%req == 0) then
        call nbor%task%mesg%request (nbor%task%id*100 + 1)
      end if
      if (nbor%task%mesg%is_complete()) then
        !-----------------------------------------------------------------------
        ! When a request is complete, unpack the task and issue a new request
        !-----------------------------------------------------------------------
        call self%unpack (nbor%task%mesg, nbor%link)
        call nbor%task%mesg%request (nbor%task%istep+1)
      end if
    end if
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Each thread has a thread-private new task allocated, ready to receive the
  ! next new virtual task, distinguished by having tag=0
  !-----------------------------------------------------------------------------
  if (.not.associated(virtual)) then
    call create_new_task
    call virtual%mesg%request (0)
  end if
  if (virtual%mesg%is_complete()) then
    virtual%id = virtual%unpack_id (virtual%mesg)
    call self%unpack (virtual%mesg, virtual%link)
    call self%add_new_link (virtual%link)
    nullify(virtual)
  end if
  call trace%end (itimer)
contains
  !=============================================================================
  ! Internal subroutine to create a new virtual task
  !=============================================================================
  subroutine create_new_task
    class(link_t), pointer:: task_link
    !---------------------------------------------------------------------------
    call trace%begin ('task_mesg_t%create_new_task')
    allocate (virtual)
    allocate (task_link)
    virtual%link => task_link
    task_link%task => virtual
    call virtual%init
    call task_link%init
    call virtual%allocate_mesg
    call virtual%set (bits%virtual)
    call trace%end ()
  end subroutine create_new_task
END SUBROUTINE recv_from_vnbors

END MODULE task_mesg_mod
