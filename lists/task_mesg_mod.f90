!===============================================================================
!> Message handling for task lists.  Some of the methods are only used by
!> dispatcher0_t, so should perhaps be moved over to that file.
!>
!> Depending on switches, incoming MPU messages are handled by one of four
!> procedures (the values of the recv_active and recv_priv switches are shown
!> after the colon).  Only the ones with recv_active==.false. can (FIXME) be
!> used in AMR runs.
!>
!> check_virtual:TT specifically checks (only) for msgs to known virtual tasks
!> check_active :TF ditto, but not split up into thread-private lists
!> check_priv   :FT puts messages in recv_list and/or unpk_list
!> check_mpi    :FF original method, also uses progressive list handling
!>
!> When doing active receives (asking actively for the package to each virtual
!> task, it is necessary that the MPI message tag equals (or is a function of)
!> the task ID, while another, unique ID should be reserved for load balance
!> messages. To allow detecting packages for new patches (AMR), and still use
!> active receives, one could reserve another unique ID for "new" patches, for
!> which then the actual ID would be extracted from the header part of the
!> buffer.  Upon receipt of such a package, the new virtual patch would be
!> created, the first direct request would be sent, and the message would be
!> added to the (thread-private) list of messages.
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
  implicit none
  private
  type, public, extends(list_t):: task_mesg_t
    integer:: method=0
  contains
    procedure:: init
    procedure:: check_mpi
    procedure, private:: unpack
    procedure:: init_virtual
  end type
  integer, save:: verbose=0
  type(task_mesg_t), save:: virtual_list, remove_list
  !$omp threadprivate (virtual_list)
  public:: task2patch
  type(task_mesg_t), public:: task_mesg
CONTAINS

!===============================================================================
!> Check for and unpack MPI messages.
!===============================================================================
SUBROUTINE init (self, name)
  class(task_mesg_t):: self
  character(len=*), optional:: name
  integer:: iostat
  type(link_t):: link
  namelist /task_mesg_params/ verbose
  character(len=120):: id = &
    '$Id$ lists/task_mesg_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%print_id (id)
  rewind (io%input)
  read (io%input, task_mesg_params, iostat=iostat)
  if (io%master) write (io%output, task_mesg_params)
  call link%init_verbose (verbose)
END SUBROUTINE init

!===============================================================================
!> Check for and unpack MPI messages.
!===============================================================================
SUBROUTINE check_mpi (self)
  class(task_mesg_t):: self
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
  ! Step 1: Prune the send list -- critical regions are handled inside.
  ! Thread 0 checks the send list if it the only thread, or if it has messages
  ! on a thread private send list  -- otherwise it goes directly to handling 
  ! incoming messages.
  !-----------------------------------------------------------------------------
  !$omp atomic read
  nq = self%nq
  if (omp%nthreads==1 .or. mpi_mesg%send_priv) then
    call mpi_mesg%sent_list%check_sent (nq)
  else if (omp%thread > 0) then
    call mpi_mesg%sent_list%check_sent (nq)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Choose method used to receive messages:
  ! - 'check_virtual' uses prepared, thread private lists of virtual patches, 
  !   and thus does not work with AMR
  ! - 'check_active' runs through all virtual patches and issues irecv requests,
  !   re-issuing requests when they are completed
  ! - 'check_priv' used mpi_mesg%get to handle all incoming MPI packages, and
  !   and is thus the most reliable method, suffering only from the randomness
  !   that appears to inflict some MPI installations, where some messages may
  !   get tucked away for a significant time before being received
  ! - the default 'check_mpi' procedure is similar to 'check_priv', but uses 
  !   lists that are common to all threads, which requires using critical regions
  !-----------------------------------------------------------------------------
  if (mpi_mesg%recv_active) then
    if (mpi_mesg%recv_priv) then
      call check_virtual (self)
    else
      call check_active (self)
    end if
    return
  end if
  if (mpi_mesg%recv_priv) then
    call check_priv (self)
    return
  end if
  call trace%begin ('task_mesg_t%check_mpi', 1, itimer=itimer)
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
  !-----------------------------------------------------------------------------
  ! Step 3: if there are more than max_recv messages in recv_list, then wait
  ! for and unpack enough messages to reach the limit.  If task_list%nq is less
  ! then min_nq, reduce the max_recv to zero, forcing a wait on all receives.
  ! By splitting the operation into one part that moves the messages to a
  ! thead-local list, and a second part that does the wait and unpack, we
  ! minimize the blocking of other threads.
  !-----------------------------------------------------------------------------
  if (self%nq >= mpi_mesg%min_nq) then
    max_recv = 0
  else
    max_recv = mpi_mesg%max_recv
  end if
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
  !-----------------------------------------------------------------------------
  ! The indications of time reversal discoverde in Step 2 have now been cleared,
  ! but other ones could possibly exist, so we reset the test
  !-----------------------------------------------------------------------------
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
    call mpi_mesg%recv_list%remove (mesg)
    call unpk_tmp%add (mesg)
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
  return
END SUBROUTINE check_mpi

!===============================================================================
!> Find the link with a specified task id.  If there isn't any, create one!
!> If the patch already existed, it's nbor list should be OK, but the status
!> bit of it, and all its nbors, need to updated.
!===============================================================================
FUNCTION find_task (self, id, new) RESULT (link)
  class(task_mesg_t):: self
  integer:: id
  logical:: new
  logical:: error
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  class(experiment_t), pointer:: exper
  !.............................................................................
  call trace%begin ('task_mesg_t%find_task')
!  call self%lock%set ('find_task')
  new = .false.
  link => self%head
  do while (associated(link))
    patch => task2patch(link%task)
    if (patch%id == id) then
      if (verbose>1) then
        write (io_unit%log,*) mpi%rank,' task_mesg_t%find_task: match', id, associated(patch%mem)
        flush(io_unit%log)
      end if
!      call self%lock%unset ('find_task')
      call trace%end()
      return
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! No task found, which means it's a new virtual task that needs to be created
  !-----------------------------------------------------------------------------
  new = .true.
  allocate (link)
  allocate (exper)
  link%task => exper
  exper%link => link
  exper%id = id
  exper%box = self%size
  !-----------------------------------------------------------------------------
  ! The virtual task needs to have all of the capabilities and procedures of the
  ! original boundary task, since it may need to play a role in extras features.
  !-----------------------------------------------------------------------------
  call exper%init
  if (mpi_mesg%recv_priv .or. mpi_mesg%recv_active) then
    call exper%allocate_mesg
    if (mpi_mesg%recv_priv) then
      !-------------------------------------------------------------------------
      ! If communication is handled by a thread-private virtual task list, then
      ! append the new task to the (random) thread that happened to find it.
      !-------------------------------------------------------------------------
      task => exper
      call virtual_list%append (task, nbor=exper%link%nbor, &
        needed=exper%link%needed, needs_me=exper%link%needs_me)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Append to the task_list.  Since we don't yet know the position of the task,
  ! we can't update the nbor list yet; the calling procedure must do that, and
  ! we notify to do it via 'new', so it's not done every time a virtual task
  ! that already exists is updated.
  !-----------------------------------------------------------------------------
  call exper%set (bits%virtual)
  call self%append_link (link)
  if (verbose>0) &
    write (io_unit%log,*) ' task_mesg_t%find_task: append', id, exper%time
!  call self%lock%unset ('find_task')
  call trace%end()
END FUNCTION find_task

!===============================================================================
!> Unpack a message, where the MPI tag is the task id.  Use that to search
!> for the task, apply its unpack method, and check if any nbors become ready.
!===============================================================================
SUBROUTINE unpack (self, mesg, link)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer, optional:: link
  class(link_t), pointer:: link1, link2
  class(experiment_t), pointer:: task
  logical:: found, failed, new
  integer:: id, n_added
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%unpack', itimer=itimer)
  if (mesg%nbuf < 40) then
    call load_balance%unpack (mesg%buffer)
    return
  end if
  call self%lock%set ('unpack')
  !-----------------------------------------------------------------------------
  ! This entire operation should be threadsafe, since no other thread should be
  ! working on the same message and the same patch. And find_task uses only the
  ! static links in the task list -- not the dynamic next_time links.
  !-----------------------------------------------------------------------------
  !!omp critical (unpack_cr)
  new = .false.
  mesg%id = mesg%id
  if (present(link)) then
    link1 => link
    if (self%verbose > 1) &
      write (io_unit%log,'(f12.6,2x,a,2i5,z12)') &
        wallclock(), 'task_mesg_t%unpack: id =', mesg%id, link1%task%id
  else
    link1 => find_task (self, mesg%id, new)
    if (self%verbose > 1) &
      write (io_unit%log,'(f12.6,2x,a,2i5,z12)') &
        wallclock(), 'task_mesg_t%unpack: found id =', mesg%id, link1%task%id
  end if
  !-----------------------------------------------------------------------------
  ! Make sure the task pointer has all experiment_t attributes
  !-----------------------------------------------------------------------------
  associate (ltask=>link1%task)
  select type (ltask)
  class is (experiment_t)
  task => ltask
  end select
  end associate
  failed = .false.
  !-----------------------------------------------------------------------------
  if (associated(link1%task)) then
   if (verbose > 1) &
    write (io_unit%log,'(a,i6,2x,5l1)') 'unpack: id, bits BVRES =', mesg%id, &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%ready), &
      task%is_set (bits%remove), &
      task%is_set (bits%swap_request)
  end if
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
  ! Measure the update cadence for this virtual task
  !---------------------------------------------------------------------------
  block
  real(8):: wc
  wc = wallclock()
  task%update_cadence = wc - task%update_last
  task%update_last = wc
  end block
  !---------------------------------------------------------------------------
  ! Unpack a patch message (which includes swapping the roles of boundary bits).
  ! Since an already existing patch may, at any one time, be under investigation
  ! by check_nbors, it must be protected by a critical region (or an OMP
  ! lock) while it is being updated here
  !---------------------------------------------------------------------------
  id = task%id
  call task%lock%set ('unpack')
  call task%unpack (mesg)
  if (verbose > 1) then
    write (io_unit%mpi,'(a,3i6,f12.6,2x,3f10.4,2x,5l1)') &
      'unpack: after task%unpack BVRES =', &
      id, mesg%id, task%id, task%time, task%position, &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%ready), &
      task%is_set (bits%remove), &
      task%is_set (bits%swap_request)
    flush (io_unit%mpi)
  end if
  !-----------------------------------------------------------------------------
  ! Check for suicide note
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%remove)) then
    if (verbose > 0) then
      write (io_unit%mpi,'(f12.6,2x,a,i6)') &
        wallclock(), 'unpack: suicide note received for id =', task%id
      flush (io_unit%mpi)
    end if
    !---------------------------------------------------------------------------
    ! If the task is on a thread-private virtual task list things are more
    ! complicated; we need to keep trying until this procedure is executed
    ! by the correct thread
    !---------------------------------------------------------------------------
    if (mpi_mesg%recv_priv) then
      call remove_list%append_link (link1)
    !---------------------------------------------------------------------------
    ! If not we need a normal remove, including new nbor lists etc.
    !---------------------------------------------------------------------------
    else
      call self%remove_and_reset (link1)
    end if
    call self%lock%unset('unpack')
    call trace%end (itimer)
    return
  end if
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
      !$omp atomic update
      self%na = self%na+1
      !$omp atomic update
      self%nb = self%nb+1
      !$omp atomic update
      self%nv = self%nv-1
      call self%init_nbors (link1)
      call task%clear (bits%swap_request+bits%ready)
      call self%update_nbor_status (link1)
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
        write (io_unit%log,*) 'task_mesg_t%unpack: new virtual patch:', task%id
      self%nv = self%nv+1
      call self%init_nbors (link1)
      call task%clear (bits%swap_request+bits%ready)
      call self%update_nbor_status (link1)
      call self%count_status
    end if
  end if
  call task%lock%unset ('unpack')
  !-----------------------------------------------------------------------------
  ! If a virtual task was just created it, it and its nbors need new nbor lists.
  ! Also, since the proper size and position were not known when a new task was
  ! created, but were set by the task%unpack, we need to redo the task%setup,
  ! which regenerates the meshes and sets up all geometric properties.  New AMR
  ! tasks may be created both directly, via calls from selective_refine(), or
  ! indirectly, via calls in check_support(); at this point this makes little
  ! or no difference.  Check_support is associated one-to-one with new AMR
  ! patches, so it should be sufficient to include the call under the if (new).
  !-----------------------------------------------------------------------------
  if (new) then
    call task%setup
    if (verbose > 0) then
      write (io_unit%mpi,'(f12.6,2x,a,i6,2i4,f12.6,2x,3f10.5,2x,"new")') &
        wallclock(), 'task_mesg_t%unpack: ', id, mesg%seq, omp%thread, &
        task%time, task%position
      flush(io_unit%mpi)
    end if
    !call self%lock%set ('unpack')
    call self%init_nbors (link1)
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    call refine%check_support (self, link1, n_added)
    call self%check_nbors (link1)
    !call self%lock%unset ('unpack')
  else
    if (verbose > 0) then
      write (io_unit%mpi,'(f12.6,2x,a,i6,2i4,f12.6,2x,3f10.5,2x)') &
        wallclock(), 'task_mesg_t%unpack: ', id, mesg%seq, omp%thread, &
        task%time, task%position
      flush(io_unit%mpi)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Patches whose nbor lists are reset by init_nbor_nbors have the init_nbors
  ! bit set, in order to replicate the behavior on the nbor ranks
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%init_nbors)) then
    call self%init_nbors (link1)
    if (verbose > 1) &
      write (io_unit%log,*) 'unpack: bits%init_nbors id =', task%id
  end if
  !!omp end critical (unpack_cr)
  !-----------------------------------------------------------------------------
  ! As the task has now been updated, we need to check if any of the neighbors
  ! have become ready to update. If the task has just been swapped into being
  ! a boundary patch, it could possibly also be ready to update.
  !-----------------------------------------------------------------------------
  if (self%method==0) then
    if (verbose > 1) &
      write (io_unit%log,*) 'unpack: check_nbors, id =', task%id
    call self%check_nbors (link1)
  end if
  call self%lock%unset ('unpack')
  call trace%end (itimer)
END SUBROUTINE unpack

!===============================================================================
!> Check for and unpack MPI messages, using threadprivate recv_list and
!> unpk_list
!===============================================================================
SUBROUTINE check_priv (self)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg, next
  integer, save:: itimer=0
  integer:: n
  !.............................................................................
  if (mpi%size <= 1) return
  call trace%begin ('task_mesg_t%check_priv', 1, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Step 1: Add incoming messages to either recv_list or unpk_list
  !-----------------------------------------------------------------------------
  call mpi_mesg%get (mesg)
  n = 0
  do while (associated(mesg))
    if (mesg%is_complete()) then
      if (verbose > 1) &
        write (io_unit%log,*) 'unpacked directly', mesg%id, n
      call self%unpack (mesg)
      call recv_list%delete (mesg, .false.)
    else
      call recv_list%add (mesg)
      n = n+1
      if (verbose > 1) &
        write (io_unit%log,*) 'addded mesg to recv_list', mesg%id, n
    end if
    call mpi_mesg%get (mesg)
  end do
  !-----------------------------------------------------------------------------
  ! Step 2: check the recv_list for completed messages
  !-----------------------------------------------------------------------------
  mesg => recv_list%head
  do while (associated(mesg))
    next => mesg%next
    if (mesg%is_complete()) then
      call self%unpack (mesg)
      call recv_list%remove (mesg)
      n = n-1
      if (verbose > 1) &
        write (io_unit%log,*) 'unpacked from list', mesg%id, n 
      call recv_list%delete (mesg, .false.)
    end if
    mesg => next
  end do
  call trace%end (itimer)
END SUBROUTINE check_priv

!===============================================================================
!> Run through all virtual patches, count them, and create separate threadprivate
!> lists that together include all of them.
!===============================================================================
SUBROUTINE init_virtual (self)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: nvirt, i, i1, i2
  integer, save:: itimer=0
  !.............................................................................
  if (mpi%size==1) return
  call trace%begin ('task_mesg_t%init_virtual', itimer=itimer)
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
  !$omp parallel shared(verbose,nvirt,mpi,bits,self,io) private(link,i,i1,i2,task) default(none)
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
        if (verbose > 0) &
          write (io_unit%output,*) 'init_virtual irecv request:', task%id, task%mesg%id, task%mesg%tag
        !if (mpi%rank==0) &
        !  write(io%output,*) mpi%rank,'MK',omp_get_thread_num(),i,task%id,'irec'
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
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: nunpk
  real(8):: wc
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('task_mesg_t%check_virtual', itimer=itimer)
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
    call set_test_time (self, task)
    if (task%mesg%is_complete()) then
      !write(stdout,*) 'recv:', task%id, task%mesg%id, task%mesg%tag, task%rank
      if (verbose > 1) &
        write (io_unit%log,*) wallclock(),' task_mesg_t%unpack: check_virtual hit ', task%id
      call self%unpack (task%mesg, link=link)
      !$omp atomic
      mpi_mesg%n_recv = mpi_mesg%n_recv+1
      !$omp atomic
      timer%bytes_recv = timer%bytes_recv + 4.0_8*task%mesg%nbuf
      !$omp atomic
      timer%n_recv = timer%n_recv + 1_8
      wc = wallclock()
      if (wc-task%wc_last > 10.) write(io%output, &
        '("WARNING: virtual patch not updated in",f5.1," sec")') wc
      task%wc_last = wc
      call task%mesg%irecv (task%rank, task%id)
      !write(stdout,*) 'ircv:', task%id, task%mesg%id, task%mesg%tag, task%rank
      nunpk = nunpk+1
    else if (verbose > 2) then
      write (io_unit%log,*) wallclock(),' task_mesg_t%unpack: check_virtual miss', task%id
    end if
    end select
    link => link%next
  end do
  if (verbose>1) write (io_unit%log,*) ' task_mesg_t%unpack: check_virtual', nunpk
  call trace%end (itimer)
END SUBROUTINE check_virtual

!===============================================================================
!> Set the mesg%test_time
!===============================================================================
SUBROUTINE set_test_time (self, task)
  class(task_mesg_t):: self
  class(patch_t):: task
  !---------------------------------------------------------------------------
  ! If TEST_TIME is specified in MPI_MESG_PARAMS, set that as the minimum time
  ! between MPI_TEST for the same task
  !---------------------------------------------------------------------------
  if (mpi_mesg%test_time > 0.0) then
    task%mesg%test_time = mpi_mesg%test_time
  !---------------------------------------------------------------------------
  ! If the parameter is zero or negative, try to estimate an optimal time.
  ! The first test is allowed after half the estimated update_cadence, and
  ! after that the time interval is halved after every fail.  This point is
  ! returned to about as often as a task update occurs on a single thread,
  ! which is about the task overload times shorter than the update_cadence.
  !---------------------------------------------------------------------------
  else
    task%update_cadence = 1e-6*product(task%n)*self%na/omp%nthreads
    task%mesg%test_time = task%update_cadence*0.5**(task%mesg%n_failed+1)
  end if
END SUBROUTINE set_test_time

!===============================================================================
!> Run through all virtual patches, and allocate message data types if not
!> allocated.  If an old request exists, test if it has completed.  If yes,
!> unpack and issue a new request.  If no request exists, issue a request.
!>
!> This procedure, which is only executed by the master thread, and is only
!> used when more than one thread is available, asks for exactly the next
!> message tag, and does not issue a new irecv until that message has arrived,
!> hence it needs only one (permanent) mesg instance per task.
!===============================================================================
SUBROUTINE check_active (self)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: req
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('task_mesg_t%check_active', itimer=itimer)
  call self%lock%set ('check_active')
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%virtual)) then
      select type (task)
      class is (patch_t)
      !-------------------------------------------------------------------------
      ! Check if this virtual task has a task message allocated.  If not, it is
      ! a new virtual task -- for example created by AMR -- and should have a
      ! message allocated, with a active receive request
      !-------------------------------------------------------------------------
      if (.not.associated(task%mesg)) then
        call task%allocate_mesg
      end if
      if (task%mesg%req == 0) then
        call irecv (task)
        task%wc_last = wallclock()
        !write(stdout,*) 'ircv:', task%id, task%mesg%id, task%mesg%tag, task%rank
      end if
      !-------------------------------------------------------------------------
      ! If a mesg was recvd, then unpack and issue a new request,
      !-------------------------------------------------------------------------
      call set_test_time (self, task)
      req = task%mesg%req
      if (task%mesg%is_complete()) then
        write(stdout,*) 'iscmp:', task%id, task%mesg%id, task%mesg%tag, task%rank
        if (verbose > 0 .or. task%id == io%id_debug) &
          write (io_unit%log,'(f12.6,2x,a,i6,i5,i9,i6,z12)') wallclock(), &
            'task_mesg%tcheck_active:  recv id, seq, mesg%id, sender, req =', &
            task%id, mod(task%istep,100), task%mesg%id, task%mesg%sender, req
        !$omp atomic
        mpi_mesg%n_recv = mpi_mesg%n_recv+1
        !$omp atomic
        timer%bytes_recv = timer%bytes_recv + 4.0_8*task%mesg%nbuf
        !$omp atomic
        timer%n_recv = timer%n_recv + 1_8
        call self%unpack (task%mesg)
        call irecv (task)
        !write(stdout,*) 'irev:', task%id, task%mesg%id, task%mesg%tag, task%rank
      else
        if (verbose > 1 .or. task%id == io%id_debug) &
          write (io_unit%log,'(f12.6,2x,a,i6,i5,i9,i6,z12)') wallclock(), &
            'task_mesg%tcheck_active:  fail id, seq, mesg%id, sender, req =', &
            task%id, mod(task%istep,100), task%mesg%id, task%rank, task%mesg%req
      end if
      end select
    end if
    link => link%next
  end do
  call self%lock%unset ('check_active')
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  ! If mpi_mesg%uniq_mesg is set, the message tag is constructed intenally,
  ! and if not, it is equal to the task%id
  !-----------------------------------------------------------------------------
  subroutine irecv (task)
  class(patch_t):: task
  call task%mesg%irecv (task%rank, task%id)
  if (verbose > 0 .or. task%id == io%id_debug) &
    write (io_unit%log,'(f12.6,2x,a,i6,i5,i9,i6,z12)') wallclock(), &
      'task_mesg%tcheck_active: irecv id, seq, mesg%id, sender, req =', &
      task%id, mod(task%istep,100), task%mesg%id, task%rank, task%mesg%req
  end subroutine
END SUBROUTINE check_active

END MODULE task_mesg_mod
