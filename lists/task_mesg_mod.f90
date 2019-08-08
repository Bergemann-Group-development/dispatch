!===============================================================================
!> Message handling for task lists.  Some of the methods are only used by
!> dispatcher0_t, so should perhaps be moved over to that file.
!>
!> Depending on switches, incoming MPU messages are handled by one of four
!> procedures (the values of the recv_active and recv_priv switches are shown
!> after the colon).  Only the ones with recv_active==.false. can (FIXME) be
!> used in AMR runs.
!>
!> recv_master  : master puts MPI_IMPROBE messages in recv_list and unpk_list
!> recv_private : master checks (only) for msgs to known virtual tasks
!> recv_private : ditto, but split up into thread-private lists
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
    procedure, private:: recv_master
    procedure:: init_virtual
    procedure, private:: recv_virtual
    procedure, private:: find_task
    procedure, private:: unpack
  end type
  integer, save:: verbose=0
  type(task_mesg_t), save:: private_list, master_list, remove_list
  !$omp threadprivate (private_list)
  character(len=8), save:: method='master'
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
  namelist /task_mesg_params/ verbose, method
  character(len=120):: id = &
    '$Id$ lists/task_mesg_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%init')
  call trace%print_id (id)
  rewind (io%input)
  read (io%input, task_mesg_params, iostat=iostat)
  call link%init_verbose (verbose)
  write (io%output, task_mesg_params)
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

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
!> Check for and unpack MPI messages.
!===============================================================================
SUBROUTINE check_mpi (self)
  class(task_mesg_t):: self
  integer:: nq
  !.............................................................................
  if (mpi%size <= 1) return
  call trace%begin ('task_mesg_t%check_mpi')
  !-----------------------------------------------------------------------------
  ! Check the list of sent messages
  !-----------------------------------------------------------------------------
  !$omp atomic read
  nq = self%nq
  call mpi_mesg%sent_list%check_sent (nq)
  !-----------------------------------------------------------------------------
  ! Choose MPI recv method -- note that we force the hand in the AMR case
  !-----------------------------------------------------------------------------
  if (refine%on) then
    method = 'master'
  end if
  !-----------------------------------------------------------------------------
  ! The 'private' method runs on all threads; the others only on master
  !-----------------------------------------------------------------------------
  select case (trim(method))
  case ('master')
    if (omp%master) &
      call self%recv_master
  case ('virtual')
    if (omp%master) &
      call self%recv_virtual (master_list)
  case ('private')
    call self%recv_virtual (private_list)
  end select
  call trace%end
END SUBROUTINE check_mpi

!===============================================================================
!> Check for and unpack MPI messages, using only the master thread, to avoid
!> the need for critical regions and locks.  The cadence of checking for
!> incoming messages is once per task update on one thread, which is suitable,
!> since each task will be "served" about once per task update time. 
!===============================================================================
SUBROUTINE recv_master (self)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg, next
  integer, save:: itimer=0
  integer:: n
  !.............................................................................
  if (omp%thread > 0) return
  call trace%begin ('task_mesg_t%recv_master', 1, itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Step 1: Add incoming messages to either recv_list or unpk_list
  !-----------------------------------------------------------------------------
  mpi_mesg%uniq_mesg = .true.
  mpi_mesg%tag_type = 2
  call mpi_mesg%get (mesg)
  n = 0
  do while (associated(mesg))
    call recv_list%add (mesg)
    n = n+1
    if (verbose > 1) &
      write (io_unit%log,*) 'recv_master: addded mesg to recv_list', &
        mesg%id, mesg%seq
    call mpi_mesg%get (mesg)
  end do
  !-----------------------------------------------------------------------------
  ! Step 2: check the recv_list for completed messages and move to unpk_list
  !-----------------------------------------------------------------------------
  n = 0
  mesg => recv_list%head
  do while (associated(mesg))
    next => mesg%next
    n = n+1
    if (mesg%is_complete('recv_master')) then
      if (verbose > 1) &
        write (io_unit%log,*) 'recv_master: moved to unpk_list:', &
          mesg%id, mesg%seq, associated(next)
      call recv_list%remove (mesg)
      call unpk_list%add (mesg)
    else
      if (verbose > 2) &
        write (io_unit%log,*) 'recv_master:   not yet complete:', &
          mesg%id, mesg%seq
    end if
    mesg => next
  end do
  !-----------------------------------------------------------------------------
  ! Step 3: unpack in sequence order
  !-----------------------------------------------------------------------------
  mesg => unpk_list%head
  do while (associated(mesg))
    next => mesg%next
    if (mesg%is_in_order()) then
      call self%unpack (mesg)
      n = n-1
      if (verbose > 1) &
        write (io_unit%log,*) 'recv_master: unpacked from list:', &
          mesg%id, mesg%seq
      call unpk_list%remove (mesg)
      call unpk_list%delete (mesg, .false.)
    else
      if (verbose > 2) &
        write (io_unit%log,*) 'recv_master:  kept in unpk_list:', &
          mesg%id, mesg%seq
    end if
    mesg => next
  end do
  call trace%end (itimer)
END SUBROUTINE recv_master

!===============================================================================
!> Run through all virtual patches, count them, and optionally create separate
!> threadprivate lists that together include all of them.
!===============================================================================
SUBROUTINE init_virtual (self)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: nvirt, per, i, i1, i2
  !.............................................................................
  if (mpi%size==1) return
  if (trim(method) == 'master') return
  call trace%begin ('task_mesg_t%init_virtual')
  mpi_mesg%uniq_mesg = .true.
  mpi_mesg%tag_type = 1
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
  per = nvirt/omp%nthreads + 1
  if (trim(method) == 'private') then
    !---------------------------------------------------------------------------
    ! Each thread copies its share of virtual patches to its private list, and
    ! initializes its private virtual task list
    !---------------------------------------------------------------------------
    !$omp parallel shared(verbose,nvirt,per,mpi,bits,self,io) private(i1,i2) default(none)
    i1 = 1  + per*omp%thread
    i2 = i1 + per
    call init_virtual_list (self, private_list, i1, i2)
    !$omp end parallel
  else if (trim(method) == 'virtual') then
    call init_virtual_list (self, master_list, 1, 1+nvirt)
  end if
  call trace%end ()
END SUBROUTINE init_virtual

!===============================================================================
!> Collect the virtual tasks in the enumerated i1..i2 interval onto reduced list
!===============================================================================
SUBROUTINE init_virtual_list (self, list, i1, i2)
  class(task_mesg_t):: self, list
  integer, optional:: i1, i2
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%init_virtual_list')
  i = 1
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%virtual)) then
      if (i >= i1 .and. i < i2) then
        ! -- append the task to the virtal list --
        call list%append (task, &
          nbor=link%nbor, needed=link%needed, needs_me=link%needs_me)
        ! -- allocate a task message --
        call task%allocate_mesg
        task%wc_last = wallclock()
        ! -- initialize the 1st receive --
        call task%mesg%irecv (task%rank, task%id)
        if (verbose > 0) &
          write (io_unit%log,*) &
            'task_mesg_t%init_virtual: id, mesg%id, mesg%tag =', &
            task%id, task%mesg%id, task%mesg%tag
      end if
      i = i+1
    end if
    link => link%next
  end do
  if (verbose > 0) &
    write (io_unit%log,'(a,2i5,2x,2i4)') &
      'task_mesg_t%init_virtual_list: n, thread, i1, i2 =', &
      list%n, omp%thread, i1, i2
  call trace%end
END SUBROUTINE init_virtual_list

!===============================================================================
!> Run through all virtual patches, unpack completed receives, and issue new
!> receive requests.  The virtual_list is threadprivate, so no critical region
!> is needed.
!===============================================================================
SUBROUTINE recv_virtual (self, list)
  class(task_mesg_t):: self, list
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: nunpk
  real(8):: wc
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('task_mesg_t%recv_virtual', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Each thread runs through its list of virtual tasks, unpacks completed
  ! messages, and issues a new request
  !-----------------------------------------------------------------------------
  nunpk = 0
  link => list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
    call set_test_time (self, task)
    if (task%mesg%is_complete('virtual')) then
      if (verbose > 1 .or. task%id == io%id_debug) &
        write (io_unit%log,'(f12.6,2x,a,i6,i5,i9,i6,z12)') wallclock(), &
          'task_mesg_t%recv_virtual:  recv id, seq, mesg%id, sender =', &
          task%id, mod(task%istep,100), task%mesg%id, task%mesg%sender
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
    else if (verbose > 2 .or. task%id == io%id_debug) then
      write (io_unit%log,'(f12.6,2x,a,i6,i5,i9,i6,z12)') wallclock(), &
        'task_mesg_t%recv_virtual:  fail id, seq, mesg%id, sender, req =', &
        task%id, mod(task%istep,100), task%mesg%id, task%rank, task%mesg%req
    end if
    end select
    link => link%next
  end do
  if (verbose > 0) &
    write (io_unit%log,*) ' task_mesg_t%recv_virtual, nunpk =', nunpk
  call trace%end (itimer)
END SUBROUTINE recv_virtual

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
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('task_mesg_t%find_task', itimer=itimer)
  if (verbose > 1) &
    write (io_unit%log,*) 'task_mesg_t%find_task: id =', id
  call self%lock%set ('find_task')
  new = .true.
  link => self%head
  do while (associated(link))
    task => link%task
    if (verbose > 2) &
      write (io_unit%log,*) 'task_mesg_t%find_task: task%id =', task%id
    if (task%id == id) then
      if (verbose > 1) then
        write (io_unit%log,*) 'task_mesg_t%find_task: match', id
        flush(io_unit%log)
      end if
      new = .false.
      exit
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! If new, a virtual task needs to be created
  !-----------------------------------------------------------------------------
  if (new) then
    allocate (link)
    allocate (exper)
    link%task => exper
    exper%link => link
    exper%id = id
    exper%box = self%size
    !---------------------------------------------------------------------------
    ! The virtual task needs to have all of the capabilities and procedures of the
    ! original boundary task, since it may need to play a role in extras features.
    !---------------------------------------------------------------------------
    call exper%init
    call exper%set (bits%virtual)
    if (verbose > 0) &
      write (io_unit%log,*) &
        'task_mesg_t%find: created new task, id =', id, exper%id
  end if
  call self%lock%unset ('find_task')
  call trace%end (itimer)
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
  !-----------------------------------------------------------------------------
  ! Measure the update cadence for this virtual task
  !-----------------------------------------------------------------------------
  block
  real(8):: wc
  wc = wallclock()
  task%update_cadence = wc - task%update_last
  task%update_last = wc
  end block
  !-----------------------------------------------------------------------------
  ! Unpack a patch message (which includes swapping the roles of boundary bits).
  ! Since an already existing patch may, at any one time, be under investigation
  ! by check_nbors, it must be protected by a critical region (or an OMP
  ! lock) while it is being updated here
  !-----------------------------------------------------------------------------
  id = task%id
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
      write (io_unit%log,'(f12.6,2x,a,i6)') &
        wallclock(), 'unpack: suicide note received for id =', task%id
    end if
    call self%remove_and_reset (link1)
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
    ! If link1 has no nbors it is a newly created virtual task. Does it need
    ! an nbor list?  At least we can use the nbor list to check that link1 is
    ! in its nbors nbor lists.  A new virtual task (where no task existed) means
    ! that some nbor of it has changed from internal to boundary, which will be
    ! checked by the test_nbor_status call below, but only if an nbor list exists.
    !---------------------------------------------------------------------------
    !else if (.not.associated(link1%nbor)) then
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
      write (io_unit%mpi,'(f12.6,2x,a,3i6,2i4,f12.6,2x,3f10.5,2x,"new")') &
        wallclock(), 'task_mesg_t%unpack: ', id, task%id, link1%task%id, &
        mesg%seq, omp%thread, task%time, task%position
      flush(io_unit%mpi)
    end if
    call self%add_new_link (link1)
    !---------------------------------------------------------------------------
    ! If the task has the %support bit set, call refine_t%check_supprt on it
    !---------------------------------------------------------------------------
    if (task%is_set (bits%support)) then
      call refine%check_support (self, link1, n_added)
    end if
  else
    !-----------------------------------------------------------------------------
    ! Patches that have the init_nbors bit set should have new nbor lists, in 
    ! order to replicate the behavior on the owner rank
    !-----------------------------------------------------------------------------
    if (task%is_set (bits%init_nbors)) then
      call self%init_nbors (link1)
      if (verbose > 1) &
        write (io_unit%log,*) 'unpack: bits%init_nbors id =', task%id
    end if
    if (verbose > 0) then
      write (io_unit%mpi,'(f12.6,2x,a,i6,2i4,f12.6,2x,3f10.5,2x)') &
        wallclock(), 'task_mesg_t%unpack: ', id, mesg%seq, omp%thread, &
        task%time, task%position
      flush(io_unit%mpi)
    end if
  end if
  !!omp end critical (unpack_cr)
  !-----------------------------------------------------------------------------
  ! As the task has now been updated, we need to check if any of the neighbors
  ! have become ready to update. If the task has just been swapped into being
  ! a boundary patch, it could possibly also be ready to update.
  ! FIXME: Check how this applies to other dispatcher methods than method=0
  !-----------------------------------------------------------------------------
  if (self%method==0) then
    if (verbose > 1) &
      write (io_unit%log,*) 'unpack: check_nbors, id =', task%id
    call self%check_nbors (link1)
  end if
  call trace%end (itimer)
END SUBROUTINE unpack

END MODULE task_mesg_mod
