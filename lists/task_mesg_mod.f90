!===============================================================================
!> Simplified task mesg handling, relying on a constant number (n_recv) of MPI
!> mesg receivers.  The check_recv() procedure runs through all of them, and
!> for each mesg received it unpacks the header to find the hask key. Then,
!> using the hash table, it checks if the task exists, and if not it creates
!> a new, skeleton task.  It then unpacks the incoming MPI message, flagging
!> new tasks witb bits%external, to choose post-processing in self%unpack().
!>
!> This is a significant simplification relative to having one recv mesg for
!> each virtual task.  Limitations on the number of open requests makes it
!> necessary to instead use a fixed number of recv messages.  Having as many
!> mesg%buffer as there are virtual tasks is an advatage, in that each MPI
!> Isend then can send the message directly into the recv buffer (if the
!> size of the buffer is less the default limit for "eager" coommication --
!> this can be controlled via I_MPI_EAGER_LIMIT or corresponding env vars).
!> Hoever, even if there are more virtual tasks than buffers, the latency
!> will still be low in most cases, with MPI instead buffering the message,
!> until the mesg%is_complete() call.
!>
!> Messages are sent with a tag = the number of one of the receivers, and
!> with sender rank set to MPI_ANY_SOURCE on the recv side.  Since the next
!> message is only requested once the previous mesg has been received, mesg
!> order reversals should not happen (will be detected if so).
!>
!> The protocal could be optimized additionally, by making it possible to
!> issue the next request immediately after receipt, before self%unpack().
!> This just requires setting aside the mesg%buffer, while allocating a new
!> buffer for receiving the next message, and dellallocating the previous one
!> after self%unpack() processing.
!===============================================================================
MODULE task_mesg_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mod
  USE mpi_lib
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
    logical:: offload=.false.
    type(mesg_list_t), pointer:: recv_list
    class(link_t), pointer:: new_link => null()
  contains
    procedure:: init
    procedure:: check_recv
    procedure:: send_to_vnbors
    procedure:: post_update
    procedure:: remove_and_reset
  end type
  !.............................................................................
  integer, save:: verbose=0
  integer, save:: max_recv=2_8**31-1          ! max number of recv messages
  integer, save:: n_recv=0                    ! number of recv messages
  integer, save:: n_boundary=0                ! number of boundary tasks
  logical, save:: lock_recv=.false.
  logical, save:: use_copy = .true.
  character(len=8), save:: method='irecv'
  public:: task2patch
  type(task_mesg_t), public:: task_mesg
CONTAINS

!===============================================================================
!> Initialize incoming message structures
!===============================================================================
SUBROUTINE init (self, name)
  class(task_mesg_t), target:: self
  character(len=*), optional:: name
  integer:: iostat, i_recv, n_virtual
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(mesg_t), pointer:: mesg
  class(experiment_t), pointer:: experiment
  character(32):: label
  logical, save:: first_time = .true.
  namelist /task_mesg_params/ verbose, method, lock_recv, use_copy, max_recv
  !-----------------------------------------------------------------------------
#ifdef MPI
  if (mpi%size == 1) &
    return
  call trace%begin ('task_mesg_t%init')
  rewind (io%input)
  read (io%input, task_mesg_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, task_mesg_params)
  end if
  max_recv = min (max_recv, mpi%max_tag)
  write (io%namelist, task_mesg_params)
  flush (io%namelist)
  !-----------------------------------------------------------------------------
  ! Determine the number of boundary and virtual tasks
  !-----------------------------------------------------------------------------
  write (stdout,*) 'number of active   tasks per rank:', self%n
  n_boundary = 0
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%boundary)) then
      n_boundary = n_boundary + 1
    end if
    link => link%next
  end do
  write (stdout,*) 'number of boundary tasks per rank:', n_boundary
  call io%namelist_warning ('task_mesg_params', iostat)
  !-----------------------------------------------------------------------------
  ! Initialize a task%mesg in all boundary and virtual tasks
  !-----------------------------------------------------------------------------
  allocate (self%recv_list)
  call self%recv_list%init ()
  call os%mem_register (storage_size(self%recv_list), 1, 'recv_list')
  n_boundary = 0
  n_virtual = 0
  link => self%head
  do while (associated(link))
    task => link%task
    if (task%is_set (bits%boundary+bits%virtual)) then
      if (.not.associated(task%mesg)) then
        !select type (task)
        !class is (experiment_t)
          call task%allocate_mesg ()
        !end select
      else
        write (stdlog,*) task%id,' mesg was allocated'
      end if
    end if
    task%tag = mod(task%id-1,mpi%max_tag) + 1
    !---------------------------------------------------------------------------
    ! Setup mesg parameters for sending in boundary tasks, based on the
    ! boundary task numbering on this rank
    !---------------------------------------------------------------------------
    if (task%is_set (bits%boundary)) then
      task%seq = -2
      n_boundary = n_boundary + 1
      task%tag = n_boundary
      if (verbose > 0) &
        write (io_unit%mpi,*) 'task_mesg%init: id,tag=', task%id, task%tag,' B'
    !---------------------------------------------------------------------------
    ! Setup mesg parameters for receiving, add mesg to the mesg_recv list,
    ! make a first request, allowing any tag and source, and allocate a new mesg
    !---------------------------------------------------------------------------
    else if (task%is_set (bits%virtual)) then
      n_virtual = n_virtual + 1
      mesg => task%mesg
      mesg%tag = mpi%max_tag - task%id
      mesg%sender = task%rank
      call self%recv_list%add (mesg)
      call new_request (mesg, task)
      if (verbose > 0) &
        write (io_unit%mpi,*) 'task_mesg%init: id,tag=', task%id, mesg%tag,' V'
      !.........................................................................
      if (io%task_logging > 2) then
        write (label,'(a,3i7)') 'initmesg', mesg%id, mesg%tag, self%recv_list%n
        call task%log (label, 3)
      end if
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Add a new task, with a mesg listening for tag==0 messages
  !-----------------------------------------------------------------------------
  call create_new (self, self%new_link)
  n_recv = n_boundary + 1
  !.............................................................................
  write (stdout,*) 'number of virtual  tasks per rank:', n_virtual
  write (stdout,*) 'number of  MPI receivers per rank:', n_recv
  write (stdout,*) 'max        MPI receivers per rank:', max_recv
#endif
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE check_ready (self, link, lock)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  logical, optional:: lock
  !.............................................................................
  if (link%task%is_set (bits%virtual)) then
    if (link%task%mesg%is_complete()) then
    end if
  else
    call self%list_t%check_ready (link, lock)
  end if
END SUBROUTINE check_ready

!===============================================================================
!> Issue a request for the next message
!===============================================================================
SUBROUTINE new_request (mesg, task)
  class(mesg_t), pointer:: mesg
  class(task_t), pointer:: task
  character(32):: label
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%new_request')
  if (associated(mesg)) then
    call mesg%request ()
    !.............................................................................
    if (io%task_logging > 2) then
      if (associated(task)) then
        write (label,'("req ",4i7)') mpi%rank, mesg%tag, task%tag
        call task%log (label, 3)
      else
        write (stdlog,*) 'new_request: task not associated'
        flush (stdlog)
      end if
    end if
  else
    write (stdlog,*) 'new_request: mesg not associated'
    flush (stdlog)
  end if
  call trace%end ()
END SUBROUTINE new_request

!===============================================================================
!> Check all instances in recv_list; each one is able to receive a package,
!> adressed to task%mesg.  Since tasks may be removed in AMR runs, the access
!> to the recv_list must be controlled by an OMP lock.
!===============================================================================
SUBROUTINE check_recv (self, n)
  class(task_mesg_t):: self
  integer, optional:: n
  !.............................................................................
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  class(mesg_t), pointer:: mesg, next
  integer, pointer:: buffer(:)
  character(32):: label
  integer:: seq, i_recv, n_unpk, wrap
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) return
  n_unpk = 0
  call trace%begin ('task_mesg_t%check_recv', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Check all ranks for messages; nornally no lock would be needed, but in case
  ! a task that belongs to the recv_list can be removed for indirect reasons,
  ! and is to be removed by some other thead, then the lock here is required.
  !-----------------------------------------------------------------------------
  call                                     self%recv_list%lock%set ('recv_list')
  mesg => self%recv_list%head
  do while (associated(mesg))
    next => mesg%next
    if (verbose > 1) then
      write (stdlog,'(a,4i8)') 'check_recv: id,sender,tag,nbuf=', &
        mesg%id, mesg%sender, mesg%tag, shape(mesg%buffer)
      flush (stdlog)
    end if
    !---------------------------------------------------------------------------
    ! Keep receiving from rank until there are no more complete messages
    !---------------------------------------------------------------------------
    if (mesg%is_complete()) then
      !-----------------------------------------------------------------------
      ! Find a link to the task, and redirect the task%mesg to point to it
      !-----------------------------------------------------------------------
      link => find_link (self, mesg)
      task => link%task
      !.......................................................................
      if (io%task_logging > 1) then
        write (label,'("recv",4i7)') mpi%rank, mesg%tag, task%tag
        call task%log (label, 2)
      end if
      !.......................................................................
      if (verbose > 0) then
        write (stdlog,'(a,4i8)') &
          'task_mesg_t%check_recv: id,rank,tag,mesg%tag=', &
          task%id, task%rank, task%tag, task%mesg%tag
        flush (stdlog)
      end if
      !-------------------------------------------------------------------------
      ! Grab a pointer to the mesg%buffer, and allocate a clone to allow request
      !-------------------------------------------------------------------------
      !buffer => mesg%buffer
      !allocate (mesg%buffer, mold=buffer)
      !-------------------------------------------------------------------------
      ! Issue a new request for the current mesg, and then unpack the buffer
      !-------------------------------------------------------------------------
      !call new_request (mesg, task)
      !call unpack (self, buffer, link)
      !deallocate (buffer)
      ! Unpack before new_request
      !-------------------------------------------------------------------------
      ! The unpack sets task%id and task%tag
      !-------------------------------------------------------------------------
      call unpack (self, mesg%buffer, link)
      n_unpk = n_unpk+1
      !-------------------------------------------------------------------------
      ! Unless the task has been deleted, or has the %remove bit set, request
      ! a new message
      !-------------------------------------------------------------------------
      if (associated(link)) then
        if (task%is_clear(bits%remove)) then
          !---------------------------------------------------------------------
          ! Make sure each virtual task has a unique tag; this is the tag used to
          ! recv the NEXT mesg; cf the similar line in send_to_vbors
          !---------------------------------------------------------------------
          wrap = mpi_mesg%tag_wrap
          mesg%tag = mod(task%seq+1,wrap) + task%tag*wrap
          mesg%sender = task%rank
          call new_request (mesg, task)
        else if (verbose > 1) then
          write (stdlog,*) 'is_set(bits%remove)'
          flush (stdlog)
        end if
      else if (verbose > 1) then
        write (stdlog,*) 'no link'
        flush (stdlog)
      end if
    end if
    mesg => next
  end do
  call                                   self%recv_list%lock%unset ('recv_list')
  if (present(n)) then
    n = n_unpk
  end if
  call trace%end (itimer)
END SUBROUTINE check_recv

!===============================================================================
!> Find the link that corresponds to an incoming message
!===============================================================================
FUNCTION find_link (self, mesg) RESULT (link)
  class(task_mesg_t):: self
  class(mesg_t), pointer:: mesg
  class(link_t), pointer:: link
  type(patch_t):: patch
  real(8):: latency
  character(32):: label
  integer:: key(4), tag, rank
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Get the task%key from the header and try to get the link from that.  If no
  ! link found, create a new task and return the link
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%find_link', itimer=itimer)
  call patch%unpack_meta (mesg, latency, key, tag, rank)
  call self%hash_table%get (key, link)
  if (io%task_logging > 3) then
    write (label,'("key",4i4,l2)') key, associated(link)
  end if
  if (associated(link)) then
    call link%task%log (label, 4)
    call link%task%clear (bits%external)
    if (verbose > 0) &
      write (stdlog,'(a,2i6)') &
      'task_mesg_t%find_link: found id=', link%task%id, mesg%id
    call link%task%log ('found', 4)
  else
    !---------------------------------------------------------------------------
    ! Return the self%new_link pointer
    !---------------------------------------------------------------------------
    link => self%new_link
    call link%task%log (label, 4)
    call link%task%set (bits%external)
    !---------------------------------------------------------------------------
    ! Create a new self%new_link as the receiver of tag==0 mesgs
    !---------------------------------------------------------------------------
    call create_new (self, self%new_link)
    if (verbose > 0) &
      write (stdlog,'(a,2i6)') &
      'task_mesg_t%find_link: created new task, id=', link%task%id
  end if
  link%task%latency = latency
  !-----------------------------------------------------------------------------
  ! Change to pointed tag and sender mesg
  !-----------------------------------------------------------------------------
  ! WRITE (stdlog,*) 'XX1', link%task%id, tag, mesg%tag
  mesg%sender = rank
  link%task%tag = tag
  call trace%end (itimer)
END FUNCTION find_link

!===============================================================================
!> Create a new skeleton virtual task, including a new task%mesg with an open
!> request.
!===============================================================================
SUBROUTINE create_new (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(halo_t), pointer:: halo
  class(task_t), pointer:: task
  class(experiment_t), pointer:: patch
  character(32):: label
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%create_new_task', itimer=itimer)
  allocate (patch)
  allocate (halo)
  link => halo
  patch%link => halo
  halo%task => patch
  patch%id = 1
  call patch%init ()
  patch%id = 0
  call halo%init ()
  patch%rank = -1
  patch%seq = 0
  patch%task_list => self%list
  patch%status = bits%virtual + bits%external
  !-----------------------------------------------------------------------------
  ! Allocate a mesg_t instance and add it to the recv_list
  !-----------------------------------------------------------------------------
  call patch%allocate_mesg ()
  task => patch
  task%mesg%tag = 0
#ifdef MPI
  task%mesg%sender = MPI_ANY_SOURCE
#endif
  call self%recv_list%add (task%mesg)
  call new_request (task%mesg, task)
  !.............................................................................
  if (io%task_logging > 2) then
    write (label,'(a,3i7)') 'create_new', &
    task%mesg%id, task%mesg%seq, self%recv_list%n
    call task%log (label, 3)
  end if
  call trace%end (itimer)
END SUBROUTINE create_new

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
  class(task_t), pointer:: task
  integer:: wrap
  !.............................................................................
  if (mpi%size == 1) return
  call trace%begin ('task_mesg_t%send_to_vnbors')
  task => link%task
  ! if (verbose > 0) call task%sent_list%check ('b')
  if (.not. associated(task%sent_list)) then
    allocate (task%sent_list)
    call task%sent_list%init
    task%sent_list%id = task%id
    call task%log ('sent_list%init')
    if (verbose > 0) &
      write (stdlog,*) &
        task%id,' has allocated a sent_list'
  end if
  !-----------------------------------------------------------------------------
  ! Allocates a new task%mesg, which is kept until all ranks have received this
  ! mesg; cf. check_sent()
  !-----------------------------------------------------------------------------
  call task%allocate_mesg ()
  task%mesg%id = task%id
  if (verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,i5,2i8,i4)') wallclock(), &
      'task_mesg_t%send2vnb: new id,sender,seq,tag,n_try =', &
      task%mesg%id, task%mesg%sender, task%mesg%seq, task%mesg%tag, &
      task%mesg%n_try
  end if
  !-----------------------------------------------------------------------------
  ! Pack task into task%mesg, making sure the mesg%tag is unique
  !-----------------------------------------------------------------------------
  if (task%seq == -2) then
    !--- tag for the first send ---
    task%mesg%tag = mpi%max_tag - task%id
    task%seq = 0
  else if (task%seq == -1) then
    !--- tag for a new AMR task ---
    task%mesg%tag = 0
    task%seq = 0
    n_boundary = n_boundary + 1
    task%tag = n_boundary
  else
    !--- tag for subsequent sends ---
    wrap = mpi_mesg%tag_wrap
    task%seq = mod(task%seq+1,wrap)
    task%mesg%tag = task%seq + task%tag*wrap
  end if
  call task%pack ()
  !-----------------------------------------------------------------------------
  ! Sending with list_t%send_to_vnbors() creates a task%mesg instance
  !-----------------------------------------------------------------------------
  call self%list_t%send_to_vnbors (link)
  call task%sent_list%add (task%mesg)
  !$omp atomic update
  mpi_mesg%nq_send = mpi_mesg%nq_send+1
  call check_sent (self, task)

  call trace%end ()
END SUBROUTINE send_to_vnbors

!===============================================================================
!> Check sent messages for this task.  This is done by the thread that owns
!> the task and just sent a message, and should thus not need any OMP locking.
!===============================================================================
SUBROUTINE check_sent (self, task)
  class(task_mesg_t):: self
  class(task_t), pointer:: task
  !.............................................................................
  class(mesg_t), pointer:: mesg, next
  character(32):: label
  !-----------------------------------------------------------------------------
  ! The function mesg%test_sent() returns true when all ranks have received
  ! the message, which can then be removed and deleted, thus reclaiming memory.
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%check_sent')
  if (io%task_logging > 2) then
    write (label,"('check_sent',i7)") task%sent_list%n
    call task%log (trim(label), 3)
    flush (io_unit%task)
  end if
  if (associated(task%sent_list)) then
    mesg => task%sent_list%head
    do while (associated(mesg))
      next => mesg%next
      if (verbose > 1) then
        write (stdlog,*) mesg%id, associated(next), task%sent_list%n
        flush (stdlog)
      end if
      if (mesg%test_sent()) then
        if (verbose > 1) then
          write (stdlog,*) mesg%id, associated(next), ' Y'
          flush (stdlog)
        end if
        call task%sent_list%remove (mesg)
        call task%sent_list%delete (mesg, send=.true.)
        nullify (mesg)
      else
        if (verbose > 1) then
          write (stdlog,*) mesg%id, associated(next), ' N'
          flush (stdlog)
        end if
      end if
      mesg => next
    end do
    !---------------------------------------------------------------------------
    ! Make sure task%mesg is not pointing to some random mesg
    !---------------------------------------------------------------------------
    task%mesg => task%sent_list%tail
  end if
  call task%log ('check_sent end', 3)
  call trace%end()
END SUBROUTINE check_sent

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
SUBROUTINE unpack (self, buffer, link)
  class(task_mesg_t):: self
  integer, pointer:: buffer(:)
  class(link_t), pointer:: link
  !.............................................................................
  class(task_t), pointer:: task
  logical:: bonafide, new
  integer:: save_bits, it, key(4), tag
  character(len=24):: label
  type(patch_t):: patch
  real(8):: wc
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%unpack', itimer=itimer)
  task => link%task
  !-----------------------------------------------------------------------------
  ! If the task about to receive the message doesn't yet have a task id, get it
  ! from the message itself, and flag the task as new
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%external)) then
    call task%clear (bits%external)
    write (label,'(a,i6,i4)') 'new  id,rank=', task%id, task%rank
    call task%log (label)
    new = .true.
  else
    new = .false.
  end if
  !-----------------------------------------------------------------------------
  bonafide = .true.
  !-----------------------------------------------------------------------------
  ! Guard against lingering extra messages to a swapped patch
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%ready)) then
    write (stdlog,'(f12.6,2x,a,i9,3x,5l1)') wallclock(), &
     'task_mesg_mod::unpack ERROR, received mpi_mesg for task with ready bit:',&
      task%id, &
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
  call task%unpack (buffer)
  task%status = ior (task%status, save_bits)
  !$omp atomic
  mpi_mesg%n_unpk = mpi_mesg%n_unpk+1
  !-----------------------------------------------------------------------------
  ! Diagnostics
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    write (stdlog,'(f12.6,3x,a,i7,":",i2.2,5i4,2x,4l1,i4,f8.4,2f11.6)') &
      wallclock(), 'unpack: id,key,nq,BVRS,nsnd,nrcv,latency,last,time =', &
      task%id, task%key, task%rank, task%nq, &
      task%is_set(bits%boundary), &
      task%is_set(bits%virtual), &
      task%is_set(bits%remove), &
      task%is_set(bits%swap_request), &
      self%recv_list%n, task%latency, task%wc_last, task%time
    if (io%verbose > 2) then
      select type (task)
      class is (patch_t)
        write (io_unit%log,*) &
          task%id, 'after unpack', minval(task%mem(:,:,:,1,task%it,1))
        end select
    end if
  end if
  !-----------------------------------------------------------------------------
  if (io%task_logging > 1) then
    write (label,'(a,i4,i8,1x,4a1)') 'unpack ', task%rank, task%istep, &
      merge('B','b',task%is_set (bits%boundary)), &
      merge('V','v',task%is_set (bits%virtual)), &
      merge('R','r',task%is_set (bits%remove)), &
      merge('S','s',task%is_set (bits%swap_request))
    call task%log (label, 3)
  end if
  !-----------------------------------------------------------------------------
  ! Since an already existing patch may, at any one time, be under investigation
  ! by check_nbors, it must be protected by a critical region (or an OMP
  ! lock) while it is being updated here
  !-----------------------------------------------------------------------------
  if (self%verbose > 0) then
    write (stdlog,'(f12.6,3x,a,i6,f12.6,i3,2x,7l1,i8,10i3)') &
      wallclock(), 'task_mesg_t%unpack: task%id, time, IBVrRSN, tag =', &
      task%id, task%time, task%it, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%ready), &
      task%is_set (bits%remove), &
      task%is_set (bits%swap_request), &
      new, task%tag, task%iit
    flush (io_unit%mpi)
  end if
  !-----------------------------------------------------------------------------
  ! Logging and verbose output
  !-----------------------------------------------------------------------------
  if (mpi_mesg%debug .or. task%id==io%id_debug) then
    !$omp critical (debug_cr)
    write (stdlog,'(f12.6,2x,a,2i9,z12,2x,5l1)') wallclock(), &
      'DBG unpk: id, sender, req =', task%id, task%rank, task%istep, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%external), &
      task%is_set (bits%swap_request)
    !$omp end critical (debug_cr)
  end if
  if (mpi_mesg%debug) &
    write (stdlog,'(f12.6,2x,a,i9,1p,e18.6)') wallclock(), &
      'unpk: id, time =', task%id, task%time
  !-----------------------------------------------------------------------------
  if (bonafide) then
    if (self%verbose > 1) write(stdlog,*) 'bonafide', new, &
      task%is_set(bits%init_nbors), task%is_set(bits%remove), &
      task%is_set(bits%swap_request)
    !---------------------------------------------------------------------------
    ! 1) Suicide note, remove the task and delete the message
    !---------------------------------------------------------------------------
    if (task%is_set(bits%remove)) then
      if (self%verbose > 0) write(stdlog,*) task%id, 'remove'
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vremove', task%id, task%time
      end if
      call self%remove_and_reset (link)
    !---------------------------------------------------------------------------
    ! 2) swap_request bit set; the task has just changed owner, and the task
    !    should become a boundary task on this rank.  Indeed, since it was
    !    sent as a virtual task, the flipping should have been done in the
    !    unpack procedure already.   Delete the message.
    !---------------------------------------------------------------------------
    else if (task%is_set(bits%swap_request)) then
      if (self%verbose > 0) write(stdlog,*) task%id, 'swap'
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vswap', task%id, task%time
      end if
      call task%log ('unpack swap')
      if (verbose > 0) then
        write(io_unit%log,*) task%id, 'swap'
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vswap', task%id, task%time
      end if
      call                             self%lock%set_read ('task_mesg_t%unpack')
      if (load_balance%is_on) &
        call                   load_balance%lock%set_read ('task_mesg_t%unpack')
      call task%log ('task_mesg_t%unpack -> list_t%swap')
      call self%swap (link, .false.)
      if (load_balance%is_on) &
        call                 load_balance%lock%unset_read ('task_mesg_t%unpack')
      call                           self%lock%unset_read ('task_mesg_t%unpack')
    !---------------------------------------------------------------------------
    ! 3) the task is a new virtual task.  It has the init_nbors bit set, and
    !    since it is a virtual task that bit should be cleared the next time
    !    the task is updated, coming back not as new task but as a virtual
    !    task update.  But first, issue a request for the next message.
    !---------------------------------------------------------------------------
    else if (new) then
      call task%log ('new', 4)
      if (self%verbose > 0) write(stdlog,*) task%id, 'new', task%it
      !-------------------------------------------------------------------------
      ! Re-initialize task, now with metadata from task%unpack(), and with
      ! optionally extra attributes added at the component level
      !-------------------------------------------------------------------------
      key = task%key
      tag = task%tag
      call task%init ()
      task%key = key
      task%tag = tag
      !-------------------------------------------------------------------------
      ! Let all time slots refer to the mem slot the new task data is in
      !-------------------------------------------------------------------------
      it = task%it
      task%t(:) = task%time
      task%dt(:) = task%dtime
      task%iit(:) = it
      !-------------------------------------------------------------------------
      ! Add the new task to the task list, and give it an nbor list, etc
      !-------------------------------------------------------------------------
      call self%add_by_level (link, count=.true.)   ! add to task list
      task%previous_status = bits%external
      call self%add_new_link (link)                 ! default actions
      !.........................................................................
      if (verbose > 0) then
        write (stdlog,'(a,i6,":",i2.2,3i4)') &
          'new task id,key:', task%id, task%key
        if (associated(link%parent)) then
          write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,3i6)') &
            wallclock(), omp%thread, 'vnew', task%id, task%time, task%istep, &
            link%parent%task%id
        else
          write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,3i6)') &
            wallclock(), omp%thread, 'vnew', task%id, task%time, task%istep
        end if
      end if
      !.........................................................................
      write (label,'(a,4l1)') 'processed BVRS=', &
        task%is_set (bits%boundary), &
        task%is_set (bits%virtual), &
        task%is_set (bits%remove), &
        task%is_set (bits%swap_request)
      call task%log (label)
    !---------------------------------------------------------------------------
    ! 4) the MPI message updated an existing task, which should be followed by
    !    exactly the same actions as on the owner rank; cf. post_update()
    !---------------------------------------------------------------------------
    else
      call task%log ('upd', 4)
      if (self%verbose > 1) write(stdlog,*) task%id, 'update'
      !-------------------------------------------------------------------------
      ! task was updated, so do the common task list post_update()
      !-------------------------------------------------------------------------
      call self%post_update (link)
      if (verbose > 0) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'vupdate', task%id, task%time
      end if
    end if
  else
    if (self%verbose > 1) write(stdlog,*) task%id, 'crap'
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

!===============================================================================
!> The common postprocessing to be done for virtual and real tasks
!===============================================================================
SUBROUTINE post_update (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%post_update', itimer=itimer)
  task => link%task
  !.............................................................................
  if (verbose > 1) &
    write (io_unit%mpi,'(a,2i6,2(2x,3l1))') & 
      'task_mesg_t%post_update: id,IBV,iu=', task%id, task%rank, &
      task%is_set (bits%internal), &
      task%is_set (bits%boundary), &
      task%is_set (bits%virtual), &
      task%is_set (bits%init_nbors), &
      task%is_set (bits%update_status)
  !-----------------------------------------------------------------------------
  ! Check if the update has cleared any of the original nbors for update
  !-----------------------------------------------------------------------------
  call self%check_nbors (link)
  !-----------------------------------------------------------------------------
  ! Set experiment properties known only at the task list level
  !-----------------------------------------------------------------------------
  call self%init_task (link%task)
  !-----------------------------------------------------------------------------
  ! Set the rank, by whatever (possibly overloaded) rule
  !-----------------------------------------------------------------------------
  call link%task%set_rank ()
  !-----------------------------------------------------------------------------
  ! If the init_nbors bit is set, make a new nbor list, which also resets the
  ! status (IBVE) of the link task
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%init_nbors)) then
    call self%init_nbors (link)
    call link%check_status ('post_update1')
    call task%clear (bits%init_nbors)
    !---------------------------------------------------------------------------
    ! The busy bit should not be cleared until after init_nbors
    !---------------------------------------------------------------------------
    call task%log ('clear busy', 3)
    call task%clear (bits%busy)
    !---------------------------------------------------------------------------
    ! Take care of tasks turned external
    !---------------------------------------------------------------------------
    if (task%is_set(bits%external)) then
      call self%remove_and_reset (link)
    else
      !-------------------------------------------------------------------------
      ! Check if the task is ready to update with the new nbor list
      !-------------------------------------------------------------------------
      if (task%is_clear (bits%virtual)) then
        call self%check_ready (link)
      end if
      !-------------------------------------------------------------------------
      ! If the task has the update_status bit set, re-evaluate status bits
      !-------------------------------------------------------------------------
      if (task%is_set (bits%update_status)) then
        call self%update_status (link)
      end if
      !-------------------------------------------------------------------------
      ! For good measure, check if any of the new nbors are ready to update
      !--------------------------------------------------------------.----------
      call self%check_nbors (link)
    end if
  else
    !---------------------------------------------------------------------------
    ! The busy bit should be cleared to allow check_ready()
    !---------------------------------------------------------------------------
    call task%log ('clear busy', 3)
    call task%clear (bits%busy)
    !---------------------------------------------------------------------------
    ! Unless the task is virtual, check if it is already OK for another update
    !---------------------------------------------------------------------------
    if (task%is_clear (bits%virtual)) then
      call self%check_ready (link)
    end if
  end if
  !call self%init_nbors2 (link)
  !-----------------------------------------------------------------------------
  ! Check if this task has changed status
  !-----------------------------------------------------------------------------
  call link%check_status ('post_update2')
  !-----------------------------------------------------------------------------
  ! Mark thread no longer on this task
  !-----------------------------------------------------------------------------
  omp%thread_task = -1
  call trace%end (itimer)
END SUBROUTINE post_update

!===============================================================================
!> If the link task is virtual and the process uses a recv_list, remove the
!> task mesg from it.  Then remove a task from consideration, making sure that
!> virtual copies do the same.  Set bits%remove, which excludes the task from
!> nbor lists in check_ready(), and signals rank nbors to repeat this
!> remove_and_reset().
!===============================================================================
SUBROUTINE remove_and_reset (self, link)
  class(task_mesg_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: count, level
  character(32):: label
  !-----------------------------------------------------------------------------
  call trace%begin ('task_mesg_t%remove_and_reset')
  task => link%task
                                        call self%qlock%set ('remove_and_reset')
  !.............................................................................
  if (verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,2l3)') &
      wallclock(), 'task_mesg_t%remove_and_reset: id,BV=', task%id, &
      task%is_set(bits%boundary), task%is_set(bits%virtual)
  end if
  !-----------------------------------------------------------------------------
  ! By setting bits%remove immediately, we make sure that the task is no longer
  ! considered in check_ready() calls, and also that the task is no longer
  ! counted in check_ready() calls to nbors of it.
  !-----------------------------------------------------------------------------
  call link%task%set (bits%remove)            ! cf. check_ready() and unpack()
  call link%task%set (bits%init_nbors)        ! trigger new nbor list
  call link%task%clear (bits%ready+bits%busy) ! just to be safe
  !-----------------------------------------------------------------------------
  ! Send a suicide note to virtual nbors, if this is a boundary patch
  !-----------------------------------------------------------------------------
  if (link%task%is_set (bits%boundary)) then
    call self%send_to_vnbors (link)
  end if
  !-----------------------------------------------------------------------------
  ! On the virtual side, remove the task%mesg from the recv_list.  On the
  ! boundary side the self%sent_list is deallocated as part of self%remove()
  !-----------------------------------------------------------------------------
  if (task%is_set(bits%virtual)) then
    if (associated(self%recv_list)) then
      if (associated(task%mesg)) then
        !.......................................................................
        if (io%task_logging > 2) then
          write (label,'(a,3i7)') 'delmesg(V)', &
            task%mesg%id, task%mesg%seq, self%recv_list%n
          call task%log (label, 3)
        end if
        !-----------------------------------------------------------------------
        ! If this procedure is executed by a non-master thread, lock/unlock
        !-----------------------------------------------------------------------
        if (omp%thread == 0) then
          call self%recv_list%remove (task%mesg)
        else
          call                      self%recv_list%lock%set ('reset_and_remove')
          call self%recv_list%remove (task%mesg)
          call                    self%recv_list%lock%unset ('reset_and_remove')
        end if
        call task%mesg%dealloc ()
        deallocate (task%mesg)
        nullify (task%mesg)
      end if
    end if
  end if
  !.............................................................................
  level = link%task%level
  !$omp atomic update
  timer%removed(level) = &
  timer%removed(level) + 1
  call self%remove (link)                  ! remove link from task list
  call self%init_nbors (link)              ! make a fresh nbor list
  call self%check_nbors (link)             ! check nbors, now link is removed
  call self%set_init_nbors (link)          ! trigger new nbor lists on nbors
  call link%garbage_collect (link)         ! add to garbage bin
                                      call self%qlock%unset ('remove_and_reset')
  call trace%end ()
END SUBROUTINE remove_and_reset

END MODULE task_mesg_mod
