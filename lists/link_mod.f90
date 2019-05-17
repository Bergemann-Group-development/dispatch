!*******************************************************************************
!> Module with list handling for generic class task_t objects.
!>
!> NOTE: The advantages of using the same links to organize a ready queue and
!> an active queue is that the process is faster and simpler when no new links
!> have to be created
!*******************************************************************************
MODULE link_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE omp_timer_mod
  USE omp_lock_mod
  implicit none
  private

  type, public:: link_t
    class(task_t), pointer :: task => null()          ! task
    class(link_t), pointer :: nbor => null()          ! next nbor; wouldn't `nbor_head` be a better name?
    class(link_t), pointer :: link => null()          ! pointer from nbor link to task list link (can be useful)
    class(link_t), pointer :: next => null()          ! next link in task list
    class(link_t), pointer :: prev => null()          ! previous link in task list
    class(link_t), pointer :: parent => null()        ! parent task link
    class(link_t), pointer :: next_time => null()     ! next link in time order
    class(link_t), pointer :: next_active => null()   ! next link to an active process, in time order
    class(link_t), pointer :: nbors_by_level => null()! decreasing level sorted list
    type(lock_t):: lock
    logical:: initialized = .false.
    logical:: check = .true.
    !---------------------------------------------------------------------------
    ! Flags making it possible to distinguish between nbors needed for link%task,
    ! and nbors that need link%task.  These flags are only used in nbor lists,
    ! so when parsing the list starting with link%nbor, one uses nbor%needed
    ! in check_ready(), and nbor%needs_me in check_nbors().
    !---------------------------------------------------------------------------
    logical:: needed   = .true.     ! marks nbors included in check_nbors()
    logical:: needs_me  = .true.    ! marks nbors included in check_ready()
    logical:: download = .true.     ! marks nbors included in download_link()
  contains
    procedure:: init
    procedure:: add_nbor_by_quality
    procedure:: add_nbor_by_rank
    procedure:: remove_nbor
    procedure:: remove_from_nbors
    procedure:: remove_nbor_list
    procedure:: copy_nbor_list
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    !procedure:: update_nbors
    procedure:: update_nbor_status
    procedure:: update_status
    procedure:: give_to
    procedure:: send_to_vnbors
    procedure:: log_nbors
    procedure:: print_nbors
    procedure:: test_nbor_status
    procedure:: check_level_sort
    procedure, nopass:: garbage_collect
    procedure, nopass:: garbage_remove
    procedure:: decrement_needed
    procedure:: increment_needed
    procedure:: delete
    procedure:: sort_nbors_by_level
  end type
  type(link_t), target, public:: garbage
  integer:: garbage_n=0, verbose=1
  logical, save:: debug=.false.
CONTAINS

!===============================================================================
!> If link is not already initialized, initialize its lock.  If an init is at
!> all made (not done for nbor links), it is done by the thread that allocated
!> the link, so there should be no need for a crtical region.
!===============================================================================
SUBROUTINE init (self)
  class(link_t):: self
  !.............................................................................
  if (.not.self%initialized) then
    self%initialized = .true.
    call self%lock%init ('link')
  end if
END SUBROUTINE init

!===============================================================================
!> Initialize the nbors chain of a new task (self), by
!> 1) adding the parent
!> 2) searching the parent nbors, adding tasks that overlap
!> Each time it adds an nbor, it should also add itself to that nbor's nbor list
!===============================================================================
SUBROUTINE make_new_nbors (self, selfp, parent)
  class(link_t):: self
  class(link_t), pointer:: selfp, parent, nbor
  class(task_t), pointer:: t1, t2
  !> ...........................................................................
  call trace_begin('list_t%make_new_nbors', 1)
  t1 => self%task
  call self%make_new_nbor (selfp, parent)       ! make me a neighbor of my parent
  call parent%lock%set ('make_new_nbors')
  nbor => parent%nbor                           ! start on parent's nbors chain
  do while (associated(nbor))
    t2 => nbor%task
    if (.not.associated(self%task,nbor%task) &
      .and. t1%overlaps (t2)) &   ! if my task overlaps with that task
      call self%make_new_nbor (selfp, nbor%nbor)! make a neighbor (and make me one)
    nbor => nbor%next                           ! continue on parents nbors chain
  end do
  call parent%lock%unset ('make_new_nbors')
  call trace_end
END SUBROUTINE make_new_nbors

!===============================================================================
!> Add a link to a new neighbor to my nbors chain, and add myself into that
!> neighbor's neighbor chain.  An nbor chain link has three pointers: a pointer
!> to its task, a pointer to the next link in the chain, and a pointer back to
!> the task list link to the task.  We need to make and add two now nbor chain
!> links; one for our own nbors chain, and one for that's nbors chain.  Both
!> are added in quality order
!===============================================================================
SUBROUTINE make_new_nbor (self, selfp, that)
  class(link_t):: self
  class(link_t), pointer:: selfp, that
  class(link_t), pointer:: link1, link2
  !> ...........................................................................
  call trace_begin('list_t%make_new_nbor', 1)
  allocate (link1)                              ! for our own nbors chain
  link1%task => that%task                       ! the link points to that task
  link1%nbor => that                            ! pointer to the nbor link
  link1%link => that                            ! pointer to the nbor link
  link1%parent => that%parent                   ! direct link to task parent
  call self%add_nbor_by_quality (link1)         ! add to that chain in quality order
  allocate (link2)                              ! for that's nbors chain
  link2%task => self%task                       ! that's new link points to my task
  link2%nbor => selfp                           ! pointer back to self
  link2%link => selfp                           ! pointer back to self
  link2%parent => selfp%parent                  ! direct link to self parent
  call that%add_nbor_by_quality (link2)         ! add to my nbors chain in quality order
  call trace_end
END SUBROUTINE make_new_nbor

!===============================================================================
!> Add a link node into the nbors list of self, in increasing quality order.
!> self must be a link whose 'nbors' points to a chain of nbors
!===============================================================================
SUBROUTINE add_nbor_by_quality (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_nbor_by_quality')
  call self%lock%set ('add_nbor_by_quality')
  nullify (prev)
  next => self%nbor
  do while (associated(next))
    if (.not.associated(next%task,this%task)) then
      if (next%task%level >= this%task%level) then
        !-------------------------------------------------------------------------
        ! Found a task we should be ahead of, so insert here
        !-------------------------------------------------------------------------
        this%next => next
        if (associated(prev)) then
          prev%next => this
        else
          self%nbor => this
        end if
        self%task%n_nbors = self%task%n_nbors+1
        call trace_end
        goto 9
      end if
    end if
    prev => next
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! This task has the highest quality, so append it
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    prev%next => this
    nullify (this%next)
  else
    self%nbor => this
  end if
  self%task%n_nbors = self%task%n_nbors+1
9 continue
  call self%lock%unset ('add_nbor_by_quality')
  call trace_end
END SUBROUTINE add_nbor_by_quality

!===============================================================================
!> Give a task to another rank, making sure to update the nbor relations on
!> this rank -- doing the same on the receiver rank is handled when unpacking,
!===============================================================================
SUBROUTINE give_to (self, link, rank)
  class(link_t):: self
  class(link_t), pointer:: link
  integer:: rank
  !-----------------------------------------------------------------------
  ! Set the bits and rank as they should be on the local rank, and mark
  ! the task with a swap_request bit
  !-----------------------------------------------------------------------
  call self%task%clear (bits%boundary)
  call self%task%set (bits%virtual+bits%swap_request)
  self%task%rank = rank
  !-----------------------------------------------------------------------
  ! With the new rank and status settings, update the status of the nbors
  ! and their nbor lists, which may all be affected by the change of rank.
  ! The ones that are new boundary tasks will be sent to nbor ranks by
  ! the update_nbor_status procedure.  The rank itself will be sent by
  ! task_list_mod::update, and should not be sent from here.
  !-----------------------------------------------------------------------
  call self%update_nbor_status (link)
  if (io%verbose>1) &
    call test_nbor_status (link)
END SUBROUTINE give_to

!===============================================================================
!> Check the status of a link and all nbor links
!===============================================================================
SUBROUTINE test_nbor_status (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  !.............................................................................
  write (io_unit%log,*) 'test_nbor_status: link', self%task%id
  call test_status (self)
  nbor => self%nbor
  do while (associated(nbor))
    call test_status (nbor%link)
    nbor => nbor%next
  end do
END SUBROUTINE test_nbor_status

!===============================================================================
!> Check the status of a link and all nbor links
!===============================================================================
SUBROUTINE check_level_sort (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  integer:: level
  !.............................................................................
  nbor => self%nbor
  level = -1
  do while (associated(nbor))
    if (nbor%task%level < level) then
       print *, self%task%id, 'WARNING: nbors not sorted by level'
       nbor => self%nbor
       do while (associated(nbor))
         print *, 'nbor, level =', nbor%task%id, nbor%task%level
         nbor => nbor%next
       end do
       return
    end if
    level = nbor%task%level
    nbor => nbor%next
  end do
END SUBROUTINE check_level_sort

!===============================================================================
!> Check the status of a link task  with respect to its nbors
!===============================================================================
SUBROUTINE test_status (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  integer:: status, test, n
  !.............................................................................
  write (io_unit%log,*) 'test_status: link', self%task%id
  if (self%task%rank == mpi%rank) then
    status = bits%internal
  else
    status = bits%external
  end if
  nbor => self%nbor
  do while (associated(nbor))
    if      (nbor%task%rank /= mpi%rank .and. status==bits%internal) then
      status = bits%boundary
    else if (nbor%task%rank == mpi%rank .and. status==bits%external) then
      status = bits%virtual
    end if
    nbor => nbor%next
  end do
  test = iand(self%task%status,bits%internal+bits%boundary+bits%virtual+bits%external)
  if (status /= test) then
    write (io_unit%log,'(i9,2x,a,i6,2(2x,4l1))') self%task%id, 'inconsistent:', &
      self%task%rank, &
      self%task%is_set(bits%internal), &
      self%task%is_set(bits%boundary), &
      self%task%is_set(bits%virtual) , &
      self%task%is_set(bits%external), &
      iand(status,bits%internal)/=0, &
      iand(status,bits%boundary)/=0, &
      iand(status,bits%virtual )/=0, &
      iand(status,bits%external)/=0
    n = 0
    nbor => self%nbor
    do while (associated(nbor))
      n = n+1
      write (io_unit%log,'(2i9,2x,a,2i6)') self%task%id, nbor%task%id, &
        'inconsistent ranks', self%task%rank, nbor%task%rank
      nbor => nbor%next
    end do
    if (n==0) write (io_unit%log,*) self%task%id, ' inconsistent bits: has no nbors'
  end if
END SUBROUTINE test_status

!===============================================================================
!> Add a link node into the nbors list of self, in increasing rank order.
!> self must be a link whose 'nbors' points to a chain of nbors
!===============================================================================
SUBROUTINE add_nbor_by_rank (self, nbors, this)
  class(link_t):: self
  class(link_t), pointer:: nbors, this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_nbor_by_rank', 3)
  call self%lock%set ('add_nbor_by_rank')
  if (io%verbose>2) &
  write (io_unit%log,*) '  adding task', this%task%id, &
    this%task%is_set(bits%boundary), this%task%is_set(bits%virtual), &
    ' to nbor list of', self%task%id
  nullify (prev)
  next => nbors
  do while (associated(next))
    if (.not.associated(next%task,this%task)) then
      if (next%task%rank >= this%task%rank) then
        !-------------------------------------------------------------------------
        ! Found a task we should be ahead of, so insert here
        !-------------------------------------------------------------------------
        this%next => next
        if (associated(prev)) then
          prev%next => this
        else
          nbors => this
        end if
        self%task%n_nbors = self%task%n_nbors+1
        goto 9
      end if
    end if
    prev => next
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! This task has the highest rank, so append it
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    prev%next => this
    nullify (this%next)
  else
    nbors => this
  end if
  self%task%n_nbors = self%task%n_nbors+1
9 continue
  call self%lock%unset ('add_nbor_by_rank')
  call trace_end
END SUBROUTINE add_nbor_by_rank

!===============================================================================
!> Remove this from the nbor list, and return the nbor pointer that was pointing
!> to it
!===============================================================================
SUBROUTINE remove_nbor (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: nbor, prev, next
  !.............................................................................
  call trace_begin ('link_mod::remove_nbor', 2)
  call self%lock%set ('remove_nbor')
  if (io%verbose>1) &
    write (io_unit%log,*) 'removing task', this%task%id,associated(this%link), &
      ' from nbor list of', self%task%id
  nullify (prev)
  nbor => self%nbor
  do while (associated(nbor))
    next => nbor%next
    if (associated(nbor%task,this%task)) then
      if (associated(prev)) then
        prev%next => nbor%next
      else
        self%nbor => nbor%next
      end if
      deallocate (nbor)
      self%task%n_nbors = self%task%n_nbors-1
      goto 9
    end if
    prev => nbor
    nbor => next
  end do
  if (io%verbose>2) &
    write (io_unit%log,*) 'WARNING: remove_nbor could not find task', this%task%id
9 continue
  call self%lock%unset ('remove_nbor')
  call trace_end
END SUBROUTINE remove_nbor

!===============================================================================
!> Remove this task from the nbor lists of the nbors of self
!===============================================================================
SUBROUTINE remove_from_nbors (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: nbor
  !.............................................................................
  call trace_begin ('link_mod::remove_from_nbors', 2)
  call self%lock%set ('remove_from_nbors')
  nbor => self%nbor
  do while (associated(nbor))
    call nbor%link%remove_nbor (this)
    nbor => nbor%next
  end do
  call self%lock%unset ('remove_from_nbors')
  call trace_end
END SUBROUTINE remove_from_nbors

!===============================================================================
!> Remove, deallocate, and nullify an existing nbor list
!===============================================================================
SUBROUTINE remove_nbor_list (self, nbors)
  class(link_t):: self
  class(link_t), pointer:: nbors, nbor, next
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%remove_nbor_list')
  nbor => nbors
  do while (associated(nbor))
    next => nbor%next
    deallocate (nbor)
    nbor => next
  end do
  nullify (nbors)
  call trace%end
END SUBROUTINE remove_nbor_list

!===============================================================================
!> Create a temporary nbor list
!===============================================================================
SUBROUTINE copy_nbor_list (self, old_head, new_head)
  class(link_t):: self
  class(link_t), pointer:: old_head, new_head
  !.............................................................................
  class(link_t), pointer:: nbor, head, tail, new
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%copy_nbors')
  nullify (head, tail)
  nbor => old_head
  do while (associated(nbor))
    allocate (new, source=nbor)
    !---------------------------------------------------------------------------
    ! If the new list is empty the new link becomes the first one, and if not,
    ! it becomes the next one after the last one
    !---------------------------------------------------------------------------
    if (associated(tail)) then
      tail%next => new
    else
      head => new
    end if
    tail => new
    nbor => nbor%next
  end do
  new_head => head
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE copy_nbor_list

!===============================================================================
!> Update the nbor list of self, adding self to the each nbor's nbor list
!===============================================================================
! SUBROUTINE update_nbors (self)
  ! class(link_t):: self
  ! class(link_t), pointer:: nbor, new
  ! !.............................................................................
  ! call trace_begin('link_t%update_nbors')
  ! nbor => self%nbor                                     ! start checking nbors
  ! do while (associated(nbor))                           ! loop over nbor list
    ! allocate (new)                                      ! new link for nbor's nbor list
    ! new%task => self%task                               ! its task is the link task
    ! new%link => self%task%link                          ! add also the link pointer
    ! call nbor%add_nbor_by_rank (self%nbor, new)         ! add a new nbor link
    ! nbor => nbor%next                                   ! continue on nbor list
  ! end do
  ! call trace_end
! END SUBROUTINE update_nbors

!===============================================================================
!> Update the nbor/vnbor status of a links nbor list, including the self link.
!> If the boundary bit gets set, send the task to its vnbors.  Before doing so,
!> the nbor list of the nbors need to be updated, to reflect the new task rank
!===============================================================================
SUBROUTINE update_nbor_status (self, link)
  class(link_t):: self
  class(link_t), pointer:: link, nbor, this
  integer:: status
  !-----------------------------------------------------------------------------
  call trace_begin ('link_mod::update_nbor_status')
  call self%lock%set ('update_nbor_status')
  !
  nbor => self%nbor                                     ! then help others
  do while (associated(nbor))
    if (io%verbose>1) &
      write (io_unit%log,*) 'update_nbor_status:', nbor%task%id, &
      nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual)
    status = nbor%task%status                           ! current status
    call nbor%link%remove_nbor (link)                   ! remove out-of-order
    allocate (this)                                     ! new nbor link
    this%link => link                                   ! attach task link
    this%task => link%task                              ! attach task
    call nbor%link%add_nbor_by_rank (self%nbor, this)   ! add in rank order
    call nbor%link%update_status                        ! new status
    if (nbor%task%is_set (bits%boundary).and.iand(status,bits%boundary)==0) then
      if (io%verbose>0) &
        write (io_unit%log,*) 'sending new bdry task',nbor%task%id,' to vnbors'
      if (.not.nbor%task%is_clear (bits%virtual)) then
        call nbor%task%clear (bits%virtual)
        write (io_unit%log,*) 'WARNING: needed to clear virtual bit(1)'
      end if
      call nbor%task%set (bits%swap_request)
      nbor%task%rank = mpi%rank                         ! since it's a boundary
      call nbor%link%send_to_vnbors
      call nbor%link%task%clear (bits%swap_request + &  ! sometimes set -- make
                                 bits%init_nbors)       ! sure clear after send
    end if
    nbor => nbor%next                                   ! check next nbor
  end do
  !
  if (io%verbose>0) &
    write (io_unit%log,*) 'update_nbor_status: id =', self%task%id, associated(link)
  status = self%task%status                             ! current status
  call self%update_status                               ! new status
  if (self%task%is_set (bits%boundary).and.iand(status,bits%boundary)==0) then
    if (io%verbose>0) &
      write (io_unit%log,*) 'sending new bdry task',self%task%id,' to vnbors'
    if (.not.nbor%task%is_clear (bits%virtual)) then
      write (io_unit%log,*) 'WARNING: needed to clear virtual bit(2)'
      call nbor%task%clear (bits%virtual)
    end if
    call nbor%task%set (bits%swap_request)
    nbor%task%rank = mpi%rank                           ! since it's a boundary
    call self%send_to_vnbors                            ! need to send ...
    call nbor%task%clear (bits%swap_request)
  end if
  call self%lock%unset ('update_nbor_status')
  call trace_end
END SUBROUTINE update_nbor_status

!===============================================================================
!> Update the nbor/vnbor status of a link.  If self belongs to the local rank
!> it is either an internal or boundary rank, if not it is virtual or external,
!> depending on what nbors it has.
!===============================================================================
SUBROUTINE update_status (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call self%lock%set ('update_status')
  nbor => self%nbor
  if      (self%task%rank == mpi%rank) then
    do while (associated(nbor))
      if (nbor%task%rank /= mpi%rank)  then
        call self%task%set   (bits%boundary)
        call self%task%clear (bits%internal+bits%external+bits%virtual)
        if (io%verbose>0) &
          write (io_unit%log,*) ' set to boundary id =', self%task%id
        return
      end if
      nbor => nbor%next
    end do
    call self%task%set   (bits%internal)
    call self%task%clear (bits%boundary+bits%external+bits%virtual)
    if (io%verbose>0) &
      write (io_unit%log,*) ' set to internal id =', self%task%id, associated(self%nbor)
  else
    do while (associated(nbor))
      if (nbor%task%rank == mpi%rank)  then
        call self%task%set   (bits%virtual)
        call self%task%clear (bits%internal+bits%external+bits%boundary)
        if (io%verbose>0) &
          write (io_unit%log,*) ' set to virtual  id =', self%task%id
        return
      end if
      nbor => nbor%next
    end do
    call self%task%set   (bits%external)
    call self%task%clear (bits%internal+bits%boundary+bits%virtual)
   !call self%remove_from_nbors
    if (io%verbose>0) &
      write (io_unit%log,*) '    set to external id =', self%task%id
  end if
  call self%lock%unset ('update_status')
END SUBROUTINE update_status

!===============================================================================
!> Print a list of nbors to the rank-thread log file
!===============================================================================
SUBROUTINE log_nbors (self, label)
  class (link_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 1) return
  if (present(label)) then
    if (trim(label) /= '') write (io_unit%log,'(a)') label
  end if
  write (io_unit%log,'(a,i5,2x,a,f14.8,2x,a,l3,2x,a,i3,2x,a,$)') &
    'task', self%task%id, 'time', self%task%time, &
    'ready', self%task%is_set(bits%ready), 'level', self%task%level, 'nbors:'
  if (io%verbose>4) write (io_unit%log,*) ''
  link => self%nbor
  do while (associated(link))
    if (io%verbose>4) then
      write (io_unit%log,'(i7,f14.8,l3)') &
        link%task%id, link%task%time, link%task%is_set(bits%ready)
    else
      write (io_unit%log,'(i4,l2,$)') link%task%id, link%task%is_set(bits%ready)
    end if
    link => link%next
  end do
  write (io_unit%log,*) ''
END SUBROUTINE log_nbors

!===============================================================================
!> Print a list of nbors to stdout
!===============================================================================
SUBROUTINE print_nbors (self, label)
  class (link_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: nbor
  integer:: idmax
  character(len=16):: fmt
  !.............................................................................
  if (present(label)) write (io_unit%output,'(a,$)') label
  write (io_unit%output, '(i6,l1," nbors:",$)') &
    self%task%id, self%task%is_set (bits%ready)
  idmax = 0
  nbor => self%nbor
  do while (associated(nbor))
    idmax = max(idmax, nbor%task%id)
    nbor => nbor%next
  end do
  idmax = floor(log10(real(idmax)))+2
  write (fmt, '("(i",i1,",2l1,$)")') idmax
  nbor => self%nbor
  do while (associated(nbor))
    write (io_unit%output,fmt) nbor%task%id, nbor%download, &
      nbor%task%is_ahead_of(self%task)
    if (.not.associated(nbor%link)) &
      print *, 'link_t%print_nbors: WARNING, nbor%link not associated for task', &
        nbor%task%id
    nbor => nbor%next
  end do
  write (io_unit%output,*) ''
END SUBROUTINE print_nbors

!===============================================================================
!> Send a package with updated time slice info to neighbor ranks.
!>
!> 1) assemble task data into a message
!> 3) send the message to the nbors rank (only once per rank)
!> 4) add the message to a sent_list
!> 5) check the sent list for completed sends
!> 6) remove the message from the sent_list and disassemble the message
!===============================================================================
SUBROUTINE send_to_vnbors (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  class(mesg_t), pointer:: mesg
  logical:: flag
  integer:: ierr, rank, nsend
  integer, dimension(:), allocatable:: ranks
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('link_mod::send_to_vnbors', itimer=itimer)
  call self%lock%set ('send_to_vnbors')
  if (self%task%id == io%id_debug) &
    write(io%output,*) 'DBG link_t%send_to_vnbors: id, rank =', &
      self%task%id, mpi%rank
  !-----------------------------------------------------------------------------
  ! Send to all unique ranks in the nbor list (which is sorted by rank).  The
  ! task%pack procedure is only called once (at the first send), and each new
  ! send adds requests to its reqs() array of send requests.
  !-----------------------------------------------------------------------------
  rank = -1                                                     ! previous rank
  nbor => self%nbor                                             ! first nbor
  flag = .false.
  nsend = 0
  allocate (ranks(mpi%size))
  do while (associated(nbor))                                   ! until end
    task => nbor%task
    if (self%task%id == io%id_debug) &
      write(io%output,*) 'DBG link_t%send_to_vnbors: nbor, rank =', &
          task%id, task%rank
    if (verbose > 1 .and. task%rank /= mpi%rank) then
      write (io_unit%log,'(f12.6,2x,a,2i6,2i4)') wallclock(), &
        'send_to_vnb: id,nbor,rank =', self%task%id, task%id, task%rank, rank
      flush (io_unit%log)
    end if
    if (task%rank/=mpi%rank .and. task%rank/=rank) then         ! new rank?
      rank = nbor%task%rank                                     ! current rank
      if (.not.flag) then                                       ! if not already packed
        task => self%task                                       ! convert to task
        call task%pack (mesg)                                   ! pack into mesg
        flag = .true.                                           ! mesg has been packed
      end if
      if (self%task%id == io%id_debug) &
        write(io%output,*) 'DBG link_t%send_to_vnbors: rank =', nbor%task%rank
      call mesg%send (rank)                                     ! send it
      nsend = nsend+1
      if (verbose > 0) then
        write (io_unit%log,'(f12.6,2x,a,i9,1p,e18.6,2x,a,i5,2x,5l1)') &
          wallclock(), 'sent', mesg%id, self%task%time, 'to', rank, &
          self%task%is_set (bits%internal), &
          self%task%is_set (bits%boundary), &
          self%task%is_set (bits%virtual), &
          self%task%is_set (bits%external), &
          self%task%is_set (bits%swap_request)
        flush (io_unit%log)
      end if
      ranks(nsend) = rank
    end if
    nbor => nbor%next                                           ! next nbor
  end do
  if (io%log_sent > 0) then
    !$omp critical (log_sent_cr)
    write (io_unit%sent,'(f16.6,i4,2x,a,i6,f16.6,8i5)') wallclock(), omp%thread, 'snt', self%task%id, &
      self%task%time, nsend, ranks(1:nsend)
    flush (io_unit%sent)
    !$omp end critical (log_sent_cr)
  end if
  if (verbose > 1) &
    write (io_unit%output,*) 'sent patch', self%task%id,' to',nsend,' other ranks'
  !-----------------------------------------------------------------------------
  ! Add the message to mesg%sent_list, if it was sent, and clear the swap bit
  !-----------------------------------------------------------------------------
  if (flag) call mpi_mesg%sent (mesg)                           ! add to sent_list
  call self%task%clear (bits%swap_request)
  call trace_end (itimer)
  call self%lock%unset ('send_to_vnbors')
END SUBROUTINE send_to_vnbors

!===============================================================================
!> Collect garbage, and check if there is garbage to remove.  This needs to be
!> done under lock protection, since many threads may be utilizing the GC at the
!> same time.
!===============================================================================
SUBROUTINE garbage_collect (link)
  class(link_t), pointer:: link, nbor
  !.............................................................................
  call trace%begin ('link_t%garbage_collect')
  !$omp atomic update
  link%task%n_needed = link%task%n_needed - 1   ! because linked task derefined
  !-----------------------------------------------------------------------------
  ! Immediately decrement the n_needed of nbors_by_level and remove that list.
  ! The task itself will be removed when the n_needed count drops to zero.
  !-----------------------------------------------------------------------------
  call link%decrement_needed (link%nbors_by_level)
  call link%remove_nbor_list (link%nbors_by_level)
  !-----------------------------------------------------------------------------
  ! If no other task needs this task remove it immediately, else add to garbage
  !-----------------------------------------------------------------------------
  if (link%task%n_needed <= 0) then
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,2x,a,i4,2l4)') &
        wallclock(), 'garbage_collect deleted: id =', link%task%id, &
        link%task%n_needed, link%task%is_set(bits%virtual)
    end if
    !---------------------------------------------------------------------------
    ! At this point, the link task has been removed from the task list, and
    ! will never be the target of a download again, and hence will not be
    ! needing the tasks on its nbor list.  Their n_needed counts should thus
    ! be decremented, and the sorted nbors list can also be removed. The
    ! actual removal of the task is delayed until no other task needs it.
    !---------------------------------------------------------------------------
    call garbage%lock%set ('garbage')
    call garbage_remove                       ! check consequences for nbors
    call garbage%lock%unset ('garbage')
    call link%delete (link)
  else
    call garbage%lock%set ('garbage')
    garbage_n = garbage_n+1
    link%next => garbage%next
    garbage%next => link
    if (verbose > 1) then
      write (io_unit%output,*) 'garbage_collect: task, n_needed, garbage_n =', &
        link%task%id, link%task%n_needed, garbage_n
      flush (io_unit%output)
    end if
    !---------------------------------------------------------------------------
    !> The tasks in the nbor list of the task pending deletion are no longer
    !> needed by this task, so their n_needed counters should be decremented
    !---------------------------------------------------------------------------
    call garbage_remove
    call garbage%lock%unset ('garbage')
  end if
  call trace%end
END SUBROUTINE garbage_collect

!===============================================================================
!> Check if there are links that may be deleted.  The GC list with one member 
!> looks like so: garbage%next -> link, link%next -> null()
!===============================================================================
SUBROUTINE garbage_remove
  class(link_t), pointer:: link, next, prev
  !.............................................................................
  call trace%begin ('link_t%garbage_remove')
  link => garbage%next
  prev => garbage                             ! in case we remove the 1st link
  do while (associated(link))
    next => link%next
    if (link%task%n_needed <= 0) then
      if (verbose > 0) then
        write (io_unit%log,'(f12.6,2x,a,i4,l4)') &
          wallclock(), 'garbage_remove deleted: id =', link%task%id, &
          link%task%n_needed, link%task%is_set(bits%virtual)
      end if
      garbage_n = garbage_n-1
      call link%delete (link)                 ! deallocate & delete task + link
      prev%next => next                       ! remove from GC list
    else
      prev => link                            ! advance prev if not removed
    end if
    link => next
  end do
  call trace%end
END SUBROUTINE garbage_remove

!===============================================================================
!> Decrement the n_needed counters of the nbor tasks
!===============================================================================
SUBROUTINE decrement_needed (self, nbors)
  class(link_t):: self
  class(link_t), pointer:: nbors, nbor
  !.............................................................................
  call trace%begin ('link_t%decrement_needed')
  nbor => nbors
  do while (associated(nbor))
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed - 1
    if (verbose > 1) then
      write (io_unit%output,*) 'decrement_needed: self, task, n_needed =', &
        self%task%id, nbor%task%id, nbor%task%n_needed
      flush (io_unit%output)
    end if
    nbor => nbor%next
  end do
  call trace%end
END SUBROUTINE decrement_needed

!===============================================================================
!> Increment the n_needed counters of the nbor tasks
!===============================================================================
SUBROUTINE increment_needed (self, nbors)
  class(link_t):: self
  class(link_t), pointer:: nbors, nbor
  !.............................................................................
  call trace%begin ('link_t%increment_needed')
  nbor => nbors
  do while (associated(nbor))
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed + 1
    if (verbose > 1) then
      write (io_unit%output,*) 'increment_needed: self, task, n_needed =', &
        self%task%id, nbor%task%id, nbor%task%n_needed
      flush (io_unit%output)
    end if
    nbor => nbor%next
  end do
  call trace%end
END SUBROUTINE increment_needed

!===============================================================================
!> Delete a link, and the corresponding task.  First make sure to lock the link
!> and unset any remaining nested locks.  Then delete the task, deallocate the
!> task itself, the link nbor list.  Finally release the lock and delete the
!> link.
!===============================================================================
SUBROUTINE delete (self, link)
  class(link_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor, next
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%delete')
  !-----------------------------------------------------------------------------
  call link%lock%set ('link_t%delete')
  do while (link%lock%level > 1)
    write (io_unit%output,*) 'delete_link unlocking lock level', link%lock%level
    call link%lock%unset
  end do
  !-----------------------------------------------------------------------------
  ! Remove / deallocate the unsorted nbor list before deleting the task
  !-----------------------------------------------------------------------------
  call link%remove_nbor_list (link%nbor)
  call link%task%dealloc()
  if (verbose > 1) then
    write (io_unit%output,'(f12.5,2x,a,2i6)') &
      wallclock(), 'delete: task, garbage_n =', link%task%id, garbage_n
    flush (io_unit%output)
  end if
  deallocate (link%task)
  nullify (link%task)
  !-----------------------------------------------------------------------------
  call link%lock%unset ('link_t%delete')
  deallocate (link)
  nullify (link)
  !-----------------------------------------------------------------------------
  call trace%end()
END SUBROUTINE delete

!===============================================================================
!> Create a temporary nbor list, sorted by decreasing level
!===============================================================================
SUBROUTINE sort_nbors_by_level (self, old_head, new_head)
  class(link_t):: self
  class(link_t), pointer:: old_head, new_head
  !.............................................................................
  class(link_t), pointer:: nbor, head, prev, find
  integer:: n_nbors
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%sort_nbors_by_level')
  if (verbose > 3) call nbor_print (old_head)
  nullify (head)
  n_nbors = 0
  nbor => old_head
  do while (associated(nbor))
    call insert
    n_nbors = n_nbors+1
    nbor => nbor%next
  end do
  new_head => head
  !-----------------------------------------------------------------------------
  if (verbose > 2) call nbor_print (head)
  if (verbose > 0) call nbor_check (new_head)
  call trace%end()
contains
  !-----------------------------------------------------------------------------
  ! Find a place to prepend a clone of nbor into chain with decreasing level.
  ! As we scan through the new, sorted list, prev points to the previous link
  !-----------------------------------------------------------------------------
  subroutine insert
    nullify(prev)
    find => head
    do while (associated(find))
      if (find%task%level <= nbor%task%level) then
        call prepend
        return
      end if
      prev => find
      find => find%next
    end do
    !---------------------------------------------------------------------------
    ! Reached the end if the new nbor list, so the nbor%task should be appended,
    ! which is the same as prepending to the previous (and thus last) link
    !---------------------------------------------------------------------------
    call prepend
    return
  end subroutine insert
  !-----------------------------------------------------------------------------
  ! Clone nbor and insert after prev, or as head
  !-----------------------------------------------------------------------------
  subroutine prepend
    class(link_t), pointer:: new
    allocate (new, source=nbor)
    new%next => find
    if (associated(prev)) then
      prev%next => new
    else
      head => new
    end if
  end subroutine prepend
  !-----------------------------------------------------------------------------
  ! Print nbor list
  !-----------------------------------------------------------------------------
  subroutine nbor_print (head)
    class(link_t), pointer:: head, nbor
    !---------------------------------------------------------------------------
    write(io_unit%output,*) 'target: ', self%task%id
    write(io_unit%output,*) 'sorted nbors:'
    nbor => head
    do while (associated(nbor))
      write (io_unit%output,*) nbor%task%id, nbor%task%level
      nbor => nbor%next
    end do
  end subroutine nbor_print
  !-----------------------------------------------------------------------------
  ! Check new nbor list for size and consistency
  !-----------------------------------------------------------------------------
  subroutine nbor_check (head)
    class(link_t), pointer:: head, nbor
    integer:: n, level
    !---------------------------------------------------------------------------
    n = 0
    level = 999
    nbor => head
    do while (associated(nbor))
      if (nbor%task%level > level) then
        verbose = 5
        call nbor_print (head)
        call io%abort ('ERROR: sorted nbor list not monotonic')
      end if
      level = nbor%task%level
      n = n+1
      nbor => nbor%next
    end do
    if (n /= n_nbors) then
      write (io_unit%output,*) 'nbor_check: n, nbors =', n, n_nbors
      call io%abort ('ERROR: wrong number of nbors in sorted nbor list')
    end if
  end subroutine nbor_check
END SUBROUTINE sort_nbors_by_level

END MODULE link_mod
