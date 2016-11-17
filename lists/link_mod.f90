!> *****************************************************************************
!> Module with list handling for generic class task_t objects
!> *****************************************************************************
MODULE link_mod
  USE io_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE omp_timer_mod
  implicit none
  private

  type, public:: link_t
    class(task_t), pointer :: task => null()         ! task
    class(link_t), pointer :: nbor => null()         ! next nbor; wouldn't `nbor_head` be a better name?
    class(link_t), pointer :: link => null()         ! pointer from nbor link to task list link (can be useful)
    class(link_t), pointer :: next => null()         ! next link in task list
    class(link_t), pointer :: prev => null()         ! previous link in task list
    class(link_t), pointer :: parent => null()       ! parent task link
    class(link_t), pointer :: next_time => null()    ! next link in time order
    logical:: check = .true.
  contains
    procedure:: add_nbor_by_quality
    procedure:: add_nbor_by_rank
    procedure:: remove_nbor
    procedure:: remove_from_nbors
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    procedure:: update_nbors
    procedure:: update_nbor_status
    procedure:: update_status
    procedure:: give_to
    procedure:: send_to_vnbors
    procedure:: print_nbors
    procedure:: test_nbor_status
  end type
  logical, save:: debug=.false.
CONTAINS

!> =============================================================================
!> Initialize the nbors chain of a new task (self), by
!> 1) adding the parent
!> 2) searching the parent nbors, adding tasks that overlap
!> Each time it adds an nbor, it should also add itself to that nbor's nbor list
!> =============================================================================
SUBROUTINE make_new_nbors (self, selfp, parent)
  class(link_t):: self
  class(link_t), pointer:: selfp, parent, nbor
  class(task_t), pointer:: t1, t2
  !> ...........................................................................
  call trace_begin('list_t%make_new_nbors', 1)
  t1 => self%task
  call self%make_new_nbor (selfp, parent)       ! make me a neighbor of my parent
  nbor => parent%nbor                           ! start on parent's nbors chain
  do while (associated(nbor))
    t2 => nbor%task
    if (.not.associated(self%task,nbor%task) &
      .and. t1%overlaps (t2)) &   ! if my task overlaps with that task
      call self%make_new_nbor (selfp, nbor%nbor)! make a neighbor (and make me one)
    nbor => nbor%next                           ! continue on parents nbors chain
  end do
  call trace_end
END SUBROUTINE make_new_nbors

!> =============================================================================
!> Add a link to a new neighbor to my nbors chain, and add myself into that
!> neighbor's neighbor chain.  An nbor chain link has three pointers: a pointer
!> to its task, a pointer to the next link in the chain, and a pointer back to
!> the task list link to the task.  We need to make and add two now nbor chain
!> links; one for our own nbors chain, and one for that's nbors chain.  Both
!> are added in quality order
!> =============================================================================
SUBROUTINE make_new_nbor (self, selfp, that)
  class(link_t):: self
  class(link_t), pointer:: selfp, that
  class(link_t), pointer:: link1, link2
  !> ...........................................................................
  call trace_begin('list_t%make_new_nbor', 1)
  allocate (link1)                              ! for our own nbors chain
  link1%task => that%task                       ! the link points to that task
  link1%nbor => that                            ! link to the nbor link
  link1%parent => that%parent                   ! direct link to task parent
  call self%add_nbor_by_quality (link1)         ! add to that chain in quality order
  allocate (link2)                              ! for that's nbors chain
  link2%task => self%task                       ! that's new link points to my task
  link2%nbor => selfp                           ! pointer back to self
  link2%parent => selfp%parent                  ! direct link to self parent
  call that%add_nbor_by_quality (link2)         ! add to my nbors chain in quality order
  call trace_end
END SUBROUTINE make_new_nbor

!> =============================================================================
!> Add a link node into the nbors list of self, in increasing quality order.
!> self must be a link whose 'nbors' points to a chain of nbors
!> =============================================================================
SUBROUTINE add_nbor_by_quality (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_nbor_by_quality')
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
  call trace_end
END SUBROUTINE add_nbor_by_quality

!> =============================================================================
!> Give a task to another rank, making sure to update the nbor relations on
!> this rank -- doing the same on the receiver rank is handled when unpacking,
!> =============================================================================
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

!> =============================================================================
!> Check the status of a link and all nbor links
!> =============================================================================
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

!> =============================================================================
!> Check the status of a link task  with respect to its nbors
!> =============================================================================
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

!> =============================================================================
!> Add a link node into the nbors list of self, in increasing rank order.
!> self must be a link whose 'nbors' points to a chain of nbors
!> =============================================================================
SUBROUTINE add_nbor_by_rank (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_nbor_by_rank', 3)
  if (io%verbose>2) &
  write (io_unit%log,*) '  adding task', this%task%id, &
    this%task%is_set(bits%boundary), this%task%is_set(bits%virtual), &
    ' to nbor list of', self%task%id
  nullify (prev)
  next => self%nbor
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
          self%nbor => this
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
    self%nbor => this
  end if
  self%task%n_nbors = self%task%n_nbors+1
9 continue
  call trace_end
END SUBROUTINE add_nbor_by_rank

!> =============================================================================
!> =============================================================================
FUNCTION remove_nbor (self, this) RESULT (nbor)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: nbor, prev
  !.............................................................................
  call trace_begin ('link_mod::remove_nbor', 2)
  if (io%verbose>1) &
    write (io_unit%log,*) 'removing task', this%task%id,associated(this%link), &
      ' from nbor list of', self%task%id
  nullify (prev)
  nbor => self%nbor
  do while (associated(nbor))
    if (associated(nbor%task,this%task)) then
      if (associated(prev)) then
        prev%next => nbor%next
      else
        self%nbor => nbor%next
      end if
      self%task%n_nbors = self%task%n_nbors-1
      goto 9
    end if
    prev => nbor
    nbor => nbor%next
  end do
  if (io%verbose>2) &
    write (io_unit%log,*) 'WARNING: remove_nbor could not find task', this%task%id
9 continue
  call trace_end
END FUNCTION remove_nbor

!> =============================================================================
!> =============================================================================
FUNCTION remove_from_nbors (self, this) RESULT (nbor)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: nbor, prev
  !.............................................................................
  call trace_begin ('link_mod::remove_from_nbors', 2)
  nbor => self%nbor
  do while (associated(nbor))
    prev => nbor%link%remove_nbor (this)
    nbor => nbor%next
  end do
  call trace_end
END FUNCTION remove_from_nbors

!> =============================================================================
!> Update the nbor list of self, adding self to the each nbor's nbor list
!> =============================================================================
SUBROUTINE update_nbors (self)
  class(link_t):: self
  class(link_t), pointer:: nbor, new
  !.............................................................................
  call trace_begin('link_mod::init_all_nbors')
  nbor => self%nbor                                     ! start checking nbors
  do while (associated(nbor))                           ! loop over nbor list
    allocate (new)                                      ! new link for nbor's nbor list
    new%task => self%task                               ! its task is the link task
    call nbor%add_nbor_by_rank (new)                    ! add a new nbor link
    nbor => nbor%next                                   ! continue on nbor list
  end do
  call trace_end
END SUBROUTINE update_nbors

!> =============================================================================
!> Update the nbor/vnbor status of a links nbor list, including the self link.
!> If the boundary bit gets set, send the task to its vnbors.  Before doing so,
!> the nbor list of the nbors need to be updated, to reflect the new task rank
!> =============================================================================
SUBROUTINE update_nbor_status (self, link)
  class(link_t):: self
  class(link_t), pointer:: link, nbor, this
  integer:: status
  !-----------------------------------------------------------------------------
  call trace_begin ('link_mod::update_nbor_status')
  !
  nbor => self%nbor                                     ! then help others
  do while (associated(nbor))
    if (io%verbose>1) &
      write (io_unit%log,*) 'update_nbor_status:', nbor%task%id, &
      nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual)
    status = nbor%task%status                           ! current status
    this => nbor%link%remove_nbor (link)                ! remove out-of-order
    if (associated(this)) then
      call nbor%link%add_nbor_by_rank (this)            ! add in rank order
      if (io%verbose>2) &
        write (io_unit%log,*) nbor%task%id, associated(nbor%link), associated(this%link)
    else
      allocate (this)                                   ! new nbor link
      this%link => link                                 ! attach task link
      this%task => link%task                            ! attach task
      call nbor%link%add_nbor_by_rank (this)            ! add in rank order
    end if
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
      call nbor%task%clear (bits%swap_request)
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
  call trace_end
END SUBROUTINE update_nbor_status

!> =============================================================================
!> Update the nbor/vnbor status of a link.  If self belongs to the local rank
!> it is either an internal or boundary rank, if not it is virtual or external,
!> depending on what nbors it has.
!> =============================================================================
SUBROUTINE update_status (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
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
END SUBROUTINE update_status

!> =============================================================================
!> Print a list of nbors
!> =============================================================================
SUBROUTINE print_nbors (self, label)
  class (link_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
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
END SUBROUTINE print_nbors

!> =============================================================================
!> Send a package with updated time slice info to neighbor ranks.
!>
!> 1) assemble task data into a message
!> 3) send the message to the nbors rank (only once per rank)
!> 4) add the message to a sent_list
!> 5) check the sent list for completed sends
!> 6) remove the message from the sent_list and disassemble the message
!> =============================================================================
SUBROUTINE send_to_vnbors (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  class(mesg_t), pointer:: mesg
  logical:: flag
  integer:: ierr, rank
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('link_mod::send_to_vnbors', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Send to all unique ranks in the nbor list (which is sorted by rank).  The
  ! task%pack procedure is only called once (at the first send), and each new
  ! send adds requests to its reqs() array of send requests.
  !-----------------------------------------------------------------------------
  rank = -1                                                     ! previous rank
  nbor => self%nbor                                             ! first nbor
  flag = .false.
  do while (associated(nbor))                                   ! until end
    task => nbor%task
    if (io%verbose>2) &
      write (io_unit%log,*) 'send_to_vnb: id,nbor,rank =', self%task%id, task%id, task%rank
    flush (io_unit%log)
    if (task%rank/=mpi%rank .and. task%rank/=rank) then         ! new rank?
      rank = nbor%task%rank                                     ! current rank
      if (.not.flag) then                                       ! if not already packed
        task => self%task                                       ! convert to task
        call task%pack (mesg)                                   ! pack into mesg
      end if
      call mesg%send (rank)                                     ! send it
      flag = .true.
      if (io%verbose>1) then
        write (io_unit%log,'(f12.6,2x,a,i9,1p,e18.6,2x,a,i5,2x,5l1)') &
          wallclock(), 'sent', mesg%id, self%task%time, 'to', rank, &
          self%task%is_set (bits%internal), &
          self%task%is_set (bits%boundary), &
          self%task%is_set (bits%virtual), &
          self%task%is_set (bits%external), &
          self%task%is_set (bits%swap_request)
        flush (io_unit%log)
      end if
    end if
    nbor => nbor%next                                           ! next nbor
  end do
  !-----------------------------------------------------------------------------
  ! Add the message to mesg%sent_list, if it was sent, and clear the swap bit
  !-----------------------------------------------------------------------------
  if (flag) call mpi_mesg%sent (mesg)                           ! add to sent_list
  call self%task%clear (bits%swap_request)
  call trace_end (itimer)
END SUBROUTINE send_to_vnbors

END MODULE link_mod
