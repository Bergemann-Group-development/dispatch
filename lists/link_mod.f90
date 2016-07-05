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
  USE patch_mod
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
  contains
    procedure:: add_nbor_by_quality
    procedure:: add_nbor_by_rank
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    procedure:: update_nbors
    procedure:: send_to_vnbors
    procedure:: print_nbors
  end type
  public task2patch
  logical, save:: debug=.false.
CONTAINS

!> =============================================================================
!> Initialize the nbors chain of a new patch (self), by
!> 1) adding the parent
!> 2) searching the parent nbors, adding patches that overlap
!> Each time it adds an nbor, it should also add itself to that nbor's nbor list
!> =============================================================================
SUBROUTINE make_new_nbors (self, selfp, parent)
  class(link_t):: self
  class(link_t), pointer:: selfp, parent, nbor
  class(patch_t), pointer:: p1, p2
  !> ...........................................................................
  call trace_begin('list_t%make_new_nbors', 1)
  p1 => task2patch (self%task)
  call self%make_new_nbor (selfp, parent)       ! make me a neighbor of my parent
  nbor => parent%nbor                           ! start on parent's nbors chain
  do while (associated(nbor))
    p2 => task2patch (nbor%task)
    if (.not.associated(self%task,nbor%task) &
      .and. p1%overlaps (p2)) &   ! if my task overlaps with that task
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
!> Add a link node into the nbors list of self, in increasing quality order.
!> self must be a link whose 'nbors' points to a chain of nbors
!> =============================================================================
SUBROUTINE add_nbor_by_rank (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_nbor_by_rank')
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
        call trace_end
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
FUNCTION task2patch (task) RESULT (patch)
  class(task_t), pointer:: task
  class(patch_t), pointer:: patch
  select type (task)
  class is (patch_t)
  patch => task
  class default
  nullify(patch)
  stop
  end select
END FUNCTION

!> =============================================================================
!> Update the nbor list of self, adding self to the each nbor's nbor list
!> =============================================================================
SUBROUTINE update_nbors (self)
  class(link_t):: self
  class(link_t), pointer:: nbor, new
  !.............................................................................
  call trace_begin('list_mod::init_all_nbors')
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
!> Print a list of nbors
!> =============================================================================
SUBROUTINE print_nbors (self, label)
  class (link_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 2) return
  if (present(label)) then
    if (trim(label) /= '') print'(a)', label
  end if
  print'(a,i5,2x,a,f14.8,2x,a,l3,2x,a,i3,2x,a,$)', &
    'task', self%task%id, 'time', self%task%time, &
    'ready', self%task%is_set(bits%ready), 'level', self%task%level, 'nbors:'
  if (io%verbose>4) print*,''
  link => self%nbor
  do while (associated(link))
    if (io%verbose>4) then
      print '(i7,f14.8,l3)', &
        link%task%id, link%task%time, link%task%is_set(bits%ready)
    else
      print '(i4,l2,$)', link%task%id, link%task%is_set(bits%ready)
    end if
    link => link%next
  end do
  print*,''
END SUBROUTINE print_nbors

!> =============================================================================
!> Send a package with updated time slice info to neighbor ranks.
!>
!> 1) assemble patch data into a message
!> 3) send the message to the nbors rank (only once per rank)
!> 4) add the message to a sent_list
!> 5) check the sent list for completed sends
!> 6) remove the message from the sent_list and disassemble the message
!> =============================================================================
SUBROUTINE send_to_vnbors (self, mesg)
  class(link_t):: self
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  class(mesg_t), pointer:: mesg, next
  logical:: flag
  integer:: ierr, rank
  !.............................................................................
  call trace_begin('list_mod::send_to_vnbors')
  !-----------------------------------------------------------------------------
  ! Send to all unique ranks in the nbor list (which is sorted by rank)
  !-----------------------------------------------------------------------------
  rank = -1                                                     ! previous rank
  nbor => self%nbor                                             ! first nbor
  do while (associated(nbor))                                   ! until end
    task => nbor%task
    if (task%rank/=mpi%rank .and. task%rank/=rank) then         ! new rank?
      rank = nbor%task%rank                                     ! current rank
      if (debug) then
        write (io_unit%log,*) ' send mesg', mesg%id, ' to', rank; flush(io_unit%log)
        flush (io%log_unit)
      end if
      call mesg%send (rank)                                     ! send it
    end if
    nbor => nbor%next                                           ! next nbor
  end do
  !-----------------------------------------------------------------------------
  ! Add the message to a list of sent messages, and check which previously sent
  ! messages that may now be deallocated.  The test is fast, so we can leave the
  ! list locked while updating.
  !-----------------------------------------------------------------------------
  call mpi_mesg%sent_list%lock%set ('send')
  if (debug) write(io%log_unit,*) 'sent_list%add mesg', mesg%id
  if (debug) flush(io%log_unit)
  call mpi_mesg%sent_list%add (mesg)                            ! add to sent_list
  mesg => mpi_mesg%sent_list%head                               ! top of list
  do while (associated(mesg))
    next => mesg%next                                           ! in case removed
    call mesg%test_all (flag)                                   ! done with?
    if (flag) then
      if (debug) write(io%log_unit,*) 'sent_list remove mesg', mesg%id
      if (debug) flush(io%log_unit)
!print *, mpi%rank, 'sent_list remove mesg', mesg%id; flush(6)
      call mpi_mesg%sent_list%remove (mesg)                     ! yes, remove
    else
!print *, mpi%rank, 'send_to_vnbors mesg', mesg%id, 'not removed'; flush(6)
    end if
    mesg => next                                                ! true next
  end do
  call mpi_mesg%sent_list%lock%unset ('send')
  call trace_end
END SUBROUTINE send_to_vnbors

END MODULE link_mod
