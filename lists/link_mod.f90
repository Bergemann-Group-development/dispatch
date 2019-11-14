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
  USE timer_mod
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
    type(lock_t):: lock, qlock
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
    procedure, nopass:: init_verbose
    procedure:: add_nbor_by_quality
    procedure:: add_nbor_by_rank
    procedure:: add_nbor_link_by_rank
    procedure:: remove_nbor
    procedure:: remove_from_nbors
    procedure:: remove_nbor_list
    procedure:: remove_nbor_list2
    procedure:: copy_nbor_list
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    !procedure:: update_nbors
    procedure:: log_nbors
    procedure:: print_nbors
    procedure:: check_level_sort
    procedure, nopass:: garbage_collect
    procedure, nopass:: garbage_remove
    procedure:: delete
    procedure:: sort_nbors_by_level
    procedure:: nbor_info
    procedure:: info
  end type
  type(link_t), target, public:: garbage
  integer:: garbage_n=0, verbose=0
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
    call self%lock%init ('nbor')
    call self%qlock%init ('queue')
    if (associated(self%task)) then
      self%lock%id  = self%task%id
      self%qlock%id = self%task%id
    end if
  end if
END SUBROUTINE init

!===============================================================================
!> Set the verbosity level (called from task_mesg_mod -- duplicates the setting
!> of verbose in the task_mesg_params namelist.
!===============================================================================
SUBROUTINE init_verbose (verb)
  integer:: verb
  !.............................................................................
  !$omp atomic write
  verbose = verb
END SUBROUTINE init_verbose

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
  !$omp atomic update
  this%task%n_needed = this%task%n_needed + 1   ! because linked task derefined
  if (verbose > 1) &
    write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
    this%task%id, 'needed by', this%task%n_needed, 'add_by_q'
  flush (io_unit%mpi)
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
!> Add a task link to an nbor list, setting also the nbor relations, based on
!> the ranks of the link%task and the self%task (which must be the "owner" of
!> the nbor_list)
!===============================================================================
SUBROUTINE add_nbor_link_by_rank (self, nbor_list, link)
  class(link_t):: self
  class(link_t), pointer:: nbor_list, link
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%add_nbor_link_by_rank')
  allocate (nbor)
  nbor%link => link
  nbor%task => link%task
  call self%add_nbor_by_rank (nbor_list, nbor)
  call self%task%nbor_relations (nbor%task, nbor%needed, nbor%needs_me, &
    nbor%download)
  call trace%end ()
END SUBROUTINE add_nbor_link_by_rank

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
  !$omp atomic update
  this%task%n_needed = this%task%n_needed + 1   ! because linked task derefined
  if (verbose > 1) &
    write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
    this%task%id, 'needed by', this%task%n_needed, 'add_by_rank'
  flush (io_unit%mpi)
  if (io%verbose>2) &
  write (io_unit%log,*) '  adding task', this%task%id, &
    this%task%is_set(bits%boundary), this%task%is_set(bits%virtual), &
    ' to nbor list of', self%task%id
  nullify (prev)
  next => nbors
  !-------------------------------------------------------------------------
  ! Check that the task is not self
  !-------------------------------------------------------------------------
  if (this%task%id == self%task%id) then
    goto 9
  end if
  !-------------------------------------------------------------------------
  ! Check that the new task is not already in the list
  !-------------------------------------------------------------------------
  do while (associated(next))
    if (next%task%id == this%task%id) then
      goto 9
    end if
    next => next%next
  end do
  !-------------------------------------------------------------------------
  ! If not, find the right place to add it
  !-------------------------------------------------------------------------
  next => nbors
  do while (associated(next))
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
  call trace_end
END SUBROUTINE add_nbor_by_rank


!===============================================================================
!> Remove this from the nbor list
!===============================================================================
SUBROUTINE remove_nbor (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: nbor, prev, next
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('link_mod::remove_nbor', 1, itimer=itimer)
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
  call trace%end (itimer)
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
  class(link_t), pointer:: nbors, nbor, next, link2, nbor2, next2, prev
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%begin ('link_t%remove_nbor_list', itimer=itimer)
  if (verbose > 1) &
    write (io_unit%mpi,'(f12.6,i4,2x,a)') wallclock(), omp%thread, 'remove_nbor_list'
    nbor => nbors
  do while (associated(nbor))
    next => nbor%next
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed - 1
    if (verbose > 1) then
      write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
        nbor%task%id, 'needed by', nbor%task%n_needed, 'remove_nb_list'
      flush (io_unit%mpi)
    end if
    deallocate (nbor)
    nbor => next
  end do
  nullify (nbors)
  if (timer%detailed) &
    call trace%end (itimer)
END SUBROUTINE remove_nbor_list

!===============================================================================
!> Remove, deallocate, and nullify an existing nbor list
!===============================================================================
SUBROUTINE remove_nbor_list2 (self, nbors)
  class(link_t):: self
  class(link_t), pointer:: nbors, nbor, next, link2, nbor2, next2, prev
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%remove_nbor_list2', itimer=itimer)
  if (verbose > 1) &
    write (io_unit%mpi,'(f12.6,i4,2x,a)') wallclock(), omp%thread, 'remove_nbor_list'
    nbor => nbors
  do while (associated(nbor))
    next => nbor%next
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed - 1
    if (verbose > 1) then
      write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
        nbor%task%id, 'needed by', nbor%task%n_needed, 'remove_nb_list'
      flush (io_unit%mpi)
    end if
    !---------------------------------------------------------------------------
    ! Also remove self from the nbor list of nbor
    !---------------------------------------------------------------------------
    nullify (prev)
    nbor2 => nbor%link%nbor
    do while (associated(nbor2))
      if (nbor2%task%id == self%task%id) then
        if (associated(prev)) then
          prev%next => nbor2%next
        else
          nbor%link%nbor => nbor2%next
        end if
      end if
      prev => nbor2
      nbor2 => nbor2%next
    end do
    deallocate (nbor)
    nbor => next
  end do
  nullify (nbors)
  call trace%end (itimer)
END SUBROUTINE remove_nbor_list2

!===============================================================================
!> Create a temporary nbor list
!===============================================================================
SUBROUTINE copy_nbor_list (self, old_head, new_head)
  class(link_t):: self
  class(link_t), pointer:: old_head, new_head
  !.............................................................................
  class(link_t), pointer:: nbor, head, tail, new
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%copy_nbors', itimer=itimer)
  nullify (head, tail)
  nbor => old_head
  do while (associated(nbor))
    allocate (new, source=nbor)
    new%task => nbor%task
    new%link => nbor%link
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed + 1
    if (verbose > 1) &
      write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
      nbor%task%id, 'needed by', nbor%task%n_needed, 'copy_nb_list'
    flush (io_unit%mpi)
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
  call trace%end (itimer)
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
!> Print a list of nbors to io_unit%mpi
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
!> Collect garbage, and check if there is garbage to remove.  This needs to be
!> done under lock protection, since many threads may be utilizing the GC at the
!> same time.
!===============================================================================
SUBROUTINE garbage_collect (link)
  class(link_t), pointer:: link, nbor
  !.............................................................................
  call trace%begin ('link_t%garbage_collect')
  flush (io_unit%mpi)
  !-----------------------------------------------------------------------------
  ! Immediately decrement the n_needed of the task, corresponding to the pending
  ! deletion, and also reduce the n_needed of nbors and nbors_by_level, and remove
  ! those list.  The task itself will be removed when the n_needed count drops
  ! to zero.
  !-----------------------------------------------------------------------------
  !$omp atomic update
  link%task%n_needed = link%task%n_needed-1
  if (verbose > 1) &
    write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
    link%task%id, 'needed by', link%task%n_needed, 'garbage_collect'
  call link%remove_nbor_list (link%nbor)
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
  call trace%end()
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
        write (io_unit%log,'(f12.6,2x,a,i6,i4,l4)') &
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
  call trace%end()
END SUBROUTINE garbage_remove

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
  call link%task%log ('delete')
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
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%begin ('link_t%sort_nbors_by_level', itimer=itimer)
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
  if (timer%detailed) &
    call trace%end (itimer)
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
    new%task => nbor%task
    new%link => nbor%link
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed + 1
    if (verbose > 1) &
      write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
      nbor%task%id, 'needed by', nbor%task%n_needed, 'sort_by_level'
    flush (io_unit%mpi)
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

!===============================================================================
!> Info for one nbor
!===============================================================================
SUBROUTINE nbor_info (self, task)
  class(link_t):: self
  class(task_t):: task
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(3x,a,i6,2i4,2x,3l1,2x,2l1,3x,3f9.3)') &
    'pa_t%nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
    self%task%id, self%task%rank, self%task%level, self%needed, self%needs_me, &
    self%download, self%task%is_set(bits%boundary), &
    self%task%is_set(bits%virtual), self%task%distance(task)/ &
    (0.5d0*(task%size+self%task%size))
END SUBROUTINE nbor_info

!===============================================================================
!> Info for all nbors
!===============================================================================
SUBROUTINE info (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call self%task%task_info()
  nbor => self%nbor
  do while (associated(nbor))
    call nbor%nbor_info (self%task)
    nbor => nbor%next
  end do
END SUBROUTINE info

END MODULE link_mod
