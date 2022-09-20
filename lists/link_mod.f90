!*******************************************************************************
!> Module with list handling for generic class task_t objects.
!>
!> NOTE: The advantages of using the same links to organize a ready queue and
!> an active queue is that the process is faster and simpler when no new links
!> have to be created
!*******************************************************************************
MODULE link_mod
  USE io_mod
  USE os_mod
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
  USE shared_mod
  implicit none
  private
  type, public:: nbor_t
    class(nbor_t), pointer :: next => null()          ! next link in nbor list
    class(link_t), pointer :: link => null()          ! pointer to task link
    class(task_t), pointer :: task => null()          ! pointer to task
    !---------------------------------------------------------------------------
    ! Flags making it possible to distinguish between nbors needed for link%task,
    ! and nbors that need link%task.  These flags are only used in nbor lists,
    ! so when parsing the list starting with link%nbor, one uses nbor%needed
    ! in check_ready(), and nbor%needs_me in check_nbors().
    !---------------------------------------------------------------------------
    logical:: needed   = .true.     ! marks nbors included in check_nbors()
    logical:: needs_me = .true.     ! marks nbors included in check_ready()
    logical:: download = .true.     ! marks nbors included in download_link()
    real:: u_max = 0.0
    integer:: check_step = 0        ! step when target cleared check_ready
    real(8):: distance(3) = 0d0
    class(task_t), pointer:: target => null()
    class(task_t), pointer:: source => null()
  contains
    procedure:: init => init_nbor
    procedure:: update => update_nbor
    procedure:: nbor_info
  end type
  type, public:: link_t
    class(task_t), pointer :: task => null()          ! task
    class(nbor_t), pointer :: nbor => null()          ! next nbor; wouldn't `nbor_head` be a better name?
    class(nbor_t), pointer :: nbors_by_level => null()! FIXME: obsolete
    class(link_t), pointer :: link => null()          ! pointer from nbor link to task list link (can be useful)
    class(link_t), pointer :: next => null()          ! next link in task list
    class(link_t), pointer :: prev => null()          ! previous link in task list
    class(link_t), pointer :: parent => null()        ! parent task link
    class(link_t), pointer :: next_time => null()     ! next link in time order
    class(link_t), pointer :: prev_time => null()     ! prev link in time order
    class(link_t), pointer :: next_active => null()   ! next link to an active process, in time order
    !class(nbor_t), pointer :: nbors_by_level => null()! decreasing level sorted list
    type(lock_t):: lock
    logical:: needs_nbors = .true.  ! flag indicating if a link is in need of adding nbors
    logical:: initialized = .false.
    logical:: check = .true.
    integer:: verbose = 0
  contains
    procedure:: init
    procedure:: reset_link
    procedure:: update
    procedure:: download
    procedure, nopass:: init_verbose
    procedure, nopass:: alloc_link
    procedure, nopass:: alloc_nbor
    procedure:: remove_nbor
    procedure:: remove_from_nbors
    procedure:: remove_nbor_list
    procedure:: remove_nbor_list2
    procedure:: copy_nbor_list
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    procedure:: check_nbor_relations
    procedure:: log_nbors
    procedure:: print_nbors
    procedure:: check_level_sort
    procedure, nopass:: garbage_collect
    procedure, nopass:: garbage_remove
    procedure:: delete
    procedure:: sort_nbors_by_level
    procedure:: reset_status
    procedure:: needs_to_add
    procedure:: add_nbor
    procedure:: add_to_nbors
    procedure:: add_nbor_by_level
    procedure:: nbor_inconsistency
    procedure:: info
    procedure:: request_support
    procedure:: has_support
    procedure:: check_status
  end type
  type(link_t), target, public:: garbage
  integer:: garbage_n=0, verbose=0
  logical, save:: debug=.false.
CONTAINS

!===============================================================================
!> For now, we set the distance from known task positions.
!> FIXME: The nbor%distance should ultimately be the only position information
!===============================================================================
SUBROUTINE init_nbor (self)
  class(nbor_t):: self
  !.............................................................................
  call trace%begin ('nbor_t%init')
  call trace%end ()
END SUBROUTINE init_nbor

!===============================================================================
!> Update the part of a halo (possibly abstract) related to an nbor
!===============================================================================
SUBROUTINE update_nbor (self, all, only)
  class(nbor_t):: self
  logical, optional:: all
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('nbor_t%update')
  call trace%end ()
END SUBROUTINE update_nbor

!===============================================================================
!> If link is not already initialized, initialize its lock.  If an init is at
!> all made (not done for nbor links), it is done by the thread that allocated
!> the link, so there should be no need for a crtical region.
!===============================================================================
SUBROUTINE init (self)
  class(link_t):: self
  !.............................................................................
  if (.not. self%initialized) then
    self%initialized = .true.
    call self%lock%init ('link')
    if (associated(self%task)) then
      self%lock%id  = self%task%id
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Link locking is mandatory
  !-----------------------------------------------------------------------------
  self%lock%on = .true.
END SUBROUTINE init

!===============================================================================
!> Reset to essentially pre-init conditions
!===============================================================================
SUBROUTINE reset_link (self)
  class(link_t):: self
  !.............................................................................
  if (self%initialized) then
    self%initialized = .false.
    call self%lock%destroy ()
  end if
  nullify (self%nbor, self%next, self%prev)
END SUBROUTINE reset_link

!===============================================================================
!> Update the guard zone, by using task%dnload followed by task%bndries
!===============================================================================
SUBROUTINE update (self, all_tgt, all_src, only)
  implicit none
  class(link_t)           :: self
  logical, optional       :: all_tgt, all_src
  integer, optional       :: only
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%update')
  call self%task%dnload
  call self%task%bndries
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Stub
!===============================================================================
SUBROUTINE download (self, nbors, only, all_tgt, all_src)
  class(link_t)          :: self
  class(nbor_t), pointer :: nbors
  integer,      optional :: only
  logical,      optional :: all_tgt, all_src
  !-----------------------------------------------------------------------------
END SUBROUTINE download

!===============================================================================
!> Set the verbosity level (called from task_mesg_mod -- duplicates the setting
!> of verbose in the task_mesg_params namelist).
!===============================================================================
SUBROUTINE init_verbose (verb)
  integer:: verb
  !.............................................................................
  !$omp atomic write
  verbose = verb
  if (verb > 0) &
  write (io_unit%log,*) 'link_t%init_verbose: verbose =', verb
END SUBROUTINE init_verbose

!===============================================================================
!> Initialize the nbors chain of a new task (self), by
!> 1) adding the parent
!> 2) searching the parent nbors, adding tasks that overlap
!> Each time it adds an nbor, it should also add itself to that nbor's nbor list
!===============================================================================
SUBROUTINE make_new_nbors (self, selfp, parent)
  class(link_t):: self
  class(link_t), pointer:: selfp, parent
  class(nbor_t), pointer:: nbor
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
      call self%make_new_nbor (selfp, nbor%link)! make a neighbor (and make me one)
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
  write (stderr,*) 'list_t%make_new_nbor is deprecated - interferes with other tasks!'
  call self%add_nbor (that)
  call that%add_nbor (selfp)
  call trace_end
END SUBROUTINE make_new_nbor

!===============================================================================
!> Add a link node into the nbors list of self, by default in increasing level
!> order and optionally in quality or rank order.  Optionally return the nbor
!> link, and/or a flag indicating whether or the task was added (is rejected
!> if it is already on the list)
!===============================================================================
SUBROUTINE add_nbor (self, link, nbors, nbor, added, by_quality, by_rank, &
                     rejected)
  class(link_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer, optional:: nbors
  class(nbor_t), pointer, optional:: nbor
  logical,                optional:: added, rejected
  logical,                optional:: by_quality, by_rank
  !.............................................................................
  logical:: ok
  class(nbor_t), pointer:: next, prev, nbor_l
  character(8):: kind
  character(20):: label
  integer, save:: itimer=0
  !.............................................................................
  if (present(rejected)) rejected = .false.
  if (present(by_quality)) then
    call trace%begin ('list_t%add_nbor_by_quality', itimer=itimer)
    kind = 'quality'
  else if (present(by_rank)) then
    call trace%begin ('list_t%add_nbor_by_rank', itimer=itimer)
    kind = 'rank'
  else
    call trace%begin ('list_t%add_nbor_by_level', itimer=itimer)
    kind = 'level'
  end if
  write (label,'(a8,i6)') 'add_nbor', link%task%id
  call self%task%log (label, 4)
  !-----------------------------------------------------------------------------
  ! Use a class method to overload the allocation
  !-----------------------------------------------------------------------------
  nullify (nbor_l)
  call self%alloc_nbor (link, nbor_l, 'add_nbor')
  !-----------------------------------------------------------------------------
  ! Consider this carefully!
  ! From Andrius: please consult me before changing this.
  !-----------------------------------------------------------------------------
  call self%check_nbor_relations(nbor_l%task, nbor_l%needed, nbor_l%needs_me, & 
                                 nbor_l%download)
  !.............................................................................
  if (link%task%id == io%id_debug) &
    write (stdlog,'(a,2(i6,":",i2.2),1x,3l1)') &
      'link_t%add_nbor: id,nbid,needed,needs_me,download=', &
      link%task%id, link%task%rank, nbor_l%task%id, nbor_l%task%rank, &
      nbor_l%needed, nbor_l%needs_me, nbor_l%download
  !.............................................................................
  if (.not. nbor_l%needed   .and. &
      .not. nbor_l%needs_me .and. &
      .not. nbor_l%download) then
    call os%mem_register (storage_size(nbor_l), -1, 'add_nbor')
    deallocate (nbor_l)
    nullify (nbor_l)
    if (present(added)) &
      added=.false.
    if (present(nbor)) &
      nullify (nbor)
    call trace%end (itimer)
    return
  end if
  !.............................................................................
  if (verbose > 1) &
    write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
    link%task%id, 'needed by', link%task%n_needed, 'add_by_'//trim(kind)
  !.............................................................................
  ok = .true.
  nullify (prev)
  if (present(nbors)) then
    next => nbors
  else
    next => self%nbor
  end if
  do while (associated(next))
    if (next%task%id == link%task%id) then
      !-------------------------------------------------------------------------
      ! Return if link task is already in the nbor list
      !-------------------------------------------------------------------------
      ok = .false.
      exit
    else if (present(by_quality)) then
      if (next%task%quality > link%task%quality) then
        !-----------------------------------------------------------------------
        ! Found a task we should be ahead of, so insert here
        !-----------------------------------------------------------------------
        nbor_l%next => next
        exit
      end if
    else if (present(by_rank)) then
      if (next%task%rank > link%task%rank) then
        !-----------------------------------------------------------------------
        ! Found a task we should be ahead of, so insert here
        !-----------------------------------------------------------------------
        nbor_l%next => next
        exit
      end if
    else
      if (next%task%level > link%task%level) then
        !-----------------------------------------------------------------------
        ! Found a task we should be ahead of, so insert here
        !-----------------------------------------------------------------------
        nbor_l%next => next
        exit
      end if
    end if
    prev => next
    next => next%next
  end do
  if (ok) then
    if (present(added)) added = .true.
    !---------------------------------------------------------------------------
    ! Found a place to insert, it link it in
    !---------------------------------------------------------------------------
    if (associated(prev)) then
      prev%next => nbor_l
    else
      !-------------------------------------------------------------------------
      ! If prev is not associated, this is the first nbor
      !-------------------------------------------------------------------------
      if (present(nbors)) then
        nbors => nbor_l
      else
        self%nbor => nbor_l
      end if
    end if
    !$omp atomic update
    link%task%n_nbors = link%task%n_nbors + 1
    !$omp atomic update
    link%task%n_needed = link%task%n_needed + 1
    if (verbose > 1) &
      write (stderr,'(2(a,i6,":",i2.2))') &
        ' link_t%add_nbor: added', link%task%id, link%task%level,&
        ' to the nbor list of task', self%task%id, self%task%level
  else
    if (present(added)) added = .false.
    if (present(rejected)) rejected = .true.
    write (stderr,'(2(a,i6,":",i2.2))') &
      'link_t%add_nbor WARNING: trying to add existing nbor', &
      next%task%id, next%task%level,&
      ' to the nbor list of task', &
      self%task%id, self%task%level
  end if
  if (present(nbor)) then
    nbor => nbor_l
  end if
  call trace%end (itimer)
END SUBROUTINE add_nbor

!===============================================================================
!> Add linked task to the nbor lists of tasks on it's owne nbor list
!===============================================================================
SUBROUTINE add_to_nbors (self, link)
  class(link_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  nbor => link%nbor
  do while (associated(nbor))
    call nbor%link%add_nbor (link)
    nbor => nbor%next
  end do
END SUBROUTINE add_to_nbors

!===============================================================================
!> Add an nbor node into the nbors list of self, in increasing rank order.
!> self must be a link whose 'nbors' points to a chain of nbors, and 'nbor'
!> must be an nbor with both a link and task pointer set.
!===============================================================================
SUBROUTINE add_nbor_by_level (self, nbors, nbor, added)
  class(link_t):: self
  class(nbor_t), pointer:: nbors, nbor
  class(nbor_t), pointer:: next, prev
  logical, optional:: added
  integer:: unit
  !.............................................................................
  call trace_begin ('link_t%add_nbor_by_level', 1)
  if (io%verbose > 2 .or. self%task%id == io%id_debug) then
    unit = merge (stdout, stdlog, self%task%id==io%id_debug)
    write (unit,'(2(a,i6,i3,3x),a,2l1)') &
      'link_t%add_nbor_by_level:  adding task', nbor%task%id, nbor%task%level, &
      'to nbor list of', self%task%id, self%task%level, &
      'BV=', nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual)
  end if
  if (present(added)) added = .false.
  nullify (prev)
  next => nbors
  !-----------------------------------------------------------------------------
  ! Check that the task is not self
  !-----------------------------------------------------------------------------
  if (nbor%task%id == self%task%id) then
    write (stderr,*) self%task%id, 'ERROR: trying to add self to nbor list of', &
      nbor%task%id
    goto 9
  end if
  !-----------------------------------------------------------------------------
  ! Consider this carefully!
  ! From Andrius: please consult me before changing this.
  !-----------------------------------------------------------------------------
  call self%check_nbor_relations(nbor%task, nbor%needed, nbor%needs_me, & 
                                 nbor%download)
  if (.not. nbor%needed   .and. &
      .not. nbor%needs_me .and. &
      .not. nbor%download) then
    !call os%mem_register (storage_size(nbor), -1, 'add_nbor')
    deallocate (nbor)
    nullify (nbor)
    if (present(added)) &
      added = .false.
    nullify (nbor)
    call trace%end ()
    return
  end if
  !-----------------------------------------------------------------------------
  ! Add it before the next higher levelm unless it is found before then
  !-----------------------------------------------------------------------------
  nullify (prev)
  next => nbors
  do while (associated(next))
    if (next%task%id == nbor%task%id) then
      if (verbose > 1) then
        write (stdlog,'(f12.6,a,2i6)') wallclock(), pad_trace(1:indent)// &
          'INFO: nbor already on nbor list, self,nbor=', &
          self%task%id, nbor%task%id
      end if
      goto 9
    else if (next%task%level > nbor%task%level) then
      !-------------------------------------------------------------------------
      ! Found a task we should be ahead of, so insert here
      !-------------------------------------------------------------------------
      nbor%next => next
      exit
    end if
    prev => next
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! This task has the highest level, so append it's nbor link
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    prev%next => nbor
  else
    nbors => nbor
  end if
  nullify (nbor%next)
  if (present(added)) added = .true.
  self%task%n_nbors = self%task%n_nbors + 1
  !$omp atomic update
  nbor%task%n_needed = nbor%task%n_needed + 1   ! because linked task derefined
9 continue
  call trace%end ()
END SUBROUTINE add_nbor_by_level

!===============================================================================
!> Overloaded procedure
!===============================================================================
SUBROUTINE alloc_nbor (link, nbor, label)
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  allocate (nbor_t:: nbor)
  call os%mem_register (storage_size(nbor), 1, shared%option(label, 'alloc_nbor'))
  call nbor%init ()
  nbor%link => link
  nbor%task => link%task
END SUBROUTINE alloc_nbor

!===============================================================================
!> Overloaded procedure
!===============================================================================
SUBROUTINE alloc_link (task, link)
  class(task_t), target:: task
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  if (.not. associated(link)) then
    allocate (link)
    call os%mem_register (storage_size(link), 1, 'alloc_link')
    link%task => task
  end if
END SUBROUTINE alloc_link

!===============================================================================
!> Ensure "reflexive" nbor relations between self%task and nbtask, while
!> allowing only self%ask to change status bits
!===============================================================================
SUBROUTINE check_nbor_relations (self, nbtask, needed, needs_me, download)
  class(link_t):: self
  class(task_t), pointer:: nbtask
  logical:: needed, needs_me, download
  !.............................................................................
  logical:: my_needed, my_needs_me, my_download
  logical:: nb_needed, nb_needs_me, nb_download
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%check_nbor_relations')
  call self%task%nbor_relations(nbtask,my_needed,my_needs_me,my_download)
  call nbtask%nbor_relations(self%task,nb_needed,nb_needs_me,nb_download)
  if (my_needed /= nb_needs_me .or. my_needs_me /= nb_needed) then
    print *, 'nbor_relations inconstistent for tasks ', self%task%id, nbtask%id
    print *, my_needed,my_needs_me,my_download
    print *, nb_needed,nb_needs_me,nb_download
  endif
  needed   = (my_needed   .or. nb_needs_me)
  needs_me = (my_needs_me .or. nb_needed)
  download = my_download
  call trace%end ()
END SUBROUTINE check_nbor_relations

!===============================================================================
!> Check the status of a link and all nbor links
!===============================================================================
SUBROUTINE check_level_sort (self)
  class(link_t):: self
  class(nbor_t), pointer:: nbor
  integer:: level
  !.............................................................................
  nbor => self%nbor
  level = -1
  do while (associated(nbor))
    if (nbor%task%level < level) then
! print *, self%task%id, 'WARNING: nbors not sorted by level'
       nbor => self%nbor
       do while (associated(nbor))
! print *, 'nbor, level =', nbor%task%id, nbor%task%level
         nbor => nbor%next
       end do
       return
    end if
    level = nbor%task%level
    nbor => nbor%next
  end do
END SUBROUTINE check_level_sort

!===============================================================================
!> Remove the this task from links nbor list 
!===============================================================================
SUBROUTINE remove_nbor (self, this, nbor_list)
  class(link_t):: self
  class(task_t), pointer:: this
  class(nbor_t), pointer:: nbor_list
  !.............................................................................
  class(nbor_t), pointer:: nbor, prev, next
  logical:: failed
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%remove_nbor', 1, itimer=itimer)
  !.............................................................................
  if (verbose > 0) &
    write (io_unit%log,*) 'link_t%remove_nbor: ', this%id, &
      ' from nbor list of', self%task%id
  !.............................................................................
  call                                             self%lock%set ('remove_nbor')
  nullify (prev)
  failed = .true.
  nbor => nbor_list
  do while (associated(nbor))
    next => nbor%next
    if (associated(nbor%task,this)) then
      if (associated(prev)) then
        prev%next => nbor%next
      else
        nbor_list => nbor%next
      end if
      self%task%n_nbors = self%task%n_nbors-1
      call os%mem_register (storage_size(nbor), -1, 'remove_nbor')
      deallocate (nbor)
      !$omp atomic update
      this%n_needed = this%n_needed - 1
      nullify (nbor)
      failed = .false.
      exit
    end if
    prev => nbor
    nbor => next
  end do
  call                                           self%lock%unset ('remove_nbor')
  if (failed) then
    write (io_unit%log,*) &
      'link_t%remove_nbor: WARNING, could not find task', this%id
  end if
  call trace%end (itimer)
END SUBROUTINE remove_nbor

!===============================================================================
!> Remove this task from the nbor lists of the nbors of self
!===============================================================================
SUBROUTINE remove_from_nbors (self, this)
  class(link_t):: self
  class(link_t), pointer:: this
  class(nbor_t), pointer:: nbor
  !.............................................................................
  call trace_begin ('link_t%remove_from_nbors', 2)
  call self%lock%set ('remove_from_nbors')
  !write (stderr,*) 'link_t%remove_from:nbors is DEPRECATED'
  nbor => self%nbor
  do while (associated(nbor))
    call nbor%link%remove_nbor (this%task, nbor%link%nbor)
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
  class(nbor_t), pointer:: nbors, nbor, next
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed) &
    call trace%begin ('link_t%remove_nbor_list', itimer=itimer)
  if (verbose > 1) &
    write (io_unit%mpi,'(f12.6,i4,2x,a)') wallclock(), omp%thread, 'remove_nbor_list'
  if (associated(nbors)) then
    nbor => nbors
    do while (associated(nbor))
      next => nbor%next
      !$omp atomic update
      nbor%task%n_needed = nbor%task%n_needed - 1
      if (verbose > 2) then
        write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
          nbor%task%id, 'needed by', nbor%task%n_needed, 'remove_nb_list'
        flush (io_unit%mpi)
      end if
      call os%mem_register (storage_size(nbor), -1, 'remove_nbor_list')
      deallocate (nbor)
      nbor => next
    end do
    nullify (nbors)
  end if
  if (timer%detailed) &
    call trace%end (itimer)
END SUBROUTINE remove_nbor_list

!===============================================================================
!> Remove, deallocate, and nullify an existing nbor list
!===============================================================================
SUBROUTINE remove_nbor_list2 (self, nbors)
  class(link_t):: self
  class(nbor_t), pointer:: nbors, nbor, next, nbor2, next2, prev
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
    if (verbose > 2) then
      write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
        nbor%task%id, 'needed by', nbor%task%n_needed, 'remove_nb_list2'
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
    call os%mem_register (storage_size(nbor), -1, 'remove_nbor_list2')
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
  class(nbor_t), pointer:: old_head, new_head
  !.............................................................................
  class(nbor_t), pointer:: nbor, head, tail, new
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%copy_nbors', itimer=itimer)
  nullify (head, tail)
  nbor => old_head
  do while (associated(nbor))
    allocate (new, source=nbor)
    call os%mem_register (storage_size(new), 1, 'remove_nbor_list2')
    new%task => nbor%task
    new%link => nbor%link
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed + 1
    if (verbose > 2) then
      write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i5,2x,a)') wallclock(), omp%thread, &
      nbor%task%id, 'needed by', nbor%task%n_needed, 'copy_nb_list'
      flush (io_unit%mpi)
    end if
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
!> Print a list of nbors to the rank-thread log file
!===============================================================================
SUBROUTINE log_nbors (self, label)
  class (link_t):: self
  character(len=*), optional:: label
  class(nbor_t) , pointer:: nbor
  !.............................................................................
  if (io%verbose < 1) return
  if (present(label)) then
    if (trim(label) /= '') write (io_unit%log,'(a)') label
  end if
  write (io_unit%log,'(a,i5,2x,a,f14.8,2x,a,l3,2x,a,i3,2x,a,$)') &
    'task', self%task%id, 'time', self%task%time, &
    'ready', self%task%is_set(bits%ready), 'level', self%task%level, 'nbors:'
  if (io%verbose>4) write (io_unit%log,*) ''
  nbor => self%nbor
  do while (associated(nbor))
    if (io%verbose>4) then
      write (io_unit%log,'(i7,f14.8,l3)') &
        nbor%task%id, nbor%task%time, nbor%task%is_set(bits%ready)
    else
      write (io_unit%log,'(i4,l2,$)') nbor%task%id, nbor%task%is_set(bits%ready)
    end if
    nbor => nbor%next
  end do
  write (io_unit%log,*) ''
END SUBROUTINE log_nbors

!===============================================================================
!> Print a list of nbors to io_unit%mpi
!===============================================================================
SUBROUTINE print_nbors (self, label, verbose, unit, nbors)
  class (link_t):: self
  character(len=*), optional:: label
  integer, optional:: verbose, unit
  class(nbor_t), optional, pointer:: nbors
  !.............................................................................
  class(nbor_t), pointer:: nbor, nbors_L
  integer:: idmax, unit_L
  logical:: ahead
  character(len=16):: fmt
  !.............................................................................
  if (present(label)) write (io_unit%output,'(a,$)') label
  unit_L = shared%option (unit, stdout)
  write (unit_L, '(i6,l1," nbors:",$)') &
    self%task%id, self%task%is_set (bits%ready)
  if (present(nbors)) then
    nbors_L => nbors
  else
    nbors_L => self%nbor
  end if
  nbor => nbors_L
  idmax = 0
  idmax = 1
  do while (associated(nbor))
    idmax = max(idmax, nbor%task%id)
    nbor => nbor%next
  end do
  idmax = floor(log10(real(idmax)))+2
  write (fmt, '("(i",i1,",l1,$)")') idmax
  nbor => nbors_L
  do while (associated(nbor))
    !ahead = nbor%task%is_ahead_of (self%task)
    write (unit_L,fmt) nbor%task%id, nbor%download
    nbor => nbor%next
  end do
  write (unit_L,*) ''
END SUBROUTINE print_nbors

!===============================================================================
!> Collect garbage, and check if there is garbage to remove.  This needs to be
!> done under lock protection, since many threads may be utilizing the GC at the
!> same time.
!===============================================================================
SUBROUTINE garbage_collect (link)
  class(link_t), pointer:: link
  !class(nbor_t), pointer:: nbor
  !.............................................................................
  call trace%begin ('link_t%garbage_collect')
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
  !-----------------------------------------------------------------------------
  ! If no other task needs this task remove it immediately, else add to garbage
  !-----------------------------------------------------------------------------
  if (link%task%n_needed <= 0) then
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,2x,a,i4,2l4)') &
        wallclock(), 'garbage_collect deleted: id =', link%task%id, &
        link%task%n_needed, link%task%is_set(bits%virtual)
      write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
        wallclock(), omp%thread, 'removed', link%task%id, link%task%time
    end if
    !---------------------------------------------------------------------------
    ! At this point, the link task has been removed from the task list, and
    ! will never be the target of a download again, and hence will not be
    ! needing the tasks on its nbor list.  Their n_needed counts should thus
    ! be decremented, and the sorted nbors list can also be removed. The
    ! actual removal of the task is delayed until no other task needs it.
    !---------------------------------------------------------------------------
    call                                            garbage%lock%set ('garbage')
    call garbage_remove                       ! check consequences for nbors
    call                                          garbage%lock%unset ('garbage')
    call link%delete (link)
  else
    call                                            garbage%lock%set ('garbage')
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
    call                                          garbage%lock%unset ('garbage')
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
      !$omp atomic update
      timer%levelpop (link%task%level) = &        ! level population count
      timer%levelpop (link%task%level) - 1
      !$omp atomic update
      timer%deleted (link%task%level) = &        ! level population count
      timer%deleted (link%task%level) + 1
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
  class(link_t), pointer:: next
  character(2):: label
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%delete ')
  !-----------------------------------------------------------------------------
  call                                           link%lock%set ('link_t%delete')
  if (io%task_logging > 1) then
    write (label,'(2a1)') &
      merge('B','b',link%task%is_set(bits%boundary)), &
      merge('V','v',link%task%is_set(bits%virtual))
    call link%task%log ('delete task '//label, 2)
  end if
  do while (link%lock%level > 1)
    write (io_unit%output,*) 'delete_link unlocking lock level', link%lock%level
    call link%lock%unset
  end do
  !-----------------------------------------------------------------------------
  ! Remove / deallocate the nbor list before deleting the task
  !-----------------------------------------------------------------------------
  call link%remove_nbor_list (link%nbor)
  call link%task%dealloc()
  if (verbose > 1) then
    write (io_unit%output,'(f12.5,2x,a,2i6)') &
      wallclock(), 'delete: task, garbage_n =', link%task%id, garbage_n
    flush (io_unit%output)
  end if
  call os%mem_register (storage_size(link%task), -1, 'link_t%delete')
  deallocate (link%task)
  nullify (link%task)
  !-----------------------------------------------------------------------------
  call                                         link%lock%unset ('link_t%delete')
  call os%mem_register (storage_size(link), -1, 'link_t%delete')
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
  class(nbor_t), pointer:: old_head, new_head
  !.............................................................................
  class(nbor_t), pointer:: nbor, head, prev, find
  integer:: n_nbors
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (timer%detailed .or. verbose > 0) &
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
  if (timer%detailed .or. verbose > 0) &
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
    class(nbor_t), pointer:: new
    allocate (new, source=nbor)
    call os%mem_register (storage_size(new), 1, 'link_t%prepend')
    new%task => nbor%task
    new%link => nbor%link
    !$omp atomic update
    nbor%task%n_needed = nbor%task%n_needed + 1
    if (verbose > 1) then
      write (io_unit%mpi,'(f12.6,i4,i6,2x,a,i4,2x,a)') wallclock(), omp%thread, &
      nbor%task%id, 'needed by', nbor%task%n_needed, 'sort_by_level'
      flush (io_unit%mpi)
    end if
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
    class(nbor_t), pointer:: head, nbor
    !---------------------------------------------------------------------------
    write(stdlog,*) 'target: ', self%task%id
    write(stdlog,*) 'sorted nbors:'
    nbor => head
    do while (associated(nbor))
      write (stdlog,*) nbor%task%id, nbor%task%level
      nbor => nbor%next
    end do
  end subroutine nbor_print
  !-----------------------------------------------------------------------------
  ! Check new nbor list for size and consistency
  !-----------------------------------------------------------------------------
  subroutine nbor_check (head)
    class(nbor_t), pointer:: head, nbor
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
!> Reset the status bits and flags on a linked task and an nbor list
!===============================================================================
SUBROUTINE reset_status (self, nbor_list, check)
  class(link_t):: self
  class(link_t), pointer:: nblink
  class(nbor_t), pointer:: nbor_list
  logical, optional:: check
  !.............................................................................
  class(nbor_t), pointer:: nb
  class(task_t), pointer:: task
  type(task_t):: new
  logical:: flags(3,2)
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%reset_status')
  !-----------------------------------------------------------------------------
  ! Set new status, as a modified copy of self%status
  !-----------------------------------------------------------------------------
  task => self%task
  new%status = task%status
  call new%clear (bits%boundary+bits%virtual)
  if (task%rank == mpi%rank) then
    call new%set (bits%internal)
  else
    call new%set (bits%external)
  end if
  !-----------------------------------------------------------------------------
  ! Loop over nbor list and modify bits; B and V tasks have only that bit
  !-----------------------------------------------------------------------------
  nb => nbor_list
  do while (associated(nb))
    if (verbose > 1) then
      write (stdlog,'(a,2(i7,i5,2x))') &
      'reset_status: nbid,rank=', task%id, task%rank, nb%task%id, nb%task%rank
    end if
    if (task%rank == mpi%rank) then
      if (nb%task%rank /= mpi%rank) then
        call new%set (bits%boundary)
        call new%clear (bits%internal+bits%external)
      end if
    else
      if (nb%task%rank == mpi%rank) then
        call new%set (bits%virtual)
        call new%clear (bits%internal+bits%external)
      end if
    end if
    nb => nb%next
  end do
  if (task%is_set(bits%internal)) then
    write(stdout, *) task%id, 'has been set internal in link_mod'
  endif
  !-----------------------------------------------------------------------------
  ! If any change, call check_status() to take proper action
  !-----------------------------------------------------------------------------
  if (new%status /= task%status) then
    !$omp atomic write
    task%status = new%status
  end if
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  if (shared%option(check)) then
    !call self%check_nbor_relations 
  end if
  !.............................................................................
  if (verbose > 0) then
    write(stdlog,'(a,i6,i4,3x,4l1)') &
      'reset_status: id,rank,IBVE=', self%task%id, self%task%rank, &
      self%task%is_set(bits%internal), &
      self%task%is_set(bits%boundary), &
      self%task%is_set(bits%virtual ), &
      self%task%is_set(bits%external)
  end if
  call trace%end ()
END SUBROUTINE reset_status

!===============================================================================
!> First, check if there is relation between self and and nbor link.  If so,
!> check if the nbor task is already in the nbor list, and return false if so.
!===============================================================================
LOGICAL FUNCTION needs_to_add (self, nblink)
  class(link_t):: self
  class(link_t), pointer:: nblink
  !.............................................................................
  class(nbor_t), pointer:: nbor
  logical:: needed, needs_me, download, done
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%needs_to_add')
  call self%check_nbor_relations (nblink%link%task, needed, needs_me, download)
  if (needed .or. needs_me) then
    needs_to_add = .true.
    nbor => self%nbor
    do while (associated(nbor))
      if (associated (nbor%task, nblink%task)) then
        needs_to_add = .false.
        exit
      end if
      nbor => nbor%next
    end do
  end if
  call trace%end ()
END FUNCTION needs_to_add

!===============================================================================
!> Check that nbor relations are cross-symmetric
!===============================================================================
LOGICAL FUNCTION nbor_inconsistency (self)
  class(link_t):: self
  !.............................................................................
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%nbor_inconsistency')
  nbor_inconsistency = .false.
  nbor => self%nbor
  do while (associated(nbor))
    if ((self%nbor%needed   .neqv. nbor%needs_me) .or. &
        (self%nbor%needs_me .neqv. nbor%needed)) then
      nbor_inconsistency = .true.
      exit
    end if
    nbor => nbor%next
  end do
  call trace%end ()
END FUNCTION nbor_inconsistency

!===============================================================================
!> Info for one nbor
!===============================================================================
SUBROUTINE nbor_info (self, task, unit)
  class(nbor_t):: self
  class(task_t):: task
  integer, optional:: unit
  real(8):: distance(3)
  integer:: unit_L
  !-----------------------------------------------------------------------------
  unit_L = shared%option (unit, io_unit%mpi)
  distance = self%task%distance(task)
  write(unit_L,'(3x,2a,i6,":",i2.2,i4,2x,3l1,2x,2l1,3x,3f9.3)') &
    'pa_t%nbor_info: kind,id,rank,level,needed,needs_me,dnload,BV,pos= ', &
    self%task%kind(1:6),self%task%id, self%task%level, self%task%rank, &
    self%needed, self%needs_me, self%download, self%task%is_set(bits%boundary), &
    self%task%is_set(bits%virtual), distance/(0.5d0*(task%size+self%task%size))
END SUBROUTINE nbor_info

!===============================================================================
!> Info for all nbors
!===============================================================================
SUBROUTINE info (self, label, unit)
  class(link_t):: self
  character(*), optional:: label
  integer, optional:: unit
  !.............................................................................
  class(nbor_t), pointer:: nbor
  integer:: unit_L
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%info')
  unit_L = shared%option (unit, io_unit%mpi)
  if (present(label)) then
    write (unit_L,'(f12.6,2x,a)') wallclock(), trim(label)
  end if
  call self%task%task_info (unit=unit_L)
  nbor => self%nbor
  do while (associated(nbor))
    call nbor%nbor_info (self%task, unit=unit_L)
    nbor => nbor%next
  end do
  flush (io_unit%mpi)
  call trace%end ()
END SUBROUTINE info

!===============================================================================
!> Stubs
!===============================================================================
LOGICAL FUNCTION has_support (self)
  class(link_t):: self
  !-----------------------------------------------------------------------------
  has_support = .false.
  call io%abort ('link_t%has_support')
END FUNCTION has_support
!===============================================================================
SUBROUTINE request_support (self, source)
  class(link_t)          :: self
  class(task_t), pointer :: source
  call io%abort ('link_t%request_support called' )
END SUBROUTINE request_support

!===============================================================================
!> Monitor change of task status and take appropriate action.  This procedure
!> is called from task_list%update(), after generating a new nbor list -- which
!> as a byproduct sets new task status bits.  The first time this routine is
!> called the task%previous_status is zero, and therefore different from the
!> task%status, but none of the tests are triggered, so no action is taken,
!> except for copying the task%status to task%previous_status, so the next
!> time this routine is called there is a valid previous status to compare with
!> and act upon.
!===============================================================================
SUBROUTINE check_status (self, label)
  class(link_t):: self
  character(*):: label
  !.............................................................................
  class(task_t), pointer:: task
  type(task_t):: previous
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%check_status')
  task => self%task
  previous%status = task%previous_status
  call task%log (trim(label)//': check_status', 4)
  if (task%status /= previous%status) then
    if (task%is_set (bits%boundary)) then
      if (previous%is_set (bits%internal)) then
        !-----------------------------------------------------------------------
        ! Change from internal to boundary: set tag for sending "new" task
        !-----------------------------------------------------------------------
        call log ('internal', 'boundary')
        !call task%set (bits%send_new)
      else if (previous%is_set (bits%virtual)) then
        !-----------------------------------------------------------------------
        ! Change from virtual to boundary; the task may already exist at the
        ! recv end, and the task already has an active sequence number
        !-----------------------------------------------------------------------
        call log ('virtual', 'boundary')
      end if
    else if (task%is_set (bits%virtual)) then
      if (previous%is_set (bits%internal)) then
        !-----------------------------------------------------------------------
        ! Should never happen
        !-----------------------------------------------------------------------
        call log ('internal', 'virtual')
        call io%warning ('link_t%check_status: virtual -> internal')
      else if (previous%is_set (bits%boundary)) then
        !-----------------------------------------------------------------------
        ! Change from boundary to virtual; automatically starts receiving
        !-----------------------------------------------------------------------
        call log ('boundary', 'virtual')
      else if (previous%is_set (bits%external)) then
        !-----------------------------------------------------------------------
        ! Change from boundary to virtual; automatically starts receiving
        !-----------------------------------------------------------------------
        call log ('external', 'virtual')
      end if
    else if (task%is_set (bits%internal)) then
      if (previous%is_set (bits%boundary)) then
        !-----------------------------------------------------------------------
        ! Change from boundary to internal: tag to remove sent_list
        !-----------------------------------------------------------------------
        call log ('boundary', 'internal')
        call task%set (bits%remove_sent)
      else if (previous%is_set (bits%virtual)) then
        !-----------------------------------------------------------------------
        ! Should never happen
        !-----------------------------------------------------------------------
        call log ('virtual', 'internal')
        call io%warning ('task_list_t%check_status: internal -> virtual')
      end if
    end if
    task%previous_status = task%status
  end if
  call trace%end ()
contains
  !=============================================================================
  subroutine log (from, to)
    character(*):: from, to
    integer:: unit
    !---------------------------------------------------------------------------
    if (verbose > 1) then
      write (stderr,'(f12.6,i5,2x,a,i6," changed from ",a," to ",a)') &
        wallclock(), mpi%rank, 'link_t%check_status: task', task%id, from, to
    end if
    if (verbose > 0) then
      write (stdlog,'(f12.6,2x,a,i6," changed from ",a," to ",a)') &
        wallclock(), 'link_t%check_status: task', task%id, from, to
      flush (stdlog)
    end if
    call task%log (from//' -> '//to, 2)
  end subroutine
END SUBROUTINE check_status

END MODULE link_mod
