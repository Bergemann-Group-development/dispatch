!> *****************************************************************************
!> Module with list handling for generic class tasks_t objects
!> *****************************************************************************
MODULE list_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE patch_mod
  implicit none
  private
  type, public:: nbor_t
    class(task_t), pointer :: task => null()          ! task
    class(link_t), pointer :: nbor => null()          ! next nbor
    class(link_t), pointer :: link => null()          ! nbor's link (maybe not needed)
  end type
  type, extends(nbor_t), public:: link_t
    type(link_t),  pointer :: next => null()          ! next link in task list
    type(link_t),  pointer :: prev => null()          ! previous link in task list
    type(link_t),  pointer :: parent => null()        ! parent task link
    type(link_t),  pointer :: next_time => null()     ! next link in time order
  contains
    procedure:: download
    procedure:: add_nbor_by_quality
    procedure:: make_new_nbors
    procedure:: make_new_nbor
    procedure:: print_nbors
  end type
  type, public:: list_t
    character(len=64):: name=''
    type(link_t), pointer:: head => null()
    type(link_t), pointer:: tail => null()
    type(link_t), pointer:: queue => null()
    integer:: n = 0
    integer:: nq = 0
  contains
    procedure:: init
    procedure:: update
    procedure:: reset
    procedure:: append
    procedure:: remove_task
    procedure:: remove_link
    procedure:: remove_head
    procedure:: queue_by_time
    procedure:: print => print_list
  end type
  public task2patch
CONTAINS

!> =============================================================================
!> Initialize list
!> =============================================================================
SUBROUTINE init (self, name)
  class(list_t):: self
  character(len=*), optional:: name
  !.............................................................................
  if (present(name)) then
    self%name = name
    print*,'init_list: ', name
  end if
END SUBROUTINE init

!> =============================================================================
!> Update simple list
!> =============================================================================
SUBROUTINE update (self)
  class(list_t):: self
END

!> =============================================================================
!> Reset list, by deallocating memory used by tasks, and deallocating tasks and
!> links
!> =============================================================================
SUBROUTINE reset (self)
  class(list_t):: self
  type(link_t), pointer:: link, old
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin ('list_t%reset '//self%name)
  link => self%head
  do while (associated(link))
    task => link%task
    if (associated(link)) then
      call task%deallocate
      deallocate (task)
    end if
    old => link
    link => link%next
    deallocate (old)
  end do
  nullify (self%head)
  nullify (self%tail)
  call trace_end
END SUBROUTINE reset

!> =============================================================================
!> =============================================================================
SUBROUTINE append (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin ('list_t%append '//self%name)
  allocate (link)
  link%task => task
  if (associated(self%tail)) then
    link%prev => self%tail
    link%prev%next => link
    self%tail => link
  else
    self%head => link
    self%tail => link
  end if
  self%n = self%n+1
  call trace_end
END SUBROUTINE append

!> =============================================================================
!> Remove a link from the list, taking special care with the head and tail
!> =============================================================================
SUBROUTINE remove_link (self, link)
  class(list_t):: self
  type(link_t), pointer:: link, next
  !.............................................................................
  call trace_begin ('list_t%remove_link '//self%name)
  next => link%next
  !-------------------------------------------------------------------------
  ! Check for tail task
  !-------------------------------------------------------------------------
  if (associated(next)) then
    next%prev => link%prev                  ! point backwards, possibly null
  else
    self%tail => link%prev                  ! new tail link
  end if
  !-------------------------------------------------------------------------
  ! Check for head task
  !-------------------------------------------------------------------------
  if (associated(link%prev)) then           ! not head task
    link%prev%next => next                  ! make prev point to next
  else                                      ! link was head
    self%head => next                       ! make head point to next
  end if
  deallocate (link)                         ! free to remove link
  self%n = self%n-1
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Remove the head link
!> =============================================================================
SUBROUTINE remove_head (self)
  class(list_t):: self
  type(link_t), pointer:: head, link, next
  !.............................................................................
  call trace_begin ('list_t%remove_head '//self%name)
  head => self%head
  next => head%next
  !-------------------------------------------------------------------------
  ! Check for tail task (only one task left)
  !-------------------------------------------------------------------------
  if (associated(next)) then
    next%prev => head%prev                  ! point backwards, possibly null
  else
    self%tail => head%prev                  ! new tail link
  end if
  !-------------------------------------------------------------------------
  ! Remove the head link
  !-------------------------------------------------------------------------
  self%head => next                         ! make head point to next
  self%n = self%n-1
  call trace_end
END SUBROUTINE remove_head

!> =============================================================================
!> Remove an task from the list, by removing the link pointing to it
!> =============================================================================
SUBROUTINE remove_task (self, task)
  class(list_t):: self
  type(link_t), pointer:: link, prev
  class(task_t), pointer:: task
  integer:: i
  !.............................................................................
  call trace_begin ('list_t%remove '//self%name)
  !-----------------------------------------------------------------------------
  ! Prevent damage to the time queue when removing the task
  !-----------------------------------------------------------------------------
!  print *, 'patching the queue when removing task', task%id, self%nq, self%n
  nullify(prev)
  link => self%queue
  i = 0
  do while (associated(link))
    i = i+1
    if (associated(link%task,task)) then
      if (associated(prev)) then
!        print*, 'patching previous link', prev%task%id, associated(link%next_time)
        prev%next_time => link%next_time
      else
!        print*, 'patching self%queue', prev%task%id, associated(link%next_time)
        self%queue => link%next_time
      end if
    end if
    prev => link
    link => link%next_time
  end do
  !-----------------------------------------------------------------------------
  ! Now find the task link, and remove it
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    if (associated(link%task,task)) then
      call self%remove_link (link)
      call trace_end
      return
    end if
    link => link%next
  end do
  call trace_end
END SUBROUTINE remove_task

!> =============================================================================
!> Add a link node into the queue chain, in ascending time order
!> =============================================================================
SUBROUTINE queue_by_time (self, this)
  class(list_t):: self
  class(link_t), pointer:: this, next, prev
  !.............................................................................
  call trace_begin ('list_t%queue_by_time '//self%name)
  nullify (this%next_time)
  nullify (prev)
!  print*, self%nq, associated(next)
  !-----------------------------------------------------------------------------
  ! Run through the queue link chain
  !-----------------------------------------------------------------------------
  next => self%queue
  do while (associated(next))
    !---------------------------------------------------------------------------
    ! Point to the next link in the time chain, which may or may not exist
    !---------------------------------------------------------------------------
    this%next_time => next
!    print *, this%task%id, this%task%time, next%task%id, next%task%time
    if (next%task%time >= this%task%time) then
      !-------------------------------------------------------------------------
      ! The next link time is larger or the same, so insert this link before
      !-------------------------------------------------------------------------
      self%nq = self%nq+1
      if (associated(prev)) then
        !-----------------------------------------------------------------------
        ! If the previous time link exist, just redirect its next_time to this
        !-----------------------------------------------------------------------
        prev%next_time => this
!        print*, 'inserted1: ', prev%task%id, this%task%id, this%next_time%task%id, self%nq
      else
        !-----------------------------------------------------------------------
        ! If the previous time link doesn't exist we should be the head of the list
        !-----------------------------------------------------------------------
        self%queue => this
!        print*, 'inserted2: ', this%task%id, self%nq
      end if
      call trace_end
      return
    end if
    prev => next
    next => next%next_time
  end do
  !-----------------------------------------------------------------------------
  ! No time was larger or equal to this time, so it should be appended
  !-----------------------------------------------------------------------------
  self%nq = self%nq+1
  if (associated(prev)) then
    prev%next_time => this
    nullify (this%next_time)
!    print*, 'inserted3: ', this%task%id, self%nq
  else
    !---------------------------------------------------------------------------
    ! If there is nothing to append to, this should be the head
    !---------------------------------------------------------------------------
    self%queue => this
!    print*, 'inserted4: ', this%task%id, self%nq
  end if
  call trace_end
END SUBROUTINE queue_by_time

!> =============================================================================
!> Insert 'new' task before 'old' task
!> =============================================================================
SUBROUTINE insert (self, old, new)
  class(list_t):: self
  class(link_t), pointer:: old, new
  !.............................................................................
  call trace_begin ('list_t%insert '//self%name)
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    self%head => new
  end if
  new%next => old
  new%next%prev => new
  self%n = self%n + 1
  call trace_end
END SUBROUTINE insert

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
!> Run through neighbor list and pick up guard cells, after which the whole
!> task is valid for download into other patches guard cells
!> =============================================================================
SUBROUTINE download (self, all_cells)
  implicit none
  class (link_t):: self
  logical, optional:: all_cells
  class (link_t),  pointer:: nbor
  class (patch_t), pointer:: patch, nbor_patch
  class (task_t), pointer:: task
  integer :: nguard, special
  real:: dmin
  !.............................................................................
  call trace_begin('link_t%download')
  patch => task2patch(self%task)
  nguard=0
  nbor => self%nbor
  do while (associated (nbor))
    task => nbor%task
    nbor_patch => task2patch(task)
    if (abs(patch%level-nbor_patch%level) < 2) then
      if (nbor_patch%level == patch%level+1) special = 6
      call nbor_patch%download_to (task, special, all_cells)
      nguard = nguard+1
    end if
    nbor => nbor%next
  end do
  if (io%verbose > 5) &
    print *, 'download: Number of neighbors with guard cells:', nguard
  call trace_end
END SUBROUTINE download

!> =============================================================================
!> Initialize the nbors chain of a new task (self), by
!> 1) adding the parent
!> 2) searching the parent nbors, adding patches that overlap
!> Each time it adds an nbor, it should also add itself to that nbor's nbor list
!> =============================================================================
SUBROUTINE make_new_nbors (self, selfp, parent)
  class(link_t):: self
  class(link_t), pointer:: selfp, parent, nbor
  class(task_t), pointer:: t1, t2
  class(patch_t), pointer:: p1, p2
  !> ...........................................................................
  call trace_begin('task_list_t%make_new_nbors', 1)
  p1 => task2patch(self%task)
  call self%make_new_nbor (selfp, parent)         ! make me a neighbor of my parent
  nbor => parent%nbor                             ! start on parent's nbors chain
  do while (associated(nbor))
    p2 => task2patch(nbor%task)
    if (.not.associated(self%task,nbor%task) &
      .and. p1%overlaps(p2)) &   ! if my task overlaps with that task
      call self%make_new_nbor (selfp, nbor%nbor)  ! make a neighbor (and make me one)
    nbor => nbor%next                             ! continue on parents nbors chain
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
  call trace_begin('task_list_t%make_new_nbor', 1)
  allocate (link1)                              ! for our own nbors chain
  link1%task => that%task                     ! the link points to that task
  link1%nbor => that                            ! link to the nbor link
  link1%parent => that%parent                   ! direct link to task parent
  call self%add_nbor_by_quality (link1)         ! add to that chain in quality order
  allocate (link2)                              ! for that's nbors chain
  link2%task => self%task                     ! that's new link points to my task
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
  type(link_t), pointer:: next, prev
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
        return
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
  call trace_end
END SUBROUTINE add_nbor_by_quality

!> =============================================================================
!> Print a list of nbors
!> =============================================================================
SUBROUTINE print_nbors (self, label)
  class (link_t):: self
  character(len=*), optional:: label
  type(link_t) , pointer:: link, nbor
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
!> =============================================================================
SUBROUTINE print_list (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  type(link_t), pointer:: link, prev
  class(task_t), pointer:: task
  character(len=32):: type
  !.............................................................................
  call trace_begin ('list_t%print')
  print*,'list_t%print: ',trim(self%name), self%n
  link => self%head
  do while (associated(link))
    task => link%task
    type = ''
    !select type (task)
    !type is (task_t)
    !  type = 'task_t'
    !type is (particle_t)
    !  type = 'particle_t'
    !type is (planet_t)
    !  type = 'planet_t'
    !type is (star_t)
    !  type = 'star_t'
    !type is (pebble_t)
    !  type = 'pebble_t'
    !class is (patch_t)
      type = 'patch_t'
    !end select
    if (io%master) print '(a,i7,f10.6,3x,3f10.6,i4,i6,l5,2x,a)', 'list_t%print:', &
      task%id, task%time, task%position, task%status, task%n_nbors, associated(link%nbor), &
      type
    link => link%next
  end do
  call trace_end
END SUBROUTINE print_list

END MODULE list_mod
