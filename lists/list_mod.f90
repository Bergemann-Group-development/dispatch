!> *****************************************************************************
!> Module with list handling for generic class tasks_t objects
!> *****************************************************************************
MODULE list_mod
  USE task_mod
  USE trace_mod
  implicit none
  type link_t
    class(task_t), pointer :: task => null()
    type(link_t),  pointer :: next => null()
    type(link_t),  pointer :: prev => null()
    type(link_t),  pointer :: nbor => null()
  end type
  type list_t
    character(len=64):: name
    type(link_t), pointer:: head => null()
    type(link_t), pointer:: tail => null()
    integer:: size = 0
    contains
    procedure:: init
    procedure:: update
    procedure:: reset
    procedure:: append
    procedure:: remove_task
    procedure:: remove_link
    procedure:: remove_head
    procedure:: add_by_time
    procedure:: print
  end type
PRIVATE
PUBLIC list_t, link_t
CONTAINS

!> =============================================================================
!> Initialize list
!> =============================================================================
SUBROUTINE init (self, name)
  class(list_t):: self
  character(len=*):: name
  !.............................................................................
  self%name = name
  print*,'init_list: ', name
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
  self%size = self%size+1
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
  self%size = self%size-1
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
  self%size = self%size-1
  call trace_end
END SUBROUTINE remove_head

!> =============================================================================
!> Remove an task from the list, by removing the link pointing to it
!> =============================================================================
SUBROUTINE remove_task (self, task)
  class(list_t):: self
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin ('list_t%remove '//self%name)
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
!> =============================================================================
SUBROUTINE add_by_time (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link, new
  !.............................................................................
  call trace_begin ('list_t%add_by_time '//self%name)
  link => self%head
  do while (associated(link))
    if (link%task%time > task%time) then
      allocate (new)
      new%task => task
      call insert (self, link, new)
      self%size = self%size+1
      call trace_end
      return
    end if
    link => link%next
  end do
  call append (self, task)
  self%size = self%size+1
  call trace_end
END SUBROUTINE add_by_time

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
  self%size = self%size + 1
  call trace_end
END SUBROUTINE insert

!> =============================================================================
!> =============================================================================
SUBROUTINE print (self)
  class(list_t):: self
  type(link_t), pointer:: link, prev
  class(task_t), pointer:: task
  !.............................................................................
  link => self%head
  print*,'print: ',trim(self%name), self%size
  do while (associated(link))
    task => link%task
    call task%print
    link => link%next
  end do
END SUBROUTINE print
END MODULE list_mod
