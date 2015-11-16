!> *****************************************************************************
!> Module with list handling for generic class tasks_t objects
!> *****************************************************************************
MODULE list_mod
  USE item_mod
  USE trace_mod
  implicit none
  type link_t
    class(item_t), pointer :: item => null()
    type(link_t),  pointer :: next => null()
    type(link_t),  pointer :: prev => null()
  end type
  type list_t
    character(len=64):: name
    type(link_t), pointer:: head => null()
    type(link_t), pointer:: tail => null()
    type(link_t), pointer:: item => null()
    integer:: size = 0
    contains
    procedure:: init
    procedure:: reset
    procedure:: append
    procedure:: remove
    procedure:: remove_link
    procedure:: add_by_time
    procedure:: print
  end type
PRIVATE
PUBLIC list_t, link_t
CONTAINS

!> =============================================================================
!> Initialize list
!> =============================================================================
SUBROUTINE init (list, name)
  class(list_t):: list
  character(len=*):: name
  !.............................................................................
  list%name = name
  print*,'init_list: ', name
END SUBROUTINE init

!> =============================================================================
!> Reset list, by deallocating memory used by items, and deallocating items and
!> links
!> =============================================================================
SUBROUTINE reset (list)
  class(list_t):: list
  type(link_t), pointer:: link, old
  class(item_t), pointer:: item
  !.............................................................................
  call trace_begin ('list_t%reset '//list%name)
  link => list%head
  do while (associated(link))
    item => link%item
    if (associated(link)) then
      call item%deallocate
      deallocate (item)
    end if
    old => link
    link => link%next
    deallocate (old)
  end do
  nullify (list%head)
  nullify (list%tail)
  nullify (list%item)
  call trace_end
END SUBROUTINE reset

!> =============================================================================
!> Remove a link from the list, taking special care with the head and tail
!> =============================================================================
SUBROUTINE remove_link (list, link)
  class(list_t):: list
  type(link_t), pointer:: link, next
  !.............................................................................
  call trace_begin ('list_t%remove_link '//list%name)
  next => link%next
  !-------------------------------------------------------------------------
  ! Check for tail item
  !-------------------------------------------------------------------------
  if (associated(next)) then
    next%prev => link%prev                  ! point backwards, possibly null
  else
    list%tail => link%prev                  ! new tail link
  end if
  !-------------------------------------------------------------------------
  ! Check for head item
  !-------------------------------------------------------------------------
  if (associated(link%prev)) then           ! not head item
    link%prev%next => next                  ! make prev point to next
  else                                      ! link was head
    list%head => next                       ! make head point to next
  end if
  deallocate (link)                         ! free to remove link
  list%size = list%size-1
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Remove a link from the list, taking special care with the head and tail
!> =============================================================================
SUBROUTINE remove (list, this)
  class(list_t):: list
  type(link_t), pointer:: link
  class(item_t), pointer:: this
  !.............................................................................
  call trace_begin ('list_t%remove '//list%name)
  link => list%head
  do while (associated(link))
    if (associated(link%item,this)) then
      call list%remove_link (link)
      call trace_end
      return
    end if
    link => link%next
  end do
  call trace_end
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE append (list, item)
  class(list_t):: list
  class(item_t), pointer:: item
  class(link_t), pointer:: link
  !.............................................................................
  !call trace_begin ('list_t%append '//list%name)
  allocate (link)
  link%item => item
  if (associated(list%tail)) then
    link%prev => list%tail
    link%prev%next => link
    list%tail => link
  else
    list%head => link
    list%tail => link
  end if
  list%size = list%size+1
  !call trace_end
END SUBROUTINE append

!> =============================================================================
!> =============================================================================
SUBROUTINE add_by_time (list, item)
  class(list_t):: list
  class(item_t), pointer:: item
  class(link_t), pointer:: link, new
  !.............................................................................
  call trace_begin ('list_t%add_by_time '//list%name)
  link => list%head
  do while (associated(link))
    if (link%item%time > item%time) then
      allocate (new)
      new%item => item
      call insert (list, link, new)
      list%size = list%size+1
      call trace_end
      return
    end if
    link => link%next
  end do
  call append (list, item)
  list%size = list%size+1
  call trace_end
END SUBROUTINE add_by_time

!> =============================================================================
!> Insert 'new' item before 'link' item
!> =============================================================================
SUBROUTINE insert (list, old, new)
  class(list_t):: list
  class(link_t), pointer:: old, prev, next, new
  !.............................................................................
  call trace_begin ('list_t%insert '//list%name)
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    list%head => new
  end if
  new%next => old
  new%next%prev => new
  list%size = list%size + 1
  call trace_end
END SUBROUTINE insert

!> =============================================================================
!> =============================================================================
SUBROUTINE print (list)
  class(list_t):: list
  type(link_t), pointer:: link, prev
  class(item_t), pointer:: item
  !.............................................................................
  link => list%head
  print*,'print: ',trim(list%name), list%size
  do while (associated(link))
    item => link%item
    call item%print
    link => link%next
  end do
END SUBROUTINE print
END MODULE list_mod
