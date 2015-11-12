!> *****************************************************************************
!> Module with list handling for generic class tasks_t objects
!> *****************************************************************************
MODULE list_mod
  USE item_mod
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
END SUBROUTINE reset

!> =============================================================================
!> =============================================================================
SUBROUTINE append (list, link)
  class(list_t):: list
  class(link_t), pointer:: link
  !.............................................................................
  if (associated(list%tail)) then
    link%prev => list%tail
    link%prev%next => link
    list%tail => link
  else
    list%head => link
    list%tail => link
  end if
  list%size = list%size+1
END SUBROUTINE append

!> =============================================================================
!> =============================================================================
SUBROUTINE add_by_time (list, new)
  class(list_t):: list
  class(link_t), pointer:: link, new
  !.............................................................................
  link => list%head
  do while (associated(link))
    if (link%item%time > new%item%time) then
      call insert (list, link, new)
      list%size = list%size+1
      return
    end if
    link => link%next
  end do
  call append (list, new)
  list%size = list%size+1
END SUBROUTINE add_by_time

!> =============================================================================
!> Insert 'new' item before 'link' item
!> =============================================================================
SUBROUTINE insert (list, old, new)
  class(list_t):: list
  class(link_t), pointer:: old, prev, next, new
  !.............................................................................
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    list%head => new
  end if
  new%next => old
  new%next%prev => new
  list%size = list%size + 1
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
