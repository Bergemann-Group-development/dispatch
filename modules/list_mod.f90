!> *****************************************************************************
!> Module with list handling for generic class task_t objects
!> *****************************************************************************
MODULE list_mod
  USE task_mod
  implicit none
  type link_t
    class(task_t), pointer :: item => null()
    type(link_t),  pointer :: next => null()
    type(link_t),  pointer :: prev => null()
  end type
  type list_t
    character(len=64):: name
    type(link_t), pointer:: head => null() 
    type(link_t), pointer:: tail => null() 
    type(link_t), pointer:: item => null()
    integer(8):: n = 0
    contains
    procedure:: init, reset, append, add_by_time, insert, print
  end type
PRIVATE
PUBLIC list_t, link_t
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (list, name)
  class(list_t):: list
  character(len=*):: name
  !.............................................................................
  list%name = name
  print*,'init_list: ', name
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE reset (list)
  class(list_t):: list
  type(link_t), pointer:: link
  class(task_t), pointer:: item
  !.............................................................................
  link => list%head
  do while (associated(link))
    item => link%item
    if (associated(link)) then
      call item%deallocate
      deallocate (item)
    end if
    link => link%next
  end do
  nullify (list%head)
  nullify (list%tail)
  nullify (list%item)
END SUBROUTINE reset

!> =============================================================================
!> =============================================================================
SUBROUTINE append (list, link)
  class(list_t):: list
  class(task_t), pointer:: item
  type(link_t), pointer:: link
  !.............................................................................
  item => link%item
  if (associated(list%tail)) then
    link%prev => list%tail
    link%prev%next => link
    list%tail => link
  else
    list%head => link
    list%tail => link
  end if
  list%n = list%n+1
END SUBROUTINE append

!> =============================================================================
!> =============================================================================
SUBROUTINE add_by_time (list, new)
  class(list_t):: list
  type(link_t), pointer:: link, new
  !.............................................................................
  link => list%head
  do while (associated(link))
    if (link%item%time > new%item%time) then
      call insert (list, link, new)
      list%n = list%n+1
      return
    end if
    link => link%next
  end do
  call append (list, new)
  list%n = list%n+1
END SUBROUTINE add_by_time

!> =============================================================================
!> Insert 'new' item before 'link' item
!> =============================================================================
SUBROUTINE insert (list, old, new)
  class(list_t):: list
  type(link_t), pointer:: old, prev, next, new
  class(task_t), pointer:: item
  !.............................................................................
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    list%head => new
  end if
  new%next => old
  new%next%prev => new
  list%n = list%n + 1
END SUBROUTINE insert

!> =============================================================================
!> =============================================================================
SUBROUTINE print (list)
  class(list_t):: list
  type(link_t), pointer:: link, prev
  class(task_t), pointer:: item
  !.............................................................................
  link => list%head
  print*,'print: ',list%name, list%n
  do while (associated(link))
    item => link%item
    call item%print
    link => link%next
  end do
END SUBROUTINE print
END MODULE list_mod
