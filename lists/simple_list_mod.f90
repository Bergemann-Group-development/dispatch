!> *****************************************************************************
!> Module for list that only have forward links
!> *****************************************************************************
MODULE simple_list_mod
  USE item_mod
  USE list_mod
  implicit none
  type simple_link_t
    class(item_t), pointer :: item => null()
    type(simple_link_t),  pointer :: next => null()
  end type
  type simple_list_t
    character(len=64):: name
    type(simple_link_t), pointer:: head => null()
    type(simple_link_t), pointer:: tail => null()
    type(simple_link_t), pointer:: item => null()
    integer(8):: n = 0
    contains
    procedure:: init
    procedure:: purge
    procedure:: append
  end type
PRIVATE
PUBLIC simple_list_t, simple_link_t
CONTAINS

!> =============================================================================
!> Initialize list
!> =============================================================================
SUBROUTINE init (self, name, list, item)
  class(simple_list_t):: self
  character(len=*):: name
  class(list_t), optional, pointer:: list
  class(item_t), optional, pointer:: item
  !.............................................................................
  self%name = name
  call self%purge
  print*,'init_list: ', name
END SUBROUTINE init

!> =============================================================================
!> Remove all links
!> =============================================================================
SUBROUTINE purge (list)
  class(simple_list_t):: list
  class(simple_link_t), pointer:: link, prev
  !.............................................................................
  link => list%head
  do while (associated(link))
    prev => link
    link => link%next
    deallocate (prev)
  end do
  nullify (list%head, list%tail)
  list%n = 0
END SUBROUTINE purge

!> =============================================================================
!> =============================================================================
SUBROUTINE append (list, item)
  class(simple_list_t):: list
  class(simple_link_t), pointer:: link
  class(item_t), optional, pointer:: item
  !.............................................................................
  allocate (link)
  link%item => item
  if (associated(list%tail)) then
    list%tail%next => link
    list%tail => link
  else
    list%head => link
    list%tail => link
  end if
  list%n = list%n+1
END SUBROUTINE append

END MODULE simple_list_mod
