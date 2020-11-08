MODULE list_mod
  USE trace_mod
  USE link_mod
  implicit none
  private
  type, public:: list_t
    class(link_t), pointer:: prev => null()
    class(link_t), pointer:: next => null()
    integer:: n
  contains
    procedure:: init
  end type
CONTAINS

SUBROUTINE init (self)
  class(list_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%init')
  call trace%end ()
END SUBROUTINE init

END MODULE list_mod
