MODULE link_mod
  USE trace_mod
  USE dll_mod
  USE task_mod
  implicit none
  private
  !.............................................................................
  type, public:: nbor_t
    class(nbor_t), pointer:: next => null()
    class(link_t), pointer:: link
    class(task_t), pointer:: task => null()
    real:: u_max
    integer:: n
  contains
    procedure:: init => init_nbor
  end type
  !.............................................................................
  type, public, extends(dll_node_t):: link_t
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: task => null()
  contains
    procedure:: init
  end type
CONTAINS

SUBROUTINE init_nbor (self)
  class(nbor_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('nbor_t%init')
  call trace%end ()
END SUBROUTINE init_nbor

SUBROUTINE init (self)
  class(link_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%init')
  call self%dll_node_t%init
  call trace%end ()
END SUBROUTINE init

END MODULE link_mod
