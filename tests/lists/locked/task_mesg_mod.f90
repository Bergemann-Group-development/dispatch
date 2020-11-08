MODULE task_mesg_mod
  USE nbor_mod
  USE link_mod
  USE list_mod
  implicit none
  private
  type, public, extends(list_t):: task_mesg_t
  contains
    procedure:: init
  end type
CONTAINS

SUBROUTINE init (self)
  class(task_mesg_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%init')
  call self%list_t%init
  call trace%end ()
END SUBROUTINE init

END MODULE task_mesg_mod