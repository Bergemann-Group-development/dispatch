MODULE task_list_mod
  USE nbor_mod
  USE link_mod
  USE list_mod
  USE task_mesg_mod
  implicit none
  private
  type, public, extends(task_mesg_t):: task_list_t
  contains
    procedure:: init
  end type
CONTAINS

SUBROUTINE init (self)
  class(task_list_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%init')
  call self%task_mesg_t%init
  call trace%end ()
END SUBROUTINE init

END MODULE task_list_mod