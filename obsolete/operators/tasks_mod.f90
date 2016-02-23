!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE tasks_mod
  USE io_mod
  USE task_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Time slot info and procedures added in tasks_t
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: tasks_t
    real(8), dimension(:), pointer:: t
    real(8), dimension(:), pointer:: dt
    integer, dimension(:), pointer:: iit
    integer:: nt
    integer:: it = 1
    integer:: new = 2
    integer:: verbose = 1
  contains
    !procedure:: overlaps
    !procedure:: send_to
  end type
CONTAINS

END MODULE
