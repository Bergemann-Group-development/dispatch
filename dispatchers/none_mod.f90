!> =============================================================================
!> Do not use a dispatcher, but call task_list%execute, which relies on threads
!> handling a ready queue.
!> =============================================================================
MODULE dispatcher_mod
  USE task_list_mod
  implicit none
  private
  type, public:: dispatcher_t
  contains
    procedure, nopass:: execute
  end type
  type(dispatcher_t), public:: dispatcher
CONTAINS

!> =============================================================================
!> Execute a task list
!> =============================================================================
SUBROUTINE execute (task_list)
  type(task_list_t):: task_list
  !-----------------------------------------------------------------------------
  call task_list%execute
END SUBROUTINE execute

END MODULE dispatcher_mod
