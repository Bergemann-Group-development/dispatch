!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE tasks_mod
  USE io_mod
  USE task_mod
  USE trace_mod
  implicit none
  private
  integer, public, parameter:: ntime=7
  !-----------------------------------------------------------------------------
  ! Time slot info and procedures added in tasks_t
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: tasks_t
    real(8):: t(ntime) = 0d0
    real(8):: dt(ntime) = 1d0
    integer:: iit(ntime) = 1
    integer:: nt = ntime
    integer:: it = 1
    integer:: verbose = 1
  contains
    procedure:: overlaps
    procedure:: send_to
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
FUNCTION overlaps (self, task)
  class(tasks_t):: self
  class(task_t):: task
  logical:: overlaps
  call trace_begin ('task_t%overlap')
  overlaps = .true.
  call trace_end
END FUNCTION overlaps

!> =============================================================================
!> =============================================================================
SUBROUTINE send_to (self, task)
  class(tasks_t):: self
  class(task_t):: task
  call trace_begin ('tasks_t%send_to')
  print *, 'tasks_t%send_to', task%rank
  call trace_end
END SUBROUTINE send_to

END MODULE
