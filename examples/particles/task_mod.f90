!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE io_mod
  USE trace_mod
  USE bits_mod
  implicit none
  private
  type, public:: task_t
    integer:: id
    integer:: status
    integer:: level
    real:: time=0.0, dtime=0.1, end_time=1.0, grace=0.0
    real:: out_time=0.1, out_next=0.0           ! controls output
  contains
    procedure:: init
    procedure:: update
    procedure:: output
    procedure:: set
    procedure:: clear
    procedure:: is_set
    procedure:: is_clear
    procedure:: ahead_of
  end type
  integer:: id = 0
CONTAINS

!> =============================================================================
!> Initialize the task.
!> =============================================================================
SUBROUTINE init (self)
  class(task_t):: self
  !.............................................................................
  call trace_begin ('task_t%init', 3)
  !$omp atomic
  id = id + 1                                   ! increment global
  if (io%verbose>2) &
    print*,'task_t%init: id=', id
  self%id = id                                  ! new task id
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Update the task
!> =============================================================================
SUBROUTINE update (self)
  class(task_t):: self
  !.............................................................................
  call trace_begin ('task_t%update', 3)
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Update the task
!> =============================================================================
SUBROUTINE output (self)
  class(task_t):: self
  !.............................................................................
  call trace_begin ('task_t%output', 3)
  call trace_end
END SUBROUTINE output

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%set')
  self%status = ior(self%status, bit)
!  call trace_end
END SUBROUTINE
!> =============================================================================
SUBROUTINE clear (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%clear')
  self%status = iand(self%status, not(bit))
!  call trace_end
END SUBROUTINE
!> =============================================================================
LOGICAL FUNCTION is_set (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%is_set')
  is_set = (iand(self%status, bit) /= 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
!  call trace_begin ('task_t%is_clear')
  is_clear = (iand(self%status, bit) == 0)
!  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION ahead_of (self, task)
  class (task_t):: self, task
  !.............................................................................
  ahead_of = self%time > task%time + task%grace
!  print*,'ahead_of: ', self%time, task%time, task%grace, ahead_of
END FUNCTION

!> =============================================================================
END MODULE task_mod
