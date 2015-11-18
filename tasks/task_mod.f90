!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for tasks.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  USE io_mod
  USE mpi_mod
  USE bits_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type, public:: task_t
    integer:: id = 0
    integer:: status = 0
    integer:: rank = 0
    integer:: n_check = 0
    real(8):: time = 0d0
    real(8):: end_time = 0d0
    real(8):: dtime = 0.1d0
    real(4):: grace = 0.0
    real(8):: position(3) = 0d0
    real(8):: velocity(3) = 0d0
    real(8):: size(3) = 0d0
  contains
    procedure:: init
    procedure:: update
    procedure:: ahead_of
    procedure:: set
    procedure:: clear
    procedure:: print
    procedure:: is_set
    procedure:: is_clear
    procedure:: send_to
    procedure:: allocate => void
    procedure:: deallocate => void
  end type
  integer:: id=0
CONTAINS

!> =============================================================================
SUBROUTINE void (self)
  class (task_t):: self
  print*, 'ERROR: task_t%void called'
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class (task_t):: self
  !.............................................................................
  call trace_begin ('task_t%init')
  id = id+1
  self%id = id
  call trace_end
END SUBROUTINE

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class (task_t):: self
  !.............................................................................
  call trace_begin ('task_t%update')
  self%time = self%time + self%dtime
  call trace_end
END SUBROUTINE

!> =============================================================================
!> =============================================================================
LOGICAL FUNCTION ahead_of (self, task)
  class (task_t):: self, task
  !.............................................................................
  ahead_of = self%time > task%time + task%grace
END FUNCTION

!> =============================================================================
!> =============================================================================
SUBROUTINE send_to (self, task)
  class(task_t):: self
  class(task_t):: task
  call trace_begin ('tasks_t%send_to')
  print *, 'tasks_t%send_to', task%rank
  call trace_end
END SUBROUTINE send_to

!> =============================================================================
SUBROUTINE print (self)
  class (task_t):: self
  if (io%master) print*, 'task_t%print:', self%id, self%time
END SUBROUTINE

!> =============================================================================
!> Set and clear status bits
!> =============================================================================
SUBROUTINE set (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('task_t%set')
  self%status = ior(self%status, bit)
  call trace_end
END SUBROUTINE
!> =============================================================================
SUBROUTINE clear (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('task_t%clear')
  self%status = iand(self%status, not(bit))
  call trace_end
END SUBROUTINE
!> =============================================================================
LOGICAL FUNCTION is_set (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('task_t%is_set')
  is_set = (ior(self%status, bit) /= 0)
  call trace_end
END FUNCTION
!> =============================================================================
LOGICAL FUNCTION is_clear (self, bit)
  class (task_t):: self
  integer:: bit
  !.............................................................................
  call trace_begin ('task_t%is_clear')
  is_clear = (iand(self%status, bit) == 0)
  call trace_end
END FUNCTION

END MODULE task_mod
