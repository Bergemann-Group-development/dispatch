!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  implicit none
  integer, parameter:: ntime=7
  type task_t
    integer:: id = 0
    integer:: it = 0
    integer:: rank = 0
    integer:: status = 0
    real(8):: time = 0d0
    real(8):: dtime = 0d0
    real(8):: pos(3) = 0d0
    real(8):: vel(3) = 0d0
    real(8):: size(3) = 0d0
    contains
    procedure init, allocate, deallocate, update, print
  end type task_t
PRIVATE
PUBLIC task_t
CONTAINS

!> =============================================================================
SUBROUTINE init (self,id)
  class(task_t):: self
  integer:: id
  !.............................................................................
  self%id = id
END SUBROUTINE init

!> =============================================================================
SUBROUTINE allocate (self)
  class(task_t):: self
  !.............................................................................
END SUBROUTINE allocate

!> =============================================================================
SUBROUTINE deallocate (self)
  class(task_t):: self
  !.............................................................................
END SUBROUTINE deallocate

!> =============================================================================
SUBROUTINE update (self)
  class(task_t):: self
  !.............................................................................
  self%time = self%time + self%dtime
END SUBROUTINE update

!> =============================================================================
SUBROUTINE print (self)
  class(task_t):: self
  !.............................................................................
  print*, self%id, self%time
END SUBROUTINE print

!> =============================================================================
END MODULE task_mod
