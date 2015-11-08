!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for planet task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE planet_mod
  USE task_mod
  USE particle_mod
  implicit none
  !-----------------------------------------------------------------------------
  type, extends(timeslots_t):: planet_t
    real:: mass, radius, surface_temperature
  end type
PRIVATE
PUBLIC planet_t
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
END MODULE
