!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for dust tasks
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pebble_mod
  USE task_mod
  USE particle_mod
  implicit none
  !-----------------------------------------------------------------------------
  type, extends(particles_t):: pebble_t
    real:: mass, radius, surface_temperature
    contains
    procedure:: init       => init_pebbles
    procedure:: allocate   => allocate_pebbles
    procedure:: deallocate => deallocate_pebbles
    procedure:: update     => update_pebbles
    procedure:: print      => print_pebbles
  end type
!PRIVATE
!PUBLIC pebble_t
CONTAINS

!> =============================================================================
SUBROUTINE init_pebbles (self)
  class(pebble_t):: self
  !.............................................................................
  call self%task_t%init
END SUBROUTINE init_pebbles

!> =============================================================================
SUBROUTINE allocate_pebbles (self)
  class(pebble_t):: self
  !.............................................................................
END SUBROUTINE allocate_pebbles

!> =============================================================================
SUBROUTINE deallocate_pebbles (self)
  class(pebble_t):: self
  !.............................................................................
END SUBROUTINE deallocate_pebbles

!> =============================================================================
SUBROUTINE update_pebbles (self)
  class(pebble_t):: self
  !.............................................................................
  self%time = self%time + self%dtime
END SUBROUTINE update_pebbles

!> =============================================================================
SUBROUTINE print_pebbles (self)
  class(pebble_t):: self
  !.............................................................................
  print*, self%id, self%time
END SUBROUTINE print_pebbles

!> =============================================================================
END MODULE
