!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for planet task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE planet_mod
  USE task_mod
  USE particle_mod
  implicit none
  !-----------------------------------------------------------------------------
  type, extends(particles_t):: planet_t
    real:: mass, radius, surface_temperature
    contains
    procedure:: init       => init_planet
    procedure:: allocate   => allocate_planet
    procedure:: deallocate => deallocate_planet
    procedure:: update     => update_planet
    procedure:: print      => print_planet
  end type
!PRIVATE
!PUBLIC planet_t
CONTAINS

!> =============================================================================
SUBROUTINE init_planet (self)
  class(planet_t):: self
  !.............................................................................
END SUBROUTINE init_planet

!> =============================================================================
SUBROUTINE allocate_planet (self)
  class(planet_t):: self
  !.............................................................................
END SUBROUTINE allocate_planet

!> =============================================================================
SUBROUTINE deallocate_planet (self)
  class(planet_t):: self
  !.............................................................................
END SUBROUTINE deallocate_planet

!> =============================================================================
SUBROUTINE update_planet (self)
  class(planet_t):: self
  !.............................................................................
  self%time = self%time + self%dtime
END SUBROUTINE update_planet

!> =============================================================================
SUBROUTINE print_planet (self)
  class(planet_t):: self
  !.............................................................................
  print*, self%id, self%time
END SUBROUTINE print_planet

!> =============================================================================
END MODULE
