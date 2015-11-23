!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for planet task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE star_mod
  USE task_mod
  USE particle_mod
  implicit none
  !-----------------------------------------------------------------------------
  type, extends(particles_t):: star_t
    real:: mass, radius, surface_temperature
    contains
    procedure:: init       => init_star
    procedure:: allocate   => allocate_star
    procedure:: deallocate => deallocate_star
    procedure:: update     => update_star
    procedure:: print      => print_star
  end type
PRIVATE
PUBLIC star_t
CONTAINS

!> =============================================================================
SUBROUTINE init_star (self)
  class(star_t):: self
  !.............................................................................
  call self%task_t%init
  print*,'star is ID', self%id
END SUBROUTINE init_star

!> =============================================================================
SUBROUTINE allocate_star (self)
  class(star_t):: self
  !.............................................................................
END SUBROUTINE allocate_star

!> =============================================================================
SUBROUTINE deallocate_star (self)
  class(star_t):: self
  !.............................................................................
END SUBROUTINE deallocate_star

!> =============================================================================
SUBROUTINE update_star (self)
  class(star_t):: self
  !.............................................................................
  self%time = self%time + self%dtime
END SUBROUTINE update_star

!> =============================================================================
SUBROUTINE print_star (self)
  class(star_t):: self
  !.............................................................................
  print*, self%id, self%time
END SUBROUTINE print_star

!> =============================================================================
END MODULE
