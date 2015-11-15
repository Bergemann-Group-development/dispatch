!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for dust tasks
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE dust_mod
  USE io_mod
  USE task_mod
  USE particle_mod
  implicit none
  !-----------------------------------------------------------------------------
  type, extends(particles_t):: dust_t
    real:: mass, radius, surface_temperature
    contains
    procedure:: init       => init_dust
    procedure:: allocate   => allocate_dust
    procedure:: deallocate => deallocate_dust
    procedure:: update     => update_dust
    procedure:: print      => print_dust
  end type
PRIVATE
PUBLIC dust_t
CONTAINS

!> =============================================================================
SUBROUTINE init_dust (self)
  class(dust_t):: self
  !.............................................................................
  call io%init
  self%verbose = io%verbose
  if (self%verbose > 0) print*,'init_dust:', self%id, self%time
END SUBROUTINE init_dust

!> =============================================================================
SUBROUTINE allocate_dust (self)
  class(dust_t):: self
  !.............................................................................
END SUBROUTINE allocate_dust

!> =============================================================================
SUBROUTINE deallocate_dust (self)
  class(dust_t):: self
  !.............................................................................
END SUBROUTINE deallocate_dust

!> =============================================================================
SUBROUTINE update_dust (self)
  class(dust_t):: self
  !.............................................................................
  self%time = self%time + self%dtime
END SUBROUTINE update_dust

!> =============================================================================
SUBROUTINE print_dust (self)
  class(dust_t):: self
  !.............................................................................
  print*, self%id, self%time
END SUBROUTINE print_dust

!> =============================================================================
END MODULE
