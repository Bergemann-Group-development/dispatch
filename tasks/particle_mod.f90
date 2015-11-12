!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE particle_mod
  USE task_mod
  implicit none
  !-----------------------------------------------------------------------------
  ! Particles.  The first dimension is over size bins, which is of suitable
  ! size for vectorization.  The position is single precision, with only the task
  ! bulk position in double precision
  !-----------------------------------------------------------------------------
  type particle_t
    real(4), dimension(:,:), pointer:: p, v     ! position, velocity
    real(4), dimension(:),   pointer:: m, s     ! mass, size
    real(4), dimension(:),   pointer:: t, dt    ! position, velocity
  end type
  type, extends(tasks_t):: particles_t
    type(particle_t), dimension(:,:), pointer:: particles
    contains
    procedure:: init
    procedure:: allocate
    procedure:: deallocate
    procedure:: update
    procedure:: print
    procedure:: check
    procedure:: send_to
    procedure:: timestep
  end type
!PRIVATE
!PUBLIC particles_t
CONTAINS

!> =============================================================================
SUBROUTINE init (self)
  class(particles_t):: self
  !.............................................................................
END SUBROUTINE

!> =============================================================================
SUBROUTINE check (self)
  class(particles_t):: self
  !.............................................................................
END SUBROUTINE

!> =============================================================================
SUBROUTINE timestep (self)
  class(particles_t):: self
  !.............................................................................
END SUBROUTINE

!> =============================================================================
SUBROUTINE send_to (self, rank)
  class(particles_t):: self
  integer:: rank
  !.............................................................................
END SUBROUTINE

!> =============================================================================
SUBROUTINE allocate (self)
  class(particles_t):: self
  !.............................................................................
END SUBROUTINE allocate

!> =============================================================================
SUBROUTINE deallocate (self)
  class(particles_t):: self
  !.............................................................................
END SUBROUTINE deallocate

!> =============================================================================
SUBROUTINE update (self)
  class(particles_t):: self
  !.............................................................................
  self%time = self%time + self%dtime
END SUBROUTINE update

!> =============================================================================
SUBROUTINE print (self)
  class(particles_t):: self
  !.............................................................................
  print*, self%id, self%time
END SUBROUTINE print

!> =============================================================================
END MODULE
