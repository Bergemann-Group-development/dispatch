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
    procedure:: init       => init_part
    procedure:: allocate   => allocate_part
    procedure:: deallocate => deallocate_part
    procedure:: update     => update_part
    procedure:: print      => print_part
  end type
!PRIVATE
!PUBLIC particles_t
CONTAINS

!> =============================================================================
SUBROUTINE init_part (self, id)
  class(particles_t):: self
  integer:: id
  !.............................................................................
  self%id = id
END SUBROUTINE init_part

!> =============================================================================
SUBROUTINE allocate_part (self)
  class(particles_t):: self
  !.............................................................................
END SUBROUTINE allocate_part

!> =============================================================================
SUBROUTINE deallocate_part (self)
  class(particles_t):: self
  !.............................................................................
END SUBROUTINE deallocate_part

!> =============================================================================
SUBROUTINE update_part (self)
  class(particles_t):: self
  !.............................................................................
  self%time = self%time + self%dtime
END SUBROUTINE update_part

!> =============================================================================
SUBROUTINE print_part (self)
  class(particles_t):: self
  !.............................................................................
  print*, self%id, self%time
END SUBROUTINE print_part

!> =============================================================================
END MODULE
