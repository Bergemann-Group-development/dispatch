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
    real(4), dimension(:), pointer:: p, v     ! position, velocity
    real(4),               pointer:: m, s     ! mass, size
    real(4),               pointer:: t, dt    ! position, velocity
  end type
  type, abstract, extends(timeslots_t):: particles_t
    type(particle_t), dimension(:,:), pointer:: particles
  end type
END MODULE
