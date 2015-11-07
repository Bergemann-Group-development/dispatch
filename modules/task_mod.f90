!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE task_mod
  implicit none
  integer, parameter:: ntime=7, nvect=16
  !-----------------------------------------------------------------------------
  ! Basic task info
  !-----------------------------------------------------------------------------
  type task_t
    integer:: id = 0
    integer:: rank = 0
    real(8):: time = 0d0
    real(8):: dtime = 0d0
    real(8):: position(3) = 0d0
    real(8):: velocity(3) = 0d0
    real(8):: size(3) = 0d0
    integer:: status = 0
    contains
    procedure init, allocate, deallocate, update, print
  end type
  !-----------------------------------------------------------------------------
  ! Timeslot info added
  !-----------------------------------------------------------------------------
  type, extends(task_t):: timeslots_t
    real(8):: t(ntime)
    real(8):: dt(ntime)
    integer:: it = 0
  end type
  !-----------------------------------------------------------------------------
  ! Mesh and memory for HD added
  !-----------------------------------------------------------------------------
  type mesh_t
    real(8):: p, s
    real(4), dimension(:), pointer:: r
    integer:: lb, ub, lbi, ubi, lbo, ubo
  end type
  type, extends(timeslots_t):: hd_t
    real, dimension(:,:,:,:,:,:), pointer:: mem
    type(mesh_t), dimension(:), pointer:: m
    real:: nu(4)
  end type
  !-----------------------------------------------------------------------------
  ! Additional MHD parameters added
  !-----------------------------------------------------------------------------
  type, extends(hd_t):: mhd_t
    real:: eta(4)
  end type
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
  type, extends(timeslots_t):: particles_t
    type(particle_t), dimension(:), pointer:: particles
  end type
  !-----------------------------------------------------------------------------
  type, extends(timeslots_t):: star_t
    real:: mass, radius, surface_temperature
  end type
  !-----------------------------------------------------------------------------
  type, extends(timeslots_t):: planet_t
    real:: mass, radius, surface_temperature
  end type
PRIVATE
PUBLIC task_t, hd_t, mhd_t, particles_t, planet_t, star_t
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
