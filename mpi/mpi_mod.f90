!> =============================================================================
!> $Id$
!> =============================================================================
MODULE mpi_mod
  USE omp_mod
  USE mpi_basic
  USE mpi_coords
  USE mpi_all
  USE mpi_reduce
  USE mpi_rma
  USE mpi_reduce
  USE mpi_send
  USE mpi_timer
  implicit none
  type box_t
    real:: size(3)
    real:: num(3)
    integer:: dim(3)
  end type
  type(box_t):: box, rank, task
  type mpi_t
    integer:: size, rank, coord(3)
  contains
    procedure:: init
  end type
  type(mpi_t):: info
  !-----------------------------------------------------------------------------
  ! Template task, with procedure to place it it MPI-space
  !-----------------------------------------------------------------------------
  type template_t
    type(box_t):: box, mpi, task
    real(8):: position(3)
    integer:: itask
    integer:: rank
  contains
    procedure:: init => init_template
  end type
PRIVATE
PUBLIC mpi_t, template_t
CONTAINS

!> =============================================================================
!> Compute a task position, based on the linear index itask (0-indexed).
!>
!> The 3-D array template%mpi%num should be set to the desired MPI-geometry =
!> the number of MPI-ranks per direction.   The product(mpi%num) should thus
!> be equal to mpi_size
!>
!> The template%task%dim should be set to the desired task dimensions; i.e.,
!> the number of task points per direction (not counting guard cells).  The
!> box dimensions need to be divisible with the task dimensions, obviously.
!> =============================================================================
SUBROUTINE init_template (self)
  class(template_t):: self
  integer:: tasks_per_dir(3), task_pos(3), ranks_per_dir(3), mpi_pos(3)
  !-----------------------------------------------------------------------------
  ! The number of tasks per direction gives the physicals size from the box size
  !-----------------------------------------------------------------------------
  tasks_per_dir = self%box%dim/self%task%dim
  self%task%size = self%box%size/tasks_per_dir
  !-----------------------------------------------------------------------------
  ! The accumulated number of tasks per direction gives the task position, first
  ! as an integer 3-D index, then as a physical position
  !-----------------------------------------------------------------------------
  tasks_per_dir(2) = tasks_per_dir(2)*tasks_per_dir(1)
  tasks_per_dir(3) = tasks_per_dir(3)*tasks_per_dir(3)
  task_pos = self%itask/tasks_per_dir
  self%position = (task_pos+0.5d0)*self%task%size
  if (product(t%mpi%num) /= mpi_size) then
    print*,'rank%num, mpi_size', t%mpi%num, mpi_size
    call abort_mpi ('wrong size')
  end if
  !-----------------------------------------------------------------------------
  ! The number of ranks per direction gives the MPI-position and the mpi rank
  !-----------------------------------------------------------------------------
  ranks_per_dir = self%box%dim/self%mpi%dim
  !-----------------------------------------------------------------------------
  ranks_per_dir(2) = ranks_per_dir(2)*ranks_per_dir(1)
  ranks_per_dir(3) = ranks_per_dir(3)*ranks_per_dir(3)
  mpi_pos = self%itask/ranks_per_dir
  self%rank = product(mpi_pos+1)-1
  !
  if (all(t%task%num*t%task%dim /= t%mpi%dim)) then
    print'(3(3i8,2x))', t%task%num, t%task%dim, t%mpi%dim
    call abort_mpi ('wrong size')
  end if
END SUBROUTINE init_template

!> =============================================================================
!> Initialize an object with the most important pieces of information about MPI:
!> mpi%size = the number of MPI ranks
!> mpi%rank = the rank of the current MPI process (0-indexed)
!> mpi%coord = the 3-D position of the rank, if a Cartesian arrangement is used
!> =============================================================================
SUBROUTINE init (mpi)
  class(mpi_t):: mpi
  !.............................................................................
  mpi%size = mpi_size
  mpi%rank = mpi_rank
  mpi%coord = mpi_coord
END SUBROUTINE init

END MODULE mpi_mod
