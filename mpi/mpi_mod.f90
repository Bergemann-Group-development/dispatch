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
  type mpi_info_t
    integer:: size, rank, coord(3)
  contains
    procedure:: init
  end type
  type(mpi_info_t):: info
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
CONTAINS

!> =============================================================================
!> Compute a task position, based on the linear index itask (0-indexed)
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
  !-----------------------------------------------------------------------------
  ! The number of ranks per direction gives the MPI-position and the mpi rank
  !-----------------------------------------------------------------------------
  ranks_per_dir = self%box%dim/self%mpi%dim
  !-----------------------------------------------------------------------------
  ranks_per_dir(2) = ranks_per_dir(2)*ranks_per_dir(1)
  ranks_per_dir(3) = ranks_per_dir(3)*ranks_per_dir(3)
  mpi_pos = self%itask/ranks_per_dir
  self%rank = product(mpi_pos+1)-1
END SUBROUTINE init_template

SUBROUTINE init (mpi)
  class(mpi_info_t):: mpi
  integer:: rank_num(3)=1, box_dim=64, patch_dim=16
  real:: box_size(3)=1.0
  type(template_t):: t
  namelist /mpi_params/ box_size, box_dim, rank_num, patch_dim
  !.............................................................................
  t%box%size  = box_size
  t%box%dim   = box_dim
  !
  t%mpi%num   = rank_num
  t%mpi%dim   = t%box%dim/t%mpi%num
  t%mpi%size  = t%box%size/t%mpi%num
  !
  t%task%dim  = patch_dim
  t%task%num  = t%mpi%dim/t%mpi%dim
  t%task%size = t%mpi%size/t%mpi%dim
  !
  if (product(t%mpi%num) /= mpi_size) then
    print*,'rank%num, mpi_size', t%mpi%num, mpi_size
    call abort_mpi ('wrong size')
  end if
  if (all(t%task%num*t%task%dim /= t%mpi%dim)) then
    print'(3(3i8,2x))', t%task%num, t%task%dim, t%mpi%dim
    call abort_mpi ('wrong size')
  end if
END SUBROUTINE init

FUNCTION get_mpi_info ()
  type(mpi_info_t):: get_mpi_info
  get_mpi_info%size  = mpi_size
  get_mpi_info%rank  = mpi_rank
  get_mpi_info%coord = mpi_coord
END FUNCTION

END MODULE mpi_mod
