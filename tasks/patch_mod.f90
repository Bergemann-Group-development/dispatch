!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for patches, which adds pointers to memory and mesh, and
!> number of dimensions and variables.  The number of dimensions is hardwired
!> to 3, since we can just choose the actual dimensions to get 2-D and 1-D.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE patch_mod
  USE mpi_mod
  USE item_mod
  USE task_mod
  USE mesh_mod
  USE io_mod
  implicit none
  integer, parameter:: ndim = 3
  type, extends(tasks_t):: patch_t
    real, dimension(:,:,:,:,:,:), pointer:: mem
    type(mesh_t), dimension(:), pointer:: m
    integer:: ndim = ndim           ! Number of dimensions
    integer:: nvar                  ! Number of variables
  contains
  end type
  !-----------------------------------------------------------------------------
  ! Template task, with procedure to place it it MPI-space
  !-----------------------------------------------------------------------------
  type box_t
    real:: size(3)
    real:: num(3)
    integer:: dim(3)
  end type
  type template_t
    type(box_t):: box, mpi, task
    real(8):: position(3)
    integer:: itask
    integer:: rank
  contains
    procedure:: init => init_template
  end type
PRIVATE
PUBLIC patch_t, template_t
CONTAINS

!> =============================================================================
!> Compute a task position, based on the linear index itask (0-indexed).
!>
!> The 3-D array template%mpi%num should be set to the desired MPI-geometry =
!> the number of MPI-ranks per direction.   The product(mpi%num) should thus
!> be equal to mpi%size
!>
!> The template%task%dim should be set to the desired task dimensions; i.e.,
!> the number of task points per direction (not counting guard cells).  The
!> box dimensions need to be divisible with the task dimensions, obviously.
!> =============================================================================
SUBROUTINE init_template (self)
  class(template_t):: self
  type(mpi_t):: mpi
  integer:: tasks_per_dir(3), task_pos(3), ranks_per_dir(3), mpi_pos(3)
  !-----------------------------------------------------------------------------
  ! The number of tasks per direction gives the physicals size from the box size
  !-----------------------------------------------------------------------------
  call mpi%init
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
  if (product(self%mpi%num) /= mpi%size) then
    print*,'rank%num, mpi%size', self%mpi%num, mpi%size
    call mpi%abort ('wrong size')
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
  if (all(self%task%num*self%task%dim /= self%mpi%dim)) then
    print'(3(3i8,2x))', self%task%num, self%task%dim, self%mpi%dim
    call mpi%abort ('wrong size')
  end if
END SUBROUTINE init_template

END MODULE
