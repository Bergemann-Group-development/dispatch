!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for unigrid patches.  This is a temporary extension, whose
!> common properties in the end are copied over to the permanent patch.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE unigrid_mod
  USE mpi_coord_mod
  USE patch_mod
  USE io_mod
  implicit none
  type box_t
    real:: size(3)
    real:: num(3)
    integer:: dim(3)
  end type
  type, extends(patch_t):: unigrid_t
    integer:: mpi_dims(3)
    integer:: mpi_coords(3)
    real:: resolution(3)
    type(box_t):: box, rank_box, task_box
  contains
    procedure:: generate
  end type
  type(unigrid_t):: unigrid
PRIVATE
PUBLIC unigrid_t, unigrid, box_t
CONTAINS

SUBROUTINE generate (self)
  class(unigrid_t):: self
END SUBROUTINE

!> =============================================================================
!> The main extra task of the init procedure here is to figure out the exact
!> position and size of the patch
!> =============================================================================
SUBROUTINE init (self)
  class(unigrid_t):: self
  type(mpi_coord_t):: mpi
  integer:: tasks_per_dir(3), task_pos(3), ranks_per_dir(3), mpi_pos(3)
END SUBROUTINE init

SUBROUTINE init2 (self)
  class(unigrid_t):: self
  type(mpi_coord_t):: mpi
  integer:: tasks_per_dir(3), task_pos(3), ranks_per_dir(3), mpi_pos(3)
  !.............................................................................
END SUBROUTINE init2

END MODULE
