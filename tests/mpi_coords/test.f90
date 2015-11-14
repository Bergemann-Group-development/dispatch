!> ============================================================================
!> $Id$
!> ============================================================================
PROGRAM mpi_coords_test
  USE mpi_basic
  USE mpi_coords
  implicit none
  type(mpi_coord_t):: mpi
  integer, dimension(3):: dims, mpi_dims
  integer:: i
  !.............................................................................
  call mpi%mpi_t%init
  mpi_dims  = [2,2,1]
  call mpi%init (mpi_dims)
  call mpi%print
  mpi_dims  = [1,3,1]
  call mpi%init (mpi_dims)
  call mpi%print
  dims  = [40,40,27]
  call mpi%init (dims=dims)
  call mpi%print
  call mpi%end
  !.............................................................................
  CONTAINS
END PROGRAM
