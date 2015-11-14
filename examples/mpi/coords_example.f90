!> ============================================================================
!> $Id$
!> Example of how to use the MPI_coords module to choose an MPI-geometry, and
!> print a summary of it.  Try it with, for example
!>   mpirun -np 2 test.x
!>   mpirun -np 3 test.x
!>   mpirun -np 4 test.x
!>   mpirun -np 9 test.x
!> ============================================================================
PROGRAM mpi_coords_test
  USE mpi_coords
  implicit none
  type(mpi_coord_t):: mpi
  integer, dimension(3):: dims, mpi_dims
  !.............................................................................
  mpi_dims  = [2,2,1]
  call mpi%init (mpi_dims)
  call mpi%print
  !.............................................................................
  mpi_dims  = [1,3,1]
  call mpi%init (mpi_dims)
  call mpi%print
  !.............................................................................
  dims  = [40,40,27]
  call mpi%init (dims=dims)
  call mpi%print
  !.............................................................................
  call mpi%end
END PROGRAM
