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
  USE mpi_coord_mod
  implicit none
  integer, dimension(3):: dims, mpi_dims
  !.............................................................................
  mpi_dims  = [2,2,1]
  call mpi_coord%init (mpi_dims)
  call mpi_coord%print
  !.............................................................................
  mpi_dims  = [1,3,1]
  call mpi_coord%init (mpi_dims)
  call mpi_coord%print
  !.............................................................................
  dims  = [40,40,27]
  call mpi_coord%init (dims=dims)
  call mpi_coord%print
  !.............................................................................
  call mpi_coord%end
END PROGRAM
