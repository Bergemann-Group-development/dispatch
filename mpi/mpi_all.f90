!> =============================================================================
!> $Id$
!> Module with all-to-all MPI calls 
!> =============================================================================
MODULE mpi_all
  USE mpi_basic
  implicit none 
PUBLIC
CONTAINS

!> =============================================================================
SUBROUTINE scatter_reals_mpi (send, recv, count, rank, nrank, comm)
  implicit none
  integer count, rank, nrank, comm
  real, dimension(count, nrank):: send
  real, dimension(count):: recv
!...............................................................................
  call MPI_Scatter (send, count, MPI_REAL, recv, count, MPI_REAL, rank, comm, mpi_err)
  call assert_mpi ('allscatter_mpi', mpi_err)
END SUBROUTINE scatter_reals_mpi

!> =============================================================================
SUBROUTINE gather_reals_mpi (send, recv, count, rank, nrank, comm)
  implicit none
  integer count, rank, nrank, comm
  real, dimension(count):: send
  real, dimension(count, nrank):: recv
!...............................................................................
  call MPI_Gather (send, count, MPI_REAL, recv, count, MPI_REAL, rank, comm, mpi_err)
  call assert_mpi ('allscatter_mpi', mpi_err)
END SUBROUTINE gather_reals_mpi

!> =============================================================================
SUBROUTINE allgather_reals_mpi (send, recv, count, rank, nrank, comm)
  implicit none
  integer count, rank, nrank, comm
  real, dimension(count, nrank):: send, recv
!...............................................................................
  call MPI_Allgather (send, count, MPI_REAL, recv, count, MPI_REAL, comm, mpi_err)
  call assert_mpi ('allgather_mpi', mpi_err)
END SUBROUTINE allgather_reals_mpi

!> =============================================================================
SUBROUTINE alltoall_reals_mpi (send, recv, count, nrank, comm)
  implicit none
  integer count, nrank, comm
  real, dimension(count, nrank):: send, recv
!...............................................................................
  call MPI_Alltoall (send, count, MPI_REAL, recv, count, MPI_REAL, comm, mpi_err)
  call assert_mpi ('alltoall_mpi', mpi_err)
END SUBROUTINE alltoall_reals_mpi

!> =============================================================================
SUBROUTINE alltoallv_reals_mpi (send, send_count, send_offset, &
                          recv, recv_count, recv_offset, nrank, comm)
  implicit none
  integer mpi_type, nrank, comm
  integer, dimension(nrank):: send_count, recv_count, send_offset, recv_offset
  real, dimension(:):: send, recv
!...............................................................................
  call MPI_Alltoallv (send, send_count, send_offset, MPI_REAL, &
                      recv, recv_count, recv_offset, MPI_REAL, comm, mpi_err)
  call assert_mpi ('alltoallv_mpi', mpi_err)
END SUBROUTINE alltoallv_reals_mpi

END MODULE mpi_all
