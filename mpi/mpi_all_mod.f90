!> =============================================================================
!> $Id$
!> Module with all-to-all MPI calls 
!> =============================================================================
MODULE mpi_all_mod
  USE mpi_mod
  implicit none 
  include 'mpif.h'
  integer, private:: mpi_err
  type, extends(mpi_t):: mpi_all_t
  contains
    procedure:: scatter_reals
    procedure:: gather_reals
    procedure:: allgather_reals
    procedure:: alltoall_reals
    procedure:: alltoallv_reals
  end type
  type(mpi_all_t):: mpi_all
PRIVATE
PUBLIC mpi_all_t, mpi_all
CONTAINS

!> =============================================================================
SUBROUTINE scatter_reals (mp, send, recv, count, rank, nrank)
  implicit none
  class(mpi_all_t):: mp
  integer count, rank, nrank
  real, dimension(count, nrank):: send
  real, dimension(count):: recv
!...............................................................................
  call MPI_Scatter (send, count, MPI_REAL, recv, count, MPI_REAL, rank, mpi%comm, mpi_err)
  call mpi%assert ('allscatter', mpi_err)
END SUBROUTINE scatter_reals

!> =============================================================================
SUBROUTINE gather_reals (mp, send, recv, count, rank, nrank)
  implicit none
  class(mpi_all_t):: mp
  integer count, rank, nrank
  real, dimension(count):: send
  real, dimension(count, nrank):: recv
!...............................................................................
  call MPI_Gather (send, count, MPI_REAL, recv, count, MPI_REAL, rank, mpi%comm, mpi_err)
  call mpi%assert ('allscatter', mpi_err)
END SUBROUTINE gather_reals

!> =============================================================================
SUBROUTINE allgather_reals (mp, send, recv, count, rank, nrank)
  implicit none
  class(mpi_all_t):: mp
  integer count, rank, nrank
  real, dimension(count, nrank):: send, recv
!...............................................................................
  call MPI_Allgather (send, count, MPI_REAL, recv, count, MPI_REAL, mpi%comm, mpi_err)
  call mpi%assert ('allgather', mpi_err)
END SUBROUTINE allgather_reals

!> =============================================================================
SUBROUTINE alltoall_reals (mp, send, recv, count, nrank)
  implicit none
  class(mpi_all_t):: mp
  integer count, nrank
  real, dimension(count, nrank):: send, recv
!...............................................................................
  call MPI_Alltoall (send, count, MPI_REAL, recv, count, MPI_REAL, mpi%comm, mpi_err)
  call mpi%assert ('alltoall', mpi_err)
END SUBROUTINE alltoall_reals

!> =============================================================================
SUBROUTINE alltoallv_reals (mp, send, send_count, send_offset, &
                                    recv, recv_count, recv_offset, nrank)
  implicit none
  class(mpi_all_t):: mp
  integer mpi_type, nrank
  integer, dimension(nrank):: send_count, recv_count, send_offset, recv_offset
  real, dimension(:):: send, recv
!...............................................................................
  call MPI_Alltoallv (send, send_count, send_offset, MPI_REAL, &
                      recv, recv_count, recv_offset, MPI_REAL, mpi%comm, mpi_err)
  call mpi%assert ('alltoallv', mpi_err)
END SUBROUTINE alltoallv_reals

END MODULE mpi_all_mod
