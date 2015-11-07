! $Id: 0e5bfbb201ad3dac127445435aa8d590b4a50a74 $
!*******************************************************************************
MODULE mpi_send
  USE mpi_base
  implicit none 
PUBLIC
CONTAINS

!*******************************************************************************
SUBROUTINE bcast_ints_mpi (a, nw, rank)
  implicit none
  integer:: nw, rank
  real:: a(nw)
  call MPI_Bcast (a, nw, MPI_INT, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_ints_mpi', mpi_err)
END SUBROUTINE bcast_ints_mpi

!*******************************************************************************
SUBROUTINE bcast_reals_mpi (a, nw, rank)
  implicit none
  integer:: nw, rank
  real:: a(nw)
  call MPI_Bcast (a, nw, MPI_REAL, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_reals_mpi', mpi_err)
END SUBROUTINE bcast_reals_mpi

!*******************************************************************************
SUBROUTINE bcast_real8_mpi (a, rank)
  implicit none
  integer:: rank
  real(8):: a
  call MPI_Bcast (a, 1, MPI_REAL8, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_reals_mpi', mpi_err)
END SUBROUTINE bcast_real8_mpi

!*******************************************************************************
SUBROUTINE bcast_real8s_mpi (a, nw, rank)
  implicit none
  integer:: nw, rank
  real(8):: a(nw)
  call MPI_Bcast (a, nw, MPI_REAL8, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_reals_mpi', mpi_err)
END SUBROUTINE bcast_real8s_mpi

!*******************************************************************************
SUBROUTINE isend_ints_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  integer:: a(nw)
  call MPI_Isend (a, nw, MPI_INT, rank, mpi_rank, mpi_comm, req, mpi_err)
  call assert_mpi ('isend_reals_mpi', mpi_err)
END SUBROUTINE isend_ints_mpi

!*******************************************************************************
SUBROUTINE irecv_ints_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  integer:: a(nw)
  call MPI_Irecv (a, nw, MPI_INT, rank, rank, mpi_comm, req, mpi_err)
  call assert_mpi ('irecv_reals_mpi', mpi_err)
END SUBROUTINE irecv_ints_mpi

!*******************************************************************************
SUBROUTINE isend_reals_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Isend (a, nw, MPI_REAL, rank, mpi_rank, mpi_comm, req, mpi_err)
  call assert_mpi ('isend_reals_mpi', mpi_err)
END SUBROUTINE isend_reals_mpi

!*******************************************************************************
SUBROUTINE issend_reals_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Issend (a, nw, MPI_REAL, rank, mpi_rank, mpi_comm, req, mpi_err)
  call assert_mpi ('isend_reals_mpi', mpi_err)
END SUBROUTINE issend_reals_mpi

!*******************************************************************************
SUBROUTINE irecv_reals_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Irecv (a, nw, MPI_REAL, rank, rank, mpi_comm, req, mpi_err)
  call assert_mpi ('irecv_reals_mpi', mpi_err)
END SUBROUTINE irecv_reals_mpi

!*******************************************************************************
SUBROUTINE wait_mpi (req)
  implicit none
  integer:: req
  integer:: status(MPI_STATUS_SIZE)
  call MPI_Wait (req, status, mpi_err)
  call assert_mpi ('wait_mpi', mpi_err)
END SUBROUTINE wait_mpi

!*******************************************************************************
SUBROUTINE waitall_mpi (req, n)
  implicit none
  integer:: n, req(n)
  call MPI_Waitall (n, req, MPI_STATUSES_IGNORE, mpi_err)
  call assert_mpi ('waitall_mpi', mpi_err)
END SUBROUTINE waitall_mpi

END MODULE mpi_send
