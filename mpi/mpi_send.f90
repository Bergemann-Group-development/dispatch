!> =============================================================================
!> $Id$
!> MPI point-to-point communication calls
!> =============================================================================
MODULE mpi_send
  USE mpi_basic
  implicit none 
PUBLIC
CONTAINS

!> =============================================================================
!> Broadcase nw ints from rank
!> =============================================================================
SUBROUTINE bcast_ints_mpi (a, nw, rank)
  implicit none
  integer:: nw, rank
  real:: a(nw)
  call MPI_Bcast (a, nw, MPI_INT, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_ints_mpi', mpi_err)
END SUBROUTINE bcast_ints_mpi

!> =============================================================================
!> Broadcase nw reals from rank
!> =============================================================================
SUBROUTINE bcast_reals_mpi (a, nw, rank)
  implicit none
  integer:: nw, rank
  real:: a(nw)
  call MPI_Bcast (a, nw, MPI_REAL, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_reals_mpi', mpi_err)
END SUBROUTINE bcast_reals_mpi

!> =============================================================================
!> Broadcase one real8 from rank
!> =============================================================================
SUBROUTINE bcast_real8_mpi (a, rank)
  implicit none
  integer:: rank
  real(8):: a
  call MPI_Bcast (a, 1, MPI_REAL8, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_reals_mpi', mpi_err)
END SUBROUTINE bcast_real8_mpi

!> =============================================================================
!> Broadcase nw real8s from rank
!> =============================================================================
SUBROUTINE bcast_real8s_mpi (a, nw, rank)
  implicit none
  integer:: nw, rank
  real(8):: a(nw)
  call MPI_Bcast (a, nw, MPI_REAL8, rank, mpi_comm, mpi_err)
  call assert_mpi ('bcast_reals_mpi', mpi_err)
END SUBROUTINE bcast_real8s_mpi

!> =============================================================================
!> Non-blocking send of nw ints to rank
!> =============================================================================
SUBROUTINE isend_ints_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  integer:: a(nw)
  call MPI_Isend (a, nw, MPI_INT, rank, mpi_rank, mpi_comm, req, mpi_err)
  call assert_mpi ('isend_reals_mpi', mpi_err)
END SUBROUTINE isend_ints_mpi

!> =============================================================================
!> Non-blocking send of nw ints to rank
!> =============================================================================
SUBROUTINE irecv_ints_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  integer:: a(nw)
  call MPI_Irecv (a, nw, MPI_INT, rank, rank, mpi_comm, req, mpi_err)
  call assert_mpi ('irecv_reals_mpi', mpi_err)
END SUBROUTINE irecv_ints_mpi

!> =============================================================================
!> Non-blocking send of nw reals to rank
!> =============================================================================
SUBROUTINE isend_reals_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Isend (a, nw, MPI_REAL, rank, mpi_rank, mpi_comm, req, mpi_err)
  call assert_mpi ('isend_reals_mpi', mpi_err)
END SUBROUTINE isend_reals_mpi

!> =============================================================================
!> Non-blocking send of nw reals to rank
!> =============================================================================
SUBROUTINE issend_reals_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Issend (a, nw, MPI_REAL, rank, mpi_rank, mpi_comm, req, mpi_err)
  call assert_mpi ('isend_reals_mpi', mpi_err)
END SUBROUTINE issend_reals_mpi

!> =============================================================================
!> Non-blocking recv of nw reals from rank
!> =============================================================================
SUBROUTINE irecv_reals_mpi (a, nw, rank, req)
  implicit none
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Irecv (a, nw, MPI_REAL, rank, rank, mpi_comm, req, mpi_err)
  call assert_mpi ('irecv_reals_mpi', mpi_err)
END SUBROUTINE irecv_reals_mpi

!> =============================================================================
!> Wait for the completion of an MPI operation
!> =============================================================================
SUBROUTINE wait_mpi (req)
  implicit none
  integer:: req
  integer:: status(MPI_STATUS_SIZE)
  call MPI_Wait (req, status, mpi_err)
  call assert_mpi ('wait_mpi', mpi_err)
END SUBROUTINE wait_mpi

!> =============================================================================
!> Wait for the completion of all MPI operations corresponding to req(1:n)
!> =============================================================================
SUBROUTINE waitall_mpi (req, n)
  implicit none
  integer:: n, req(n)
  call MPI_Waitall (n, req, MPI_STATUSES_IGNORE, mpi_err)
  call assert_mpi ('waitall_mpi', mpi_err)
END SUBROUTINE waitall_mpi

END MODULE mpi_send
