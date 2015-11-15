!> =============================================================================
!> $Id$
!> MPI point-to-point communication calls
!> =============================================================================
MODULE mpi_send_mod
  USE mpi_mod
  implicit none 
  include "mpif.h"
  integer, private:: mpi_err
  type, extends(mpi_t):: mpi_send_t
  contains
    procedure:: bcast_ints
    procedure:: bcast_reals
    procedure:: bcast_real8
    procedure:: bcast_real8s
    procedure:: isend_ints
    procedure:: irecv_ints
    procedure:: isend_reals
    procedure:: issend_reals
    procedure:: irecv_reals
    procedure:: wait
    procedure:: waitall
  end type
  type(mpi_send_t):: mpi_send
PRIVATE
PUBLIC mpi_send_t, mpi_send
CONTAINS

!> =============================================================================
!> Broadcase nw ints from rank
!> =============================================================================
SUBROUTINE bcast_ints (mp, a, nw, rank)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank
  real:: a(nw)
  call MPI_Bcast (a, nw, MPI_INT, rank, mp%comm, mpi_err)
  call mp%assert ('bcast_ints', mpi_err)
END SUBROUTINE bcast_ints

!> =============================================================================
!> Broadcase nw reals from rank
!> =============================================================================
SUBROUTINE bcast_reals (mp, a, nw, rank)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank
  real:: a(nw)
  call MPI_Bcast (a, nw, MPI_REAL, rank, mp%comm, mpi_err)
  call mp%assert ('bcast_reals', mpi_err)
END SUBROUTINE bcast_reals

!> =============================================================================
!> Broadcase one real8 from rank
!> =============================================================================
SUBROUTINE bcast_real8 (mp, a, rank)
  implicit none
  class(mpi_send_t):: mp
  integer:: rank
  real(8):: a
  call MPI_Bcast (a, 1, MPI_REAL8, rank, mp%comm, mpi_err)
  call mp%assert ('bcast_reals', mpi_err)
END SUBROUTINE bcast_real8

!> =============================================================================
!> Broadcase nw real8s from rank
!> =============================================================================
SUBROUTINE bcast_real8s (mp, a, nw, rank)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank
  real(8):: a(nw)
  call MPI_Bcast (a, nw, MPI_REAL8, rank, mp%comm, mpi_err)
  call mp%assert ('bcast_reals', mpi_err)
END SUBROUTINE bcast_real8s

!> =============================================================================
!> Non-blocking send of nw ints to rank
!> =============================================================================
SUBROUTINE isend_ints (mp, a, nw, rank, req)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank, req
  integer:: a(nw)
  call MPI_Isend (a, nw, MPI_INT, rank, mp%rank, mp%comm, req, mpi_err)
  call mp%assert ('isend_reals', mpi_err)
END SUBROUTINE isend_ints

!> =============================================================================
!> Non-blocking send of nw ints to rank
!> =============================================================================
SUBROUTINE irecv_ints (mp, a, nw, rank, req)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank, req
  integer:: a(nw)
  call MPI_Irecv (a, nw, MPI_INT, rank, rank, mp%comm, req, mpi_err)
  call mp%assert ('irecv_reals', mpi_err)
END SUBROUTINE irecv_ints

!> =============================================================================
!> Non-blocking send of nw reals to rank
!> =============================================================================
SUBROUTINE isend_reals (mp, a, nw, rank, req)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Isend (a, nw, MPI_REAL, rank, mp%rank, mp%comm, req, mpi_err)
  call mp%assert ('isend_reals', mpi_err)
END SUBROUTINE isend_reals

!> =============================================================================
!> Non-blocking send of nw reals to rank
!> =============================================================================
SUBROUTINE issend_reals (mp, a, nw, rank, req)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Issend (a, nw, MPI_REAL, rank, mp%rank, mp%comm, req, mpi_err)
  call mp%assert ('isend_reals', mpi_err)
END SUBROUTINE issend_reals

!> =============================================================================
!> Non-blocking recv of nw reals from rank
!> =============================================================================
SUBROUTINE irecv_reals (mp, a, nw, rank, req)
  implicit none
  class(mpi_send_t):: mp
  integer:: nw, rank, req
  real:: a(nw)
  call MPI_Irecv (a, nw, MPI_REAL, rank, rank, mp%comm, req, mpi_err)
  call mp%assert ('irecv_reals', mpi_err)
END SUBROUTINE irecv_reals

!> =============================================================================
!> Wait for the completion of an MPI operation
!> =============================================================================
SUBROUTINE wait (mp, req)
  implicit none
  class(mpi_send_t):: mp
  integer:: req
  integer:: status(MPI_STATUS_SIZE)
  call MPI_Wait (req, status, mpi_err)
  call mp%assert ('wait', mpi_err)
END SUBROUTINE wait

!> =============================================================================
!> Wait for the completion of all MPI operations corresponding to req(1:n)
!> =============================================================================
SUBROUTINE waitall (mp, req, n)
  implicit none
  class(mpi_send_t):: mp
  integer:: n, req(n)
  call MPI_Waitall (n, req, MPI_STATUSES_IGNORE, mpi_err)
  call mp%assert ('waitall', mpi_err)
END SUBROUTINE waitall

END MODULE mpi_send_mod
