!> =============================================================================
!> $Id$
!> MPI point-to-point communication calls
!> 
!> To make use of this module, make sure it is compiled together with your code,
!> using either the Makefile and the 'make' command, or compiling manually, with
!> 
!> mpi = ../../../../../mpi                     # or wherever
!> mpifort -c $mpi/mpi_mod.f90                  # compiler mpi_mod.f90
!> mpifort -c $mpi/mpi_coords.f90               # add any other module your need
!> mpifort -c $mpi/mpi_send.f90                 # add any other module your need
!> mpifort -c your_code.f90                     # compile your code
!> mpifort *.o -o your_code.x                   # link together into your_code.x 
!>
!> In your_code.f90 you add lines such as:
!>
!> USE mpi_coords_mod                           ! makes the mpi_coords% object available
!> USE mpi_send_mod                             ! makes the mpi_send% object available
!> implicit none
!> real, dimension(:), allocatable: send        ! send buffer
!> real, dimension(:), allocatable: recv        ! receive buffer
!> integer:: rq(2)                              ! requests
!> integer:: nw                                 ! number of words in buffers
!> ...
!> call mpi_coords%init(dims=[mx,my,mz])        ! initializes the mpi_coords%object
!> call mpi_coords%print                        ! prints a summary of the MPI-geometry
!> ...
!> allocate(send(nw),recv(nw))                  ! allocate buffers
!> ...
!> call mpi_send%isend_reals(nw,send,to,rq(1))  ! send nw floating point values to MPI rank 'to'
!> call mpi_send%irecv_reals(nw,recv,fr,rq(2))  ! send nw floating point values from MPI rank 'fr'
!> call mpi_send%waitall(2,rq)                  ! wait for the send and recvs to complete
!> ...
!> call mpi_coords%end                          ! closes MPI
!>
!> Too see which variables and procedures are inside mpi_send%, just look below!
!>
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
SUBROUTINE bcast_ints (mp, nw, a, rank)
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
SUBROUTINE bcast_reals (mp, nw, a, rank)
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
SUBROUTINE bcast_real8s (mp, nw, a, rank)
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
SUBROUTINE isend_ints (mp, nw, a, rank, req)
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
SUBROUTINE irecv_ints (mp, nw, a, rank, req)
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
SUBROUTINE isend_reals (mp, nw, a, rank, req)
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
SUBROUTINE issend_reals (mp, nw, a, rank, req)
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
SUBROUTINE irecv_reals (mp, nw, a, rank, req)
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
SUBROUTINE waitall (mp, n, req)
  implicit none
  class(mpi_send_t):: mp
  integer:: n, req(n)
  call MPI_Waitall (n, req, MPI_STATUSES_IGNORE, mpi_err)
  call mp%assert ('waitall', mpi_err)
END SUBROUTINE waitall

END MODULE mpi_send_mod
