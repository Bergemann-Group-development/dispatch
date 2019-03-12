!===============================================================================
!> $Id$
!===============================================================================
MODULE upwind_mod
  USE mpi_send_mod
  USE mpi_coords_mod
  USE derivatives_mod
  implicit none
  private
  type, public, extends(derivatives_t):: upwind_t
    real, pointer :: velocity(:)
  contains
    procedure:: dx, dy, dz
  end type
  logical:: first_time=.true.
CONTAINS

!===============================================================================
!> Derivatives evaluated with central differences.  When using MPI domain
!> decomposition, one needs to get a single slice of "ghost cells" on each of
!> the 6 faces of the local domain.
!===============================================================================
FUNCTION dx (self, f)
  class(upwind_t):: self
  integer:: m, rank, rq(2), k
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dx
  real, dimension(self%n(2),self%n(3)):: send, recv
  real:: c
  !.............................................................................
  if (first_time) then
    first_time = .false.
    if (mpi_coords%rank==0) print *,'using upwind differences'
  end if
  m = self%n(1)
  c = 1./self%ds(1)

  ! velocity is positive, we need stuff from "downstream"
  if (self%velocity(1) > 0.0) then
    dx(2:m,  :,:) = (f(2:m,:,:)-f(1:m-1,:,:))*c  ! this part is local

    send = f(m,:,:)
    call mpi_send%isend_reals (self%n(2)*self%n(3), send, mpi_coords%up(1), rq(1))
    call mpi_send%irecv_reals (self%n(2)*self%n(3), recv, mpi_coords%dn(1), rq(2))
    call mpi_send%waitall (2, rq)
    dx(1    ,:,:) = (f(1  ,:,:)-recv)*c          ! remote (n,:,:) slice needed!

  ! velocity is negative, we need stuff from "upstream"
  else
    dx(1:m-1,:,:) = (f(2:m,:,:)-f(1:m-1,:,:))*c  ! this part is local

    send = f(1,:,:)
    call mpi_send%isend_reals (self%n(2)*self%n(3), send, mpi_coords%dn(1), rq(1))
    call mpi_send%irecv_reals (self%n(2)*self%n(3), recv, mpi_coords%up(1), rq(2))
    call mpi_send%waitall (2, rq)
    dx(  m  ,:,:) = (recv-f(  m,:,:))*c          ! remote (1,:,:) slice needed!
  endif
END FUNCTION dx
!===============================================================================
FUNCTION dy (self, f)
  class(upwind_t):: self
  integer:: n(3), m, rq(2)
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dy
  real, dimension(self%n(1),self%n(3)):: send, recv
  real:: c
  !.............................................................................
  m = self%n(2)
  c = 1./self%ds(2)

  ! velocity is positive, we need stuff from "downstream"
  if (self%velocity(2) > 0.0) then
    dy(:,2:m,:) = (f(:,2:m,:)-f(:,1:m-1,:))*c   ! this part is local

    send = f(:,m,:)
    call mpi_send%isend_reals (self%n(1)*self%n(3), send, mpi_coords%up(2), rq(1))
    call mpi_send%irecv_reals (self%n(1)*self%n(3), recv, mpi_coords%dn(2), rq(2))
    call mpi_send%waitall (2, rq)
    dy(:,1  ,:) = (f(:,1  ,:)-recv)*c           ! remote (:,n,:) slice needed!

  ! velocity is negative, we need stuff from "upstream"
  else
    dy(:,1:m-1,:) = (f(:,2:m,:)-f(:,1:m-1,:))*c ! this part is local

    send = f(:,1,:)
    call mpi_send%isend_reals (self%n(1)*self%n(3), send, mpi_coords%dn(2), rq(1))
    call mpi_send%irecv_reals (self%n(1)*self%n(3), recv, mpi_coords%up(2), rq(2))
    call mpi_send%waitall (2, rq)
    dy(:,  m  ,:) = (recv-f(:,  m,:))*c         ! remote (:,1,:) slice needed!
  endif
END FUNCTION dy
!===============================================================================
FUNCTION dz (self, f)
  class(upwind_t):: self
  integer:: n(3), m, rq(2)
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dz
  real, dimension(self%n(1),self%n(2)):: send, recv
  real:: c
  !.............................................................................
  m = self%n(3)
  c = 1./self%ds(3)

  ! velocity is positive, we need stuff from "downstream"
  if (self%velocity(2) > 0.0) then
    dz(:,:,2:m) = (f(:,:,2:m)-f(:,:,1:m-1))*c   ! this part is local

    send = f(:,:,m)
    call mpi_send%isend_reals (self%n(1)*self%n(2), send, mpi_coords%up(3), rq(1))
    call mpi_send%irecv_reals (self%n(1)*self%n(2), recv, mpi_coords%dn(3), rq(2))
    call mpi_send%waitall (2, rq)
    dz(:,:,1  ) = (f(:,:,1  )-recv)*c           ! remote (:,:,n) slice needed!

  ! velocity is negative, we need stuff from "upstream"
  else
    dz(:,:,1:m-1) = (f(:,:,2:m)-f(:,:,1:m-1))*c ! this part is local

    send = f(:,:,1)
    call mpi_send%isend_reals (self%n(1)*self%n(2), send, mpi_coords%dn(3), rq(1))
    call mpi_send%irecv_reals (self%n(1)*self%n(2), recv, mpi_coords%up(3), rq(2))
    call mpi_send%waitall (2, rq)
    dz(:,:,  m  ) = (recv-f(:,:,  m))*c         ! remote (:,:,1) slice needed!
  endif
END FUNCTION dz

END MODULE upwind_mod
