!> =============================================================================
!> $Id$
!> =============================================================================
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

!> =============================================================================
!> Derivatives evaluated with central differences.  When using MPI domain
!> decomposition, one needs to get a single slice of "ghost cells" on each of
!> the 6 faces of the local domain.
!> =============================================================================
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
  c = 0.5/self%ds(1)

  !***************
  ! HINT1: To make a correct UPWIND operator you have to check
  ! the sign of self%velocity(1) and create the corresponding operator
  !
  ! HINT2: Since the upwind only uses information to one side of the
  ! grid-point, you only need one of the 
  !    send = ...
  !    call mpi_send%isend...
  !    call mpi_send%irecv...
  !    call mpi_send%waitall
  !    dx(...
  ! blocks below in EACH case.
  !
  ! HINT3: Take care with array indices. Sketch on a piece of paper first.
  ! They are slightly different for Upwind and central, and easy to get wrong.
  ! Specially iff you copy-paste from x to yo to z. What is the size of "c"?
  !***************

  dx(2:m-1,:,:) = (f(3:m,:,:)-f(1:m-2,:,:))*c  ! this part is local

  send = f(m,:,:)
  call mpi_send%isend_reals (self%n(2)*self%n(3), send, mpi_coords%up(1), rq(1))
  call mpi_send%irecv_reals (self%n(2)*self%n(3), recv, mpi_coords%dn(1), rq(2))
  call mpi_send%waitall (2, rq)
  dx(1    ,:,:) = (f(2  ,:,:)-recv        )*c  ! remote (n,:,:) slice needed!

  send = f(1,:,:)
  call mpi_send%isend_reals (self%n(2)*self%n(3), send, mpi_coords%dn(1), rq(1))
  call mpi_send%irecv_reals (self%n(2)*self%n(3), recv, mpi_coords%up(1), rq(2))
  call mpi_send%waitall (2, rq)
  dx(  m  ,:,:) = (recv      -f(  m-1,:,:))*c  ! remote (1,:,:) slice needed!
END FUNCTION dx
!> =============================================================================
FUNCTION dy (self, f)
  class(upwind_t):: self
  integer:: n(3), m, rq(2)
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dy
  real, dimension(self%n(1),self%n(3)):: send, recv
  real:: c
  !.............................................................................
  m = self%n(2)
  c = 0.5/self%ds(2)
  dy(:,2:m-1,:) = (f(:,3:m,:)-f(:,1:m-2,:))*c  ! this part is local
  send = f(:,m,:)
  call mpi_send%isend_reals (self%n(1)*self%n(3), send, mpi_coords%up(2), rq(1))
  call mpi_send%irecv_reals (self%n(1)*self%n(3), recv, mpi_coords%dn(2), rq(2))
  call mpi_send%waitall (2, rq)
  dy(:,1    ,:) = (f(:,2  ,:)-recv        )*c  ! remote (:,n,:) slice needed!
  send = f(:,1,:)
  call mpi_send%isend_reals (self%n(1)*self%n(3), send, mpi_coords%dn(2), rq(1))
  call mpi_send%irecv_reals (self%n(1)*self%n(3), recv, mpi_coords%up(2), rq(2))
  call mpi_send%waitall (2, rq)
  dy(:,  m  ,:) = (recv      -f(:,  m-1,:))*c  ! remote (:,1,:) slice needed!
END FUNCTION dy
!> =============================================================================
FUNCTION dz (self, f)
  class(upwind_t):: self
  integer:: n(3), m, rq(2)
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dz
  real, dimension(self%n(1),self%n(2)):: send, recv
  real:: c
  !.............................................................................
  m = self%n(3)
  c = 0.5/self%ds(3)
  dz(:,:,2:m-1) = (f(:,:,3:m)-f(:,:,1:m-2))*c  ! this part is local
  send = f(:,:,m)
  call mpi_send%isend_reals (self%n(1)*self%n(2), send, mpi_coords%up(3), rq(1))
  call mpi_send%irecv_reals (self%n(1)*self%n(2), recv, mpi_coords%dn(3), rq(2))
  call mpi_send%waitall (2, rq)
  dz(:,:,1    ) = (f(:,:,2  )-recv        )*c  ! remote (:,:,n) slice needed!
  send = f(:,:,1)
  call mpi_send%isend_reals (self%n(1)*self%n(2), send, mpi_coords%dn(3), rq(1))
  call mpi_send%irecv_reals (self%n(1)*self%n(2), recv, mpi_coords%up(3), rq(2))
  call mpi_send%waitall (2, rq)
  dz(:,:,  m  ) = (recv      -f(:,:,  m-1))*c  ! remote (:,:,1) slice needed!
END FUNCTION dz

END MODULE upwind_mod
