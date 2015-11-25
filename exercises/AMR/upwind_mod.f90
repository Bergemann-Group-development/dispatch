!> =============================================================================
!> $Id$
!> =============================================================================
MODULE upwind_mod
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
  integer:: m, rank, rq(2), k, i
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dx
  real, dimension(self%n(2),self%n(3)):: send, recv, dx1, dx2
  real:: c, delta
  !.............................................................................
  if (first_time) then
    first_time = .false.
  end if
  m = self%n(1)
  c = 0.5/self%ds(1)

  delta = self%velocity(1)*self%dt*self%factor/self%ds(1)
  do i=2,m-1
    dx1 = f(i  ,:,:)-f(i-1,:,:)
    dx2 = f(i+1,:,:)-f(i  ,:,:)
    dx(i,:,:) = c*(dx1*(1.0+delta)+dx2*(1.0-delta))
  end do
  dx(2:m-1,:,:) = (f(3:m,:,:)-f(1:m-2,:,:))*c  ! this part is local

  dx1 = f(1,:,:)-recv
  dx2 = f(2,:,:)-f(1 ,:,:)
  dx(1,:,:) = c*(dx1*(1.0+delta)+dx2*(1.0-delta))

  dx1 = f(m,:,:)-f(m-1 ,:,:)
  dx2 = recv-f(m,:,:)
  dx(m,:,:) = c*(dx1*(1.0+delta)+dx2*(1.0-delta))
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
  dy(:,1    ,:) = (f(:,2  ,:)-recv        )*c  ! remote (:,n,:) slice needed!
  send = f(:,1,:)
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
  dz(:,:,1    ) = (f(:,:,2  )-recv        )*c  ! remote (:,:,n) slice needed!
  send = f(:,:,1)
  dz(:,:,  m  ) = (recv      -f(:,:,  m-1))*c  ! remote (:,:,1) slice needed!
END FUNCTION dz

END MODULE upwind_mod
