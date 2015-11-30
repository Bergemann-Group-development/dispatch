!> =============================================================================
!> $Id$
!> =============================================================================
MODULE central_mod
  USE trace_mod
  USE derivatives_mod
  implicit none
  private
  type, public, extends(derivatives_t):: central_t
  contains
    procedure:: dx, dy, dz
  end type
  logical:: first_time=.true.
CONTAINS

!> =============================================================================
!> Derivatives evaluated with central differences.  When using neighbor lists,
!> guard zones are updated by the neighbor tasks
!> =============================================================================
FUNCTION dx (self, f)
  class(central_t):: self
  integer:: m, rank, rq(2), k
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dx
  real, dimension(self%n(2),self%n(3)):: send, recv
  real:: c
  !.............................................................................
  call trace_begin ('central_t%dx', 3)
  print*,'central', self%ds
  print*,'central', self%n
  if (first_time) then
    first_time = .false.
  end if
  m = self%n(1)
  c = 0.5/self%ds(1)
  print*,'central', c
  dx(2:m-1,:,:) = (f(3:m,:,:)-f(1:m-2,:,:))*c  ! this part is local
  recv = f(m,:,:)
  dx(1    ,:,:) = (f(2  ,:,:)-recv        )*c  ! remote (n,:,:) slice needed!
  recv = f(1,:,:)
  dx(  m  ,:,:) = (recv      -f(  m-1,:,:))*c  ! remote (1,:,:) slice needed!
  call trace_end
END FUNCTION dx
!> =============================================================================
FUNCTION dy (self, f)
  class(central_t):: self
  integer:: n(3), m, rq(2)
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dy
  real, dimension(self%n(1),self%n(3)):: send, recv
  real:: c
  !.............................................................................
  call trace_begin ('central_t%dy', 3)
  m = self%n(2)
  c = 0.5/self%ds(2)
  dy(:,2:m-1,:) = (f(:,3:m,:)-f(:,1:m-2,:))*c  ! this part is local
  recv = f(:,m,:)
  dy(:,1    ,:) = (f(:,2  ,:)-recv        )*c  ! remote (:,n,:) slice needed!
  recv = f(:,1,:)
  dy(:,  m  ,:) = (recv      -f(:,  m-1,:))*c  ! remote (:,1,:) slice needed!
  call trace_end
END FUNCTION dy
!> =============================================================================
FUNCTION dz (self, f)
  class(central_t):: self
  integer:: n(3), m, rq(2)
  real, dimension(self%n(1),self%n(2),self%n(3)):: f, dz
  real, dimension(self%n(1),self%n(2)):: send, recv
  real:: c
  !.............................................................................
  call trace_begin ('central_t%dz', 3)
  m = self%n(3)
  c = 0.5/self%ds(3)
  dz(:,:,2:m-1) = (f(:,:,3:m)-f(:,:,1:m-2))*c  ! this part is local
  recv = f(:,:,m)
  dz(:,:,1    ) = (f(:,:,2  )-recv        )*c  ! remote (:,:,n) slice needed!
  recv = f(:,:,1)
  dz(:,:,  m  ) = (recv      -f(:,:,  m-1))*c  ! remote (:,:,1) slice needed!
  call trace_end
END FUNCTION dz
!> =============================================================================
END MODULE central_mod

