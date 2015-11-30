!> =============================================================================
!> $Id$
!> =============================================================================
MODULE central_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: central_t
    real:: ds(3)=1.0, dt=0.0
    integer:: n(3)
    logical:: do_upwind=.false.
  contains
    procedure:: init
    procedure:: dx, dy, dz
  end type
  logical:: first_time=.true.
CONTAINS

!> =============================================================================
SUBROUTINE init (self, dims, ds)
  class(central_t):: self
  integer:: dims(3)
  real:: ds(3)
  !.............................................................................
  call trace_begin('central_t%init', 1)
  self%n = dims
  self%ds = ds
  call trace_end
END SUBROUTINE init

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
  if (self%n(1)==1) then
    dx = 0.0
  else
    m = self%n(1)
    c = 0.5/self%ds(1)
    print*,'central', c
    dx(2:m-1,:,:) = (f(3:m,:,:)-f(1:m-2,:,:))*c  ! this part is local
    recv = f(m,:,:)
    dx(1    ,:,:) = (f(2  ,:,:)-recv        )*c  ! remote (n,:,:) slice needed!
    recv = f(1,:,:)
    dx(  m  ,:,:) = (recv      -f(  m-1,:,:))*c  ! remote (1,:,:) slice needed!
  end if
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
  if (self%n(2)==1) then
    dy = 0.0
  else
    m = self%n(2)
    c = 0.5/self%ds(2)
    dy(:,2:m-1,:) = (f(:,3:m,:)-f(:,1:m-2,:))*c  ! this part is local
    recv = f(:,m,:)
    dy(:,1    ,:) = (f(:,2  ,:)-recv        )*c  ! remote (:,n,:) slice needed!
    recv = f(:,1,:)
    dy(:,  m  ,:) = (recv      -f(:,  m-1,:))*c  ! remote (:,1,:) slice needed!
  end if
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
  if (self%n(3)==1) then
    dz = 0.0
  else
    m = self%n(3)
    c = 0.5/self%ds(3)
    dz(:,:,2:m-1) = (f(:,:,3:m)-f(:,:,1:m-2))*c  ! this part is local
    recv = f(:,:,m)
    dz(:,:,1    ) = (f(:,:,2  )-recv        )*c  ! remote (:,:,n) slice needed!
    recv = f(:,:,1)
    dz(:,:,  m  ) = (recv      -f(:,:,  m-1))*c  ! remote (:,:,1) slice needed!
  end if
  call trace_end
END FUNCTION dz
!> =============================================================================
END MODULE central_mod

