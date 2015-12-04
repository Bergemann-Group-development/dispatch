!> =============================================================================
!> $Id$
!> =============================================================================
MODULE upwind_mod
  USE io_mod
  USE trace_mod
  USE central_mod
  implicit none
  private
  type, public, extends(central_t):: upwind_t
    real, pointer, dimension(:,:,:) :: vx,vy,vz
    real:: factor=1.5
  contains
    procedure:: init
    procedure:: dx
    procedure:: dy
    procedure:: dz
    procedure:: init_upwind
  end type
  logical:: first_time=.true.
CONTAINS

SUBROUTINE init_upwind (self, vx, vy, vz)
  class(upwind_t):: self
  real, dimension(:,:,:), pointer:: vx, vy, vz
  self%vx => vx
  self%vy => vy
  self%vz => vz
END SUBROUTINE

!> =============================================================================
SUBROUTINE init (self, dims, ds)
  class(upwind_t):: self
  integer:: dims(3)
  real:: ds(3)
  !.............................................................................
  call trace_begin('central_t%init', 1)
  self%n = dims
  self%ds = ds
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Derivatives evaluated with upwind differences. We do not care about guard-cells.
!> They are supposed to be updated elsewhere.
!> =============================================================================
FUNCTION dx (self, f, n)
  class(upwind_t):: self
  integer:: m, n(3), i
  real, dimension(n(1),n(2),n(3)):: f, dx
  real, dimension(n(2),n(3)):: dx1, dx2, delta
  real :: c
  !.............................................................................
  if (first_time) then
    first_time = .false.
  end if
  if (n(1)==1) then
    dx = 0.0
  else
    m = n(1)
    c = 0.5/self%ds(1)
    do i=2,m-1
      delta = self%vx(i,:,:)*self%dt*self%factor/self%ds(1)
      dx1 = f(i  ,:,:)-f(i-1,:,:)
      dx2 = f(i+1,:,:)-f(i  ,:,:)
      dx(i,:,:) = c*(dx1*(1.0+delta)+dx2*(1.0-delta))
    end do
  end if
END FUNCTION dx
!> =============================================================================
FUNCTION dy (self, f, n)
  class(upwind_t):: self
  integer:: n(3), m, i
  real, dimension(n(1),n(2),n(3)):: f, dy
  real, dimension(n(1),n(3)):: dy1, dy2, delta
  real:: c
  !.............................................................................
  if (n(2)==1) then
    dy = 0.0
  else
    m = n(2)
    c = 0.5/self%ds(2)
    do i=2,m-1
      delta = self%vy(:,i,:)*self%dt*self%factor/self%ds(2)
      dy1 = f(:,i  ,:)-f(:,i-1,:)
      dy2 = f(:,i+1,:)-f(:,i  ,:)
      dy(:,i,:) = c*(dy1*(1.0+delta)+dy2*(1.0-delta))
    end do
  end if
END FUNCTION dy
!> =============================================================================
FUNCTION dz (self, f, n)
  class(upwind_t):: self
  integer:: n(3), m, i
  real, dimension(n(1),n(2),n(3)):: f, dz
  real, dimension(n(1),n(2)):: dz1, dz2, delta
  real:: c
  !.............................................................................
  if (n(3)==1) then
    dz = 0.0
  else
    m = n(3)
    c = 0.5/self%ds(3)
    do i=2,m-1
      delta = self%vz(:,:,i)*self%dt*self%factor/self%ds(3)
      dz1 = f(:,:,i  )-f(:,:,i-1)
      dz2 = f(:,:,i+1)-f(:,:,i  )
      dz(:,:,i) = c*(dz1*(1.0+delta)+dz2*(1.0-delta))
    end do
  end if
END FUNCTION dz

END MODULE upwind_mod
