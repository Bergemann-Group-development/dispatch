!> =============================================================================
!> $Id$
!> =============================================================================
MODULE stagger_mod
  implicit none
  public
  integer, private:: verbose=0
  logical, save:: hardwire=.false.
  type, public:: stagger_t
  contains
    procedure, nopass:: ddxdn1
    procedure, nopass:: ddydn1
    procedure, nopass:: ddzdn1
    procedure, nopass:: ddxup1
    procedure, nopass:: ddyup1
    procedure, nopass:: ddzup1
    procedure, nopass:: xdn1
    procedure, nopass:: ydn1
    procedure, nopass:: zdn1
    procedure, nopass:: xup1
    procedure, nopass:: yup1
    procedure, nopass:: zup1
    procedure, nopass:: ddxdn
    procedure, nopass:: ddydn
    procedure, nopass:: ddzdn
    procedure, nopass:: ddxup
    procedure, nopass:: ddyup
    procedure, nopass:: ddzup
    procedure, nopass:: xdn
    procedure, nopass:: ydn
    procedure, nopass:: zdn
    procedure, nopass:: xup
    procedure, nopass:: yup
    procedure, nopass:: zup
  end type
CONTAINS

!> =============================================================================
FUNCTION ddxdn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: ds(3)
  b = ddxdn (ds, a)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: ds(3)
  b = ddydn (ds, a)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: ds(3)
  b = ddzdn (ds, a)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: ds(3)
  b = ddxup (ds, a)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: ds(3)
  b = ddyup (ds, a)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8) :: ds(3)
  b = ddzup (ds, a)
END FUNCTION ddzup1

!*******************************************************************************
FUNCTION xdn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = xdn (a)
END FUNCTION xdn1

!*******************************************************************************
FUNCTION ydn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ydn (a)
END FUNCTION ydn1

!*******************************************************************************
FUNCTION zdn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = zdn (a)
END FUNCTION zdn1

!*******************************************************************************
FUNCTION xup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = xup (a)
END FUNCTION xup1

!*******************************************************************************
FUNCTION yup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = yup (a)
END FUNCTION yup1

!*******************************************************************************
FUNCTION zup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = zup (a)
END FUNCTION zup1

!*******************************************************************************
FUNCTION ddxdn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  real(8) :: ds(3)
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,1) > 1) then
    b(1:2,:,:) = 0.0
    b(size(a,1),:,:) = 0.0
    c1 = -1./(24.*ds(1))
    c2 = 1./ds(1)-3.*c1
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=3,size(a,1)-1
      b(ix,iy,iz) = c2*(a(ix  ,iy,iz) - a(ix-1,iy,iz)) &
                  + c1*(a(ix+1,iy,iz) - a(ix-2,iy,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddxdn

!*******************************************************************************
FUNCTION ddydn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz
  real    :: c1, c2
  real(8) :: ds(3)
!...............................................................................
  if (size(a,2) > 1) then
    b(:,1:2,:) = 0.0
    b(:,size(a,2),:) = 0.0
    c1 = -1./(24.*ds(2))
    c2 = 1./ds(2)-3.*c1
    do iz=1,size(a,3)
    do iy=3,size(a,2)-1
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy  ,iz) - a(ix,iy-1,iz)) &
                  + c1*(a(ix,iy+1,iz) - a(ix,iy-2,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddydn

!*******************************************************************************
FUNCTION ddzdn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  real(8) :: ds(3)
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,3) > 1) then
    b(:,:,1:2) = 0.0
    b(:,:,size(a,3)) = 0.0
    c1 = -1./(24.*ds(3))
    c2 = 1./ds(3)-3.*c1
    do iz=3,size(a,3)-1
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy,iz  ) - a(ix,iy,iz-1)) &
                  + c1*(a(ix,iy,iz+1) - a(ix,iy,iz-2))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddzdn

!*******************************************************************************
FUNCTION ddxup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz
  real    :: c1, c2
  real(8) :: ds(3)
!...............................................................................
  if (size(a,1) > 1) then
    b(1,:,:) = 0.0
    b(size(a,1)-1:,:,:) = 0.0
    c1 = -1./(24.*ds(1))
    c2 = 1./ds(1)-3.*c1
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=2,size(a,1)-2
      b(ix,iy,iz) = c2*(a(ix+1,iy,iz) - a(ix  ,iy,iz)) &
                  + c1*(a(ix+2,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddxup

!*******************************************************************************
FUNCTION ddyup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz
  real    :: c1, c2
  real(8) :: ds(3)
!...............................................................................
  if (size(a,2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,1,:) = 0.0
    b(:,size(a,2)-1:,:) = 0.0
    c1 = -1./(24.*ds(2))
    c2 = 1./ds(2)-3.*c1
    do iz=1,size(a,3)
    do iy=2,size(a,2)-2
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy+1,iz) - a(ix,iy  ,iz)) &
                  + c1*(a(ix,iy+2,iz) - a(ix,iy-1,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddyup

!*******************************************************************************
FUNCTION ddzup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz
  real    :: c1, c2
  real(8) :: ds(3)
!...............................................................................
  if (size(a,3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,1) = 0.0
    b(:,:,size(a,3)-1:) = 0.0
    c1 = -1./(24.*ds(3))
    c2 = 1./ds(3)-3.*c1
    do iz=2,size(a,3)-2
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy,iz+1) - a(ix,iy,iz  )) &
                  + c1*(a(ix,iy,iz+2) - a(ix,iy,iz-1))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddzup

!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c1 = -1./16.
    c2 = 0.5-c1
    b(1:2,:,:) = a(1:2,:,:)
    b(size(a,1),:,:) = a(size(a,1),:,:)
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=3,size(a,1)-1
      b(ix,iy,iz) = c2*(a(ix  ,iy,iz) + a(ix-1,iy,iz)) &
                  + c1*(a(ix+1,iy,iz) + a(ix-2,iy,iz))
    end do
    end do
    end do
  else
    b = a
    end if
END FUNCTION xdn

!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c1 = -1./16.
    c2 = 0.5-c1
    b(:,1:2,:) = a(:,1:2,:)
    b(:,size(a,2),:) = a(:,size(a,2),:)
    do iz=1,size(a,3)
    do iy=3,size(a,2)-1
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy  ,iz) + a(ix,iy-1,iz)) &
                  + c1*(a(ix,iy+1,iz) + a(ix,iy-2,iz))
    end do
    end do
    end do
  else
    b = a
    end if
END FUNCTION ydn

!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c1 = -1./16.
    c2 = 0.5-c1
    b(:,:,1:2) = a(:,:,1:2)
    b(:,:,size(a,3)) = a(:,:,size(a,3))
    do iz=3,size(a,3)-1
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy,iz  ) + a(ix,iy,iz-1)) &
                  + c1*(a(ix,iy,iz+1) + a(ix,iy,iz-2))
    end do
    end do
    end do
  else
    b = a
    end if
END FUNCTION zdn

!*******************************************************************************
FUNCTION xup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c1 = -1./16.
    c2 = 0.5-c1
    b(1,:,:) = a(1,:,:)
    b(size(a,1)-1:,:,:) = a(size(a,1)-1:,:,:)
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=2,size(a,1)-2
      b(ix,iy,iz) = c2*(a(ix+1,iy,iz) + a(ix  ,iy,iz)) &
                  + c1*(a(ix+2,iy,iz) + a(ix-1,iy,iz))
    end do
    end do
    end do
  else
    b = a
  endif
END FUNCTION

!*******************************************************************************
FUNCTION yup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c1 = -1./16.
    c2 = 0.5-c1
    b(:,1,:) = a(:,1,:)
    b(:,size(a,2)-1:,:) = a(:,size(a,2)-1:,:)
    do iz=1,size(a,3)
    do iy=2,size(a,2)-2
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy+1,iz) + a(ix,iy  ,iz)) &
                  + c1*(a(ix,iy+2,iz) + a(ix,iy-1,iz))
    end do
    end do
    end do
  else
    b = a
  endif
END FUNCTION

!*******************************************************************************
FUNCTION zup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz
!...............................................................................
  if (size(a,3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c1 = -1./16.
    c2 = 0.5-c1
    b(:,:,1) = a(:,:,1)
    b(:,:,size(a,3)-1:) = a(:,:,size(a,3)-1:)
    do iz=2,size(a,3)-2
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c2*(a(ix,iy,iz+1) + a(ix,iy,iz  )) &
                  + c1*(a(ix,iy,iz+2) + a(ix,iy,iz-1))
    end do
    end do
    end do
  else
    b = a
  endif
END FUNCTION

!> =============================================================================
!> Unit test of 4th/3th order operators
!> =============================================================================
SUBROUTINE stagger_test
  integer       :: m(3)
  real, pointer :: a(:,:,:)
  real, pointer :: b(:,:,:)
  real, pointer :: c(:,:,:)
  real(8)       :: ds(3), x, y, z, d
  integer       :: i
  logical       :: ok
  real, parameter :: eps = 2e-6
  !-----------------------------------------------------------------------------
  m = [10,11,12]
  ds = 1d0/m
  allocate (a(m(1),m(2),m(3)))
  allocate (b(m(1),m(2),m(3)))
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,m(1)
    x = i*ds(1)
    a(i,:,:) = x**3
  end do
  b = xdn(a)
  ok = .true.
  do i=4,m(1)-2
    x = (i-0.5d0)*ds(1)
    d = x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'xdn:', ok
  b = xup(a)
  ok = .true.
  do i=3,m(1)-3
    x = (i+0.5d0)*ds(1)
    d = x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'xup:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,m(2)
    x = i*ds(2)
    a(:,i,:) = x**3
  end do
  b = ydn(a)
  ok = .true.
  do i=4,m(2)-2
    y = (i-0.5d0)*ds(2)
    d = y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'ydn:', ok
  b = yup(a)
  ok = .true.
  do i=3,m(2)-3
    y = (i+.5d0)*ds(2)
    d = y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'yup:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,m(3)
    x = i*ds(3)
    a(:,:,i) = x**3
  end do
  b = zdn(a)
  ok = .true.
  do i=4,m(3)-2
    z = (i-0.5d0)*ds(3)
    d = z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'zdn:', ok
  b = zup(a)
  ok = .true.
  do i=3,m(3)-3
    z = (i+0.5d0)*ds(3)
    d = z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'zup:', ok
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,m(1)
    x = i*ds(1)
    a(i,:,:) = 1d0 + x**4
  end do
  b = ddxdn(ds,a)
  ok = .true.
  do i=4,m(1)-2
    x = (i-0.5d0)*ds(1)
    d = 4d0*x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'ddxdn:', ok
  b = ddxup(ds,a)
  ok = .true.
  do i=3,m(1)-3
    x = (i+0.5d0)*ds(1)
    d = 4d0*x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'ddxup:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,m(2)
    x = i*ds(2)
    a(:,i,:) = 1d0 + x**4
  end do
  b = ddydn(ds,a)
  ok = .true.
  do i=4,m(2)-2
    y = (i-0.5d0)*ds(2)
    d = 4d0*y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'ddydn:', ok
  b = ddyup(ds,a)
  ok = .true.
  do i=3,m(2)-3
    y = (i+.5d0)*ds(2)
    d = 4d0*y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'ddyup:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,m(3)
    x = i*ds(3)
    a(:,:,i) = 1d0 + x**4
  end do
  b = ddzdn(ds,a)
  ok = .true.
  do i=4,m(3)-2
    z = (i-0.5d0)*ds(3)
    d = 4d0*z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'ddzdn:', ok
  b = ddzup(ds,a)
  ok = .true.
  do i=3,m(3)-3
    z = (i+0.5d0)*ds(3)
    d = 4d0*z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'ddzup:', ok
  !-----------------------------------------------------------------------------
  deallocate (a, b)
END SUBROUTINE stagger_test

END MODULE
