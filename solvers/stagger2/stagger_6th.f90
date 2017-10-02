!> =============================================================================
!> 6th order stagger operators, with self-test procedure
!> =============================================================================
MODULE stagger_mod
  implicit none
  public
  integer, private:: verbose=0
  logical, save:: hardwire=.false.
  type, public:: stagger_t
  contains
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
    procedure, nopass:: stagger_test
  end type
CONTAINS

!*******************************************************************************
FUNCTION ddxdn (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(1) > 5) then
    d = ds(1)
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(    1:3,:,:) = 0.0
    b(m(1)-1:,:,:) = 0.0
    do k=1,m(3)
    do j=1,m(2)
    do i=4,m(1)-2
      b(i  ,j,k) = (c3*(a(i+2,j,k)-a(i-3,j,k)) &
                  + c2*(a(i+1,j,k)-a(i-2,j,k))) &
                  + c1*(a(i  ,j,k)-a(i-1,j,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddxdn

!*******************************************************************************
FUNCTION ddydn (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(2) > 5) then
    d = ds(2)
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,    1:3,:) = 0.0
    b(:,m(2)-1:,:) = 0.0
    do k=1,m(3)
    do j=4,m(2)-2
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j+2,k)-a(i,j-3,k)) &
                + c2*(a(i,j+1,k)-a(i,j-2,k))) &
                + c1*(a(i,j  ,k)-a(i,j-1,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddydn

!*******************************************************************************
FUNCTION ddzdn (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(3) > 5) then
    d = ds(3)
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,:,    1:3) = 0.0
    b(:,:,m(3)-1:) = 0.0
    do k=4,m(3)-2
    do j=1,m(2)
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j,k+2)-a(i,j,k-3)) &
                + c2*(a(i,j,k+1)-a(i,j,k-2))) &
                + c1*(a(i,j,k  )-a(i,j,k-1))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddzdn

!*******************************************************************************
FUNCTION ddxup (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  d = ds(1)
  if (m(1) > 5) then
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(    1:2,:,:) = 0.0
    b(m(1)-2:,:,:) = 0.0
    do k=1,m(3)
    do j=1,m(2)
    do i=3,m(1)-3
      b(i  ,j,k) = (c3*(a(i+3,j,k)-a(i-2,j,k)) &
                  + c2*(a(i+2,j,k)-a(i-1,j,k))) &
                  + c1*(a(i+1,j,k)-a(i  ,j,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddxup

!*******************************************************************************
FUNCTION ddyup (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  d = ds(2)
  if (m(1) > 5) then
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,    1:2,:) = 0.0
    b(:,m(2)-2:,:) = 0.0
    do k=1,m(3)
    do j=3,m(2)-3
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j+3,k)-a(i,j-2,k)) &
                + c2*(a(i,j+2,k)-a(i,j-1,k))) &
                + c1*(a(i,j+1,k)-a(i,j  ,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddyup

!*******************************************************************************
FUNCTION ddzup (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  d = ds(3)
  if (m(1) > 5) then
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,:,    1:2) = 0.0
    b(:,:,m(1)-2:) = 0.0
    do k=3,m(3)-3
    do j=1,m(2)
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j,k+3)-a(i,j,k-2)) &
                + c2*(a(i,j,k+2)-a(i,j,k-1))) &
                + c1*(a(i,j,k+1)-a(i,j,k  ))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ddzup

!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(    1:2,:,:) = a(    1:2,:,:)
    b(m(1)-2:,:,:) = a(m(1)-2:,:,:)
    do k=1,m(3)
    do j=1,m(2)
    do i=4,m(1)-2
      b(i  ,j,k) = (c3*(a(i+2,j,k)+a(i-3,j,k)) &
                  + c2*(a(i+1,j,k)+a(i-2,j,k))) &
                  + c1*(a(i  ,j,k)+a(i-1,j,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION xdn

!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,    1:2,:) = a(:,    1:2,:)
    b(:,m(2)-2:,:) = a(:,m(2)-2:,:)
    do k=1,m(3)
    do j=4,m(2)-2
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j+2,k)+a(i,j-3,k)) &
                + c2*(a(i,j+1,k)+a(i,j-2,k))) &
                + c1*(a(i,j  ,k)+a(i,j-1,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION ydn

!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,:,    1:2) = a(:,:,    1:2)
    b(:,:,m(1)-2:) = a(:,:,m(1)-2:)
    do k=4,m(3)-2
    do j=1,m(2)
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j,k+2)+a(i,j,k-3)) &
                + c2*(a(i,j,k+1)+a(i,j,k-2))) &
                + c1*(a(i,j,k  )+a(i,j,k-1))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION zdn

!*******************************************************************************
FUNCTION xup (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(    1:2,:,:) = a(    1:2,:,:)
    b(m(1)-2:,:,:) = a(m(1)-2:,:,:)
    do k=1,m(3)
    do j=1,m(2)
    do i=3,m(1)-3
      b(i  ,j,k) = (c3*(a(i+3,j,k)+a(i-2,j,k)) &
                  + c2*(a(i+2,j,k)+a(i-1,j,k))) &
                  + c1*(a(i+1,j,k)+a(i  ,j,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION xup

!*******************************************************************************
FUNCTION yup (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,    1:2,:) = a(:,    1:2,:)
    b(:,m(2)-2:,:) = a(:,m(2)-2:,:) 
    do k=1,m(3)
    do j=3,m(2)-3
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j+3,k)+a(i,j-2,k)) &
                + c2*(a(i,j+2,k)+a(i,j-1,k))) &
                + c1*(a(i,j+1,k)+a(i,j  ,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION yup

!*******************************************************************************
FUNCTION zup (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, m(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  m = shape(a)
  if (m(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,:,    1:2) = a(:,:,    1:2)
    b(:,:,m(1)-2:) = a(:,:,m(1)-2:)
    do k=3,m(3)-3
    do j=1,m(2)
    do i=1,m(1)
      b(i,j,k) = (c3*(a(i,j,k+3)+a(i,j,k-2)) &
                + c2*(a(i,j,k+2)+a(i,j,k-1))) &
                + c1*(a(i,j,k+1)+a(i,j,k  ))
    end do
    end do
    end do
  else
    b = 0.0
  endif
END FUNCTION zup

!> =============================================================================
!> Unit test of 6th/5th order operators
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
    a(i,:,:) = x**5
  end do
  b = xdn(a)
  ok = .true.
  do i=4,m(1)-2
    x = (i-0.5d0)*ds(1)
    d = x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'xdn:', ok
  b = xup(a)
  ok = .true.
  do i=3,m(1)-3
    x = (i+0.5d0)*ds(1)
    d = x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'xup:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,m(2)
    x = i*ds(2)
    a(:,i,:) = x**5
  end do
  b = ydn(a)
  ok = .true.
  do i=4,m(2)-2
    y = (i-0.5d0)*ds(2)
    d = y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'ydn:', ok
  b = yup(a)
  ok = .true.
  do i=3,m(2)-3
    y = (i+.5d0)*ds(2)
    d = y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'yup:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,m(3)
    x = i*ds(3)
    a(:,:,i) = x**5
  end do
  b = zdn(a)
  ok = .true.
  do i=4,m(3)-2
    z = (i-0.5d0)*ds(3)
    d = z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'zdn:', ok
  b = zup(a)
  ok = .true.
  do i=3,m(3)-3
    z = (i+0.5d0)*ds(3)
    d = z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'zup:', ok
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,m(1)
    x = i*ds(1)
    a(i,:,:) = 1d0 + x**6
  end do
  b = ddxdn(ds,a)
  ok = .true.
  do i=4,m(1)-2
    x = (i-0.5d0)*ds(1)
    d = 6d0*x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'ddxdn:', ok
  b = ddxup(ds,a)
  ok = .true.
  do i=3,m(1)-3
    x = (i+0.5d0)*ds(1)
    d = 6d0*x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  print *,'ddxup:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,m(2)
    x = i*ds(2)
    a(:,i,:) = 1d0 + x**6
  end do
  b = ddydn(ds,a)
  ok = .true.
  do i=4,m(2)-2
    y = (i-0.5d0)*ds(2)
    d = 6d0*y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'ddydn:', ok
  b = ddyup(ds,a)
  ok = .true.
  do i=3,m(2)-3
    y = (i+.5d0)*ds(2)
    d = 6d0*y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  print *,'ddyup:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,m(3)
    x = i*ds(3)
    a(:,:,i) = 1d0 + x**6
  end do
  b = ddzdn(ds,a)
  ok = .true.
  do i=4,m(3)-2
    z = (i-0.5d0)*ds(3)
    d = 6d0*z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'ddzdn:', ok
  b = ddzup(ds,a)
  ok = .true.
  do i=3,m(3)-3
    z = (i+0.5d0)*ds(3)
    d = 6d0*z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  print *,'ddzup:', ok
  !-----------------------------------------------------------------------------
  deallocate (a, b)
END SUBROUTINE stagger_test

END MODULE stagger_mod
