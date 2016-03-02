!> -----------------------------------------------------------------------------
!> Module for Lagrange interpolation
!> -----------------------------------------------------------------------------
MODULE lagrange_mod
  implicit none
  private
  type, public:: lagrange_t
    integer:: verbose=0
  contains
    procedure lagrange_interpolation
    procedure lagrange_weights
    procedure lagrange_interval
    procedure lagrange_sequence
    procedure lagrange_sequence_weights
  end type
CONTAINS

!> -----------------------------------------------------------------------------
!> Classical Lagrange interpolation
!> -----------------------------------------------------------------------------
FUNCTION lagrange_interpolation (self, x0, x, y, n) RESULT(f)
  class(lagrange_t):: self
  real(kind=8), dimension(n):: x, y
  real(kind=8):: x0, f, g
  integer:: n, i, j
  integer:: ii(n)
  !.............................................................................
  f = 0.0
  do i=1,n
    g = y(i)
    do j=1,n
      g = merge(g,g*(x0-x(j))/(x(i)-x(j)),j==i)
    end do
    f = f + g
  end do
END FUNCTION lagrange_interpolation

!> -----------------------------------------------------------------------------
!> Classical Lagrange interpolation weights.  w(i) is the contribution factor
!> of a value y(i) to the function value at x0, given that there are mesh points
!> at x = [x(j), j=1,n], so the function value f(x0) = sum(y(i)*w(i))
!> -----------------------------------------------------------------------------
FUNCTION lagrange_weights (self, x0, x, n) RESULT(w)
  class(lagrange_t):: self
  real(kind=8), dimension(n):: x, w
  real(kind=8):: x0, g
  integer:: n, i, j
  !.............................................................................
  do i=1,n
    w(i) = 1.0
    do j=1,n
      w(i) = merge(w(i),w(i)*(x0-x(j))/(x(i)-x(j)),j==i)
    end do
  end do
END FUNCTION lagrange_weights

!> -----------------------------------------------------------------------------
!> Prepare to use Lagrange interpolation of given order to interpolate to
!> coordinate x0 in a sequence that may be longer than needed
!> -----------------------------------------------------------------------------
SUBROUTINE lagrange_interval (self, x0, x, n, i0, i2, o, order)
  class(lagrange_t):: self
  integer:: n
  integer, optional:: order
  real(kind=8), dimension(n):: x, y
  real(kind=8):: x0, f
  integer:: i0, i1, i2, o
  !.............................................................................
  ! Default order = 3, or n-1 if that is smaller
  o = 3
  if (present(order)) o = order
  o = min(o,n-1)
  ! Search for first point beyond x0
  i1 = 1
  do while (x(i1)<x0 .and. i1<n)
    i1 = i1+1
  end do
  if (i1 < 3 .or. i1 > n-1) then
    o = min(o,2)
  end if

  if (i1 < 3+o/2) then
    ! We are at the start of the interval, and need 1 point before i1
    i0 = i1 - (o+1)/2
    if (i0 < 1) then
      i0 = 1
      i1 = i0 + (o+1)/2
    end if
    i2 = i0+o
    if (i2 > n) then
      i2 = n
      o = min(o,i2-i0)
    end if
  else
    ! We are in the middle or end of the interval
    i2 = i1 + (o-1)/2
    ! Shift left if not enough points
    if (i2 > n) then
      i2 = n
      i1 = n-1
    end if
    i0 = i2 - o
    if (i0 < 1) then
      i0 = 1
      i2 = min(i2,n)
      o = min(o,i2-i0)
    end if
  end if
  if (self%verbose>0) then
    print*,'lagrange_interval: i0, i2, o =', i0, i2, o
  end if
END SUBROUTINE lagrange_interval

!> -----------------------------------------------------------------------------
!> Use Lagrange interpolation of given order to interpolate to coordinate x0
!> in a sequence that may be longer than needed
!> -----------------------------------------------------------------------------
FUNCTION lagrange_sequence (self, x0, x, y, n, order, correct) RESULT(f)
  class(lagrange_t):: self
  integer:: n
  integer, optional:: order
  real(kind=8), optional:: correct
  real(kind=8), dimension(n):: x, y
  real(kind=8):: x0, f
  integer:: i0, i2, o
  !.............................................................................
  call self%lagrange_interval (x0, x, n, i0, i2, o, order)
  if (self%verbose>0) then
    print'(" i0,i2,order=",3i3,4x,"x0,x,y  =",1p,18g12.4)', &
      i0, i2, o, x0, x(i0:i2), y(i0:i2)
  end if
  f = self%lagrange_interpolation (x0, x(i0:i2), y(i0:i2), i2-i0+1)
  if (self%verbose>0) then
    if (present(correct)) then
      print'(" i0,i2,order=",3i3,4x,"x0,x,y,f,correct=",1p,18g12.4)', &
        i0, i2, o, x0, x(i0:i2), y(i0:i2), f, correct
    else
      print'(" i0,i2,order=",3i3,4x,"x0,x,y,f=",1p,18g12.4)', &
        i0, i2, o, x0, x(i0:i2), y(i0:i2), f
    end if
  end if
END FUNCTION lagrange_sequence

!> -----------------------------------------------------------------------------
!> Use Lagrange interpolation of given order to get weights for interpolate to
!> coordinate x0 in a sequence of increasing x that may be longer than needed
!> -----------------------------------------------------------------------------
SUBROUTINE lagrange_sequence_weights (self, x0, x, n, i0, i2, w, order)
  class(lagrange_t):: self
  integer:: n
  integer, optional:: order
  real(kind=8), dimension(n):: x, w
  real(kind=8):: x0
  integer:: i0, i2, o
  !.............................................................................
  call self%lagrange_interval (x0, x, n, i0, i2, o, order)
  w = 0d0
  w(1:i2-i0+1) = self%lagrange_weights (x0, x(i0:i2), i2-i0+1)
END SUBROUTINE lagrange_sequence_weights

END MODULE
