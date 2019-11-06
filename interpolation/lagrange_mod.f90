!===============================================================================
!> Module for Lagrange interpolation
!===============================================================================
MODULE lagrange_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: lagrange_t
    integer:: verbose=0
    integer:: order=2
    logical:: initialized=.false.
  contains
    procedure interpolation4
    procedure interpolation8
    generic:: interpolation => interpolation4, interpolation8
    procedure sequence4
    procedure sequence8
    generic:: sequence => sequence4, sequence8
    procedure interval4
    procedure interval8
    generic:: interval => interval4, interval8
    procedure weights4
    procedure weights8
    generic:: weights => weights4, weights8
    procedure deriv_weights4
    procedure deriv_weights8
    generic:: deriv_weights => deriv_weights4, deriv_weights8
    procedure deriv2_weights4
    procedure deriv2_weights8
    generic:: deriv2_weights => deriv2_weights4, deriv2_weights8
    procedure sequence_weights4
    procedure sequence_weights8
    generic:: sequence_weights => sequence_weights4, sequence_weights8
    procedure deriv_sequence_weights4
    procedure deriv_sequence_weights8
    generic:: deriv_sequence_weights => deriv_sequence_weights4, deriv_sequence_weights8
    procedure deriv2_sequence_weights4
    procedure deriv2_sequence_weights8
    generic:: deriv2_sequence_weights => deriv2_sequence_weights4, deriv2_sequence_weights8
    procedure interpolate3d
    procedure trilinear3d
    procedure test
    procedure debug
    procedure init
    procedure interpolate1d4
    procedure interpolate1d8
    generic:: interpolate1d => interpolate1d4, interpolate1d8
  end type
  type(lagrange_t), public:: lagrange
CONTAINS

!> -----------------------------------------------------------------------------
!> Classical Lagrange interpolation
!> -----------------------------------------------------------------------------
FUNCTION interpolation4 (self, x0, x, y) RESULT(f)
  class(lagrange_t):: self
  real(kind=4), dimension(:):: x, y
  real(kind=4):: x0, f, g
  integer:: n, i, j
  !.............................................................................
  f = 0.0
  n = size(x)
  do i=1,n
    g = y(i)
    do j=1,n
      g = merge(g,g*(x0-x(j))/(x(i)-x(j)),j==i)
    end do
    f = f + g
  end do
END FUNCTION interpolation4
FUNCTION interpolation8 (self, x0, x, y) RESULT(f)
  class(lagrange_t):: self
  real(kind=8), dimension(:):: x, y
  real(kind=8):: x0, f, g
  integer:: n, i, j
  !.............................................................................
  f = 0.0
  n = size(x)
  do i=1,n
    g = y(i)
    do j=1,n
      g = merge(g,g*(x0-x(j))/(x(i)-x(j)),j==i)
    end do
    f = f + g
  end do
END FUNCTION interpolation8

!> -----------------------------------------------------------------------------
!> Classical Lagrange interpolation weights.  w(i) is the contribution factor
!> of a value y(i) to the function value at x0, given that there are mesh points
!> at x = [x(j), j=1,n], so the function value f(x0) = sum(y(i)*w(i))
!> -----------------------------------------------------------------------------
FUNCTION weights4 (self, x0, x) RESULT(w)
  class(lagrange_t):: self
  real(kind=4), dimension(:):: x
  real(kind=4), dimension(size(x)):: w
  real(kind=4):: x0, g
  integer:: n, i, j
  !.............................................................................
  n = size(x)
  do i=1,n
    w(i) = 1.0
    do j=1,n
      if (j==i) cycle
      if (x(j)==x(i)) then
        w = 0.0
        w(1) = 1.0
        return
      end if
      w(i) = w(i)*(x0-x(j))/(x(i)-x(j))
    end do
  end do
END FUNCTION weights4
FUNCTION weights8 (self, x0, x) RESULT(w)
  class(lagrange_t):: self
  real(kind=8), dimension(:):: x
  real(kind=8), dimension(size(x)):: w
  real(kind=8):: x0, g
  integer:: n, i, j
  !.............................................................................
  n = size(x)
  do i=1,n
    w(i) = 1.0
    do j=1,n
      if (j==i) cycle
      if (x(j)==x(i)) then
        w = 0d0
        w(1) = 1d0
        return
      end if
      w(i) = w(i)*(x0-x(j))/(x(i)-x(j))
    end do
  end do
END FUNCTION weights8

!> -----------------------------------------------------------------------------
!> Classical Lagrange interpolation weights.  w(i) is the contribution factor
!> of a value y(i) to the derivative value at x0, given that there are mesh
!> points at x = [x(j), j=1,n], so the function value f(x0) = sum(y(i)*w(i))
!> -----------------------------------------------------------------------------
FUNCTION deriv_weights4 (self, x0, x) RESULT(w)
  class(lagrange_t):: self
  real(kind=4), dimension(:):: x
  real(kind=4), dimension(size(x)):: w
  real(kind=4):: x0, nomin, denom, g
  integer:: n, i, j, k
  !.............................................................................
  n = size(x)
  do i=1,n
    denom = 1.0_4
    nomin = 0.0_4
    do j=1,n
      if (j==i) cycle
      if (x(i)==x(j)) then
        w = 0.0
        return
      end if
      g = 1.
      do k=1,n
        if (k==j .or. k==i) cycle
        g = g*(x0-x(k))
      end do
      nomin = nomin + g
      denom = denom*(x(i)-x(j))
    end do
    w(i) = nomin/denom
  end do
END FUNCTION deriv_weights4
FUNCTION deriv_weights8 (self, x0, x) RESULT(w)
  class(lagrange_t):: self
  real(kind=8), dimension(:):: x
  real(kind=8), dimension(size(x)):: w
  real(kind=8):: x0, nomin, denom, g
  integer:: n, i, j, k
  !.............................................................................
  n = size(x)
  do i=1,n
    denom = 1.0_8
    nomin = 0.0_8
    do j=1,n
      if (j==i) cycle
      if (x(i)==x(j)) then
        w = 0.0_8
        return
      end if
      g = 1.0_8
      do k=1,n
        if (k==j .or. k==i) cycle
        g = g*(x0-x(k))
      end do
      nomin = nomin + g
      denom = denom*(x(i)-x(j))
    end do
    w(i) = nomin/denom
  end do
END FUNCTION deriv_weights8

!> -----------------------------------------------------------------------------
!> Classical Lagrange interpolation weights.  w(i) is the contribution factor
!> of a value y(i) to the 2nd derivative value at x0, given that there are mesh
!> points at x = [x(j), j=1,n], so the function value f(x0) = sum(y(i)*w(i))
!> -----------------------------------------------------------------------------
FUNCTION deriv2_weights4 (self, x0, x) RESULT(w)
  class(lagrange_t):: self
  real(kind=4), dimension(:):: x
  real(kind=4), dimension(size(x)):: w
  real(kind=4):: x0, nomin, denom, g
  integer:: n, i, j, k, l
  !.............................................................................
  n = size(x)
  do i=1,n
    denom = 1.0_4
    nomin = 0.0_4
    do j=1,n
      if (j==i) cycle
      if (x(i)==x(j)) then
        w = 0.0
        return
      end if
      do k=1,n
        if (k==j .or. k==i) cycle
        g = 1.
        do l=1,n
          if (l==k .or. l==j .or. l==i) cycle
          g = g*(x0-x(l))
        end do
        nomin = nomin + g
      end do
      denom = denom*(x(i)-x(j))
    end do
    w(i) = nomin/denom
  end do
END FUNCTION deriv2_weights4
FUNCTION deriv2_weights8 (self, x0, x) RESULT(w)
  class(lagrange_t):: self
  real(kind=8), dimension(:):: x
  real(kind=8), dimension(size(x)):: w
  real(kind=8):: x0, nomin, denom, g
  integer:: n, i, j, k, l
  !.............................................................................
  n = size(x)
  do i=1,n
    denom = 1.0_8
    nomin = 0.0_8
    do j=1,n
      if (j==i) cycle
      if (x(i)==x(j)) then
        w = 0.0_8
        return
      end if
      do k=1,n
        if (k == j .or. k == i) cycle
        g = 1.0_8
        do l=1,n
          if (l == k .or. l == j .or. l == i) cycle
          g = g*(x0-x(l))
        end do
        nomin = nomin + g
      end do
      denom = denom*(x(i)-x(j))
    end do
    w(i) = nomin/denom
  end do
END FUNCTION deriv2_weights8

!> -----------------------------------------------------------------------------
!> Prepare to use Lagrange interpolation of given order to interpolate to
!> coordinate x0 in a sequence that may be longer than needed
!> -----------------------------------------------------------------------------
SUBROUTINE interval4 (self, x0, x, i0, i2, o, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=4), dimension(:):: x
  real(kind=4), dimension(size(x)):: y
  real(kind=4):: x0, f
  integer:: i0, i1, i2, o, n
  !.............................................................................
  ! Default order = 3, or n-1 if that is smaller
  n = size(x)
  o = 3
  if (present(order)) o = order
  o = min(o,n-1)
  ! Search for identical initial points
  i1 = 1
  do while (x(i1)==x(1) .and. i1<n)
    i1 = i1+1
  end do
  ! All points are equal to x(1)
  if (x(i1)==x(1)) then
    i0 = 1
    i2 = 1
    o = 0
    return
  ! Only the last point is not equal to x(1)
  else if (i1==n) then
    i0 = n-1
    i2 = n
    o = 1
    return
  end if
  ! Normal search
  do while (x(i1)<x0 .and. i1<n)
    i1 = i1+1
  end do
  if (i1 < 3 .or. i1 > n-1) then
    o = min(o,2)
  end if

  if (i1 < 2+o/2) then
    ! We are at the start of the interval, and need 1 point before i1
    i0 = i1 - (o+1)/2
    if (i0 < 1) then
      i0 = 1
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
END SUBROUTINE interval4
SUBROUTINE interval8 (self, x0, x, i0, i2, o, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=8), dimension(:):: x
  real(kind=8), dimension(size(x)):: y
  real(kind=8):: x0, f
  integer:: i0, i1, i2, o, n
  !.............................................................................
  ! Default order = 3, or n-1 if that is smaller
  n = size(x)
  o = 3
  if (present(order)) o = order
  o = min(o,n-1)
  ! Search for identical initial points
  i1 = 1
  do while (x(i1)==x(1) .and. i1<n)
    i1 = i1+1
  end do
  ! All points are equal to x(1)
  if (x(i1)==x(1)) then
    i0 = 1
    i2 = 1
    o = 0
    return
  ! Only the last point is not equal to x(1)
  else if (i1==n) then
    i0 = n-1
    i2 = n
    o = 1
    return
  end if
  ! Normal search
  do while (x(i1)<x0 .and. i1<n)
    i1 = i1+1
  end do
  if (i1 < 3 .or. i1 > n-1) then
    o = min(o,2)
  end if

  if (i1 < 2+o/2) then
    ! We are at the start of the interval, and need 1 point before i1
    i0 = i1 - (o+1)/2
    if (i0 < 1) then
      i0 = 1
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
    end if
    i0 = i2 - o
    if (i0 < 1) then
      i0 = 1
      i2 = min(i2,n)
      o = min(o,i2-i0)
    end if
  end if
  if (self%debug(1)) then
    print '("lagrange_interval: i0,i2:",2i3,3x,"o:",i2,3x,"x0:",f10.6,3x,"x:",10f10.6)', &
                                i0,i2,          o,         x0,           x
  end if
END SUBROUTINE interval8

!> -----------------------------------------------------------------------------
!> Use Lagrange interpolation of given order to interpolate to coordinate x0
!> in a sequence that may be longer than needed
!> -----------------------------------------------------------------------------
FUNCTION sequence4 (self, x0, x, y, order, correct) RESULT(f)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=4), optional:: correct
  real(kind=4), dimension(:):: x, y
  real(kind=4):: x0, f
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  if (self%verbose>0) then
    print'(" i0,i2,order=",3i3,4x,"x0,x,y  =",1p,18g12.4)', &
      i0, i2, o, x0, x(i0:i2), y(i0:i2)
  end if
  f = self%interpolation (x0, x(i0:i2), y(i0:i2))
  if (self%verbose>0) then
    if (present(correct)) then
      print'(" i0,i2,order=",3i3,4x,"x0,x,y,f,correct=",1p,18g12.4)', &
        i0, i2, o, x0, x(i0:i2), y(i0:i2), f, correct
    else
      print'(" i0,i2,order=",3i3,4x,"x0,x,y,f=",1p,18g12.4)', &
        i0, i2, o, x0, x(i0:i2), y(i0:i2), f
    end if
  end if
END FUNCTION sequence4
FUNCTION sequence8 (self, x0, x, y, order, correct) RESULT(f)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=8), optional:: correct
  real(kind=8), dimension(:):: x, y
  real(kind=8):: x0, f
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  if (self%verbose>0) then
    print'(" i0,i2,order=",3i3,4x,"x0,x,y  =",1p,18g12.4)', &
      i0, i2, o, x0, x(i0:i2), y(i0:i2)
  end if
  f = self%interpolation (x0, x(i0:i2), y(i0:i2))
  if (self%verbose>0) then
    if (present(correct)) then
      print'(" i0,i2,order=",3i3,4x,"x0,x,y,f,correct=",1p,18g12.4)', &
        i0, i2, o, x0, x(i0:i2), y(i0:i2), f, correct
    else
      print'(" i0,i2,order=",3i3,4x,"x0,x,y,f=",1p,18g12.4)', &
        i0, i2, o, x0, x(i0:i2), y(i0:i2), f
    end if
  end if
END FUNCTION sequence8

!> -----------------------------------------------------------------------------
!> Use Lagrange interpolation of given order to get weights for interpolate to
!> coordinate x0 in a sequence of increasing x that may be longer than needed
!> -----------------------------------------------------------------------------
SUBROUTINE sequence_weights4 (self, x0, x, i0, i2, w, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=4), dimension(:):: x, w
  real(kind=4):: x0
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  w = 0d0
  if (o==0) then
    w(1) = 1.0
  else
    w(1:i2-i0+1) = self%weights (x0, x(i0:i2))
  end if
END SUBROUTINE sequence_weights4
SUBROUTINE sequence_weights8 (self, x0, x, i0, i2, w, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=8), dimension(:):: x, w
  real(kind=8):: x0
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  w = 0d0
  if (o==0) then
    w(1) = 1.0
  else
    w(1:i2-i0+1) = self%weights (x0, x(i0:i2))
  end if
END SUBROUTINE sequence_weights8

!> -----------------------------------------------------------------------------
!> Use Lagrange interpolation of given order to get weights for interpolate to
!> coordinate x0 in a sequence of increasing x that may be longer than needed
!> -----------------------------------------------------------------------------
SUBROUTINE deriv_sequence_weights4 (self, x0, x, i0, i2, w, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=4), dimension(:):: x, w
  real(kind=4):: x0
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  w = 0d0
  w(1:i2-i0+1) = self%deriv_weights (x0, x(i0:i2))
END SUBROUTINE deriv_sequence_weights4
SUBROUTINE deriv_sequence_weights8 (self, x0, x, i0, i2, w, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=8), dimension(:):: x, w
  real(kind=8):: x0
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  w = 0d0
  w(1:i2-i0+1) = self%deriv_weights (x0, x(i0:i2))
END SUBROUTINE deriv_sequence_weights8

!> -----------------------------------------------------------------------------
!> Use Lagrange interpolation of given order to get weights for interpolate to
!> coordinate x0 in a sequence of increasing x that may be longer than needed
!> -----------------------------------------------------------------------------
SUBROUTINE deriv2_sequence_weights4 (self, x0, x, i0, i2, w, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=4), dimension(:):: x, w
  real(kind=4):: x0
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  w = 0d0
  w(1:i2-i0+1) = self%deriv2_weights (x0, x(i0:i2))
END SUBROUTINE deriv2_sequence_weights4
SUBROUTINE deriv2_sequence_weights8 (self, x0, x, i0, i2, w, order)
  class(lagrange_t):: self
  integer, optional:: order
  real(kind=8), dimension(:):: x, w
  real(kind=8):: x0
  integer:: i0, i2, o, n
  !.............................................................................
  n = size(x)
  call self%interval (x0, x, i0, i2, o, order)
  w = 0d0
  w(1:i2-i0+1) = self%deriv2_weights (x0, x(i0:i2))
END SUBROUTINE deriv2_sequence_weights8

!> -----------------------------------------------------------------------------
!> Interpolate in 3D from fi to fo
!> -----------------------------------------------------------------------------
SUBROUTINE interpolate3d (self, xi, yi, zi, fi, xo, yo, zo, fo, order)
  class(lagrange_t):: self
  real, dimension(:):: xi, yi, zi, xo, yo, zo
  real, dimension(:,:,:):: fi, fo
  integer:: order
  !.............................................................................
  real, dimension(:,:,:), pointer:: f1, f2
  integer:: nix, niy, niz, nox, noy, noz, ix, iy, iz
  !-----------------------------------------------------------------------------
  nix = size(xi)
  niy = size(yi)
  niz = size(zi)
  nox = size(xo)
  noy = size(yo)
  noz = size(zo)
  allocate (f1(nix,niy,noz))
  do iz=1,noz
  do iy=1,niy
  do ix=1,nix
    f1(ix,iy,iz) = self%sequence (zo(iz), zi, fi(ix,iy,:), order)
  end do
  end do
  end do
  allocate (f2(nix,noy,noz))
  do iz=1,noz
  do iy=1,noy
  do ix=1,nix
    f2(ix,iy,iz) = self%sequence (yo(iy), yi, f1(ix,:,iz), order)
  end do
  end do
  end do
  deallocate (f1)
  do iz=1,noz
  do iy=1,noy
  do ix=1,nox
    fo(ix,iy,iz) = self%sequence (xo(ix), xi, f2(:,iy,iz), order)
  end do
  end do
  end do
  deallocate (f2)
END SUBROUTINE interpolate3d

!> -----------------------------------------------------------------------------
!> Interpolate in 3D from fi to fo
!> -----------------------------------------------------------------------------
SUBROUTINE trilinear3d (self, xi, yi, zi, fi, xo, yo, zo, fo, order)
  class(lagrange_t):: self
  real, dimension(:):: xi, yi, zi, xo, yo, zo
  real, dimension(:,:,:):: fi, fo
  integer:: order
  !.............................................................................
  integer:: nix, niy, niz, nox, noy, noz, ix, iy, iz, jx ,jy, jz
  real:: px, py, pz, qx, qy, qz, dx, dy, dz
  !-----------------------------------------------------------------------------
  nix = size(xi)
  niy = size(yi)
  niz = size(zi)
  nox = size(xo)
  noy = size(yo)
  noz = size(zo)
  dx = (xi(nix)-xi(1))/(nix-1)
  dy = (yi(niy)-yi(1))/(niy-1)
  dz = (zi(niz)-zi(1))/(niz-1)
  do iz=1,noz
    do jz=1,niz-2
      if (zi(jz+1)>zo(iz)) exit
    end do
    pz = (zo(iz)-zi(jz))/(zi(jz+1)-zi(jz))
    pz = max(0.0,min(1.0,pz))
    qz = 1.0-pz
    do iy=1,noy
      py = (yo(iy)-yi(1))/dy + 1.
      jy = py
      jy = max(1,min(niy-1,jy))
      py = max(0.0,min(1.0,py -jy))
      qy = 1.0-py
      do ix=1,nox
        px = (xo(ix)-xi(1))/dx + 1.
        jx = px
        jx = max(1,min(nix-1,jx))
        px = max(0.0,min(1.0,px -jx))
        qx = 1.0-px
        fo(ix,iy,iz) = &
          qz*(qy*(qx*fi(jx  ,jy  , jz  ) + px*fi(jx+1,jy  , jz  )) + &
              py*(qx*fi(jx  ,jy+1, jz  ) + px*fi(jx+1,jy+1, jz  ))) + &
          pz*(qy*(qx*fi(jx  ,jy  , jz+1) + px*fi(jx+1,jy  , jz+1)) + &
              py*(qx*fi(jx  ,jy+1, jz+1) + px*fi(jx+1,jy+1, jz+1)))
      end do
    end do
  end do
END SUBROUTINE trilinear3d

!> -----------------------------------------------------------------------------
!> Test for a polynomial of order 2
!> -----------------------------------------------------------------------------
SUBROUTINE test (self)
  class(lagrange_t):: self
  integer, parameter:: n=5
  integer:: order, i, i0, i1
  real(8), dimension(n):: x, y, w
  real(8):: x0, y0
  !-----------------------------------------------------------------------------
  if (.not. io_unit%master) return
  print *,'--------------- lagrange_t%test ------------------'
  x = 0.0
  y = 1.0
  order = 1
  print '(5x,a,11x,a,11x,a,4x,a)', 'x', 'y', 'y==1', 'order'
  do i=1,n+1
    x0 = i-0.5
    call self%sequence_weights (x0, x, i0, i1, w, order)
    y0 = sum(y(i0:i1)*w(1:1+i1-i0))
    print '(3g12.4,i6)', x0, y0, x0**order, i1-i0
  end do
  y = 0.0
  x(n) = 1.0
  y(n) = 1.0
  print '(5x,a,11x,a,11x,a,4x,a)', 'x', 'y', 'x-4', 'order'
  do i=1,n+1
    x0 = i-0.5
    call self%sequence_weights (x0, x, i0, i1, w, order)
    y0 = sum(y(i0:i1)*w(1:1+i1-i0))
    print '(3g12.4,i6)', x0, y0, x0**order, i1-i0
  end do
  order = 2
  do i=1,n
    x(i) = i
    y(i) = x(i)**order
  end do
  print '(5x,a,11x,a,11x,a,4x,a)', 'x', 'y', 'x**2', 'order'
  do i=1,n+1
    x0 = i-0.5
    call self%sequence_weights (x0, x, i0, i1, w, order)
    y0 = sum(y(i0:i1)*w(1:1+i1-i0))
    print '(3g12.4,i6)', x0, y0, x0**order, i1-i0
  end do
  print *,'--------------------------------------------------'
END SUBROUTINE test

!> -----------------------------------------------------------------------------
!> -----------------------------------------------------------------------------
FUNCTION debug (self, verbose)
  class(lagrange_t):: self
  integer:: verbose
  logical:: debug
  !-----------------------------------------------------------------------------
  if (.not.self%initialized) call self%init
  debug = (self%verbose >= verbose)
END FUNCTION debug

!> -----------------------------------------------------------------------------
!> Initialize (only) this instance of lagrange
!> -----------------------------------------------------------------------------
SUBROUTINE init (self)
  class(lagrange_t):: self
  integer:: iostat
  integer:: verbose=0, order=2
  namelist /lagrange_params/ verbose, order
  !-----------------------------------------------------------------------------
  call trace%begin ('lagrange_t%init')
  !$omp critical (input_cr)
  if (.not.self%initialized) then
    self%initialized = .true.
    rewind (io_unit%input)
    read (io_unit%input, lagrange_params, iostat=iostat)
    if (io_unit%master) write (io_unit%output, lagrange_params)
    self%verbose = verbose
    self%order = order
  end if
  !$omp end critical (input_cr)
  call trace%end ()
END SUBROUTINE init

!> =============================================================================
!> Use default 1-D Lagrange interpolation
!> =============================================================================
FUNCTION interpolate1d4 (self,x0, n, x, y) RESULT(f)
  class(lagrange_t):: self
  integer:: n
  real(kind=4):: x0, x(n), y(n), f
  !.............................................................................
  if (x0 < x(1)) then
    f = y(1)
  else if (x0 > x(n)) then
    f = y(n)
  else
    f = self%sequence (x0, x, y)
  end if
END FUNCTION interpolate1d4

FUNCTION interpolate1d8 (self,x0, n, x, y) RESULT(f)
  class(lagrange_t):: self
  integer:: n
  real(kind=8):: x0, x(n), y(n), f
  !.............................................................................
  if (x0 < x(1)) then
    f = y(1)
  else if (x0 > x(n)) then
    f = y(n)
  else
    f = self%sequence (x0, x, y)
  end if
END FUNCTION interpolate1d8

END MODULE
