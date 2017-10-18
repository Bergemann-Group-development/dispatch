!> -----------------------------------------------------------------------------
!> Module for Lagrange interpolation
!> -----------------------------------------------------------------------------
MODULE lagrange_mod
  USE io_mod
  implicit none
  private
  type, public:: lagrange_t
    integer:: verbose=0
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
    procedure sequence_weights4
    procedure sequence_weights8
    generic:: sequence_weights => sequence_weights4, sequence_weights8
    procedure interpolate3d
    procedure trilinear3d
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
      w(i) = merge(w(i),w(i)*(x0-x(j))/(x(i)-x(j)),j==i)
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
      w(i) = merge(w(i),w(i)*(x0-x(j))/(x(i)-x(j)),j==i)
    end do
  end do
END FUNCTION weights8

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
  w(1:i2-i0+1) = self%weights (x0, x(i0:i2))
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
  w(1:i2-i0+1) = self%weights (x0, x(i0:i2))
END SUBROUTINE sequence_weights8

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
    pz = (zo(iz)-zi(1))/dz
    jz = pz
    jz = max(1,min(niz-1,jz))
    pz = max(0.0,min(1.0,pz -jz))
    qz = 1.0-pz
    do iy=1,noy
      py = (yo(iy)-yi(1))/dy
      jy = py
      jy = max(1,min(niy-1,jy))
      py = max(0.0,min(1.0,py -jy))
      qy = 1.0-py
      do ix=1,nox
        px = (xo(ix)-xi(1))/dx
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

END MODULE
