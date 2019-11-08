!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_mod
  USE trace_mod
  USE omp_timer_mod
  USE mesh_mod
  implicit none
  public
  integer, private:: verbose=0
  type, public:: stagger_t
    real(8):: ds(3)
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
    procedure, nopass:: h2c
    procedure, nopass:: h2f
    procedure, nopass:: h31c
    procedure, nopass:: h31f
    procedure, nopass:: h32c
    procedure, nopass:: h32f
  end type
CONTAINS

!===============================================================================
FUNCTION ddxdn1 (mesh, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real, dimension(size(a,1),size(a,2),size(a,3)):: d1
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = shape(b)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d1(i,j,k) = mesh(1)%d
      end do
    end do
  end do
  b= ddxdn (d1, a)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (mesh, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real, dimension(size(a,1),size(a,2),size(a,3)):: d2
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = shape(b)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d2(i,j,k) = mesh(2)%d
      end do
    end do
  end do
  b= ddydn (d2, a)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (mesh, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real, dimension(size(a,1),size(a,2),size(a,3)):: d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = shape(b)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d3(i,j,k) = mesh(3)%d
      end do
    end do
  end do
  b= ddzdn (d3, a)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (mesh, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real, dimension(size(a,1),size(a,2),size(a,3)):: d1
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = shape(b)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d1(i,j,k) = mesh(1)%d
      end do
    end do
  end do
  b= ddxup (d1, a)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (mesh, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real, dimension(size(a,1),size(a,2),size(a,3)):: d2
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = shape(b)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d2(i,j,k) = mesh(2)%d*mesh(1)%h2c(i)
      end do
    end do
  end do
  b= ddyup (d2, a)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (mesh, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real, dimension(size(a,1),size(a,2),size(a,3)):: d3
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = shape(b)
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        d3(i,j,k) = mesh(3)%d*mesh(1)%h31c(i)*mesh(2)%h32c(j)
      end do
    end do
  end do
  b= ddzup (d3, a)
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
FUNCTION ddxdn (s1, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a, s1
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddxdn', itimer=itimer)
  n = shape(a)
  if (n(1) > 1) then
    b(1,:,:) = 0.0
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=2,n(1)
      b(ix,iy,iz) = (a(ix  ,iy,iz) - a(ix-1,iy,iz))/s1(ix,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
END FUNCTION ddxdn

!*******************************************************************************
FUNCTION ddydn (s2, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a, s2
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddydn', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then
    b(:,1,:) = 0.0
    do iz=1,n(3)
    do iy=2,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = (a(ix,iy  ,iz) - a(ix,iy-1,iz))/s2(ix,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
END FUNCTION ddydn

!*******************************************************************************
FUNCTION ddzdn (s3, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a, s3
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddzdn', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then
    b(:,:,1) = 0.0
    do iz=2,n(3)
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = (a(ix,iy,iz  ) - a(ix,iy,iz-1))/s3(ix,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
END FUNCTION ddzdn

!*******************************************************************************
FUNCTION ddxup (s1, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a, s1
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  n = shape(a)
  if (itimer==0) print *,'stagger: n =', n
  call trace%begin ('stagger_mod::ddxup', itimer=itimer)
  if (n(1) > 1) then
    b(n(1),:,:) = 0.0
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=1,n(1)-1
      b(ix,iy,iz) = (a(ix+1,iy,iz) - a(ix  ,iy,iz))/s1(ix,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
END FUNCTION ddxup

!*******************************************************************************
FUNCTION ddyup (s2, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a, s2
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real, dimension(size(a,1)):: a1
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddyup', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,n(2),:) = 0.0
    do iz=1,n(3)
    do iy=1,n(2)-1
    do ix=1,n(1)
      b(ix,iy,iz) = (a(ix,iy+1,iz) - a(ix,iy,iz))/s2(ix,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
END FUNCTION ddyup

!*******************************************************************************
FUNCTION ddzup (s3, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a, s3
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddzup', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,n(3)) = 0.0
    do iz=1,n(3)-1
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = (a(ix,iy,iz+1) - a(ix,iy,iz  ))/s3(ix,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
END FUNCTION ddzup

!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xdn', itimer=itimer)
  n = shape(a)
  if (n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(1,:,:) = a(1,:,:)
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=2,n(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
    end do
    end do
    end do
  else
    b = a
  end if
!  call trace%end (itimer)
END FUNCTION xdn

!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ydn', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,1,:) = a(:,1,:)
    do iz=1,n(3)
    do iy=2,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
    end do
    end do
    end do
  else
    b = a
  end if
!  call trace%end (itimer)
END FUNCTION ydn

!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::zdn', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,:,1) = a(:,:,1)
    do iz=2,n(3)
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
    end do
    end do
    end do
  else
    b = a
  end if
!  call trace%end (itimer)
END FUNCTION zdn

!*******************************************************************************
FUNCTION xup (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xup', itimer=itimer)
  n = shape(a)
  if (n(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(n(1),:,:) = a(n(1),:,:)
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=1,n(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
    end do
    end do
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
END FUNCTION xup

!*******************************************************************************
FUNCTION yup (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::yup', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,n(2),:) = a(:,n(2),:)
    do iz=1,n(3)
    do iy=1,n(2)-1
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
    end do
    end do
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
END FUNCTION yup

!*******************************************************************************
FUNCTION zup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::zup', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    b(:,:,n(3)) = a(:,:,n(3))
    do iz=1,n(3)-1
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
    end do
    end do
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
END FUNCTION zup

!*******************************************************************************
FUNCTION h2c (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k) = in(i,j,k)*mesh(1)%h2c(i)
      end do
    end do
  end do
END FUNCTION h2c
!===============================================================================
FUNCTION h2f (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k) = in(i,j,k)*mesh(1)%h2f(i)
      end do
    end do
  end do
END FUNCTION h2f
!===============================================================================
FUNCTION h31c (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k) = in(i,j,k)*mesh(1)%h31c(i)
      end do
    end do
  end do
END FUNCTION h31c
!===============================================================================
FUNCTION h31f (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k) = in(i,j,k)*mesh(1)%h31f(i)
      end do
    end do
  end do
END FUNCTION h31f
!===============================================================================
FUNCTION h32c (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k) = in(i,j,k)*mesh(2)%h32c(j)
      end do
    end do
  end do
END FUNCTION h32c
!===============================================================================
FUNCTION h32f (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k) = in(i,j,k)*mesh(2)%h32f(j)
      end do
    end do
  end do
END FUNCTION h32f
!===============================================================================
FUNCTION h3c (mesh, in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  class(mesh_t), dimension(:):: mesh
  integer:: i, j, k, n(3)
  !-----------------------------------------------------------------------------
  n = [size(in,1),size(in,2),size(in,3)]
  do k=1,n(3)
    do j=1,n(2)
      do i=1,n(1)
        out(i,j,k) = in(i,j,k)*mesh(1)%h31c(i)*mesh(2)%h32c(j)
      end do
    end do
  end do
END FUNCTION h3c
!*******************************************************************************
END MODULE
