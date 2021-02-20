!===============================================================================
!> 6th order stagger operators, with self-test procedure
!===============================================================================
MODULE stagger_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  implicit none
  public
  integer, private:: verbose=0
  logical, save:: hardwire=.false.
  type, public:: stagger_t
    integer:: flops=0
    integer:: count=0
    integer:: up(2)=[2,-3]
    integer:: dn(2)=[3,-2]
    integer:: order=6
  contains
    procedure, nopass:: ddxdn1
    procedure, nopass:: ddydn1
    procedure, nopass:: ddzdn1
    procedure, nopass:: ddx => ddx1
    procedure, nopass:: ddy => ddy1
    procedure, nopass:: ddz => ddz1
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
    procedure, nopass:: ddnx
    procedure, nopass:: ddny
    procedure, nopass:: ddnz
    procedure, nopass:: dupx
    procedure, nopass:: dupy
    procedure, nopass:: dupz
    procedure, nopass:: xyzsm
    procedure, nopass:: xmax3
    procedure, nopass:: ymax3
    procedure, nopass:: zmax3
    procedure, nopass:: xmax5
    procedure, nopass:: ymax5
    procedure, nopass:: zmax5
    procedure, nopass:: xsm
    procedure, nopass:: ysm
    procedure, nopass:: zsm
    procedure, nopass:: max3
    procedure, nopass:: stagger_test
  end type
  type(stagger_t):: stagger
CONTAINS

!*******************************************************************************
FUNCTION ddxdn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxdn', itimer=itimer)
  n = shape(a)
  if (n(1) > 1) then
    c = 1./ds(1)
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=2,n(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
    b(1,:,:) = 0.0
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddydn', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then
    c = 1./ds(2)
    do iz=1,n(3)
    do iy=2,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
    end do
    end do
    b(:,1,:) = 0.0
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(3)*(n(2)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzdn', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then
    c = 1./ds(3)
    do iz=2,n(3)
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
    end do
    end do
    end do
    b(:,:,1) = 0.0
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(2)*(n(3)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxup', itimer=itimer)
  n = shape(a)
  if (n(1) > 1) then
    b(n(1),:,:) = 0.0
    c = 1./ds(1)
    !print*,'mk'
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=1,n(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddyup', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,n(2),:) = 0.0
    c = 1./ds(2)
    do iz=1,n(3)
    do iy=1,n(2)-1
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(3)*(n(2)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real(8):: ds(3), c
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzup', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,n(3)) = 0.0
    c = 1./ds(3)
    do iz=1,n(3)-1
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(2)*(n(3)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddzup1

FUNCTION ddx1 (ds, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddx', itimer=itimer)
  n = shape(a)
!  if (hardwire) then
!    if (all(n==[16,16,16])) then
!      call dxdn_16 (ds,a,b)
!    else if (all(n==[20,20,20])) then
!      call dxdn_20 (ds,a,b)
!    else if (all(n==[24,24,24])) then
!      call dxdn_24 (ds,a,b)
!    else if (all(n==[32,32,32])) then
!      call dxdn_32 (ds,a,b)
!    else if (all(n==[36,36,36])) then
!      call dxdn_36 (ds,a,b)
!    end if
!  else if (n(1) > 1) then
  if (n(1) > 1) then
    b(   1,:,:) = 0.0
    b(n(1),:,:) = 0.0
    c = 0.5/ds(1)
    do iz=1,n(3)
    do iy=1,n(2)
    !$vector always assert
    do ix=2,n(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddx1

!*******************************************************************************
FUNCTION ddy1 (ds, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddy', itimer=itimer)
  n = shape(a)
!  if (hardwire) then
!    if (all(n==[16,16,16])) then
!      call dxdn_16 (ds,a,b)
!    else if (all(n==[20,20,20])) then
!      call dxdn_20 (ds,a,b)
!    else if (all(n==[24,24,24])) then
!      call dxdn_24 (ds,a,b)
!    else if (all(n==[32,32,32])) then
!      call dxdn_32 (ds,a,b)
!    else if (all(n==[36,36,36])) then
!      call dxdn_36 (ds,a,b)
!    end if
!  else if (n(1) > 1) then
  if (n(2) > 1) then
    b(:,   1,:) = 0.0
    b(:,n(2),:) = 0.0
    c = 0.5/ds(2)
    do iz=1,n(3)
    do iy=2,n(2)-1
    !$vector always assert
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy-1,iz))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddy1

!*******************************************************************************
FUNCTION ddz1 (ds, a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real(8):: ds(3), c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  call trace%begin ('stagger_mod::ddz', itimer=itimer)
  n = shape(a)
!  if (hardwire) then
!    if (all(n==[16,16,16])) then
!      call dxdn_16 (ds,a,b)
!    else if (all(n==[20,20,20])) then
!      call dxdn_20 (ds,a,b)
!    else if (all(n==[24,24,24])) then
!      call dxdn_24 (ds,a,b)
!    else if (all(n==[32,32,32])) then
!      call dxdn_32 (ds,a,b)
!    else if (all(n==[36,36,36])) then
!      call dxdn_36 (ds,a,b)
!    end if
!  else if (n(1) > 1) then
  if (n(3) > 1) then
    b(:,:,   1) = 0.0
    b(:,:,n(3)) = 0.0
    c = 0.5/ds(3)
    do iz=2,n(3)-1
    do iy=1,n(2)
    !$vector always assert
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz-1))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddz1

!*******************************************************************************
FUNCTION xdn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
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
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION xdn1

!*******************************************************************************
FUNCTION ydn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
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
  stagger%flops = stagger%flops + 2*n(1)*n(3)*(n(2)-1)
  stagger%count = stagger%count + 1
END FUNCTION ydn1

!*******************************************************************************
FUNCTION zdn1 (a) RESULT (b)
  real, dimension(:,:,:):: a
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
  stagger%flops = stagger%flops + 2*n(1)*n(2)*(n(3)-1)
  stagger%count = stagger%count + 1
END FUNCTION zdn1

!*******************************************************************************
FUNCTION xup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
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
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION xup1

!*******************************************************************************
FUNCTION yup1 (a) RESULT (b)
  real, dimension(:,:,:):: a
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
  stagger%flops = stagger%flops + 2*n(1)*n(3)*(n(2)-1)
  stagger%count = stagger%count + 1
END FUNCTION yup1

!*******************************************************************************
FUNCTION zup1 (a) RESULT (b)
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
  stagger%flops = stagger%flops + 2*n(1)*n(2)*(n(3)-1)
  stagger%count = stagger%count + 1
END FUNCTION zup1

!*******************************************************************************
FUNCTION ddxdn (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(1) > 5) then
    d = ds(1)
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(    1:3,:,:) = 0.0
    b(n(1)-1:,:,:) = 0.0
    do k=1,n(3)
    do j=1,n(2)
    do i=4,n(1)-2
      b(i  ,j,k) = (c3*(a(i+2,j,k)-a(i-3,j,k)) &
                  + c2*(a(i+1,j,k)-a(i-2,j,k))) &
                  + c1*(a(i  ,j,k)-a(i-1,j,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  stagger%flops = stagger%flops + 8*n(2)*n(3)*(n(1)-5)
  stagger%count = stagger%count + 1
END FUNCTION ddxdn

!*******************************************************************************
FUNCTION ddydn (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(2) > 5) then
    d = ds(2)
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,    1:3,:) = 0.0
    b(:,n(2)-1:,:) = 0.0
    do k=1,n(3)
    do j=4,n(2)-2
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j+2,k)-a(i,j-3,k)) &
                + c2*(a(i,j+1,k)-a(i,j-2,k))) &
                + c1*(a(i,j  ,k)-a(i,j-1,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(3)*(n(2)-5)
  stagger%count = stagger%count + 1
END FUNCTION ddydn

!*******************************************************************************
FUNCTION ddzdn (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(3) > 5) then
    d = ds(3)
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,:,    1:3) = 0.0
    b(:,:,n(3)-1:) = 0.0
    do k=4,n(3)-2
    do j=1,n(2)
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j,k+2)-a(i,j,k-3)) &
                + c2*(a(i,j,k+1)-a(i,j,k-2))) &
                + c1*(a(i,j,k  )-a(i,j,k-1))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(2)*(n(3)-5)
  stagger%count = stagger%count + 1
END FUNCTION ddzdn

!*******************************************************************************
FUNCTION ddxup (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  d = ds(1)
  if (n(1) > 5) then
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(    1:2,:,:) = 0.0
    b(n(1)-2:,:,:) = 0.0
    do k=1,n(3)
    do j=1,n(2)
    do i=3,n(1)-3
      b(i  ,j,k) = (c3*(a(i+3,j,k)-a(i-2,j,k)) &
                  + c2*(a(i+2,j,k)-a(i-1,j,k))) &
                  + c1*(a(i+1,j,k)-a(i  ,j,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  stagger%flops = stagger%flops + 8*n(2)*n(3)*(n(1)-5)
  stagger%count = stagger%count + 1
END FUNCTION ddxup

!*******************************************************************************
FUNCTION ddyup (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  d = ds(2)
  if (n(2) > 5) then
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,    1:2,:) = 0.0
    b(:,n(2)-2:,:) = 0.0
    do k=1,n(3)
    do j=3,n(2)-3
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j+3,k)-a(i,j-2,k)) &
                + c2*(a(i,j+2,k)-a(i,j-1,k))) &
                + c1*(a(i,j+1,k)-a(i,j  ,k))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(3)*(n(2)-5)
  stagger%count = stagger%count + 1
END FUNCTION ddyup

!*******************************************************************************
FUNCTION ddzup (ds, a) RESULT (b)
  real(8) :: ds(3)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  d = ds(3)
  if (n(3) > 5) then
    c3 = (-1.+(3.**5-3.)/(3.**3-3.))/(5.**5-5.-5.*(3.**5-3))
    c2 = (-1.-120.*c3)/24.
    c1 = (1.-3.*c2-5.*c3)/d
    c2 = c2/d
    c3 = c3/d
    b(:,:,    1:2) = 0.0
    b(:,:,n(3)-2:) = 0.0
    do k=3,n(3)-3
    do j=1,n(2)
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j,k+3)-a(i,j,k-2)) &
                + c2*(a(i,j,k+2)-a(i,j,k-1))) &
                + c1*(a(i,j,k+1)-a(i,j,k  ))
    end do
    end do
    end do
  else
    b = 0.0
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(2)*(n(3)-5)
  stagger%count = stagger%count + 1
END FUNCTION ddzup

!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(    1:3,:,:) = a(    1:3,:,:)
    b(n(1)-2:,:,:) = a(n(1)-2:,:,:)
    do k=1,n(3)
    do j=1,n(2)
    do i=4,n(1)-2
      b(i  ,j,k) = (c3*(a(i+2,j,k)+a(i-3,j,k)) &
                  + c2*(a(i+1,j,k)+a(i-2,j,k))) &
                  + c1*(a(i  ,j,k)+a(i-1,j,k))
    end do
    end do
    end do
  else
    b = a
  endif
  stagger%flops = stagger%flops + 8*n(2)*n(3)*(n(1)-5)
  stagger%count = stagger%count + 1
END FUNCTION xdn

!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(2) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,    1:3,:) = a(:,    1:3,:)
    b(:,n(2)-2:,:) = a(:,n(2)-2:,:)
    do k=1,n(3)
    do j=4,n(2)-2
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j+2,k)+a(i,j-3,k)) &
                + c2*(a(i,j+1,k)+a(i,j-2,k))) &
                + c1*(a(i,j  ,k)+a(i,j-1,k))
    end do
    end do
    end do
  else
    b = a
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(3)*(n(2)-5)
  stagger%count = stagger%count + 1
END FUNCTION ydn

!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(3) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,:,    1:3) = a(:,:,    1:3)
    b(:,:,n(3)-2:) = a(:,:,n(3)-2:)
    do k=4,n(3)-2
    do j=1,n(2)
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j,k+2)+a(i,j,k-3)) &
                + c2*(a(i,j,k+1)+a(i,j,k-2))) &
                + c1*(a(i,j,k  )+a(i,j,k-1))
    end do
    end do
    end do
  else
    b = a
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(2)*(n(3)-5)
  stagger%count = stagger%count + 1
END FUNCTION zdn

!*******************************************************************************
FUNCTION xup (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(1) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(    1:2,:,:) = a(    1:2,:,:)
    b(n(1)-2:,:,:) = a(n(1)-2:,:,:)
    do k=1,n(3)
    do j=1,n(2)
    do i=3,n(1)-3
      b(i  ,j,k) = (c3*(a(i+3,j,k)+a(i-2,j,k)) &
                  + c2*(a(i+2,j,k)+a(i-1,j,k))) &
                  + c1*(a(i+1,j,k)+a(i  ,j,k))
    end do
    end do
    end do
  else
    b = a
  endif
  stagger%flops = stagger%flops + 8*n(2)*n(3)*(n(1)-5)
  stagger%count = stagger%count + 1
END FUNCTION xup

!*******************************************************************************
FUNCTION yup (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(2) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,    1:2,:) = a(:,    1:2,:)
    b(:,n(2)-2:,:) = a(:,n(2)-2:,:) 
    do k=1,n(3)
    do j=3,n(2)-3
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j+3,k)+a(i,j-2,k)) &
                + c2*(a(i,j+2,k)+a(i,j-1,k))) &
                + c1*(a(i,j+1,k)+a(i,j  ,k))
    end do
    end do
    end do
  else
    b = a
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(3)*(n(2)-5)
  stagger%count = stagger%count + 1
END FUNCTION yup

!*******************************************************************************
FUNCTION zup (a) RESULT (b)
  real    :: a(:,:,:)
  real    :: b(size(a,1),size(a,2),size(a,3))
  !.............................................................................
  integer :: i, j, k, n(3)
  real    :: c1, c2, c3, d
  !-----------------------------------------------------------------------------
  n = shape(a)
  if (n(3) > 5) then
    c3 = 3./256.
    c2 = -25./256.
    c1 = .5-c2-c3
    b(:,:,    1:2) = a(:,:,    1:2)
    b(:,:,n(3)-2:) = a(:,:,n(3)-2:)
    do k=3,n(3)-3
    do j=1,n(2)
    do i=1,n(1)
      b(i,j,k) = (c3*(a(i,j,k+3)+a(i,j,k-2)) &
                + c2*(a(i,j,k+2)+a(i,j,k-1))) &
                + c1*(a(i,j,k+1)+a(i,j,k  ))
    end do
    end do
    end do
  else
    b = a
  endif
  stagger%flops = stagger%flops + 8*n(1)*n(2)*(n(3)-5)
  stagger%count = stagger%count + 1
END FUNCTION zup

!*******************************************************************************
FUNCTION xsm (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xsm', itimer=itimer)
  n = size(a,1)
  if (n > 2) then
    c = 0.25
    do i=2,n-1
      b(i,:,:) = (a(i-1,:,:) + a(i+1,:,:))*c +  a(i,:,:)*(1.-2.*c)
    end do
    do i=1,n,n-1
      b(i,:,:) = a(i,:,:)
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
  !stagger%flops = stagger%flops + 4*n*n*(n-2)
  !stagger%count = stagger%count + 1
END FUNCTION

!*******************************************************************************
FUNCTION ysm (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ysm', itimer=itimer)
  n = size(a,2)
  if (n > 2) then
    c = 0.25
    do i=2,n-1
      b(:,i,:) = (a(:,i-1,:) + a(:,i+1,:))*c +  a(:,i,:)*(1.-2.*c)
    end do
    do i=1,n,n-1
      b(:,i,:)  = a(:,i,:) 
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
  !stagger%flops = stagger%flops + 4*n*n*(n-2)
  !stagger%count = stagger%count + 1
END FUNCTION

!*******************************************************************************
FUNCTION zsm (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::zsm', itimer=itimer)
  n = size(a,3)
  if (n > 2) then
    c = 0.25
    do i=2,n-1
      b(:,:,i) = (a(:,:,i-1) + a(:,:,i+1))*c +  a(:,:,i)*(1.-2.*c)
    end do
    do i=1,n,n-1
      b(:,:,i) = a(:,:,i)
    end do
  else
    b = a
  endif
!  call trace%end (itimer)
  !stagger%flops = stagger%flops + 4*n*n*(n-2)
  !stagger%count = stagger%count + 1
END FUNCTION

!*******************************************************************************
FUNCTION xyzsm(a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer, save:: itimer=0
  call trace%begin ('stagger_mod::xyzsm', itimer=itimer)
  b = zsm(ysm(xsm(a)))
  call trace%end (itimer)
END FUNCTION


!*******************************************************************************
FUNCTION ddnx (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxdn', itimer=itimer)
  n = shape(a)
  if (n(1) > 1) then
    b(1,:,:) = 0.0
    do iz=1,n(3)
    do iy=1,n(2)
    !$vector always assert
    do ix=2,n(1)
      b(ix,iy,iz) = a(ix  ,iy,iz) - a(ix-1,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddnx

!*******************************************************************************
FUNCTION dupx (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddxup', itimer=itimer)
  n = shape(a)
  if (n(1) > 1) then
    b(n(1),:,:) = 0.0
    do iz=1,n(3)
    do iy=1,n(2)
    !$vector always assert
    do ix=1,n(1)-1
      b(ix,iy,iz) = a(ix+1,iy,iz) - a(ix  ,iy,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(3)*n(2)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION dupx

!*******************************************************************************
FUNCTION dupy (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddyup', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,n(2),:) = 0.0
    do iz=1,n(3)
    do iy=1,n(2)-1
    !$vector always assert
    do ix=1,n(1)
      b(ix,iy,iz) = a(ix,iy+1,iz) - a(ix,iy  ,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(3)*(n(2)-1)
  stagger%count = stagger%count + 1
END FUNCTION dupy

!*******************************************************************************
FUNCTION dupz (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzup', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    b(:,:,n(3)) = 0.0
    do iz=1,n(3)-1
    do iy=1,n(2)
    !$vector always assert
    do ix=1,n(1)
      b(ix,iy,iz) = a(ix,iy,iz+1) - a(ix,iy,iz  )
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(2)*(n(3)-1)
  stagger%count = stagger%count + 1
END FUNCTION dupz

!*******************************************************************************
FUNCTION ddny (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddydn', itimer=itimer)
  n = shape(a)
  if (n(2) > 1) then
    b(:,1,:) = 0.0
    do iz=1,n(3)
    do iy=2,n(2)
    !$vector always assert
    do ix=1,n(1)
      b(ix,iy,iz) = a(ix,iy  ,iz) - a(ix,iy-1,iz)
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(3)*(n(2)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddny

!*******************************************************************************
FUNCTION ddnz (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::ddzdn', itimer=itimer)
  n = shape(a)
  if (n(3) > 1) then
    b(:,:,1) = 0.0
    do iz=2,n(3)
    do iy=1,n(2)
    !$vector always assert
    do ix=1,n(1)
      b(ix,iy,iz) = a(ix,iy,iz  ) - a(ix,iy,iz-1)
    end do
    end do
    end do
  else
    b = 0.0
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 2*n(1)*n(2)*(n(3)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddnz

!*******************************************************************************
FUNCTION max3(a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer, save:: itimer=0
  call trace%begin ('stagger_mod::max3', itimer=itimer)
  b = zmax3(ymax3(xmax3(a)))
  call trace%end (itimer)
END FUNCTION

!*******************************************************************************
FUNCTION xmax3 (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xsm', itimer=itimer)
  n = size(a,1)
  if (n > 2) then
    do i=2,n-1
      b(i,:,:) = max(a(i-1,:,:),a(i,:,:),a(i+1,:,:))
    end do
    b(1,:,:)  = b(2,:,:)
    b(n,:,:)  = b(n-1,:,:)
  else
    b = a
  endif

!  call trace%end (itimer)
  !stagger%flops = stagger%flops + 4*n*n*(n-2)
  !stagger%count = stagger%count + 1
END FUNCTION

!*******************************************************************************
FUNCTION ymax3 (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xsm', itimer=itimer)
  n = size(a,2)
  if (n > 2) then
    do i=2,n-1
      b(:,i,:) = max(a(:,i-1,:),a(:,i,:),a(:,i+1,:))
    end do
    b(:,1,:)  = b(:,2,:)
    b(:,n,:)  = b(:,n-1,:)
  else
    b = a
  endif
!  call trace%end (itimer)
  !stagger%flops = stagger%flops + 4*n*n*(n-2)
  !stagger%count = stagger%count + 1
END FUNCTION ymax3

!*******************************************************************************
FUNCTION zmax3 (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xsm', itimer=itimer)
  n = size(a,3)
  if (n > 2) then
    do i=2,n-1
      b(:,:,i) = max(a(:,:,i-1),a(:,:,i),a(:,:,i+1))
    end do
    b(:,:,1)  = b(:,:,2)
    b(:,:,n)  = b(:,:,n-1)
  else
    b = a
  endif
!  call trace%end (itimer)
  !stagger%flops = stagger%flops + 4*n*n*(n-2)
  !stagger%count = stagger%count + 1
END FUNCTION zmax3

!*******************************************************************************
FUNCTION max5(a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer, save:: itimer=0
  call trace%begin ('stagger_mod::max3', itimer=itimer)
  b = zmax5(ymax5(xmax5(a)))
  call trace%end (itimer)
END FUNCTION

!*******************************************************************************
FUNCTION xmax5 (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xsm', itimer=itimer)
  n = size(a,1)
  if (n > 5) then
    do i=3,n-2
      b(i,:,:) = max(a(i-2,:,:),a(i-1,:,:),a(i,:,:),a(i+1,:,:),a(i+2,:,:))
    end do
    b(2,:,:) = max(a(2,:,:),a(3,:,:),a(4,:,:))
    b(1,:,:) = max(a(1,:,:),a(2,:,:),a(3,:,:))
    b(n-1,:,:) = max(a(n-1,:,:),a(n-2,:,:),a(n-3,:,:))
    b(n  ,:,:) = max(a(n  ,:,:),a(n-1,:,:),a(n-2,:,:))
  else
    b = a
  endif

!  call trace%end (itimer)
  stagger%flops = stagger%flops + 4*n*n*(n-2)
  stagger%count = stagger%count + 1
END FUNCTION

!*******************************************************************************
FUNCTION ymax5 (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xsm', itimer=itimer)
  n = size(a,2)
  if (n > 5) then
    do i=3,n-2
      b(:,i,:) = max(a(:,i-2,:),a(:,i-1,:),a(:,i,:),a(:,i+1,:),a(:,i+2,:))
    end do
    b(:,2,:) = max(a(:,2,:),a(:,3,:),a(:,4,:))
    b(:,1,:) = max(a(:,1,:),a(:,2,:),a(:,3,:))
    b(:,n-1,:) = max(a(:,n-1,:),a(:,n-2,:),a(:,n-3,:))
    b(:,n  ,:) = max(a(:,n  ,:),a(:,n-1,:),a(:,n-2,:))
  else
    b = a
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 4*n*n*(n-2)
  stagger%count = stagger%count + 1
END FUNCTION ymax5

!*******************************************************************************
FUNCTION zmax5 (a) RESULT (b)
  real, dimension(:,:,:), intent(in):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: i, n
  integer, save:: itimer=0
!...............................................................................
!  call trace%begin ('stagger_mod::xsm', itimer=itimer)
  n = size(a,3)
  if (n > 5) then
    do i=3,n-2
      b(:,:,i) = max(a(:,:,i-2),a(:,:,i-1),a(:,:,i),a(:,:,i+1),a(:,:,i+2))
    end do
    b(:,:,2) = max(a(:,:,2),a(:,:,3),a(:,:,4))
    b(:,:,1) = max(a(:,:,1),a(:,:,2),a(:,:,3))
    b(:,:,n-1) = max(a(:,:,n-1),a(:,:,n-2),a(:,:,n-3))
    b(:,:,n  ) = max(a(:,:,n ),a(:,:,n-1),a(:,:,n-2))
  else
    b = a
  endif
!  call trace%end (itimer)
  stagger%flops = stagger%flops + 4*n*n*(n-2)
  stagger%count = stagger%count + 1
END FUNCTION zmax5

!===============================================================================
!> Unit test of 6th/5th order operators
!===============================================================================
SUBROUTINE stagger_test
  USE mpi_mod
  ! --------------------------------------------------------------------
  integer       :: n(3)
  real, pointer :: a(:,:,:)
  real, pointer :: b(:,:,:)
  real, pointer :: c(:,:,:)
  real(8)       :: ds(3), x, y, z, d
  integer       :: i
  logical       :: ok, allok
  real, parameter :: eps = 2e-6
  !-----------------------------------------------------------------------------
  n = [10,11,12]
  ds = 1d0/n
  allocate (a(n(1),n(2),n(3)))
  allocate (b(n(1),n(2),n(3)))
  allok = .true.
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,n(1)
    x = i*ds(1)
    a(i,:,:) = x**2
  end do
  b = ddxdn1(ds,a)
  ok = .true.
  do i=4,n(1)-2
    x = (i-0.5d0)*ds(1)
    d = 2.*x
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddxdn1:', ok
  b = ddxup1(ds,a)
  ok = .true.
  do i=3,n(1)-3
    x = (i+0.5d0)*ds(1)
    d = 2.*x
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddxup1:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,n(2)
    y = i*ds(2)
    a(:,i,:) = y**2
  end do
  b = ddydn1(ds,a)
  ok = .true.
  do i=4,n(2)-2
    y = (i-0.5d0)*ds(2)
    d = 2.*y
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddydn1:', ok
  b = ddyup1(ds,a)
  ok = .true.
  do i=3,n(2)-3
    y = (i+.5d0)*ds(2)
    d = 2.*y
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddyup1:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,n(3)
    z = i*ds(3)
    a(:,:,i) = z**2
  end do
  b = ddzdn1(ds,a)
  ok = .true.
  do i=4,n(3)-2
    z = (i-0.5d0)*ds(3)
    d = 2.*z
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddzdn1:', ok
  b = ddzup1(ds,a)
  ok = .true.
  do i=3,n(3)-3
    z = (i+0.5d0)*ds(3)
    d = 2.*z
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddzup1:', ok
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,n(1)
    x = i*ds(1)
    a(i,:,:) = x
  end do
  b = xdn1(a)
  ok = .true.
  do i=4,n(1)-2
    x = (i-0.5d0)*ds(1)
    d = x
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'xdn1:', ok
  b = xup1(a)
  ok = .true.
  do i=3,n(1)-3
    x = (i+0.5d0)*ds(1)
    d = x
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'xup1:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,n(2)
    y = i*ds(2)
    a(:,i,:) = y
  end do
  b = ydn1(a)
  ok = .true.
  do i=4,n(2)-2
    y = (i-0.5d0)*ds(2)
    d = y
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ydn1:', ok
  b = yup1(a)
  ok = .true.
  do i=3,n(2)-3
    y = (i+.5d0)*ds(2)
    d = y
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'yup1:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,n(3)
    z = i*ds(3)
    a(:,:,i) = z
  end do
  b = zdn1(a)
  ok = .true.
  do i=4,n(3)-2
    z = (i-0.5d0)*ds(3)
    d = z
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'zdn1:', ok
  b = zup1(a)
  ok = .true.
  do i=3,n(3)-3
    z = (i+0.5d0)*ds(3)
    d = z
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'zup1:', ok
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,n(1)
    x = i*ds(1)
    a(i,:,:) = x**5
  end do
  b = xdn(a)
  ok = .true.
  do i=4,n(1)-2
    x = (i-0.5d0)*ds(1)
    d = x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'xdn:', ok
  b = xup(a)
  ok = .true.
  do i=3,n(1)-3
    x = (i+0.5d0)*ds(1)
    d = x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'xup:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,n(2)
    x = i*ds(2)
    a(:,i,:) = x**5
  end do
  b = ydn(a)
  ok = .true.
  do i=4,n(2)-2
    y = (i-0.5d0)*ds(2)
    d = y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ydn:', ok
  b = yup(a)
  ok = .true.
  do i=3,n(2)-3
    y = (i+.5d0)*ds(2)
    d = y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'yup:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,n(3)
    x = i*ds(3)
    a(:,:,i) = x**5
  end do
  b = zdn(a)
  ok = .true.
  do i=4,n(3)-2
    z = (i-0.5d0)*ds(3)
    d = z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'zdn:', ok
  b = zup(a)
  ok = .true.
  do i=3,n(3)-3
    z = (i+0.5d0)*ds(3)
    d = z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'zup:', ok
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,n(1)
    x = i*ds(1)
    a(i,:,:) = 1d0 + x**6
  end do
  b = ddxdn(ds,a)
  ok = .true.
  do i=4,n(1)-2
    x = (i-0.5d0)*ds(1)
    d = 6d0*x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddxdn:', ok
  b = ddxup(ds,a)
  ok = .true.
  do i=3,n(1)-3
    x = (i+0.5d0)*ds(1)
    d = 6d0*x**5
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddxup:', ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,n(2)
    x = i*ds(2)
    a(:,i,:) = 1d0 + x**6
  end do
  b = ddydn(ds,a)
  ok = .true.
  do i=4,n(2)-2
    y = (i-0.5d0)*ds(2)
    d = 6d0*y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddydn:', ok
  b = ddyup(ds,a)
  ok = .true.
  do i=3,n(2)-3
    y = (i+.5d0)*ds(2)
    d = 6d0*y**5
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddyup:', ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,n(3)
    x = i*ds(3)
    a(:,:,i) = 1d0 + x**6
  end do
  b = ddzdn(ds,a)
  ok = .true.
  do i=4,n(3)-2
    z = (i-0.5d0)*ds(3)
    d = 6d0*z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddzdn:', ok
  b = ddzup(ds,a)
  ok = .true.
  do i=3,n(3)-3
    z = (i+0.5d0)*ds(3)
    d = 6d0*z**5
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  allok = allok .and. ok
  if (.not.ok) write (stderr,*) 'ddzup:', ok
  write (stderr,*) io%hl
  if (omp%master) then
    if (allok) then
      write (stdout,*) '6th/5th order stagger operator test passed'
      write (stdout,*) io%hl
    else
      write (stdout,*) '6th/5th order stagger operator test failed!'
      write (stdout,*)
      call mpi%abort('stagger_test')
    end if
  end if
  !-----------------------------------------------------------------------------
  deallocate (a, b)
END SUBROUTINE stagger_test

END MODULE stagger_mod
