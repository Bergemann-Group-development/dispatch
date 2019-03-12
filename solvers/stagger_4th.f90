!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_mod
  USE io_mod
  USE trace_mod
  implicit none
  public
  integer, private:: verbose=0
  logical, save:: hardwire=.false.
  type, public:: stagger_t
    integer:: flops=0
    integer:: count=0
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
    procedure, nopass:: xyzsm
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
    b(1,:,:) = 0.0
    c = 1./ds(1)
    do iz=1,n(3)
    do iy=1,n(2)
    do ix=2,n(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
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
    b(:,1,:) = 0.0
    c = 1./ds(2)
    do iz=1,n(3)
    do iy=2,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
    end do
    end do
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
    b(:,:,1) = 0.0
    c = 1./ds(3)
    do iz=2,n(3)
    do iy=1,n(2)
    do ix=1,n(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
    end do
    end do
    end do
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
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  real(8) :: ds(3)
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddxdn

!*******************************************************************************
FUNCTION ddydn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real    :: c1, c2
  real(8) :: ds(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(3)*(n(2)-3)
  stagger%count = stagger%count + 1
END FUNCTION ddydn

!*******************************************************************************
FUNCTION ddzdn (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  real(8) :: ds(3)
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(2)*(n(3)-3)
  stagger%count = stagger%count + 1
END FUNCTION ddzdn

!*******************************************************************************
FUNCTION ddxup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
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
  stagger%flops = stagger%flops + 5*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION ddxup

!*******************************************************************************
FUNCTION ddyup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real    :: c1, c2
  real(8) :: ds(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(3)*(n(2)-3)
  stagger%count = stagger%count + 1
END FUNCTION ddyup

!*******************************************************************************
FUNCTION ddzup (ds, a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, n(3)
  real    :: c1, c2
  real(8) :: ds(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(2)*(n(3)-3)
  stagger%count = stagger%count + 1
END FUNCTION ddzup

!*******************************************************************************
FUNCTION xdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION xdn

!*******************************************************************************
FUNCTION ydn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(3)*(n(2)-3)
  stagger%count = stagger%count + 1
END FUNCTION ydn

!*******************************************************************************
FUNCTION zdn (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(2)*(n(3)-3)
  stagger%count = stagger%count + 1
END FUNCTION zdn

!*******************************************************************************
FUNCTION xup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(2)*n(3)*(n(1)-1)
  stagger%count = stagger%count + 1
END FUNCTION

!*******************************************************************************
FUNCTION yup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(3)*(n(2)-3)
  stagger%count = stagger%count + 1
END FUNCTION

!*******************************************************************************
FUNCTION zup (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c1, c2
  integer :: ix, iy, iz, n(3)
!...............................................................................
  n = shape(a)
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
  stagger%flops = stagger%flops + 5*n(1)*n(2)*(n(3)-3)
  stagger%count = stagger%count + 1
END FUNCTION

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
  if (n > 1) then
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
  stagger%flops = stagger%flops + 4*n*n*(n-2)
  stagger%count = stagger%count + 1
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
  if (n > 1) then
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
  stagger%flops = stagger%flops + 4*n*n*(n-2)
  stagger%count = stagger%count + 1
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
  if (n > 1) then
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
  stagger%flops = stagger%flops + 4*n*n*(n-2)
  stagger%count = stagger%count + 1
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

!===============================================================================
!> Unit test of 4th/3th order operators
!===============================================================================
SUBROUTINE stagger_test
  integer       :: m(3)
  real, pointer :: a(:,:,:)
  real, pointer :: b(:,:,:)
  real, pointer :: c(:,:,:)
  real(8)       :: ds(3), x, y, z, d
  integer       :: i
  logical       :: ok, allok
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
  allok = .true.
  do i=3,m(1)-1
    x = (i-0.5d0)*ds(1)
    d = x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  if (.not.ok) print *,'xdn:', ok
  allok = allok.and.ok
  b = xup(a)
  ok = .true.
  do i=2,m(1)-2
    x = (i+0.5d0)*ds(1)
    d = x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  if (.not.ok) print *,'xup:', ok
  allok = allok.and.ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,m(2)
    x = i*ds(2)
    a(:,i,:) = x**3
  end do
  b = ydn(a)
  ok = .true.
  do i=3,m(2)-1
    y = (i-0.5d0)*ds(2)
    d = y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  if (.not.ok) print *,'ydn:', ok
  allok = allok.and.ok
  b = yup(a)
  ok = .true.
  do i=2,m(2)-2
    y = (i+.5d0)*ds(2)
    d = y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  if (.not.ok) print *,'yup:', ok
  allok = allok.and.ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,m(3)
    x = i*ds(3)
    a(:,:,i) = x**3
  end do
  b = zdn(a)
  ok = .true.
  do i=3,m(3)-1
    z = (i-0.5d0)*ds(3)
    d = z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  if (.not.ok) print *,'zdn:', ok
  allok = allok.and.ok
  b = zup(a)
  ok = .true.
  do i=2,m(3)-2
    z = (i+0.5d0)*ds(3)
    d = z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  if (.not.ok) print *,'zup:', ok
  allok = allok.and.ok
  !-----------------------------------------------------------------------------
  ! x-operators
  !-----------------------------------------------------------------------------
  do i=1,m(1)
    x = i*ds(1)
    a(i,:,:) = 1d0 + x**4
  end do
  b = ddxdn(ds,a)
  ok = .true.
  do i=3,m(1)-1
    x = (i-0.5d0)*ds(1)
    d = 4d0*x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  if (.not.ok) print *,'ddxdn:', ok
  allok = allok.and.ok
  b = ddxup(ds,a)
  ok = .true.
  do i=2,m(1)-2
    x = (i+0.5d0)*ds(1)
    d = 4d0*x**3
    ok = ok .and. (maxval(abs(b(i,:,:)-d)) < eps)
  end do
  if (.not.ok) print *,'ddxup:', ok
  allok = allok.and.ok
  !-----------------------------------------------------------------------------
  ! y-operators
  !-----------------------------------------------------------------------------
  do i=1,m(2)
    x = i*ds(2)
    a(:,i,:) = 1d0 + x**4
  end do
  b = ddydn(ds,a)
  ok = .true.
  do i=3,m(2)-1
    y = (i-0.5d0)*ds(2)
    d = 4d0*y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  if (.not.ok) print *,'ddydn:', ok
  allok = allok.and.ok
  b = ddyup(ds,a)
  ok = .true.
  do i=2,m(2)-2
    y = (i+.5d0)*ds(2)
    d = 4d0*y**3
    ok = ok .and. (maxval(abs(b(:,i,:)-d)) < eps)
  end do
  if (.not.ok) print *,'ddyup:', ok
  allok = allok.and.ok
  !-----------------------------------------------------------------------------
  ! z-operators
  !-----------------------------------------------------------------------------
  do i=1,m(3)
    x = i*ds(3)
    a(:,:,i) = 1d0 + x**4
  end do
  b = ddzdn(ds,a)
  ok = .true.
  do i=3,m(3)-1
    z = (i-0.5d0)*ds(3)
    d = 4d0*z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  if (.not.ok) print *,'ddzdn:', ok
  allok = allok.and.ok
  b = ddzup(ds,a)
  ok = .true.
  do i=2,m(3)-2
    z = (i+0.5d0)*ds(3)
    d = 4d0*z**3
    ok = ok .and. (maxval(abs(b(:,:,i)-d)) < eps)
  end do
  if (.not.ok) print *,'ddzup:', ok
  allok = allok.and.ok
  print *,io%hl
  if (allok) then
    print *,'4th/3rd order stagger operator test passed'
    print *,io%hl
  else
    print *,'4th/3rd order stagger operator test failed!'
    print *,io%hl
    error stop 'stagger_test'
  end if
  !-----------------------------------------------------------------------------
  deallocate (a, b)
END SUBROUTINE stagger_test

END MODULE
