!> =============================================================================
!> $Id$
!> =============================================================================
MODULE stagger_20
  implicit none
  public
CONTAINS

!*******************************************************************************
FUNCTION ddxdn_20 (ds, a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
    b(1,:,:) = 0.0
    c = 1./ds(1)
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=2,size(a,1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
    b = 0.0
END FUNCTION ddxdn_20

!*******************************************************************************
FUNCTION ddydn_20 (ds, a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
    b(:,1,:) = 0.0
    c = 1./ds(2)
    do iz=1,size(a,3)
    do iy=2,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
    end do
    end do
    end do
    b = 0.0
END FUNCTION ddydn_20

!*******************************************************************************
FUNCTION ddzdn_20 (ds, a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
    b(:,:,1) = 0.0
    c = 1./ds(3)
    do iz=2,size(a,3)
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
    end do
    end do
    end do
    b = 0.0
END FUNCTION ddzdn_20

!*******************************************************************************
FUNCTION ddxup_20 (ds, a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
    b(size(a,1),:,:) = 0.0
    c = 1./ds(1)
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=1,size(a,1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
    end do
    end do
    end do
    b = 0.0
END FUNCTION ddxup_20

!*******************************************************************************
FUNCTION ddyup_20 (ds, a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
    b(:,size(a,2),:) = 0.0
    c = 1./ds(2)
    do iz=1,size(a,3)
    do iy=1,size(a,2)-1
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
    end do
    end do
    end do
    b = 0.0
END FUNCTION ddyup_20

!*******************************************************************************
FUNCTION ddzup_20 (ds, a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
    b(:,:,size(a,3)) = 0.0
    c = 1./ds(3)
    do iz=1,size(a,3)-1
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
    end do
    end do
    end do
    b = 0.0
END FUNCTION ddzup_20

!*******************************************************************************
FUNCTION xdn_20 (a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
    c = 0.5
    b(1,:,:) = a(1,:,:)
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=2,size(a,1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
    end do
    end do
    end do
END FUNCTION xdn_20

!*******************************************************************************
FUNCTION ydn_20 (a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
    c = 0.5
    b(:,1,:) = a(:,1,:)
    do iz=1,size(a,3)
    do iy=2,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
    end do
    end do
    end do
END FUNCTION ydn_20

!*******************************************************************************
FUNCTION zdn_20 (a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
    c = 0.5
    b(:,:,1) = a(:,:,1)
    do iz=2,size(a,3)
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
    end do
    end do
    end do
END FUNCTION zdn_20

!*******************************************************************************
FUNCTION xup_20 (a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
    c = 0.5
    b(size(a,1),:,:) = a(size(a,1),:,:)
    do iz=1,size(a,3)
    do iy=1,size(a,2)
    do ix=1,size(a,1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
    end do
    end do
    end do
END FUNCTION

!*******************************************************************************
FUNCTION yup_20 (a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
    c = 0.5
    b(:,size(a,2),:) = a(:,size(a,2),:)
    do iz=1,size(a,3)
    do iy=1,size(a,2)-1
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
    end do
    end do
    end do
END FUNCTION

!*******************************************************************************
FUNCTION zup_20 (a) RESULT (b)
  real, dimension(20,20,20):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
    c = 0.5
    b(:,:,size(a,3)) = a(:,:,size(a,3))
    do iz=1,size(a,3)-1
    do iy=1,size(a,2)
    do ix=1,size(a,1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
    end do
    end do
    end do
END FUNCTION

END MODULE
