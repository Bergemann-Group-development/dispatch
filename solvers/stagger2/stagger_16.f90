!> =============================================================================
!> $Id$
!> =============================================================================
MODULE stagger_16
  implicit none
  public
CONTAINS

!*******************************************************************************
FUNCTION ddxdn_16 (ds, a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION ddxdn_16

!*******************************************************************************
FUNCTION ddydn_16 (ds, a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION ddydn_16

!*******************************************************************************
FUNCTION ddzdn_16 (ds, a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION ddzdn_16

!*******************************************************************************
FUNCTION ddxup_16 (ds, a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION ddxup_16

!*******************************************************************************
FUNCTION ddyup_16 (ds, a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION ddyup_16

!*******************************************************************************
FUNCTION ddzup_16 (ds, a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION ddzup_16

!*******************************************************************************
FUNCTION xdn_16 (a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION xdn_16

!*******************************************************************************
FUNCTION ydn_16 (a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION ydn_16

!*******************************************************************************
FUNCTION zdn_16 (a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
END FUNCTION zdn_16

!*******************************************************************************
FUNCTION xup_16 (a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
FUNCTION yup_16 (a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
FUNCTION zup_16 (a) RESULT (b)
  real, dimension(16,16,16):: a
  real, dimension(16,16,16):: b
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
