!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_20
  implicit none
  public
CONTAINS

!*******************************************************************************
SUBROUTINE ddxdn_20 (ds, a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(1,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,20
  do iy=1,20
  do ix=2,20
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxdn_20

!*******************************************************************************
SUBROUTINE ddydn_20 (ds, a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,1,:) = 0.0
  c = 1./ds(2)
  do iz=1,20
  do iy=2,20
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ddydn_20

!*******************************************************************************
SUBROUTINE ddzdn_20 (ds, a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(:,:,1) = 0.0
  c = 1./ds(3)
  do iz=2,20
  do iy=1,20
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE ddzdn_20

!*******************************************************************************
SUBROUTINE ddxup_20 (ds, a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(20,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,20
  do iy=1,20
  do ix=1,20-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxup_20

!*******************************************************************************
SUBROUTINE ddyup_20 (ds, a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,20,:) = 0.0
  c = 1./ds(2)
  do iz=1,20
  do iy=1,20-1
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE ddyup_20

!*******************************************************************************
SUBROUTINE ddzup_20 (ds, a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,:,20) = 0.0
  c = 1./ds(3)
  do iz=1,20-1
  do iy=1,20
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE ddzup_20

!*******************************************************************************
SUBROUTINE xdn_20 (a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(1,:,:) = a(1,:,:)
  do iz=1,20
  do iy=1,20
  do ix=2,20
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE xdn_20

!*******************************************************************************
SUBROUTINE ydn_20 (a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,1,:) = a(:,1,:)
  do iz=1,20
  do iy=2,20
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ydn_20

!*******************************************************************************
SUBROUTINE zdn_20 (a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,1) = a(:,:,1)
  do iz=2,20
  do iy=1,20
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE zdn_20

!*******************************************************************************
SUBROUTINE xup_20 (a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(20,:,:) = a(20,:,:)
  do iz=1,20
  do iy=1,20
  do ix=1,20-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE yup_20 (a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,20,:) = a(:,20,:)
  do iz=1,20
  do iy=1,20-1
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE zup_20 (a, b)
  real, dimension(20,20,20), intent(in):: a
  real, dimension(20,20,20):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,20) = a(:,:,20)
  do iz=1,20-1
  do iy=1,20
  do ix=1,20
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE

END MODULE
