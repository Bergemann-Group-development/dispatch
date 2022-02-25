!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_28
  implicit none
  public
CONTAINS

!*******************************************************************************
SUBROUTINE ddxdn_28 (ds, a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(1,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,28
  do iy=1,28
  do ix=2,28
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxdn_28

!*******************************************************************************
SUBROUTINE ddydn_28 (ds, a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,1,:) = 0.0
  c = 1./ds(2)
  do iz=1,28
  do iy=2,28
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ddydn_28

!*******************************************************************************
SUBROUTINE ddzdn_28 (ds, a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(:,:,1) = 0.0
  c = 1./ds(3)
  do iz=2,28
  do iy=1,28
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE ddzdn_28

!*******************************************************************************
SUBROUTINE ddxup_28 (ds, a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(28,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,28
  do iy=1,28
  do ix=1,28-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxup_28

!*******************************************************************************
SUBROUTINE ddyup_28 (ds, a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,28,:) = 0.0
  c = 1./ds(2)
  do iz=1,28
  do iy=1,28-1
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE ddyup_28

!*******************************************************************************
SUBROUTINE ddzup_28 (ds, a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,:,28) = 0.0
  c = 1./ds(3)
  do iz=1,28-1
  do iy=1,28
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE ddzup_28

!*******************************************************************************
SUBROUTINE xdn_28 (a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(1,:,:) = a(1,:,:)
  do iz=1,28
  do iy=1,28
  do ix=2,28
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE xdn_28

!*******************************************************************************
SUBROUTINE ydn_28 (a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,1,:) = a(:,1,:)
  do iz=1,28
  do iy=2,28
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ydn_28

!*******************************************************************************
SUBROUTINE zdn_28 (a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,1) = a(:,:,1)
  do iz=2,28
  do iy=1,28
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE zdn_28

!*******************************************************************************
SUBROUTINE xup_28 (a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(28,:,:) = a(28,:,:)
  do iz=1,28
  do iy=1,28
  do ix=1,28-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE yup_28 (a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,28,:) = a(:,28,:)
  do iz=1,28
  do iy=1,28-1
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE zup_28 (a, b)
  real, dimension(28,28,28), intent(in):: a
  real, dimension(28,28,28):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,28) = a(:,:,28)
  do iz=1,28-1
  do iy=1,28
  do ix=1,28
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE

END MODULE
