!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_16
  implicit none
  public
CONTAINS

!*******************************************************************************
SUBROUTINE ddxdn_16 (ds, a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(1,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,16
  do iy=1,16
  do ix=2,16
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxdn_16

!*******************************************************************************
SUBROUTINE ddydn_16 (ds, a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,1,:) = 0.0
  c = 1./ds(2)
  do iz=1,16
  do iy=2,16
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ddydn_16

!*******************************************************************************
SUBROUTINE ddzdn_16 (ds, a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(:,:,1) = 0.0
  c = 1./ds(3)
  do iz=2,16
  do iy=1,16
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE ddzdn_16

!*******************************************************************************
SUBROUTINE ddxup_16 (ds, a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(16,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,16
  do iy=1,16
  do ix=1,16-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxup_16

!*******************************************************************************
SUBROUTINE ddyup_16 (ds, a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,16,:) = 0.0
  c = 1./ds(2)
  do iz=1,16
  do iy=1,16-1
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE ddyup_16

!*******************************************************************************
SUBROUTINE ddzup_16 (ds, a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,:,16) = 0.0
  c = 1./ds(3)
  do iz=1,16-1
  do iy=1,16
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE ddzup_16

!*******************************************************************************
SUBROUTINE xdn_16 (a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(1,:,:) = a(1,:,:)
  do iz=1,16
  do iy=1,16
  do ix=2,16
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE xdn_16

!*******************************************************************************
SUBROUTINE ydn_16 (a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,1,:) = a(:,1,:)
  do iz=1,16
  do iy=2,16
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ydn_16

!*******************************************************************************
SUBROUTINE zdn_16 (a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,1) = a(:,:,1)
  do iz=2,16
  do iy=1,16
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE zdn_16

!*******************************************************************************
SUBROUTINE xup_16 (a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(16,:,:) = a(16,:,:)
  do iz=1,16
  do iy=1,16
  do ix=1,16-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE yup_16 (a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,16,:) = a(:,16,:)
  do iz=1,16
  do iy=1,16-1
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE zup_16 (a, b)
  real, dimension(16,16,16), intent(in):: a
  real, dimension(16,16,16):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,16) = a(:,:,16)
  do iz=1,16-1
  do iy=1,16
  do ix=1,16
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE

END MODULE
