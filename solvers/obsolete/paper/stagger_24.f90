!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_24
  implicit none
  public
CONTAINS

!*******************************************************************************
SUBROUTINE ddxdn_24 (ds, a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(1,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,24
  do iy=1,24
  do ix=2,24
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxdn_24

!*******************************************************************************
SUBROUTINE ddydn_24 (ds, a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,1,:) = 0.0
  c = 1./ds(2)
  do iz=1,24
  do iy=2,24
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ddydn_24

!*******************************************************************************
SUBROUTINE ddzdn_24 (ds, a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(:,:,1) = 0.0
  c = 1./ds(3)
  do iz=2,24
  do iy=1,24
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE ddzdn_24

!*******************************************************************************
SUBROUTINE ddxup_24 (ds, a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(24,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,24
  do iy=1,24
  do ix=1,24-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxup_24

!*******************************************************************************
SUBROUTINE ddyup_24 (ds, a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,24,:) = 0.0
  c = 1./ds(2)
  do iz=1,24
  do iy=1,24-1
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE ddyup_24

!*******************************************************************************
SUBROUTINE ddzup_24 (ds, a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,:,24) = 0.0
  c = 1./ds(3)
  do iz=1,24-1
  do iy=1,24
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE ddzup_24

!*******************************************************************************
SUBROUTINE xdn_24 (a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(1,:,:) = a(1,:,:)
  do iz=1,24
  do iy=1,24
  do ix=2,24
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE xdn_24

!*******************************************************************************
SUBROUTINE ydn_24 (a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,1,:) = a(:,1,:)
  do iz=1,24
  do iy=2,24
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ydn_24

!*******************************************************************************
SUBROUTINE zdn_24 (a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,1) = a(:,:,1)
  do iz=2,24
  do iy=1,24
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE zdn_24

!*******************************************************************************
SUBROUTINE xup_24 (a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(24,:,:) = a(24,:,:)
  do iz=1,24
  do iy=1,24
  do ix=1,24-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE yup_24 (a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,24,:) = a(:,24,:)
  do iz=1,24
  do iy=1,24-1
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE zup_24 (a, b)
  real, dimension(24,24,24), intent(in):: a
  real, dimension(24,24,24):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,24) = a(:,:,24)
  do iz=1,24-1
  do iy=1,24
  do ix=1,24
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE

END MODULE
