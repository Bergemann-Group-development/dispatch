!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_32
  implicit none
  public
CONTAINS

!*******************************************************************************
SUBROUTINE ddxdn_32 (ds, a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(1,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,32
  do iy=1,32
  do ix=2,32
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxdn_32

!*******************************************************************************
SUBROUTINE ddydn_32 (ds, a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,1,:) = 0.0
  c = 1./ds(2)
  do iz=1,32
  do iy=2,32
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ddydn_32

!*******************************************************************************
SUBROUTINE ddzdn_32 (ds, a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(:,:,1) = 0.0
  c = 1./ds(3)
  do iz=2,32
  do iy=1,32
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE ddzdn_32

!*******************************************************************************
SUBROUTINE ddxup_32 (ds, a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(32,:,:) = 0.0
  c = 1./ds(1)
  do iz=1,32
  do iy=1,32
  do ix=1,32-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE ddxup_32

!*******************************************************************************
SUBROUTINE ddyup_32 (ds, a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,32,:) = 0.0
  c = 1./ds(2)
  do iz=1,32
  do iy=1,32-1
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE ddyup_32

!*******************************************************************************
SUBROUTINE ddzup_32 (ds, a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,:,32) = 0.0
  c = 1./ds(3)
  do iz=1,32-1
  do iy=1,32
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE ddzup_32

!*******************************************************************************
SUBROUTINE xdn_32 (a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(1,:,:) = a(1,:,:)
  do iz=1,32
  do iy=1,32
  do ix=2,32
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE xdn_32

!*******************************************************************************
SUBROUTINE ydn_32 (a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,1,:) = a(:,1,:)
  do iz=1,32
  do iy=2,32
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ydn_32

!*******************************************************************************
SUBROUTINE zdn_32 (a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,1) = a(:,:,1)
  do iz=2,32
  do iy=1,32
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE zdn_32

!*******************************************************************************
SUBROUTINE xup_32 (a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(32,:,:) = a(32,:,:)
  do iz=1,32
  do iy=1,32
  do ix=1,32-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE yup_32 (a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,32,:) = a(:,32,:)
  do iz=1,32
  do iy=1,32-1
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE zup_32 (a, b)
  real, dimension(32,32,32), intent(in):: a
  real, dimension(32,32,32):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,32) = a(:,:,32)
  do iz=1,32-1
  do iy=1,32
  do ix=1,32
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE

END MODULE
