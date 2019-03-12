!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_36
  implicit none
  public
CONTAINS

!*******************************************************************************
SUBROUTINE ddxdn_36 (ds, a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  c = 1./ds(1)
  do iz=1,36
  do iy=1,36
    do ix=2,36
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
    end do
    b(1,iy,iz) = 0.0
  end do
  end do
END SUBROUTINE ddxdn_36

!*******************************************************************************
SUBROUTINE ddydn_36 (ds, a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,1,:) = 0.0
  c = 1./ds(2)
  do iz=1,36
  do iy=2,36
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ddydn_36

!*******************************************************************************
SUBROUTINE ddzdn_36 (ds, a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: ds(3), c
  integer :: ix, iy, iz
!...............................................................................
  b(:,:,1) = 0.0
  c = 1./ds(3)
  do iz=2,36
  do iy=1,36
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE ddzdn_36

!*******************************************************************************
SUBROUTINE ddxup_36 (ds, a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  c = 1./ds(1)
  do iz=1,36
  do iy=1,36
    do ix=1,36-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
    end do
    b(36,iy,iz) = 0.0
  end do
  end do
END SUBROUTINE ddxup_36

!*******************************************************************************
SUBROUTINE ddyup_36 (ds, a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,36,:) = 0.0
  c = 1./ds(2)
  do iz=1,36
  do iy=1,36-1
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE ddyup_36

!*******************************************************************************
SUBROUTINE ddzup_36 (ds, a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  integer :: ix, iy, iz
  real    :: ds(3), c
!...............................................................................
  b(:,:,36) = 0.0
  c = 1./ds(3)
  do iz=1,36-1
  do iy=1,36
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE ddzup_36

!*******************************************************************************
SUBROUTINE xdn_36 (a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(1,:,:) = a(1,:,:)
  do iz=1,36
  do iy=1,36
  do ix=2,36
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
END SUBROUTINE xdn_36

!*******************************************************************************
SUBROUTINE ydn_36 (a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,1,:) = a(:,1,:)
  do iz=1,36
  do iy=2,36
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
END SUBROUTINE ydn_36

!*******************************************************************************
SUBROUTINE zdn_36 (a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,1) = a(:,:,1)
  do iz=2,36
  do iy=1,36
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
  end do
  end do
  end do
END SUBROUTINE zdn_36

!*******************************************************************************
SUBROUTINE xup_36 (a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(36,:,:) = a(36,:,:)
  do iz=1,36
  do iy=1,36
  do ix=1,36-1
    b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE yup_36 (a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,36,:) = a(36,:,:)
  do iz=1,36
  do iy=1,36-1
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
  end do
  end do
  end do
END SUBROUTINE

!*******************************************************************************
SUBROUTINE zup_36 (a, b)
  real, dimension(36,36,36), intent(in):: a
  real, dimension(36,36,36):: b
  real    :: c
  integer :: ix, iy, iz
!...............................................................................
  c = 0.5
  b(:,:,36) = a(:,:,36)
  do iz=1,36-1
  do iy=1,36
  do ix=1,36
    b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy,iz  ))
  end do
  end do
  end do
END SUBROUTINE

END MODULE
