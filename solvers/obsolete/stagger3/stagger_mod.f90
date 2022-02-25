!===============================================================================
!> $Id$
!===============================================================================
MODULE stagger_mod
  USE io_mod
  USE mesh_mod
  implicit none
  private
  type, public:: stagger_t
    class(mesh_t), dimension(:), pointer:: mesh
  contains
    procedure ddxdn1
    procedure ddydn1
    procedure ddzdn1
    procedure ddxup1
    procedure ddyup1
    procedure ddzup1
    procedure xdn1
    procedure ydn1
    procedure zdn1
    procedure xup1
    procedure yup1
    procedure zup1
    procedure ddxdn
    procedure ddydn
    procedure ddzdn
    procedure ddxup
    procedure ddyup
    procedure ddzup
    procedure xdn
    procedure ydn
    procedure zdn
    procedure xup
    procedure yup
    procedure zup
  end type
CONTAINS

!*******************************************************************************
FUNCTION ddxdn1 (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddxdn (self, mesh, a)
END FUNCTION ddxdn1

!*******************************************************************************
FUNCTION ddydn1 (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddydn (self, mesh, a)
END FUNCTION ddydn1

!*******************************************************************************
FUNCTION ddzdn1 (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddzdn (self, mesh, a)
END FUNCTION ddzdn1

!*******************************************************************************
FUNCTION ddxup1 (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddxup (self, mesh, a)
END FUNCTION ddxup1

!*******************************************************************************
FUNCTION ddyup1 (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddyup (self, mesh, a)
END FUNCTION ddyup1

!*******************************************************************************
FUNCTION ddzup1 (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ddzup (self, mesh, a)
END FUNCTION ddzup1

!*******************************************************************************
FUNCTION xdn1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = xdn (self, a)
END FUNCTION xdn1

!*******************************************************************************
FUNCTION ydn1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = ydn (self, a)
END FUNCTION ydn1

!*******************************************************************************
FUNCTION zdn1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = zdn (self, a)
END FUNCTION zdn1

!*******************************************************************************
FUNCTION xup1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = xup (self, a)
END FUNCTION xup1

!*******************************************************************************
FUNCTION yup1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = yup (self, a)
END FUNCTION yup1

!*******************************************************************************
FUNCTION zup1 (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  b = zup (self, a)
END FUNCTION zup1

!*******************************************************************************
FUNCTION ddxdn (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = mesh%lb
  u = mesh%ub
  if (u(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./mesh(1)%d
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1)+1,u(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) - a(ix-1,iy,iz))
    end do
    end do
    end do
  endif
  b(l(1),:,:) = 0.0
END FUNCTION ddxdn

!*******************************************************************************
FUNCTION ddydn (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = mesh%lb
  u = mesh%ub
  if (u(2) > 1) then
    c = 1./mesh(2)%d
    do iz=l(3),u(3)
    do iy=l(2)+1,u(2)
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy  ,iz) - a(ix,iy-1,iz))
    end do
    end do
    end do
  endif
  b(:,l(2),:) = 0.0
END FUNCTION ddydn

!*******************************************************************************
FUNCTION ddzdn (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = mesh%lb
  u = mesh%ub
  if (u(3) > 1) then
    c = 1./mesh(3)%d
    do iz=l(3)+1,u(3)
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz  ) - a(ix,iy,iz-1))
    end do
    end do
    end do
  endif
  b(:,:,l(3)) = 0.0
END FUNCTION ddzdn

!*******************************************************************************
FUNCTION ddxup (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = mesh%lb
  u = mesh%ub
  if (io%verbose>2) print *,l,u
  if (u(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 1./mesh(1)%d
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1),u(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) - a(ix  ,iy,iz))
    end do
    end do
    end do
  endif
  b(u(1),:,:) = 0.0
END FUNCTION ddxup

!*******************************************************************************
FUNCTION ddyup (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = mesh%lb
  u = mesh%ub
  if (u(2) > 1) then
    c = 1./mesh(2)%d
    do iz=l(3),u(3)
    do iy=l(2),u(2)-1
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) - a(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  b(:,u(2),:) = 0.0
END FUNCTION ddyup

!*******************************************************************************
FUNCTION ddzup (self, mesh, a) RESULT (b)
  class(stagger_t):: self
  class(mesh_t), dimension(:), pointer:: mesh
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = mesh%lb
  u = mesh%ub
  if (u(3) > 1) then
    c = 1./mesh(3)%d
    do iz=l(3),u(3)-1
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) - a(ix,iy,iz  ))
    end do
    end do
    end do
  endif
  b(:,:,u(3)) = 0.0
END FUNCTION ddzup

!*******************************************************************************
FUNCTION xdn (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer:: ix, iy, iz, l(3), u(3)
  real:: c
!...............................................................................
  l = 1
  u = shape(a)
  if (u(1) > 1) then
    c = 0.5
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1)+1,u(1)
      b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
    end do
    end do
    end do
  end if
  b(l(1),:,:) = a(l(1),:,:)
END FUNCTION xdn

!*******************************************************************************
FUNCTION ydn (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real:: c
  integer:: ix, iy, iz, l(3), u(3)
!...............................................................................
  l = 1
  u = shape(a)
  if (u(2) > 1) then
    c = 0.5
    do iz=l(3),u(3)
    do iy=l(2)+1,u(2)
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
    end do
    end do
    end do
  end if
  b(:,l(2),:) = a(:,l(2),:)
END FUNCTION ydn

!*******************************************************************************
FUNCTION zdn (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real:: c
  integer:: ix, iy, iz, l(3), u(3)
!...............................................................................
  l = 1
  u = shape(a)
  if (u(3) > 1) then
    c = 0.5
    do iz=l(3)+1,u(3)
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy,iz-1))
    end do
    end do
    end do
  end if
  b(:,:,l(3)) = a(:,:,l(3))
END FUNCTION zdn

!*******************************************************************************
FUNCTION xup (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = 1
  u = shape(a)
  if (u(1) > 1) then ! x-derivative only non-zero iff x-dimension larger than 1
    c = 0.5
    do iz=l(3),u(3)
    do iy=l(2),u(2)
    do ix=l(1),u(1)-1
      b(ix,iy,iz) = c*(a(ix+1,iy,iz) + a(ix  ,iy,iz))
    end do
    end do
    end do
  endif
  b(u(1),:,:) = a(u(1),:,:)
END FUNCTION

!*******************************************************************************
FUNCTION yup (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = 1
  u = shape(a)
  if (u(2) > 1) then
    c = 0.5
    do iz=l(3),u(3)
    do iy=l(2),u(2)-1
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy+1,iz) + a(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  b(:,u(2),:) = a(:,u(2),:)
END FUNCTION

!*******************************************************************************
FUNCTION zup (self, a) RESULT (b)
  class(stagger_t):: self
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  integer :: ix, iy, iz, l(3), u(3)
  real    :: c
!...............................................................................
  l = 1
  u = shape(a)
  if (u(3) > 1) then
    c = 0.5
    do iz=l(3),u(3)-1
    do iy=l(2),u(2)
    do ix=l(1),u(1)
      b(ix,iy,iz) = c*(a(ix,iy,iz+1) + a(ix,iy  ,iz))
    end do
    end do
    end do
  endif
  b(:,:,u(3)) = a(:,:,u(3))
END FUNCTION

!*******************************************************************************
END MODULE
