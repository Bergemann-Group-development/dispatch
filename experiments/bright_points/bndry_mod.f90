! $Id$ 
!===============================================================================
! Variable centering in the PP-code is such that for face centered variables
! the physical boundaries are at index lb-1/2 and ub-1/2.  Thus the ranges for
! which ghost zones should be set are
!                 lower        upper
! face centered:  1,lb-1       ub  ,n
! cell centered:  1,lb         ub+1,n
!===============================================================================
MODULE bndry_mod
  USE mesh_mod,              ONLY: mesh_t
  implicit none
  type, public :: bndry_t
  contains
    procedure,nopass:: symmetric_x_lb_center
    procedure,nopass:: symmetric_x_ub_center
    procedure,nopass:: symmetric_x_lb_face
    procedure,nopass:: symmetric_x_ub_face
    procedure,nopass:: extrapolate_x_lb_center
    procedure,nopass:: extrapolate_x_ub_center
    procedure,nopass:: extrapolate_x_lb_face
    procedure,nopass:: extrapolate_x_ub_face
    procedure,nopass:: log_extrapolate_x_lb_center
    procedure,nopass:: log_extrapolate_x_ub_center
    procedure,nopass:: antisymmetric_x_lb_center
    procedure,nopass:: antisymmetric_x_ub_center
    procedure,nopass:: antisymmetric_x_lb_face
    procedure,nopass:: antisymmetric_x_ub_face
    procedure,nopass:: divB_correction_x_lb
    procedure,nopass:: divB_correction_x_ub
  end type

CONTAINS
!> =============================================================================
SUBROUTINE symmetric_x_lb_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,lb
  lb = m(1)%lo
  do ix=1,lb-1
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z-cell center (and determined on the boundary by physics)
    v(lb-ix,iy,iz) = v(lb+ix,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE symmetric_x_lb_center
!===============================================================================
SUBROUTINE extrapolate_x_lb_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,lb
  lb = m(1)%lo
  do ix=1,lb-1
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be known at the boundary, with symmetric derivative
    v(lb-ix,iy,iz) = 2.*v(lb,iy,iz)-v(lb+ix,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE extrapolate_x_lb_center
!===============================================================================
SUBROUTINE log_extrapolate_x_lb_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,lb
  lb = m(1)%lo
  do ix=1,lb-1
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be known at the boundary, with symmetric derivative
    v(lb-ix,iy,iz) = exp(2.*log(v(lb,iy,iz))-log(v(lb+ix,iy,iz)))
  enddo
  enddo
  enddo
END SUBROUTINE log_extrapolate_x_lb_center
!===============================================================================
SUBROUTINE extrapolate_x_ub_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,ub
  ub = m(1)%uo
  do ix=ub+1,m(1)%gn
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be known at the boundary, with symmetric derivative
    v(ix,iy,iz) = 2.*v(ub,iy,iz)-v(ub-(ix-ub),iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE extrapolate_x_ub_center
!===============================================================================
SUBROUTINE log_extrapolate_x_ub_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,ub
  ub = m(1)%uo
  do ix=ub+1,m(1)%gn
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be known at the boundary, with symmetric derivative
    v(ix,iy,iz) = exp(2.*log(v(ub,iy,iz))-log(v(ub-(ix-ub),iy,iz)))
  enddo
  enddo
  enddo
END SUBROUTINE log_extrapolate_x_ub_center
!===============================================================================
SUBROUTINE antisymmetric_x_lb_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,lb
  lb = m(1)%lo
  do ix=1,lb-1
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z cell center (and zero on the boundary)
    v(lb-ix,iy,iz) = - v(lb+ix,iy,iz)
  enddo
  enddo
  enddo
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z cell center (and zero on the boundary)
    v(lb,iy,iz) = 0.
  enddo
  enddo
END SUBROUTINE antisymmetric_x_lb_center
!===============================================================================
SUBROUTINE symmetric_x_lb_face(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,lb
  lb = m(1)%lo
  do ix=1,lb-1
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z face center, with physical boundary at lb+1/2
    v(lb-ix,iy,iz) = v(lb+ix-1,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE symmetric_x_lb_face
!===============================================================================
SUBROUTINE antisymmetric_x_lb_face(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,lb
  lb = m(1)%lo
  do ix=1,lb-1
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z face center (and zero on the boundary)
    v(lb-ix,iy,iz) = - v(lb+ix-1,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE antisymmetric_x_lb_face
!===============================================================================
SUBROUTINE extrapolate_x_lb_face(v, m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,lb
  lb = m(1)%lo
  do ix=1,lb-1
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z face center (and zero on the boundary)
    v(lb-ix,iy,iz) = 2.0*v(lb,iy,iz) - v(lb+ix-1,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE extrapolate_x_lb_face
!===============================================================================
SUBROUTINE symmetric_x_ub_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,ub
  ub = m(1)%uo
  do ix=1,m(1)%gn-ub
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z cell center (and determined on the boundary by physics)
    v(ub+ix,iy,iz) = v(ub-ix,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE symmetric_x_ub_center
!===============================================================================
SUBROUTINE antisymmetric_x_ub_center(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,ub
  ub = m(1)%uo
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    v(ub,iy,iz) = 0.
  enddo
  enddo
  do ix=1,m(1)%gn-ub
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be centered at z-faces (and zero on the boundary)
    v(ub+ix,iy,iz) = - v(ub-ix,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE antisymmetric_x_ub_center
!===============================================================================
SUBROUTINE symmetric_x_ub_face(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,ub
  ub = m(1)%uo
  do ix=0,m(1)%gn-ub
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z faces
    v(ub+ix,iy,iz) = v(ub-ix-1,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE symmetric_x_ub_face
!===============================================================================
SUBROUTINE antisymmetric_x_ub_face(v,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,ub
  ub = m(1)%uo
  do ix=0,m(1)%gn-ub
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z faces
    v(ub+ix,iy,iz) = - v(ub-ix-1,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE antisymmetric_x_ub_face
!===============================================================================
SUBROUTINE extrapolate_x_ub_face(v, m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real,    dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: v
  integer                                  :: ix,iy,iz,ub
  ub = m(1)%uo
  do ix=0,m(1)%gn-ub
  do iy=1,m(2)%gn
  do iz=1,m(3)%gn
    ! variable is supposed to be at z faces
    v(ub+ix,iy,iz) = 2.*v(ub,iy,iz) - v(ub-ix-1,iy,iz)
  enddo
  enddo
  enddo
END SUBROUTINE extrapolate_x_ub_face
!===============================================================================
SUBROUTINE divB_correction_x_lb (bx,by,bz,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: bx,by,bz
  integer                                  :: iy,iz,lb
  lb = m(1)%lo
  do iy=1,m(2)%gn-1
  do iz=1,m(3)%gn-1
    bx(lb,iy,iz) = bx(lb+1,iy,iz) + (m(1)%d/m(2)%d)*(by(lb,iy+1,iz)-by(lb,iy,iz)) &
                                  + (m(1)%d/m(3)%d)*(bz(lb,iy,iz+1)-bz(lb,iy,iz))
  enddo
  enddo
END SUBROUTINE divB_correction_x_lb
!===============================================================================
SUBROUTINE divB_correction_x_ub (bx,by,bz,m)
  implicit none
  class(mesh_t), pointer                   :: m(:)
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: bx,by,bz
  integer                                  :: iy,iz,ub
  ub = m(1)%uo
  do iy=1,m(2)%gn-1
  do iz=1,m(3)%gn-1
    bx(ub+1,iy,iz) = bx(ub,iy,iz) - (m(1)%d/m(2)%d)*(by(ub,iy+1,iz)-by(ub,iy,iz)) &
                                  - (m(1)%d/m(3)%d)*(bz(ub,iy,iz+1)-bz(ub,iy,iz))
  enddo
  enddo
END SUBROUTINE divB_correction_x_ub
!===============================================================================
END MODULE bndry_mod
