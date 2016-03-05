!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE io_mod
  USE mesh_mod
  implicit none
  private
  type, public:: boundary_t
  contains
    procedure condition
  end type
  type(boundary_t), public:: boundary
CONTAINS

!> =============================================================================
!> Impose boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE condition (self, mem, m, h, nv)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  type(mesh_t)           :: m(3)
  real(8)                :: h(3)
  integer                :: nv
  !-----------------------------------------------------------------------------
END SUBROUTINE condition

END MODULE boundary_mod
