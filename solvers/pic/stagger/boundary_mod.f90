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
SUBROUTINE condition (self, mem, m, h)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  class(mesh_t), pointer :: m(:)
  real(8)                :: h(3)
  !-----------------------------------------------------------------------------
END SUBROUTINE condition

END MODULE boundary_mod
