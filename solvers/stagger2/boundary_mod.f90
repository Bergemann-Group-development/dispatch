!> =============================================================================
!> Boundary conditions for a planet surface
!> =============================================================================
MODULE boundary_mod
  USE trace_mod
  USE io_mod
  USE mesh_mod
  implicit none
  private
  type, public:: boundary_t
    integer:: status = 0
    logical:: unset = .true.
    class(mesh_t), pointer:: m(:)
  contains
    procedure:: init
    procedure:: set
    procedure condition
  end type
  type(boundary_t), public:: boundary
CONTAINS

!> =============================================================================
!> Initialize the boundary data type
!> =============================================================================
SUBROUTINE init (self, mesh, position, radius)
  class(boundary_t):: self
  class(mesh_t), pointer, optional:: mesh(:)
  real(8), optional:: position(3), radius
  !-----------------------------------------------------------------------------
  call trace%begin ('boundary_t%init')

  self%unset = .false.
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Set bits in the status word
!> =============================================================================
SUBROUTINE set (self, bits)
  class(boundary_t):: self
  integer:: bits
  !-----------------------------------------------------------------------------
  self%status = ior (self%status, bits)
END SUBROUTINE set

!> =============================================================================
!> Check bits in the status word
!> =============================================================================
FUNCTION is_set (self, bits) RESULT (out)
  class(boundary_t):: self
  integer:: bits
  logical:: out
  !-----------------------------------------------------------------------------
  out = iand (self%status, bits)
END FUNCTION is_set


!> =============================================================================
!> Impose boundary conditions, taking staggering into account
!> =============================================================================
SUBROUTINE condition (self, mem, m, h)
  class(boundary_t)      :: self
  real, dimension(:,:,:) :: mem
  class(mesh_t)          :: m(3)
  real(8)                :: h(3)
  !-----------------------------------------------------------------------------
END SUBROUTINE condition

END MODULE boundary_mod
