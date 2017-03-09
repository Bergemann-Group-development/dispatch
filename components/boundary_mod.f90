!> =============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!> =============================================================================
MODULE boundary_mod
  USE trace_mod
  USE mpi_mod
  USE bits_mod
  implicit none
  private
  type, public:: boundary_t
    integer:: status = 0
    logical:: unset = .true.
    real(8):: cool_temp, hot_temp
  contains
    procedure:: init
    procedure:: set
    procedure:: apply
  end type
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

!> =====================================================================
!> Set the boundary condition, whatever it means
!> =====================================================================
SUBROUTINE apply(self)
  class(boundary_t):: self
  
END SUBROUTINE apply

END MODULE boundary_mod
