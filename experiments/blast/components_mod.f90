!===============================================================================
!> Initializiation for the MHD blast experiment
!> The set up draws heavily from the "blast" problem in AZEuS/ZEUS-3D.
!===============================================================================
MODULE components_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE mhd_mod
  USE rubiks_mod
  implicit none
  private

  type, public, extends(rubiks_t):: components_t
  contains
    procedure:: init => init_components
  end type

CONTAINS
!===============================================================================
SUBROUTINE init_components (self)
  class(components_t):: self
  !-----------------------------------------------------------------------------
  call trace_begin ('components_t%init')

  call self%rubiks_t%init_hierarchy

  call trace_end
END SUBROUTINE init_components

END MODULE components_mod
