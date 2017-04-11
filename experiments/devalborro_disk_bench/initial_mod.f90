!> =============================================================================
!> A dummy `initial_mod` to override the default one.
!> =============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  implicit none
  private

  type, public:: initial_t
  contains
    procedure:: init
    procedure:: selected
  end type

CONTAINS

!> =============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!> =============================================================================
SUBROUTINE init (self)
  class(initial_t):: self
END SUBROUTINE init

SUBROUTINE selected (self, m, ff)
  class(initial_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: ff
  class(mesh_t), pointer, dimension(:):: m
END SUBROUTINE selected

END MODULE initial_mod
