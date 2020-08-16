!===============================================================================
!> Boundary conditions for centered variables, which have the physical boundary
!> in-between mesh points lo and li, or/and in-between mesh points ui and uo.
!===============================================================================
MODULE boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE index_mod
  USE standard_boundaries_mod
  USE patch_mod
  implicit none
  private
  type, extends(standard_boundaries_t), public:: boundaries_t
    ! --------------------------------------------------------------------------
    ! Experiment - specific variables
    ! --------------------------------------------------------------------------
  contains
    ! --------------------------------------------------------------------------
    ! Default procedures
    ! --------------------------------------------------------------------------
    procedure:: init
    procedure:: primary
    procedure:: derivative
    ! --------------------------------------------------------------------------
    ! Experiment - specific procedures
    ! --------------------------------------------------------------------------
  end type
CONTAINS

!===============================================================================
!> Initialize the boundary data type, while testing also if the patch contains
!> a spherical boundary, or possibly is completely immersed in such a boundary.
!===============================================================================
SUBROUTINE init (self, patch)
  class(boundaries_t)     :: self
  class(patch_t)          :: patch
  !.............................................................................
  class(mesh_t) , pointer :: mesh(:)
  real(8)                 :: position(3), radius
  integer                 :: id
  !.............................................................................
  id = patch%id
  mesh => patch%mesh
  position = patch%position
  radius = patch%radius
  call self%standard_boundaries_t%init (patch)
END SUBROUTINE init

SUBROUTINE primary(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer :: mem
  ! ----------------------------------------------------------------------------
  !if (self%mesh(3)%lower_boundary) call self%lower(mem)
  !if (self%mesh(3)%upper_boundary) call self%upper(mem)

END SUBROUTINE primary

SUBROUTINE derivative(self, mem)
  class(boundaries_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:), pointer :: mem
  !
  !if (self%mesh(3)%lower_boundary) call self%lower_damp(mem)

END SUBROUTINE derivative


END MODULE boundaries_mod
