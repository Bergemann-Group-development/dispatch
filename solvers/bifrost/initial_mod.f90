!===============================================================================
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  type, public:: initial_t
    integer:: id=0
    character(len=64):: solver
    procedure(void), pointer:: condition=>null()
  contains
    procedure:: init
    procedure, nopass:: cleanup
  end type
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma)
   class(initial_t):: self
   character(len=64):: solver
   real, optional:: gamma
   !............................................................................
END SUBROUTINE init

!===============================================================================
SUBROUTINE cleanup
END SUBROUTINE cleanup

!===============================================================================
!> Initialize a single 3-D wavenumber, with random amplitude and phase
!===============================================================================
SUBROUTINE void (self, m, ff, idx)
  class(initial_t)      :: self
  class(mesh_t), pointer, dimension(:):: m
  real, dimension(:,:,:,:) :: ff
  type(index_t):: idx
  !.............................................................................
END SUBROUTINE void

END MODULE
