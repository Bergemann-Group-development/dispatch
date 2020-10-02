!===============================================================================
!> Simple initical condition module for tests
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE mesh_mod
  USE index_mod
  implicit none
  private
  type, public:: initial_t
    logical:: mhd=.false.
  contains
    procedure:: init
    procedure:: condition
  end type
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma1)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma1
  !.............................................................................
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
END SUBROUTINE init

!===============================================================================
!> Initialize a single 3-D wavenumber
!===============================================================================
SUBROUTINE condition (self, m, ff, idx)
  class(initial_t) :: self
  class(mesh_t)    :: m(3)
  real             :: ff(:,:,:,:)
  type(index_t)    :: idx
  !-----------------------------------------------------------------------------
END SUBROUTINE condition

END MODULE
