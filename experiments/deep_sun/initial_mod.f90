!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_mod
  !USE const
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE omp_mod
  USE index_mod
  USE io_unit_mod
  USE math_mod
  USE kinds_mod
  implicit none
  private
  type, public:: initial_t
    integer:: seed
    logical:: mhd=.true.
    real(8):: time=0d0
    character(len=64):: solver
  contains
    procedure:: init
    procedure:: condition
  end type
  real, save:: csound=1.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real(8), optional:: gamma
  !---------------------------------------------------------------------
  call trace_begin('initial_t%init')
  self%solver = trim(solver)
  call trace%end()
END subroutine init


SUBROUTINE condition (self, m, ff, idx)
  class(initial_t) :: self
  type(index_t):: idx
  class(mesh_t)    :: m(3)
  real, dimension(:,:,:,:), pointer :: ff

END SUBROUTINE condition

END MODULE initial_mod
