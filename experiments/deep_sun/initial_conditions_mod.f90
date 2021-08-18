!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE initial_conditions_mod
  !USE const
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE omp_mod
  USE index_mod
  USE patch_mod
  USE gpatch_mod
  USE io_unit_mod
  USE math_mod
  USE kinds_mod
  USE eos_mod
  implicit none
  private
  type, public:: initial_conditions_t

  contains
    procedure:: pre_init
    procedure:: init
  end type
  real, save:: csound=1.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace_begin('initial_t%init')

  call trace%end()
END subroutine pre_init

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t)            :: self

  call trace%end()
END SUBROUTINE init

END MODULE initial_conditions_mod
