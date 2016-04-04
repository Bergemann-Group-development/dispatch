!> =============================================================================
!> Initializiation and update of each RAMSES patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE mhd_mod
  USE scaling_mod
  USE timestep_mod
  implicit none
  private
  type, public, extends(mhd_t):: experiment_t
  contains
    procedure init
    procedure update
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call self%mhd_t%init
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  ! time evolution for this experiment has been intentionally disabled.
  call trace_end
END SUBROUTINE update

END MODULE experiment_mod
