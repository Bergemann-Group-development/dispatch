!> =============================================================================
!> Initializiation and update of each RAMSES patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE mhd_mod
  USE scaling_mod
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
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  ! the update procedure for this experiment has been intentionally disabled;
  ! it is not designed to be evolved on its own, but rather used as a set up for
  ! another experiment.
  self%dtime = 0.001 ! dummy time step
  call trace_end
END SUBROUTINE update

END MODULE experiment_mod
