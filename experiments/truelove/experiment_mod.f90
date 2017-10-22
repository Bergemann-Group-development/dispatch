!> =============================================================================
!> Initializiation and update for the Truelove grav. collapse experiment
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE solver_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
  end type

CONTAINS

!> =============================================================================
!> Initialise this experiment.
!> =============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace%begin('experiment_t%init')

  self%nv = 5
  self%periodic = .true. ! not periodic in the y-direction
  self%mhd = .false.

  call self%solver_t%init()

  call trace%end
END SUBROUTINE init_experiment

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')

  call self%solver_t%update()

  call trace%end
END SUBROUTINE update_experiment

END MODULE experiment_mod
