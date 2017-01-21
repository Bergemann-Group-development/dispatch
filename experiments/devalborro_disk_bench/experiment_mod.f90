!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE mhd_mod
  USE boundary_mod
  USE initial_mod
  USE potential_mod
  implicit none
  private

  type, public, extends(mhd_t):: experiment_t
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
  call trace_begin('experiment_t%init')

  call self%mhd_t%init()
  call initial_conditions(self)
  call potential%init()
  call boundary%conditions(self)

  call trace_end
END SUBROUTINE init_experiment

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call boundary%conditions(self)
  ! update the potential (for a growing or moving planet)
  potential%omega = self%grid_vel(3)
  self%gp%q = potential%selected(self%mesh, self%time, self%dtime)
  call self%mhd_t%update()
  call boundary%conditions(self)

  call trace_end
END SUBROUTINE update_experiment

END MODULE experiment_mod
