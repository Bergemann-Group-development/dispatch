!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE timestep_mod
  USE mesh_mod
  USE mhd_mod
  USE ic_mod
  USE pboundary_mod
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
  call ics%init()
  call ics%init_mhd_vars(self)
  call pboundary%conditions(ics, self%mhd_t, self%it)

  call trace_end
END SUBROUTINE init_experiment

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  type(timestep_t):: timestep
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call pboundary%conditions(ics, self%mhd_t, self%it)
  call self%mhd_t%update()
  call pboundary%conditions(ics, self%mhd_t, self%new)

  call trace_end
END SUBROUTINE update_experiment

END MODULE experiment_mod
