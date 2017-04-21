!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE boundary_mod
  USE mhd_mod
  USE initial_mod
  USE pboundary_mod
  implicit none
  private

  type, public, extends(mhd_t):: experiment_t
    type(boundary_t):: boundary
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

  self%nv = 8
  call self%mhd_t%init()
  call self%initial%init(self)
  call pboundaries%conditions(self, self%it)

  call self%output

  call trace_end
END SUBROUTINE init_experiment

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')

  call pboundaries%conditions(self, self%it)
  call self%mhd_t%update()
  call pboundaries%conditions(self, self%new)

  call trace_end
END SUBROUTINE update_experiment

END MODULE experiment_mod
