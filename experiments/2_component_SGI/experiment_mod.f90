!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE solver_mod
  USE pboundary_mod
  USE potential_mod
  USE parameters_mod
  USE scaling_mod
  USE index_mod
  USE download_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
  end type

  real(8), parameter :: pi = asin(1.0d0) * 2.0d0
CONTAINS

!> =============================================================================
!> Initialise this experiment.
!> =============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')
  self%nv = experiment_params%nv                        ! Number of total variables 
  self%nc = experiment_params%nc                        ! Number of fluid components
  call self%mhd_t%init()
  if (self%nc==2) self%no_pressure = 2
  call potential%init()
  call pboundaries%conditions(self, self%it)

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
