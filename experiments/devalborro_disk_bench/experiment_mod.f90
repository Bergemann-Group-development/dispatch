!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE boundary_mod
  USE mhd_mod
  USE pboundary_mod
  USE initial_mod
  USE potential_mod
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

  call self%mhd_t%init()
  call initial_conditions(self)
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
  ! update the potential (for a growing or moving planet)
  potential%omega = self%grid_vel(3)
  self%mem(:,:,:,self%i_phi,self%it,1) = potential%selected(self%mesh, self%time, self%dtime)
  call self%mhd_t%update()
  call pboundaries%conditions(self, self%new)

  call trace_end
END SUBROUTINE update_experiment

END MODULE experiment_mod
