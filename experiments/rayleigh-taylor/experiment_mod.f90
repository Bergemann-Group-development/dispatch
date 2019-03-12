!===============================================================================
!> Initializiation and update for the Rayleigh-Taylor experiment
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE mesh_mod
  USE pboundary_mod
  USE mhd_mod
  implicit none
  private

  type, public, extends(mhd_t):: experiment_t
  contains
    procedure :: init => init_experiment
    procedure :: update => update_experiment
  end type

  logical, save :: first_time = .true.

CONTAINS

!===============================================================================
!> Initialise this experiment.
!===============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace_begin('experiment_t%init')

  self%nv = 5
  self%periodic(2) = .false. ! not periodic in the y-direction
  self%mhd = .false.

  call self%force%init(self%kind)
  call self%mhd_t%init()

  call trace_end
END SUBROUTINE init_experiment

!===============================================================================
!> Update procedure for this experiment.
!===============================================================================
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
