!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE initial_mod
  USE pboundary_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type
  logical, save:: first_time = .true.

CONTAINS

!> =============================================================================
!> Initialise the experiment data type.
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')

  self%nv = 2
  self%periodic = .false.
  
  call self%gpatch_t%init()

  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  
  call pboundaries%conditions(self, self%it)
  call self%gravity()
  call self%dnload  (only=self%i_phi)
  call self%gravity (only=self%i_phi)
  call pboundaries%conditions(self, self%new)

  call trace%end
END SUBROUTINE update

END MODULE experiment_mod
