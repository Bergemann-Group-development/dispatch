!> =============================================================================
!> Initializiation and update of each ZEUS patch
!> =============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE pboundaries_mod
  !!!USE pebble_patch_mod
  USE potential_mod
  implicit none
  private

  !type, public, extends(pebble_patch_t):: experiment_t
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init => init_experiment
    procedure:: update => update_experiment
  end type

CONTAINS

!> =============================================================================
!> Initialise this experiment.
!> Things that need to be initialised per-patch happens here.
!> =============================================================================
SUBROUTINE init_experiment (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')

  self%mhd = .true.
  self%initial%mhd = .true.
  call self%initial%init(self%kind, real(self%gamma))
  !call self%pebble_patch_t%init()
  call self%solver_t%init()
  ! intialize time-independent gravitational potential of the central star
  call potential%init(self%mesh, self%idx, self%initial%star_mass, self%mem(:,:,:,:,1,1))
  call pboundaries%init_bdries(self%mesh,self%position,self%idx)

  call trace%end
END SUBROUTINE init_experiment

!> =============================================================================
!> Update procedure for this experiment.
!> =============================================================================
SUBROUTINE update_experiment (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')

  call pboundaries%apply(self, self%it)
  call self%solver_t%update()
  call pboundaries%apply(self, self%new)
  !if (self%do_pebbles) call self%pebble_patch_t%update_pebbles

  call trace%end
END SUBROUTINE update_experiment

END MODULE experiment_mod
