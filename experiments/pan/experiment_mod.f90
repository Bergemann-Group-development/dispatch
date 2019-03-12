!===============================================================================
!> Experiment initialization and update
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  !USE particle_patch_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type
  !type (particle_patch_t):: particle_patch
CONTAINS

!===============================================================================
!> Setup a uniform initial state
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, pointer:: d(:,:,:), p(:,:,:,:)
  integer:: i
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call self%solver_t%init
  !----------------------------------------------------------------------------
  ! Support (disabled for now) for particle patches (carrying sink particles)
  !----------------------------------------------------------------------------
  !call particle_patch%init
  !call particle_patch%test
  call trace_end
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%solver_t%update
  call trace_end
END SUBROUTINE update

END MODULE
