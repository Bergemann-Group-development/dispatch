!===============================================================================
!> Minimal experiment initialization -- all other calls fall through to solver
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
  end type
CONTAINS

!===============================================================================
!> Generic init procedure
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  real, pointer:: d(:,:,:), p(:,:,:,:)
  integer:: i
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call self%solver_t%init
  !----------------------------------------------------------------------------
  ! Set type of variables, for interpolations
  !----------------------------------------------------------------------------
  self%unsigned(self%idx%d) = .true.
  self%unsigned(self%idx%e) = .true.
  call trace_end
END SUBROUTINE init

END MODULE
