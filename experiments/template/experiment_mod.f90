!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!=======================================================================
!> Template experiment
!=======================================================================
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE trace_mod
  USE solver_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!===============================================================================
!> Initialize
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin ('experiment_t%init')
  call self%solver_t%init ()
  self%mem = 1.0
  call trace_end ()
END SUBROUTINE init

!===============================================================================
!> Update
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: iv,i1,i2
  logical :: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin ('experiment_t%update')
  call self%solver_t%update ()
  call trace_end ()
END SUBROUTINE update

END MODULE experiment_mod
