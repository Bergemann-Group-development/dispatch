!===============================================================================
!> Initializiation and update of each RAMSES patch
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    integer:: key(4)
  contains
    procedure init
    procedure update
  end type
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init', itimer=itimer)
  call self%solver_t%init
  call trace%end (itimer)
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%output
  call self%solver_t%update
  call trace_end
END SUBROUTINE update

END MODULE experiment_mod
