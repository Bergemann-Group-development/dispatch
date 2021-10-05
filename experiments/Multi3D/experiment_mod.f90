!===============================================================================
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE multi3d_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  call trace_end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call multi3d%update
  call trace_end
END SUBROUTINE update

END MODULE
