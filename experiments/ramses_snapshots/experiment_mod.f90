!===============================================================================
!> Initializiation and update of each RAMSES patch
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE bits_mod
  USE shared_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure init
    procedure dnload
    procedure update
    procedure info
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
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%dnload')
  call self%solver_t%dnload
  call trace%end()
END SUBROUTINE dnload

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  call self%output
  if (self%is_clear (bits%frozen)) then
    call self%solver_t%update
  end if
  call trace%end()
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(experiment_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%info')
  if (self%is_clear (bits%frozen)) then
    call self%solver_t%info (nq, ntask, experiment_name)
  end if
  call trace%end()
END SUBROUTINE info

END MODULE experiment_mod
