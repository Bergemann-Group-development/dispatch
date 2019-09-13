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
    integer:: key(4)
  contains
    procedure init
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
  print *,'experiment_t%init: unisigned =', &
    self%unsigned(self%idx%d), self%unsigned(self%idx%e)
  self%unsigned(self%idx%d) = .true.
  self%unsigned(self%idx%e) = .true.
  call trace%end (itimer)
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  call self%output
  if (self%level >= shared%levelmin) then
    call self%solver_t%update
  else
    call self%set (bits%frozen)
  end if
  call trace_end
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE info (self, nq, ntask, experiment_name)
  class(experiment_t):: self
  integer, optional:: nq, ntask
  character(len=64), optional:: experiment_name
  !----------------------------------------------------------------------------
  if (self%level >= shared%levelmin) &
    call self%solver_t%info (nq, ntask, experiment_name)
END SUBROUTINE info

END MODULE experiment_mod
