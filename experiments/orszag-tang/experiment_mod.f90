!===============================================================================
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE initial_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    type(initial_t):: initial
  contains
    procedure init
    procedure update
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  !----------------------------------------------------------------------------
  ! Allocate memory and set default values
  !----------------------------------------------------------------------------
  self%mhd = .true.
  call self%solver_t%init
  call self%initial%init (self%kind)
  call self%initial%condition (self%mesh, self%mem(:,:,:,:,1,1), self%idx)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Solver update; the default should be built into the solver
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  !
  call self%solver_t%update
  !
  call trace_end
END SUBROUTINE update

END MODULE
