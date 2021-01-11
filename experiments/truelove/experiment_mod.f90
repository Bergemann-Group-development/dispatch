!===============================================================================
!> Initializiation and update for the Truelove grav. collapse experiment
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE bits_mod
  USE solver_mod
  USE pboundary_mod
  implicit none
  private

  type, public, extends(solver_t):: experiment_t
  contains
    procedure :: init
    procedure :: update
  end type

CONTAINS

!===============================================================================
!> Initialise this experiment.  The experiment can be run both as periodic and
!> non-periodic, depending on the setting 'fixedbcs' in pboundary_mod.f90
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t) :: self
  !----------------------------------------------------------------------------
  if (self%mem_allocated) return
  call trace%begin('experiment_t%init')
  call self%solver_t%init()
  call pboundaries%init(self)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Update procedure for this experiment. In this case the call falls through to
!> mhd_t%update, which organizes the preliminary and final selfgravity solutions
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer:: m
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  call pboundaries%conditions(self, self%it)
  call self%solver_t%update()
  call pboundaries%conditions(self, self%new)
  call trace%end()
END SUBROUTINE update

END MODULE experiment_mod
