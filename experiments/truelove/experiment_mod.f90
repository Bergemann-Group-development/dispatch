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
  USE validate_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: update
    procedure:: bndries
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
  call validate%init
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Update procedure for this experiment. In this case the call falls through to
!> mhd_t%update, which organizes the preliminary and final selfgravity solutions
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer:: m, it
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  ! call self%stats ('before upd', stdout)
  call self%solver_t%update()
  it = self%it
  self%it = self%new
  ! call self%stats ('after  upd', stdout)
  self%it = it
  call trace%end()
END SUBROUTINE update

!===============================================================================
!> Hook that catches the halo_t%update call to task%bndries
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  integer:: m
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%bndries')
  call pboundaries%conditions(self, self%it)
  call trace%end()
END SUBROUTINE bndries

END MODULE experiment_mod
