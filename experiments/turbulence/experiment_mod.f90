!===============================================================================
!> This module contains all experiment specific information.  The experiment_t
!> object extends the generic MHD object mhd_t, which is again an extension of
!> the basic patch_t object for tasks (object task_t) with meshes.
!>
!> The MHD object is chosen with the SOLVER macro in the Makefile in the curren
!> directory, while the generic task_t and patch_t objects are defined in the
!> $(TOP)/tasks/ directory, in task_mod.f90 and patch_mod.f90.
!>
!> The init procedure here (experiment_t::init) is called from the cartesian_mod
!> module, which is specified as the patch distributor in the main program
!> dispatch.f90.   The update procedure here is called by the task_list_mod
!> update procedure in $(TOP)/lists/task_list_mod.f90.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE initial_mod
  USE timer_mod
  USE link_mod
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
!> Setup an initial state.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin('experiment_t%init')
  call self%solver_t%init
  self%d0 = self%initial%d0
  call self%set (bits%static)
  if (self%initial%time > 0d0) then
    self%time = self%initial%time
    self%t(self%it) = self%initial%time
  end if
  !call self%periodic_grid
  !call particle_patch%init
  !call particle_patch%test
  call trace%end
END SUBROUTINE init

!===============================================================================
!> The update procedure calls a (generic name) MHD update procedure, whose
!> actual implementation is in $(TOP)/solvers/$(SOLVER)/mhd_mod.f90.  It then
!> calls a courant_condition procedure, which if not overloaded is the one in
!> the patch_t object (but one can replace it at the mhd_t or experiment_t
!> level, by simply adding a courant_condition procedure there).
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  class(link_t), pointer:: link
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%update')
  if (first_time) then
    !$omp critical (cleanup_cr)
    if (first_time) then
      first_time = .false.
      call self%initial%cleanup
    end if
    !$omp end critical (cleanup_cr)
  end if
  call self%solver_t%update
  call trace%end
END SUBROUTINE update

END MODULE
