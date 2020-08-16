!===============================================================================
!> Template version of a module intended for adding extra features, at a level
!> between the basic patch and the layer of solvers.
!>
!> To add extra features to the code, copy this file to experiments/your_dir/,
!> and add (only) the extra dependencies to the Makefile there.
!>
!>                    experiment_t
!>                     solver_t
!>                       mhd_t
!>                      extras_t
!>                      /    \
!>                     /      \
!>                extras_1_t extras_2_t
!>                       \    /
!>                      gpatch_t
!>                       patch_t
!>                       /
!>                    task_t
!>
!> The extra modules needed in an experiment should be chosen and coordinated
!> in this module, and the ones that are "below" (USEd by) this module should
!> generally refer to array memory in patch_t, and/or to local memory allocated
!> inside these extra modules.  Only key physics effects, such as forcing,
!> cooling or heating, ionization, etc., should be added to the patch_t module,
!>
!> Modules that prepare input to the solvers should generally be called from the
!> pre_upate() procedure, while modules that rely on results from the solver
!> should be called from the post_update() procedure.
!>
!> This module should generally NOT use the h5_mod (HDF5) module itself. It
!> should instead be USEd from the actual extra modules, to assist in saving
!> ad hoc results in a standardized format.
!===============================================================================

MODULE extras_mod
  !-----------------------------------------------------------------------------
  ! Generic modules, which are always included in the code in any case
  !-----------------------------------------------------------------------------
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE patch_mod
  USE task_mod
  USE list_mod
  USE link_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.
  !-----------------------------------------------------------------------------
  USE forces_mod                                                ! forces
  USE selfgravity_mod                                           ! selfg
  USE cooling_mod                                               ! cool
  USE particle_solver_mod                                       ! dust
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(forces_t):: forces                                     ! forces
    type(selfgravity_t):: selfgravity                           ! selfg
    type(cooling_t):: cooling                                   ! cool
    class(gpatch_t), pointer:: gpatch => null()                 ! keep
    type(particle_solver_t):: particle_solver                   ! dust
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: dealloc
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    procedure:: check_refine
    procedure:: courant_condition
    procedure:: statistics
  end type
  type(extras_t), public:: extras
CONTAINS

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_init')
  self%gpatch => self
  call self%forces%pre_init (self%gpatch)                       ! forces
  call self%cooling%pre_init (self%gpatch)                      ! cool
  call self%selfgravity%pre_init (self%gpatch)                  ! selfg
  call self%particle_solver%pre_init (self%gpatch)              ! dust
  call trace%end()
END SUBROUTINE pre_init

 !===============================================================================
 !> Initialize the modules that are to be included and used
 !===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%forces%init                                         ! forces
  call self%selfgravity%init                                    ! selfg
  call self%particle_solver%init                                ! dust
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Deallocate allocated arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%dealloc')
  call self%forces%dealloc                                      ! forces
  call self%gpatch_t%dealloc
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> The pre_update hierarchy is called before the solver updates
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target:: self
  class(patch_t), pointer:: patch
  !.............................................................................
  call trace%begin ('extras_t%pre_update')
  call self%gpatch_t%pre_update                                 ! keep
  self%gpatch => self                                           ! keep
  patch => self                                                 ! keep
  call self%forces%pre_update                                   ! forces
  call self%selfgravity%pre_update                              ! selfg
  call self%cooling%pre_update                                  ! cool
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%post_update')
  call self%selfgravity%post_update                             ! selfg
  call self%particle_solver%update                              ! dust
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Optionally call specific output modules for the extra features
!===============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
  call self%gpatch_t%output                                     ! keep
END SUBROUTINE output

!===============================================================================
!> Interface to refinement.
!===============================================================================
INTEGER function check_refine (self, patch)
  class(extras_t):: self
  class(extras_t), pointer:: patch
  class(gpatch_t), pointer:: gpatch
  !.............................................................................
  gpatch => patch
  check_refine = -1
END FUNCTION check_refine

!===============================================================================
!> Intercept courant_condition() calls; make sure the gpatch_t call is the last
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(extras_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%courant_condition')
  call self%gpatch_t%courant_condition (detailed_timer)          ! keep
  call self%cooling%courant_condition                            ! cool
  call trace%end()
END SUBROUTINE courant_condition

!===============================================================================
!> Called to allow statistics collection
!===============================================================================
SUBROUTINE statistics (self)
  class(extras_t):: self
  call self%particle_solver%statistics
END SUBROUTINE statistics

END MODULE extras_mod
