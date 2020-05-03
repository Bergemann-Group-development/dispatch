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
!>                       \    /   |
!>                      gpatch_t  |
!>                       patch_t  |
!>                       /    \   |
!>                    task_t  connect_t
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
!>
!> The call hierarchy that initializes things runs in this order (using as an
!> example patch creation by the cartesian_mod module)
!>
!> cartesian%init
!>   allocate ...
!>   experiment_t%init
!>     solver_t%init
!>       mhd_t%pre_init
!>         extras_t%pre_init
!>       mhd_t%init
!>         extras_t%init
!>
!> The corresponding update call hierarchy is
!>
!> task_list%update
!>   experiment_t%update
!>     solver_t%update
!>       mhd_t%pre_update
!>         extras_t%pre_update
!>       mhd_t%update
!>       mhd_t%post_update
!>         extras_t%post_update
!>
!> The general convention is that the solver_t layer (which is there basically
!> to allow the actual solver layer to have any module name) calls the actual
!> solver layers' procedures, but that layer may choose to not have for example
!> a pre_ and post_update, so those calls "fall through" to the extras_t layer.
!>
!> If the solver layer implements pre_ and post_ calls, it has the responsibility
!> to call the the extras_ layer.
!>
!> Step 1: We refer to the static mesa instance (not the sinks instance), and
!>         check that only one sink is created
!> Step 2: We refer to the mesa instance for creation, but to the mesa nbor task
!>         for anything that requires access to the individual sink particle
!===============================================================================

MODULE extras_mod
  !-----------------------------------------------------------------------------
  ! Generic modules, which are always included in the code in any case
  !-----------------------------------------------------------------------------
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE gpatch_mod
  USE patch_mod
  USE list_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.  
  !-----------------------------------------------------------------------------
  USE forces_mod                                                ! forces
  USE selfgravity_mod                                           ! selfg
  USE cooling_mod                                               ! cool
  USE sinks_mod                                                 ! sinks
  USE mesa_mod                                                  ! mesa
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
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure:: dealloc
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    procedure:: check_refine
    procedure:: courant_condition
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
  call self%selfgravity%pre_init(self%gpatch)                           ! selfg
  call self%cooling%pre_init(self%gpatch)                               ! cool
  call self%forces%pre_init(self%gpatch)
  call mesa%pre_init (self%gpatch)                              ! sinks
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
  call mesa%init                                                ! sinks
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> This intercept is the correct place to add calls to extras 'underlings' that
!> need information both about the corresponding MHD patch (self), AND the task
!> list.  The procedure is called for every task on the task list, just before
!> task list updates start.
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(extras_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('extras_t%init_task_list')
  call self%gpatch_t%init_task_list (task_list)                 ! keep
  call trace%end ()
END SUBROUTINE init_task_list

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
  call self%forces%pre_update                                   ! forces
  call self%selfgravity%pre_update                              ! selfg
  call self%cooling%pre_update                                  ! cool
  call mesa%pre_update (self%gpatch)                            ! sinks
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%post_update')
  call self%forces%post_update                                  ! forces
  call self%selfgravity%post_update                             ! elfg
  call mesa%post_update (self%gpatch)                           ! sinks
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
INTEGER FUNCTION check_refine (self, patch)
  class(extras_t):: self
  class(extras_t), pointer:: patch
  class(gpatch_t), pointer:: gpatch
  !.............................................................................
  gpatch => patch
  call sinks%check_refine (gpatch)                              ! sinks
  check_refine = -1                                             ! keep
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

END MODULE extras_mod
