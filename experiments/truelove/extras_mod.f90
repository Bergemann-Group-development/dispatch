!===============================================================================
!> Template version of a module intended for adding extra features, at a level
!> between the basic task and patch layer, and the layer of solvers.
!>
!> To add extra features to the code, copy this file to experiments/your_dir/,
!> and add the extra dependencies to the Makefile (only) there.
!>
!> When doing to, the sketch of module dependencies below may be useful.
!> In particular, it shows that extras underlings, via gpatch_t, may use the
!> list_t data type (and hence they can access the copy of the task_list
!> pointer made available there).
!>
!> The sketch also illustrates that extras underlings can at most cast a task
!> pointer to a gpatch_t pointer, but not to an extras_t pointer (which is why
!> generic array data needs to be stored in data types below extras_t).
!>
!>                               task_list_t
!>                                / |  |
!>                    experiment_t  |  |
!>                     solver_t     |  |
!>                       mhd_t      |  |
!>                         | refine_t  |
!>                         |  /        |
!>                      extras_t       |
!>                      /    \         |
!>                     /      \        |
!>              extras_1_t extras_2_t  |
!>                 |     \    /        |
!>                 |     gpatch_t      |
!>                 |       |     \     |
!>                 |       |      list_t
!>                 |       |     /  |
!>                 |      patch_t   |
!>                 |      /   \     |
!>                 |     /     link_t
!>                 |    /       |
!>                connect_t   task_t
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
!> This module should generally NOT use the h5_mod (HDF5) module itself. h5_mod
!> should instead be USEd from the actual extra modules that wish to use it to
!> save ad hoc results in a standardized format.
!===============================================================================

MODULE extras_mod
  !-----------------------------------------------------------------------------
  ! Generic modules, which are always included in the code in any case
  !-----------------------------------------------------------------------------
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.  
  !-----------------------------------------------------------------------------
  USE selfgravity_mod
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(selfgravity_t):: selfgravity
  contains
    procedure, nopass:: cast2extras
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    procedure:: check_refine
  end type
  type(extras_t), public:: extras
CONTAINS

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call self%selfgravity%init_extras (self)
END SUBROUTINE init

!===============================================================================
!> The pre_update hierarchy is called before the solver updates
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%selfgravity%pre_update (self)
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%selfgravity%post_update (self)
END SUBROUTINE post_update

!===============================================================================
!> Call the specific output modules of the extra features, as well as the
!> generic output
!===============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
  call self%gpatch_t%output                                     ! keep
END SUBROUTINE output

!===============================================================================
!> Interface to refinement
!===============================================================================
INTEGER function check_refine (self, patch)
  class(extras_t):: self
  class(extras_t), pointer:: patch
  class(gpatch_t), pointer:: gpatch
  !.............................................................................
  gpatch => patch
  check_refine = -1
  !check_refine = max (check_refine, &
  !  sinkparticles%check_refine (gpatch))                        ! sinkp
END FUNCTION check_refine

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2extras (task) RESULT(extras)
  USE task_mod
  class(task_t), pointer:: task
  class(extras_t), pointer:: extras
  !.............................................................................
  select type (task)
  class is (extras_t)
  extras => task
  class default
  nullify(extras)
  call io%abort ('extras_t%cast: failed to cast a task to extras_t')
  end select
END FUNCTION cast2extras

END MODULE extras_mod
