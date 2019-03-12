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
!>                      gpatch_t
!>                      extras_t
!>                        /
!>                     extras_12_t
!>                     /      \
!>                extras_1_t extras_2_t
!>                       \  /
!>                      patch_t
!>                       task_t
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
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE gpatch_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.  
  !-----------------------------------------------------------------------------
  USE selfgravity_mod
  !-----------------------------------------------------------------------------
  implicit none
  private
  type, public, extends(gpatch_t):: extras_t
    type(selfgravity_t):: selfgravity
  contains
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
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
END SUBROUTINE output

END MODULE extras_mod
