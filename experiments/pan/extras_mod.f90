!> =============================================================================
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
!> =============================================================================

MODULE extras_mod
  !-----------------------------------------------------------------------------
  ! Generic modules, which are always included in the code in any case
  !-----------------------------------------------------------------------------
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE connect_mod
  USE patch_mod
  USE list_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.  
  !-----------------------------------------------------------------------------
  USE forces_mod                                        ! forces
  USE trace_particles_mod                               ! tracep
  USE sinkparticles_mod                                 ! sinkp
  !USE selfgravity_mod                                   ! self-g
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(forces_t):: forces                             ! forces
    type(sinkparticles_t):: sinkparticles               ! sinkp
    !type(selfgravity_t) :: selfgravity                  ! self-g
    class(trace_particles_t), pointer:: trace_particles ! tracep
    real, dimension(:,:,:,:), pointer:: v               ! tracep
    class(list_t), pointer:: task_list
  contains
    procedure:: init
    procedure:: init_task_list
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    procedure:: check_refine
  end type
  type(extras_t), public:: extras
CONTAINS

!> =============================================================================
!> Initialize the modules that are to be included and used
!> =============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call self%forces%init (self%link)                     ! forces
  allocate (self%trace_particles)                       ! tracep
  self%connect%trace_particles => self%trace_particles  ! tracep
  call self%trace_particles%add (self%link)             ! tracep
END SUBROUTINE init

!> =============================================================================
!> =============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(extras_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  self%task_list => task_list
  call self%sinkparticles%init_task_list (task_list)    ! sinkp
END SUBROUTINE init_task_list

!> =============================================================================
!> The pre_update hierarchy is called before the solver updates
!> =============================================================================
SUBROUTINE pre_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%forces%pre_update                           ! forces
  call self%trace_particles%update                      ! tracep
END SUBROUTINE pre_update

!> =============================================================================
!> The post_update hierarchy is called after the solver updates
!> =============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%trace_particles%post_update                 ! tracep
END SUBROUTINE post_update

!> =============================================================================
!> Optionally call specific output modules for the extra features
!> =============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
END SUBROUTINE output

!> =============================================================================
!> Interface to refinement
!> =============================================================================
INTEGER function check_refine (self, patch)
  class(extras_t):: self
  class(patch_t), pointer:: patch
  !.............................................................................
  check_refine = -1
  check_refine = max (check_refine, &
    sinkparticles%check_refine (patch))                 ! sinkp
END FUNCTION check_refine

END MODULE extras_mod
