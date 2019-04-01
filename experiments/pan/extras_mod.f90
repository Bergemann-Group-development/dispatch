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
!===============================================================================

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
  USE forces_mod                                                ! forces
  USE trace_particles_mod                                       ! tracep
  USE sink_patch_mod                                            ! sinkp
  USE particle_solver_mod                                       ! sinkp
  !USE selfgravity_mod                                           ! selfg
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(forces_t):: forces                                     ! forces
    class(trace_particles_t), pointer:: trace_particles         ! tracep
    real, dimension(:,:,:,:), pointer:: v                       ! tracep
    !type(selfgravity_t) :: selfgravity                         ! selfg
  contains
    procedure:: init
    procedure:: dealloc
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
  call trace%begin ('extras_t%init')
  call self%forces%init (self%link)                             ! forces
  allocate (self%trace_particles)                               ! tracep
  self%connect%trace_particles => self%trace_particles          ! tracep
  call self%trace_particles%add (self%link)                     ! tracep
  call sink_patch%init                                          ! sinkp
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Deallocate allocated arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%dealloc')
  call self%forces%dealloc                                      ! forces
  call self%trace_particles%dealloc                             ! sinkp
  deallocate (self%trace_particles)                             ! tracep
  call sink_patch%dealloc                                       ! sinkp
  call self%gpatch_t%dealloc
  call trace%end
END SUBROUTINE dealloc

!===============================================================================
!> The pre_update hierarchy is called before the solver updates
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%forces%pre_update                                   ! forces
  call self%trace_particles%update                              ! tracep
  call particle_solver%force_field (self)                       ! sinkp
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%forces%post_update                                  ! forces
  call self%trace_particles%post_update                         ! tracep
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
!> Interface to refinement. This uses the static instance of sink_patch, where
!> the check_refine may cause creation of dynamic sink_patch_t instances.
!===============================================================================
INTEGER function check_refine (self, patch)
  class(extras_t):: self
  class(extras_t), pointer:: patch
  class(gpatch_t), pointer:: gpatch
  !.............................................................................
  gpatch => patch
  check_refine = -1
  check_refine = max (check_refine, &
    sink_patch%check_refine (gpatch))                            ! sinkp
END FUNCTION check_refine

END MODULE extras_mod
