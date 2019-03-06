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
  !-----------------------------------------------------------------------------
  ! Optional modules start here.  
  !-----------------------------------------------------------------------------
  !USE forces_mod                                        ! forces
  !USE trace_particles_mod                               ! trace_particles
  !USE spitzer_mod
  !USE gravity_mod
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    !type(forces_t):: forces                             ! forces
    !type (spitzer_t) :: spitzer
    !type(gravity_t) :: gravity
    !class(trace_particles_t), pointer:: trace_particles ! trace_particles
    !real, dimension(:,:,:,:), pointer:: v               ! trace_particles
  contains
    procedure:: init
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
  end type
  type(extras_t), public:: extras
CONTAINS

!> =============================================================================
!> Initialize the modules that are to be included and used
!> =============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  !call self%forces%init (self)                          ! forces
  !-----------------------------------------------------------------------------
  ! Allocate a trace_particles data type, connect a pointer, and initialize
  !-----------------------------------------------------------------------------
  !allocate (self%trace_particles)                       ! trace_particles
  !self%connect%trace_particles => self%trace_particles  ! trace_particles
  !call self%trace_particles%add (self%link%task)        ! trace_particles
  !--------- Spitzer conductivity ----------
  ! call self%spitzer%init(self%link)
  !--------- gravitational acceleration ----
  !call self%gravity%init(self%link)
END SUBROUTINE init

!> =============================================================================
!> The pre_update hierarchy is called before the solver updates
!> =============================================================================
SUBROUTINE pre_update (self)
  class(extras_t):: self
  !.............................................................................
  !call self%forces%pre_update (self)                    ! forces
  !call self%trace_particles%update (self)               ! trace_particles
  !call self%spitzer%update                              ! spitzer conductivity
  !call self%gravity%update
END SUBROUTINE pre_update

!> =============================================================================
!> The post_update hierarchy is called after the solver updates
!> =============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  !call self%trace_particles%post_update (self)          ! trace_particles
END SUBROUTINE post_update

!> =============================================================================
!> Optionally call specific output modules for the extra features
!> =============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
END SUBROUTINE output

END MODULE extras_mod
