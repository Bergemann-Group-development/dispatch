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
  USE patch_mod
  USE list_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.  
  !-----------------------------------------------------------------------------
  USE initial_conditions_mod                                    ! initial
  USE forces_mod                                                ! forces
  USE trace_particles_mod                                       ! tracep
  USE selfgravity_mod                                           ! selfg
  USE cooling_mod                                               ! cool
  USE lagrangian_mod                                            ! move
  USE sinks_mod                                                 ! sinks
  USE SN_mod                                                    ! SN
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(initial_conditions_t):: initial_conditions             ! initial
    type(forces_t):: forces                                     ! forces
    class(trace_particles_t), pointer:: trace_particles         ! tracep
    real, dimension(:,:,:,:), pointer:: v                       ! tracep
    type(selfgravity_t):: selfgravity                           ! selfg
    type(cooling_t):: cooling                                   ! cool
    type(lagrangian_t):: lagrangian                             ! move
    type(SN_t):: SN                                             ! SN
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
  class(sinks_t), pointer:: sinks => null()                     ! sinks
  type(extras_t), public:: extras
CONTAINS

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_init')
  self%gpatch => self                                           ! keep
  call self%gpatch_t%pre_init                                   ! keep
  call self%forces%pre_init             (self%gpatch)           ! forces
  call self%initial_conditions%pre_init (self%gpatch)           ! initial
  call self%selfgravity%pre_init        (self%gpatch)           ! selfg
  call self%cooling%pre_init            (self%gpatch)           ! cool
  call io%lock%set ('extras_t%pre_init')
  if (.not. associated(sinks)) then                             ! sinks
    allocate (sinks)                                            ! sinks
    call sinks%pre_init                                         ! sinks
  end if                                                        ! sinks
  call io%lock%unset ('extras_t%pre_init')
  call self%SN%pre_init                 (self%gpatch, sinks)    ! sinks
  !allocate (self%trace_particles)                               ! tracep
  !call self%trace_particles%pre_init (self%gpatch)              ! tracep
   !call self%lagrangian%pre_init     (self%gpatch)              ! lagr
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that need the allocated patch to initialize
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%gpatch_t%init                                       ! keep
  call self%initial_conditions%init                             ! initial
  call self%forces%init                                         ! forces
  call self%selfgravity%init                                    ! selfg
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
  !call self%trace_particles%dealloc                             ! tracep
  !deallocate (self%trace_particles)                             ! tracep
  call self%gpatch_t%dealloc
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> The pre_update hierarchy is called before the solver updates
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_update')
  self%gpatch => self                                           ! keep
  call self%gpatch_t%pre_update                                 ! keep
  call self%forces%pre_update                                   ! forces
  call sinks%pre_update                                         ! sinks
  !call self%trace_particles%update                              ! tracep
  call self%selfgravity%pre_update                              ! selfg
  call self%cooling%pre_update                                  ! cool
  call self%SN%pre_update                                       ! SN
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%forces%post_update                                  ! forces
  !call self%trace_particles%post_update                         ! tracep
  call self%selfgravity%post_update                             ! selfg
  !call self%lagrangian%post_update                              ! move
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
INTEGER FUNCTION check_refine (self)
  class(extras_t), target:: self
  !.............................................................................
  gpatch => patch
  check_refine = -1
  !check_refine = max (check_refine, &                            ! sinkp
  !  sink_particles%check_refine ())                              ! sinkp
END FUNCTION check_refine

!===============================================================================
!> Intercept courant_condition() calls; make sure the gpatch_t call is the last
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(extras_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%courant_condition')
  call self%cooling%courant_condition                            ! cool
  call self%SN%courant_condition                                 ! SN
  call self%gpatch_t%courant_condition (detailed_timer)          ! keep (last!)
  call trace%end()
END SUBROUTINE courant_condition

END MODULE extras_mod
