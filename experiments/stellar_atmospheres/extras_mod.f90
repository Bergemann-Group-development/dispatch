!===============================================================================
!> Template version of a module intended for adding extra features, at a level
!> between the basic task and patch layer, and the layer of solvers.
!>
!> To add extra features to the code, copy this file to experiments/your_dir/,
!> and add the extra dependencies to the Makefile (only) there.
!>
!> When doing so, the sketch of module dependencies below may be useful.
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
  USE gpatch_mod
  USE list_mod
  USE patch_mod
  USE task_mod
  USE io_mod
  USE trace_mod
  USE connect_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.
  !-----------------------------------------------------------------------------
  USE rt_mod                                                           ! RT
  USE newton_mod                                                       ! Newton
  !USE forces_mod                                                       ! forces
  !USE trace_particles_mod                                              ! tracep
  !USE spitzer_mod                                                      ! spitzer
  USE gravity_mod                                                      ! gravity
  !USE sink_patch_mod                                                   ! sinkp
  !USE particle_solver_mod                                              ! sinkp
  !USE pic_mod                                                          ! PIC
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    !type(forces_t):: forces                                            ! forces
    !class(pic_t), pointer:: pic                                        ! pic
    !type (spitzer_t):: spitzer                                         ! spitzer
    type(gravity_t):: gravity                                          ! gravity
    class(rt_t), pointer:: rt => null()                                ! RT
    type(newton_t):: newton                                            ! Newton
    !type(sinkparticles_t):: sinkparticles                              ! sinkp
    !class(patch_force_t), pointer:: patch_force                        ! sinkp
    !class(trace_particles_t), pointer:: trace_particles                ! tracep
    !real, dimension(:,:,:,:), pointer:: v                              ! tracep
    logical :: do_mhd = .true.
  contains
    procedure, nopass:: cast2extras
    procedure:: pre_init
    procedure:: init
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
!> Pre-initialize the modules that are to be included and used. These routines
!> are used just to get an indication if more memory slots should be used or any
!> specific variables ought to be modified, e.g. index, number of guard zones,
!> periodicity, etc.
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  !call self%forces%init (self)                                         ! forces
  !-----------------------------------------------------------------------------
  ! Allocate a trace_particles data type, connect a pointer, and initialize
  !-----------------------------------------------------------------------------
  allocate(self%rt)                                                    ! RT
  call self%rt%pre_init(self%link)                                     ! RT
  !allocate (self%trace_particles)                                      ! tracep
  !self%connect%trace_particles => self%trace_particles                 ! tracep
  !allocate(self%pic)                                                   ! pic
  !self%connect%pic_patch => self%pic                                   ! pic
  !call self%pic%pre_init(self%link)                                    ! pic
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  !call self%forces%init (self)                                         ! forces
  !-----------------------------------------------------------------------------
  ! Allocate a trace_particles data type, connect a pointer, and initialize
  !-----------------------------------------------------------------------------
  call self%newton%init (self%link)                                    ! Newton
  !allocate (self%trace_particles)                                      ! tracep
  !self%connect%trace_particles => self%trace_particles                 ! tracep
  !call self%trace_particles%add (self%link%task)                       ! tracep
  !call self%pic%init                                                   ! pic
  !if (self%pic%on) self%do_mhd = self%pic%do_mhd                       ! pic
  !call self%spitzer%init(self%link)                                    ! spitzer
  !call self%spitzer%init (self%link)                                   ! spitzer
  call self%gravity%init (self%link)                                   ! gravity
  call self%rt%init_task (self%rt)                                     ! RT
  ! -- ensure consistency (a value is already set via rt%init)         ! experiment
  self%rt%rt_boundaries%grav = self%gravity%const_grav                 ! experiment
  !call sink_patch%init (self%link)                                     ! sinkp
  !allocate (self%patch_force)                                          ! sinkp
  !call self%patch_force%init                                           ! sinkp 
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Deallocate allocated arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%dealloc')
  !call self%forces%dealloc                                             ! forces
  !call trace_particles%dealloc                                         ! sinkp
  !deallocate (self%trace_particles)                                    ! tracep
  !call sink_patch%dealloc                                              ! sinkp
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> The pre_update hierarchy is called before the solver updates. To allow for
!> a general case, with several force data types adding their contributions,
!> the generic patch_t%arrays must (if allocated), be set to zero here, so each
!> force can be added separately, and to prevent these contributions to pile up.
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_update')
  if (allocated (self%force_per_unit_mass)) &                          ! keep
    self%force_per_unit_mass = 0.0                                     ! keep
  if (allocated (self%force_per_unit_volume)) &                        ! keep
    self%force_per_unit_volume = 0.0                                   ! keep
  if (allocated (self%heating_per_unit_mass)) &                        ! keep
    self%heating_per_unit_mass = 0.0                                   ! keep
  if (allocated (self%heating_per_unit_volume)) &                      ! keep
    self%heating_per_unit_volume = 0.0                                 ! keep
  !
  call self%rt%pre_update()                                            ! RT
  call self%newton%pre_update()                                        ! Newton
  !call self%forces%pre_update                                          ! forces
  !call self%trace_particles%update                                     ! tracep
  !call self%spitzer%pre_update()                                       ! spitzer
  call self%gravity%pre_update()                                       ! gravity
  !call particle_solver%force_field (self)                              ! sinkp
  !call self%pic%pre_update()                                           ! pic
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates, and might be
!> called from either solver_t%update, or an intercepting from mhd_t%post_update
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%rt%post_update                                             ! RT
  !call self%forces%post_update                                         ! forces
  !call self%trace_particles%post_update (self)                         ! tracep
  !call self%pic%post_update()                                          ! pic
END SUBROUTINE post_update

!===============================================================================
!> Optionally call specific output modules for the extra features
!===============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
  call self%gpatch_t%output                                            ! keep
  call self%rt%output                                                  ! RT
END SUBROUTINE output

!===============================================================================
!> Interface to refinement
!===============================================================================
INTEGER function check_refine (self, patch)
  class(extras_t):: self
  class(extras_t), pointer:: patch
  !.............................................................................
  check_refine = -1                                                    ! keep
  !check_refine = max (check_refine, &                                  ! sinkp
  !  sinkparticles%check_refine (gpatch))                               ! sinkp
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

!===============================================================================
!> Optional intercept of courant_condition() call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(extras_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%courant_condition')
  call self%rt%courant_condition (detailed_timer)                      ! RT
  call self%gpatch_t%courant_condition (detailed_timer)                ! keep
  call trace%end()
END SUBROUTINE courant_condition

END MODULE extras_mod
