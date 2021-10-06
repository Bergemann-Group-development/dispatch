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
!>                       \    /        |
!>                       gpatch_t      |
!>                         |     \     |
!>                         |      list_t
!>                         |     /  |
!>                        patch_t   |
!>                            \     |
!>                             link_t
!>                              |
!>                            task_t
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
  USE initial_conditions_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.
  !-----------------------------------------------------------------------------
  USE m3d_rt_mod                                                       ! RT
  USE newton_mod                                                       ! Newton
  USE gravity_mod                                                      ! gravity
  !USE forces_mod                                                       ! forces
  !USE trace_particles_mod                                              ! tracep
  !USE spitzer_mod                                                      ! spitzer
  !USE pic_mod                                                          ! PIC
  !USE cooling_mod                                                      ! PIC
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    class(gpatch_t), pointer:: gpatch => null()                         ! keep
    type(initial_conditions_t):: initial_conditions                     ! keep
    type(gravity_t):: gravity                                          ! gravity
    type(newton_t):: newton                                            ! Newton
    class(m3d_rt_t), pointer:: rt => null()                            ! RT
    !type(forces_t):: forces                                            ! forces
    !class(pic_t), pointer:: pic                                        ! pic
    !type (spitzer_t):: spitzer                                         ! spitzer
    !real, dimension(:,:,:,:), pointer:: v                              ! tracep
    !type(cooling_t):: cooling                                          ! cool
  contains
    procedure, nopass:: cast2extras
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure:: dealloc
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    procedure:: check_refine
    procedure:: nbor_relations
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
  self%gpatch => self                                                   ! keep
  call self%gpatch_t%pre_init                                           ! keep
  allocate(self%rt)                                                     ! RT
  self%rt%patch => self
  call self%rt%pre_init                                                 ! RT
  call self%newton%pre_init  (self)                                     ! Newton
  call self%gravity%pre_init (self)                                      ! gravity
  !call self%forces%init (self)                                          ! forces
  !call self%spitzer%pre_init (self)                                     ! spitzer
  !allocate (self%pic)                                                   ! pic
  !call self%pic%pre_init (self%link)                                    ! pic
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%gpatch_t%init                                               ! keep
  !-----------------------------------------------------------------------------
  ! Optional inits
  !-----------------------------------------------------------------------------
  call self%newton%init                                                 ! Newton
  call self%gravity%init                                                ! gravity
  call self%rt%init                                                     ! RT
  !call self%forces%init (self)                                         ! forces
  !call self%pic%init                                                   ! pic
  !if (self%pic%on) self%do_mhd = self%pic%do_mhd                       ! pic
  !call self%spitzer%init                                                ! spitzer
  ! -- ensure consistency (a value is already set via rt%init)         ! experiment
  !if (associated(self%rt)) then
  !  if (self%rt%on) self%rt%boundaries%gravity = self%gravity%const_grav                 ! experiment
  !end if
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
  call self%gpatch_t%init_task_list (task_list)                        ! keep
  call self%rt%init_task_list (task_list)                              ! RT
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Deallocate allocated arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%dealloc')
  !if (associated(self%rt)) then
  !  call self%rt%dealloc()
  !  deallocate(self%rt)
  !end if
  !call self%forces%dealloc                                             ! forces
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
  call self%newton%pre_update()                                        ! Newton
  call self%gravity%pre_update()                                       ! gravity
  !call self%forces%pre_update                                          ! forces
  !call self%spitzer%pre_update()                                       ! spitzer
  !call self%pic%pre_update()                                           ! pic
  !call self%selfgravity%post_update (self)                             ! selfg
  !call self%cooling%pre_update (self)                                  ! cool
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates, and might be
!> called from either solver_t%update, or an intercepting from mhd_t%post_update
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  !call self%spitzer%post_update
  !call self%forces%post_update                                         ! forces
  !call self%pic%post_update()                                          ! pic
  !call self%selfgravity%post_update (self)                             ! selfg
END SUBROUTINE post_update

!===============================================================================
!> Optionally call specific output modules for the extra features
!===============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
  call self%gpatch_t%output                                            ! keep
  call self%rt%output()                                                ! RT
  !call self%rt%output                                                  ! RT
END SUBROUTINE output

!===============================================================================
!> Interface to refinement
!===============================================================================
INTEGER FUNCTION check_refine (self)
  class(extras_t), target:: self
  !.............................................................................
  check_refine = -1                                                    ! keep
END FUNCTION check_refine

!===============================================================================
!> Cast a generic task_t to extras_t
!===============================================================================
FUNCTION cast2extras (task) RESULT (extras)
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
!> Set default nbor relations, and allow underlings to modify them.  This is
!> called from list_t%init_nbors and from list_t%reset_status()
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(extras_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  !-----------------------------------------------------------------------------
  call trace%begin ('extraswx_t%nbor_relations')
  !-----------------------------------------------------------------------------
  ! First check the main relations, then specific ones
  !-----------------------------------------------------------------------------
  call self%task_t%nbor_relations(another, needed, needs_me, download) ! keep
  select type (another)
  class is (extras_t)
    ! do nothing, all is done in the base task_t layer
  class default
    needed=.false.; needs_me=.false.; download=.false.
  end select
  call trace%end()
END SUBROUTINE nbor_relations

!===============================================================================
!> Optional intercept of courant_condition() call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(extras_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%courant_condition')
  !call self%spitzer%courant_condition (detailed_timer)                 ! spitzer
  if (associated(self%rt)) then
    self%u_max=max(self%u_max,self%rt%u_max/self%rt%rt_freq)
  end if
  call self%gpatch_t%courant_condition (detailed_timer)                ! keep
  call trace%end()
END SUBROUTINE courant_condition

END MODULE extras_mod
