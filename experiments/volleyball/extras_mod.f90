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
  USE newton_mod                                                       ! Newton
  USE forces_mod                                                       ! forces
  USE region_mod
  USE gravity_mod                                                      ! gravity
  USE rad_diffusion_mod
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(newton_t):: newton                                            ! Newton
    type(gravity_t):: gravity                                          ! gravity
    type(rad_diffusion_t):: rad_diffusion
    class(gpatch_t), pointer:: patch => null()
    logical :: do_mhd = .true.
  contains
    procedure, nopass:: cast2extras
    procedure:: pre_init
    procedure:: init
    procedure:: pre_update
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
  class(gpatch_t), pointer:: patch
  patch => self
  self%patch => patch                                                  ! keep
  call self%gpatch_t%pre_init                                          ! keep
  !................................Æ.............................................
  call self%newton%pre_init (patch)                                    ! Newton
  call self%rad_diffusion%pre_init(patch)
  call self%gravity%pre_init (self%patch)                              ! gravity
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%gpatch_t%init                                              ! keep
  !call self%rad_diffusion%init
  call self%gravity%init                                               ! gravity
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> The pre_update hierarchy is called before the solver updates. To allow for
!> a general case, with several force data types adding their contributions,
!> the generic patch_t%arrays must (if associated), be set to zero here, so each
!> force can be added separately, and to prevent these contributions to pile up.
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target :: self
  !.............................................................................
  call trace%begin ('extras_t%pre_update')
  call self%newton%pre_update()                                        ! Newton
  call self%gravity%pre_update()                                       ! gravity
  call self%rad_diffusion%pre_update
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> Interface to refinement
!===============================================================================
INTEGER function check_refine (self)
  class(extras_t), target:: self
  !.............................................................................
  check_refine = -1                                                    ! keep
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
!> Set default nbor relations, and allow underlings to modify them.  This is
!> called from list_t%init_nbors and from list_t%reset_status()
!===============================================================================
SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(extras_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%nbor_relations')
  !
  !-----------------------------------------------------------------------------
  ! First check the main relations, then specific ones
  !-----------------------------------------------------------------------------
  call self%task_t%nbor_relations(another, needed, needs_me, download) ! keep
  !
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
  !call self%rt%courant_condition (detailed_timer)                      ! RT
  call self%rad_diffusion%courant_condition()
  call self%gpatch_t%courant_condition (detailed_timer)                ! keep
  call trace%end()
END SUBROUTINE courant_condition

END MODULE extras_mod
