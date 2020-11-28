!===============================================================================
!> Extras module stripped down to handle only sink particles
!>
!> SINKS:
!> There may be a large number of sinks, which are tasks that occur in the
!> task list, and whose positions and properties are updated by their own
!> update procedure.  As part of the sinks_t%update, their nbor lists are
!> kept up-to-date, which means each patch that has overlapping sink tasks
!> have them in their nbor linsts.
!>
!> Call hierarchy:
!>
!> solver_t%init
!>   extras_t%pre_init
!>   extras_t%init
!>
!> task_list%update
!>   sinks_t%update
!>   solver_t%update
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
  !-----------------------------------------------------------------------------
  ! Optional modules start here.
  !-----------------------------------------------------------------------------
  USE sinks_mod                                                         ! sinks
  USE cooling_mod                                                       ! cool
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(cooling_t):: cooling                                           ! cool
    class(patch_t), pointer:: patch => null()
  contains
    procedure, nopass:: cast2extras
    procedure:: pre_init
    procedure:: init
    !procedure:: dealloc
    procedure:: pre_update
    procedure:: post_update
    procedure:: output
    procedure:: courant_condition
    procedure:: check_refine
    procedure:: nbor_relations
  end type
  type(extras_t), public:: extras
CONTAINS

!===============================================================================
!> Pre-initialize the modules that are to be included and used. The %pre_init
!> hierarchy stops here; there is no gpatch_t%pre_init() procedure
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_init')
  self%patch => self
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%cooling%init                                                ! cool
  call trace%end ()
END SUBROUTINE init

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
  call self%cooling%pre_update (self)                                   ! cool
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates, and might be
!> called from either solver_t%update, or an intercepting from mhd_t%post_update
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%post_update')
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Optionally call specific output modules for the extra features
!===============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%output')
  call self%gpatch_t%output                                             ! keep
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Optional intercept of courant_condition() call
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(extras_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%courant_condition')
  call self%gpatch_t%courant_condition (detailed_timer)                 ! keep
  call self%cooling%courant_condition (self)                            ! cool
  call trace%end ()
END SUBROUTINE courant_condition

!===============================================================================
!> Interface to refinement
!===============================================================================
INTEGER function check_refine (self, patch)
  class(extras_t):: self
  class(extras_t), pointer:: patch
  !.............................................................................
  check_refine = -1                                                     ! keep
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
  call trace%begin ('extras_t%nbor_relations')
  call self%task_t%nbor_relations(another, needed, needs_me, download)  ! keep
  call trace%end ()
END SUBROUTINE nbor_relations

END MODULE extras_mod
