!===============================================================================
!> Example of an extras_mod that creates new tasks, which start operating
!> independently.
!>
!> In this case, the refinement criterion for creating sink particles is
!> checked before each update of the parent patch task.  When the criterion
!> is met, the refinement procedure calls a procedure that creates a new
!> task and adds it to the task list.  The new task should be created with
!> an nbor list that includes the patches that overlap with the accretion
!> sphere.
!>
!> Since the new task gets called by the dispatcher, its update procedure
!> should not be called from here, while other minions to the patch task,
!> such as forcing and cooling, are called directly from here, in pre_update()
!> or post_update() procedures.
!>
!> The call hierarchy creating the sink task and accreting to it is:
!>   task_list%update
!>     refine_t%check_current
!>      refine_t%refine_needed
!>        extras_t%pre_update
!>          sink_t%pre_update
!>            sink_t%check
!>              sink_t%create
!>                sink_t%accrete
!>                sink_t%add_task
!>                  task_list_t%append
!>              sink_t%accrete
!>
!> The extras_mod role in the patch update call hierarchy is:
!>   experiment_t%update
!>     solver_t%update
!>       mhd_t%pre_update
!>         extras_t%pre_update
!>       mhd_t%update
!>       mhd_t%post_update
!>         extras_t%post_update
!>
!> The call hierarchy involved in updating the sink task is:
!>   task_list_t%update
!>     sink_t%update
!===============================================================================

MODULE extras_mod
  !-----------------------------------------------------------------------------
  ! Generic modules, which are always included in the code in any case
  !-----------------------------------------------------------------------------
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE gpatch_mod
  USE patch_mod
  USE list_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.  
  !-----------------------------------------------------------------------------
  USE forces_mod                                                ! forces
  USE selfgravity_mod                                           ! selfg
  USE cooling_mod                                               ! cool
  USE sink_mod                                                 ! sink
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    type(forces_t):: forces                                     ! forces
    type(selfgravity_t):: selfgravity                           ! selfg
    type(cooling_t):: cooling                                   ! cool
    type(sink_t):: sink                                       ! sink
    class(gpatch_t), pointer:: patch => null()                  ! keep
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
  type(extras_t), public:: extras
CONTAINS

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_init')
  self%patch => self
  call self%gpatch_t%pre_init                                   ! keep
  call self%selfgravity%pre_init(self%patch)                    ! selfg
  call self%forces%pre_init(self%patch)                         ! forces
  call self%sink%pre_init (self%patch)                          ! sink
  call trace%end()
END SUBROUTINE pre_init

 !===============================================================================
 !> Initialize the modules that are to be included and used
 !===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%gpatch_t%init                                       ! keep
  call self%forces%init                                         ! forces
  call self%sink%init                                           ! sink
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
  call self%sink%init_task_list (task_list)
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
  call self%patch_t%dealloc
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> The pre_update hierarchy is called before the solver updates
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target:: self
  class(patch_t), pointer:: patch
  !.............................................................................
  call trace%begin ('extras_t%pre_update')
  call self%forces%pre_update                                   ! forces
  call self%selfgravity%pre_update                              ! selfg
  call self%cooling%pre_update                                  ! cool
  call self%sink%pre_update                                     ! sink
  !call self%sink%pre_update                                    ! sink
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%post_update')
  call self%forces%post_update                                  ! forces
  call self%selfgravity%post_update                             ! elfg
  !call self%sink%update                                        ! sink
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Optionally call specific output modules for the extra features
!===============================================================================
SUBROUTINE output (self)
  class(extras_t):: self
  !.............................................................................
  !call self%sink%output                                        ! sink
END SUBROUTINE output

!===============================================================================
!> Interface to refinement.
!===============================================================================
INTEGER FUNCTION check_refine (self)
  class(extras_t), target:: self
  !.............................................................................
  gpatch => patch
  check_refine = -1                                             ! keep
END FUNCTION check_refine

!===============================================================================
!> Intercept courant_condition() calls; make sure the gpatch_t call is the last
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(extras_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%courant_condition')
  call self%patch_t%courant_condition (detailed_timer)           ! keep
  call self%cooling%courant_condition                            ! cool
  call trace%end()
END SUBROUTINE courant_condition

END MODULE extras_mod
