!===============================================================================
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
  USE rt_mod                                                           ! RT
  USE newton_mod                                                       ! Newton
  USE gravity_mod                                                      ! gravity
  !-----------------------------------------------------------------------------
  implicit none
  private
  type, public, extends(gpatch_t):: extras_t
    type(gravity_t):: gravity                                          ! gravity
    type(newton_t):: newton                                            ! Newton
    class(rt_t), pointer:: rt => null()                                ! RT
    class(gpatch_t), pointer:: patch => null()                         ! keep
    logical :: do_mhd = .true.                                         ! keep
  contains
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
!> Pre-initialize modules, reading parameters, and picking up self%patch
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_init')
  self%patch => self
  allocate(self%rt)                                                    ! RT
  call self%rt%pre_init (self%patch)                                   ! RT
  call self%newton%pre_init (self%patch)                               ! Newton
  call self%gravity%pre_init (self%patch)                              ! gravity
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%gravity%init                                               ! gravity
  call self%newton%init                                                ! Newton
  call self%rt%init_task (self%rt)                                     ! RT
  self%rt%rt_boundaries%grav = self%gravity%const_grav                 ! experiment
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Deallocate allocated arrays
!===============================================================================
SUBROUTINE dealloc (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%dealloc')
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> The pre_update hierarchy is called before the solver updates.
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_update')
  call self%rt%pre_update()                                            ! RT
  call self%newton%pre_update()                                        ! Newton
  call self%gravity%pre_update()                                       ! gravity
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update hierarchy is called after the solver updates
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call self%rt%post_update                                             ! RT
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
END FUNCTION check_refine

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
