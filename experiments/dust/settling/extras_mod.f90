!===============================================================================
!===============================================================================

MODULE extras_mod
  !-----------------------------------------------------------------------------
  ! Generic modules, which are always included in the code in any case
  !-----------------------------------------------------------------------------
  USE io_mod
  USE trace_mod
  USE task_mod
  USE gpatch_mod
  USE link_mod
  USE list_mod
  !-----------------------------------------------------------------------------
  ! Optional modules start here.
  !-----------------------------------------------------------------------------
  USE particle_solver_mod
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Note that modules that contain patch-private data, such as random number
  ! states, must be accessed via a private instance in the extras_t data type
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: extras_t
    class(gpatch_t), pointer:: gpatch => null()                         ! keep
    type(particle_solver_t):: particle_solver                           ! part
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure:: pre_update
    procedure:: post_update
    procedure:: check_refine
  end type
CONTAINS

!===============================================================================
!> Pick up a pointer to the gpatch the underlings are used inside, to give them
!> access to all the gpatch properties (including nbor lists, task list, etc)
!===============================================================================
SUBROUTINE pre_init (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_init')
  self%gpatch => self                                                   ! keep
  call self%particle_solver%pre_init (self%gpatch)                      ! part
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the modules that are to be included and used
!===============================================================================
SUBROUTINE init (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%init')
  call self%particle_solver%init                                        ! part
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize particle solver nbor list, after task list is created
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(extras_t):: self
  class(list_t), pointer:: task_list
  class(nbor_t), pointer:: nbor, nbor2
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('extras_t%init_task_list')
  call self%particle_solver%init_task_list (task_list)                  ! part
  call self%gpatch_t%init_task_list (task_list)                         ! keep
  nbor => self%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (extras_t)
      !-------------------------------------------------------------------------
      ! Allocate an nbor, point to the particle_solver_t part of it, and prepend
      !-------------------------------------------------------------------------
      allocate (nbor2)
      nbor2%task => task%particle_solver
      nbor2%next => self%particle_solver%nbor
      self%particle_solver%nbor => nbor2
      !print *, 'extras_t: id, nbor', self%id, task%particle_solver%patch%id
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Do whatever needs to be done before the HD/MHD udpate
!===============================================================================
SUBROUTINE pre_update (self)
  class(extras_t), target:: self
  !.............................................................................
  call trace%begin ('extras_t%pre_update')                             ! keep
  call trace%end()
END SUBROUTINE pre_update

!===============================================================================
!> Do whatever needs to be done after the HD/MHD udpate
!===============================================================================
SUBROUTINE post_update (self)
  class(extras_t):: self
  !.............................................................................
  call trace%begin ('extras_t%post_update')                            ! keep
  call self%particle_solver%update                                     ! part
  call trace%end()
END SUBROUTINE post_update

!===============================================================================
!> Interface to refinement
!===============================================================================
INTEGER function check_refine (self, patch)
  class(extras_t):: self
  class(extras_t), pointer:: patch
  !.............................................................................
  call trace%begin ('extras_t%check_refine')                           ! keep
  check_refine = -1                                                    ! keep
  call trace%end()
END FUNCTION check_refine

END MODULE extras_mod
