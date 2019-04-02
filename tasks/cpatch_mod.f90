!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Cartesian patch for cylinder_t arrangement
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE cpatch_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE solver_mod
  implicit none
  private
  type, public, extends(solver_t):: cpatch_t
    real(8):: z                 ! vertical position
    real(8):: r                 ! radial position
    real(8):: omega(3)          ! angular velocity of mid point
  contains
    procedure:: init
    procedure, nopass:: cast2cpatch
  end type
  type(cpatch_t), public:: cpatch
CONTAINS

!===============================================================================
!> Initialize Cartesian cylinder patch
!===============================================================================
SUBROUTINE init (self)
  class(cpatch_t):: self
  !.............................................................................
  !-----------------------------------------------------------------------------
  call trace%begin ('cpatch_t%init')
  call self%solver_t%init                       ! basic patch initialization
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Cast a generic task_t to patch_t
!===============================================================================
FUNCTION cast2cpatch (task) RESULT(patch)
  class(task_t), pointer:: task
  class(cpatch_t), pointer:: patch
  !.............................................................................
  select type (task)
  class is (cpatch_t)
  patch => task
  class default
  nullify(patch)
  call io%abort ('patch_t%cast: failed to cast a task to patch_t')
  end select
END FUNCTION cast2cpatch

END MODULE cpatch_mod
