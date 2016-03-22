
MODULE patch_mod
  type, extends(tasks_t), public:: patch_t
  contains
    procedure:: init
  end type
CONTAINS

SUBROUTINE init (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  call self%mhd_t%init (self)
  call init_rt_internal (self)
  call init_rt_forward (self)
  call init_rt_reverse (self)
END SUBROUTINE init

SUBROUTINE init_rt_internal (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  call self%task_t%init
  !...
END SUBROUTINE init_rt_internal

SUBROUTINE init_rt_forward (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  call self%task_t%init
  !...
END SUBROUTINE init_rt_forward

SUBROUTINE init_rt_reverse (self)
  class(patch_t):: self
  !-----------------------------------------------------------------------------
  call self%task_t%init
  !...
END SUBROUTINE init_rt_reverse

END MODULE patch_mod
