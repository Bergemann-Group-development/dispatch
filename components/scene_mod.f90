!> =============================================================================
!> This module is a placeholder for scene information, solving the problem of
!> making a piece of information available btw any two points in the code, without
!> having to worry about module dependency.  Since this module does not depend
!> on anything else, we can put shared information here
!> =============================================================================
MODULE scene_mod
  implicit none
  public
  type, public:: scene_t
    integer:: n_repeat = 1
    integer:: n_patch
  end type
  type(scene_t):: scene
END MODULE scene_mod
