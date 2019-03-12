!===============================================================================
!> This module is a placeholder for scene information, solving the problem of
!> making a piece of information available btw any two points in the code, without
!> having to worry about module dependency.  Since this module does not depend
!> on anything else, we can put shared information here
!===============================================================================
MODULE scene_mod
  implicit none
  public
  type, public:: scene_t
    integer:: n_repeat(3) = 1
    integer:: n_patch
    integer:: n_initial = 0
    integer:: n_add = 0
    integer:: n_patch_total = 0
    integer:: patch_n(3) = 0
    real:: f_zoom(3)=1.0
  end type
  type(scene_t):: scene
END MODULE scene_mod
