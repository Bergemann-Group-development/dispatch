!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for HD task
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE hd_mod
  USE patch_mod
  implicit none
  type, extends(patch_t):: hd_t
    real:: nu(4)
  end type
END MODULE
