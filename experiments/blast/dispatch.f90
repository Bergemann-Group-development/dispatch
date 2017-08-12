PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE setup_mod
  USE components_mod
  type (components_t):: components
  !.............................................................................
  call setup%init
  call components%init                   ! Initialize the patch list and patches
  call components%task_list%execute      ! update the list = update member patches
  call setup%end
!> =============================================================================
END PROGRAM dispatch
