PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE setup_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init
  call cartesian%init                   ! Initialize the patch list and patches
  call cartesian%task_list%execute      ! update the list = update member patches
  call setup%end
!> =============================================================================
END PROGRAM dispatch
