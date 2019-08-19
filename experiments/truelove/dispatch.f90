PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE setup_mod
  USE cartesian_mod
  USE dispatcher_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init
  call dispatcher%init
  call cartesian%init                   ! Initialize the patch list and patches
  call dispatcher%execute(cartesian%task_list)      ! update the list = update member patches
  call setup%end
!===============================================================================
END PROGRAM dispatch
