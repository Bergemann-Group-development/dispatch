PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE dispatcher_mod
  USE cartesian_mod
  USE setup_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init
  call dispatcher%init
  call cartesian%init     
  call dispatcher%execute (cartesian%task_list) ! Call dispatcher
  call setup%end
!===============================================================================
END PROGRAM dispatch
