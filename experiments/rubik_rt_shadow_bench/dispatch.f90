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
  call rt_setup (cartesian%task_list)              ! add RT tasks to the task list
  call dispatcher%execute (cartesian%task_list) ! Call dispatcher
  print*,'stopping: all tasks finished'
  call setup%end
!===============================================================================
END PROGRAM dispatch
