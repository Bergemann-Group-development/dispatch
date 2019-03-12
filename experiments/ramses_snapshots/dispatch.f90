PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE setup_mod
  USE ramses_mod
  USE dispatcher_mod
  !.............................................................................
  call setup%init                               ! initialize MPI, I/O, scaling ...
  call dispatcher%init                          ! initialize dispatcher
  call ramses%init('ramses_read')               ! initialize the task list
  call dispatcher%execute (ramses%task_list)    ! execute the task list
  call setup%end                                ! finalize
!===============================================================================
END PROGRAM dispatch
