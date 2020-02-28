PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE setup_mod
  USE cartesian_mod
  USE dispatcher_mod
  implicit none
  type (cartesian_t):: cartesian
  !.....................................................................
  call setup%init                                                       ! Standard setup (MPI, I/O, scaling, ...)
  call dispatcher%init                                                  ! Set dispatcher method
  call cartesian%init                                                   ! Initialize the patch list and patches
  call dispatcher%execute (cartesian%task_list)                         ! Call dispatcher
  call setup%end
END PROGRAM dispatch
