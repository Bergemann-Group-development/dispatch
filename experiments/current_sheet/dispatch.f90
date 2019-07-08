PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the Poisson solver
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE setup_mod
  USE cartesian_mod
  USE dispatcher_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init
  call dispatcher%init                          ! Set dispatcher method
  call cartesian%init('orszag-tang')    ! Initialize the patch list and patches
  call dispatcher%execute (cartesian%task_list) ! Cáll dispatcher
  call setup%end
  print*,'stopping: all tasks finished'
!===============================================================================
END PROGRAM dispatch
