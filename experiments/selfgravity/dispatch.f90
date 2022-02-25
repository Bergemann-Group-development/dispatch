PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the Poisson solver
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE task_list_mod
  class (task_list_t):: task_list
  real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  call task_list%init ('selfgravity')   ! Initialize the patch list and patches
  do while (task_list%n > 0)            ! As long as there are patches in the list
    call task_list%update               ! update the list = update member patches
  end do
  print*,'stopping: all tasks finished'
!===============================================================================
END PROGRAM dispatch
