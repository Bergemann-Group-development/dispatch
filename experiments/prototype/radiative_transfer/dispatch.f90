PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Rubiks patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE setup_mod
  !.............................................................................
  call io%init ('radiative_transfer')
 !call mpi%init
  call setup%init                       ! initialize the patch list and patches
  call setup%task_list%execute          ! execute the list = update member patches
  call setup%close                      ! close the experiment
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
