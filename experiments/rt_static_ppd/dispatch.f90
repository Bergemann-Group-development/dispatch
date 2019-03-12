PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Setup a set of patches on a task list, and update until finished.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE setup_mod
  !.............................................................................
  call mpi%init
  call io%init('static RT PPD benchmark')
  call setup%init                   ! Initialize the patch list and patches
  call setup%task_list%execute      ! update the list = update member patches
  call mpi%end
!===============================================================================
END PROGRAM dispatch
