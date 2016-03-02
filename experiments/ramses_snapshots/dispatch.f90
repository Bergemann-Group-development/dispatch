PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE ramses_mod
  type(ramses_t):: ramses
  !.............................................................................
  call io%init
  !call mpi%init
  call ramses%init('ramses_read')       ! initialize the patch list and patches
  call ramses%task_list%execute         ! update the list = update member patches
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
