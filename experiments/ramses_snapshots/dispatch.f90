PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE mpi_mod
  USE ramses_mod
  type(ramses_t):: ramses
  !.............................................................................
  call mpi%init
  call ramses%init('ramses_read')       ! initialize the patch list and patches
  call ramses%evaluate
 !call ramses%task_list%execute         ! update the list = update member patches
  call mpi%end
!> =============================================================================
END PROGRAM dispatch
