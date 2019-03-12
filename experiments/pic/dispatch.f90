PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE cartesian_mod
  implicit none
  type (cartesian_t):: cartesian        ! An instance of a Cartesian patch setup
  !.............................................................................
  call mpi%init                         ! Initialize the DISPATCH part of MPI
  call io%init('pic')                   ! Initialize I/O units and log files
  call init_ppcode                      ! Read PPcode parameters and initialize
  call cartesian%init                   ! Initialize the patch list and patches
  call cartesian%task_list%execute      ! update the list = update member patches
  call mpi%end                          ! Finalize MPI
!===============================================================================
END PROGRAM dispatch
