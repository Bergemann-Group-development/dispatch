PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE eos_mod
  USE opacity_mod
  USE setup_mod
  !.............................................................................
  call mpi%init
  call io%init('shadow benchmark')
  call eos%init
  call opacity%init
  call setup%init
  call setup%task_list%execute          ! execute the list = update member patches
  call setup%close                      ! close the experiment
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
