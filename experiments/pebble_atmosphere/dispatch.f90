PROGRAM dispatch
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE omp_mod
  USE mpi_mod
  USE scaling_mod
  USE eos_mod
  USE download_mod
  USE setup_mod
  USE dispatcher_mod
  !.............................................................................
  call omp%init
  call mpi%init
  call io%init ('pebble_atmosphere')
  call scaling%init
  call eos%init
  call download%init
  call dispatcher%init
  call setup%init                       ! initialize the patch list and patches
  call dispatcher%execute (setup%task_list) ! Run the dispatcher
  call setup%close                      ! close the experiment
  print*,'stopping: all tasks finished'
!===============================================================================
END PROGRAM dispatch
