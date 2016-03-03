PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a set of Cartesian patches on a task list, and update until finished
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE setup_mod
  USE experiment_mod
  USE omp_timer_mod
  !.............................................................................
  call io%init
 !call mpi%init
  call setup%init('planet_atmosphere')  ! initialize the patch list and patches
  call tic
  call setup%task_list%execute          ! execute the list = update member patches
  call toc ('update', n_update)
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
