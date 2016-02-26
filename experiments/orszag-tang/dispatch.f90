PROGRAM dispatch
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the Poisson solver
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE mpi_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  type (task_list_t), pointer:: task_list
  real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  call io%init
  !call mpi%init
  call cartesian%init('orszag-tang')     ! Initialize the patch list and patches
  task_list => cartesian%task_list      ! for now, just short hand -- could assemble
  task_list%n_tasks = product(cartesian%dims)
  !$omp parallel
  do while (task_list%n > 0)            ! As long as there are patches in the list
    call task_list%update               ! update the list = update member patches
  end do
  !$omp end parallel
  print*,'stopping: all tasks finished'
!> =============================================================================
END PROGRAM dispatch
