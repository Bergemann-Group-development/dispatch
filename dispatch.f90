PROGRAM dispatch
  USE experiment_mod
  USE mpi_mod
  USE io_mod
  implicit none
  class(experiment_t), pointer:: experiment
  !.............................................................................
  call mpi%init
  call io%init
  allocate (experiment)
  call experiment%init
  do while (experiment%task_list%n > 0)
    call experiment%task_list%update
  end do
  call experiment%end
  call mpi%end
END PROGRAM
