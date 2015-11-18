PROGRAM dispatch
  USE experiment_mod
  USE mpi_mod
  USE io_mod
  implicit none
  class(experiment_list_t), pointer:: experiment_list
  !.............................................................................
  call mpi%init
  call io%init
  allocate (experiment_list)
  call experiment_list%init ('experiment_list')
  do while (experiment_list%n > 0)
    call experiment_list%update
  end do
  call experiment_list%output
  call mpi%end
END PROGRAM
