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
  do while (experiment_list%size > 0)
    print*,experiment_list%name
    print*,experiment_list%ready_queue%name
    call experiment_list%update
  end do
  call mpi%end
END PROGRAM
