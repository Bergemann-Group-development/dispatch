PROGRAM dispatch
  !USE task_list_mod
  USE mpi_mod
  USE cartesian_mod
  USE io_mod
  type (cartesian_t):: cartesian
  type (task_list_t), pointer:: task_list
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  !real:: position(3)=0.0, size(3)=1.0
  !.............................................................................
  print*, 'OMP_GET_NUM_THREADS',OMP_GET_NUM_THREADS()
  call mpi%init
  call io%init
  call cartesian%init ('C-shock')
  call cartesian%task_list%execute	     ! update the list = update member patches
  call mpi%end
!> =============================================================================
END PROGRAM dispatch
