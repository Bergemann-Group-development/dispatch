PROGRAM dispatch
  USE mpi_mod
  USE cartesian_mod
  USE io_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call mpi%init
  call io%init
  call cartesian%init ('C-shock')
  call cartesian%task_list%execute     ! update the list = update member patches
  call mpi%end
!> =============================================================================
END PROGRAM dispatch
