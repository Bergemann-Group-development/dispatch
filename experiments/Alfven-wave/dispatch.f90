PROGRAM dispatch
  USE setup_mod
!  USE mpi_mod
  USE cartesian_mod
  USE io_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                       ! Standard setup (MPI, I/O, scaling, ...)
!  call mpi%init
!  call io%init
  call cartesian%init ('Alfven-wave')
  call cartesian%task_list%execute     ! update the list = update member patches
  call setup%end
!  call mpi%end
!> =============================================================================
END PROGRAM dispatch
