!===============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!===============================================================================
PROGRAM dispatch
  USE io_unit_mod
  USE setup_mod
  USE cartesian_mod
  USE dispatcher_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                               ! Standard setup (MPI, I/O, scaling, ...)
  call dispatcher%init                          ! Set dispatcher method
  call cartesian%init                           ! Initialize the patch list and patches
  call dispatcher%execute (cartesian%task_list) ! call dispatcher
  write (io_unit%mpi,*) 'calling setup%end'
  flush (io_unit%mpi)
  call setup%end
!===============================================================================
END PROGRAM dispatch
