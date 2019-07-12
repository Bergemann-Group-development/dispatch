!===============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!===============================================================================
PROGRAM dispatch
  USE setup_mod
  USE dispatcher_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian                ! Use Cartesian patch arrangement
  !.............................................................................
  call setup%init                               ! Standard setup (MPI, I/O, scaling, ...)
  call dispatcher%init                          ! Initialize the dispatcher
  call cartesian%init                           ! Initialize the task list
  call dispatcher%execute (cartesian%task_list) ! Run dispatcher on the task_list
  call setup%end                                ! End setup
!===============================================================================
END PROGRAM dispatch
