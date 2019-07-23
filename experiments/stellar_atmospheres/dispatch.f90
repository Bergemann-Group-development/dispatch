!===============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!===============================================================================
PROGRAM dispatch
  USE read_snapshot_mod
  USE setup_mod
  USE cartesian_mod
  USE dispatcher_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                               ! Standard setup (MPI, I/O, scaling, ...)
  call dispatcher%init                          ! Set dispatcher method
  call cartesian%init                           ! Initialize the patch list and patches
  call read_snapshot%init (cartesian%task_list) ! Construct ICs
  call dispatcher%execute (cartesian%task_list) ! Call dispatcher
  call setup%end
!===============================================================================
END PROGRAM dispatch
