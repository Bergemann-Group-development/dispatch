!===============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!===============================================================================
PROGRAM dispatch
  USE read_snapshot_mod
  USE setup_mod
  USE restart_mod
  USE cartesian_mod
  USE dispatcher_mod
  USE refine_mod
  type(cartesian_t):: cartesian                 ! Use Cartesian patch arrangement
  class(task_list_t), pointer:: task_list => null()
  !.............................................................................
  call setup%init                               ! Standard setup (MPI, I/O, scaling, ...)
  call dispatcher%init                          ! Set dispatcher method
  call restart%init (task_list)                 ! Possibly initialize task_list
  call cartesian%init (task_list)               ! Initialize the patch list and patches
  !call refine%initial_refine(cartesian%task_list)
  call read_snapshot%init (task_list)           ! Construct ICs
  call dispatcher%execute (task_list)           ! Call dispatcher
  call setup%end
!===============================================================================
END PROGRAM dispatch
