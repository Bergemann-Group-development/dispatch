!> =============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!> =============================================================================
PROGRAM dispatch
  USE setup_mod
  USE rubiks_mod
  type (rubiks_t):: rubiks
  !.............................................................................
  call setup%init                       ! Standard setup (MPI, I/O, scaling, ...)
  call rubiks%init_hierarchy            ! Initialize the patch list and patches
  call rubiks%task_list%execute         ! Update task list = update all tasks
  call setup%end
!> =============================================================================
END PROGRAM dispatch
