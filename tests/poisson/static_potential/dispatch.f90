!> =============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!> =============================================================================
PROGRAM dispatch
  USE setup_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                       ! Standard setup (MPI, I/O, scaling, ...)
  call cartesian%init                   ! Initialize the patch list and patches
  call cartesian%task_list%execute      ! Update task list = update all tasks
  call setup%end
!> =============================================================================
END PROGRAM dispatch
