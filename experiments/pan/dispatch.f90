!> =============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!> =============================================================================
PROGRAM dispatch
  USE setup_mod
  USE cartesian_mod
  USE dispatcher_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                               ! Standard setup (MPI, I/O, scaling, ...)
  call cartesian%init                           ! Initialize the patch list and patches
  call dispatcher%init                          ! Set dispatcher method
  call dispatcher%method (cartesian%task_list)  ! Cáll dispatcher
  call setup%end
!> =============================================================================
END PROGRAM dispatch
