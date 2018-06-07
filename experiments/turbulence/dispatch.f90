!> =============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!> =============================================================================
PROGRAM dispatch
  USE setup_mod
  USE dispatcher_mod
  USE cartesian_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                               ! Standard setup
  call dispatcher%init                          ! Initialize the dispatcher
  call cartesian%init                           ! Initialize the patch list
  call dispatcher%method (cartesian%task_list)  ! Run the dispatcher
  call setup%end                                ! End setup
!> =============================================================================
END PROGRAM dispatch
