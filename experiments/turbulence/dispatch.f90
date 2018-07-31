!> =============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!> =============================================================================
PROGRAM dispatch
  USE setup_mod
  USE dispatcher_mod
  USE cartesian_mod
  USE lagrange_mod
  type (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                               ! Standard setup
  call lagrange%test                            ! Test time interpolation
  call dispatcher%init                          ! Initialize the dispatcher
  call cartesian%init                           ! Initialize the patch list
  call dispatcher%method (cartesian%task_list)  ! Run the dispatcher
  call setup%end                                ! End setup
!> =============================================================================
END PROGRAM dispatch
