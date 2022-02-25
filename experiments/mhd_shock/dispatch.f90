!===============================================================================
!> Generic DISPATCH main program for Cartesian meshes
!===============================================================================
PROGRAM dispatch
  USE setup_mod
  USE cartesian_mod
  USE dispatcher_mod
  class (cartesian_t):: cartesian
  !.............................................................................
  call setup%init                               ! Standard setup (MPI, I/O, scaling, ...)
  call dispatcher%init                          ! Set dispatcher method
  call cartesian%init                           ! Initialize the patch list and patches
  call dispatcher%execute (cartesian%task_list) ! Cáll dispatcher
  call setup%end
!===============================================================================
END PROGRAM dispatch
