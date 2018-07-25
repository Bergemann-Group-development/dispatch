!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Setup a cylindrical annulus in cylindrical coords. and add a gravitating point
!> mass that represents a massive planet (but has no real surface). A rotating
!> coord. system, centred at the planet's orbital radius, is possible.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
PROGRAM dispatch
  USE setup_mod
  USE mpi_mod
  USE dispatcher_mod
  implicit none
  !.............................................................................
  call mpi%init
  call setup%init
  call dispatcher%init
  call dispatcher%method(setup%task_list)
  call mpi%end
  print*,'All tasks finished. Goodbye.'
!> =============================================================================
END PROGRAM dispatch
