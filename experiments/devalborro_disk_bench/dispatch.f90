!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Setup a cylindrical annulus in cylindrical coords. and add a gravitating point
!> mass that represents a massive planet (but has no real surface). A rotating
!> coord. system, centred at the planet's orbital radius, is possible.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
PROGRAM dispatch
  USE setup_mod
  USE mpi_mod
  USE io_mod
  implicit none
  !.............................................................................
  call mpi%init
  call setup%init
  call setup%task_list%execute
  call mpi%end
  print*,'All tasks finished. Goodbye.'
!> =============================================================================
END PROGRAM dispatch
