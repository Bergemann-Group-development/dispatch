PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Make list of links with pointers to patches, and with procedures to initialize
!> and update the patches.  Keep updating until the list is empty
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE io_mod
  USE pde_mod
  type(pde_t):: poisson
  !.............................................................................
  call poisson%init
  call poisson%update
  call poisson%output
!> =============================================================================
END PROGRAM experiment
