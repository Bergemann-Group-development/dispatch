PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the Poisson solver
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE pde_mod
  type(pde_t):: poisson
  !.............................................................................
  call poisson%init
  call poisson%solve
  call poisson%output
!> =============================================================================
END PROGRAM experiment
