PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the radiative transfer solver
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE pde_mod
  type(pde_t):: radiation
  !.............................................................................
  call radiation%init
  call radiation%solve
  call radiation%output
!> =============================================================================
END PROGRAM experiment
