PROGRAM experiment
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Setup a single patch with a density tophat and call the Poisson solver
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
  USE pde_mod
  type(pde_t):: BE
  !.............................................................................
  call BE%init
  call BE%solve
  call BE%output
!> =============================================================================
END PROGRAM experiment
