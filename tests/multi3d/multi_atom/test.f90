!===============================================================================
!> Test program for running multiple copies of a code
!===============================================================================
PROGRAM multiple_atoms
  USE omp_mod
  call omp%init ()
  !$omp parallel
  call code
  !$omp end parallel
END PROGRAM
