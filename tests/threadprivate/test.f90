!*******************************************************************************
!> Module with list handling for generic class task_t objects
!*******************************************************************************
MODULE test_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  real:: test=0.0
  !$omp threadprivate (test)
END MODULE

!===============================================================================
!===============================================================================
PROGRAM test_omp_io
  USE test_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  real:: b

  call omp%init
  b = 2.
  !$omp parallel firstprivate (b, test)
  test=1.0
  print*, 1, omp%thread, test, b
  test = omp%thread+2.0
  print*, 2, omp%thread, test, b
  !$omp end parallel
END PROGRAM
