!==============================================================================
PROGRAM test
  USE dll_mod
  USE omp_mod
  USE io_unit_mod
  USE omp_lock_mod
  type(dll_t):: dll
  !------------------------------------------------------------------------------
  call omp%init                         ! OpenMP startup
  call omp_lock%init ('setup')          ! needs omp, io%init
  io_unit%do_trace = .false.
  io_unit%omp_trace = .false.
  call dll%test1
  call dll%test2
  call dll%test3
END PROGRAM test
