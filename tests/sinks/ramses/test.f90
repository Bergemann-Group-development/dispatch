PROGRAM test
  USE omp_mod
  USE io_mod
  USE mpi_mod
  USE scaling_mod
  USE sinks_mod
  type(sinks_t) sinks
  !-----------------------------------------------------------------------------
  call omp%init
  call mpi%init
  call io%init
  call scaling%init
  call sinks%init
END PROGRAM test
