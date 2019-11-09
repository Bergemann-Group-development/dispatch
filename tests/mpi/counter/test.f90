!> =============================================================================
!> Test program for remote memory access via MPI_Put, MPI_Get
!> =============================================================================
PROGRAM test_mpi_io
  USE io_mod
  USE mpi_mod
  USE omp_timer_mod
  implicit none
  type(mpi_counter_t):: counter
  integer(4):: n, i
  !-----------------------------------------------------------------------------
  ! Initialize MPI and create a counter
  !-----------------------------------------------------------------------------
  call mpi%init
  n = mpi%size
  call counter%init (n)
  i = counter%update(0)
  call mpi%barrier
  print 1,  mpi%rank, i, wallclock(), counter%i, counter%value, 'initial'
1 format("rank:",i4,3x,"value:",i4,3x,"time:",f10.6,3x,"%i:",i4,3x,"%v:",i4,4x,a)
  call mpi%barrier
  i = counter%update(-1)
  print 1,  mpi%rank, i, wallclock(), counter%i, counter%value
  call mpi%barrier
  i = counter%update(0)
  print 1,  mpi%rank, i, wallclock(), counter%i, counter%value, 'final'
  call mpi%barrier
  call mpi%end
END PROGRAM
