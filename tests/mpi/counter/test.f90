!> =============================================================================
!> Test program for counters
!> =============================================================================
PROGRAM test_counters
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE counters_mod
  implicit none
  type(mpi_counter_t):: mpi_counter
  integer(4):: ni, nj, i, j, k
  type(counter_t):: counter
  !-----------------------------------------------------------------------------
  ! Initialize MPI and create a counter
  !-----------------------------------------------------------------------------
  call mpi%init
  call io%init
  nj = 10
  ni = 4
  !$omp parallel shared(i,j)
  call mpi_counter%init (0)
  !$omp single
  k = mpi_counter%update(0)
  !$omp end single
  !-----------------------------------------------------------------------------
  ! The intent is that each rank should perform the loop ni times, and that
  ! all threads in each process should contribute to increasing "j".
  !-----------------------------------------------------------------------------
  i = 0
  do while (i < ni)
    !$omp single
    call counter%init (0)
    j = counter%update(0)
    write (io_unit%log,*) wallclock(), omp%thread, 'do', mpi%size, i, j
    !$omp end single nowait
    do while (j < nj)
      write (io_unit%log,1)  mpi%rank, omp%thread, i, j, wallclock()
      1 format("rank, thread:",2i4,3x,"values:",2i4,3x,"time:",f10.6)
      !$omp critical
      if (j < nj) then
        j = counter%update(1)
        if (mod(j,nj) == 0) then
          call mpi%barrier ('before')
          i = i+1
          write (io_unit%log,*) 'i, rank, thread =', &
            i, mpi%rank, omp%thread, wallclock()
          call mpi%barrier ('after')
        end if
      end if
      !$omp end critical
    end do
    write (io_unit%log,1) wallclock(), mpi%rank, omp%thread, i, j
  end do
  !$omp barrier
  !$omp end parallel
  write (io_unit%log,*) wallclock(), mpi%rank, 'arrived'
  call mpi%barrier ('end')
  call mpi%end
END PROGRAM
