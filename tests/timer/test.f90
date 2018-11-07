!> =============================================================================
!> Measure granularity of MPI_Wtime(), system_clock, and omp_get_wtime()
!> =============================================================================
PROGRAM test
  USE mpi
  USE omp_lib
  implicit none
  real(8):: a, start, time1, time2, time3
  integer:: istart, itime, rate, i, j, k, n
  !-----------------------------------------------------------------------------
  call system_clock (count_rate=rate)
  a = 0.0
  print *, '      n    MPI_Wtime     system_clock  omp_get_wtime'
  do i=1,14
    n = 2**i
    do k=1,3
      start = MPI_Wtime()
      do j=1,n
        a=a+sin(j*0.1)
      end do
      time1 = MPI_Wtime()-start

      call system_clock (count=istart)
      start = real(istart,kind=8)/real(rate,kind=8)
      do j=1,n
        a=a+sin(j*0.1)
      end do
      call system_clock (count=itime)
      time2 = real(itime-istart,kind=8)/real(rate,kind=8)

      start = omp_get_wtime()
      do j=1,n
        a=a+sin(j*0.1)
      end do
      time3 = omp_get_wtime()-start

      print '(i8,3f14.8)', n, time1, time2, time3
    end do
  end do
  print *,a
END PROGRAM
