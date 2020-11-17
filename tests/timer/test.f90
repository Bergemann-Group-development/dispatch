!===============================================================================
!> Measure granularity of MPI_Wtime(), system_clock, and omp_get_wtime()
!===============================================================================
PROGRAM test
  USE mpi
  USE io_unit_mod
  USE omp_lib
  USE omp_timer_mod
  USE omp_lock_mod
  implicit none
  real(8):: a, start, time1, time2, time3, end_time, test_time
  integer:: istart, itime, rate, i, j, k, n
  type(lock_t):: lock
  !-----------------------------------------------------------------------------
  call system_clock (count_rate=rate)
  io_unit%input = 1
  call omp_lock%init
  call lock%init
  a = 0.0
  print *, '      n    MPI_Wtime     system_clock  omp_get_wtime'
  do i=1,14
    n = 2**i
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
9 print *,a

  n=0
  test_time=2d0
  end_time = omp_get_wtime() + test_time
  do while (omp_get_wtime() < end_time)
    n=n+1
  end do
  print '(6p,f7.3,a)', test_time/n, ' microsec per call to omp_get_wtime()'

  n=0
  end_time = wallclock() + test_time
  do while (wallclock() < end_time)
    n=n+1
  end do
  print '(6p,f7.3,a)', test_time/n, ' microsec per call to wallclock()'

  n=0
  end_time = omp_get_wtime() + test_time
  do while (omp_get_wtime() < end_time)
    n=n+1
  end do
  print '(6p,f7.3,a)', test_time/n, ' microsec per call to omp_get_wtime()'
 
  omp_lock%on = .true.
  omp_lock%verbose = 0
  call omp_set_num_threads (2)
  start = omp_get_wtime()
  !$omp parallel
  !$omp master
  do i=1,n
    continue
  end do
  !$omp end master
  !$omp barrier
  !$omp end parallel
  test_time = omp_get_wtime() - start
  print '(6p,f7.3,a)', test_time/n, ' microsec per empty loop'

  j = 0
  start = omp_get_wtime()
  !$omp parallel shared(j,n) private(i)
  !$omp master
  do i=1,n
    !$omp atomic update
    j = j+1
  end do
  !$omp end master
  !$omp barrier
  !$omp end parallel
  test_time = omp_get_wtime() - start
  print '(6p,f7.3,a)', test_time/n, ' microsec per omp atomic update'
 
  omp_lock%on = .true.
  omp_lock%verbose = 0
  start = omp_get_wtime()
  !$omp parallel
  !$omp master
  do i=1,n
    call lock%set
    call lock%unset
  end do
  !$omp end master
  !$omp barrier
  !$omp end parallel
  test_time = omp_get_wtime() - start
  print '(6p,f7.3,a)', test_time/n, ' microsec per call lock%set/unset'

END PROGRAM
