!> =============================================================================
!> This test shows that it costs about 45 microseconds to create a task in 
!> gfortran on the T430s. The best performance (on two cores) is with 15-25
!> threads.
!> =============================================================================
PROGRAM test
  USE omp_lib
  implicit none
  integer:: i, j
  real(8):: t0, t1, a
  !$omp parallel
  !$omp master
  t0 = omp_get_wtime()
  do i=1,10
    !$omp task firstprivate(i), private(j,a,t1)
      print *, omp_get_wtime()-t0, i, omp_get_thread_num()
      t1 = omp_get_wtime()
      a = 0.0
      do j=1,100000
        a = a + sin(0.1*j)
      end do
      print *, omp_get_wtime()-t0, i, omp_get_thread_num(), a, omp_get_wtime()-t1
    !$omp end task
    print *, omp_get_wtime()-t0, i
  end do
  !$omp end master
  !$omp end parallel
END PROGRAM test
