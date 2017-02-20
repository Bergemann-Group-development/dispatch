!> =============================================================================
!> This test shows how to parallelize over a loop.
!> The best performance (on two cores) is with 5-15 threads.
!> =============================================================================
PROGRAM test
  USE omp_lib
  implicit none
  integer:: i, j
  real(8):: t0, t1, a
  !$omp parallel
  t0 = omp_get_wtime()
  !$omp do private(i,j,a,t1)
  do i=1,10
    print *, omp_get_wtime()-t0, i, omp_get_thread_num()
    t1 = omp_get_wtime()
    a = 0.0
    do j=1,100000
      a = a + sin(0.1*j)
    end do
    print *, omp_get_wtime()-t0, i, omp_get_thread_num(), a, omp_get_wtime()-t1
  end do
  !$omp end do
  !$omp end parallel
END PROGRAM test
