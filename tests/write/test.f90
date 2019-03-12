!===============================================================================
!> Test if printing to a unit outher than stdout is thread-safe -- it appear to
!> be, with both gfortran and ifort.
!===============================================================================
PROGRAM test
  USE mpi
  USE omp_lib
  implicit none
  integer:: i, j
  !-----------------------------------------------------------------------------
  !$omp parallel
  !$omp master
  print *, 'testing with', omp_get_num_threads(), 'threads'
  !$omp end master
  !$omp end parallel
  do j=1,1000
    if (mod(j,100)==0) print *, j/10, '%'
    open (unit=20, file='test.log', form='formatted', status='unknown')
    !$omp parallel
    !$omp do private(i)
    do i=1,10000
      write (20,*) i, omp_get_thread_num(), i+omp_get_thread_num()
    end do
    !$omp end parallel
    close (20)
  end do
  print *, 'Check test.log for confusion'
END PROGRAM
