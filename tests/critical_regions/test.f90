
!===============================================================================
PROGRAM test
  USE io_mod
  USE mpi_mod
  USE omp_timer_mod
  implicit none
  integer, parameter:: n=10000000
  real:: a(n), b(n)
  integer:: i
  !-----------------------------------------------------------------------------
  call io%init
  call mpi%init
  !-----------------------------------------------------------------------------
  call tic
  do i=1,n
    a(i)=i+1
  end do
  call timer%toc("no omp critical",n)
  !-----------------------------------------------------------------------------
  call tic
  do i=1,n
    !$omp critical (test_cr)
    a(i)=i+2
    !$omp end critical (test_cr)
  end do
  call timer%toc("omp critical",n)
  !-----------------------------------------------------------------------------
END PROGRAM
