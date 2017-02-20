!> *****************************************************************************
!> Module with list handling for generic class task_t objects
!> *****************************************************************************
MODULE test_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  implicit none
  private
  integer, public:: unit
  integer, public, parameter:: nbuf=100000
  !$omp threadprivate (unit)
  real, public, dimension(nbuf):: buffer, bufs
  !$omp threadprivate (bufs)
END MODULE

!> =============================================================================
PROGRAM test_locks
  USE test_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  implicit none
  integer:: i, n
  real(8):: wc, a

  call omp%init
  call lock%init

  n = 2**22

  a = 0
  wc = wallclock()
  do i=1,n
    a = a+1d0
  end do
  wc = wallclock()-wc
  print*,1e6*wc/n, a

  !$omp parallel private(i,wc)
  wc = wallclock()
  do i=1,n
    !$omp critical (list_cr)
    !$omp end critical (list_cr)
  end do
  wc = wallclock()-wc
  if (omp%thread==0) print*,1e6*wc/n
  !$omp end parallel

  !$omp parallel private(i,wc)
  a = 0.
  wc = wallclock()
  do i=1,n
    !$omp critical (list_cr)
    a = a+1d0
    !$omp end critical (list_cr)
  end do
  wc = wallclock()-wc
  if (omp%thread==0) print*,1e6*wc/n
  !$omp end parallel

  !$omp parallel private(i,wc)
  a = 0.
  wc = wallclock()
  do i=1,n
    call lock%set
    a = a+1d0
    call lock%unset
  end do
  wc = wallclock()-wc
  if (omp%thread==0) print*,1e6*wc/n
  !$omp end parallel

  !$omp parallel private(i,wc)
  a = 0.
  wc = wallclock()
  do i=1,n
    call lock%set
    !$omp flush
    a = a+1d0
    !$omp flush
    call lock%unset
  end do
  wc = wallclock()-wc
  if (omp%thread==0) print*,1e6*wc/n
  !$omp end parallel

END PROGRAM
