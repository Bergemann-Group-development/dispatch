!*******************************************************************************
!> Module with list handling for generic class task_t objects
!*******************************************************************************
MODULE test_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  private
  integer, public:: unit
  integer, public, parameter:: nbuf=100000
  !$omp threadprivate (unit)
  real, public, dimension(nbuf):: buffer, bufs
  !$omp threadprivate (bufs)
  public create
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE create
  implicit none
  integer:: i
  open (unit=10, file='test.tmp', status='unknown')
  close (unit=10, status='delete')
  open (unit=10, file='test.tmp', form='unformatted', status='unknown', &
        access='direct', recl=sizeof(buffer))
  print '(/a)', 'write:'
  call tic
  do i=1,omp%nthreads
    buffer = i
    !print *,i
    write (10,rec=i) buffer
  end do
  close (10)
  call timer%toc ('write', nbuf)
END SUBROUTINE
END MODULE

!===============================================================================
!===============================================================================
PROGRAM test_omp_io
  USE test_mod
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: i

  call omp%init

  call create
  print '(/a)', 'parallel read:'
  call tic
  !$omp parallel
  unit = 100 + omp_mythread
  open (unit=unit, file='test.tmp', form='unformatted', status='old', &
        access='direct', recl=sizeof(bufs))
  read (unit,rec=1+omp%thread) bufs
  print*,unit, omp%thread, omp%mythread(), bufs(1)
  close (unit)
  !$omp end parallel
  call timer%toc ('read', nbuf)

  call create
  print '(/a)', 'critical read:'
  call tic
  !$omp parallel
  !$omp critical
  unit = 100 + omp_mythread
  open (unit=unit, file='test.tmp', form='unformatted', status='old', &
        access='direct', recl=sizeof(bufs))
  read (unit,rec=1+omp%thread) bufs
  print*,unit, omp%thread, omp%mythread(), bufs(1)
  close (unit)
  !$omp end critical
  !$omp end parallel
  call timer%toc ('read', nbuf)

  call create
  print '(/a)', 'serial read:'
  call tic
  do i=1,omp%nthreads
    unit = 100
    open (unit=unit, file='test.tmp', form='unformatted', status='old', &
          access='direct', recl=sizeof(bufs))
    read (unit,rec=i) bufs
    print*, unit, i, omp%mythread(), bufs(1)
    close (unit)
  end do
  call timer%toc ('read', nbuf)

END PROGRAM
