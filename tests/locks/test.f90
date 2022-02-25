PROGRAM test
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  implicit none
  type(lock_t):: lock
  integer(8), parameter:: n=2**12
  integer:: i, j, k

  call init_omp
  call lock%init ('test')

  call tic
  do i=1,n
  do j=1,n
    call omp_set_lock (lock%lock)
    call omp_unset_lock (lock%lock)
  end do
  end do
  call timer%toc ('direct', n**2)

  call tic
  do i=1,n
  do j=1,n
    call lock%set ('test')
    call lock%unset ('test')
  end do
  end do
  call timer%toc ('call', n**2)

  call tic
  do i=1,n
  do j=1,n
    call omp_set_lock (lock%lock)
    !$omp flush
    call omp_unset_lock (lock%lock)
  end do
  end do
  call timer%toc ('direct+flush', n**2)

  call tic
  do i=1,n
  do j=1,n
    call lock%set ('test')
    !$omp flush
    call lock%unset ('test')
  end do
  end do
  call timer%toc ('call+flush', n**2)

  call tic
  do i=1,n
  do j=1,n
    !$omp critical
    k = j
    !$omp end critical
  end do
  end do
  call timer%toc ('critical', n**2)

END PROGRAM
