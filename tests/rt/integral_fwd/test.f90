!===============================================================================
!===============================================================================
PROGRAM test
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: i, j, k, m, n, iostat, it, nit
  real, dimension(:,:), pointer:: dtau, s, qfw
  !-----------------------------------------------------------------------------

  call omp%init
  n = 28
  do k=1,3
    m = 400
    nit = 20000*20/n
    do j=1,5
      print*,'n,m,nit =', n, m, nit
      allocate (dtau(m,n), s(m,n), qfw(m,n))

      call tic
      do it=1,nit
        dtau = 1.0
        s = it
        qfw = 0.0
        call transfer0 (m, n, dtau, s, qfw)
      end do
      call timer%toc ('transfer             ', nit*n*m)

      call tic
      do it=1,nit
        dtau = 1.0
        s = it
        qfw = 0.0
        call transfer1 (m, n, dtau, s, qfw)
      end do
      call timer%toc ('transfer (vectorized)', nit*n*m)

      call tic
      do it=1,nit
        dtau = 1.0
        s = it
        qfw = 0.0
        call transfer2 (m, n, dtau, s, qfw)
      end do
      call timer%toc ('transfer  (optimized)', nit*n*m)

      deallocate (dtau, s, qfw)
      call tic
      do it=1,nit
        allocate (dtau(m,n), s(m,n), qfw(m,n))
        dtau = 1.0
        s = it
        qfw = 0.0
        call transfer2 (m, n, dtau, s, qfw)
        deallocate (dtau, s, qfw)
      end do
      call timer%toc ('transfer   (allocate)', nit*n*m)
      m = m/2
      nit = nit*2
    end do
    n = n/2
    nit = nit*2
  end do
  print*,''
  print*,'Conclusion: Vectorization over an independent index is a must,'
  print*,'and one needs of the order of 100 cases to vectorize over.'
  print*,''
  print*,'Tiling over a few depth zones is a good thing, because then'
  print*,'the previous depth values of ds1 and qfw are still in cache.'
END PROGRAM

!===============================================================================
SUBROUTINE transfer0 (m, n, dtau, s, qfw)
  implicit none
  integer m, n
  real, dimension(m,n) :: dtau, s, qfw
  real, dimension(n) :: ds1
  real :: ex0, ex1, ex2, d1s, d2s, ctau
  integer :: i, j
  !-----------------------------------------------------------------------------
  do j=1,m
  do i=2,n
    ds1(i) = (s(j,i)-s(j,i-1))/dtau(j,i)                            ! 1d    1a      ! half-point centered dS/dtau
  end do
  !-----------------------------------------------------------------------------
  ! Forward rays (2a+3m)
  !-----------------------------------------------------------------------------
  do i=2,n-1
    ctau = 1./(dtau(j,i+1)+dtau(j,i))                               ! 1d    1a
    d1s = (ds1(i+1)*dtau(j,i)+ds1(i)*dtau(j,i+1))*ctau              !    3m 1a      ! dS/dtau, 2nd. order accurate
    d2s = (ds1(i+1)-ds1(i))*2.*ctau                                 !    2m 1a      ! d2S/dtau2, 2nd. order accurate
    ex0 = exp(-dtau(j,i))                                           !          1e
    ex1 = 1.-ex0                                                    !       1a
    qfw(j,i) = qfw(j,i-1)*ex0 - d1s*ex1 + d2s*(ex1-ex0*dtau(j,i))   !    4m 3a      ! formal solution forward
  end do
  end do                                                            ! 2d 9m 8a 1e
END

!===============================================================================
SUBROUTINE transfer1 (m, n, dtau, s, qfw)
  implicit none
  integer m, n
  real, dimension(m,n) :: dtau, s, qfw
  real, dimension(m,n) :: ds1
  real :: ex0, ex1, ex2, d1s, d2s, ctau
  integer :: i, j
  !-----------------------------------------------------------------------------
  do i=2,n
  do j=1,m
    ds1(j,i) = (s(j,i)-s(j,i-1))/dtau(j,i)                          ! 1d    1a      ! half-point centered dS/dtau
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Forward rays (2a+3m)
  !-----------------------------------------------------------------------------
  do i=2,n-1
  do j=1,m                                                          ! explicit loop
    ctau = 1./(dtau(j,i+1)+dtau(j,i))                               ! 1d    1a
    d1s = (ds1(j,i+1)*dtau(j,i)+ds1(j,i)*dtau(j,i+1))*ctau          !    3m 1a      ! dS/dtau, 2nd. order accurate
    d2s = (ds1(j,i+1)-ds1(j,i))*2.*ctau                             !    2m 1a      ! d2S/dtau2, 2nd. order accurate
    ex0 = exp(-dtau(j,i))                                           !          1e
    ex1 = 1.-ex0                                                    !       1a
    qfw(j,i) = qfw(j,i-1)*ex0 - d1s*ex1 + d2s*(ex1-ex0*dtau(j,i))   !    4m 3a      ! formal solution forward
  end do
  end do                                                            ! 2d 9m 8a 1e
END

!===============================================================================
SUBROUTINE transfer2 (m, n, dtau, s, qfw)
  implicit none
  integer m, n
  real, dimension(m,n) :: dtau, s, qfw
  real, dimension(m,n) :: ds1
  real :: ex0, ex1, ex2, d1s, ctau
  integer :: i, j, k
  !-----------------------------------------------------------------------------
  do i=2,n
  do j=1,m
    ds1(j,i) = (s(j,i)-s(j,i-1))/dtau(j,i)                          ! 1d    1a      ! half-point centered dS/dtau
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Forward rays (2a+3m)
  !-----------------------------------------------------------------------------
  do k=2,n-5,5
  do i=k,k+4
  do j=1,m                                                          ! explicit loop
    ctau = 1./(dtau(j,i+1)+dtau(j,i))                               ! 1d    1a
    ex0 = exp(-dtau(j,i))                                           !          1e
    ex1 = 1.-ex0                                                    !       1a
    qfw(j,i) = qfw(j,i-1)*ex0 &
             - ex1*(ds1(j,i+1)*dtau(j,i)+ds1(j,i)*dtau(j,i+1))*ctau &
             + (ds1(j,i+1)-ds1(j,i))*2.*ctau*(ex1-ex0*dtau(j,i))
  end do
  end do                                                            ! 2d 9m 8a 1e
  end do                                                            ! 2d 9m 8a 1e
END
