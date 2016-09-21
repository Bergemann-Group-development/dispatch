!> =============================================================================
!> =============================================================================
PROGRAM test
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer:: i, m, n, iostat, it, nit
  real, dimension(:,:), pointer:: dtau, s, qfw
  !-----------------------------------------------------------------------------

  call omp%init
  do while (.true.)
    read (*,*,iostat=iostat) n,m,nit
    if (iostat /= 0) exit
    allocate (dtau(m,n), s(m,n), qfw(m,n))
    dtau = 1.0
    s = 1.0
    qfw = 0.0

    call tic
    do it=1,nit
      call transfer (m, n, dtau, s, qfw)
    end do
    call toc ('transfer', nit*n*m)

    call tic
    do it=1,nit
      call transfer2 (m, n, dtau, s, qfw)
    end do
    call toc ('transfer2', nit*n*m)

    deallocate (dtau, s, qfw)
  end do
END PROGRAM

!> =============================================================================
SUBROUTINE transfer (m, n, dtau, s, qfw)
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

!> =============================================================================
SUBROUTINE transfer2 (m, n, dtau, s, qfw)
  implicit none
  integer m, n
  real, dimension(m,n) :: dtau, s, qfw
  real, dimension(m,n) :: ds1
  real :: ex0, ex1, ex2, d1s, ctau
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
    ex0 = exp(-dtau(j,i))                                           !          1e
    ex1 = 1.-ex0                                                    !       1a
    qfw(j,i) = qfw(j,i-1)*ex0 &
             - ex1*(ds1(j,i+1)*dtau(j,i)+ds1(j,i)*dtau(j,i+1))*ctau &
             + (ds1(j,i+1)-ds1(j,i))*2.*ctau*(ex1-ex0*dtau(j,i))
  end do
  end do                                                            ! 2d 9m 8a 1e
END
