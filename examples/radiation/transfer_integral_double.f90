! $Id$
!***********************************************************************
SUBROUTINE transfer (n, dtau, s, Ifw0, Ire1, q, Ire0, Ifw1)
  implicit none
  integer n
  real, dimension(n) :: dtau, s, q, qfw, qre
  real Ifw0, Ifw1, Ire0, Ire1
!
!  Input: n, dtau(2:n), s(1:n), qfw(1), qre(n)
! Output: q(2:n-1), qfw(n-1), qre(2)
!
! The input is the source functions and opacities at all points, the net
! forward intensity (I-S) at the 1st point and the net reverse intensity
! at the last point.  As a result, we get the net bi-directional intensity
! everywhere, and also the outgoing forward and reverse intensities at 
! points number 1 and np.
!
  real, dimension(n) :: ex0, ex1, ex2, ds1dtau, dsdtau, d2sdtau2
  real ctau, ex
  real(kind=8) dtaui
  integer i
!-----------------------------------------------------------------------
!
  do i=2,n
      dtaui = dtau(i)                                                   ! compute in double prec.
      ex = dexp(-dtaui)                                                 ! compute exp in double prec.
      ex0(i) = ex                                                       ! store in single -- no need for dbl
      ex1(i) = 1d0-ex
      ex2(i) = 1d0-ex-dtaui*ex
      ds1dtau(i) = (s(i)-s(i-1))/dtau(i)                                ! half-point centered dS/dtau
  end do
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    ctau = 1./(dtau(i+1)+dtau(i))
    dsdtau(i) = (ds1dtau(i+1)*dtau(i)+ds1dtau(i)*dtau(i+1))*ctau        ! dS/dtau, 2nd. order accurate
    d2sdtau2(i) = (ds1dtau(i+1)-ds1dtau(i))*2.*ctau                     ! d2S/dtau2, 2nd. order accurate
  end do
!
!  Forward rays (2a+3m)
!
  qfw(1) = Ifw0 - s(1)                                                  ! initial forward Q
  do i=2,n
    qfw(i) = qfw(i-1)*ex0(i) - dsdtau(i)*ex1(i) + d2sdtau2(i)*ex2(i)    ! formal solution forward
  end do
  Ifw1 = qfw(n) + s(n)                                                  ! forward boundary intensity
!
!  Reverse rays (2a+3m)
!
  qre(n) = Ire1 - s(n)                                                  ! initial reverse Q
  do i=n-1,1,-1
    qre(i) = &
       qre(i+1)*ex0(i+1) + dsdtau(i)*ex1(i+1) + d2sdtau2(i)*ex2(i+1)    ! formal solution backwards
    q(i) = 0.5*(qre(i)+qfw(i))                                          ! mean Q
  end do
  Ire0 = qre(2  ) + s(2  )                                              ! reverse boundary intensity
END
