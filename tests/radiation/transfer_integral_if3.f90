! $Id$
!***********************************************************************
SUBROUTINE transfer (N, dtau, S, Ifw0, Irv1, Q, Irv0, Ifw1)
  implicit none
  integer n
  real, dimension(n) :: dtau, s, q, qfw, qre
  real Ifw0, Ifw1, Irv0, Irv1
!
!  Input: N, dtau(2:n), S(1:n), Ifw(1), Irv(n)
! Output: Q(1:n-1), qfw(n-1), qre(2)
!
! The input is the source functions and opacities at all points, the net
! forward intensity (I-S) at the 1st point and the net reverse intensity
! at the last point.  As a result, we get the net bi-directional intensity
! everywhere, and also the outgoing forward and reverse intensities at 
! points number 1 and np.
!
! The equations solved are Q' = Q - dS/dtau in the forward direction and
! Q' = -Q + dS/dtau in the reverse direction.  S is assumed to be quadratic,
! do dS/dtau is assumed to be linear in tau.  The analytic solution may
! be written as a linear combination of the formal solution for constant
! and linear source functions.
!
  real, dimension(n) :: ex0, ex1, ex2, ds1dtau, dsdtau, d2sdtau2
  real ctau, ex
  integer i
!-----------------------------------------------------------------------
!
  do i=2,n
    if (dtau(i) < 0.01) then
      ex0(i) = 1.-dtau(i)
      ex1(i) = dtau(i)
      ex2(i) = ex1(i)-dtau(i)*ex0(i)
    else if (dtau(i) > 10.) then
      ex0(i) = 0.
      ex1(i) = 1.
      ex2(i) = 1.
    else
      ex0(i) = exp(-dtau(i))
      ex1(i) = 1.-ex0(i)
      ex2(i) = ex1(i)-dtau(i)*ex0(i)
    end if
    ds1dtau(i) = (s(i)-s(i-1))/dtau(i)                                  ! half-point centered dS/dtau
  end do
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    ctau = 1./(dtau(i+1)+dtau(i))
    dsdtau(i) = (ds1dtau(i+1)*dtau(i)+ds1dtau(i)*dtau(i+1))*ctau        ! dS/dtau, 2nd. order accurate
    d2sdtau2(i) = (ds1dtau(i+1)-ds1dtau(i))*2.*ctau                     ! d2S/dtau2, 2nd. order accurate
  end do
  d2sdtau2(1) = d2sdtau2(2)                                             ! assume constant d2S/dtau2
  dsdtau(1) = dsdtau(2) - dtau(2)*d2sdtau2(2)                           ! extrapolate to 1
  d2sdtau2(n) = d2sdtau2(n-1)                                           ! assume constant d2S/dtau2
  dsdtau(n) = dsdtau(n-1) + dtau(n)*d2sdtau2(n)                         ! extrapolate to n
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
  qre(n) = Irv1 - s(n)                                                  ! initial reverse Q
  q(n) = 0.5*(qre(n)+qfw(n))                                            ! mean Q
  do i=n-1,1,-1
    qre(i) = &
       qre(i+1)*ex0(i+1) + dsdtau(i)*ex1(i+1) + d2sdtau2(i)*ex2(i+1)    ! formal solution backwards
    q(i) = 0.5*(qre(i)+qfw(i))                                          ! mean Q
  end do
  Irv0 = qre(1) + s(1)                                                  ! reverse boundary intensity
END
