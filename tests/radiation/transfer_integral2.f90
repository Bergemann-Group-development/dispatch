! $Id$
!***********************************************************************
SUBROUTINE transfer2 (M, N, dtau, S, Ifw0, Irv1, Q, Irv0, Ifw1)
  implicit none
  integer m, n
  real, dimension(m,n) :: dtau, s, q, qfw, qre
  real, dimension(m)   :: Ifw0, Ifw1, Irv0, Irv1
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
! so dS/dtau is assumed to be linear in tau.  The analytic solution may
! be written as a linear combination of the formal solution for constant
! and linear source functions.
!
  real, dimension(m,n) :: ex0, ex1, ex2, ds1, d1s, d2s
  real, dimension(m)   :: ctau
  integer i
!-----------------------------------------------------------------------
!
  do i=2,n
      ex0(:,i) = exp(-dtau(:,i))                        !          1e
      ex1(:,i) = 1.-ex0(:,i)                            !       1a
      ex2(:,i) = ex1(:,i)-dtau(:,i)*ex0(:,i)            !    1m 1a
      ds1(:,i) = (s(:,i)-s(:,i-1))/dtau(:,i)            ! 1d    1a      ! half-point centered dS/dtau
  end do
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    ctau = 1./(dtau(:,i+1)+dtau(:,i))                   ! 1d    1a
    d1s(:,i) = &
      (ds1(:,i+1)*dtau(:,i)+ds1(:,i)*dtau(:,i+1))*ctau  !    3m 1a      ! dS/dtau, 2nd. order accurate
    d2s(:,i) = (ds1(:,i+1)-ds1(:,i))*2.*ctau            !    2m 1a      ! d2S/dtau2, 2nd. order accurate
  end do

    d2s(:,1) = d2s(:,2)                                                 ! assume constant d2S/dtau2
    d1s(:,1) = d1s(:,2) - dtau(:,2)*d2s(:,2)                            ! extrapolate to 1
    d2s(:,n) = d2s(:,n-1)                                               ! assume constant d2S/dtau2
    d1s(:,n) = d1s(:,n-1) + dtau(:,n)*d2s(:,n)                          ! extrapolate to n
    qfw(:,1) = Ifw0(:) - s(:,1)                                         ! initial forward Q
!
!  Forward rays (2a+3m)
!
  do i=2,n
    qfw(:,i) = qfw(:,i-1)*ex0(:,i) &
           - d1s(:,i)*ex1(:,i) + d2s(:,i)*ex2(:,i)      !    3m 2a      ! formal solution forward
  end do
    Ifw1(:) = qfw(:,n) + s(:,n)                                         ! forward boundary intensity
    qre(:,n) = Irv1(:) - s(:,n)                                         ! initial reverse Q
    q(:,n) = 0.5*(qre(:,n)+qfw(:,n))                                    ! mean Q
!
!  Reverse rays (2a+3m)
!
  do i=n-1,1,-1
    qre(:,i) = qre(:,i+1)*ex0(:,i+1) &
           + d1s(:,i)*ex1(:,i+1) + d2s(:,i)*ex2(:,i+1)  !    3m 2a      ! formal solution backwards
    q(:,i) = 0.5*(qre(:,i)+qfw(:,i))                    !    1m 1a      ! mean Q
  end do
    Irv0(:) = qre(:,1) + s(:,1)                                         ! reverse boundary intensity

END                                                     ! 2d 13m 11a 1e ! 26f + 1e
