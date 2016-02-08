!***********************************************************************
SUBROUTINE transfer (N, dtau, S, Ifw1, IrvN, Q, Irv1, IfwN)
!
!  Solve the transfer equation, given optical depth increments dtau, 
!  source function S, infalling intensity Ifw1 at the first point and
!  IrvN at the last point.
!
!  This version works with Q = P - S, where P is Feautriers P (the average 
!  of in- and out-going specific intensity) and S is the source function.
!  The equation solved is Q" = Q - S", with boundary conditions 
!  Q' - Q = S - S' - Ifw1 at the first point and Q' + Q = -(S' + S) + Irv1
!  at the last point
!
!  The choice of equations produces a high precision in Q at large optical
!  depths (where P - S would have round off), while at the same time
!  retaining the precision at small optical depths (see below).
!
!  The routine returns the net heating contribution Q at each point, and 
!  the outgoing intensity at the first point.
!
!  This is a second order version.  For simulations, with rapidly
!  varying absorption coefficients and source functions, this is to be
!  preferred over spline and Hermitean versions because it is positive
!  definite, in the sense that a positive source function is guaranteed
!  to result in a positive average intensity p.  Also, the flux
!  divergence is exactly equal to q, for the conventional definition
!  of the flux.
!
!  "Steins trick" is used; storing the sum of the three elements in each
!  row of the tridiagonal equation system, instead of the diagonal ele-
!  ment.  This allows accurate solutions to be otained for arbitrarily
!  small optical depths.  Without this trick, round off errors become
!  noticeable already when the smallest optical depth is less than the
!  square root of the machine precisions.
!
!  Operation count: 3 divide + 9 multiply + 8 add = 20 flops
!
!***********************************************************************
!
  implicit none
  integer, intent(in)              :: n
  real, intent(in) , dimension(n)  :: dtau, s
  real, intent(in)                 :: Ifw1, IrvN
  real, intent(out), dimension(n)  :: q
  real, intent(out)                :: Irv1, IfwN

  real, dimension(n)               :: a1, a2, a3, d
  integer                          :: k
  real                             :: a, dsdtau_1, dsdtau_2, d2sdtau2
!
!  Boundary condition at k=n; infalling forward intensity = Ifw1
!
  a2(1) = +1. + dtau(2) + 0.5*dtau(2)**2
  a3(1) = -1.
  dsdtau_1 = (s(2)-s(1))/dtau(2)
  dsdtau_2 = (s(3)-s(2))/dtau(3)
  d2sdtau2 = (dsdtau_2-dsdtau_1)*2./(dtau(2)+dtau(3))
  q(1)  = - dtau(2)*(s(1) - dsdtau_1 - Ifw1) + 0.5*dtau(2)**2*d2sdtau2  ! 
  a2(1) = + dtau(2) + 0.5*dtau(2)**2                                    ! store row sum
!
!  Compute matrix elements and eliminate the sub-diagonal in one sweep
!
  d(2) = -1./dtau(2)
  do k=2,n-1
    a       = 2./(dtau(k+1)+dtau(k))                                    ! 1d    1a
    d(k+1)  = -1./dtau(k+1)                                             ! 1d
    a3(k)   = a*d(k+1)                                                  !    1m
    a1(k)   = a*d(k)                                                    !    1m
    a2(k)   = 1.
    q(k)    = a1(k)*(s(k)-s(k-1)) + a3(k)*(s(k)-s(k+1))                 !    2m 3a
    a       = 1./(a2(k-1)-a3(k-1))                                      ! 1d    1a
    a1(k)   = a1(k)*a                                                   !    1m
    q(k)    = q(k) - a1(k)*q(k-1)                                       !    1m 1a
    a2(k)   = a2(k) - a1(k)*a2(k-1)                                     !    1m 1a
    a2(k-1) = a
  end do
!
!  Boundary condition at k=n; infalling reverse intensity = IrvN
!
  a2(n) = +1. + dtau(n) + 0.5*dtau(n)**2
  a1(n) = -1.
  dsdtau_1 = (s(n  )-s(n-1))/dtau(n  )
  dsdtau_2 = (s(n-1)-s(n-2))/dtau(n-1)
  d2sdtau2 = (dsdtau_1-dsdtau_2)*2./(dtau(n)+dtau(n-1))
  q(n)  = - dtau(n)*(s(n) + dsdtau_1 - IrvN) - 0.5*dtau(n)**2*d2sdtau2
  a2(n) = + dtau(n) + 0.5*dtau(n)**2
!
  a       = 1./(a2(n-1)-a3(n-1))
  a1(n)   = a1(n)*a
  q(n)    = q(n) - a1(n)*q(n-1)
  a2(n)   = a2(n) - a1(n)*a2(n-1)
  a2(n-1) = a
!
!  Backsubstitute
!
  q(n) = q(n)/a2(n)                                                     ! solution at n
  do k=n-1,1,-1
    q(k) = (q(k)-a3(k)*q(k+1))*a2(k)                                    !    2m 1a
  end do
!                                                                       ! 3d 9m 8a = 20f
!  Surface intensities
!
  Irv1 = 2.*(q(1) + s(1))-Ifw1                                          ! reverse surface intensity
  IfwN = 2.*(q(n) + s(n))-IrvN                                          ! forward surface intensity
END
