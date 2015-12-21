!***********************************************************************
SUBROUTINE feautrier8 (M, N, dtau, S, Ifw1, IrvN, qq, Irv1, IfwN)
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
  integer, intent(in)               :: m, n
  real, intent(in) , dimension(m,n) :: dtau, s
  real, intent(in) , dimension(m)   :: Ifw1, IrvN
  real, intent(out), dimension(m,n) :: qq
  real, intent(out), dimension(m)   :: Irv1, IfwN

  real(8), dimension(m,n)           :: a1, a2, a3, d, q
  integer                           :: j, k
  real(8)                           :: a, dsdtau_1, dsdtau_2, d2sdtau2
  logical                           :: debug=.false., first_time=.true.
  !.............................................................................
  if (debug) then
    if (first_time) open(20,file='rt8.dat',form='unformatted',status='unknown')
    first_time = .false.
    write(20) m, n
    write(20) dtau, s, Ifw1, IrvN
  end if
  !-----------------------------------------------------------------------------
  ! Boundary condition at k=1; infalling forward intensity = Ifw1
  !-----------------------------------------------------------------------------
  do j=1,m
    a2(j,1)  = + dtau(j,2) + 0.5*dtau(j,2)**2
    a3(j,1)  = -1.
    dsdtau_1 = (s(j,2)-s(j,1))/dtau(j,2)
    dsdtau_2 = (s(j,3)-s(j,2))/dtau(j,3)
    d2sdtau2 = (dsdtau_2-dsdtau_1)*2./(dtau(j,2)+dtau(j,3))
    q(j,1)   = - dtau(j,2)*(s(j,1) - dsdtau_1 - Ifw1(j)) + 0.5*dtau(j,2)**2*d2sdtau2
    d(j,2)   = -1./dtau(j,2)
  end do
  !-----------------------------------------------------------------------------
  ! Compute matrix elements and eliminate the sub-diagonal in one sweep
  !-----------------------------------------------------------------------------
  do k=2,n-1
  do j=1,m
    a        = 2./(dtau(j,k+1)+dtau(j,k))                               ! 1d    1a
    d(j,k+1) = -1./dtau(j,k+1)                                          ! 1d
    a3(j,k)  = a*d(j,k+1)                                               !    1m
    a1(j,k)  = a*d(j,k)                                                 !    1m
    a2(j,k)  = 1.
    q(j,k)   = a1(j,k)*(s(j,k)-s(j,k-1)) + a3(j,k)*(s(j,k)-s(j,k+1))    !    2m 3a
    a        = 1./(a2(j,k-1)-a3(j,k-1))                                 ! 1d    1a
    a1(j,k)  = a1(j,k)*a                                                !    1m
    q(j,k)   = q(j,k) - a1(j,k)*q(j,k-1)                                !    1m 1a
    a2(j,k)  = a2(j,k) - a1(j,k)*a2(j,k-1)                              !    1m 1a
    a2(j,k-1) = a
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Boundary condition at k=n; infalling reverse intensity = IrvN
  !-----------------------------------------------------------------------------
  do j=1,m
    a2(j,n) = dtau(j,n) + 0.5*dtau(j,n)**2
    a1(j,n) = -1.
    dsdtau_1 = (s(j,n  )-s(j,n-1))/dtau(j,n  )
    dsdtau_2 = (s(j,n-1)-s(j,n-2))/dtau(j,n-1)
    d2sdtau2 = (dsdtau_1-dsdtau_2)*2./(dtau(j,n)+dtau(j,n-1))
    q(j,n)  = - dtau(j,n)*(s(j,n) + dsdtau_1 - IrvN(j)) - 0.5*dtau(j,n)**2*d2sdtau2
!
    a       = 1./(a2(j,n-1)-a3(j,n-1))                                  ! 1d    1a
    a1(j,n) = a1(j,n)*a                                                 !    1m
    q(j,n)    = q(j,n) - a1(j,n)*q(j,n-1)                               !    1m 1a
    a2(j,n)   = a2(j,n) - a1(j,n)*a2(j,n-1)                             !    1m 1a
    a2(j,n-1) = a
    q(j,n) = q(j,n)/a2(j,n)                                             ! solution at n
  end do
  !-----------------------------------------------------------------------------
  ! Backsubstitute
  !-----------------------------------------------------------------------------
  do k=n-1,1,-1
  do j=1,m
    q(j,k) = (q(j,k)-a3(j,k)*q(j,k+1))*a2(j,k)                          !    2m 1a
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Surface intensities
  !-----------------------------------------------------------------------------
  Irv1 = 2.*(q(:,1) + s(:,1))-Ifw1                                      ! reverse surface intensity
  IfwN = 2.*(q(:,n) + s(:,n))-IrvN                                      ! forward surface intensity
  qq = q
  if (debug) write(20) qq, Irv1, IfwN
END
