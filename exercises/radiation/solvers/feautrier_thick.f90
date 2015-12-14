!***********************************************************************
SUBROUTINE transfer (n, dtau, s, Ifw0, Ire1, q, Ire0, Ifw1)
!
!  Solve the transfer equation using an optimized Feautrier method, given 
!  optical depth increments dtau, source function S, and infalling intensity
!  Ifw0.
!
!  The reverse infalling intensity Ire1 at the other boundary is not used 
!  in this version, which assumes Q=0 at the right-most bdry.
!
!  This version works with Q = P - S, where P is Feautriers P (the average 
!  of in- and out-going specific intensity) and S is the source function.
!  The equation solved is Q" = Q - S", with boundary conditions 
!  Q' - Q = S - S' - I0 at the first point and Q = 0 at the last point.
!
!  It is straightforward to change the boundary condition in the last
!  point to one with a given infalling intensity and a computed outgoing
!  intensity.
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
  real, intent(in)                 :: Ifw0, Ire1
  real, intent(out), dimension(n)  :: q
  real, intent(out)                :: Ire0, Ifw1

  real, dimension(n)               :: a2, a3, d
  integer                          :: k
  real                             :: a, a1
!
!  Boundary condition at k=1; infalling intensity = i0
!
  a2(1) = -1. - 0.5*dtau(2)
  a3(1) = +1. - 0.5*dtau(2)
  q(1)  = -a2(1)*s(1) - a3(1)*s(2) + dtau(2)*Ifw0
  a2(1) = -dtau(2)
!
!  Boundary condition at k=n; intensity = source function
!
  a2(n)=1.
  q(n)=0.0
!
!  Compute matrix elements and eliminate the sub-diagonal in one sweep
!
  d(2) = -1./dtau(2)
  do k=2,n-1
    a       = 2./(dtau(k+1)+dtau(k))                                    ! 1d    1a
    d(k+1)  = -1./dtau(k+1)                                             ! 1d
    a3(k)   = a*d(k+1)                                                  !    1m
    a1      = a*d(k)                                                    !    1m
    q(k)    = a1*(s(k)-s(k-1)) + a3(k)*(s(k)-s(k+1))                    !    2m 3a
    a       = 1./(a2(k-1)-a3(k-1))                                      ! 1d    1a
    a1      = a1*a                                                      !    1m
    q(k)    = q(k) - a1*q(k-1)                                          !    1m 1a
    a2(k)   = 1. - a1*a2(k-1)                                           !    1m 1a
    a2(k-1) = a
  end do
!
!  Backsubstitute
!
  do k=n-1,1,-1
    q(k) = (q(k)-a3(k)*q(k+1))*a2(k)                                    !    2m 1a
  end do
!                                                                       ! 3d 9m 8a = 20f
!  Surface intensity.
!
  Ire0 = 2.*(q(1) + s(1))-Ifw0
  Ifw1 = s(n)
END
