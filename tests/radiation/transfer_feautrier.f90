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
!  Operation count: 4 divide + 7 multiply + 9 add = 20 flops
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
!  Boundary condition at k=1; infalling intensity = Ifw1
!
  a2(1) = +1. + dtau(2) + 0.5*dtau(2)**2                                ! matrix element a(1,1)
  a3(1) = -1.                                                           ! matrix element a(2,1)
  dsdtau_1 = (s(2)-s(1))/dtau(2)                                        ! dS/dtau at i=1+1/2
  dsdtau_2 = (s(3)-s(2))/dtau(3)                                        ! dS/dtau at i=2+1/2
  d2sdtau2 = (dsdtau_2-dsdtau_1)*2./(dtau(2)+dtau(3))                   ! d2S/dtau2
  q(1)  = - dtau(2)*(s(1) - dsdtau_1 - Ifw1) + 0.5*dtau(2)**2*d2sdtau2  ! right hand side
!
!  Compute matrix elements in the interior
!
  d(2) = -1./dtau(2)                                      ! flop count
  do k=2,n-1
    a       = 2./(dtau(k+1)+dtau(k))                      ! 1d    1a    ! store and reuse
    d(k+1)  = -1./dtau(k+1)                               ! 1d          ! save a divide 
    a3(k)   = a*d(k+1)                                    !    1m       ! matrix element a(i-1,i)
    a1(k)   = a*d(k)                                      !    1m       ! matrix element a(i+1,i)
    a2(k)   = 1.-a1(k)-a3(k)                              !       2a    ! matrix element a(i,i)
    q(k)    = a1(k)*(s(k)-s(k-1)) + a3(k)*(s(k)-s(k+1))   !    2m 3a    ! right hand sides
  end do
!
!  Boundary condition at k=n; infalling intensity = IrvN
!
  a2(n) = +1. + dtau(n) + 0.5*dtau(n)**2                                ! matrix element a(n,n)
  a1(n) = -1.                                                           ! matrix element a(n-1,n)
  dsdtau_1 = (s(n  )-s(n-1))/dtau(n  )
  dsdtau_2 = (s(n-1)-s(n-2))/dtau(n-1)
  d2sdtau2 = (dsdtau_1-dsdtau_2)*2./(dtau(n)+dtau(n-1))
  q(n)  = - dtau(n)*(s(n) + dsdtau_1 - IrvN) - 0.5*dtau(n)**2*d2sdtau2  ! right hand side
!
!  Eliminate subdiagonal
!
  do k=2,n
    a1(k)   = a1(k)/a2(k-1)                               ! 1d          ! fraction of k-1 to subtract
    q(k)    = q(k) - a1(k)*q(k-1)                         !    1m 1a    ! subtract from row k RHS
    a2(k)   = a2(k) - a1(k)*a3(k-1)                       !    1m 1a    ! subtract from row k LHS
  end do
!
!  Backsubstitute
!
  q(n) = q(n)/a2(n)                                                     ! solution at n
  do k=n-1,1,-1
    q(k) = (q(k)-a3(k)*q(k+1))/a2(k)                      ! 1d 1m 1a    ! remaining super-diagonal
  end do
!                                                         ! 4d 9m 9a    ! 20 flops
!  Surface intensities
!
  Irv1 = 2.*(q(1) + s(1))-Ifw1                                          ! reverse surface intensity
  IfwN = 2.*(q(n) + s(n))-IrvN                                          ! forward surface intensity
END
