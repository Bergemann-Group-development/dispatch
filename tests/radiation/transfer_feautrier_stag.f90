SUBROUTINE transfer (n,dtau,s,q,xi,doxi)
!> =====================================================================
!> Solve the transfer equation, given optical depth and source function.
!> "Steins trick" is used; storing the sum of the three elements in each
!> row of the tridiagonal equation system, instead of the diagonal ele-
!> ment.  This allows accurate solutions to be otained for arbitrarily
!> small optical depths.  Without this trick, round off errors become
!> noticeable already when the smallest optical depth is less than the
!> square root of the machine precisions.
!
!> This version works with P as a variable for k <= n1, and with
!> Q = P - S for k > n1.  P is Feautriers P and S is the source function.
!> The equation solved is  Q" = P - S for k<=n1, and Q" = Q - S" for k>n1.
!> The answer returned is Q = P - S for all k.
!
!> This choice of equations produces a high precision in Q at large optical
!> depths (where P - S would have round off), while at the same time
!> retaining the precision at small optical depths (where S" becomes very
!> large and eventually would cause round off errors in the back
!> substitution).
!
!> This is a second order version.  For simulations, with rapidly
!> varying absorption coefficients and source functions, this is to be
!> preferred over spline and Hermitean versions because it is positive
!> definite, in the sense that a positive source function is guaranteed
!> to result in a positive average intensity p.  Also, the flux
!> divergence is exactly equal to q, for the conventional definition
!> of the flux.
!
!> Operation count: 5d+5m+11a = 21 flops
!
!> Timings:
!>           Alliant: 0.58 * 21 / 1.25 = 9.7 Mfl @ 31*31*31
!
!> Update history:  
!> 12/01/2017 This routine is based on transfer_feautrier_omp from
!> stagger. 
!
!> =====================================================================
  implicit none
  integer, intent(in)              :: n
  real, intent(in), dimension(n) :: dtau, s
  real, intent(out), dimension(n)::q
  real, intent(out)::xi
  real, dimension(n) :: sp1, sp2, sp3
  real   :: ex
  real :: taum, dinv, a
  integer :: nx, ny, l, k, n1p
  logical :: doxi
  !---------------------------------------------------------------------
  ny = n
  
  
  !---------------------------------------------------------------------
  ! Calculate 1-exp(-tau(l,1,iz)) straight or by series expansion,
  ! based on max tau(l,1,iz).
  !---------------------------------------------------------------------
    taum=amax1(taum,dtau(1))
    
  if (taum.gt.0.1) then
      ex=1.-exp(-dtau(1))
  else
      ex=dtau(1)*(1.-.5*dtau(1)*(1.-.3333*dtau(1)))
  endif

  !---------------------------------------------------------------------
  ! Boundary condition at l=1.
  !---------------------------------------------------------------------
    sp2(1)=dtau(2)*(1.+.5*dtau(2))
    sp3(1)=-1.
    q  (1)=s(1)*dtau(2)*(.5*dtau(2)+ex)
  
  !---------------------------------------------------------------------
  ! Matrix elements for k>2, k<np: [3d+3s]
  !---------------------------------------------------------------------
  do k=2,ny-1
      sp2(k)=1.
      dinv=2./(dtau(k+1)+dtau(k))
      sp1(k)=-dinv/dtau(k)
      sp3(k)=-dinv/dtau(k+1)
  end do
    
  !---------------------------------------------------------------------
  ! Right hand sides, for k>2, k<n1
  !---------------------------------------------------------------------
  n1p=max0(2,ny-3)
  do k=2,n1p-1
      q(k)=s(k)
  end do
  
  !---------------------------------------------------------------------
  ! RHS for k=n1p.  The equation is still in p(k), but q(k+1) is now
  ! q(k+1)=p(k+1)-s(k+1), so to get p(k+1) we have to add s(k+1) to the
  ! LHS, that is subtract it from the RHS.
  !---------------------------------------------------------------------
  if (n1p.le.ny) then
      q(n1p)=s(n1p)-sp3(n1p)*s(n1p+1)
  endif
  !---------------------------------------------------------------------
  ! RHS for k=n1p+1.  The equation is now in Q: Q" = Q - S", but q(k-1)
  ! is now p(k-1)=q(k-1)+s(k-1), so we must add sp1*s(l,k-1,iz) to the RHS,
  ! relative to the case below.
  !---------------------------------------------------------------------
  if (n1p+1.le.ny) then
      q(n1p+1)=sp1(n1p+1)* s(n1p+1) &
                +sp3(n1p+1)*(s(n1p+1)-s(n1p+2))
  endif
  
  !---------------------------------------------------------------------
  ! k>n1p+1,<np [2m+2a]
  !---------------------------------------------------------------------
  do k=n1p+2,ny-1
      q(k)=sp1(k)*(s(k)-s(k-1)) &
            +sp3(k)*(s(k)-s(k+1))
  end do
  
  !---------------------------------------------------------------------
  ! k=np
  !---------------------------------------------------------------------
      sp2(ny)=1.
      q(ny)=0.0
    
  !---------------------------------------------------------------------
  !  Eliminate subdiagonal, save factors in sp1 [1d+2m+4a]
  !---------------------------------------------------------------------
  do k=1,ny-2
      a=-sp1(k+1)/(sp2(k)-sp3(k))
      q(k+1)=q(k+1)+a*q(k)
      sp2(k+1)=sp2(k+1)+a*sp2(k)
      sp2(k)=sp2(k)-sp3(k)
  end do
  
  sp2(ny-1)=sp2(ny-1)-sp3(ny-1)
  
  
  !---------------------------------------------------------------------
  ! Backsubstitute [1d+1m+1a]
  !---------------------------------------------------------------------
  do k=ny-1,1,-1
      q(k)=(q(k)-sp3(k)*q(k+1))/sp2(k)
  end do
  
  !---------------------------------------------------------------------
  ! Subtract S in the region where the equation is in P
  !---------------------------------------------------------------------
  do k=1,n1p
      q(k)=q(k)-s(k)
  end do
  
  !---------------------------------------------------------------------
  ! Surface intensity.
  !---------------------------------------------------------------------
  if (doxi) then
      xi=2.*(1.-ex)*(q(1)+s(1))+s(1)*ex**2 
  

END SUBROUTINE transfer
