!> =============================================================================
!> Module with RT solvers -- for now just an integral solver
!> =============================================================================
MODULE rt_solver_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: rt_solver_t
    procedure(feautrier), pointer:: transfer => null()
    procedure(integral_forward), pointer:: forward  => null()
    procedure(integral_forward), pointer:: reverse  => null()
  contains
    procedure:: init
    procedure:: integral_forward
    procedure:: integral_reverse
    procedure:: feautrier
    procedure:: transfer2
    procedure:: feautrier_stag
  end type
  type (rt_solver_t), public:: rt_solver
CONTAINS

!> =============================================================================
!> Choose RT solver method
!> =============================================================================
SUBROUTINE init (self)
  class (rt_solver_t):: self
  self%forward   => integral_forward
  self%reverse   => integral_reverse
  self%transfer  => feautrier
END SUBROUTINE init

!> =============================================================================
!> Integral solution method
!> =============================================================================
SUBROUTINE integral_forward (self, dtau, ex0, S, I_in, q, I_out, one_p)
  implicit none
  class (rt_solver_t):: self
  real, dimension(:,:,:), pointer :: dtau, ex0, S, q
  real, dimension(:,:)   :: I_in, I_out
  logical :: one_p
!-----------------------------------------------------------------------
!
!  Input: N, dtau(2:n), S(1:n), Ifw(1), Irv(n)
! Output: Q(1:n-1), qfw(n-1)
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
!-----------------------------------------------------------------------
  real(kind=8), dimension(size(s,1),size(s,2),size(s,3)) :: ex1,ex2, ds1, d1s, d2s
  real(kind=8), dimension(size(s,1),size(s,2))           :: ctau
  integer i, n
!-----------------------------------------------------------------------
!

  n = size(s,3)
  do i=2,n
    ex1(:,:,i) = 1.-ex0(:,:,i)
    ex2(:,:,i) = ex1(:,:,i) - dtau(:,:,i)*ex0(:,:,i)
    ds1(:,:,i) = (s(:,:,i) - s(:,:,i-1)) / dtau(:,:,i)                  ! half-point centered dS/dtau
  end do 
  
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    ctau = 1. / (dtau(:,:,i+1) + dtau(:,:,i))
    d1s(:,:,i) =(ds1(:,:,i+1) * dtau(:,:,i) &
               + ds1(:,:,i  ) * dtau(:,:,i+1)) * ctau                   ! dS/dtau, 2nd. order accurate
    d2s(:,:,i)= (ds1(:,:,i+1) - ds1 (:,:,i  )) * 2. * ctau              ! d2S/dtau2, 2nd. order accurate
  end do
  d2s(:,:,1) = d2s(:,:,2)                                               ! assume constant d2S/dtau2
  d1s(:,:,1) = d1s(:,:,2) - dtau(:,:,2)*d2s(:,:,2)                      ! extrapolate to 1
  d2s(:,:,n) = d2s(:,:,n-1)
  d1s(:,:,n) = d1s(:,:,n-1) - dtau(:,:,n)*d2s(:,:,n)

  q(:,:,1) = I_in - s(:,:,1)                                            ! initial forward Q
  do i=2,n
    q(:,:,i) = q(:,:,i-1) *ex0(:,:,i) &
             - d1s(:,:,i) *ex1(:,:,i) &
             + d2s(:,:,i) *ex2(:,:,i)                                   ! formal solution forward
  end do
  if (one_p) then
    I_out = q(:,:,n-1) + s(:,:,n-1)                                     ! forward boundary intensity if a single/two point ray is formulated
  else
    I_out = q(:,:,n) + s(:,:,n)                                         ! forward boundary intensity
  end if

END SUBROUTINE integral_forward                         ! 2d 13m 11a 1e ! 26f + 1e

SUBROUTINE integral_reverse (self, dtau, ex0, S, I_in, q, I_out, one_p)
  implicit none
  class (rt_solver_t):: self
  real, dimension(:,:,:), pointer :: dtau, ex0, S, q
  real, dimension(:,:)   :: I_in, I_out
  logical :: one_p
!-----------------------------------------------------------------------
!
!  Input: N, dtau(2:n), S(1:n), Ifw(1), Irv(n)
! Output: Q(1:n-1), qfw(n-1)
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
!-----------------------------------------------------------------------
  real(kind=8), dimension(size(s,1),size(s,2),size(s,3)) :: ex1,ex2, ds1, d1s, d2s, st
  real(kind=8), dimension(size(s,1),size(s,2))           :: ctau
  integer i, n
!-----------------------------------------------------------------------
!

  n = size(s,3)
  do i=2,n
    ex1(:,:,i) = 1.-ex0(:,:,i)
    ex2(:,:,i) = ex1(:,:,i) - dtau(:,:,i)*ex0(:,:,i)
    ds1(:,:,i) = (s(:,:,i) - s(:,:,i-1)) / dtau(:,:,i)                  ! half-point centered dS/dtau
    st(:,:,i) = s(:,:,i) - s(:,:,i-1)
  end do 
  
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    ctau = 1. / (dtau(:,:,i+1) + dtau(:,:,i))
    d1s(:,:,i) =(ds1(:,:,i+1) * dtau(:,:,i) &
               + ds1(:,:,i  ) * dtau(:,:,i+1)) * ctau                   ! dS/dtau, 2nd. order accurate
    d2s(:,:,i)= (ds1(:,:,i+1) - ds1 (:,:,i  )) * 2. * ctau              ! d2S/dtau2, 2nd. order accurate
  end do
  d2s(:,:,1) = d2s(:,:,2)                                               ! assume constant d2S/dtau2
  d1s(:,:,1) = d1s(:,:,2) - dtau(:,:,2)*d2s(:,:,2)                      ! extrapolate to 1
  d2s(:,:,n) = d2s(:,:,n-1)
  d1s(:,:,n) = d1s(:,:,n-1) - dtau(:,:,n)*d2s(:,:,n)
  
  
  

  q(:,:,1) = I_in - s(:,:,n)                                            ! initial reverse Q
  
  do i=n-1,1,-1
    q(:,:,i) = q(:,:,i+1) *ex0(:,:,i+1) &
             + d1s(:,:,i) *ex1(:,:,i+1) &
             + d2s(:,:,i) *ex2(:,:,i+1)                                  ! formal solution reverse
  end do
    
  if (one_p) then
    I_out = q(:,:,2) + s(:,:,2)                                         ! reverse boundary intensity if a single/two point ray is formulated
  else
    I_out = q(:,:,1) + s(:,:,1)                                         ! reverse boundary intensity
  end if

END SUBROUTINE integral_reverse                         ! 2d 13m 11a 1e ! 26f + 1e

SUBROUTINE feautrier (self, dtau, S, I_in1, I_in2, Q, I_out1, I_out2)
  implicit none
  class (rt_solver_t):: self
  real, dimension(:,:,:), pointer :: dtau, S, q
  real, dimension(:,:) :: I_in1, I_in2, I_out1, I_out2

  real, dimension(size(s,1),size(s,2),size(s,3)) :: a1, a2, a3, d
  real, dimension(size(s,2)) :: a, dsdtau_1, dsdtau_2, d2sdtau2
  integer :: k, m, n, j
  
!
  m = size(s,1)
  n = size(s,3)
  
  do j=1,m
    a2(j,:,1)  = +1. + dtau(j,:,2) + 0.5*dtau(j,:,2)**2
    a3(j,:,1)  = -1.
    dsdtau_1 = (s(j,:,2)-s(j,:,1))/dtau(j,:,2)
    dsdtau_2 = (s(j,:,3)-s(j,:,2))/dtau(j,:,3)
    d2sdtau2 = (dsdtau_2-dsdtau_1)*2./(dtau(j,:,2)+dtau(j,:,3))
    q(j,:,1)   = I_in1(j,:)
    a2(j,:,1)  = + dtau(j,:,2) + 0.5*dtau(j,:,2)**2
    d(j,:,2)   = -1./dtau(j,:,2)
  end do
!
  !
!  Compute matrix elements and eliminate the sub-diagonal in one sweep
!
  do k=2,n-1
  do j=1,m
    a(:)        = 2./(dtau(j,:,k+1)+dtau(j,:,k))                               ! 1d    1a
    d(j,:,k+1) = -1./dtau(j,:,k+1)                                          ! 1d
    a3(j,:,k)  = a*d(j,:,k+1)                                               !    1m
    a1(j,:,k)  = a*d(j,:,k)                                                 !    1m
    a2(j,:,k)  = 1.
    q(j,:,k)   = a1(j,:,k)*(s(j,:,k)-s(j,:,k-1)) + a3(j,:,k)*(s(j,:,k)-s(j,:,k+1))    !    2m 3a
    a        = 1./(a2(j,:,k-1)-a3(j,:,k-1))                                 ! 1d    1a
    a1(j,:,k)  = a1(j,:,k)*a                                                !    1m
    q(j,:,k)   = q(j,:,k) - a1(j,:,k)*q(j,:,k-1)                                !    1m 1a
    a2(j,:,k)  = a2(j,:,k) - a1(j,:,k)*a2(j,:,k-1)                              !    1m 1a
    a2(j,:,k-1) = a
  end do
  end do
  
  !
!  Boundary condition at k=n; infalling reverse intensity = IrvN
!
  do j=1,m
   a2(j,:,n) = +1. + dtau(j,:,n) + 0.5*dtau(j,:,n)**2
   a1(j,:,n) = -1.
   dsdtau_1 = (s(j,:,n  )-s(j,:,n-1))/dtau(j,:,n  )
   dsdtau_2 = (s(j,:,n-1)-s(j,:,n-2))/dtau(j,:,n-1)
   d2sdtau2 = (dsdtau_1-dsdtau_2)*2./(dtau(j,:,n)+dtau(j,:,n-1))
   q(j,:,n)  = I_In2(j,:) - 0.5*dtau(j,:,n)**2*d2sdtau2
   a2(j,:,n) = + dtau(j,:,n) + 0.5*dtau(j,:,n)**2
!
   a       = 1./(a2(j,:,n-1)-a3(j,:,n-1))                                   ! 1d    1a
   a1(j,:,n) = a1(j,:,n)*a                                                  !    1m
   q(j,:,n)    = q(j,:,n) - a1(j,:,n)*q(j,:,n-1)                                !    1m 1a
   a2(j,:,n)   = a2(j,:,n) - a1(j,:,n)*a2(j,:,n-1)                              !    1m 1a
   a2(j,:,n-1) = a
   q(j,:,n) = q(j,:,n)/a2(j,:,n)                                              ! solution at n
  end do
!
!  Backsubstitute
!
  do k=n-1,1,-1
  do j=1,m
    q(j,:,k) = (q(j,:,k)-a3(j,:,k)*q(j,:,k+1))*a2(j,:,k)                          !    2m 1a
  end do
  end do
!                                                                       ! 3d 9m 8a = 20f
!  Surface intensities
!
  I_out2 = 2.*q(:,:,1) + s(:,:,1)                                           ! reverse surface intensity
  I_out1 = 2.*q(:,:,n) + s(:,:,n)                                           ! forward surface intensity
  


END SUBROUTINE feautrier

SUBROUTINE transfer2 (self,M, N, dtau, S, Ifw1, IrvN, Q, Irv1, IfwN)
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
class (rt_solver_t):: self

  integer, intent(in)               :: m, n
  real, intent(in) , dimension(m,n) :: dtau, s
  real, intent(in) , dimension(m)   :: Ifw1, IrvN
  real, intent(out), dimension(m,n) :: q
  real, intent(out), dimension(m)   :: Irv1, IfwN

  real, dimension(m,n)              :: a1, a2, a3, d
  integer                           :: j, k
  real                              :: a, dsdtau_1, dsdtau_2, d2sdtau2
!
!  Boundary condition at k=1; infalling forward intensity = Ifw1
!
  do j=1,m
    a2(j,1)  = +1. + dtau(j,2) + 0.5*dtau(j,2)**2
    a3(j,1)  = -1.
    dsdtau_1 = (s(j,2)-s(j,1))/dtau(j,2)
    dsdtau_2 = (s(j,3)-s(j,2))/dtau(j,3)
    d2sdtau2 = (dsdtau_2-dsdtau_1)*2./(dtau(j,2)+dtau(j,3))
    q(j,1)   = - dtau(j,2)*(s(j,1) - dsdtau_1 - Ifw1(j)) + 0.5*dtau(j,2)**2*d2sdtau2
    a2(j,1)  = + dtau(j,2) + 0.5*dtau(j,2)**2
    d(j,2)   = -1./dtau(j,2)
  end do
!
!  Compute matrix elements and eliminate the sub-diagonal in one sweep
!
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
!
!  Boundary condition at k=n; infalling reverse intensity = IrvN
!
  do j=1,m
   a2(j,n) = +1. + dtau(j,n) + 0.5*dtau(j,n)**2
   a1(j,n) = -1.
   dsdtau_1 = (s(j,n  )-s(j,n-1))/dtau(j,n  )
   dsdtau_2 = (s(j,n-1)-s(j,n-2))/dtau(j,n-1)
   d2sdtau2 = (dsdtau_1-dsdtau_2)*2./(dtau(j,n)+dtau(j,n-1))
   q(j,n)  = - dtau(j,n)*(s(j,n) + dsdtau_1 - IrvN(j)) - 0.5*dtau(j,n)**2*d2sdtau2
   a2(j,n) = + dtau(j,n) + 0.5*dtau(j,n)**2
!
   a       = 1./(a2(j,n-1)-a3(j,n-1))                                   ! 1d    1a
   a1(j,n) = a1(j,n)*a                                                  !    1m
   q(j,n)    = q(j,n) - a1(j,n)*q(j,n-1)                                !    1m 1a
   a2(j,n)   = a2(j,n) - a1(j,n)*a2(j,n-1)                              !    1m 1a
   a2(j,n-1) = a
   q(j,n) = q(j,n)/a2(j,n)                                              ! solution at n
  end do
!
!  Backsubstitute
!
  do k=n-1,1,-1
  do j=1,m
    q(j,k) = (q(j,k)-a3(j,k)*q(j,k+1))*a2(j,k)                          !    2m 1a
  end do
  end do
!                                                                       ! 3d 9m 8a = 20f
!  Surface intensities
!
  Irv1 = 2.*(q(:,1) + s(:,1))-Ifw1                                      ! reverse surface intensity
  IfwN = 2.*(q(:,n) + s(:,n))-IrvN                                      ! forward surface intensity
END SUBROUTINE transfer2

SUBROUTINE feautrier_stag (self,dtau,s,q,xi,doxi)
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
  class (rt_solver_t):: self
  real, dimension(:,:), pointer :: dtau, s, q
  real, dimension(:,:), pointer :: sp1, sp2, sp3
  real, dimension(:), pointer   :: ex, xi
  real :: taum, dinv, a
  integer :: nx, ny, l, k, n1p
  logical :: doxi
  !---------------------------------------------------------------------
  nx = size(s,1)
  ny = size(s,2)
  
  !---------------------------------------------------------------------
  ! Allocate the necessary arrays
  !---------------------------------------------------------------------
  allocate(sp1(nx,ny),sp2(nx,ny), sp3(nx,ny), ex(nx))
  if (doxi) allocate(xi(nx))
  
  !---------------------------------------------------------------------
  ! Calculate 1-exp(-tau(l,1,iz)) straight or by series expansion,
  ! based on max tau(l,1,iz).
  !---------------------------------------------------------------------
  taum=dtau(1,1)
  do l=2,nx
    taum=amax1(taum,dtau(l,1))
  end do
  if (taum.gt.0.1) then
    do l=1,nx
      ex(l)=1.-exp(-dtau(l,1))
    end do
  else
    do l=1,nx
      ex(l)=dtau(l,1)*(1.-.5*dtau(l,1)*(1.-.3333*dtau(l,1)))
    end do
  endif

  !---------------------------------------------------------------------
  ! Boundary condition at l=1.
  !---------------------------------------------------------------------
  do l=1,nx
    sp2(l,1)=dtau(l,2)*(1.+.5*dtau(l,2))
    sp3(l,1)=-1.
    q  (l,1)=s(l,1)*dtau(l,2)*(.5*dtau(l,2)+ex(l))
  end do
  
  !---------------------------------------------------------------------
  ! Matrix elements for k>2, k<np: [3d+3s]
  !---------------------------------------------------------------------
  do k=2,ny-1
    do l=1,nx
      sp2(l,k)=1.
      dinv=2./(dtau(l,k+1)+dtau(l,k))
      sp1(l,k)=-dinv/dtau(l,k)
      sp3(l,k)=-dinv/dtau(l,k+1)
    end do
  end do
    
  !---------------------------------------------------------------------
  ! Right hand sides, for k>2, k<n1
  !---------------------------------------------------------------------
  n1p=max0(2,ny-3)
  do k=2,n1p-1
    do l=1,nx
      q(l,k)=s(l,k)
    end do
  end do
  
  !---------------------------------------------------------------------
  ! RHS for k=n1p.  The equation is still in p(k), but q(k+1) is now
  ! q(k+1)=p(k+1)-s(k+1), so to get p(k+1) we have to add s(k+1) to the
  ! LHS, that is subtract it from the RHS.
  !---------------------------------------------------------------------
  if (n1p.le.ny) then
    do l=1,nx
      q(l,n1p)=s(l,n1p)-sp3(l,n1p)*s(l,n1p+1)
    end do
  endif
  !---------------------------------------------------------------------
  ! RHS for k=n1p+1.  The equation is now in Q: Q" = Q - S", but q(k-1)
  ! is now p(k-1)=q(k-1)+s(k-1), so we must add sp1*s(l,k-1,iz) to the RHS,
  ! relative to the case below.
  !---------------------------------------------------------------------
  if (n1p+1.le.ny) then
    do l=1,nx
      q(l,n1p+1)=sp1(l,n1p+1)* s(l,n1p+1) &
                +sp3(l,n1p+1)*(s(l,n1p+1)-s(l,n1p+2))
    end do
  endif
  
  !---------------------------------------------------------------------
  ! k>n1p+1,<np [2m+2a]
  !---------------------------------------------------------------------
  do k=n1p+2,ny-1
    do l=1,nx
      q(l,k)=sp1(l,k)*(s(l,k)-s(l,k-1)) &
            +sp3(l,k)*(s(l,k)-s(l,k+1))
    end do
  end do
  
  !---------------------------------------------------------------------
  ! k=np
  !---------------------------------------------------------------------
    do l=1,nx
      sp2(l,ny)=1.
      q(l,ny)=0.0
    end do
    
  !---------------------------------------------------------------------
  !  Eliminate subdiagonal, save factors in sp1 [1d+2m+4a]
  !---------------------------------------------------------------------
  do k=1,ny-2
    do l=1,nx
      a=-sp1(l,k+1)/(sp2(l,k)-sp3(l,k))
      q(l,k+1)=q(l,k+1)+a*q(l,k)
      sp2(l,k+1)=sp2(l,k+1)+a*sp2(l,k)
      sp2(l,k)=sp2(l,k)-sp3(l,k)
    end do
  end do
  
  do l=1,nx
    sp2(l,ny-1)=sp2(l,ny-1)-sp3(l,ny-1)
  end do
  
  !---------------------------------------------------------------------
  ! Backsubstitute [1d+1m+1a]
  !---------------------------------------------------------------------
  do k=ny-1,1,-1
    do l=1,nx
      q(l,k)=(q(l,k)-sp3(l,k)*q(l,k+1))/sp2(l,k)
    end do
  end do
  
  !---------------------------------------------------------------------
  ! Subtract S in the region where the equation is in P
  !---------------------------------------------------------------------
  do k=1,n1p
    do l=1,nx
      q(l,k)=q(l,k)-s(l,k)
    end do
  end do
  
  !---------------------------------------------------------------------
  ! Surface intensity.
  !---------------------------------------------------------------------
  if (doxi) then
    do l=1,nx
      xi(l)=2.*(1.-ex(l))*(q(l,1)+s(l,1))+s(l,1)*ex(l)**2 !??? what?
    end do
    deallocate(xi)
  end if
  
  deallocate(sp1,sp2, sp3, ex)

END SUBROUTINE feautrier_stag

END MODULE rt_solver_mod
