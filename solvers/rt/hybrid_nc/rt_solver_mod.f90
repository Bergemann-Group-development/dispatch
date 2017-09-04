!> =============================================================================
!> Module with RT solvers -- for now just an integral solver
!> =============================================================================
MODULE rt_solver_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  implicit none
  private
  type, public:: rt_solver_t
  contains
    procedure:: init
    procedure:: integral_forward
    procedure:: integral_reverse
    procedure:: integral_double
    procedure:: feautrier
  end type
  type (rt_solver_t), public:: rt_solver
CONTAINS

!> =====================================================================
!> Choose RT solver method
!> =====================================================================
SUBROUTINE init (self)
  class (rt_solver_t):: self
END SUBROUTINE init


!> =====================================================================
!> Integral solution method in forward direction
!>  Input: N, dtau(2:n), S(1:n), qfw(1)
!> Output: Q(1:n-1), qfw(n-1)
!>
!> The input is the source functions and optical depths at all points,
!> the net forward intensity (I-S) at the 1st point. As a result, we get
!> the net heating rates in forward direction everywhere, and also the
!> outgoing net forward intensity at the last point.
!>
!> The equation solved is Q' = Q - dS/dtau in the forward direction.
!> S is assumed to be quadratic, so dS/dtau is assumed to be linear in tau.
!> The analytic solution may be written as a linear combination of the
!> formal solution for constant and linear source functions.
!>
!> Two additional points for source function are necessary to avoid the
!> extrapolation of d1s and d2s.
!> =====================================================================
SUBROUTINE integral_forward (self, dtau, ex0, S, q_in, q)
  implicit none
  class (rt_solver_t):: self
  real(kind=KindRTVar), dimension(:,:,:), pointer :: dtau, ex0, S, q
  real(kind=KindRTVar), dimension(:,:)   :: q_in!, q_out
!-----------------------------------------------------------------------
  real(kind=KindRTVar), dimension(size(s,1),size(s,2),size(s,3)) :: ex1,ex2, ds1, d1s, d2s
  real(kind=KindRTVar) :: ctau
  integer i, n, m, w, j, k
!-----------------------------------------------------------------------
!
  m = size(s,1)
  w = size(s,2)
  n = size(s,3)
  do i=2,n
    do j = 1, w
      do k = 1, m
    ex1(k,j,i) = 1.-ex0(k,j,i)
    ex2(k,j,i) = ex1(k,j,i) - dtau(k,j,i)*ex0(k,j,i)
    ds1(k,j,i) = (s(k,j,i) - s(k,j,i-1)) / dtau(k,j,i)                  ! half-point centered dS/dtau
      end do
    end do
  end do 
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    do j=1,w
      do k=1,m
        ctau = 1. / (dtau(k,j,i+1) + dtau(k,j,i))
        d1s(k,j,i) =(ds1(k,j,i+1) * dtau(k,j,i) &
                   + ds1(k,j,i  ) * dtau(k,j,i+1)) * ctau               ! dS/dtau, 2nd. order accurate
        d2s(k,j,i)= (ds1(k,j,i+1) - ds1 (k,j,i  )) * 2. * ctau          ! d2S/dtau2, 2nd. order accurate
      end do
    end do  
  end do
!
  do j = 1, w
    do k = 1, m
      q(k,j,1) = q_in(k,j)                                              ! initial forward Q
    end do
  end do
!
  do i=2,n-1
    do j = 1, w
      do k = 1, m
        q(k,j,i) = q(k,j,i-1) *ex0(k,j,i) &
                 - d1s(k,j,i) *ex1(k,j,i) &
                 + d2s(k,j,i) *ex2(k,j,i)                               ! formal solution forward
      end do
    end do  
  end do
!
!  do j = 1, w
!    do k = 1, m
!      q_out(k,j) = q(k,j,n-1)                                           ! forward boundary intensity
!    end do
!  end do  
END SUBROUTINE integral_forward                                         ! 2d 13m 11a 1e ! 26f + 1e

!> =============================================================================
!> Integral solution method in forward direction
!>  Input: N, dtau(2:n), S(1:n), qrv(n)
!> Output: Q(1:n-1), qrv(2)
!>
!> The input is the source functions and optical depths at all points,
!> the net forward intensity (I-S) at the last point. As a result, we get
!> the net heating rates in reverse direction everywhere, and also the
!> outgoing net reverse intensity at the first point.
!>
!> The equation solved is Q' = -Q + dS/dtau in the reverse direction.
!> S is assumed to be quadratic, so dS/dtau is assumed to be linear in tau.
!> The analytic solution may be written as a linear combination of the
!> formal solution for constant and linear source functions.
!>
!> Two additional points for source function are necessary to avoid the
!> extrapolation of d1s and d2s.
!> =============================================================================
SUBROUTINE integral_reverse (self, dtau, ex0, S, q_in, q)
  implicit none
  class (rt_solver_t):: self
  real(kind=KindRTVar), dimension(:,:,:), pointer :: dtau, ex0, S, q
  real(kind=KindRTVar), dimension(:,:)   :: q_in!, q_out
!-----------------------------------------------------------------------
!
!  Input: N, dtau(2:n), S(1:n), Ifw(1), Irv(n)
! Output: Q(1:n-1), qrv(n-1)
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
  real(kind=KindRTVar), dimension(size(s,1),size(s,2),size(s,3)) :: ex1,ex2, ds1, d1s, d2s, st
  real(kind=KindRTVar) :: ctau
  integer i, n, m, w, j, k
!-----------------------------------------------------------------------
!
  m = size(s,1)
  w = size(s,2)
  n = size(s,3)
  do i=2,n
    do j = 1, w
      do k = 1, m
    ex1(k,j,i) = 1.-ex0(k,j,i)
    ex2(k,j,i) = ex1(k,j,i) - dtau(k,j,i)*ex0(k,j,i)
    ds1(k,j,i) = (s(k,j,i) - s(k,j,i-1)) / dtau(k,j,i)                  ! half-point centered dS/dtau
      end do
    end do
  end do 
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    do j=1,w
      do k=1,m
        ctau = 1. / (dtau(k,j,i+1) + dtau(k,j,i))
        d1s(k,j,i) =(ds1(k,j,i+1) * dtau(k,j,i) &
                   + ds1(k,j,i  ) * dtau(k,j,i+1)) * ctau               ! dS/dtau, 2nd. order accurate
        d2s(k,j,i)= (ds1(k,j,i+1) - ds1 (k,j,i  )) * 2. * ctau          ! d2S/dtau2, 2nd. order accurate
      end do
    end do  
  end do
!
  do j = 1, w
    do k = 1, m
      q(k,j,n) = q_in(k,j)                                              ! initial reverse Q
    end do
  end do
!
  do i=n-1,2,-1
    do j = 1, w
      do k = 1, m
        q(k,j,i) = q(k,j,i+1) *ex0(k,j,i+1) &
                 + d1s(k,j,i) *ex1(k,j,i+1) &
                 + d2s(k,j,i) *ex2(k,j,i+1)                             ! formal solution reverse
      end do
    end do  
  end do
!
!  do j = 1, w
!    do k = 1, m
!      q_out(k,j) = q(k,j,2)                                             ! reverse boundary intensity
!    end do
!  end do  
END SUBROUTINE integral_reverse                                         ! 2d 13m 11a 1e ! 26f + 1e

!> =====================================================================
!> Integral solution method in both directions
!>  Input: N, dtau(2:n), S(1:n), q_in_fw(1), q_in_re(n)
!> Output: Q(1:n-1), q_fw(2:n), q_re(1:n-1)
!>
!> The input is the source functions and optical depths at all points,
!> the net forward intensity (I-S) at the 1st point and the net reverse
!> intensity (I-S) at the Nth point . As a result, we get
!> the net heating rates in forward and reverse directions everywhere
!> and also the outgoing net forward intensity at the first/last points.
!>
!> The equation solved is Q' = Q - dS/dtau in the forward direction.
!> S is assumed to be quadratic, so dS/dtau is assumed to be linear in tau.
!> The analytic solution may be written as a linear combination of the
!> formal solution for constant and linear source functions.
!>
!> Two additional points for source function are necessary to avoid the
!> extrapolation of d1s and d2s.
!> =====================================================================
SUBROUTINE integral_double (self, dtau, ex0, S, q_in_fw, q_in_re, q_fw, q_re)
  implicit none
  class (rt_solver_t):: self
  real(kind=KindRTVar), dimension(:,:,:), pointer :: dtau, ex0, S, q_fw, q_re
  real(kind=KindRTVar), dimension(:,:)   :: q_in_fw, q_in_re
!-----------------------------------------------------------------------
  real(kind=KindRTVar), dimension(size(s,1),size(s,2),size(s,3)) :: ex1,ex2, ds1, d1s, d2s
  real(kind=KindRTVar) :: ctau
  integer i, n, m, w, j, k
!-----------------------------------------------------------------------
!
  m = size(s,1)
  w = size(s,2)
  n = size(s,3)
  do i=2,n
    do j = 1, w
      do k = 1, m
    ex1(k,j,i) = 1.-ex0(k,j,i)
    ex2(k,j,i) = ex1(k,j,i) - dtau(k,j,i)*ex0(k,j,i)
    ds1(k,j,i) = (s(k,j,i) - s(k,j,i-1)) / dtau(k,j,i)                  ! half-point centered dS/dtau
      end do
    end do
  end do 
!
!  dS/dtau and d2S/dtau2 (3a+5m+1d)
!
  do i=2,n-1
    do j=1,w
      do k=1,m
        ctau = 1. / (dtau(k,j,i+1) + dtau(k,j,i))
        d1s(k,j,i) =(ds1(k,j,i+1) * dtau(k,j,i) &
                   + ds1(k,j,i  ) * dtau(k,j,i+1)) * ctau               ! dS/dtau, 2nd. order accurate
        print *, 'ds1+1, ds1,ctau', ds1(k,j,i+1), ds1 (k,j,i  ), ctau
        d2s(k,j,i)= (ds1(k,j,i+1) - ds1 (k,j,i  )) * 2d0 * ctau          ! d2S/dtau2, 2nd. order accurate
      end do
    end do  
  end do
!
  do j = 1, w
    do k = 1, m
      q_fw(k,j,1) = q_in_fw(k,j)                                        ! initial forward Q
      q_re(k,j,n) = q_in_re(k,j)                                        ! initial reverse Q
    end do
  end do
!
  do i=2,n-1
    do j = 1, w
      do k = 1, m
        q_fw(k,j,i) = q_fw(k,j,i-1) *ex0(k,j,i) &
                       - d1s(k,j,i) *ex1(k,j,i) &
                       + d2s(k,j,i) *ex2(k,j,i)                         ! formal solution forward
      end do
    end do  
  end do
!
  do i=n-1,2,-1
    do j = 1, w
      do k = 1, m
        q_re(k,j,i) = q_re(k,j,i+1) *ex0(k,j,i+1) &
                       + d1s(k,j,i) *ex1(k,j,i+1) &
                       + d2s(k,j,i) *ex2(k,j,i+1)                       ! formal solution reverse
      end do
    end do  
  end do
END SUBROUTINE integral_double                                          ! 2d 13m 11a 1e ! 26f + 1e

SUBROUTINE feautrier (self,dtau,s,q,qfw,qre,doxi)
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
  real(kind=KindRTVar), dimension(:,:,:), pointer :: dtau, s, q
  real(kind=KindRTVar), dimension(:,:), pointer :: qfw, qre
  logical :: doxi
  !.....................................................................
  real(kind=KindRTVar), dimension(:,:,:), pointer :: sp1, sp2, sp3
  real(kind=KindRTVar), dimension(:,:), pointer   :: ex
  real :: taum, dinv, a
  integer :: nx, ny, nz, l, k, j, n1p
  !---------------------------------------------------------------------
  
  nx = size(s,1)
  ny = size(s,2)
  nz = size(s,3)
  
  !---------------------------------------------------------------------
  ! Allocate the necessary arrays
  !---------------------------------------------------------------------
  allocate(sp1(nx,ny,nz),sp2(nx,ny,nz), sp3(nx,ny,nz), ex(nx,ny))
  
  !---------------------------------------------------------------------
  ! Calculate 1-exp(-tau(l,k,1)) straight or by series expansion,
  ! based on max tau(l,k,1).
  !---------------------------------------------------------------------
  taum=dtau(1,1,1)
  do k=1,ny
    do l=1,nx
      taum=amax1(taum,dtau(l,k,1))
    end do
  end do
  if (taum.gt.0.1) then
    ex=1.-exp(-dtau(:,:,1))
  else
    ex=dtau(:,:,1)*(1.-.5*dtau(:,:,1)*(1.-.3333*dtau(:,:,1)))
  endif

  !---------------------------------------------------------------------
  ! Boundary condition at j=1.
  !---------------------------------------------------------------------
  sp2(:,:,1)=dtau(:,:,2)*(1.+.5*dtau(:,:,2))
  sp3(:,:,1)=-1.
  q  (:,:,1)=s(:,:,1)*dtau(:,:,2)*(.5*dtau(:,:,2)+ex)
  
  !---------------------------------------------------------------------
  ! Matrix elements for j>2, j<np: [3d+3s]
  !---------------------------------------------------------------------
  do j=2, nz-1
    do k=1,ny
      do l=1,nx
        sp2(l,k,j)=1.
        dinv=2./(dtau(l,k,j+1)+dtau(l,k,j))
        sp1(l,k,j)=-dinv/dtau(l,k,j)
        sp3(l,k,j)=-dinv/dtau(l,k,j+1)
      end do
    end do
  end do
    
  !---------------------------------------------------------------------
  ! Right hand sides, for j>2, j<n1
  !---------------------------------------------------------------------
  n1p=max0(2,nz-3)
  do j=2,n1p-1
    do k=1,ny
      do l=1,nx
        q(l,k,j)=s(l,k,j)
      end do
    end do
  end do
  
  !---------------------------------------------------------------------
  ! RHS for j=n1p.  The equation is still in p(j), but q(j+1) is now
  ! q(j+1)=p(j+1)-s(j+1), so to get p(j+1) we have to add s(j+1) to the
  ! LHS, that is subtract it from the RHS.
  !---------------------------------------------------------------------
  if (n1p.le.nz) then
    do l=1,nx
      do k=1,ny
        q(l,k,n1p)=s(l,k,n1p)-sp3(l,k,n1p)*s(l,k,n1p+1)
      end do
    end do
  endif
  !---------------------------------------------------------------------
  ! RHS for j=n1p+1.  The equation is now in Q: Q" = Q - S", but q(j-1)
  ! is now p(j-1)=q(j-1)+s(j-1), so we must add sp1*s(l,k,j-1) to the RHS,
  ! relative to the case below.
  !---------------------------------------------------------------------
  if (n1p+1.le.nz) then
    do l=1,nx
      do k=1,ny
        q(l,k,n1p+1)=sp1(l,k,n1p+1)* s(l,k,n1p+1) &
                  +sp3(l,k,n1p+1)*(s(l,k,n1p+1)-s(l,k,n1p+2))
      end do
    end do
  endif
  
  !---------------------------------------------------------------------
  ! j>n1p+1,<np [2m+2a]
  !---------------------------------------------------------------------
  do j=n1p+2,nz-1
    do k=1,ny
      do l=1,nx
        q(l,k,j)=sp1(l,k,j)*(s(l,k,j)-s(l,k,j-1)) &
                +sp3(l,k,j)*(s(l,k,j)-s(l,k,j+1))
      end do
    end do
  end do
  
  !---------------------------------------------------------------------
  ! j=np
  !---------------------------------------------------------------------
  sp2(:,:,nz)=1.
    q(:,:,nz)=0.0
    
  !---------------------------------------------------------------------
  !  Eliminate subdiagonal, save factors in sp1 [1d+2m+4a]
  !---------------------------------------------------------------------
  do j=1,nz-2
    do k=1,ny
      do l=1,nx
        a=-sp1(l,k,j+1)/(sp2(l,k,j)-sp3(l,k,j))
        q(l,k,j+1)=q(l,k,j+1)+a*q(l,k,j)
        sp2(l,k,j+1)=sp2(l,k,j+1)+a*sp2(l,k,j)
        sp2(l,k,j)=sp2(l,k,j)-sp3(l,k,j)
      end do
    end do
  end do
  
  sp2(:,:,nz-1)=sp2(:,:,nz-1)-sp3(:,:,nz-1)
  
  !---------------------------------------------------------------------
  ! Backsubstitute [1d+1m+1a]
  !---------------------------------------------------------------------
  do j=nz-1,1,-1
    do k=1,ny
      do l=1,nx
        q(l,k,j)=(q(l,k,j)-sp3(l,k,j)*q(l,k,j+1))/sp2(l,k,j)
      end do
    end do
  end do
  
  !---------------------------------------------------------------------
  ! Subtract S in the region where the equation is in P
  !---------------------------------------------------------------------
  do j=1,n1p
    do k=1,ny
      do l=1,nx
        q(l,k,j)=q(l,k,j)-s(l,k,j)
      end do
    end do
  end do
  
  !---------------------------------------------------------------------
  ! Surface intensity.
  !---------------------------------------------------------------------
  if (doxi) then
    do k=1,ny
      do l=1,nx
        qre(l,k)=2.*(q(l,k,2)+s(l,k,2))-q(l,k,1)
        qfw(l,k)=2.*(q(l,k,nz-1)+s(l,k,nz-1))-q(l,k,nz)
      end do
    end do
  end if
  deallocate(sp1,sp2, sp3, ex)

END SUBROUTINE feautrier

END MODULE rt_solver_mod
