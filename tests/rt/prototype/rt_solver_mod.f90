!> =============================================================================
!> Module with RT solvers -- for now just an integral solver
!> =============================================================================
MODULE rt_solver_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: rt_solver_t
    procedure(integral_forward), pointer:: transfer => null()
    procedure(integral_forward), pointer:: forward  => null()
    procedure(integral_forward), pointer:: reverse  => null()
  contains
    procedure:: init
    procedure:: integral_forward
    procedure:: integral_reverse
  end type
  type (rt_solver_t), public:: rt_solver
CONTAINS

!> =============================================================================
!> Choose RT solver method
!> =============================================================================
SUBROUTINE init (self)
  class (rt_solver_t):: self
  self%forward  => integral_forward
  self%reverse  => integral_reverse
END SUBROUTINE init

!> =============================================================================
!> Integral solution method
!> =============================================================================
SUBROUTINE integral_forward (self, dtau, exp0, S, qfw0, qfw, qfw1)
  implicit none
  class (rt_solver_t):: self
  real, dimension(:,:) :: dtau, exp0, S, qfw
  real, dimension(:)   :: qfw0, qfw1
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
  real, dimension(size(s,1),size(s,2)) :: ex0, ds1, d1s, d2s
  real, dimension(size(s,1))           :: ctau
  integer i, m, n
!-----------------------------------------------------------------------
!
  m = size(s,1)
  n = size(s,2)
  do i=2,n
    ex0(:,i) = exp0(:,i)                                !          1e
    ds1(:,i) = (s(:,i)-s(:,i-1))/dtau(:,i)              ! 1d    1a      ! half-point centered dS/dtau
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
  d2s(:,1) = d2s(:,2)                                                   ! assume constant d2S/dtau2
  d1s(:,1) = d1s(:,2) - dtau(:,2)*d2s(:,2)                              ! extrapolate to 1
  d2s(:,n) = d2s(:,n-1)                                                 ! assume constant d2S/dtau2
  d1s(:,n) = d1s(:,n-1) + dtau(:,n)*d2s(:,n)                            ! extrapolate to n
!
!  Forward rays (2a+3m)
!
  qfw(:,1) = qfw0(:)                                                    ! initial forward Q
  do i=2,n
    qfw(:,i) = qfw(:,i-1)*ex0(:,i) &
           - d1s(:,i)*(1.-ex0(:,i)) &
           + d2s(:,i)*(1.-ex0(:,i)*(1.+dtau(:,i)))      !    3m 2a      ! formal solution forward
  end do
  qfw1(:) = qfw(:,n)                                                    ! forward boundary intensity
END SUBROUTINE integral_forward                         ! 2d 13m 11a 1e ! 26f + 1e

!> =============================================================================
!> Integral solution method
!> =============================================================================
SUBROUTINE integral_reverse (self, dtau, exp0, S, qrv1, qrv, qrv0)
  implicit none
  class (rt_solver_t):: self
  real, dimension(:,:) :: dtau, exp0, S, qrv
  real, dimension(:)   :: qrv0, qrv1
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
  real, dimension(size(s,1),size(s,2)) :: ex0, ds1, d1s, d2s
  real, dimension(size(s,1))           :: ctau
  integer i, m, n
!-----------------------------------------------------------------------
!
  m = size(s,1)
  n = size(s,2)
  do i=2,n
    ex0(:,i) = exp0(:,i)                                !          1e
    ds1(:,i) = (s(:,i)-s(:,i-1))/dtau(:,i)              ! 1d    1a      ! half-point centered dS/dtau
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
  d2s(:,1) = d2s(:,2)                                                   ! assume constant d2S/dtau2
  d1s(:,1) = d1s(:,2) - dtau(:,2)*d2s(:,2)                              ! extrapolate to 1
  d2s(:,n) = d2s(:,n-1)                                                 ! assume constant d2S/dtau2
  d1s(:,n) = d1s(:,n-1) + dtau(:,n)*d2s(:,n)                            ! extrapolate to n
!
!  Forward rays (2a+3m)
!
  qrv(:,n) = qrv1(:)                                                    ! initial reverse Q
  do i=n,2,-1
    qrv(:,i-1) = qrv(:,i)*ex0(:,i) &
           - d1s(:,i)*(1.-ex0(:,i)) &
           + d2s(:,i)*(1.-ex0(:,i)*(1.+dtau(:,i)))      !    3m 2a      ! formal reverse solution
  end do
  qrv0(:) = qrv(:,1)                                                    ! reverse boundary intensity
END SUBROUTINE integral_reverse                         ! 2d 13m 11a 1e ! 26f + 1e

END MODULE rt_solver_mod
