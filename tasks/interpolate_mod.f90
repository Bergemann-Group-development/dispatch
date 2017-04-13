!> =============================================================================
!> =============================================================================
MODULE interpolate_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE kinds_mod
  implicit none
  private

  ! a generic interface for interpolation functions.
  abstract interface
    function interpolator_interface (q,i,j,k,iv,aa,bb,cc) result(qtilde)
      import KindScalarVar
      real(kind=KindScalarVar), intent(in) :: q(:,:,:,:)
      integer, intent(in) :: i, j, k, iv
      real(kind=KindScalarVar), intent(in) :: aa, bb, cc
      real(kind=KindScalarVar) :: qtilde
    end function
  end interface

PUBLIC SelectInterpolator, interpolator_interface
CONTAINS

!> =============================================================================
!> Create and return a procedure pointer for an interpolation function.
!> =============================================================================
FUNCTION SelectInterpolator (iorder) result(interp)
  procedure(interpolator_interface), pointer :: interp
  integer, intent(in) :: iorder
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_t%SelectInterpolator')

  select case (iorder)
  case (0) ! donor cell
    interp => interpolator_donor_cell
  case (1) ! linear van Leer
    interp => interpolator_van_Leer_3d
  case (-1) ! straight, un-limited trilinear interpolation
    interp => interpolator_trilinear
  case default
    call mpi%abort("The order you have selected has not been implemented. Abort!")
  end select

  call trace_end
END FUNCTION SelectInterpolator

!> =============================================================================
!> Donor cell (direct injection) interpolation.
!> =============================================================================
FUNCTION interpolator_donor_cell (q,i,j,k,iv,aa,bb,cc) result(qtilde)
  real(kind=KindScalarVar), intent(in) :: q(:,:,:,:)
  integer, intent(in) :: i, j, k, iv
  real(kind=KindScalarVar), intent(in) :: aa, bb, cc
  real(kind=KindScalarVar) :: qtilde
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_mod::donor_cell', 2)
  
  qtilde = q(i,j,k,iv)
  
  call trace_end
END FUNCTION interpolator_donor_cell

!> =============================================================================
!> Piece-wise monotonic, conservative linear interpolation (PLI; van Leer 1977).
!> aa, bb, cc are the linear interpolation weights in each direction.
!> p is the patch we are interpolating *FROM*.
!>
!> Note that, for an interpolation in 3-D, three 1-D, linear van Leer
!> interpolations are *not* guaranteed to be monotonic!
!> =============================================================================
FUNCTION interpolator_van_Leer_3d (q,i,j,k,iv,aa,bb,cc) result(qtilde)
  real(kind=KindScalarVar), intent(in) :: q(:,:,:,:)
  integer, intent(in) :: i, j, k, iv
  real(kind=KindScalarVar), intent(in) :: aa, bb, cc
  real(kind=KindScalarVar) :: qtilde
  !-----------------------------------------------------------------------------
  logical :: ldim(3), lzc=.false., lfc(3)=.false.
  real(8) :: qijk, dq1, dq2, dq3, qp, qm, dp, dm, aych1, aych2, aych3
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_mod::van_Leer', 2)

  ! determine what type of variable we are interpolating
!  select type (p)
!  type is (mhd_t) ! hardwired stagger mapping; FIXME!
!    if (iv <= 2) lzc = .true.
!    if (iv == 3 .or. iv == 6) lfc(1) = .true.
!    if (iv == 4 .or. iv == 7) lfc(2) = .true.
!    if (iv == 5 .or. iv == 8) lfc(3) = .true.
!  end select
  lzc = .true.  ! FIXME!

  ldim(:) = .false.
  if (lzc) ldim(:) = .true.
  if (lfc(2) .or. lfc(3)) ldim(1) = .true.
  if (lfc(3) .or. lfc(1)) ldim(2) = .true.
  if (lfc(1) .or. lfc(2)) ldim(3) = .true.

  qijk = q(i,j,k,iv)

  ! --- First, determine the van Leer slopes in each direction (as needed) ---
  dq1 = 0.0d0
  dq2 = 0.0d0
  dq3 = 0.0d0

  ! 1st dimension
  if (ldim(1)) then
    qp = q(i+1,j,k,iv)
    qm = q(i-1,j,k,iv)
    dp = qp - qijk
    dm = qijk - qm
    if ( dp * dm > 0.0d0) then
      dq1 = (dp * dm) / (qp - qm)
    else
      dq1 = 0.0d0
    end if
  end if

  ! 2nd dimension
  if (ldim(2)) then
    qp = q(i,j+1,k,iv)
    qm = q(i,j-1,k,iv)
    dp = qp - qijk
    dm = qijk - qm
    if ( dp * dm > 0.0d0) then
      dq2 = (dp * dm) / (qp - qm)
    else
      dq2 = 0.0d0
    end if
  end if

  ! 3rd dimension
  if (ldim(3)) then
    qp = q(i,j,k+1,iv)
    qm = q(i,j,k-1,iv)
    dp = qp - qijk
    dm = qijk - qm
    if ( dp * dm .gt. 0.0d0) then
      dq3 = (dp * dm) / (qp - qm)
    else
      dq3 = 0.0d0
    end if
  end if

  ! --- Second, apply a correction (if any) for curvilinear coords. ---
  ! **Currently disabled**
  aych1 = 0.0d0
  aych2 = 0.0d0
  aych3 = 0.0d0
  
!  f31zi = fthreedv ( irr, dvl1a , ifu )
!  f31fi = fthreedv ( irr, h2b * dx1a, ifu )
!  f32zj = fthreedv ( irr, dvl2a , jfu )
!  f32fj = fthreedv ( irr, dx2a  , jfu )
!  f33k  = fthreedv ( irr, dvl3a, kfu )
  
!  ! Determine the correction factor.
!  if ( lzc .or. lfc2 .or. lfc3 ) aych1 = dvl1ai(i) * f31z(i)
!  if ( iv .eq. b2idx .or. iv .eq. b3idx ) aych1 = h31bi(i) * dx1ai(i) * f31f(i)
!  if ( lzc .or. lfc1 .or. lfc3 ) aych2 = dvl2ai(j) * f32z(j)
!  if ( iv .eq. b3idx ) aych2 = dx2ai(j) * f32f(j)
!  if ( lzc .or. lfc1 .or. lfc2 ) aych3 = dvl3ai(k) * f33 (k)

!  ! Correct when montonicity does not hold.
!  if ( lzc .or. lfc3 .or. lfc2 ) then
!    q1 = abs ( dq1 * ( 1.0d0 + aych1 ) )
!    q2 = abs ( dq1 * ( 1.0d0 - aych1 ) )
!    q3 = abs (q(i  ,j,k) - q(i-1,j,k))
!    q4 = abs (q(i+1,j,k) - q(i  ,j,k))
!    if ( q1 > q3 ) then
!      dq1 = sign (1.0d0, dq1) * q3 / ( 1.0d0 + aych1 )
!    else if ( q2 > q4 ) then
!      dq1 = sign (1.0d0, dq1) * q4 / ( 1.0d0 - aych1 )
!    end if
!  end if

!  if ( lzc .or. lfc3 .or. lfc1 ) then
!    q1 = abs ( dq2 * ( 1.0d0 + aych2 ) )
!    q2 = abs ( dq2 * ( 1.0d0 - aych2 ) )
!    q3 = abs (q(i,j  ,k) - q(i,j-1,k))
!    q4 = abs (q(i,j+1,k) - q(i,j  ,k))
!    if ( q1 > q3 ) then
!      dq2 = sign (1.0d0, dq2) * q3 / ( 1.0d0 + aych2 )
!    else if ( q2 > q4 ) then
!      dq2 = sign (1.0d0, dq2) * q4 / ( 1.0d0 - aych2 )
!    end if
!  end if

!  if ( lzc .or. lfc1 .or. lfc2 ) then
!    q1 = abs ( dq3 * ( 1.0d0 + aych3 ) )
!    q2 = abs ( dq3 * ( 1.0d0 - aych3 ) )
!    q3 = abs (q(i,j,k  ) - q(i,j,k-1))
!    q4 = abs (q(i,j,k+1) - q(i,j,k  ))
!    if ( q1 > q3 ) then
!      dq3 = sign (1.0d0, dq3) * q3 / ( 1.0d0 + aych3 )
!    else if ( q2 > q4 ) then
!      dq3 = sign (1.0d0, dq3) * q4 / ( 1.0d0 - aych3 )
!    end if
!  end if

  ! Finally, determine the interpolated value.
  qtilde = qijk + dq1 * ( aa - aych1 ) &
                + dq2 * ( bb - aych2 ) &
                + dq3 * ( cc - aych3 )

  call trace_end
END FUNCTION interpolator_van_Leer_3d

!> =============================================================================
!> Plain, un-limited, trilinear interpolation.
!> =============================================================================
FUNCTION interpolator_trilinear (q,i,j,k,iv,aa,bb,cc) result(qtilde)
  real(kind=KindScalarVar), intent(in) :: q(:,:,:,:)
  integer, intent(in) :: i, j, k, iv
  real(kind=KindScalarVar), intent(in) :: aa, bb, cc
  real(kind=KindScalarVar) :: qtilde, p(3), m(3)
  integer :: i1, j1, k1
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_mod::trilinear', 2)
  i1 = merge(1,0,size(q,1) > 1)
  j1 = merge(1,0,size(q,2) > 1)
  k1 = merge(1,0,size(q,3) > 1)

  p(1) = aa
  p(2) = bb
  p(3) = cc
  m = 1.0 - p
  qtilde = m(3) * (m(2) * (m(1) * q(i  ,j   ,k   ,iv) + p(1) * q(i+i1,j   ,k   ,iv))  + &
                   p(2) * (m(1) * q(i  ,j+j1,k   ,iv) + p(1) * q(i+i1,j+j1,k   ,iv))) + &
           p(3) * (m(2) * (m(1) * q(i  ,j   ,k+k1,iv) + p(1) * q(i+i1,j   ,k+k1,iv))  + &
                   p(2) * (m(1) * q(i  ,j+j1,k+k1,iv) + p(1) * q(i+i1,j+j1,k+k1,iv)))

  call trace_end
END FUNCTION interpolator_trilinear

END MODULE interpolate_mod
