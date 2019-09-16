!===============================================================================
!===============================================================================
MODULE interpolate_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE kinds_mod
  USE index_mod
  implicit none
  private

  ! a generic interface for interpolation functions.
  abstract interface
    function interpolator_interface (q,i,j,k,w) result(qtilde)
      import KindScalarVar
      real(kind=KindScalarVar), intent(in):: q(:,:,:)
      integer, intent(in):: i, j, k
      real, intent(in):: w(3)
      real(kind=KindScalarVar):: qtilde
    end function
  end interface

  type:: interpolator_t
    integer:: order_interpolator
    procedure(interpolator_interface), pointer, nopass:: interpolator => null()
  contains
    procedure, nopass:: trilinear_log
    procedure, nopass:: trilinear_pv
    procedure, nopass:: four_d
    procedure, nopass:: four_d_log
    procedure, nopass:: four_d_pv
  end type
  type(interpolator_t):: interpolator

PUBLIC SelectInterpolator, interpolator, interpolator_interface, interpolator_unsigned
CONTAINS

!===============================================================================
!> Create and return a procedure pointer for an interpolation function.
!===============================================================================
FUNCTION SelectInterpolator (iorder) result(interp)
  procedure(interpolator_interface), pointer :: interp
  integer, intent(in) :: iorder
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_t%SelectInterpolator')

  select case (iorder)
  case (0) ! donor cell
    interp => interpolator_donor_cell
    interpolator%interpolator => interpolator_donor_cell
  case (1) ! linear van Leer
    interp => interpolator_van_Leer_3d
    interpolator%interpolator => interpolator_van_Leer_3d
  case (-1) ! straight, un-limited trilinear interpolation
    interp => interpolator_trilinear
    interpolator%interpolator => interpolator_trilinear
  case default
    call mpi%abort("The order you have selected has not been implemented. Abort!")
  end select
  interpolator%order_interpolator = iorder

  call trace_end
END FUNCTION SelectInterpolator

!===============================================================================
!> Donor cell (direct injection) interpolation.
!===============================================================================
FUNCTION interpolator_donor_cell (q,i,j,k,w) result(qtilde)
  real(kind=KindScalarVar), intent(in):: q(:,:,:)
  integer, intent(in):: i, j, k
  real, intent(in):: w(3)
  real(kind=KindScalarVar):: qtilde
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_mod::donor_cell', 2)
  
  qtilde = q(i,j,k)
  
  call trace_end
END FUNCTION interpolator_donor_cell

!===============================================================================
!> Piece-wise monotonic, conservative linear interpolation (PLI; van Leer 1977).
!> aa, bb, cc are the linear interpolation weights in each direction.
!> p is the patch we are interpolating *FROM*.
!>
!> Note that, for an interpolation in 3-D, three 1-D, linear van Leer
!> interpolations are *not* guaranteed to be monotonic!
!===============================================================================
FUNCTION interpolator_van_Leer_3d (q,i,j,k,w) result(qtilde)
  real(kind=KindScalarVar), intent(in):: q(:,:,:)
  integer, intent(in):: i, j, k
  real, intent(in):: w(3)
  real(kind=KindScalarVar):: qtilde
  !-----------------------------------------------------------------------------
  logical :: ldim(3), lzc=.false., lfc(3)=.false.
  real(8) :: qijk, dq1, dq2, dq3, qp, qm, dp, dm, aych1, aych2, aych3
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_mod::van_Leer', 1, itimer=itimer)
  ldim(:) = .false.
  if (lzc) ldim(:) = .true.
  if (lfc(2) .or. lfc(3)) ldim(1) = .true.
  if (lfc(3) .or. lfc(1)) ldim(2) = .true.
  if (lfc(1) .or. lfc(2)) ldim(3) = .true.
  if (size(q,1) <= 1) ldim(1) = .false.
  if (size(q,2) <= 1) ldim(2) = .false.
  if (size(q,3) <= 1) ldim(3) = .false.

  qijk = q(i,j,k)

  ! --- First, determine the van Leer slopes in each direction (as needed) ---
  dq1 = 0.0d0
  dq2 = 0.0d0
  dq3 = 0.0d0

  ! 1st dimension
  if (ldim(1)) then
    qp = q(i+1,j,k)
    qm = q(i-1,j,k)
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
    qp = q(i,j+1,k)
    qm = q(i,j-1,k)
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
    qp = q(i,j,k+1)
    qm = q(i,j,k-1)
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
  
  ! Finally, determine the interpolated value.
  qtilde = qijk + dq1 * ( w(1) - aych1 ) &
                + dq2 * ( w(2) - aych2 ) &
                + dq3 * ( w(3) - aych3 )

  call trace_end (itimer)
END FUNCTION interpolator_van_Leer_3d

!===============================================================================
!> Plain, un-limited, trilinear interpolation.
!===============================================================================
FUNCTION interpolator_trilinear (q,i,j,k,w) result(qtilde)
  real(kind=KindScalarVar), intent(in):: q(:,:,:)
  integer, intent(in):: i, j, k
  real, intent(in):: w(3)
  real(kind=KindScalarVar):: qtilde, m(3)
  integer:: i1, j1, k1
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace_begin('interpolate_mod::trilinear', 2, itimer=itimer)
  i1 = merge(1,0,size(q,1) > 1)
  j1 = merge(1,0,size(q,2) > 1)
  k1 = merge(1,0,size(q,3) > 1)
  m = 1.0 - w
  qtilde = m(3) * (m(2) * (m(1) * q(i  ,j   ,k   ) + w(1) * q(i+i1,j   ,k   ))  + &
                   w(2) * (m(1) * q(i  ,j+j1,k   ) + w(1) * q(i+i1,j+j1,k   ))) + &
           w(3) * (m(2) * (m(1) * q(i  ,j   ,k+k1) + w(1) * q(i+i1,j   ,k+k1))  + &
                   w(2) * (m(1) * q(i  ,j+j1,k+k1) + w(1) * q(i+i1,j+j1,k+k1)))
  !call trace_end (itimer)
END FUNCTION interpolator_trilinear

!===============================================================================
!> Plain, un-limited, trilinear interpolation.
!===============================================================================
FUNCTION trilinear_log (q,i,j,k,w) result(qtilde)
  real(kind=KindScalarVar), intent(in):: q(:,:,:)
  integer, intent(in):: i, j, k
  real, intent(in):: w(3)
  real(kind=KindScalarVar):: qtilde, m(3)
  integer:: i1, j1, k1
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace_begin('interpolate_mod::trilinear', 2, itimer=itimer)
  i1 = merge(1,0,size(q,1) > 1)
  j1 = merge(1,0,size(q,2) > 1)
  k1 = merge(1,0,size(q,3) > 1)
  m = 1.0 - w
  qtilde = m(3)*(m(2)*(m(1)*log(q(i  ,j   ,k   )) + w(1)*log(q(i+i1,j   ,k   )))  + &
                 w(2)*(m(1)*log(q(i  ,j+j1,k   )) + w(1)*log(q(i+i1,j+j1,k   )))) + &
           w(3)*(m(2)*(m(1)*log(q(i  ,j   ,k+k1)) + w(1)*log(q(i+i1,j   ,k+k1)))  + &
                 w(2)*(m(1)*log(q(i  ,j+j1,k+k1)) + w(1)*log(q(i+i1,j+j1,k+k1))))
  !call trace_end (itimer)
END FUNCTION trilinear_log

!===============================================================================
!> Plain, un-limited, trilinear interpolation.
!===============================================================================
FUNCTION trilinear_pv (q,d,i,j,k,w) result(qtilde)
  real(kind=KindScalarVar), intent(in):: d(:,:,:), q(:,:,:)
  integer, intent(in):: i, j, k
  real, intent(in):: w(3)
  real(kind=KindScalarVar):: qtilde, m(3)
  integer:: i1, j1, k1
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace_begin('interpolate_mod::trilinear', 2, itimer=itimer)
  i1 = merge(1,0,size(q,1) > 1)
  j1 = merge(1,0,size(q,2) > 1)
  k1 = merge(1,0,size(q,3) > 1)
  m = 1.0 - w
  qtilde = m(3)*(m(2)*(m(1)*q(i   ,j   ,k   )/d(i   ,j   ,k   )   + &
                       w(1)*q(i+i1,j   ,k   )/d(i+i1,j   ,k   ))  + &
                 w(2)*(m(1)*q(i   ,j+j1,k   )/d(i   ,j+j1,k   )   + &
                       w(1)*q(i+i1,j+j1,k   )/d(i+i1,j+j1,k   ))) + &
           w(3)*(m(2)*(m(1)*q(i   ,j   ,k+k1)/d(i   ,j   ,k+k1)   + &
                       w(1)*q(i+i1,j   ,k+k1)/d(i+i1,j   ,k+k1))  + &
                 w(2)*(m(1)*q(i   ,j+j1,k+k1)/d(i   ,j+j1,k+k1)   + &
                       w(1)*q(i+i1,j+j1,k+k1)/d(i+i1,j+j1,k+k1)))
  !call trace_end (itimer)
END FUNCTION trilinear_pv

!===============================================================================
!> Plain, un-limited, trilinear interpolation.
!===============================================================================
FUNCTION interpolator_trilinear_3d (q,i,j,k,w) result(qtilde)
  real(kind=KindScalarVar), intent(in):: q(:,:,:)
  integer, intent(in):: i, j, k
  real, intent(in):: w(3)
  real(kind=KindScalarVar):: qtilde, m(3)
  integer:: i1, j1, k1
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace_begin('interpolate_mod::trilinear', 2, itimer=itimer)
  m = 1.0 - w
  qtilde = m(3) * (m(2) * (m(1) * q(i  ,j  ,k  ) + w(1) * q(i+1,j  ,k  ))  + &
                   w(2) * (m(1) * q(i  ,j+1,k  ) + w(1) * q(i+1,j+1,k  ))) + &
           w(3) * (m(2) * (m(1) * q(i  ,j  ,k+1) + w(1) * q(i+1,j  ,k+1))  + &
                   w(2) * (m(1) * q(i  ,j+1,k+1) + w(1) * q(i+1,j+1,k+1)))
  !call trace_end (itimer)
END FUNCTION interpolator_trilinear_3d

!===============================================================================
!> Do the interpolation under the assumption that the variable is
!> positive-definite (i.e. `patch%unsigned(iv) = .true.`).
!===============================================================================
FUNCTION interpolator_unsigned (q,i,j,k,w) result(qtilde)
  real(kind=KindScalarVar), intent(in):: q(:,:,:)
  integer, intent(in):: i, j, k
  real, intent(in):: w(3)
  real(kind=KindScalarVar):: qtilde
  real(kind=KindScalarVar), allocatable:: logq(:,:,:)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('interpolate_mod::unsigned', 1, itimer=itimer)
  !
  allocate(logq(size(q,1),size(q,2),size(q,3)))
  logq(:,:,:) = log(q(:,:,:))
  qtilde = exp(interpolator%interpolator(logq,i,j,k,w))
  !
  deallocate (logq)
  call trace_end (itimer)
END FUNCTION interpolator_unsigned

!===============================================================================
!> Plain, un-limited, trilinear interpolation.
!===============================================================================
FUNCTION four_d (q, i, j, k, l, w) RESULT (qtilde)
  real(kind=KindScalarVar), intent(in)  :: q(:,:,:,:)
  real(kind=KindScalarVar)              :: qtilde
  integer, intent(in)                   :: i, j, k, l(2)
  real, intent(in)                      :: w(4)
  !.............................................................................
  integer:: l1, l2
  real(kind=KindScalarVar):: m(4)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace_begin('interpolate_mod::trilinear', 2, itimer=itimer)
  l1 = l(1)
  l2 = l(2)
  m = 1.0 - w
   qtilde = &
    m(4)*(m(3)*(m(2)*(m(1)*q(i  ,j  ,k  ,l1) + w(1)*q(i+1,j  ,k  ,l1))   + &
                w(2)*(m(1)*q(i  ,j+1,k  ,l1) + w(1)*q(i+1,j+1,k  ,l1)))  + &
          w(3)*(m(2)*(m(1)*q(i  ,j  ,k+1,l1) + w(1)*q(i+1,j  ,k+1,l1))   + &
                w(2)*(m(1)*q(i  ,j+1,k+1,l1) + w(1)*q(i+1,j+1,k+1,l1)))) + &
    w(4)*(m(3)*(m(2)*(m(1)*q(i  ,j  ,k  ,l2) + w(1)*q(i+1,j  ,k  ,l2))   + &
                w(2)*(m(1)*q(i  ,j+1,k  ,l2) + w(1)*q(i+1,j+1,k  ,l2)))  + &
          w(3)*(m(2)*(m(1)*q(i  ,j  ,k+1,l2) + w(1)*q(i+1,j  ,k+1,l2))   + &
                w(2)*(m(1)*q(i  ,j+1,k+1,l2) + w(1)*q(i+1,j+1,k+1,l2))))
  !call trace_end (itimer)
END FUNCTION four_d

!===============================================================================
!> Plain, un-limited, trilinear interpolation.
!===============================================================================
FUNCTION four_d_pv (q, d1, d2, i, j, k, l, w) RESULT (qtilde)
  real(kind=KindScalarVar), intent(in)  :: q(:,:,:,:), d1(:,:,:), d2(:,:,:)
  real(kind=KindScalarVar)              :: qtilde
  integer, intent(in)                   :: i, j, k, l(2)
  real, intent(in)                      :: w(4)
  !.............................................................................
  integer:: l1, l2
  real(kind=KindScalarVar):: m(4)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace_begin('interpolate_mod::trilinear', 2, itimer=itimer)
  l1 = l(1)
  l2 = l(2)
  m = 1.0 - w
   qtilde = &
    m(4)*(m(3)*(m(2)*(m(1)*q(i  ,j  ,k  ,l1)/d1(i  ,j  ,k  )    + &
                      w(1)*q(i+1,j  ,k  ,l1)/d1(i+1,j  ,k  ))   + &
                w(2)*(m(1)*q(i  ,j+1,k  ,l1)/d1(i  ,j+1,k  )    + &
                      w(1)*q(i+1,j+1,k  ,l1)/d1(i+1,j+1,k  )))  + &
          w(3)*(m(2)*(m(1)*q(i  ,j  ,k+1,l1)/d1(i  ,j  ,k+1)    + &
                      w(1)*q(i+1,j  ,k+1,l1)/d1(i+1,j  ,k+1))   + &
                w(2)*(m(1)*q(i  ,j+1,k+1,l1)/d1(i  ,j+1,k+1)    + &
                      w(1)*q(i+1,j+1,k+1,l1)/d1(i+1,j+1,k+1)))) + &
    w(4)*(m(3)*(m(2)*(m(1)*q(i  ,j  ,k  ,l2)/d2(i  ,j  ,k  )    + &
                      w(1)*q(i+1,j  ,k  ,l2)/d2(i+1,j  ,k  ))   + &
                w(2)*(m(1)*q(i  ,j+1,k  ,l2)/d2(i  ,j+1,k  )    + &
                      w(1)*q(i+1,j+1,k  ,l2)/d2(i+1,j+1,k  )))  + &
          w(3)*(m(2)*(m(1)*q(i  ,j  ,k+1,l2)/d2(i  ,j  ,k+1)    + &
                      w(1)*q(i+1,j  ,k+1,l2)/d2(i+1,j  ,k+1))   + &
                w(2)*(m(1)*q(i  ,j+1,k+1,l2)/d2(i  ,j+1,k+1)    + &
                      w(1)*q(i+1,j+1,k+1,l2)/d2(i+1,j+1,k+1))))
  !call trace_end (itimer)
END FUNCTION four_d_pv
  
!===============================================================================
!> Plain, un-limited, trilinear interpolation.
!===============================================================================
FUNCTION four_d_log (q, i, j, k, l, w) RESULT (qtilde)
  real(kind=KindScalarVar), intent(in)  :: q(:,:,:,:)
  real(kind=KindScalarVar)              :: qtilde
  integer, intent(in)                   :: i, j, k, l(2)
  real, intent(in)                      :: w(4)
  !.............................................................................
  integer:: l1, l2
  real(kind=KindScalarVar):: m(4)
  integer:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace_begin('interpolate_mod::trilinear', 2, itimer=itimer)
  l1 = l(1)
  l2 = l(2)
  m = 1.0 - w
   qtilde = &
    m(4)*(m(3)*(m(2)*(m(1)*log(q(i  ,j  ,k  ,l1)) + w(1)*log(q(i+1,j  ,k  ,l1)))   + &
                w(2)*(m(1)*log(q(i  ,j+1,k  ,l1)) + w(1)*log(q(i+1,j+1,k  ,l1))))  + &
          w(3)*(m(2)*(m(1)*log(q(i  ,j  ,k+1,l1)) + w(1)*log(q(i+1,j  ,k+1,l1)))   + &
                w(2)*(m(1)*log(q(i  ,j+1,k+1,l1)) + w(1)*log(q(i+1,j+1,k+1,l1))))) + &
    w(4)*(m(3)*(m(2)*(m(1)*log(q(i  ,j  ,k  ,l2)) + w(1)*log(q(i+1,j  ,k  ,l2)))   + &
                w(2)*(m(1)*log(q(i  ,j+1,k  ,l2)) + w(1)*log(q(i+1,j+1,k  ,l2))))  + &
          w(3)*(m(2)*(m(1)*log(q(i  ,j  ,k+1,l2)) + w(1)*log(q(i+1,j  ,k+1,l2)))   + &
                w(2)*(m(1)*log(q(i  ,j+1,k+1,l2)) + w(1)*log(q(i+1,j+1,k+1,l2)))))
  !call trace_end (itimer)
END FUNCTION four_d_log

END MODULE interpolate_mod
