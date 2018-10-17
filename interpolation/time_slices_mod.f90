!> =============================================================================
!> Compute time derivative from a sequence of time slices, using Lagrange
!> interpolation
!> =============================================================================
MODULE time_slices_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE lagrange_mod
  implicit none
  private
  type, public:: time_slices_t
  contains
    procedure, nopass:: test
    procedure, nopass:: interpolate
    procedure, nopass:: derivative1
    procedure, nopass:: derivative2
  end type
  integer, parameter:: verbose=0
  integer, parameter:: order=4
  type(time_slices_t), public:: time_slices
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE test()
  call lagrange%test
END SUBROUTINE test

!> =============================================================================
!> Time interpolation
!> =============================================================================
SUBROUTINE interpolate (patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=KindScalarVar), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  real(8), dimension(order+1):: w
  integer:: i1, i2
  !-----------------------------------------------------------------------------
  if (verbose > 0) print *, patch%id, 'time_slices_t%interpolate'
  times = patch%t(patch%iit(1:patch%nt-1))
  call lagrange%sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE interpolate

!> =============================================================================
!> First time derivate
!> =============================================================================
SUBROUTINE derivative1 (patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=KindScalarVar), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  real(8), dimension(order+1):: w
  integer:: i1, i2
  !-----------------------------------------------------------------------------
  if (verbose > 0) print *, patch%id, 'time_slices_t%derivative1'
  times = patch%t(patch%iit(1:patch%nt-1))
  call lagrange%deriv_sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE derivative1

!> =============================================================================
!> Second time derivate
!> =============================================================================
SUBROUTINE derivative2 (patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=KindScalarVar), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  real(8), dimension(order+1):: w
  integer:: i1, i2
  !-----------------------------------------------------------------------------
  if (verbose > 0) print *, patch%id, 'time_slices_t%derivative2'
  times = patch%t(patch%iit(1:patch%nt-1))
  call lagrange%deriv2_sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE derivative2

!> =============================================================================
!> Weighted sum of time slice values
!> =============================================================================
SUBROUTINE weighted_sum (patch, times, w, i1, i2, mem, buffer)
  class(patch_t):: patch
  real(8), dimension(:):: times, w
  integer:: i1, i2
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=KindScalarVar), dimension(:,:,:):: buffer
  !.............................................................................
  integer:: i, li(3), ui(3)
  !-----------------------------------------------------------------------------
  if (verbose>0) then
    print '(a,2i4,2x,10i4)', 'indices:', i1, i2, patch%iit(i1:i2)
    print '(a,f12.6,2x,10f12.6)'   , '  times:', patch%out_next, times(i1:i2)
    print '(a,f12.6,2x,10f12.3)'   , 'weights:', sum(w), w
  end if
  buffer = 0.0
  if (io%guard_zones) then
    do i=i1,i2
      buffer = buffer + w(1+i-i1)*mem(:,:,:,patch%iit(i))
    end do
  else
    li = patch%mesh%li
    ui = patch%mesh%ui
    do i=i1,i2
      buffer = buffer + w(1+i-i1)*mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),patch%iit(i))
    end do
  end if
END SUBROUTINE weighted_sum

END MODULE time_slices_mod
