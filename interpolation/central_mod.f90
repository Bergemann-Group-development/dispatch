!> =============================================================================
!> Compute time derivative from a sequence of time slices, using Lagrange
!> interpolation
!> =============================================================================
MODULE central_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE lagrange_mod
  implicit none
  private
  type, public:: central_t
  contains
    procedure, nopass:: time_derivative
    procedure, nopass:: time_interpolate
  end type
CONTAINS

!> =============================================================================
!> =============================================================================
SUBROUTINE time_derivative (patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=KindScalarVar), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  integer, parameter:: order=2
  real(8), dimension(order+1):: w
  integer:: i, i1, i2
  !-----------------------------------------------------------------------------
  times = patch%t(patch%iit(1:patch%nt-1))
  call lagrange%deriv_sequence_weights (patch%out_next, times, i1, i2, w, order)
  if (io%verbose>2) then
    print '(a,10f12.6)', '  times:', times, patch%out_next
    print '(a,10f12.6)', 'weights:', w, sum(w)
  end if
  buffer = 0.0
  do i=i1,i2
    buffer = buffer + w(1+i-i1)*mem(:,:,:,patch%iit(i))
  end do
END SUBROUTINE time_derivative

!> =============================================================================
!> =============================================================================
SUBROUTINE time_interpolate (patch, iv, buffer, order)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:):: buffer
  integer, optional:: order
  !.............................................................................
  real(8):: times(patch%nt-1)
  integer:: order_l=2
  real(8), dimension(:), allocatable:: w
  integer:: i, i1, i2
  !-----------------------------------------------------------------------------
  if (present(order)) order_l = order
  times = patch%t(patch%iit(1:patch%nt-1))
  allocate (w(order_l+1))
  call lagrange%sequence_weights (patch%out_next, times, i1, i2, w, order)
  if (io%verbose>2) then
    print '(a,10f12.6)', '  times:', times, patch%out_next
    print '(a,10f12.6)', 'weights:', w, sum(w)
  end if
  do i=i1,i2
    if (i==i1) then
      buffer =          w(1+i-i1)*patch%mem(:,:,:,iv,patch%iit(i),1)
    else
      buffer = buffer + w(1+i-i1)*patch%mem(:,:,:,iv,patch%iit(i),1)
    end if
  end do
  deallocate (w)
END SUBROUTINE time_interpolate

END MODULE central_mod
