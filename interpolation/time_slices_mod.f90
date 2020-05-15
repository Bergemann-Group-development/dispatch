!===============================================================================
!> Compute time derivative from a sequence of time slices, using Lagrange
!> interpolation
!===============================================================================
  MODULE time_slices_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE patch_mod
  USE kinds_mod
  USE lagrange_mod
  implicit none
  private
  type, public:: time_slices_t
    integer:: order=3
  contains
    procedure:: init
    procedure, nopass:: interpolate4
    procedure, nopass:: interpolate8
    generic:: interpolate => interpolate4, interpolate8
    procedure, nopass:: derivative1
    procedure, nopass:: derivative2
  end type
  integer, save:: verbose=0
  integer, save:: order=3
  integer, save:: id_debug=-1
  type(time_slices_t), public:: time_slices
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, nt)
  class(time_slices_t):: self
  integer:: nt
  !.............................................................................
  integer:: iostat
  namelist /time_slices_params/ verbose, order, id_debug
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('time_slices_t%init')
  call lagrange%init
  call io%lock%set ('time_slices_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, time_slices_params, iostat=iostat)
    call io%namelist_warning ('time_slices_params', iostat)
    order = min(order,nt-2)
    write (io%namelist, time_slices_params)
    call lagrange%test
    first_time = .false.
  end if
  call io%lock%unset ('time_slices_t%init')
  self%order = order
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Time interpolation.  mem() is assumed to contain nt time slices of a 3-D
!> variable, with times consistent with the times produced by patch%timeslots()
!===============================================================================
!> Single precision buffer version
!===============================================================================
SUBROUTINE interpolate4 (patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=4), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  integer:: iit(patch%nt-1)
  real(8), dimension(order+1):: w
  integer:: i1, i2
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(io_unit%log,*) patch%id, 'time_slices_t%interpolate4'
  call patch%timeslots (iit, times)
  call lagrange%sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum4 (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE interpolate4

SUBROUTINE interpolate8(patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=8), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  integer:: iit(patch%nt-1)
  real(8), dimension(order+1):: w
  integer:: i1, i2
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(io_unit%log,*) patch%id, 'time_slices_t%interpolate8'
  call patch%timeslots (iit, times)
  call lagrange%sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum8 (patch, times, w, i1, i2, mem, buffer)
  call trace%end ()
END SUBROUTINE interpolate8

!===============================================================================
!> First time derivate
!===============================================================================
SUBROUTINE derivative1 (patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=4), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  real(8), dimension(order+1):: w
  integer:: i1, i2
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(io_unit%log,*) patch%id, 'time_slices_t%derivative1'
  times = patch%t(patch%iit(1:patch%nt-1))
  call lagrange%deriv_sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum4 (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE derivative1

!===============================================================================
!> Second time derivate
!===============================================================================
SUBROUTINE derivative2 (patch, mem, buffer)
  class(patch_t):: patch
  integer:: iv
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=4), dimension(:,:,:):: buffer
  !.............................................................................
  real(8):: times(patch%nt-1)
  real(8), dimension(order+1):: w
  integer:: i1, i2
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(io_unit%log,*) patch%id, 'time_slices_t%derivative2'
  times = patch%t(patch%iit(1:patch%nt-1))
  call lagrange%deriv2_sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum4 (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE derivative2

!===============================================================================
!> Weighted sum of time slice values
!===============================================================================
SUBROUTINE weighted_sum4 (patch, times, w, i1, i2, mem, buffer)
  class(patch_t):: patch
  real(8), dimension(:):: times, w
  real(8):: t(patch%nt)
  integer:: i1, i2, iit(patch%nt)
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=4), dimension(:,:,:):: buffer
  !.............................................................................
  integer:: i
  !-----------------------------------------------------------------------------
  if (verbose>0 .or. patch%id==id_debug   &
                .or. patch%id==id_debug-1 &
                .or. patch%id==id_debug+1) then
    t = patch%t(patch%iit)
    write (io_unit%output,'(a,i6,2i4)')  'id,i1,i2:', patch%id, i1, i2
    write (io_unit%output,'(a,9i4)')     '     iit:', patch%iit
    write (io_unit%output,'(a,9f14.6)')  '       t:', patch%t
    write (io_unit%output,'(a,11f14.6)') '  t(iit):', patch%out_next, t(i1:i2)
    write (io_unit%output,'(a,11f14.6)') '   times:', patch%out_next, times(i1:i2)
    write (io_unit%output,'(a,11f14.3)') ' weights:', sum(w), w
  end if
  buffer = 0.0
  do i=i1,i2
    buffer = buffer + w(1+i-i1)*mem(:,:,:,i)
  end do
END SUBROUTINE weighted_sum4

SUBROUTINE weighted_sum8 (patch, times, w, i1, i2, mem, buffer)
  class(patch_t):: patch
  real(8), dimension(:):: times, w
  integer:: i1, i2
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=8), dimension(:,:,:):: buffer
  !.............................................................................
  integer:: i
  !-----------------------------------------------------------------------------
  if (verbose>0 .or. patch%id==id_debug) then
    write (io_unit%output,'(a,2i4,2x,10i4)')      'indices:', i1, i2, patch%iit(i1:i2)
    write (io_unit%output,'(a,f12.6,2x,10f14.6)') '  times:', patch%out_next, times(i1:i2)
    write (io_unit%output,'(a,f12.6,2x,10f14.3)') 'weights:', sum(w), w
  end if
  buffer = 0.0
  do i=i1,i2
    buffer = buffer + w(1+i-i1)*mem(:,:,:,i)
  end do
END SUBROUTINE weighted_sum8

END MODULE time_slices_mod
