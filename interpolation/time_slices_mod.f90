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
    procedure, nopass:: interpolate
    procedure, nopass:: derivative1
    procedure, nopass:: derivative2
  end type
  integer, save:: verbose=0
  integer, save:: order=3
  integer, save:: id_debug=0
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
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, time_slices_params, iostat=iostat)
    order = min(order,nt-2)
    write (io%output, time_slices_params)
    write (io_unit%nml, time_slices_params)
    flush (io_unit%nml)
    call lagrange%test
  end if
  !$omp end critical (input_cr)
  self%order = order
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Time interpolation.  mem() is assumed to contain nt time slices of a 3-D
!> variable, with times consistent with the times produced by patch%timeslots()
!===============================================================================
SUBROUTINE interpolate (patch, mem, buffer)
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
  if (verbose > 0) write(io_unit%log,*) patch%id, 'time_slices_t%interpolate'
  call patch%timeslots (iit, times)
  call lagrange%sequence_weights (patch%out_next, times, i1, i2, w, order)
  call weighted_sum (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE interpolate

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
  call weighted_sum (patch, times, w, i1, i2, mem, buffer)
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
  call weighted_sum (patch, times, w, i1, i2, mem, buffer)
END SUBROUTINE derivative2

!===============================================================================
!> Weighted sum of time slice values
!===============================================================================
SUBROUTINE weighted_sum (patch, times, w, i1, i2, mem, buffer)
  class(patch_t):: patch
  real(8), dimension(:):: times, w
  integer:: i1, i2
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem
  real(kind=4), dimension(:,:,:):: buffer
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
END SUBROUTINE weighted_sum

END MODULE time_slices_mod
