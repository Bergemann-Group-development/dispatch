!> =============================================================================
!> Generic validation module. The general idea is to be able to compare two runs
!> at critical points in the sequence of evaluation, by having the first run 
!> write to files, which the 2nd job reads from, and then compare the contents.
!>
!> This is primarily intended for testing and validation of small runs, where 
!> having one file per patch is not a problem.  This simplifies things, while
!> still allowing different kinds of tests:
!>
!> 1. Runs that use a single core, and which are therefore updating patches in
!>    exactly the same order, allowing bit-wise match for unchanged code.
!> 2. Runs that use several core on a single MPI rank, where results are not
!>    expected to be bit-wise identical, but where one can get useful information
!>    about the magnitudes of differences, separately in guard zones and the
!>    interior.
!> 3. The same, but withe several ranks.  Assuming patch IDs uniquely specify
!>    the location in space, by for example being the sum of integers that give
!>    the flattened position in each level.
!> =============================================================================
MODULE validate_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE patch_mod
  implicit none
  private
  type, public:: validate_t
    integer:: id
    integer:: irec=0
    integer:: iout=0
    logical:: ok=.true.
    character(len=8):: state='off'
  contains
    procedure:: init
    procedure:: set
    procedure, private:: dump0
    procedure, private:: dump1
    procedure, private:: dump2
    procedure, private:: dump3
    generic, public:: dump => dump0, dump1, dump2, dump3
  end type
  type(validate_t), public:: validate
CONTAINS

!> =============================================================================
!> Initialize the dump & compare process
!> =============================================================================
SUBROUTINE init (self)
  class(validate_t):: self
  !.............................................................................
  integer:: iostat
  logical, save:: first_time=.true.
  character(len=8), save:: state='off'
  namelist /validate_params/ state
  !.............................................................................
  call trace%begin ('validate_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, validate_params, iostat=iostat)
    write (io%output, validate_params)
  end if
  !$omp end critical (input_cr)
  self%state = state
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Change state of the dump process
!> =============================================================================
SUBROUTINE set (self, state1)
  class(validate_t):: self
  character(len=*):: state1
  !.............................................................................
  self%state = state1
END SUBROUTINE set

!> =============================================================================
!> Switch
!> =============================================================================
SUBROUTINE dump0 (self, id1, f1, label1)
  class(validate_t)    :: self
  integer              :: id1
  real                 :: f1
  character(len=*)     :: label1
  !.............................................................................
  real                 :: f0
  character(len=64)    :: label0
  integer              :: id0, sz, dims
  real                 :: diff
  !-----------------------------------------------------------------------------
  if (self%state=='write') then
    label0 = label1
    write (io_unit%validate) id1, 0, label0
    write (io_unit%validate) 1
    write (io_unit%validate) f1
  else if (self%state=='compare') then
    read (io_unit%validate) id0, sz, label0
    read (io_unit%validate) dims
    read (io_unit%validate) f0
    diff  = f0-f1
    !---------------------------------------------------------------------------
    ! Write out summary of differences
    !---------------------------------------------------------------------------
    write (io_unit%output,1) '    patch ID =',      id0, trim(label0)
    write (io_unit%output,1) '    patch ID =',      id1, trim(label1)
    write (io_unit%output,2) '       aver1 =',      f0, 'aver2 =', f1
    write (io_unit%output,2) '        diff =',    diff
    1 format(a,i12,3x,a)
    2 format(a,1pe12.3,3x,a,1pe12.3)
  end if
END SUBROUTINE dump0

!> =============================================================================
!> Switch
!> =============================================================================
SUBROUTINE dump1 (self, patch, f1, label1)
  class(validate_t)    :: self
  class(patch_t)       :: patch
  integer              :: id1
  real                 :: f1(:)
  character(len=*)     :: label1
  !.............................................................................
  real, allocatable    :: f0(:), diff1(:)
  logical, allocatable :: inside(:)
  character(len=64)    :: label0
  integer              :: id0, sz, dims, iz, l(3), u(3)
  real, dimension(2)   :: aver, averx, avery
  real                 :: di_min, di_max, do_min, do_max, rms, diff
  !-----------------------------------------------------------------------------
  if (self%state=='write') then
    write (io_unit%validate) patch%id, 1, label1
    write (io_unit%validate) shape(f1)
    write (io_unit%validate) f1
  else if (self%state=='compare') then
    read (io_unit%validate) id0, sz, label0
    read (io_unit%validate) dims
    allocate (    f0(dims))
    allocate (inside(dims))
    allocate ( diff1(dims))
    read (io_unit%validate) f0
    !---------------------------------------------------------------------------
    ! Compute averages and differences
    !---------------------------------------------------------------------------
    l = patch%mesh%li
    u = patch%mesh%ui
    aver = 0.0
    do iz=1,dims
      aver(1) = aver(1) + f0(iz)
      aver(2) = aver(2) + f1(iz)
      diff1(iz) = f1(iz)-f0(iz)
      inside(iz) = (iz >= l(3) .and. iz <= u(3))
    end do
    aver = aver/dims
    diff  = aver(2)-aver(1)
    rms   = sqrt(sum(diff1)/dims)
    di_min = minval(diff1,mask=inside)
    di_max = maxval(diff1,mask=inside)
    do_min = minval(diff1,mask=.not.inside)
    do_max = maxval(diff1,mask=.not.inside)
    deallocate (f0, inside, diff1)
    !---------------------------------------------------------------------------
    ! Write out summary of differences
    !---------------------------------------------------------------------------
    write (io_unit%output,1) '    patch ID =',      id0, trim(label0)
    write (io_unit%output,1) '    patch ID =', patch%id, trim(label1)
    write (io_unit%output,2) '       aver1 =', aver(1), 'aver2 =', aver(2)
    write (io_unit%output,2) '        diff =',    diff, '  RMS =', rms
    write (io_unit%output,2) ' inside: min =',  di_min, '  max =', di_max
    write (io_unit%output,2) 'outside: min =',  do_min, '  max =', do_max
    1 format(a,i12,3x,a)
    2 format(a,1pe12.3,3x,a,1pe12.3)
  end if
END SUBROUTINE dump1

!> =============================================================================
!> Switch
!> =============================================================================
SUBROUTINE dump2 (self, patch, f1, label1)
  class(validate_t)    :: self
  class(patch_t)       :: patch
  integer              :: id1
  real                 :: f1(:,:)
  character(len=*)     :: label1
  !.............................................................................
  real, allocatable    :: f0(:,:), diff1(:), diff2(:,:)
  logical, allocatable :: inside(:,:)
  character(len=64)    :: label0
  integer              :: id0, sz, dims(2), ix, iy, l(3), u(3)
  real, dimension(2)   :: aver, averx, avery
  real                 :: di_min, di_max, do_min, do_max, rms, diff
  !-----------------------------------------------------------------------------
  call trace%begin ('validate_t%dump3')
  if (self%state=='write') then
    label0 = label1
    write (io_unit%validate) patch%id, 2, label0
    write (io_unit%validate) shape(f1)
    write (io_unit%validate) f1
  else if (self%state=='compare') then
    read (io_unit%validate) id0, sz, label0
    read (io_unit%validate) dims
    allocate (    f0(dims(1),dims(2)))
    allocate (inside(dims(1),dims(2)))
    allocate ( diff2(dims(1),dims(2)))
    allocate ( diff1(dims(2)))
    read (io_unit%validate) f0
    !---------------------------------------------------------------------------
    ! Compute averages and differences
    !---------------------------------------------------------------------------
    l = patch%mesh%li
    u = patch%mesh%ui
    aver = 0.0
    do iy=1,dims(2)
      averx = 0.0
      do ix=1,dims(1)
        averx(1) = averx(1) + f0(ix,iy)
        averx(2) = averx(2) + f1(ix,iy)
        diff2(ix,iy) = f1(ix,iy)-f0(ix,iy)
        inside(ix,iy) = (ix >= l(1) .and. ix <= u(1)) .and. &
                        (iy >= l(2) .and. iy <= u(2))
      end do
      aver = aver + averx
    end do
    aver = aver/product(dims)
    diff  = aver(2)-aver(1)
    diff1 = sum(diff2**2,1)
    rms   = sqrt(sum(diff1)/product(dims))
    di_min = minval(diff2,mask=inside)
    di_max = maxval(diff2,mask=inside)
    do_min = minval(diff2,mask=.not.inside)
    do_max = maxval(diff2,mask=.not.inside)
    deallocate (f0, inside, diff1, diff2)
    !---------------------------------------------------------------------------
    ! Write out summary of differences
    !---------------------------------------------------------------------------
    write (io_unit%output,1) '    patch ID =',      id0, trim(label0)
    write (io_unit%output,1) '    patch ID =', patch%id, trim(label1)
    write (io_unit%output,2) '       aver1 =', aver(1), 'aver2 =', aver(2)
    write (io_unit%output,2) '        diff =',    diff, '  RMS =', rms
    write (io_unit%output,2) ' inside: min =',  di_min, '  max =', di_max
    write (io_unit%output,2) 'outside: min =',  do_min, '  max =', do_max
    1 format(a,i12,3x,a)
    2 format(a,1pe12.3,3x,a,1pe12.3)
  end if
  call trace%end
END SUBROUTINE dump2

!> =============================================================================
!> 3-D comparison
!> =============================================================================
SUBROUTINE dump3 (self, patch, f1, label1)
  class(validate_t)    :: self
  class(patch_t)       :: patch
  integer              :: id1
  real                 :: f1(:,:,:)
  character(len=*)     :: label1
  !.............................................................................
  real, allocatable    :: f0(:,:,:), diff1(:), diff2(:,:), diff3(:,:,:)
  logical, allocatable :: inside(:,:,:)
  character(len=64)    :: label0
  integer              :: id0, sz, dims(3), ix, iy, iz, l(3), u(3)
  real, dimension(2)   :: aver, averx, avery
  real                 :: di_min, di_max, do_min, do_max, rms, diff
  !-----------------------------------------------------------------------------
  if (self%state=='write') then
    label0 = label1
    write (io_unit%validate) patch%id, 3, label0
    write (io_unit%validate) shape(f1)
    write (io_unit%validate) f1
  else if (self%state=='compare') then
    read (io_unit%validate) id0, sz, label0
    read (io_unit%validate) dims
    allocate (    f0(dims(1),dims(2),dims(3)))
    allocate (inside(dims(1),dims(2),dims(3)))
    allocate ( diff3(dims(1),dims(2),dims(3)))
    allocate ( diff2(dims(2),dims(3)))
    allocate ( diff1(dims(3)))
    read (io_unit%validate) f0
    !---------------------------------------------------------------------------
    ! Compute averages and differences
    !---------------------------------------------------------------------------
    l = patch%mesh%li
    u = patch%mesh%ui
    aver = 0.0
    do iz=1,dims(3)
      avery = 0.0
      do iy=1,dims(2)
        averx = 0.0
        do ix=1,dims(1)
          averx(1) = averx(1) + f0(ix,iy,iz)
          averx(2) = averx(2) + f1(ix,iy,iz)
          diff3(ix,iy,iz) = f1(ix,iy,iz)-f0(ix,iy,iz)
          inside(ix,iy,iz) = (ix >= l(1) .and. ix <= u(1)) .and. &
                             (iy >= l(2) .and. iy <= u(2)) .and. &
                             (iz >= l(3) .and. iz <= u(3))
        end do
        avery = avery + averx
      end do
      aver = aver + avery
    end do
    aver = aver/product(dims)
    diff2 = sum(diff3**2,1)
    diff1 = sum(diff2,1)
    diff  = aver(2)-aver(1)
    rms   = sqrt(sum(diff1)/product(dims))
    di_min = minval(diff3,mask=inside)
    di_max = maxval(diff3,mask=inside)
    do_min = minval(diff3,mask=.not.inside)
    do_max = maxval(diff3,mask=.not.inside)
    deallocate (f0, inside, diff1, diff2, diff3)
    !---------------------------------------------------------------------------
    ! Write out summary of differences
    !---------------------------------------------------------------------------
    if (rms > 0.0) then
      write (io_unit%output,1) '    patch ID =',      id0, trim(label0)
      write (io_unit%output,1) '    patch ID =', patch%id, trim(label1)
      write (io_unit%output,2) '       aver1 =', aver(1), 'aver2 =', aver(2)
      write (io_unit%output,2) '        diff =',    diff, '  RMS =', rms
      write (io_unit%output,2) ' inside: min =',  di_min, '  max =', di_max
      write (io_unit%output,2) 'outside: min =',  do_min, '  max =', do_max
      self%ok = .false.
    end if
    1 format(a,i12,3x,a)
    2 format(a,1pe12.3,3x,a,1pe12.3)
  end if
END SUBROUTINE dump3

END MODULE validate_mod
