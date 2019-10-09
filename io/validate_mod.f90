!===============================================================================
!> Generic validation module. The general idea is to be able to compare two runs
!> at critical points in the sequence of evaluation, by having the first run 
!> write to files, which the 2nd job reads from, and then compare the contents.
!>
!> For now, the data are stored in a single file per rank, which implies that
!> comparisons only work as long as the updates occur in exactly the same order,
!> allowing
!>
!> 1. Runs that use a single core, and which are therefore updating patches in
!>    exactly the same order, allowing bit-wise match for unchanged code.
!>
!> However, since this functionality is primarily intended for testing and
!> validation of small runs, having one file per patch would not be a problem.
!> This would allow also the kinds of tests below:
!>
!> 2. Runs that use several core on a single MPI rank, where results are not
!>    expected to be bit-wise identical, but where one can get useful information
!>    about the magnitudes of differences, separately in guard zones and the
!>    interior.
!> 3. The same, but withe several ranks.  Assuming patch IDs are identical in
!>    two runs, the same type of quantitative comparisons as in item 2 can be
!>    made.
!>
!> Syntax:
!>   call validate%check (patch%link, scalar, 'label string')
!>   call validate%check (patch%link,  array, 'label string')
!>
!> where 'array' may be 1-D, 2-D, 3-D, or 4-D.
!===============================================================================
MODULE validate_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE patch_mod
  USE link_mod
  USE kinds_mod
  implicit none
  private
  type, public:: validate_t
    integer:: id
    integer:: irec=0
    integer:: iout=0
    integer:: verbose=0
    logical:: ok=.true.
    character(len=8):: mode='off'
  contains
    procedure:: init
    procedure:: set
    procedure, private:: check0
    procedure, private:: check1
    procedure, private:: check2
    procedure, private:: check3
    procedure, private:: check4
    procedure, private:: check0t
    procedure, private:: check1t
    procedure, private:: check2t
    procedure, private:: check3t
    procedure, private:: check4t
    generic, public:: check => check0 , check1 , check2 , check3 , check4, &
                               check0t, check1t, check2t, check3t, check4t
  end type
  type(validate_t), public:: validate
CONTAINS

!===============================================================================
!> Initialize the dump & compare process
!===============================================================================
SUBROUTINE init (self)
  class(validate_t):: self
  !.............................................................................
  integer:: iostat
  integer, save:: verbose=0
  logical, save:: first_time=.true.
  character(len=8), save:: mode='off'
  namelist /validate_params/ verbose, mode
  character(len=120):: id = &
   '$Id$ io/validate_mod.f90'
  !.............................................................................
  call trace%begin ('validate_t%init')
  call trace%print_id (id)
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, validate_params, iostat=iostat)
    write (io%output, validate_params)
  end if
  !$omp end critical (input_cr)
  self%mode = mode
  self%verbose = verbose
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Change mode of the dump process
!===============================================================================
SUBROUTINE set (self, mode1, verbose)
  class(validate_t):: self
  character(len=*):: mode1
  integer, optional:: verbose
  !.............................................................................
  self%mode = mode1
  if (present(verbose)) then
    self%verbose = verbose
  end if
END SUBROUTINE set

!===============================================================================
!> Switch
!===============================================================================
SUBROUTINE check0 (self, patch, f1, label1)
  class(validate_t)    :: self
  class(patch_t)       :: patch
  real                 :: f1
  character(len=*)     :: label1
  !.............................................................................
  real                 :: f0
  character(len=64)    :: label0
  integer              :: id0, id1, sz, dims
  real                 :: diff
  !-----------------------------------------------------------------------------
  if (self%verbose > 0) &
    write (io_unit%output,'(a,i6,2x,i1,20x,f12.6,2x,a)') &
      ' validate:', patch%id, 1, patch%time, label1
  if (self%mode=='write') then
    id1 = patch%id
    label0 = label1
    write (io_unit%validate) id1, 0, label0
    write (io_unit%validate) 1
    write (io_unit%validate) f1
  else if (self%mode=='compare') then
    read (io_unit%validate) id0, sz, label0
    id1 = patch%id
    if (id0 /= id1) write (io_unit%output,*) 'WARNING: id differs', id0, id1
    read (io_unit%validate) dims
    read (io_unit%validate) f0
    diff  = f0-f1
    !---------------------------------------------------------------------------
    ! Write out summary of differences, if there are any
    !---------------------------------------------------------------------------
    if (diff /= 0.0 .or. self%verbose > 0) then
      write (io_unit%output,1) '    patch ID =',      id0, trim(label0)
      write (io_unit%output,1) '    patch ID =',      id1, trim(label1)
      write (io_unit%output,2) '      value1 =',      f0, 'value2 =', f1
      write (io_unit%output,2) '        diff =',    diff
    end if
    1 format(a,i12,3x,a)
    2 format(a,1pe12.3,2x,a,1pe12.3)
  end if
END SUBROUTINE check0

!===============================================================================
!> Scalar comparison
!===============================================================================
SUBROUTINE check0t (self, link, f1, label1)
  class(validate_t)      :: self
  class(link_t), pointer :: link
  real                   :: f1
  character(len=*)       :: label1
  !.............................................................................
  associate (patch => link%task)
  select type (patch)
  class is (patch_t)
  call check0 (self, patch, f1, label1)
  end select
  end associate
END SUBROUTINE check0t

!===============================================================================
!> Switch
!===============================================================================
SUBROUTINE check1 (self, patch, f1, label1)
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
  if (self%verbose > 0) &
    write (io_unit%output,'(a,i6,2x,i1,2x,i4,14x,f12.6,2x,a)') &
      ' validate:', patch%id, 1, shape(f1), patch%time, label1
  if (self%mode=='write') then
    write (io_unit%validate) patch%id, 1, label1
    write (io_unit%validate) shape(f1)
    write (io_unit%validate) f1
  else if (self%mode=='compare') then
    read (io_unit%validate) id0, sz, label0
    id1 = patch%id
    if (id0 /= id1) write (io_unit%output,*) 'WARNING: id differs', id0, id1
    if (sz  /=   1) write (io_unit%output,*) 'WARNING: sz differs', sz, 1
    read (io_unit%validate) dims
    if (any(dims /= shape(f1))) &
      write (io_unit%output,'(1x,a,2(2x,3i4))') 'WARNING: dims differ', dims, shape(f1)
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
    ! Write out summary of differences, if there are any
    !---------------------------------------------------------------------------
    if (rms > 0.0 .or. self%verbose > 1) then
      write (io_unit%output,1) '    patch ID =',      id0, trim(label0)
      write (io_unit%output,1) '    patch ID =', patch%id, trim(label1)
      write (io_unit%output,2) '       aver1 =', aver(1), 'aver2 =', aver(2)
      write (io_unit%output,2) '        diff =',    diff, '  RMS =', rms
      write (io_unit%output,2) ' inside: min =',  di_min, '  max =', di_max
      write (io_unit%output,2) 'outside: min =',  do_min, '  max =', do_max
    end if
    1 format(a,i12,3x,a)
    2 format(a,1pe12.3,3x,a,1pe12.3)
  end if
END SUBROUTINE check1

!===============================================================================
!> 3-D comparison
!===============================================================================
SUBROUTINE check1t (self, link, f1, label1)
  class(validate_t)      :: self
  class(link_t), pointer :: link
  real                   :: f1(:)
  character(len=*)       :: label1
  !.............................................................................
  associate (patch => link%task)
  select type (patch)
  class is (patch_t)
  call check1 (self, patch, f1, label1)
  end select
  end associate
END SUBROUTINE check1t

!===============================================================================
!> Switch
!===============================================================================
SUBROUTINE check2 (self, patch, f1, label1)
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
  call trace%begin ('validate_t%check2')
  if (self%verbose > 0) &
    write (io_unit%output,'(a,i6,2x,i1,2x,2i4,10x,f12.6,2x,a)') &
      ' validate:', patch%id, 2, shape(f1), patch%time, label1
  if (self%mode=='write') then
    label0 = label1
    write (io_unit%validate) patch%id, 2, label0
    write (io_unit%validate) shape(f1)
    write (io_unit%validate) f1
  else if (self%mode=='compare') then
    read (io_unit%validate) id0, sz, label0
    id1 = patch%id
    if (id0 /= id1) write (io_unit%output,*) 'WARNING: id differs', id0, id1
    if (sz  /=   2) write (io_unit%output,*) 'WARNING: sz differs', sz, 2
    read (io_unit%validate) dims
    if (any(dims /= shape(f1))) &
      write (io_unit%output,'(1x,a,2(2x,3i4))') 'WARNING: dims differ', dims, shape(f1)
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
    ! Write out summary of differences, if there are any
    !---------------------------------------------------------------------------
    if (rms > 0.0 .or. self%verbose > 1) then
      write (io_unit%output,1) '    patch ID =',      id0, trim(label0)
      write (io_unit%output,1) '    patch ID =', patch%id, trim(label1)
      write (io_unit%output,2) '       aver1 =', aver(1), 'aver2 =', aver(2)
      write (io_unit%output,2) '        diff =',    diff, '  RMS =', rms
      write (io_unit%output,2) ' inside: min =',  di_min, '  max =', di_max
      write (io_unit%output,2) 'outside: min =',  do_min, '  max =', do_max
    end if
    1 format(a,i12,3x,a)
    2 format(a,1pe12.3,3x,a,1pe12.3)
  end if
  call trace%end()
END SUBROUTINE check2

!===============================================================================
!> 2-D comparison
!===============================================================================
SUBROUTINE check2t (self, link, f1, label1)
  class(validate_t)      :: self
  class(link_t), pointer :: link
  real                   :: f1(:,:)
  character(len=*)       :: label1
  !.............................................................................
  associate (patch => link%task)
  select type (patch)
  class is (patch_t)
  call check2 (self, patch, f1, label1)
  end select
  end associate
END SUBROUTINE check2t

!===============================================================================
!> 3-D comparison
!===============================================================================
SUBROUTINE check3 (self, patch, f1, label1)
  class(validate_t)        :: self
  class(patch_t)           :: patch
  real(kind=KindScalarVar) :: f1(:,:,:)
  character(len=*)         :: label1
  !.............................................................................
  real, allocatable    :: f0(:,:,:), diff1(:), diff2(:,:), diff3(:,:,:)
  logical, allocatable :: inside(:,:,:)
  character(len=64)    :: label0
  integer              :: id0, id1, sz, dims(3), ix, iy, iz, l(3), u(3)
  real, dimension(2)   :: aver, averx, avery
  real                 :: di_min, di_max, do_min, do_max, rms, diff
  !-----------------------------------------------------------------------------
  call trace%begin ('validate_t%check3')
  if (self%verbose > 0) &
    write (io_unit%output,'(a,i6,2x,i1,2x,3i4,6x,f12.6,2x,a)') &
      ' validate:', patch%id, 3, shape(f1), patch%time, label1
  if (self%mode=='write') then
    label0 = label1
    write (io_unit%validate) patch%id, 3, label0
    write (io_unit%validate) shape(f1)
    write (io_unit%validate) f1
  else if (self%mode=='compare') then
    read (io_unit%validate) id0, sz, label0
    id1 = patch%id
    if (id0 /= id1) write (io_unit%output,*) 'WARNING: id differs', id0, id1
    if (sz  /=   3) write (io_unit%output,*) 'WARNING: sz differs', sz, 3
    read (io_unit%validate) dims
    if (any(dims /= shape(f1))) &
      write (io_unit%output,'(1x,a,2(2x,3i4))') 'WARNING: dims differ', dims, shape(f1)
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
    ! Write out summary of differences, if there are any
    !---------------------------------------------------------------------------
    if (rms > 0.0 .or. self%verbose > 1) then
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
  call trace%end()
END SUBROUTINE check3

!===============================================================================
!> 3-D comparison
!===============================================================================
SUBROUTINE check3t (self, link, f1, label1)
  class(validate_t)        :: self
  class(link_t), pointer   :: link
  real(kind=KindScalarVar) :: f1(:,:,:)
  character(len=*)         :: label1
  !.............................................................................
  associate (patch => link%task)
  select type (patch)
  class is (patch_t)
  call check3 (self, patch, f1, label1)
  end select
  end associate
END SUBROUTINE check3t

!===============================================================================
!> 3-D comparison
!===============================================================================
SUBROUTINE check4 (self, patch, f1, label1)
  class(validate_t)    :: self
  class(patch_t)       :: patch
  real                 :: f1(:,:,:,:)
  character(len=*)     :: label1
  !.............................................................................
  real, allocatable    :: f0(:,:,:,:), diff1(:), diff2(:,:), diff3(:,:,:)
  logical, allocatable :: inside(:,:,:)
  character(len=64)    :: label0
  integer              :: id0, id1, sz, dims(4), ix, iy, iz, iv, l(3), u(3)
  real, dimension(2)   :: aver, averx, avery
  real                 :: di_min, di_max, do_min, do_max, rms, diff
  !-----------------------------------------------------------------------------
  call trace%begin ('validate_t%check4')
  if (self%verbose > 0) &
    write (io_unit%output,'(a,i6,2x,i1,2x,4i4,2x,f12.6,2x,a)') &
      ' validate:', patch%id, 4, shape(f1), patch%time, label1
  if (self%mode=='write') then
    label0 = label1
    write (io_unit%validate) patch%id, 4, label0
    write (io_unit%validate) shape(f1)
    write (io_unit%validate) f1
  else if (self%mode=='compare') then
    read (io_unit%validate) id0, sz, label0
    id1 = patch%id
    if (id0 /= id1) write (io_unit%output,*) 'WARNING: id differs', id0, id1
    if (sz  /=   4) write (io_unit%output,*) 'WARNING: sz differs', sz, 4
    read (io_unit%validate) dims
    if (any(dims /= shape(f1))) &
      write (io_unit%output,'(1x,a,2(2x,4i4))') 'WARNING: dims differ', dims, shape(f1)
    allocate     (f0(dims(1),dims(2),dims(3),dims(4)))
    allocate (inside(dims(1),dims(2),dims(3)))
    allocate  (diff3(dims(1),dims(2),dims(3)))
    allocate  (diff2(dims(2),dims(3)))
    allocate  (diff1(dims(3)))
    read (io_unit%validate) f0
    !---------------------------------------------------------------------------
    ! Compute averages and differences
    !---------------------------------------------------------------------------
    l = patch%mesh%li
    u = patch%mesh%ui
    do iv=1,dims(4)
      aver = 0.0
      do iz=1,dims(3)
        avery = 0.0
        do iy=1,dims(2)
          averx = 0.0
          do ix=1,dims(1)
            averx(1) = averx(1) + f0(ix,iy,iz,iv)
            averx(2) = averx(2) + f1(ix,iy,iz,iv)
            diff3(ix,iy,iz) = f1(ix,iy,iz,iv)-f0(ix,iy,iz,iv)
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
      !---------------------------------------------------------------------------
      ! Write out summary of differences, if there are any
      !---------------------------------------------------------------------------
      if (rms > 0.0 .or. self%verbose > 1) then
        write (io_unit%output,1) '    patch ID =',      id0, iv, trim(label0)
        write (io_unit%output,1) '    patch ID =', patch%id, iv, trim(label1)
        write (io_unit%output,2) '       aver1 =', aver(1), 'aver2 =', aver(2)
        write (io_unit%output,2) '        diff =',    diff, '  RMS =', rms
        write (io_unit%output,2) ' inside: min =',  di_min, '  max =', di_max
        write (io_unit%output,2) 'outside: min =',  do_min, '  max =', do_max
        self%ok = .false.
      end if
    end do
    deallocate (f0, inside, diff1, diff2, diff3)
    1 format(a,i12,i4,3x,a)
    2 format(a,1pe12.3,3x,a,1pe12.3)
  end if
  call trace%end()
END SUBROUTINE check4

!===============================================================================
!> 3-D comparison
!===============================================================================
SUBROUTINE check4t (self, link, f1, label1)
  class(validate_t)      :: self
  class(link_t), pointer :: link
  real                   :: f1(:,:,:,:)
  character(len=*)       :: label1
  !.............................................................................
  associate (patch => link%task)
  select type (patch)
  class is (patch_t)
  call check4 (self, patch, f1, label1)
  end select
  end associate
END SUBROUTINE check4t

END MODULE validate_mod
