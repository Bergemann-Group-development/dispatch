!===============================================================================
!> Use MPI parallel I/O to write everything to a single file.  No critical regions
!> should be needed here; all thread protection is done inside mpi_file_mod and
!> mpi_io_mod.
!>
!> There are two pairs of procedures; output/input and output_single/input_single.
!> The latter group supports writing each snapshot in a separate file, but only
!> for MPI libraries that support multiple thread, simultenous MPI calls.
!>
!> In each of the pairs, the io%format value chooses between an older storage
!> pattern, where the position of a patch in the disk image is always the same,
!> and independent of MPI configuration, and one where each rank writes a
!> contiguous part in the file, allowing a much faster output, and also a much
!> faster reading from IDL and Python.
!>
!> A compromise storage pattern would be one where the data is both stored in an
!> MPI-independent way, but still can be read in large chunks. A good arrangement
!> would have all patches in a variable in a contiguous disk space, arranged in
!> an order that does not depend on the MPI-configuration, nor depends on the
!> order of writing (there should be no need to get the order of writes from the
!> order of patches in the patch_rrrrr.nml files, for example).
!>
!> In the current task-ID scheme, task IDs are guaranteed to be unique, without
!> requiring MPI-communication when getting a new ID (which would have to be
!> arranged via MPI-RMA -- with indications of not being reliable). Hence, in
!> the current scheme, task-IDs are already MPI-dependent, and needs to be
!> regenerated when restarting with a different MPI-configuration.  It would
!> thus be best if the position of a patch in the file only depende on the
!> position of the patch in space, and in such a way that a consequence was to
!> have most (if not all) patches written from a single rank written to the
!> same part of the file.
!>
!> A better approach, much easier to implement, is that the reader uses the
!> meta-data saved with the snapshot to deduce the arrangement in the file,
!> and reads the data based on that.  It is then free to choose the pattern to
!> write data with based only on efficiency of writing, and speed of reading.
!===============================================================================
MODULE parallel_io_mod
  USE patch_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_coords_mod
  USE mpi_io_mod
  USE mpi_file_mod
  USE trace_mod
  USE time_slices_mod
  USE bits_mod
  USE kinds_mod
  implicit none
  private
  type, public:: parallel_io_t
    real, dimension(:,:,:,:,:,:), pointer:: buffer => null()
  contains
    procedure, nopass:: init
    procedure, nopass:: input
    procedure, nopass:: output
    procedure, nopass:: close
    procedure, nopass:: open_and_close
  end type
  integer, save:: verbose=0
  logical:: check_files=.false.
  type(mpi_file_t), pointer:: file_in=>null(), file_out=>null()
  type(mpi_file_t), pointer, dimension(:):: files_out=>null()
  type(parallel_io_t), public:: parallel_io
  integer, save:: io0, io1
  integer, save:: opened=-1, closed=-1
  integer, save:: id_debug=0
CONTAINS

!===============================================================================
!> Initialize for MPI parallel I/O.  This is an MPI collective operation, and
!> needs to be done at the start of the run.
!===============================================================================
SUBROUTINE init (patch)
  class(patch_t):: patch
  !.............................................................................
  character(len=5):: snapname
  character(len=120):: filename
  logical:: exists
  integer:: i, iostat
  logical:: first_time=.true.
  namelist /parallel_io_params/ verbose, id_debug
  character(len=120):: id = &
  '$Id$ io/parallel_io_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%print_id (id)
  call trace%begin ('parallel_io_t%init')
  call io%lock%set ('parallel_io_t%init')
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read(io%input, parallel_io_params, iostat=iostat)
    call io%namelist_warning ('parallel_io_params', iostat)
    write (io%namelist, parallel_io_params)
    if (io%ntask==0) &
      call io%abort ('parallel_io_t%init: io%ntask == 0')
    call mpi_io%init
    call mpi_io%set (nwrite=io%nwrite*(io%time_derivs+1))
    if (io%method=='parallel') then
      allocate (file_out)
      write (filename,'(a,"/snapshots.dat")') trim(io%outputname)
      call file_out%openw (filename)
    end if
  end if
  call io%lock%unset ('parallel_io_t%init')
  call time_slices%init (patch%nt)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Close the data I/O method
!===============================================================================
SUBROUTINE close
  integer:: i
  !.............................................................................
  call trace%begin ('parallel_io_t%close')
  if (io%method == 'snapshot') then
    do i=io0,io1
      call files_out(i)%close
    end do
  else
    call file_out%close
  end if
  call file_in%close
  call trace%end()
END SUBROUTINE close

!===============================================================================
!> Adjust the open / closed status of the data I/O files
!===============================================================================
SUBROUTINE open_and_close ()
  integer:: i
  !-----------------------------------------------------------------------------
  ! Only check files when a request to open or close has just been made. Only
  ! one thread should do the checking.
  !-----------------------------------------------------------------------------
  if (check_files) then
    !$omp critical (files_out_cr)
    if (check_files) then
      !$omp atomic write
      check_files = .false.
      !$omp end atomic
      call trace%begin ('parallel_io_t%open_and_close')
      if (associated(files_out)) then
        do i=io0,io1
          call files_out(i)%open_and_close()
        end do
      end if
      call trace%end()
    end if
    !$omp end critical (files_out_cr)
    if (verbose > 1) then
      write (io_unit%log,*) wallclock(), 'parallel_io_t%open_and_close: exiting'
      flush (io_unit%log)
    end if
  end if
END SUBROUTINE open_and_close

!===============================================================================
!> Write output in increments of one patch data cube
!===============================================================================
SUBROUTINE output (patch, count)
  class(patch_t):: patch
  integer:: count
  !.............................................................................
  integer(8):: snapshot_words, buffer_words, pos
  integer:: iv, n(3), jt(2), iw, nw
  real(kind=KindScalarVar), dimension(:,:,:,:,:), allocatable, target:: mem
  real(kind=4),             dimension(:,:,:,:)  , allocatable, target:: buffer
  character(len=120):: filename
  real:: pt(2)
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output', itimer=itimer)
  if (io%format /= 14) then
    io%format = 15
  end if
  !-----------------------------------------------------------------------------
  ! Initialize file_out
  !-----------------------------------------------------------------------------
  if (.not.associated(file_out)) then
    write (filename,'(a,"/snapshots.dat")') trim(io%outputname)
    allocate (file_out)
    call file_out%openw (filename)
  end if
  !-----------------------------------------------------------------------------
  ! Set size parameters; a buffer is all variables in the patch, and the
  ! snspshot offset is that times the number of patches per snapshot.
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    n = patch%gn
  else if (patch%no_mans_land) then
    n = patch%n
  else
    n = patch%n-1
  end if
  !-----------------------------------------------------------------------------
  ! io%format=14 then slots in the file are ordered (id,iv,iout), while for 15
  ! they are ordered (iv,id,iout).  The former is much faster to read, while the
  ! latter is faster to write.
  !-----------------------------------------------------------------------------
  snapshot_words = product(n)*io%nv*io%ntotal
  if (io%format==14) then
    buffer_words = product(n)
  else if (io%format==15) then
    buffer_words = product(n)*io%nv
  end if
  if (io%method == 'snapshot') then
    call mpi_io%use (files_out(patch%iout))
  else
    call mpi_io%use (file_out)
  end if
  call mpi_io%set (snapshot_words, buffer_words, io%nwrite*(io%time_derivs+1))
  !-----------------------------------------------------------------------------
  ! Info
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    write(io%output,'(a,f8.3,a,i10,a)') &
      ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
  end if
  if (verbose > 0) then
    write(io%output,'(1x,a,6i6,2i12,l3)') &
      'parallel_io_t%output: id, ip, iout, n, snapshot_words, buffer_words =', &
      patch%id, patch%ip, patch%iout, n, snapshot_words, buffer_words, &
      associated(parallel_io%buffer)
    flush (io%output)
  end if
  !-----------------------------------------------------------------------------
  ! Time interpolation to exact out_next, by first transferring the values in
  ! patch%mem to the mem(:,:,:,:) temporary, resorting into increasing time
  ! order and optionally also converting to a different representation.
  !-----------------------------------------------------------------------------
  if (verbose > 0) then
    write (io_unit%log,'(a,i6,2i4)') 'parallel_io_t%output: id, nv, io%nv =', &
      patch%id, patch%nv, io%nv
  end if
  allocate (mem(n(1),n(2),n(3),patch%nt-1,io%nv))
  allocate (buffer(n(1),n(2),n(3),io%nv))
  call prepare_time_slices
  do iv=1,io%nv
    call time_slices%interpolate (patch, mem(:,:,:,:,iv), buffer(:,:,:,iv))
    if (patch%id==id_debug) then
      write (io_unit%output,'(1x,a,i6,i5,i4,1p,g14.6,3(2x,2e11.2))') &
        'parallel_io_t%output: id, io, iv, time, minval, maxval =', &
        patch%id, patch%iout, iv, patch%time, &
        minval(patch%mem(:,:,:,iv,patch%it,1)), &
        maxval(patch%mem(:,:,:,iv,patch%it,1)), &
        minval(buffer(:,:,:,iv)), maxval(buffer(:,:,:,iv)), &
        minval(mem(:,:,:,patch%nt-1,iv)), maxval(mem(:,:,:,patch%nt-1,iv))
    end if
  end do
  call io_write ((io%time_derivs+1)*patch%iout, 0)
  !-----------------------------------------------------------------------------
  ! On-the-fly time derivatives output
  !-----------------------------------------------------------------------------
  if (io%time_derivs>0) then
    do iv=1,io%nv
      call time_slices%derivative1 (patch, mem(:,:,:,:,iv), buffer(:,:,:,iv))
    end do
    call io_write ((io%time_derivs+1)*patch%iout, 1)
    if (io%time_derivs==2) then
      do iv=1,io%nv
        call time_slices%derivative2 (patch, mem(:,:,:,:,iv), buffer(:,:,:,iv))
      end do
      call io_write ((io%time_derivs+1)*patch%iout, 2)
    end if
  end if
  if (count == 0) then
    call mpi_io%sync ()
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate mem, and buffer.  If io_write uses non-blocking MPI_File_iwrite
  ! calls, the buffer must remain allocated, until a later call to check_io
  ! finds that it is no longer needed, and deallocates it
  !-----------------------------------------------------------------------------
  deallocate (mem, buffer)
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  ! Copy variables from `patch%mem` to `mem` array, converting to unsigned (log)
  ! or pervolume representation, if requested
  !-----------------------------------------------------------------------------
  subroutine prepare_time_slices
    integer:: l1(3), u1(3), l2(3), u2(3), n(3), it, iit
    real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e, p1, p2, p3, phi
    real(kind=KindScalarVar), dimension(:,:,:), pointer:: b1d, b2d, b3d, b1u, b2u, b3u
    !...........................................................................
    call trace%begin ('parallel_io-t%output::prepare_time_slices')
    if (io%guard_zones) then
      l1 = patch%mesh%lb
      u1 = max(patch%mesh%ub,patch%mesh%gn)
    else if (patch%no_mans_land) then
      l1 = patch%mesh%li
      u1 = patch%mesh%ui
    else
      l1 = patch%mesh%li
      u1 = patch%mesh%ui-1
    end if
    n = u1-l1+1
    l2 = l1+1
    u2 = u1+1
    associate (idx=>patch%idx)
    if (verbose > 0) then
      write (io_unit%log,'(a,i6,4i4)') &
       'parallel_io_t%prepare_time_slices: nv, nw, bx, phi =', &
        patch%id, patch%nv, patch%nw, patch%idx%bx, patch%idx%bx
    end if
    do it=1,patch%nt-1
      iit = patch%iit(it)
      d  => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%d ,iit,1)
      e  => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%e ,iit,1)
      p1 => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%px,iit,1)
      p2 => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%py,iit,1)
      p3 => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%pz,iit,1)
      mem(:,:,:,it,idx%d ) = log(d)
      !-------------------------------------------------------------------------
      ! For RAMSES cases, convert to velocity and energy per unit mass,
      ! excluding kinetic energy
      !-------------------------------------------------------------------------
      if (patch%kind(1:6) == 'ramses') then
        mem(:,:,:,it,idx%px) = p1/d
        mem(:,:,:,it,idx%py) = p2/d
        mem(:,:,:,it,idx%pz) = p3/d
        mem(:,:,:,it,idx%e ) = (e-0.5*(p1**2+p2**2+p3**2)/d)/d
      else
        mem(:,:,:,it,idx%px) = p1
        mem(:,:,:,it,idx%py) = p2
        mem(:,:,:,it,idx%pz) = p3
        mem(:,:,:,it,idx%e ) = e/d
      end if
      !-------------------------------------------------------------------------
      ! MHD variables and potential, if present
      !-------------------------------------------------------------------------
      if (patch%idx%bx > 0) then
        mem(:,:,:,it,idx%e) = patch%ee_min
        b1d => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%bx,iit,1)
        b2d => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%by,iit,1)
        b3d => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%bz,iit,1)
        mem(:,:,:,it,idx%bx) = b1d
        mem(:,:,:,it,idx%by) = b2d
        mem(:,:,:,it,idx%bz) = b3d
      end if
      if (patch%idx%phi > 0) then
        phi => patch%mem(l1(1):u1(1),l1(2):u1(2),l1(3):u1(3),idx%phi,iit,1)
        mem(:,:,:,it,idx%phi) = phi
      end if
    end do
    end associate
    call trace%end()
  end subroutine prepare_time_slices
  !-----------------------------------------------------------------------------
  subroutine io_write (slot, ider)
    integer:: slot, ider, irec
    if (io%format==14) then
      !-------------------------------------------------------------------------
      ! Ordering with (ip,iv,ider,io)
      !-------------------------------------------------------------------------
      do iv=1,io%nv
        irec = patch%id + io%ntotal*((iv-1) + io%nv*ider)
        call mpi_io%write (buffer(:,:,:,iv), 1+slot, irec)
      end do
    else
      !-------------------------------------------------------------------------
      ! This ordering has the derivative order added to the patch%id, which is
      ! in turn causing an increment by (io%time_derivs+1), so effectively adding
      ! nv additional variables for each derivative order, making it easy to
      ! address the derivatives, while using snapshot numbers incremented by 1.
      !-------------------------------------------------------------------------
      call mpi_io%write (buffer, 1+slot, 1+ider+(patch%id-1)*(io%time_derivs+1))
    end if
  end subroutine
END SUBROUTINE output

!===============================================================================
!> Read input in increments of one variable data cube
!===============================================================================
SUBROUTINE input (patch, ok)
  class(patch_t):: patch
  logical:: ok
  !.............................................................................
  integer(8):: snapshot_words, buffer_words
  integer:: iv, n(3), slot
  real, dimension(:,:,:,:), pointer:: buffer
  character(len=120):: filename
  logical:: single_patch
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%input', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If the patch is a virtual patch, signal that it should be filled with IC
  ! data at first -- it will become updated via MPI later. By setting time < 0
  ! we allow restarting from the IC snapshot, for tests
  !-----------------------------------------------------------------------------
  if (patch%is_set(bits%virtual)) then
    ok = .false.
    patch%time = -1d0
    call trace%end (itimer)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Open the appropriate input file
  !-----------------------------------------------------------------------------
  if (io%method == 'snapshot') then
    write (filename,'(a,"/",i5.5,"/snapshot.dat")') &
      trim(io%outputname), io%restart
  else
    filename=trim(io%inputdir)//'snapshots.dat'
  end if
  allocate (file_in)
  call file_in%openr (filename)
  !-----------------------------------------------------------------------------
  ! Choose input reader
  !-----------------------------------------------------------------------------
  !if (io%format > 9 .and. io%format < 14) then
    !call input_single (patch, ok)
    !call trace%end (itimer)
    !return
  !end if
  !-----------------------------------------------------------------------------
  ! Set size parameters
  !-----------------------------------------------------------------------------
  if (file_in%handle==0) then
    call init (patch)
  end if
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  !-----------------------------------------------------------------------------
  ! Catch reading from a snapshot with different mpi_size
  !-----------------------------------------------------------------------------
  single_patch = any(io%mpi_odims/=io%mpi_dims)
  !-----------------------------------------------------------------------------
  ! Compute offset size constants for the mpi_io reader
  !-----------------------------------------------------------------------------
  snapshot_words = product(n)*io%ntotal*io%nv*(io%time_derivs+1)
  if (single_patch) then
    buffer_words = product(n)
  else
    buffer_words = product(n)*io%ntask
  end if
  !-----------------------------------------------------------------------------
  ! Info
  !-----------------------------------------------------------------------------
  !$omp critical (first_cr)
  if (first_time) then
    first_time = .false.
    write(io%output,'(a,f8.3,a,i10,a)') &
      ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
    write(io%output,'(a,4i6,2i12)') &
      ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
      patch%iout, n, snapshot_words, buffer_words
  end if
  !$omp end critical (first_cr)
  !-----------------------------------------------------------------------------
  call mpi_io%use (file_in)
  call mpi_io%set (snapshot_words, buffer_words, io%nwrite*(io%time_derivs+1))
  !-----------------------------------------------------------------------------
  ! Read via buffer, in case guard zones not included. FIXME: read snapshot
  ! guard zone parameter, to allow it to differ
  !-----------------------------------------------------------------------------
  allocate (buffer(n(1),n(2),n(3),io%nv))
  patch%iout = io%restart
  slot = (io%time_derivs+1)*patch%iout
  call mpi_io%read (buffer, 1+slot, 1+(patch%id-1)*(io%time_derivs+1))
  if (mpi_io%err /= 0) then
    ok = .false.
    call trace%end (itimer)
    return
  end if
  if (io%guard_zones) then
    patch%mem(:,:,:,1:io%nv,patch%it,1) = buffer
  else
    patch%mem(patch%li(1):patch%ui(1), &
              patch%li(2):patch%ui(2), &
              patch%li(3):patch%ui(3),1:io%nv,patch%it,1) = buffer
  end if
  if (verbose > 1) then
    do iv=1,io%nv
      write (io_unit%output,*) 'parallel_io_t%input: id, iv, iout, min, max =', &
        patch%id, iv, patch%iout, patch%fminval(iv), patch%fmaxval(iv)
    end do
  else if (verbose > 0) then
    write (io_unit%output,*) 'parallel_io_t%input: id, iout =', &
      patch%id, patch%iout
  end if
  deallocate (buffer)
  do iv=1,io%nv
    call convert_variables ! copy `patch%mem` into `mem`, converting if desired
  end do
  !-----------------------------------------------------------------------------
  ! Make sure all time slots have the correct time set
  !-----------------------------------------------------------------------------
  patch%t(patch%iit) = patch%time
  !-----------------------------------------------------------------------------
  ! Signal success
  !-----------------------------------------------------------------------------
  ok = .true.
  if (verbose==1) &
    write(io%output,'(a,2i6,1p,2e16.6)') ' parallel_io_t%input: iout, id, time, dtime =', &
    patch%iout, patch%id, patch%time, patch%dtime
  if (verbose>1) &
    write(io%output,'(a,4i6,2i7)') &
    ' parallel_io_t%input: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  ! Copy variables from `patch%mem` to `mem` array, converting variables based
  ! on the I/O format parameter.
  !-----------------------------------------------------------------------------
  subroutine convert_variables
    !-----------------------------------------------------------------------------
    ! Convert log(d) to d
    !-----------------------------------------------------------------------------
    if (iv == patch%idx%d) then
      patch%mem(:,:,:,iv,:,1) = exp(patch%mem(:,:,:,iv,:,1))
    end if
    !-----------------------------------------------------------------------------
    ! Convert velocity to momentum in RAMSES cases
    !-----------------------------------------------------------------------------
    if (patch%kind(1:6) == 'ramses') then
      if (iv >= patch%idx%px .and. iv <= patch%idx%pz) then
        patch%mem(:,:,:,iv,:,1) = patch%mem(:,:,:,iv,:,1)*patch%mem(:,:,:,patch%idx%d,:,1)
      !-----------------------------------------------------------------------------
      ! Convert per-unit-mass energy to per-unit-volume, adding kinetic energy
      !-----------------------------------------------------------------------------
      else if (iv == patch%idx%e) then
        patch%mem(:,:,:,iv,:,1) = patch%mem(:,:,:,patch%idx%d,:,1) * &
                                  patch%mem(:,:,:,patch%idx%e,:,1) + &
                             0.5*(patch%mem(:,:,:,patch%idx%px,:,1)**2 + &
                                  patch%mem(:,:,:,patch%idx%py,:,1)**2 + &
                                  patch%mem(:,:,:,patch%idx%pz,:,1)**2) &
                                 /patch%mem(:,:,:,patch%idx%d ,:,1)
      end if
    end if
  end subroutine
END SUBROUTINE input

END MODULE parallel_io_mod
