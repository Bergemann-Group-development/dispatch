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
  '$Id$'
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
  integer:: iv, n(3), jt(2), iw, nw, l(3), u(3)
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: mem
  real(kind=4),             dimension(:,:,:,:), allocatable, target:: buffer
  character(len=120):: filename
  real:: pt(2)
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Prefer format 14, which is much faster to read
  !-----------------------------------------------------------------------------
  if (io%format /= 15) then
    io%format = 14
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
  if (io%ntotal==0) then
    io%ntotal = io%nwrite*mpi%size
  end if
  snapshot_words = product(n)*io%nv*io%ntotal
  if (io%format==15) then
    buffer_words = product(n)*io%nv
  else
    buffer_words = product(n)
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
    write(io%output,'(1x,a,5i6,2i12,l3)') &
      'parallel_io_t%output: id,iout,n,snapshot_words,buffer_words=', &
      patch%id, patch%iout, n, snapshot_words, buffer_words, &
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
  allocate (buffer(n(1),n(2),n(3),io%nv))
  do iv=1,io%nv
    if (io%guard_zones) then
      mem => patch%mem(:,:,:,iv,:,1)
    else
      l = patch%mesh%li
      u = patch%mesh%ui
      mem => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,:,1)
    end if
    call time_slices%interpolate (patch, mem, buffer(:,:,:,iv))
    if (patch%id==id_debug) then
      write (io_unit%output,'(1x,a,i6,i5,i4,1p,g14.6,3(2x,2e11.2))') &
        'parallel_io_t%output: id, io, iv, time, minval, maxval =', &
        patch%id, patch%iout, iv, patch%time, &
        minval(patch%mem(:,:,:,iv,patch%it,1)), &
        maxval(patch%mem(:,:,:,iv,patch%it,1)), &
        minval(buffer(:,:,:,iv)), maxval(buffer(:,:,:,iv))
    end if
  end do
  call io_write ((io%time_derivs+1)*patch%iout, 0)
  !-----------------------------------------------------------------------------
  ! On-the-fly time derivatives output
  !-----------------------------------------------------------------------------
  if (io%time_derivs>0) then
    do iv=1,io%nv
      if (io%guard_zones) then
        mem => patch%mem(:,:,:,iv,:,1)
      else
        l = patch%mesh%li
        u = patch%mesh%ui
        mem => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,:,1)
      end if
      call time_slices%derivative1 (patch, mem, buffer(:,:,:,iv))
    end do
    call io_write ((io%time_derivs+1)*patch%iout, 1)
    if (io%time_derivs==2) then
      if (io%guard_zones) then
        mem => patch%mem(:,:,:,iv,:,1)
      else
        l = patch%mesh%li
        u = patch%mesh%ui
        mem => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,:,1)
      end if
      do iv=1,io%nv
        call time_slices%derivative2 (patch, mem, buffer(:,:,:,iv))
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
  deallocate (buffer)
  call trace%end (itimer)
contains
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
  integer:: iv, n(3), slot, irec
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
  !if (io%in%format > 9 .and. io%format < 14) then
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
  single_patch = any(io%in%mpi_dims/=io%mpi_dims)
  if (single_patch) then
    write (stderr,'(2(3i4,2x))') io%mpi_dims, io%in%mpi_dims
    call io%abort("changing MPI geometry not supported with method='parallel'")
  end if
  !-----------------------------------------------------------------------------
  ! Compute offset size constants for the mpi_io reader
  !-----------------------------------------------------------------------------
  snapshot_words = product(n)*io%in%ntotal*io%nv*(io%in%time_derivs+1)
  if (io%in%format==15) then
    buffer_words = product(n)*io%nv
  else
    buffer_words = product(n)
  end if
  !-----------------------------------------------------------------------------
  ! Info
  !-----------------------------------------------------------------------------
  !$omp critical (first_cr)
  if (first_time) then
    first_time = .false.
    write(io%output,'(a,f8.3,a,i10,a)') &
      ' input snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
    !write(io%output,'(a,4i6,2i12)') &
      !' parallel_io_t%input: iout, n, snapshot_words, buffer_words =', &
      !patch%iout, n, snapshot_words, buffer_words
  end if
  !$omp end critical (first_cr)
  !-----------------------------------------------------------------------------
  call mpi_io%use (file_in)
  call mpi_io%set (snapshot_words, buffer_words, io%nwrite*(io%in%time_derivs+1))
  !-----------------------------------------------------------------------------
  ! Read via buffer, in case guard zones not included. FIXME: read snapshot
  ! guard zone parameter, to allow it to differ
  !-----------------------------------------------------------------------------
  if (verbose > 0) then
    write(io%output,'(1x,a,5i6,2i12,l3)') &
      'parallel_io_t%input: id,iout,n,snapshot_words,buffer_words=', &
      patch%id, patch%iout, n, snapshot_words, buffer_words
    flush (io%output)
  end if
  !.............................................................................
  allocate (buffer(n(1),n(2),n(3),io%nv))
  patch%iout = io%restart
  slot = (io%in%time_derivs+1)*patch%iout
  !-----------------------------------------------------------------------------
  if (io%in%format==15) then
    !---------------------------------------------------------------------------
    ! This ordering has the derivative order added to the patch%id, which is
    ! in turn causing an increment by (io%time_derivs+1), so effectively adding
    ! nv additional variables for each derivative order, making it easy to
    ! address the derivatives, while using snapshot numbers incremented by 1.
    !---------------------------------------------------------------------------
    call mpi_io%read (buffer, 1+slot, 1+(patch%id-1)*(io%in%time_derivs+1))
  else
    !---------------------------------------------------------------------------
    ! Ordering with (ip,iv,ider,io)
    !---------------------------------------------------------------------------
    do iv=1,io%nv
      irec = patch%id + io%in%ntotal*((iv-1))
      call mpi_io%read (buffer(:,:,:,iv), 1+slot, irec)
    end do
  end if

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
  deallocate (buffer)
  !-----------------------------------------------------------------------------
  ! Make sure all time slots have the correct time set
  !-----------------------------------------------------------------------------
  patch%t(patch%iit) = patch%time
  !-----------------------------------------------------------------------------
  ! Signal success
  !-----------------------------------------------------------------------------
  ok = .true.
  if (verbose == 2) &
    write(io%output,'(a,2i6,1p,2e16.6)') &
    ' parallel_io_t%input: iout, id, time, dtime =', &
    patch%iout, patch%id, patch%time, patch%dtime
  if (verbose > 2) &
    write(io%output,'(a,4i6,2i7)') &
    ' parallel_io_t%input: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call trace%end (itimer)
END SUBROUTINE input

END MODULE parallel_io_mod
