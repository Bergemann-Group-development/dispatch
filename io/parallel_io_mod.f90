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
  namelist /parallel_io_params/ verbose
  character(len=120):: id = &
  '$Id$ io/parallel_io_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%print_id (id)
  call trace%begin ('parallel_io_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input); read(io%input, parallel_io_params, iostat=iostat)
     write (io%output, parallel_io_params)
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
  !$omp end critical (input_cr)
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
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: mem
  real(kind=4), dimension(:,:,:,:), pointer:: buffer
  character(len=120):: filename
  real:: pt(2)
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output', itimer=itimer)
  if (io%format > 9 .and. io%format < 14) then
    call output_single (patch, count)
    call trace%end (itimer)
    return
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
  else
    n = patch%n
  end if
  !-----------------------------------------------------------------------------
  ! If 9 < io%format < 14 then slots in the file are ordered (id,iv,iout),
  ! while otherwise they are ordered (iv,id,iout).  The former is much faster
  ! to read, while the latter is faster to write.
  !-----------------------------------------------------------------------------
  snapshot_words = product(n)*io%nv*io%ntotal
  buffer_words = product(n)*io%nv
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
    write(io%output,'(a,4i6,2i12)') &
      ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
      patch%iout, n, snapshot_words, buffer_words
  end if
  !-----------------------------------------------------------------------------
  ! Time interpolation to exact out_next
  !-----------------------------------------------------------------------------
  allocate (mem(n(1),n(2),n(3),patch%nt))
  allocate (buffer(n(1),n(2),n(3),io%nv))
  do iv=1,io%nv
    call convert_variables ! copy `patch%mem` into `mem`, converting if desired
    call time_slices%interpolate (patch, mem, buffer(:,:,:,iv))
  end do
  call io_write ((io%time_derivs+1)*patch%iout, 0)
  !-----------------------------------------------------------------------------
  ! On-the-fly time derivatives output
  !-----------------------------------------------------------------------------
  if (io%time_derivs>0) then
    do iv=1,io%nv
      call convert_variables ! copy `patch%mem` into `mem`, converting if desired
      call time_slices%derivative1 (patch, mem, buffer(:,:,:,iv))
    end do
    call io_write ((io%time_derivs+1)*patch%iout, 1)
    if (io%time_derivs==2) then
      do iv=1,io%nv
        call convert_variables ! copy `patch%mem` into `mem`, converting if desired
        call time_slices%derivative2 (patch, mem, buffer(:,:,:,iv))
      end do
      call io_write ((io%time_derivs+1)*patch%iout, 2)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate mem, but do NOT deallocate buffer.  Since io_write uses non-
  ! blocking MPI_File_iwrite calls, the buffer must remain allocated, until a
  ! later call to check_io finds that it is no longer needed, and deallocates it
  !-----------------------------------------------------------------------------
  deallocate (mem)
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  ! Copy variables from `patch%mem` to `mem` array, converting variables based
  ! on the I/O format parameter.  This procedure is internal to the 'output'
  ! procedure, and as such only supports io%format in the range 6-9.
  !-----------------------------------------------------------------------------
  subroutine convert_variables
    integer:: nt1, l(3), u(3), i, j, k
    !...........................................................................
    call trace%begin ('parallel_io-t%output::convert_variables')
    if (io%guard_zones) then
      l = patch%mesh%lb
      u = max(patch%mesh%ub,patch%mesh%gn)
    else
      l = patch%mesh%li
      u = patch%mesh%ui
    end if
    nt1 = patch%nt-1
    if (io%format >= 8 .and. iv==patch%idx%d) then
      ! --- Convert density to log density ---
      do k=l(3),u(3); do j=l(2),u(2); do i=l(1),u(1)
        mem(i-l(1)+1,j-l(2)+1,k-l(3)+1,1:nt1) &
        = log(patch%mem(i,j,k,iv,patch%iit(1:nt1),1))
      end do; end do; end do
    else if ((io%format == 6 .or. io%format==8) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert momentum to velocity ---
      do k=l(3),u(3); do j=l(2),u(2); do i=l(1),u(1)
        mem(i-l(1)+1,j-l(2)+1,k-l(3)+1,1:nt1) &
          = patch%mem(i,j,k,iv,patch%iit(1:nt1),1) &
          / patch%mem(i,j,k,patch%idx%d,patch%iit(1:nt1),1)
      end do; end do; end do
    else
      mem = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,patch%iit(1:nt1),1)
    end if
    call trace%end()
  end subroutine
  !-----------------------------------------------------------------------------
  subroutine io_write (slot, ider)
    integer:: slot, ider
    !---------------------------------------------------------------------------
    ! This ordering has the derivative order added to the snapshot number, which
    ! is incremented by (io%time_derivs+1), with patch%id starting with 1, and
    ! ider and slot starting with 0. This make the time derivative results 
    ! appear to be extra snapshots (thus making accessing the snapshot metadata
    ! difficult). Introduced for compatibility with Paolo's analysis scripts.
    !---------------------------------------------------------------------------
    if (io%method == 'pan') then
      call mpi_io%iwrite (buffer, 1+ider+slot, patch%id)
    !---------------------------------------------------------------------------
    ! This ordering has the derivative order added to the patch%id, which is
    ! in turn causing an increment by (io%time_derivs+1), so effectively adding
    ! nv additional variables for each derivative order, making it easy to
    ! address the derivatives, while using snapshot numbers incremented by 1.
    !---------------------------------------------------------------------------
    else
      call mpi_io%iwrite (buffer, 1+slot, 1+ider+(patch%id-1)*(io%time_derivs+1))
    end if
  end subroutine
END SUBROUTINE output

!===============================================================================
!> Write output as one single buffer
!===============================================================================
SUBROUTINE output_single (patch, count)
  class(patch_t):: patch
  integer:: count
  !.............................................................................
  integer(8):: snapshot_words, buffer_words, pos
  integer:: iv, n(3), jt(2), iw, nw, nder, irec, ider, ip, i, iout
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: mem
  real(kind=4), dimension(:,:,:), pointer:: buffer
  real:: pt(2)
  real(8):: start
  logical:: first_time=.true.
  character(len=120):: filename
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output_single', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! method='snapshot' => One file per snapshot
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    write (io_unit%log,*) wallclock(), patch%id, patch%time, '  begin'
    flush (io_unit%log)
  end if
  if (io%method == 'snapshot') then
    if (.not. associated(files_out)) then
      io0 = patch%iout
      io1 = io0 + nint((io%end_time-patch%out_next)/io%out_time) + 1
      closed = io0-1
      if (verbose > 0 .and. omp%master) then
        write (io_unit%log,'(a,i4," -- ",i4)') &
          'parallel_io_t%output: initializing for snapshots', io0, io1
      end if
      allocate (files_out(io0:io1))
      do i=io0,io0+1
        write (filename,'(a,"/",i5.5,"/snapshot.dat")') trim(io%outputname), i
        call files_out(i)%openw (filename)
      end do
    end if
    !---------------------------------------------------------------------------
    ! If new snapshot number, open a new file
    !---------------------------------------------------------------------------
    if (patch%iout+1 > opened) then
      !$omp atomic write
      opened = patch%iout+1
      write (filename,'(a,"/",i5.5,"/snapshot.dat")') trim(io%outputname), opened
      if (verbose > 1) then
        write (io_unit%log,*) wallclock(), patch%id, patch%time, '  request_open'
        flush (io_unit%log)
      end if
      call files_out(opened)%request_open (filename)
      !$omp atomic write
      check_files = .true.
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Set size parameters; a buffer is all variables in the patch, and the
  ! snspshot offset is that times the number of patches per snapshot.
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    write (io_unit%log,*) wallclock(), patch%id, patch%time, '  set size'
    flush (io_unit%log)
  end if
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  !-----------------------------------------------------------------------------
  ! Slots in the file are ordered (id,iv,iout).
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    write (io_unit%log,*) wallclock(), patch%id, patch%time, '  call mpi_io'
    flush (io_unit%log)
  end if
  if (io%ntask==0) &
    call io%abort ('parallel_io_t%output_single: io%ntask == 0')
  !-----------------------------------------------------------------------------
  nder = io%time_derivs+1
  snapshot_words = product(n)*io%ntotal*io%nv*nder
  buffer_words = product(n)*io%ntask
  if (io%method == 'snapshot') then
    call mpi_io%use (files_out(patch%iout))
  else
    call mpi_io%use (file_out)
  end if
  !-----------------------------------------------------------------------------
  ! Info
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    write(io%output,'(a,f8.3,a,i10,a)') &
      ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
    write(io%output,'(a,4i6,2i12)') &
      ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
      patch%iout, n, snapshot_words, buffer_words
  end if
  !-----------------------------------------------------------------------------
  call mpi_io%set (snapshot_words, buffer_words, 99)
  !-----------------------------------------------------------------------------
  if (.not.associated(parallel_io%buffer)) then
    allocate (parallel_io%buffer(n(1),n(2),n(3),io%ntask,io%nv,nder))
    if (verbose > 0) &
      write (io%output,*) 'shape io%buffer =', shape(parallel_io%buffer)
  end if
  ! --- rank-local patch 1D index; cf. cartesian_mod ---
  ip = patch%ip
  if (ip < 0 .or. ip > io%ntask .and. patch%time < io%end_time) then
      write (io_unit%output,*) &
        'parallel_io_t%output_single: WARNING io, ntask =', ip, io%ntask
    ip = max(1,min(ip,io%ntask))
  end if
  if (verbose > 0) then
    write(io%output,'(a,6i6,2i12,l3)') &
      ' parallel_io_t%output_single: id, ip, iout, n, snapshot_words, buffer_words =', &
      patch%id, ip, patch%iout, n, snapshot_words, buffer_words, associated(parallel_io%buffer)
    flush (io%output)
  end if
  if (verbose > 1) then
    write (io_unit%log,*) wallclock(), patch%id, patch%time, '  allocate mem'
    flush (io_unit%log)
  end if
  allocate (mem(n(1),n(2),n(3),patch%nt-1))
  !-----------------------------------------------------------------------------
  ! Time interpolation to exact time = out_next.
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    write (io_unit%log,*) wallclock(), patch%id, patch%time, '  convert var, interpolate'
    flush (io_unit%log)
  end if
  do iv=1,io%nv
    call convert_variables ! copy `patch%mem` into `mem`, converting if desired
    buffer => parallel_io%buffer(:,:,:,ip,iv,1)
    call time_slices%interpolate (patch, mem, buffer)
  end do
  !-----------------------------------------------------------------------------
  ! On-the-fly time derivatives output
  !-----------------------------------------------------------------------------
  if (io%time_derivs>0) then
    do iv=1,io%nv
      call convert_variables ! copy `patch%mem` into `mem`, converting if desired
      buffer => parallel_io%buffer(:,:,:,ip,iv,2)
      call time_slices%derivative1 (patch, mem, buffer)
    end do
    if (io%time_derivs==2) then
      do iv=1,io%nv
        call convert_variables ! copy `patch%mem` into `mem`, converting if desired
        buffer => parallel_io%buffer(:,:,:,ip,iv,3)
        call time_slices%derivative2 (patch, mem, buffer)
      end do
    end if
  end if
  deallocate (mem)
  if (verbose > 0) then
    write (io_unit%log,*) wallclock(), patch%id, patch%time, count, ' buffer'
    flush (io_unit%log)
  end if
  if (count == 0) then
    if (io%method == 'snapshot') then
      iout = 0
      file_out => files_out(patch%iout)
    else
      iout = patch%iout
    end if
    start = wallclock()
    if (verbose > 0) then
      write (io_unit%log,*) wallclock(), patch%id, patch%time, '  start:'//trim(file_out%filename)
      flush (io_unit%log)
    end if
    do ider=1,nder
    do iv=1,io%nv
      !-------------------------------------------------------------------------
      ! Each "record" is product(n)*io%ntasks words long, and contains values
      ! for one variable, for all the patches in a rank.  The ranks are written
      ! out one after the other for one variable, followed by the next var
      !-------------------------------------------------------------------------
      irec = 1 + mpi%rank + mpi%size*(iv-1 + (ider-1)*io%nv)
      if (verbose > 1) &
        write (io_unit%output,*) 'write(1): irec, iout =', irec, iout
      write (io_unit%log,*) 'mpi_io%write: iout, irec =', 1+iout, irec
      call mpi_io%write (parallel_io%buffer(:,:,:,:,iv,ider), 1+iout, irec)
    end do
    end do
    if (verbose > 0) then
      write (io_unit%log,*) wallclock(), patch%id, patch%time, '    end:'//trim(file_out%filename)
      flush (io_unit%log)
    end if
    write (io_unit%log,'(a,f7.3,1x,a)') &
      'parallel_io_t%output_single: MPI_File_write_at I/O', wallclock()-start, 'sec'
    if (io%method == 'snapshot') then
      if (verbose > 0) &
        write (io_unit%log,*) wallclock(), patch%id, patch%time, '  request_close'
      if (patch%iout-1 > closed) then
        !$omp atomic write
        closed = patch%iout-1
        call files_out(closed)%request_close ()
        !$omp atomic write
        check_files = .true.
      end if
    end if
  end if
  if (verbose > 1) then
    write (io_unit%log,*) wallclock(), patch%id, patch%time, ' END'
    flush (io_unit%log)
  end if
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  ! Copy variables from `patch%mem` to `mem` array, converting variables based
  ! on the I/O format parameter.  This procedure is internal to the 'output_single'
  ! procedure, and as such only supports io%format in the range 10-13.  In the
  ! process, the time slot ordering is rearranged into increasing time.
  !-----------------------------------------------------------------------------
  subroutine convert_variables
    integer:: nt1, l(3), u(3), i, j, k
    !...........................................................................
    call trace%begin ('parallel_io-t%output_single::convert_variables')
    if (io%guard_zones) then
      l = patch%mesh%lb
      u = max(patch%mesh%ub,patch%mesh%gn)
    else
      l = patch%mesh%li
      u = patch%mesh%ui
    end if
    nt1 = patch%nt-1
    if (io%format >= 12 .and. iv==patch%idx%d) then
      ! --- Convert density to log density ---
      do k=l(3),u(3); do j=l(2),u(2); do i=l(1),u(1)
        mem(i-l(1)+1,j-l(2)+1,k-l(3)+1,1:nt1) &
        = log(patch%mem(i,j,k,iv,patch%iit(1:nt1),1))
      end do; end do; end do
    else if ((io%format == 10 .or. io%format==12) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert momentum to velocity ---
      do k=l(3),u(3); do j=l(2),u(2); do i=l(1),u(1)
        mem(i-l(1)+1,j-l(2)+1,k-l(3)+1,1:nt1) &
          = patch%mem(i,j,k,iv,patch%iit(1:nt1),1) &
          / patch%mem(i,j,k,patch%idx%d,patch%iit(1:nt1),1)
      end do; end do; end do
    else
      mem = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,patch%iit(1:nt1),1)
    end if
    call trace%end()
  end subroutine
END SUBROUTINE output_single

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
  if (io%format > 9 .and. io%format < 14) then
    call input_single (patch, ok)
    call trace%end (itimer)
    return
  end if
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
    ! Using log(d) if io%format==8 or 9
    !-----------------------------------------------------------------------------
    if (io%format >= 8 .and. iv==patch%idx%d) then
      ! --- Convert log density to density ---
      patch%mem(:,:,:,iv,:,1) = exp(patch%mem(:,:,:,iv,:,1))
    end if
    !-----------------------------------------------------------------------------
    ! Restore to momentum variables for RAMSES snapshots, if format=6,8
    !-----------------------------------------------------------------------------
    if ((io%format == 6 .or. io%format==8) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert velocity to momentum ---
      patch%mem(:,:,:,iv,:,1) = patch%mem(:,:,:,iv,:,1)*patch%mem(:,:,:,patch%idx%d,:,1)
    end if
  end subroutine
END SUBROUTINE input

!===============================================================================
!> Read input in just one read per variable, reading all patches for this rank
!> in one go.  Use the same parallel_io%buffer as for writes.
!===============================================================================
SUBROUTINE input_single (patch, ok)
  class(patch_t):: patch
  logical:: ok
  !.............................................................................
  integer(8):: snapshot_words, buffer_words
  integer:: iv, n(3), slot, nder, irec, ip, iout, rank, nrank
  logical, save:: read_failed=.false., first_time=.true.
  integer, save:: itimer=0
  real, allocatable:: buffer(:,:,:)
  logical:: single_patch
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%input_single', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If the patch is a virtual patch, read nothing, but return ok=.true. never-
  ! theless, relying on the negative time to force nbor tasks to wait for the
  ! first MPI refresh of the task.
  !-----------------------------------------------------------------------------
  if (patch%is_set(bits%virtual)) then
    ok = .true.
    patch%time = -1.0
    call trace%end (itimer)
    return
  end if
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
  ! Slots in the file are ordered (ip,rank,iv,iout), so all patches for one
  ! variable that belong to this rank may be read at once
  !-----------------------------------------------------------------------------
  if (io%ntask==0) &
    call io%abort ('parallel_io_t%input_single: io%ntask == 0')
  !-----------------------------------------------------------------------------
  nder = io%time_derivs+1
  snapshot_words = product(n)*io%ntotal*io%nv*nder
  if (single_patch) then
    buffer_words = product(n)
  else
    buffer_words = product(n)*io%ntask
  end if
  !-----------------------------------------------------------------------------
  !$omp critical (allocate_cr)
  if (.not.associated(parallel_io%buffer)) then
    allocate (parallel_io%buffer(n(1),n(2),n(3),io%ntask,io%nv,nder))
    if (verbose > 0) &
      write (io%output,*) 'shape io%buffer =', shape(parallel_io%buffer)
  end if
  !$omp end critical (allocate_cr)
  !-----------------------------------------------------------------------------
  call mpi_io%use (file_in)
  call mpi_io%set (snapshot_words, buffer_words, 99)
  !-----------------------------------------------------------------------------
  ! Make sure only the first thread to come here allocates the buffer, if needed
  !-----------------------------------------------------------------------------
  patch%iout = io%restart
  ip = patch%ip
  if (single_patch) then
    !-----------------------------------------------------------------------------
    ! Info
    !-----------------------------------------------------------------------------
    if (first_time) then
    !$omp critical (first_cr)
    if (first_time) then
      write(io%output,'(a,f8.3,a,i10,a)') &
        ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
        io%ntotal, ' patches'
        write(io%output,'(a,4i6,2i12)') &
      ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
        patch%iout, n, snapshot_words, buffer_words
      first_time = .false.
    end if
    !$omp end critical (first_cr)
    end if
    allocate (buffer(n(1),n(2),n(3)))
    !-----------------------------------------------------------------------------
    ! The previous number of ranks is equal to the product of the previous mpi_dims.
    ! The previous rank is computable from what the old %ipos would have been.
    !-----------------------------------------------------------------------------
    nrank = product(io%mpi_odims)
    rank = mpi_coords%coords_to_rank (patch%ipos*io%mpi_odims/io%dims)
    if (verbose > 0) then
      irec = 1 + rank
      write (io_unit%output,'(a,i6,4(2x,3i5))') &
        'single_patch: id, ip, rank, nrank, per_rank, ipos =', &
          patch%id, ip, rank, nrank, io%dims/io%mpi_odims, patch%ipos, &
          irec, patch%iout
    end if
    do iv=1,io%nv
      irec = 1 + rank + nrank*(iv-1)
      if (verbose > 1) &
        write (io_unit%output,*) 'read(1): irec, iout =', irec, patch%iout
      call mpi_io%read (buffer, 1+patch%iout, irec)
      if (io%guard_zones) then
        patch%mem(:,:,:,iv,patch%it,1) = buffer
      else
        patch%mem(patch%li(1):patch%ui(1), &
                  patch%li(2):patch%ui(2), &
                  patch%li(3):patch%ui(3),iv,patch%it,1) = buffer
      end if
    end do
    deallocate (buffer)
  else
    if (first_time) then
    !$omp critical (read_first_cr)
    if (first_time) then
      !---------------------------------------------------------------------------
      ! Info
      !---------------------------------------------------------------------------
      write(io%output,'(a,f8.3,a,i10,a)') &
        ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
        io%ntotal, ' patches'
      write(io%output,'(a,4i6,2i12)') &
      ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
        patch%iout, n, snapshot_words, buffer_words
      !---------------------------------------------------------------------------
      ! Read via buffer, in case guard zones not included. FIXME: read snapshot
      ! guard zone parameter, to allow it to differ
      !---------------------------------------------------------------------------
      if (io%method == 'snapshot') then
        iout = 0
      else
        iout = patch%iout
      end if
      do iv=1,io%nv
        irec = 1 + mpi%rank + mpi%size*(iv-1)
        if (verbose > 1) &
          write (stdout,*) 'read(2): irec, iout =', irec, iout
        call mpi_io%read (parallel_io%buffer(:,:,:,:,iv,1), 1+iout, irec)
        if (mpi_io%err /= 0) then
          read_failed = .true.
          write(io%output,*) 'io%input_single: read failed', io%restart
          exit
        end if
      end do
      if (.not.read_failed) &
        write(io%output,*) 'io%input_single: read succeeded, snapshot', io%restart
      first_time = .false.
    end if
    !$omp end critical (read_first_cr)
    end if
    if (read_failed) then
      ok = .false.
      call trace%end (itimer)
      return
    end if
    if (ip < 1 .or. ip > io%ntask) then
      write (io_unit%output,*) &
        'parallel_io_t%input_single: WARNING io, ntask =', ip, io%ntask
      flush (io_unit%output)
      ip = max(1,min(ip,io%ntask))
    end if
    if (io%guard_zones) then
      patch%mem(:,:,:,:,patch%it,1) = parallel_io%buffer(:,:,:,ip,:,1)
    else
      patch%mem(patch%li(1):patch%ui(1), &
                patch%li(2):patch%ui(2), &
                patch%li(3):patch%ui(3),:,patch%it,1) = parallel_io%buffer(:,:,:,ip,:,1)
    end if
  end if
  if (verbose > 1) then
    do iv=1,io%nv
      write (io_unit%output,'(a,i6,3i4,1p,2g14.6)') &
        'parallel_io_t%input: id, ip, iv, iout, min, max =', &
        patch%id, ip, iv, patch%iout, patch%fminval(iv), patch%fmaxval(iv)
    end do
  else if (verbose > 0) then
    write (io_unit%output,*) 'parallel_io_t%input: id, ip, iout =', &
      patch%id, ip, patch%iout
  end if
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
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  ! Copy variables from `patch%mem` to `mem` array, converting variables based
  ! on the I/O format parameter.
  !-----------------------------------------------------------------------------
  subroutine convert_variables
    integer:: l(3), u(3)
    call trace%begin('input_single::convert_variables')
    if (io%guard_zones) then
      l = patch%mesh%lb
      u = max(patch%mesh%ub,patch%mesh%gn)
    else
      l = patch%mesh%li
      u = patch%mesh%ui
    end if
    !-----------------------------------------------------------------------------
    ! Using log(d) if io%format==12 or 13
    !-----------------------------------------------------------------------------
    if (io%format >= 12 .and. iv==patch%idx%d) then
      ! --- Convert log density to density ---
          patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,:,1) = &
      exp(patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,:,1))
    !-----------------------------------------------------------------------------
    ! Restore to momentum variables for RAMSES snapshots, if format=10,12
    !-----------------------------------------------------------------------------
    else if ((io%format == 10 .or. io%format==12) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert velocity to momentum ---
      patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,:,1) = &
      patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,:,1) * &
      patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,:,1)
    end if
    call trace%end()
  end subroutine
END SUBROUTINE input_single

END MODULE parallel_io_mod
