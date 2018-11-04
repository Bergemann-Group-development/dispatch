!> =============================================================================
!> Use MPI parallel I/O to write everythin to a single file.  No critical regions
!> are needed here; all thread protection is one inside mpi_io_file_mod and
!> mpi_io_mod.
!> =============================================================================
MODULE parallel_io_mod
  USE patch_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_io_mod
  USE mpi_io_file_mod
  USE trace_mod
  USE time_slices_mod
  USE counters_mod
  implicit none
  private
  type, public:: parallel_io_t
    real, dimension(:,:,:,:,:,:), pointer:: buffer => null()
  contains
    procedure, nopass:: init
    procedure, nopass:: input
    procedure, nopass:: output
    procedure, nopass:: close
  end type
  integer, parameter:: verbose=0
  type(mpi_io_file_t), save:: file_in, file_out
  type(parallel_io_t), public:: parallel_io
CONTAINS

!> =============================================================================
!> Initialize for MPI paraallel I/O.  This is an MPI collective operation, and
!> needs to be done at the start of the run.
!> =============================================================================
SUBROUTINE init (patch)
  class(patch_t):: patch
  !.............................................................................
  character(len=5):: snapname
  character(len=120):: filename
  logical:: exists
  integer:: recl
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%init')
  call mpi_io%init
  if (io%restart > 0) then
    filename=trim(io%inputdir)//'snapshots.dat'
    inquire (file=filename, exist=exists)
    if (.not.exists) then
      filename=trim(io%inputdir)//'00000/snapshot.dat'
      inquire (file=filename, exist=exists)
      if (.not.exists) then
        if (io%master) write(io%output,*) 'neither '//trim(io%inputdir) &
          //'snapshots.dat, nor '//trim(filename)//' exist'
        call mpi%abort()
      end if
    end if
    if (mpi%master) write(io%output,*) ' input file: ',filename
    call file_in%openr (filename)
  end if
  if (io%format==5) then
    filename=trim(io%outputname)//'00000/snapshot.dat'
  else
    filename=trim(io%outputname)//'snapshots.dat'
  end if
  if (mpi%master) write(io%output,*) 'output file: ',filename
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    recl = io%nv*product(patch%gn)*4
  else
    recl = io%nv*product(patch%n)*4
  end if
  call file_out%openw (filename, recl)
  io%nwrite = io%ntask
  call mpi_io%set (nwrite=io%nwrite*(io%time_derivs+1))
  call time_slices%init (patch)
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Close the data I/O method
!> =============================================================================
SUBROUTINE close
  !.............................................................................
  call trace%begin ('parallel_io_t%close')
  call file_out%close
  call file_in%close
  call trace%end
END SUBROUTINE close

!> =============================================================================
!> Write output in increments of one patch data cube
!> =============================================================================
SUBROUTINE output (patch, count)
  class(patch_t):: patch
  integer:: count
  !.............................................................................
  integer(8):: snapshot_words, buffer_words, pos
  integer:: iv, n(3), jt(2), iw, nw
  real, dimension(:,:,:,:), pointer:: buffer, mem
  real:: pt(2)
  logical:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (io%format > 9 .and. io%format < 14) then
    call output_single (patch, count)
    return
  end if
  call trace%begin ('parallel_io_t%output', itimer=itimer)
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
  if (io%format > 9 .and. io%format < 14) then
    buffer_words = snapshot_words
  else
    buffer_words = product(n)*io%nv
  end if
  call mpi_io%use (file_out)
  call mpi_io%set (snapshot_words, buffer_words, io%nwrite*(io%time_derivs+1))
  !-----------------------------------------------------------------------------
  ! Info
  !-----------------------------------------------------------------------------
  if (io%verbose==1) write(io%output,'(a,2i6)') ' parallel_io_t%output: iout, id =', &
    patch%iout, patch%id
  if (first_time) then
    first_time = .false.
    if (io%master) write(io%output,'(a,f8.3,a,i10,a)') &
      ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
  end if
  if (io%verbose>1) write(io%output,'(a,4i6,2i12)') &
    ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  allocate (mem(patch%gn(1),patch%gn(2),patch%gn(3),patch%nt))
  !-----------------------------------------------------------------------------
  ! Time interpolation to exaxt out_next
  !-----------------------------------------------------------------------------
  allocate (buffer(n(1),n(2),n(3),io%nv))
  do iv=1,io%nv
    call convert_variables
    call time_slices%interpolate (patch, mem, buffer(:,:,:,iv))
  end do
  call io_write ((io%time_derivs+1)*patch%iout, 0)
  !-----------------------------------------------------------------------------
  ! On-the-fly time derivatives output
  !-----------------------------------------------------------------------------
  if (io%time_derivs>0) then
    allocate (buffer(n(1),n(2),n(3),io%nv))
    do iv=1,io%nv
      call convert_variables
      call time_slices%derivative1 (patch, mem, buffer(:,:,:,iv))
    end do
    call io_write ((io%time_derivs+1)*patch%iout, 1)
    if (io%time_derivs==2) then
      allocate (buffer(n(1),n(2),n(3),io%nv))
      do iv=1,io%nv
        call convert_variables
        call time_slices%derivative2 (patch, mem, buffer(:,:,:,iv))
      end do
      call io_write ((io%time_derivs+1)*patch%iout, 2)
    end if
  end if
  deallocate (mem)
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  subroutine convert_variables
    integer:: ioformat
    !...........................................................................
    if (io%format > 9 .and. io%format < 14) then
      ioformat = io%format - 4
    else
      ioformat = io%format
    end if
    if (ioformat >= 8 .and. iv==patch%idx%d) then
      ! --- Convert density to log density ---
      mem = log(patch%mem(:,:,:,iv,:,1))
    else if ((ioformat == 6 .or. ioformat==8) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert momentum to velocity ---
      mem = patch%mem(:,:,:,iv,:,1)/patch%mem(:,:,:,patch%idx%d,:,1)
    else
      mem = patch%mem(:,:,:,iv,:,1)
    end if
  end subroutine
  !-----------------------------------------------------------------------------
  subroutine io_write (slot, ider)
    integer:: slot, ider, iv, islot, m(4)
    real, dimension(:,:,:), pointer:: buf
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
    ! This ordering has all patch data for one variable inside the same N**3
    ! section of the disk, so it can be read directly, and reordered in memory
    !---------------------------------------------------------------------------
    else if (io%format > 9 .and. io%format < 14) then
      m = rank(buffer)
      do iv=1,io%nv
        islot = 1 + iv + io%nv*(ider + patch%iout*(io%time_derivs+1))
        allocate (buf(m(1),m(2),m(3)))
        buf = buffer(:,:,:,iv)
        call mpi_io%iwrite (buf, islot, patch%id)
      end do
      deallocate (buffer)
    !---------------------------------------------------------------------------
    ! This ordering has the derivative order added to the patch%id, which is
    ! in turn causing an increment by (io%time_derivs+1), so effectively adding
    ! nv additional variables for each derivative ordre, making it easy to
    ! address the derivatives, while using snapshot numbers incremented by 1.
    !---------------------------------------------------------------------------
    else
      call mpi_io%iwrite (buffer, 1+slot, 1+ider+(patch%id-1)*(io%time_derivs+1))
    end if
  end subroutine
END SUBROUTINE output

!> =============================================================================
!> Write output as one single buffer
!> =============================================================================
SUBROUTINE output_single (patch, count)
  class(patch_t):: patch
  integer:: count
  !.............................................................................
  integer(8):: snapshot_words, buffer_words, pos
  integer:: iv, n(3), jt(2), iw, nw, nder, irec, ider, ip
  real, dimension(:,:,:,:), pointer:: mem
  real, dimension(:,:,:), pointer:: buffer
  real:: pt(2)
  real(8):: start
  logical:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output_single', itimer=itimer)
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
  nder = io%time_derivs+1
  snapshot_words = product(n)*io%ntotal*io%nv*nder
  buffer_words = product(n)*io%ntask
  call mpi_io%use (file_out)
  call mpi_io%set (snapshot_words, buffer_words, 99)
  if (.not.associated(parallel_io%buffer)) then
    nder = io%time_derivs+1
    allocate (parallel_io%buffer(n(1),n(2),n(3),io%ntask,io%nv,nder))
    write (io%output,*) 'shape io%buffer =', shape(parallel_io%buffer)
  end if
  !-----------------------------------------------------------------------------
  ! Info
  !-----------------------------------------------------------------------------
  if (first_time .and. io%master) then
    first_time = .false.
    write(io%output,'(a,f8.3,a,i10,a)') &
      ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
  end if
  ! --- rank-local patch 1D index; cf. cartesian_mod ---
  ip = patch%ip
  if (ip < 0 .or. ip > io%ntask) then
    write(io%output,*) 'WARNING', ip
    ip = max(1,min(ip,io%ntask))
  end if
  if (io%verbose>0) then
    write(io%output,'(a,6i6,2i12,l3)') &
      ' parallel_io_t%output_single: id, ip, iout, n, snapshot_words, buffer_words =', &
      patch%id, ip, patch%iout, n, snapshot_words, buffer_words, associated(parallel_io%buffer)
    flush (io%output)
  end if
  allocate (mem(n(1),n(2),n(3),patch%nt-1))
  !-----------------------------------------------------------------------------
  ! Time interpolation to exaxt out_next
  !-----------------------------------------------------------------------------
  do iv=1,io%nv
    call convert_variables
    buffer => parallel_io%buffer(:,:,:,ip,iv,1)
    call time_slices%interpolate (patch, mem, buffer)
  end do
  !-----------------------------------------------------------------------------
  ! On-the-fly time derivatives output
  !-----------------------------------------------------------------------------
  if (io%time_derivs>0) then
    do iv=1,io%nv
      call convert_variables
      buffer => parallel_io%buffer(:,:,:,ip,iv,2)
      call time_slices%derivative1 (patch, mem, buffer)
    end do
    if (io%time_derivs==2) then
      do iv=1,io%nv
        call convert_variables
        buffer => parallel_io%buffer(:,:,:,ip,iv,3)
        call time_slices%derivative2 (patch, mem, buffer)
      end do
    end if
  end if
  deallocate (mem)
  if (count == 0) then
    start = wallclock()
    do ider=1,nder
    do iv=1,io%nv
      !-------------------------------------------------------------------------
      ! Each "record" is product(n)*io%ntasks words long, and contains values
      ! for one variable, for all the patches in a rank.  The ranks are written
      ! out one after the other for one variable, followed by then next var
      !-------------------------------------------------------------------------
      irec = 1 + mpi%rank + mpi%size*(iv-1 + (ider-1)*io%nv)
      call mpi_io%write (parallel_io%buffer(:,:,:,:,iv,ider), 1+patch%iout, irec)
    end do
    end do
    write (io%output,'(a,f7.3,1x,a)') &
      'parallel_io%output_single: MPI_File_write_at I/O', wallclock()-start, 'sec'
  end if
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  subroutine convert_variables
    integer:: nt1, l(3), u(3)
    !...........................................................................
    if (io%guard_zones) then
      l = patch%mesh%lb
      u = patch%mesh%ub
    else
      l = patch%mesh%li
      u = patch%mesh%ui
    end if
    nt1 = patch%nt-1
    if (io%format >= 8 .and. iv==patch%idx%d) then
      ! --- Convert density to log density ---
      mem = log(patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,patch%iit(1:nt1),1))
    else if ((io%format == 6 .or. io%format==8) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert momentum to velocity ---
      mem = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,patch%iit(1:nt1),1) &
          / patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%iit(1:nt1),1)
    else
      mem = patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,patch%iit(1:nt1),1)
    end if
  end subroutine
END SUBROUTINE output_single

!> =============================================================================
!> Read input in increments of one variable data cube
!> =============================================================================
SUBROUTINE input (patch, ok)
  class(patch_t):: patch
  logical:: ok
  !.............................................................................
  integer(8):: snapshot_words, buffer_words
  integer:: iv, n(3), slot
  real, dimension(:,:,:,:), pointer:: buffer
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (io%format > 9 .and. io%format < 14) then
    !!omp critical (parallel_in_cr)
    call input_single (patch, ok)
    !!omp end critical (parallel_in_cr)
    return
  end if
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%input', itimer=itimer)
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
  ! Compute offset size constants for the mpi_io reader
  !-----------------------------------------------------------------------------
  if (io%format > 9 .and. io%format < 14) then
    buffer_words = product(n)
  else
    buffer_words = product(n)*io%nv
  end if
  snapshot_words = buffer_words*io%ntotal
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
    patch%mem(:,:,:,:,patch%it,1) = buffer
  else
    patch%mem(patch%li(1):patch%ui(1), &
              patch%li(2):patch%ui(2), &
              patch%li(3):patch%ui(3),:,patch%it,1) = buffer
  end if
  if (io%verbose>1) then
    do iv=1,io%nv
      write (io_unit%output,*) 'parallel_io_t%input: id, iv, iout, min, max =', &
        patch%id, iv, patch%iout, patch%fminval(iv), patch%fmaxval(iv)
    end do
  end if
  deallocate (buffer)
  do iv=1,io%nv
    call convert_variables
  end do
  !-----------------------------------------------------------------------------
  ! Make sure all time slots have the correct time set
  !-----------------------------------------------------------------------------
  patch%time = io%restart*io%out_time
  patch%t(patch%iit) = patch%time
  !-----------------------------------------------------------------------------
  ! Signal success
  !-----------------------------------------------------------------------------
  ok = .true.
  if (io%verbose==1) &
    write(io%output,'(a,2i6,1p,e16.6)') ' parallel_io_t%input: iout, id, time =', &
    patch%iout, patch%id, patch%time
  if (io%verbose>1) &
    write(io%output,'(a,4i6,2i7)') &
    ' parallel_io_t%input: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call trace%end (itimer)
contains
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
    ! Restore to momentum variables for RAMSES snapshots
    !-----------------------------------------------------------------------------
    if ((io%format == 6 .or. io%format==8) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert velocity to momentum ---
      patch%mem(:,:,:,iv,:,1) = patch%mem(:,:,:,iv,:,1)*patch%mem(:,:,:,patch%idx%d,:,1)
    end if
  end subroutine
END SUBROUTINE input

!> =============================================================================
!> Read input in just one read per variable, reading all patches for this rank
!> in one go.  Use the same parallel_io%buffer as for writes.
!> =============================================================================
SUBROUTINE input_single (patch, ok)
  class(patch_t):: patch
  logical:: ok
  !.............................................................................
  integer(8):: snapshot_words, buffer_words
  integer:: iv, n(3), slot, nder, irec
  logical, save:: read_failed=.false., first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%input_single', itimer=itimer)
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
  ! Slots in the file are ordered (ip,rank,iv,iout), so all patches for one
  ! variable that belong to this rank may be read at once
  !-----------------------------------------------------------------------------
  nder = io%time_derivs+1
  snapshot_words = product(n)*io%ntotal*io%nv*nder
  buffer_words = product(n)*io%ntask
  call mpi_io%use (file_in)
  call mpi_io%set (snapshot_words, buffer_words, 99)
  !-----------------------------------------------------------------------------
  ! Make sure only the first thread to come here allocates the buffer, if needed
  !-----------------------------------------------------------------------------
  patch%iout = io%restart
  !$omp critical (parallel_in_cr)
  if (first_time) then
    first_time = .false.
    if (.not.associated(parallel_io%buffer)) then
      nder = io%time_derivs+1
      allocate (parallel_io%buffer(n(1),n(2),n(3),io%ntask,io%nv,nder))
      write (io%output,*) 'shape io%buffer =', shape(parallel_io%buffer)
    end if
    !---------------------------------------------------------------------------
    ! Read via buffer, in case guard zones not included. FIXME: read snapshot
    ! guard zone parameter, to allow it to differ
    !---------------------------------------------------------------------------
    do iv=1,io%nv
      irec = 1 + mpi%rank + mpi%size*(iv-1)
      call mpi_io%read (parallel_io%buffer(:,:,:,:,iv,1), 1+patch%iout, irec)
    end do
    if (mpi_io%err /= 0) then
      read_failed = .true.
      write(io%output,*) 'io%parallel_io: read failed', io%restart
    else
      write(io%output,*) 'io%parallel_io: read succeeded', io%restart
    end if
  end if
  !$omp end critical (parallel_in_cr)
  if (read_failed) then
    ok = .false.
    call trace%end (itimer)
    return
  end if
  if (io%guard_zones) then
    patch%mem(:,:,:,:,patch%it,1) = parallel_io%buffer(:,:,:,patch%ip,:,1)
  else
    patch%mem(patch%li(1):patch%ui(1), &
              patch%li(2):patch%ui(2), &
              patch%li(3):patch%ui(3),:,patch%it,1) = parallel_io%buffer(:,:,:,patch%ip,:,1)
  end if
  if (io%verbose>1) then
    do iv=1,io%nv
      write (io_unit%output,*) 'parallel_io_t%input: id, iv, iout, min, max =', &
        patch%id, iv, patch%iout, patch%fminval(iv), patch%fmaxval(iv)
    end do
  end if
  do iv=1,io%nv
    call convert_variables
  end do
  !-----------------------------------------------------------------------------
  ! Make sure all time slots have the correct time set
  !-----------------------------------------------------------------------------
  patch%time = io%restart*io%out_time
  patch%t(patch%iit) = patch%time
  !-----------------------------------------------------------------------------
  ! Signal success
  !-----------------------------------------------------------------------------
  ok = .true.
  call trace%end (itimer)
contains
  !-----------------------------------------------------------------------------
  subroutine convert_variables
    !-----------------------------------------------------------------------------
    ! Using log(d) if io%format==12 or 13
    !-----------------------------------------------------------------------------
    if (io%format >= 12 .and. iv==patch%idx%d) then
      ! --- Convert log density to density ---
      patch%mem(:,:,:,iv,:,1) = exp(patch%mem(:,:,:,iv,:,1))
    end if
    !-----------------------------------------------------------------------------
    ! Restore to momentum variables for RAMSES snapshots
    !-----------------------------------------------------------------------------
    if ((io%format == 10 .or. io%format==12) .and. patch%solver_is('ramses') &
        .and. iv >= patch%idx%px .and. iv <= patch%idx%pz) then
      ! --- Convert velocity to momentum ---
      patch%mem(:,:,:,iv,:,1) = patch%mem(:,:,:,iv,:,1)*patch%mem(:,:,:,patch%idx%d,:,1)
    end if
  end subroutine
END SUBROUTINE input_single

END MODULE parallel_io_mod
