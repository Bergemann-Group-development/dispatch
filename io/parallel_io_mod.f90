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
  USE mpi_mod
  USE mpi_io_mod
  USE mpi_io_file_mod
  USE trace_mod
  implicit none
  private
  type, public:: parallel_io_t
  contains
    procedure, nopass:: init
    procedure, nopass:: input
    procedure, nopass:: output
    procedure, nopass:: close
  end type
  integer, parameter:: order=1
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
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%init')
  if (patch%restart > 0) then
    filename=trim(io%inputname)//'snapshots.dat'
    inquire (file=filename, exist=exists)
    if (.not.exists) then
      filename=trim(io%inputname)//'00000/snapshot.dat'
      inquire (file=filename, exist=exists)
      if (.not.exists) then
        call mpi%abort('input file '//trim(filename)//' not found')
      end if
    end if
    if (mpi%master) print *,' input file: ',filename
    call file_in%openr (filename)
  end if
  if (io%format==5) then
    filename=trim(io%outputname)//'00000/snapshot.dat'
  else
    filename=trim(io%outputname)//'snapshots.dat'
  end if
  if (mpi%master) print *,'output file: ',filename
  call file_out%openw (filename)
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
!> Write output in increments of one variable data cube
!> =============================================================================
SUBROUTINE output (patch)
  class(patch_t):: patch
  if (io%format==5) then
    call output5 (patch)
  else if (io%format==6) then
    call output6 (patch)
  end if
END SUBROUTINE output

!> =============================================================================
!> Write output in increments of one variable data cube
!> =============================================================================
SUBROUTINE output5 (patch)
  class(patch_t):: patch
  !.............................................................................
  integer(8):: snapshot_words, buffer_words
  integer:: iv, n(3), jt(2), nv
  real, dimension(:,:,:), allocatable:: buffer
  type(mpi_io_t):: thread_io
  real:: pt(2)
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output')
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call thread_io%use (file_out)
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  buffer_words = product(n)
  snapshot_words = buffer_words*patch%nv*io%ntotal
  if (io%verbose==1) print '(a,2i6)', ' parallel_io_t%output: iout, id =', &
    patch%iout, patch%id
  if (io%verbose>1) print '(a,4i6,2i7)', &
    ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call thread_io%init (snapshot_words, buffer_words)
  allocate (buffer(n(1),n(2),n(3)))
  call patch%time_interval (patch%out_next, jt, pt)
  nv = merge (io%nv, patch%nv, io%nv>0)
  do iv=1,nv
    if (io%guard_zones) then
      buffer = patch%mem(:,:,:,iv,jt(1),1)*pt(1) &
             + patch%mem(:,:,:,iv,jt(2),1)*pt(2)
    else
      buffer = patch%mem(patch%li(1):patch%ui(1), &
                         patch%li(2):patch%ui(2), &
                         patch%li(3):patch%ui(3),iv,jt(1),1)*pt(1) &
             + patch%mem(patch%li(1):patch%ui(1), &
                         patch%li(2):patch%ui(2), &
                         patch%li(3):patch%ui(3),iv,jt(2),1)*pt(2)
    end if
    if (io%verbose>2) print *, 'parallel_io_t%output: iv =', &
      patch%id, iv, minval(buffer), maxval(buffer)
    if (order==0) then
      call thread_io%write (buffer, patch%iout+1, (patch%id-1)*patch%nv + iv)
    else
      call thread_io%write (buffer, patch%iout+1, patch%id + (iv-1)*io%ntotal)
    end if
    if (verbose>0) &
      write (io_unit%mpi,*) 'WR5', patch%id, patch%iout, iv, minval(buffer), maxval(buffer)
  end do
  deallocate (buffer)
  call trace%end
END SUBROUTINE output5

!> =============================================================================
!> Write output in increments of one variable data cube
!> =============================================================================
SUBROUTINE output6 (patch)
  class(patch_t):: patch
  !.............................................................................
  integer(8):: snapshot_words, buffer_words, pos
  integer:: iv, n(3), jt(2), nv
  real, dimension(:,:,:,:), allocatable:: buffer
  type(mpi_io_t):: thread_io
  real:: pt(2)
  logical:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%output')
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call thread_io%use (file_out)
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  buffer_words = product(n)*patch%nv
  snapshot_words = buffer_words*io%ntotal
  if (io%verbose==1) print '(a,2i6)', ' parallel_io_t%output: iout, id =', &
    patch%iout, patch%id
  if (first_time) then
    first_time = .false.
    if (io%master) print '(a,f8.3,a,i10,a)', &
      ' output snapshot size =', snapshot_words*4d0/1024d0**3, ' GB, in', &
      io%ntotal, ' patches'
  end if
  if (io%verbose>1) print '(a,4i6,2i7)', &
    ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call thread_io%init (snapshot_words, buffer_words)
  allocate (buffer(n(1),n(2),n(3),patch%nv))
  call patch%time_interval (patch%out_next, jt, pt)
  nv = merge (io%nv, patch%nv, io%nv>0)
  do iv=1,nv
    if (io%guard_zones) then
      buffer(:,:,:,iv) = patch%mem(:,:,:,iv,jt(1),1)*pt(1) &
                       + patch%mem(:,:,:,iv,jt(2),1)*pt(2)
    else
      buffer(:,:,:,iv) = patch%mem(patch%li(1):patch%ui(1), &
                                   patch%li(2):patch%ui(2), &
                                   patch%li(3):patch%ui(3),iv,jt(1),1)*pt(1) &
                       + patch%mem(patch%li(1):patch%ui(1), &
                                   patch%li(2):patch%ui(2), &
                                   patch%li(3):patch%ui(3),iv,jt(2),1)*pt(2)
    end if
    if (verbose>0) &
      write (io_unit%mpi,'(a,i6,i4,i3,1p,6e12.4,2i3,0p,2f7.3)') &
       'id,iout,iv,vmin,vmax:', patch%id, patch%iout, iv, &
        minval(buffer(:,:,:,iv)), maxval(buffer(:,:,:,iv)), &
        patch%fminval(patch%mem(:,:,:,iv,jt(1),1)), &
        patch%fmaxval(patch%mem(:,:,:,iv,jt(1),1)), &
        patch%fminval(patch%mem(:,:,:,iv,jt(2),1)), &
        patch%fmaxval(patch%mem(:,:,:,iv,jt(2),1)), &
        jt, pt
  end do
  if (verbose > 0) then
    pos = 4_8*(patch%iout*snapshot_words + (patch%id-1_8)*buffer_words)
    write (io_unit%mpi,*) 'WR6', patch%id, patch%iout, pos
  end if
  call thread_io%write (buffer, patch%iout+1, patch%id)
  deallocate (buffer)
  call trace%end
END SUBROUTINE output6

!> =============================================================================
!> Read input in increments of one variable data cube
!> =============================================================================
SUBROUTINE input (patch, ok)
  class(patch_t):: patch
  logical:: ok
  !.............................................................................
  integer(8):: snapshot_words, buffer_words
  integer:: iv, n(3), nv
  real, dimension(:,:,:), allocatable:: buffer
  type(mpi_io_t):: thread_io
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%input')
  if (file_in%handle==0) then
    call init (patch)
  end if
  call thread_io%use (file_in)
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  buffer_words = product(n)
  snapshot_words = buffer_words*patch%nv*io%ntotal
  call thread_io%init (snapshot_words, buffer_words)
  allocate (buffer(n(1),n(2),n(3)))
  nv = merge (io%nv, patch%nv, io%nv>0)
  do iv=1,nv
    if (order==0) then
      call thread_io%read (buffer, patch%restart+1, (patch%id-1)*patch%nv + iv)
    else
      call thread_io%read (buffer, patch%restart+1, patch%id + (iv-1)*io%ntotal)
    end if
    if (thread_io%err /= 0) then
      ok = .false.
      call trace%end
      return
    end if
    if (io%guard_zones) then
      patch%mem(:,:,:,iv,patch%it,1) = buffer
    else
      patch%mem(patch%li(1):patch%ui(1), &
                patch%li(2):patch%ui(2), &
                patch%li(3):patch%ui(3),iv,patch%it,1) = buffer
    end if
    if (verbose>0) &
      write (io_unit%mpi,*) 'RD', patch%id, patch%iout, iv, minval(buffer), maxval(buffer)
  end do
  deallocate (buffer)
  ok = .true.
  patch%time = patch%restart*io%out_time
  if (io%verbose==1) print '(a,2i6,1p,e16.6)', ' parallel_io_t%input: iout, id, time =', &
    patch%iout, patch%id, patch%time
  if (io%verbose>1) print '(a,4i6,2i7)', &
    ' parallel_io_t%input: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call trace%end
END SUBROUTINE input

END MODULE parallel_io_mod
