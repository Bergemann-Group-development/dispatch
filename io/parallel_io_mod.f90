!> =============================================================================
!> Use MPI parallel I/O to write everythin to a single file.  No critical regions
!> are needed here; all thread protection is one inside mpi_io_file_mod and
!> mpi_io_mod.
!> =============================================================================
MODULE parallel_io_mod
  USE patch_mod
  USE io_mod
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
  integer, save:: verbose=0
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
  character(len=5):: filename
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%init')
  write (filename,'(i5.5)') patch%iout
  call file_out%openw (trim(io%outputname)//filename//'/snapshot.dat')
  if (patch%restart > 0) then
    write (filename,'(i5.5)') patch%restart
    call file_in%openr (trim(io%outputname)//filename//'/snapshot.dat')
  end if
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
  !.............................................................................
  integer:: snapshot_words, buffer_words, iv, n(3), jt(2)
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
  snapshot_words = buffer_words*patch%nv*io%ntask
  if (io%verbose==1) print '(a,2i6)', ' parallel_io_t%output: iout, id =', &
    patch%iout, patch%id
  if (io%verbose>1) print '(a,4i6,2i7)', &
    ' parallel_io_t%output: iout, n, snapshot_words, buffer_words =', &
    patch%iout, n, snapshot_words, buffer_words
  call thread_io%init (snapshot_words, buffer_words)
  allocate (buffer(n(1),n(2),n(3)))
  call patch%time_interval (patch%out_next, jt, pt)
  do iv=1,patch%nv
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
    call thread_io%write (buffer, patch%iout+1, (patch%id-1)*patch%nv + iv)
  end do
  deallocate (buffer)
  call trace%end
END SUBROUTINE output

!> =============================================================================
!> Read input in increments of one variable data cube
!> =============================================================================
SUBROUTINE input (patch, ok)
  class(patch_t):: patch
  logical:: ok
  !.............................................................................
  integer:: snapshot_words, buffer_words, iv, n(3)
  real, dimension(:,:,:), allocatable:: buffer
  type(mpi_io_t):: thread_io
  !-----------------------------------------------------------------------------
  call trace%begin ('parallel_io_t%input')
  call thread_io%use (file_in)
  if (io%guard_zones) then
    n = patch%gn
  else
    n = patch%n
  end if
  buffer_words = product(n)
  snapshot_words = buffer_words*patch%nv
  call thread_io%init (snapshot_words, buffer_words)
  allocate (buffer(n(1),n(2),n(3)))
  do iv=1,patch%nv
    call thread_io%read (buffer, patch%restart+1, (patch%id-1)*patch%nv + iv)
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
  end do
  deallocate (buffer)
  ok = .true.
  call trace%end
END SUBROUTINE input

END MODULE parallel_io_mod
