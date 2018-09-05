!> *****************************************************************************
!> Module for handling blocking and non-blocking MPI parallel I/O to a single file.
!> The module is initialized by specifying the chunk size (which can be 1- or
!> 3-dimensional). A chunk could be just a patch, or a single chunk of data per
!> rank.
!>
!> The module can be used to handle any collection of patches, by using the
!> patch%id number to map each patch into the file.  If the patch sizes are fixed,
!> this is a unique mapping, and if each snapshot is written into a separate file
!> this is also sufficiently general.  The header part of the file may be written
!> to a separate file.
!>
!> Typical multi-threaded use:
!>
!>   type(mpi_io_file_t):: file                         ! holds file info
!>   type(mpi_io_t):: thread_io                         ! holds thread info
!>   !$omp threadprivate (thread_io)                    ! threadprivate data type
!>   ...
!>   call file%openw ('some_file.dat')                  ! implicit MPI barrier
!>   call thread_io%use (file)                          ! all threads use the same file
!>   ...
!>   call thread_io%init (snapshot_size, patch_size)    ! define offset factors
!>   call thread_io%write (patch_data, iout1, id1)      ! write snapshot iout1, patch id1
!>   ...
!>   call thread_io%read (patch_data, iout2, id2)       ! read snapshot iout2, patch id2
!>   ...
!>   !$omp barrier                                      ! all threads must agree
!>   call file%close                                    ! implicit MPI barrier
!>
!> *****************************************************************************
MODULE mpi_io_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE trace_mod
  USE mpi_mod, only: mp=>mpi
  USE mpi_io_file_mod
  USE dll_mod
#ifdef MPI
  USE mpi
#endif MPI
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! dll_node extension with a request and buffer pointer
  !-----------------------------------------------------------------------------
  type, public, extends(dll_node_t):: iwrite_t
    integer:: req, handle, words
    integer(kind=MPI_OFFSET_KIND):: pos
    real, dimension(:,:,:,:), pointer:: buffer
  end type
  !-----------------------------------------------------------------------------
  ! dll_t extension with a procedure for testing request completion
  !-----------------------------------------------------------------------------
  type, public, extends(dll_t):: iwrite_list_t
    logical:: first_time=.true.
    integer:: thread=-1
  contains
    procedure:: check
  end type
  !-----------------------------------------------------------------------------
  ! MPI I/O data type
  !-----------------------------------------------------------------------------
  type, public:: mpi_io_t
    type(mpi_io_file_t):: file
    integer:: req
    integer(8):: rec_words
    integer:: chunk_words
    integer:: err
    real, dimension(:,:,:), allocatable:: buf
#ifdef MPI
    integer:: status(MPI_STATUS_SIZE)
#else
    integer:: status
#endif
  contains
    procedure:: use
    procedure:: init
    procedure:: set
    procedure:: check_init
    procedure, private:: write3
    procedure, private:: write4
    generic, public :: write => write3, write4
    procedure:: iwrite
    procedure, private:: read3
    procedure, private:: read4
    generic, public :: read => read3, read4
    procedure:: assert
  end type
  character(len=32), save:: fmt='(1x,a,4i8,1p,5e12.3)'
  integer, save:: verbose=0
  type(iwrite_list_t), public:: iwrite_list
CONTAINS

!> =============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!> =============================================================================
SUBROUTINE use (self, file)
  class(mpi_io_t):: self
  class(mpi_io_file_t):: file
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%use')
  self%file = file
  call trace%end
END SUBROUTINE use

!> =============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!> =============================================================================
SUBROUTINE init (self)
  class(mpi_io_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%init')
  !$omp critical (iwrite_cr)
  if (iwrite_list%first_time) then
    iwrite_list%first_time = .false.
    call iwrite_list%init ('iwrite_list')
    print *, 'mpi_io_t%init done'
  end if
  !$omp end critical (iwrite_cr)
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!> =============================================================================
SUBROUTINE set (self, rec_words, chunk_words)
  class(mpi_io_t):: self
  integer(8), optional:: rec_words, chunk_words
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%init')
  if (present(rec_words)) then
    self%rec_words = rec_words
    if (verbose > 2) print *, 'mpi_io_t%init: setting rec_words =', rec_words
  end if
  if (present(chunk_words)) then
    self%chunk_words = chunk_words
    if (verbose > 2) print *, 'mpi_io_t%init: setting chunk_words =', chunk_words
  end if
  call trace%end
END SUBROUTINE set

!> =============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!> =============================================================================
SUBROUTINE check_init (self)
  class(mpi_io_t):: self
  !-----------------------------------------------------------------------------
  call self%assert ('mpi_io_t:rec_words is not positive', err=self%rec_words <= 0)
  call self%assert ('mpi_io_t:chunk_words is not positive', err=self%chunk_words <= 0)
END SUBROUTINE check_init

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE write3 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer:: rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%write')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose>1) &
    write (io_unit%mpi, '(a,3i8,i15,1p,2e14.6)') ' mpi_io_t%write3: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !!omp critical (mpi_io_cr)
  call MPI_file_write_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                          self%status, self%err)
  !!omp end critical (mpi_io_cr)
  call self%assert ('mpi_io_mod::write MPI_File_write_at '//trim(self%file%filename))
  call trace%end
#endif
END SUBROUTINE write3

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE write4 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:,:):: f
  integer:: rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%write')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose>1) &
    write (io_unit%mpi, '(a,3i8,i15,1p,2e14.6)') ' mpi_io_t%write4: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !!omp critical (mpi_io_cr)
  call MPI_file_write_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                          self%status, self%err)
  !!omp end critical (mpi_io_cr)
  call self%assert ('mpi_io_mod::write MPI_File_write_at '//trim(self%file%filename))
  call trace%end
#endif
END SUBROUTINE write4

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE iwrite (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:,:), pointer:: f
  integer:: rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  class(iwrite_t), pointer:: item
  class(dll_node_t), pointer:: node
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%iwrite', itimer=itimer)
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose>1) &
    print '(a,3i8,i15,2f12.0)', ' mpi_io_t%iwrite: rec, id, words, pos =', &
    rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !-----------------------------------------------------------------------------
  ! Add details we need in order to check for buffers to deallocate
  !-----------------------------------------------------------------------------
  allocate (item)
  item%handle = self%file%handle
  item%words  = self%chunk_words
  item%pos    = pos
  item%buffer => f
  !$omp critical (iwrite_cr)
  node => item
  call iwrite_list%append (node)
  if (verbose>1) &
    print '(1x,a,z12,i4)', 'appended req, n =', item%req, iwrite_list%n
  !$omp end critical (iwrite_cr)
  call trace%end (itimer)
#endif
END SUBROUTINE iwrite

!> =============================================================================
!> Read a chunk (self%lsize) from the open file
!> =============================================================================
SUBROUTINE read3 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  if (rec < 1) return
  call trace%begin ('mpi_io_t%read3')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  !$omp critical (mpi_io_cr)
  call MPI_File_read_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                         self%status, self%err)
  if (verbose>1) &
    write(io_unit%mpi, '(a,3i8,i15,1p,2e14.6)') ' mpi_io_t%read: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !$omp end critical (mpi_io_cr)
  call trace%end
#endif
END SUBROUTINE read3

!> =============================================================================
!> Read a chunk (self%lsize) from the open file
!> =============================================================================
SUBROUTINE read4 (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:,:):: f
  integer rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  if (rec < 1) return
  call trace%begin ('mpi_io_t%read4')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  !$omp critical (mpi_io_cr)
  call MPI_File_read_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                         self%status, self%err)
  if (verbose>1) &
    write(io_unit%mpi, '(a,3i8,i15,1p,2e14.6)') ' mpi_io_t%read: rec, id, words, pos =', &
      rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !$omp end critical (mpi_io_cr)
  call trace%end
#endif
END SUBROUTINE read4

!> =============================================================================
!> Assert no error
!> =============================================================================
SUBROUTINE assert (self, label, err)
  class(mpi_io_t):: self
  character(len=*):: label
  logical, optional:: err
  !-----------------------------------------------------------------------------
  if (present(err)) then
    if (err) then
      call mp%abort (label)
    end if
  else if (self%err /= 0) then
    print*,mp%rank,omp%thread,' error =',self%err
    call mp%abort (label)
  end if
END SUBROUTINE assert

!> =============================================================================
!> Multi-thread interface
!> =============================================================================
SUBROUTINE check (self)
  class(iwrite_list_t):: self
  if (omp%nthreads > 1) then
    !$omp task default(shared)
    call check_real (self)
    !$omp end task
  else
    call check_real (self)
  end if
END SUBROUTINE check

!> =============================================================================
!> Check for queued I/O, and do it all with one single thread, to avoid
!> slowing down the other threads.  Other threads are free to append more
!> items at the end of the iwrite_list in the mean time -- this does not
!> disturb the 
!> =============================================================================
SUBROUTINE check_real (self)
  class(iwrite_list_t):: self
  class(dll_node_t), pointer:: node, next
  class(iwrite_t), pointer:: item
  integer:: status(MPI_STATUS_SIZE), err
  integer, save:: itimer=0
  integer:: n, done, thread, omp_get_thread_num
  logical:: ok
  !-----------------------------------------------------------------------------
  !$omp atomic read
  thread = iwrite_list%thread
  if (iwrite_list%n == 0 .or. thread >=0) return
  !$omp critical (iwrite_cr)
  ok = (iwrite_list%thread == -1)
  iwrite_list%thread = omp_get_thread_num()
  call trace%begin ('iwrite_list_t%check_real', itimer=itimer)
  n = iwrite_list%n
  done = 0
  node => iwrite_list%head
  !$omp end critical (iwrite_cr)
  if (ok) then
    do while (n>0 .and. associated(node))
      err = 0
      select type (node)
      class is (iwrite_t)
        item => node
        call MPI_file_write_at (item%handle, item%pos, item%buffer, item%words, &
                                MPI_REAL, status, err)
      end select
      deallocate (item%buffer)
      !$omp critical (iwrite_cr)
      next => node%next
      call iwrite_list%remove (node)
      n = iwrite_list%n
      done = done+1
      if (verbose > 0 .and. mod(done,10000) == 0) &
        print '(1p,g14.6,3x,a,i7)', wallclock(), 'I/O done =', done
      deallocate (node)
      node => next
      !$omp end critical (iwrite_cr)
    end do
    !$omp critical (iwrite_cr)
    if (verbose > 0) &
      print '(1p,g14.6,3x,a,i7)', wallclock(), 'I/O done =', done, verbose
    iwrite_list%thread = -1
    !$omp end critical (iwrite_cr)
  end if
  call trace%end (itimer)
END SUBROUTINE check_real

END MODULE mpi_io_mod
