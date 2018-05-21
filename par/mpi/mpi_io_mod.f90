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
  USE omp_mod
  USE trace_mod
  USE mpi_mod, only: mp=>mpi
  USE mpi_io_file_mod
  USE mpi
  implicit none
  private
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
    procedure:: check_init
    procedure:: write
    procedure:: iwrite
    procedure:: read
    procedure:: test
    procedure:: wait
    procedure:: assert
  end type
  character(len=32), save:: fmt='(1x,a,4i8,1p,5e12.3)'
  integer, save:: verbose=0
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
SUBROUTINE init (self, rec_words, chunk_words)
  class(mpi_io_t):: self
  integer(8), optional:: rec_words, chunk_words
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%init')
  if (present(rec_words)) then
    self%rec_words = rec_words
    if (verbose > 0) print *, 'mpi_io_t%init: setting rec_words =', rec_words
  end if
  if (present(chunk_words)) then
    self%chunk_words = chunk_words
    if (verbose > 0) print *, 'mpi_io_t%init: setting chunk_words =', chunk_words
  end if
  call trace%end
END SUBROUTINE init

!> =============================================================================
!> Initialize the local size self%lsize, global size self%gsize
!> =============================================================================
SUBROUTINE check_init (self)
  class(mpi_io_t):: self
  !-----------------------------------------------------------------------------
  call self%assert ('rec_words is not positive', err=self%rec_words <= 0)
  call self%assert ('chunk_words is not positive', err=self%chunk_words <= 0)
END SUBROUTINE check_init

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE write (self, f, rec, id)
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
    print '(a,3i8,i15,2f12.0)', ' mpi_io_t%write: rec, id, words, pos =', &
    rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !$omp critical (mpi_io_cr)
  call MPI_file_write_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                          self%status, self%err)
  !$omp end critical (mpi_io_cr)
  call self%assert ('mpi_io_mod::write MPI_File_write_at '//trim(self%file%filename))
  call trace%end
#endif
END SUBROUTINE write

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE iwrite (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer:: rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%iwrite')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  if (verbose>1) &
    print '(a,3i8,i15,2f12.0)', ' mpi_io_t%iwrite: rec, id, words, pos =', &
    rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !$omp critical (mpi_io_cr)
  if (allocated(self%buf)) deallocate (self%buf)
  allocate (self%buf(size(f,1),size(f,2),size(f,3)))
  self%buf = f
  call MPI_file_iwrite_at (self%file%handle, pos, self%buf, self%chunk_words, &
                           MPI_REAL, self%req, self%err)
  !$omp end critical (mpi_io_cr)
  call self%assert ('mpi_io_mod::write MPI_File_iwrite_at '//trim(self%file%filename))
  call trace%end
#endif
END SUBROUTINE iwrite

!> =============================================================================
!> Read a chunk (self%lsize) from the open file
!> =============================================================================
SUBROUTINE read (self, f, rec, id)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer rec, id
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  if (rec < 1) return
  call trace%begin ('mpi_io_t%read')
  call self%check_init
  pos = 4_8*((rec-1_8)*self%rec_words + (id-1_8)*self%chunk_words)
  !$omp critical (mpi_io_cr)
  call MPI_File_read_at (self%file%handle, pos, f, self%chunk_words, MPI_REAL, &
                         self%status, self%err)
  if (verbose>1) &
    print '(a,3i8,i15,2f12.0)', ' mpi_io_t%read: rec, id, words, pos =', &
    rec, id, self%chunk_words, pos, minval(f), maxval(f)
  !$omp end critical (mpi_io_cr)
  call trace%end
#endif
END SUBROUTINE read

!> =============================================================================
!> Test if an I/O operation is complete
!> =============================================================================
FUNCTION test (self) RESULT (flag)
  class(mpi_io_t):: self
  logical:: flag
  !-----------------------------------------------------------------------------
#ifdef MPI
  call trace%begin ('mpi_io_t%test')
  call MPI_Test (self%req, flag, self%status, self%err)
  call self%assert ('mpi_io_mod::read MPI_Test')
  if (flag) then
    if (allocated(self%buf)) deallocate (self%buf)
    self%req = 0
  end if
  call trace%end
#endif
END FUNCTION test

!> =============================================================================
!> Wait for an MPI I/O operation to complete
!> =============================================================================
SUBROUTINE wait (self)
  class(mpi_io_t):: self
  !-----------------------------------------------------------------------------
#ifdef MPI
  if (self%req==0) return
  call trace%begin ('mpi_io_t%wait')
  call MPI_Wait (self%req, self%status, self%err)
  call self%assert ('mpi_io_mod::read MPI_Wait')
  deallocate (self%buf)
  self%req = 0
  call trace%end
#endif
END SUBROUTINE wait

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

END MODULE mpi_io_mod
