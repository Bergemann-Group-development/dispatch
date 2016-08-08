!> *****************************************************************************
!> Module for handling non-blocking MPI parallel I/O to a single file.  We have
!> no interest in making this general, since we do not plan to write to several
!> files at the same time.  The module is initialized by specifying the size
!> and relative position of the local chunk relative to the whole file. A chunk
!> could be just a patch, or a single chunk of data per rank (but this would
!> require that a rank held a cubic region, which would be too restrictive).
!>
!> As is, this module is oriented towards handling unigrid files (e.g. the ones
!> used by the original STAGGER code), but the method can be generalized to
!> handle any collection of patches, by using the patch%id number to map each
!> patch into the file.  If the patch sizes are fixed, this is a unique mapping,
!> and if each snapshot is written into a separate file this is also sufficiently
!> general.
!> *****************************************************************************
MODULE mpi_io_mod
  USE mpi_mod, only: mp=>mpi
  USE io_mod
  USE mpi
  implicit none
  !include 'mpif.h'
  private
  type, public:: mpi_io_t
    integer:: gsize(3)
    integer:: lsize(3)
    integer:: start(3)
    integer:: handle
    integer:: filetype
    integer:: req
    integer:: count
    integer:: err
    integer:: status(MPI_STATUS_SIZE)
    character(len=64):: filename
  contains
    procedure:: init
    procedure:: open
    procedure:: openw
    procedure:: openr
    procedure:: close
    procedure:: write
    procedure:: read
    procedure:: test
    procedure:: wait
    procedure:: assert
  end type
  character(len=32), save:: fmt='(1x,a,4i8,1p,5e12.3)'
CONTAINS

!> =============================================================================
!> Initialize the local size self%lsize, global size self%gsize, and start
!> address self%start for MPI parallel I/O
!> =============================================================================
SUBROUTINE init (self, gsize, lsize, start)
  class(mpi_io_t):: self
  integer, dimension(3):: gsize, lsize, start
  !-----------------------------------------------------------------------------
  self%gsize = gsize
  self%lsize = lsize
  self%start = start
  self%count = product(self%lsize)
END SUBROUTINE init

!> =============================================================================
!> Open a file with given mode, and given local size self%lsize, global size
!> self%gsize, and start address self%start; these must be initialized first
!> =============================================================================
SUBROUTINE open (self, filename, mode)
  class(mpi_io_t):: self
  character(len=*) filename
  integer mode
  !-----------------------------------------------------------------------------
  self%filename = filename
  CALL MPI_FILE_OPEN (mp%comm, filename, mode, MPI_INFO_NULL, self%handle, &
                      self%err)
  if (io%verbose>0 .and. io%master) then
    print'(a)',' MPI_FILE_OPEN: file='//trim(filename)//'  mode='// &
      trim(merge(   'MPI_MODE_CREATE', &
                    '               ',iand(mode,MPI_MODE_CREATE) .ne. 0)) // &
      trim(merge('+MPI_MODE_RDWR'  , &
                 '              '    ,iand(mode,MPI_MODE_RDWR)   .ne. 0)) // &
      trim(merge('MPI_MODE_RDONLY'    , &
                 '               '   ,iand(mode,MPI_MODE_RDONLY) .ne. 0))
  endif
  call self%assert ('mpi_io_mod::open MPI_FILE_OPEN')
  !
  CALL MPI_TYPE_CREATE_SUBARRAY (3, self%gsize, self%lsize, self%start, &
    MPI_ORDER_FORTRAN, MPI_REAL, self%filetype, self%err)
  call self%assert ('mpi_io_mod::open MPI_TYPE_CREATE_SUBARRAY')
  !
  CALL MPI_TYPE_COMMIT (self%filetype, self%err)
  call self%assert ('mpi_io_mod::open MPI_TYPE_COMMIT')
  !
  if (io%verbose > 0 .and. io%master) then
    print'(1x,a,2i8)','mpi_io_mod::open: '//trim(filename), mp%rank, self%err
  end if
END SUBROUTINE open

!> =============================================================================
!> Open a file for writing, releaving the user from having to remember (and
!> having access to) the appropriate mode integers
!> =============================================================================
SUBROUTINE openw (self, filename)
  class(mpi_io_t):: self
  character(len=*) filename
  integer mode
  !-----------------------------------------------------------------------------
  call self%open (filename, MPI_MODE_CREATE + MPI_MODE_RDWR)
END SUBROUTINE openw

!> =============================================================================
!> Open a file for reading, releaving the user from having to remember (and
!> having access to) the appropriate mode integers
!> =============================================================================
SUBROUTINE openr (self, filename)
  class(mpi_io_t):: self
  character(len=*) filename
  integer mode
  !-----------------------------------------------------------------------------
  call self%open (filename, MPI_MODE_RDONLY)
END SUBROUTINE openr

!> =============================================================================
!> Close a file that is open for MPI parallel I/O
!> =============================================================================
SUBROUTINE close (self)
  class(mpi_io_t):: self
  !-----------------------------------------------------------------------------
  CALL MPI_FILE_CLOSE(self%handle, self%err)
  if (io%verbose > 0 .and. io%master) then
    print'(1x,a,2i8)','mpi_io_mod::close: '//trim(self%filename)
  end if
END SUBROUTINE close

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE write (self, f, rec)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer:: rec
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  pos = 4*(rec-1)*product(self%gsize)
  CALL MPI_FILE_SET_VIEW (self%handle, pos, MPI_REAL, self%filetype, 'native', &
                          MPI_INFO_NULL, self%err)
  call self%assert ('mpi_io_mod::write MPI_FILE_SET_VIEW')
  !
  CALL MPI_FILE_IWRITE (self%handle, f, self%count, MPI_REAL, self%req, self%err)
  call self%assert ('mpi_io_mod::write MPI_FILE_IWRITE')
  !
  if (io%verbose > 0 .and. io%master) then
    print fmt,'mpi_io_mod::write:', &
      mp%rank, self%count, rec, pos, f(3:5,3,3)
  end if
END SUBROUTINE write

!> =============================================================================
!> Read a chunk (self%lsize) from the open file
!> =============================================================================
SUBROUTINE read (self, f, rec)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer rec
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  pos = 4*(rec-1)*product(self%gsize)
  CALL MPI_FILE_SET_VIEW (self%handle, pos, MPI_REAL, self%filetype, 'native', &
                          MPI_INFO_NULL, self%err)
  call self%assert ('mpi_io_mod::read MPI_FILE_SET_VIEW')
  !
  CALL MPI_FILE_IREAD (self%handle, f , self%count, MPI_REAL, self%req, self%err)
  call self%assert ('mpi_io_mod::read MPI_FILE_IREAD')
  !
  if (io%verbose > 0 .and. io%master) then
    print fmt,'mpi_io_mod::read:', mp%rank, self%count, rec, pos, f(3:5,3,3)
  end if
END SUBROUTINE read

!> =============================================================================
!> Test if an I/O operation is complete
!> =============================================================================
FUNCTION test (self) RESULT (flag)
  class(mpi_io_t):: self
  logical:: flag
  !-----------------------------------------------------------------------------
  call MPI_TEST (self%req, flag, self%status, self%err)
  call self%assert ('mpi_io_mod::read MPI_TEST')
END FUNCTION test

!> =============================================================================
!> Wait for an MPI I/O operation to complete
!> =============================================================================
SUBROUTINE wait (self)
  class(mpi_io_t):: self
  call MPI_WAIT (self%req, self%status, self%err)
  call self%assert ('mpi_io_mod::read MPI_WAIT')
END SUBROUTINE wait

!> =============================================================================
!> Assert no error
!> =============================================================================
SUBROUTINE assert (self, label)
  class(mpi_io_t):: self
  character(len=*):: label
  if (self%err /= 0) then
    print*,mp%rank,'ERROR: '//trim(label)//' error =',self%err
    call mp%abort (label)
  end if
END SUBROUTINE assert
END MODULE mpi_io_mod
