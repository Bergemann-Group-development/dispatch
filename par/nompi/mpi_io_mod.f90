!> *****************************************************************************
!> Module for handling non-blocking MPI parallel I/O to a single file. The module
!> is initialized by specifying the size (which can be 1- or 3-dimensional)
!> and relative position of the local chunk relative to the whole file. A chunk
!> could be just a patch, or a single chunk of data per rank (but this would
!> require that a rank held a cubic region, which would be too restrictive).
!>
!> The module can be used to handle any collection of patches, by using the
!> patch%id number to map each patch into the file.  If the patch sizes are fixed,
!> this is a unique mapping, and if each snapshot is written into a separate file
!> this is also sufficiently general.  The header part of the file may be written
!> to a separate file, as demonstrated in tests/mpi_io/test.f90.
!> *****************************************************************************
MODULE mpi_io_mod
  USE mpi_mod, only: mp=>mpi
  USE omp_mod
  USE io_mod
  USE mpi
  implicit none
  private
  type, public:: mpi_io_t
    integer:: n
    integer:: gsize(3) = 1
    integer:: lsize(3) = 1
    integer:: start(3) = 0
    integer:: handle
    integer:: filetype
    integer:: req
    integer:: count
    integer:: err
    integer:: status
    character(len=64):: filename
  contains
    procedure:: init
    procedure:: open
    procedure:: openw
    procedure:: openr
    procedure:: close
    generic, public:: write => write1, write3
    procedure, private::       write1
    procedure, private::       write3
    generic, public:: read => read1, read3
    procedure, private::      read1
    procedure, private::      read3
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
  integer, dimension(:):: gsize, lsize, start
  integer:: n
  !-----------------------------------------------------------------------------
  n = size(gsize,1)
  self%n = n
  self%gsize(1:n) = gsize
  self%lsize(1:n) = lsize
  self%start(1:n) = start
  self%count = product(lsize)
END SUBROUTINE init

!> =============================================================================
!> Open a file with given mode, and given local size self%lsize, global size
!> self%gsize, and start address self%start; these must be initialized first
!> =============================================================================
SUBROUTINE open (self, filename, mode)
  class(mpi_io_t):: self
  character(len=*) filename
  integer:: mode
  integer:: n
  !-----------------------------------------------------------------------------
  self%filename = filename
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
  call self%open (filename, 0)
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
  call self%open (filename, 0)
END SUBROUTINE openr

!> =============================================================================
!> Close a file that is open for MPI parallel I/O
!> =============================================================================
SUBROUTINE close (self)
  class(mpi_io_t):: self
  !-----------------------------------------------------------------------------

  if (io%verbose > 0 .and. io%master) then
    print'(1x,a,2i8)','mpi_io_mod::close: '//trim(self%filename)
  end if
END SUBROUTINE close

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE write3 (self, f, rec)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer:: rec
  integer(kind=4):: pos
  !-----------------------------------------------------------------------------
  pos = 4*(rec-1)*product(self%gsize)

  !
  if (io%verbose==1 .and. io%master) then
    print fmt,'mpi_io_mod::write: '//trim(self%filename), &
      mp%rank, self%count, rec, pos, f(3:5,3,3)
  else if (io%verbose>1) then
    print fmt,'mpi_io_mod::read: '//trim(self%filename), &
      mp%rank, omp%thread, rec, pos, f(3:5,3,3)
  end if
END SUBROUTINE write3

!> =============================================================================
!> Write a chunk (self%lsize) out to the open file
!> =============================================================================
SUBROUTINE write1 (self, f, rec)
  class(mpi_io_t):: self
  real(kind=4), dimension(:):: f
  integer:: rec
  integer(kind=4):: pos
  !-----------------------------------------------------------------------------
  pos = 4*(rec-1)*self%gsize(1)

  !
  if (io%verbose==1 .and. io%master) then
    print fmt,'mpi_io_mod::write: '//trim(self%filename), &
      mp%rank, self%count, rec, pos, f(3:5)
  else if (io%verbose>1) then
    print fmt,'mpi_io_mod::read: '//trim(self%filename), &
      mp%rank, omp%thread, rec, pos, f(3:5)
  end if
END SUBROUTINE write1

!> =============================================================================
!> Read a chunk (self%lsize) from the open file
!> =============================================================================
SUBROUTINE read3 (self, f, rec)
  class(mpi_io_t):: self
  real(kind=4), dimension(:,:,:):: f
  integer rec
  integer(kind=4):: pos
  !-----------------------------------------------------------------------------
  pos = 4*(rec-1)*product(self%gsize)

  !
  if (io%verbose==1 .and. io%master) then
    print fmt,'mpi_io_mod::read: '//trim(self%filename), &
      mp%rank, self%count, rec, pos, f(3:5,3,3)
  else if (io%verbose>1) then
    print fmt,'mpi_io_mod::read: '//trim(self%filename), &
      mp%rank, omp%thread, rec, pos, f(3:5,3,3)
  end if
END SUBROUTINE read3

!> =============================================================================
!> Read a chunk (self%lsize) from the open file
!> =============================================================================
SUBROUTINE read1 (self, f, rec)
  class(mpi_io_t):: self
  real(kind=4), dimension(:):: f
  integer rec
  integer(kind=4):: pos
  !-----------------------------------------------------------------------------
  pos = 4*(rec-1)*self%gsize(1)

  !
  if (io%verbose==1 .and. io%master) then
    print fmt,'mpi_io_mod::read: '//trim(self%filename), &
      mp%rank, self%count, rec, pos, f(3:5)
  else if (io%verbose>1) then
    print fmt,'mpi_io_mod::read: '//trim(self%filename), &
      mp%rank, omp%thread, rec, pos, f(3:5)
  end if
END SUBROUTINE read1

!> =============================================================================
!> Test if an I/O operation is complete
!> =============================================================================
FUNCTION test (self) RESULT (flag)
  class(mpi_io_t):: self
  logical:: flag
  !-----------------------------------------------------------------------------

END FUNCTION test

!> =============================================================================
!> Wait for an MPI I/O operation to complete
!> =============================================================================
SUBROUTINE wait (self)
  class(mpi_io_t):: self

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
