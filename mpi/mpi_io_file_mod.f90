!> *****************************************************************************
!> Module for handling blocking and non-blocking MPI parallel I/O to a single file.
!> *****************************************************************************
MODULE mpi_io_file_mod
  USE mpi_mod, only: mp=>mpi
  USE omp_mod
  USE io_mod
  USE trace_mod
#ifdef MPI
  USE mpi
#endif
  implicit none
  private
  type, public:: mpi_io_file_t
    integer:: handle = 0
    integer:: comm = 0
    integer:: err = 0
    character(len=64):: filename
  contains
    procedure:: open
    procedure:: openw
    procedure:: openr
    procedure:: close
    procedure:: assert
  end type
  integer, save:: verbose=2
CONTAINS

!> =============================================================================
!> Open a file with given mode, and given local size self%lsize.  It is not
!> necessary that all threads execute this code, but at least one thread must
!> do it, before other threads can use the file.
!> =============================================================================
SUBROUTINE open (self, filename, mode)
  class(mpi_io_file_t):: self
  character(len=*) filename
  integer:: mode
  !.............................................................................
  integer:: n
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
#endif
  !-----------------------------------------------------------------------------
#ifdef MPI
  !$omp critical (mpi_io_cr)
  if (self%handle == 0) then
    self%filename = filename
    call MPI_Comm_dup (mp%comm, self%comm, self%err)
    call MPI_File_open (self%comm, filename, mode, MPI_INFO_NULL, self%handle, self%err)
    if (io%verbose > 0) then
      write(io_unit%mpi,'(a," status:",i4)') &
        ' MPI_File_open: file='//trim(filename)//'  mode='// &
        trim(merge('MPI_MODE_CREATE', &
                   '               ' ,iand(mode,MPI_MODE_CREATE) .ne. 0)) // &
        trim(merge('+MPI_MODE_RDWR' , &
                   '              '  ,iand(mode,MPI_MODE_RDWR)   .ne. 0)) // &
        trim(merge('MPI_MODE_RDONLY', &
                   '               ' ,iand(mode,MPI_MODE_RDONLY) .ne. 0)), self%err
      flush (io_unit%mpi)
    end if
  end if
  !$omp end critical (mpi_io_cr)
#endif
END SUBROUTINE open

!> =============================================================================
!> Open a file for writing, relieving the user from having to remember (and
!> having access to) the appropriate mode integers
!> =============================================================================
SUBROUTINE openw (self, filename, recl)
  class(mpi_io_file_t):: self
  character(len=*) filename
  integer, optional:: recl
  integer mode
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%openw')
#ifdef MPI
  call self%open (filename, MPI_MODE_CREATE + MPI_MODE_RDWR)
#endif
  if (present(recl)) then
    write (io%output,*) 'open:', trim(filename), recl
    open (unit=io_unit%direct, file=filename, access='direct', recl=recl, &
      status='unknown')
  end if
  call trace%end
END SUBROUTINE openw

!> =============================================================================
!> Open a file for reading, releaving the user from having to remember (and
!> having access to) the appropriate mode integers
!> =============================================================================
SUBROUTINE openr (self, filename)
  class(mpi_io_file_t):: self
  character(len=*) filename
  integer mode
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%openr')
#ifdef MPI
  call self%open (filename, MPI_MODE_RDONLY)
#endif
  call trace%end
END SUBROUTINE openr

!> =============================================================================
!> Close a file that is open for MPI parallel I/O.  It is not necessary that all
!> threads execute this code, but all threads must have finalized I/O before one
!> thread closes the file here. This may require an !$omp barrier in the calling
!> code.
!> =============================================================================
SUBROUTINE close (self)
  class(mpi_io_file_t):: self
  !-----------------------------------------------------------------------------
#ifdef MPI
  !$omp critical (close_cr)
  if (self%handle /= 0) then
    call MPI_File_closE(self%handle, self%err)
    if (io%verbose > 0 .and. io%master) then
      print'(1x,a,2i8)','mpi_io_mod::close: '//trim(self%filename)
    end if
    self%handle = 0
  end if
  !$omp end critical (close_cr)
#endif MPI
END SUBROUTINE close

!> =============================================================================
!> Assert no error
!> =============================================================================
SUBROUTINE assert (self, label, err)
  class(mpi_io_file_t):: self
  character(len=*):: label
  logical, optional:: err
  !-----------------------------------------------------------------------------
  if (present(err)) then
    if (err) then
      call mp%abort (label)
    end if
  else if (self%err /= 0) then
    print*,mp%rank,'error =',self%err
    call mp%abort (label)
  end if
END SUBROUTINE assert

END MODULE mpi_io_file_mod
