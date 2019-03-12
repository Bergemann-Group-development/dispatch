!*******************************************************************************
!> Module for handling blocking and non-blocking MPI parallel I/O to a single file.
!*******************************************************************************
MODULE mpi_file_mod
  USE mpi_mod, only: mp=>mpi
  USE omp_mod
  USE omp_timer_mod
  USE io_mod
  USE trace_mod
#ifdef MPI
  USE mpi
#endif
  implicit none
  private
  type, public:: mpi_file_t
    integer:: handle = 0
    integer:: comm = 0
    integer:: mode = 0
    integer:: err = 0
    logical:: closed=.true.
    character(len=64):: filename
  contains
    procedure:: request_open
    procedure:: request_close
    procedure:: open_and_close
    procedure:: open
    procedure:: openw
    procedure:: openr
    procedure:: close
    procedure:: assert
  end type
  integer, save:: verbose=2
CONTAINS

!===============================================================================
!> Request that a file is opened, when convenient
!===============================================================================
SUBROUTINE request_open (self, filename)
  class(mpi_file_t):: self
  character(len=*) filename
  character(len=120):: id = &
  '$Id$ mpi/mpi_file_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_file_t%request_open')
  call trace%print_id (id)
  if (self%closed) then
    self%mode = MPI_MODE_CREATE + MPI_MODE_RDWR
    self%filename = filename
    self%closed = .false.
  end if
  call trace%end
END SUBROUTINE request_open

!===============================================================================
!> Request that a file is closed, when convenient
!===============================================================================
SUBROUTINE request_close (self)
  class(mpi_file_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_file_t%request_close')
  self%closed = .true.
  call trace%end
END SUBROUTINE request_close

!===============================================================================
!> Make a files actual status match the requested status
!===============================================================================
SUBROUTINE open_and_close (self)
  class(mpi_file_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_file_t%open_and_close',2)
  if (self%closed .and. self%handle/=0) then
    call self%close ()
    write(io_unit%mpi,'(g12.5,i4,2x,a)') wallclock(), omp%thread, &
      'mpi_file_t%open_and_close: closed '//trim(self%filename)
    flush (io_unit%mpi)
  else if (.not. self%closed .and. self%handle==0) then
    call self%openw (self%filename)
    write(io_unit%mpi,'(g12.5,i4,2x,a)') wallclock(), omp%thread, &
      'mpi_file_t%open_and_close: opened '//trim(self%filename)
    flush (io_unit%mpi)
  end if
  call trace%end
END SUBROUTINE open_and_close

!===============================================================================
!> Open a file with given mode, and given local size self%lsize.  It is not
!> necessary that all threads execute this code, but at least one thread must
!> do it, before other threads can use the file.
!===============================================================================
SUBROUTINE open (self, filename, mode)
  class(mpi_file_t):: self
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
    if (mp%mode == MPI_THREAD_MULTIPLE) then
      call MPI_File_open (mp%comm, filename, mode, MPI_INFO_NULL, self%handle, self%err)
    else
      !$omp critical (mpi_cr)
      call MPI_File_open (mp%comm, filename, mode, MPI_INFO_NULL, self%handle, self%err)
      !$omp end critical (mpi_cr)
    end if
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
  self%closed = .false.
  !$omp end critical (mpi_io_cr)
#endif
END SUBROUTINE open

!===============================================================================
!> Open a file for writing, relieving the user from having to remember (and
!> having access to) the appropriate mode integers
!===============================================================================
SUBROUTINE openw (self, filename, recl)
  class(mpi_file_t):: self
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
  call trace%end()
END SUBROUTINE openw

!===============================================================================
!> Open a file for reading, releaving the user from having to remember (and
!> having access to) the appropriate mode integers
!===============================================================================
SUBROUTINE openr (self, filename)
  class(mpi_file_t):: self
  character(len=*) filename
  integer mode
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%openr')
#ifdef MPI
  call self%open (filename, MPI_MODE_RDONLY)
#endif
  call trace%end()
END SUBROUTINE openr

!===============================================================================
!> Close a file that is open for MPI parallel I/O.  It is not necessary that all
!> threads execute this code, but all threads must have finalized I/O before one
!> thread closes the file here. This may require an !$omp barrier in the calling
!> code.
!===============================================================================
SUBROUTINE close (self)
  class(mpi_file_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_file_t%close')
#ifdef MPI
  !$omp critical (close_cr)
  if (self%handle /= 0) then
    if (mp%mode == MPI_THREAD_MULTIPLE) then
      call MPI_File_close (self%handle, self%err)
    else
      !$omp critical (mpi_cr)
      call MPI_File_close (self%handle, self%err)
      !$omp end critical (mpi_cr)
    end if
    if (io%verbose > 0) then
      write (io_unit%mpi,'(1x,a,i8)') &
        'mpi_io_t%close: '//trim(self%filename)//'  err:', self%err
      flush (io_unit%mpi)
    end if
    self%handle = 0
  end if
  !$omp end critical (close_cr)
#endif MPI
  call trace%end
END SUBROUTINE close

!===============================================================================
!> Assert no error
!===============================================================================
SUBROUTINE assert (self, label, err)
  class(mpi_file_t):: self
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

END MODULE mpi_file_mod
