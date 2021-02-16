!*******************************************************************************
!> Module for handling blocking and non-blocking MPI parallel I/O to a single file.
!*******************************************************************************
MODULE mpi_file_mod
  USE mpi_mod, only: mp=>mpi
  USE omp_mod
  USE omp_timer_mod
  USE io_mod
  USE trace_mod
  USE mpi_lib
  implicit none
  private
  type, public:: mpi_file_t
    integer:: handle = 0
    integer:: comm = 0
    integer:: mode = 0
    integer:: err = 0
    integer(8):: size = 0_8
    integer(8):: offset = 0_8
    logical:: closed=.true.
    character(len=64):: name
  contains
    procedure:: request_open
    procedure:: request_close
    procedure:: open_and_close
    procedure:: open
    procedure:: openw
    procedure:: openr
    procedure:: close
    procedure, private:: write1
    procedure, private:: write5
    procedure, private:: writes
    generic, public :: write => write1, write5, writes
    procedure, private:: read1
    generic, public :: read => read1
    procedure:: assert
  end type
#ifndef MPI
  integer:: MPI_REAL=1, MPI_STATUS_SIZE=8
#endif
  integer, save:: verbose=2
CONTAINS

!===============================================================================
!> Request that a file is opened, when convenient
!===============================================================================
SUBROUTINE request_open (self, filename)
  class(mpi_file_t):: self
  character(len=*) filename
  character(len=120):: id = &
  '$Id$'
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_file_t%request_open')
  call trace%print_id (id)
#ifdef MPI
  if (self%closed) then
    self%mode = MPI_MODE_CREATE + MPI_MODE_RDWR
    self%name = filename
    self%closed = .false.
  end if
#endif
  call trace%end()
END SUBROUTINE request_open

!===============================================================================
!> Request that a file is closed, when convenient
!===============================================================================
SUBROUTINE request_close (self)
  class(mpi_file_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_file_t%request_close')
  self%closed = .true.
  call trace%end()
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
      'mpi_file_t%open_and_close: closed '//trim(self%name)
    flush (io_unit%mpi)
  else if (.not. self%closed .and. self%handle==0) then
    call self%openw (self%name)
    write(io_unit%mpi,'(g12.5,i4,2x,a)') wallclock(), omp%thread, &
      'mpi_file_t%open_and_close: opened '//trim(self%name)
    flush (io_unit%mpi)
  end if
  call trace%end()
END SUBROUTINE open_and_close

!===============================================================================
!> Open a file with given mode, and given local size self%lsize.  It is not
!> necessary that all threads execute this code, but at least one thread must
!> do it, before other threads can use the file.
!===============================================================================
SUBROUTINE open (self, filename, mode)
  class(mpi_file_t):: self
  character(len=*), optional:: filename
  integer, optional:: mode
  !.............................................................................
  integer:: n, mode_l
#ifdef MPI
  integer(kind=MPI_OFFSET_KIND):: pos
  !-----------------------------------------------------------------------------
  if (present(filename)) then
    self%name = filename
  end if
  if (present(mode)) then
    mode_l = mode
  else
    mode_l = MPI_MODE_RDWR + MPI_MODE_CREATE
  end if
  !call mp%barrier ('mpi_file_t%open')
  if (self%handle == 0) then
    if (mp%mode == MPI_THREAD_MULTIPLE) then
      call MPI_File_open &
        (mp%comm, self%name, mode_l, MPI_INFO_NULL, self%handle, self%err)
      if (self%err /= 0) then
        write (stderr,*) self%err, ' error when trying to open ', trim(self%name)
        call io%abort
      end if
    else
      call MPI_File_open &
        (mp%comm, self%name, mode_l, MPI_INFO_NULL, self%handle, self%err)
      if (self%err /= 0) then
        write (stderr,*) self%err, ' error when trying to open ', trim(self%name)
        call io%abort
      end if
    end if
      if (io%verbose > 0) then
        write(io_unit%mpi,'(a," status:",i4)') &
          ' MPI_File_open: file='//trim(self%name)//'  mode='// &
          trim(merge('MPI_MODE_CREATE', &
            '               ' ,iand(mode_l,MPI_MODE_CREATE) .ne. 0)) // &
          trim(merge('+MPI_MODE_RDWR' , &
            '              '  ,iand(mode_l,MPI_MODE_RDWR)   .ne. 0)) // &
          trim(merge('MPI_MODE_RDONLY', &
            '               ' ,iand(mode_l,MPI_MODE_RDONLY) .ne. 0)), &
          self%err
      flush (io_unit%mpi)
    end if
  end if
  self%closed = .false.
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
  call trace%begin ('mpi_file_t%openw')
#ifdef MPI
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call self%open (filename, MPI_MODE_CREATE + MPI_MODE_RDWR)
  else
    !$omp critical (mpi_cr)
    call self%open (filename, MPI_MODE_CREATE + MPI_MODE_RDWR)
    !$omp end critical (mpi_cr)
  end if
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
  call trace%begin ('mpi_file_t%openr')
#ifdef MPI
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call self%open (filename, MPI_MODE_RDONLY)
  else
    !$omp critical (mpi_cr)
    call self%open (filename, MPI_MODE_RDONLY)
    !$omp end critical (mpi_cr)
  end if
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
        'mpi_file_t%close: '//trim(self%name)//'  err:', self%err
      flush (io_unit%mpi)
    end if
    self%handle = 0
  end if
#endif MPI
  call trace%end()
END SUBROUTINE close

!===============================================================================
!> Write out a buffer at a given offset into the file
!===============================================================================
SUBROUTINE write1 (self, off, size, buffer)
  class(mpi_file_t):: self
  integer(8):: off
  integer:: size, st(MPI_STATUS_SIZE), err
  integer:: buffer(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%write1')
#ifdef MPI
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_File_write_at (self%handle, off, buffer, size, MPI_REAL, st, err)
    call io%assert (err==0, 'mpi_file_t%read1 non-zero status')
  else
    !$omp critical (mpi_cr)
    call MPI_File_write_at (self%handle, off, buffer, size, MPI_REAL, st, err)
    call io%assert (err==0, 'mpi_file_t%read1 non-zero status')
    !$omp end critical (mpi_cr)
  end if
#endif MPI
  call trace%end ()
END SUBROUTINE write1

!===============================================================================
!> Read a buffer at a given offset into the file
!===============================================================================
SUBROUTINE read1 (self, off, size, buffer)
  class(mpi_file_t):: self
  integer(8):: off
  integer:: size, st(MPI_STATUS_SIZE), err
  integer:: buffer(:)
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%read1')
#ifdef MPI
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_File_read_at (self%handle, off, buffer, size, MPI_REAL, st, err)
    call io%assert (err==0, 'mpi_file_t%read1 non-zero status')
  else
    !$omp critical (mpi_cr)
    call MPI_File_read_at (self%handle, off, buffer, size, MPI_REAL, st, err)
    call io%assert (err==0, 'mpi_file_t%read1 non-zero status')
    !$omp end critical (mpi_cr)
  end if
#endif MPI
  call trace%end ()
END SUBROUTINE read1

!===============================================================================
!> Write out a buffer at a given offset into the file
!===============================================================================
SUBROUTINE write5 (self, off, size, buffer)
  class(mpi_file_t):: self
  integer(8):: off
  integer:: size, st(MPI_STATUS_SIZE), err
  real:: buffer(:,:,:,:,:)
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%write5')
#ifdef MPI
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_File_write_at (self%handle, off, buffer, size, MPI_REAL, st, err)
  else
    !$omp critical (mpi_cr)
    call MPI_File_write_at (self%handle, off, buffer, size, MPI_REAL, st, err)
    !$omp end critical (mpi_cr)
  end if
#endif MPI
  call trace%end ()
END SUBROUTINE write5

!===============================================================================
!> Write out a buffer at a given offset into the file
!===============================================================================
SUBROUTINE writes (self, buffer, off, size)
  class(mpi_file_t):: self
  integer(8), optional:: off
  integer, optional:: size
  integer:: err, st(MPI_STATUS_SIZE)
  class(*), pointer:: buffer
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_io_t%writes')
  if (present(off)) then
    self%offset = off
  end if
  if (present(size)) then
    self%size = size
  end if
#ifdef MPI
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_File_write_at (self%handle, off, buffer, size, MPI_REAL, st, err)
  else
    !$omp critical (mpi_cr)
    call MPI_File_write_at (self%handle, off, buffer, size, MPI_REAL, st, err)
    !$omp end critical (mpi_cr)
  end if
#endif MPI
  call trace%end ()
END SUBROUTINE writes

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
