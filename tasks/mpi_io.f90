! $Id$
!*******************************************************************************
MODULE io
  USE mpi_base
  USE omp_mod, only: omp_mythread, omp_trace
  implicit none
PUBLIC
  integer, private:: handle                                                     ! file handle
  integer, private:: filetype                                                   ! file type
  integer, private:: status(MPI_STATUS_SIZE)                                    ! file status
  integer(kind=MPI_OFFSET_KIND), private:: pos                                  ! position offset into file
  integer(kind=MPI_OFFSET_KIND), private:: bytesize                             ! bytes per record
  character(len=mch), dimension(30), save:: tracing                             ! active routine
  integer, private, save:: indent=4                                             ! trace indent
  integer, private, save:: level=1                                              ! trace level
  !$omp threadprivate(indent,level,tracing)
!PUBCLIC PART
  character(len=mch):: inputfile                                                ! trace format
  logical:: do_trace=.false.                                                    ! trace execution
  logical:: do_debug=.false.                                                    ! debug execution
  real(8), save:: wc0
  integer:: verbose=0
  integer, save:: dispatch_thread, ndisp=0, nwork=0
  integer::       remove_bit=1
  integer::        retry_bit=2
  integer::        ready_bit=2**2
  integer::     internal_bit=2**3
  integer::      virtual_bit=2**4
  integer::     boundary_bit=2**5
  integer::      foreign_bit=2**6
  integer::       frozen_bit=2**7
  integer:: swap_request_bit=2**8
  integer::    root_grid_bit=2**9
  integer::         star_bit=2**10
  integer::       planet_bit=2**11
  integer::         disk_bit=2**12
  integer::     particle_bit=2**13
  integer::         leaf_bit=2**14
CONTAINS

!*******************************************************************************
SUBROUTINE init_io
  call getarg(1,inputfile)                                                      ! 1st command line argument
  if (inputfile==' ') inputfile = 'input.nml'                                   ! default = input.nml
  open (input, file=inputfile, form='formatted', status='old')                  ! open for reading
END SUBROUTINE init_io

!*******************************************************************************
SUBROUTINE print_id (id)
  implicit none
  character(len=mch) id
  character(len=mch), save:: hl= &
    '--------------------------------------------------------------------------------'
!..............................................................................
  if (id /= ' ') then
    if (master) then
      print'(a)',trim(hl)
      print 1, trim(id), wallclock()
      write(output,1) trim(id), wallclock()
    1 format(1x,a,f12.3)
    end if
    if (.not.do_debug) id = ' '
  end if
END SUBROUTINE print_id

!*******************************************************************************
SUBROUTINE trace_begin (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  !.............................................................................
  call trace_omp (id//' begin')
  tracing(level) = trim(id)
  indent = indent + 3
  level = level+1
  !write(output,*) 'indent+:', indent, level, omp_mythread()
  if (do_trace .or. do_debug) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    if (master.and.verbose>0 .or. verbose>1) then
      write(*,fmt) omp_mythread(), wallclock(), trim(id)//' begin'
    end if
    write(output,fmt) omp_mythread(), wallclock(), trim(id)//' begin'
    call flush (output)
  end if
END SUBROUTINE trace_begin

!*******************************************************************************
SUBROUTINE trace_omp (id)
  implicit none
  character(len=*) id
  logical omp_in_parallel
  !.............................................................................
  if (omp_trace) then
    write (77,'(f12.6,4x,a16,i4,2i8)') &
      wallclock()-wc0, id, omp_mythread(), ndisp, nwork
  end if
END SUBROUTINE trace_omp

!*******************************************************************************
SUBROUTINE trace (id)
  implicit none
  character(len=*) id
  character(len=mch):: fmt
  call trace_omp (id)
  if (do_trace .or. do_debug) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    if (master.and.verbose>0 .or. verbose>1) then
      write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level-1))//' '//trim(id)
    end if
    write(output,fmt) omp_mythread(), wallclock(), trim(tracing(level-1))//' '//trim(id)
    call flush (output)
  end if
END SUBROUTINE trace

!*******************************************************************************
SUBROUTINE trace_end
  implicit none
  integer i
  character(len=mch):: fmt
  level  = max(level-1,1)
  call trace_omp (trim(tracing(level))//' end')
  indent = max(indent,4)
  if (do_trace .or. do_debug) then
    write(fmt,'(a,i3.3,a)') '("trace:",i3,f12.3,',indent,'x,a)'
    if (master.and.verbose>0 .or. verbose>1) then
      write(*,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' end'
    end if
    write(output,fmt) omp_mythread(), wallclock(), trim(tracing(level))//' end'
    call flush (output)
  end if
  indent = max(indent-3,4)
  !write(output,*) 'indent-:', indent, level, omp_mythread()
END SUBROUTINE trace_end

!*******************************************************************************
SUBROUTINE file_open_mpi (filename, dd, d, offset, mode)
  implicit none
  integer:: mode
  character(len=*):: filename
  integer:: d(3)                                                                ! local dimensions
  integer:: dd(3)                                                               ! global dimensions
  integer:: offset(3)                                                           ! offset into array
!...............................................................................
  bytesize = 4_8*product(int(dd,kind=8))                                        ! record size in bytes
  CALL MPI_File_open (mpi_comm, filename, mode, MPI_INFO_NULL, handle, mpi_err) ! open file
  CALL MPI_Type_create_subarray (3, dd, d, offset, MPI_ORDER_FORTRAN, &         ! make filetype
                                 MPI_REAL, filetype, mpi_err)
  CALL MPI_Type_commit (filetype, mpi_err)                                      ! commit filetype
END SUBROUTINE file_open_mpi

!*******************************************************************************
SUBROUTINE file_openw_mpi (filename, dd, d, offset)
  implicit none
  integer, dimension(3):: dd, d, offset
  character(len=*) filename
!...............................................................................
  call file_open_mpi (filename, dd, d, offset, MPI_MODE_CREATE+MPI_MODE_RDWR)   ! open read-write
END SUBROUTINE file_openw_mpi

!*******************************************************************************
SUBROUTINE file_openr_mpi (filename, dd, d, offset)
  implicit none
  integer, dimension(3):: dd, d, offset
  character(len=*) filename
!...............................................................................
  call file_open_mpi (filename, dd, d, offset, MPI_MODE_RDONLY)                 ! open read-only
END SUBROUTINE file_openr_mpi

!*******************************************************************************
SUBROUTINE file_write_mpi (f, d, rec)
  implicit none
  integer:: d(3), rec
  real, dimension(d(1),d(2),d(3)):: f
!...............................................................................
  pos = rec*bytesize                                                            ! position in file
  CALL MPI_File_set_view (handle, pos, MPI_REAL, filetype, 'native', MPI_INFO_NULL, mpi_err)
  CALL MPI_File_write_all (handle, f , product(d), MPI_REAL, status, mpi_err)   ! collective write
END SUBROUTINE file_write_mpi

!*******************************************************************************
SUBROUTINE file_read_mpi (f, d, rec)
  implicit none
  integer:: d(3), rec
  real, dimension(d(1),d(2),d(3)):: f
!...............................................................................
  pos = rec*bytesize                                                            ! position in file
  CALL MPI_FILE_SET_VIEW (handle, pos, MPI_REAL, filetype, 'native', MPI_INFO_NULL, mpi_err)
  CALL MPI_FILE_READ_ALL(handle, f , product(d), MPI_REAL, status, mpi_err)     ! collective read
END SUBROUTINE file_read_mpi

!*******************************************************************************
SUBROUTINE file_close_mpi
  implicit none
!...............................................................................
  CALL MPI_File_close (handle, mpi_err)                                         ! close the file
  CALL MPI_Type_free  (filetype, mpi_err)                                       ! commit filetype
END SUBROUTINE file_close_mpi

END MODULE io
