!> =============================================================================
!> $Id$
!> MPI calls related to input/output
!> =============================================================================
MODULE mpi_par_io
  USE omp_mod, only: omp_mythread, omp_trace
  USE mpi_basic
  USE mpi_timer
  USE mpi_coords
  USE mpi
  implicit none
PUBLIC
  integer, private:: mpi_err
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
  real(8), save:: wc0
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
  type, extends(mpi_coord_t):: mpi_par_io_t
  contains
    procedure:: file_open
    procedure:: file_openr
    procedure:: file_openw
    procedure:: file_read
    procedure:: file_write
    procedure:: file_close
  end type
CONTAINS

!> =============================================================================
!> Open a file for MPI parallel I/O
!> =============================================================================
SUBROUTINE file_open (mp, filename, dd, d, offset, mode)
  implicit none
  class(mpi_par_io_t):: mp
  integer:: mode
  character(len=*):: filename
  integer:: d(3)                                                                ! local dimensions
  integer:: dd(3)                                                               ! global dimensions
  integer:: offset(3)                                                           ! offset into array
!...............................................................................
  bytesize = 4_8*product(int(dd,kind=8))                                        ! record size in bytes
  CALL MPI_File_open (mpi_comm_world, filename, mode, MPI_INFO_NULL, handle, mpi_err) ! open file
  CALL MPI_Type_create_subarray (3, dd, d, offset, MPI_ORDER_FORTRAN, &         ! make filetype
                                 MPI_REAL, filetype, mpi_err)
  CALL MPI_Type_commit (filetype, mpi_err)                                      ! commit filetype
END SUBROUTINE file_open

!> =============================================================================
SUBROUTINE file_openw (mp, filename, dd, d, offset)
  implicit none
  class(mpi_par_io_t):: mp
  integer, dimension(3):: dd, d, offset
  character(len=*) filename
!...............................................................................
  call mp%file_open (filename, dd, d, offset, MPI_MODE_CREATE+MPI_MODE_RDWR)   ! open read-write
END SUBROUTINE file_openw

!> =============================================================================
!> Open for reading
!> =============================================================================
SUBROUTINE file_openr (mp, filename, dd, d, offset)
  implicit none
  class(mpi_par_io_t):: mp
  integer, dimension(3):: dd, d, offset
  character(len=*) filename
!...............................................................................
  call mp%file_open (filename, dd, d, offset, MPI_MODE_RDONLY)                 ! open read-only
END SUBROUTINE file_openr

!> =============================================================================
!> Open for writing
!> =============================================================================
SUBROUTINE file_write (mp, f, d, rec)
  implicit none
  class(mpi_par_io_t):: mp
  integer:: d(3), rec
  real, dimension(d(1),d(2),d(3)):: f
!...............................................................................
  pos = rec*bytesize                                                            ! position in file
  CALL MPI_File_set_view (handle, pos, MPI_REAL, filetype, 'native', MPI_INFO_NULL, mpi_err)
  CALL MPI_File_write_all (handle, f , product(d), MPI_REAL, status, mpi_err)   ! collective write
END SUBROUTINE file_write

!> =============================================================================
!> Read with MPI parallell I/O
!> =============================================================================
SUBROUTINE file_read (mp, f, d, rec)
  implicit none
  class(mpi_par_io_t):: mp
  integer:: d(3), rec
  real, dimension(d(1),d(2),d(3)):: f
!...............................................................................
  pos = rec*bytesize                                                            ! position in file
  CALL MPI_FILE_SET_VIEW (handle, pos, MPI_REAL, filetype, 'native', MPI_INFO_NULL, mpi_err)
  CALL MPI_FILE_READ_ALL(handle, f , product(d), MPI_REAL, status, mpi_err)     ! collective read
END SUBROUTINE file_read

!> =============================================================================
!> Close a file that was opened for MPI parallel I/O
!> =============================================================================
SUBROUTINE file_close (mp)
  implicit none
  class(mpi_par_io_t):: mp
!...............................................................................
  CALL MPI_File_close (handle, mpi_err)                                         ! close the file
  CALL MPI_Type_free  (filetype, mpi_err)                                       ! commit filetype
END SUBROUTINE file_close

END MODULE mpi_par_io
