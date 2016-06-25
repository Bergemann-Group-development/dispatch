!> =============================================================================
!> $Id$
!> MPI calls related to input/output
!> =============================================================================
MODULE mpi_par_io_mod
  USE omp_mod, only: omp_mythread, omp_trace
  USE mpi_mod
  USE mpi_coords_mod
  USE mpi_timer
  implicit none
  private
  include 'mpif.h'
  integer:: mpi_err
  integer:: handle                                                     ! file handle
  integer:: filetype                                                   ! file type
  integer:: status(MPI_STATUS_SIZE)                                    ! file status
  integer(kind=MPI_OFFSET_KIND):: pos                                  ! position offset into file
  type, public, extends(mpi_t):: mpi_par_io_t
    integer:: n(3)
    integer(kind=MPI_OFFSET_KIND):: bytesize                           ! bytes per record
  contains
    procedure:: open
    procedure:: openr
    procedure:: openw
    procedure:: read
    procedure:: write
    procedure:: close
  end type
  type(mpi_par_io_t), public:: mpi_par_io
CONTAINS

!> =============================================================================
!> Open a file for MPI parallel I/O
!> =============================================================================
SUBROUTINE open (mp, filename, gn, n, offset, mode)
  implicit none
  class(mpi_par_io_t):: mp
  integer:: mode
  character(len=*):: filename
  integer, dimension(3):: offset, gn, n
!...............................................................................
  mp%bytesize = 4_8*product(int(gn,kind=8))                                  ! record size in bytes
  mp%n  =  n
  CALL MPI_File_open (mpi_comm_world, filename, mode, MPI_INFO_NULL, handle, mpi_err) ! open file
  CALL MPI_Type_create_subarray (3, gn, n, offset, MPI_ORDER_FORTRAN, &
                                 MPI_REAL, filetype, mpi_err)
  CALL MPI_Type_commit (filetype, mpi_err)                                      ! commit filetype
END SUBROUTINE open

!> =============================================================================
!> Open for reading
!> =============================================================================
SUBROUTINE openr (mp, filename, gn, n, offset)
  implicit none
  class(mpi_par_io_t):: mp
  integer, dimension(3):: offset, gn, n
  character(len=*) filename
!...............................................................................
  call mp%open (filename, gn, n, offset, MPI_MODE_RDONLY)                 ! open read-only
END SUBROUTINE openr

!> =============================================================================
SUBROUTINE openw (mp, filename, gn, n, offset)
  implicit none
  class(mpi_par_io_t):: mp
  integer, dimension(3):: offset, gn, n
  character(len=*) filename
!...............................................................................
  call mp%open (filename, gn, n, offset, MPI_MODE_CREATE+MPI_MODE_RDWR)   ! open read-write
END SUBROUTINE openw

!> =============================================================================
!> Open for writing
!> =============================================================================
SUBROUTINE write (mp, f, rec)
  implicit none
  class(mpi_par_io_t):: mp
  integer:: rec
  real, dimension(mp%n(1),mp%n(2),mp%n(3)):: f
!...............................................................................
  pos = rec*mp%bytesize                                                            ! position in file
  CALL MPI_File_set_view (handle, pos, MPI_REAL, filetype, 'native', MPI_INFO_NULL, mpi_err)
  CALL MPI_File_write_all (handle, f , product(mp%n), MPI_REAL, status, mpi_err)   ! collective write
END SUBROUTINE write

!> =============================================================================
!> Read with MPI parallell I/O
!> =============================================================================
SUBROUTINE read (mp, f, rec)
  implicit none
  class(mpi_par_io_t):: mp
  integer:: rec
  real, dimension(mp%n(1),mp%n(2),mp%n(3)):: f
!...............................................................................
  pos = rec*mp%bytesize                                                            ! position in file
  CALL MPI_FILE_SET_VIEW (handle, pos, MPI_REAL, filetype, 'native', MPI_INFO_NULL, mpi_err)
  CALL MPI_FILE_READ_ALL(handle, f , product(mp%n), MPI_REAL, status, mpi_err)     ! collective read
END SUBROUTINE read

!> =============================================================================
!> Close a file that was opened for MPI parallel I/O
!> =============================================================================
SUBROUTINE close (mp)
  implicit none
  class(mpi_par_io_t):: mp
!...............................................................................
  CALL MPI_File_close (handle, mpi_err)                                         ! close the file
  CALL MPI_Type_free  (filetype, mpi_err)                                       ! commit filetype
END SUBROUTINE close

END MODULE mpi_par_io_mod
