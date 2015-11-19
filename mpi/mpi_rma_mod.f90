!> =============================================================================
!> $Id$
!> MPI remote-memory-access (RMA) calls, also referred to as one-sided MPI
!> 
!> To make use of this module, make sure it is compiled together with your code,
!> using either the Makefile and the 'make' command, or compiling manually, with
!> 
!> mpi = ../../../../../mpi                     # or wherever
!> mpifort -c $mpi/mpi_mod.f90                  # compiler mpi_mod.f90
!> mpifort -c $mpi/mpi_coords.f90               # add any other module your need
!> mpifort -c $mpi/mpi_rma.f90                  # add any other module your need
!> mpifort -c your_code.f90                     # compile your code
!> mpifort *.o -o your_code.x                   # link together into your_code.x 
!>
!> In your_code.f90 you add lines such as:
!>
!> USE mpi_coords_mod                           ! makes the mpi_coords% object available
!> USE mpi_rma_mod                              ! makes the mpi_rma% object available
!> implicit none
!> real, dimension(:), allocatable: send        ! send buffer
!> real, dimension(:), allocatable: recv        ! receive buffer
!> integer:: wins, winr                         ! windws into send and receive buffers
!> ...
!> call mpi_coords%init (dims=[mx,my,mz])       ! initializes the mpi_coords%object
!> call mpi_coords%print                        ! prints a summary of the MPI-geometry
!> ...
!> allocate (send(nw), recv(nw))                ! allocate buffers
!> call mpi_rma%win_reals (nw, send, wins)      ! open window into send buffer
!> call mpi_rma%win_reals (nw, recv, winr)      ! open window into recv buffer
!> ...
!> call mpi_rma%put_reals (nw, send, to, winr)  ! put nw floating point values to recv on MPI rank 'to'
!> call mpi_rma%get_reals (nw, recv, fr, wins)  ! get nw floating point values from send on MPI rank 'fr'
!> ...
!> call mpi_coords%end                          ! closes MPI
!>
!> Too see which variables and procedures are inside mpi_rma%, just look below!
!>
!> =============================================================================
MODULE mpi_rma_mod
  USE mpi_mod
  implicit none
  include "mpif.h"
  integer, private:: mpi_err
  type, extends(mpi_t):: mpi_rma_t
  contains
    procedure:: win_ints
    procedure:: win_reals
    procedure:: win_free
    procedure:: get_reals
    procedure:: put_reals
    procedure:: get_ints
    procedure:: put_ints
    procedure:: lock
    procedure:: lock_x
    procedure:: unlock
    procedure:: fence
  end type
  type(mpi_rma_t):: mpi_rma
PRIVATE
PUBLIC mpi_rma_t, mpi_rma
CONTAINS

!> =============================================================================
!> Create an MPI window of size "nw" words (reals)
!> =============================================================================
SUBROUTINE win_reals (mp, win, nw, b)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, flag, mode, mpi_err
  real(4), dimension(nw):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  ! Create MPI remote memory window, of size nw words.
  !----------------------------------------------------------------------------
  nbytes = nw*4
  call MPI_Win_create (b, nbytes, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call mp%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE win_reals

!> =============================================================================
!> Create an MPI window of size "nw" words (ints)
!> =============================================================================
SUBROUTINE win_ints (mp, win, nw, b)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, flag, mode, mpi_err
  integer(4), dimension(nw):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  ! Create MPI remote memory window, of size nw words.
  !----------------------------------------------------------------------------
  nbytes = nw*4
  call MPI_Win_create (b, nbytes, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call mp%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE win_ints

!> =============================================================================
!> Free a window opened with win_ints or win_reals
!> =============================================================================
SUBROUTINE win_free (mp, win)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win
  !----------------------------------------------------------------------------
  ! Assert there is no more window activity, and free the window
  !----------------------------------------------------------------------------
  call MPI_Win_fence (MPI_MODE_NOSUCCEED, win, mpi_err)
  call MPI_Win_free (win, mpi_err)
  call mp%assert ('MPI_Win_free', mpi_err)
END SUBROUTINE win_free

!> =============================================================================
!> Put nb bytes from a window on rank
!> =============================================================================
SUBROUTINE get_reals (mp, win, w, nw, offset, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank, mpi_err
  real:: w(nw)
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_get (w, nw, MPI_CHAR, rank, offset2, nw, MPI_CHAR, win, mpi_err)
  call mp%assert ('MPI_Get', mpi_err)
END SUBROUTINE get_reals

!> =============================================================================
!> Put nb bytes from a window on rank
!> =============================================================================
SUBROUTINE put_reals (mp, win, w, nw, offset, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank, mpi_err
  real:: w(nw)
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_put (w, nw, MPI_CHAR, rank, offset2, nw, MPI_CHAR, win, mpi_err)
  call mp%assert ('MPI_Put', mpi_err)
END SUBROUTINE put_reals

!> =============================================================================
!> Get nb bytes from a window on rank
!> =============================================================================
SUBROUTINE get_ints (mp, win, nb, b, offset, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nb, offset, rank, mpi_err
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_Get (b, nb, MPI_CHAR, rank, offset2, nb, MPI_CHAR, win, mpi_err)
  !call MPI_Win_sync (win, mpi_err)
  call mp%assert ('MPI_Get', mpi_err)
END SUBROUTINE get_ints

!> =============================================================================
!> Put nb bytes from a window on rank
!> =============================================================================
SUBROUTINE put_ints (mp, win, nb, b, offset, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nb, offset, rank, mpi_err
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_put (b, nb, MPI_CHAR, rank, offset2, nb, MPI_CHAR, win, mpi_err)
  call mp%assert ('MPI_Put', mpi_err)
END SUBROUTINE put_ints

!> =============================================================================
!> Lock a window, preparing for access from rank
!> =============================================================================
SUBROUTINE lock (mp, win, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_SHARED, rank, 0, win, mpi_err)
  call mp%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE lock

!> =============================================================================
!> Lock a window, preparing for exclusive access from rank
!> =============================================================================
SUBROUTINE lock_x (mp, win, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, rank, 0, win, mpi_err)
  call mp%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE lock_x

!> =============================================================================
!> Remove a lock
!> =============================================================================
SUBROUTINE unlock (mp, win, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Unlock the access to the window
  !----------------------------------------------------------------------------
  call MPI_Win_unlock (rank, win, mpi_err)
  call mp%assert ('MPI_Win_unlock', mpi_err)
END SUBROUTINE unlock

!> =============================================================================
!> Synchronize the operations on MPI open windows
!> =============================================================================
SUBROUTINE fence (mp, win)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call MPI_Win_fence (0, win, mpi_err)
  call mp%assert ('MPI_Win_fence', mpi_err)
END SUBROUTINE fence

!> =============================================================================
!> Add a value to an integer "hosted" on a certain rank (e.g. 0 = master)
!> =============================================================================
SUBROUTINE add_int (mp, win, rank, add, val)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, rank, n=1, val, add, org, tmp
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  call mp%lock_x (win, rank)
  org = val
  call MPI_Accumulate (add, n, MPI_INTEGER, rank, offset, &
                            n, MPI_INTEGER, MPI_SUM, win, mpi_err)
  call mp%lock (win, rank)
  call mp%assert ('MPI_Accumulate', mpi_err)
END SUBROUTINE add_int

!> =============================================================================
!> Find the smallest value among windows (?)
!> =============================================================================
SUBROUTINE min_real (mp, win, val)
  implicit none
  class(mpi_rma_t):: mp
  real(8):: val
  integer:: win, n=1, rank=0, sz
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !----------------------------------------------------------------------------
  call mp%lock (win, rank)
  call MPI_Accumulate (val, n, MPI_REAL8, rank, offset, &
                            n, MPI_REAL8, MPI_MIN, win, mpi_err)
  call mp%assert ('MPI_Accumulate', mpi_err)
  call mp%unlock (win, rank)
END SUBROUTINE min_real

END MODULE mpi_rma_mod
