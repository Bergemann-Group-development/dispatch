!> =============================================================================
!> $Id$
!> MPI remote-memory-access (RMA) calls, also referred to as one-sided MPI
!> =============================================================================
MODULE mpi_rma
  USE mpi_basic
  implicit none
  include "mpif.h"
  integer, private:: mpi_err
  type, extends(mpi_t):: mpi_rma_t
  contains
    procedure:: win_create
    procedure:: win_create_real
    procedure:: win_free
    procedure:: win_lock
    procedure:: win_lock_x
    procedure:: win_unlock
    procedure:: win_fence
    procedure:: win_sync
    procedure:: win_create_dynamic
    procedure:: win_attach
    procedure:: win_detach
    procedure:: win_get
    procedure:: win_put
  end type
PUBLIC
CONTAINS

!> =============================================================================
!> Free a window opened with win_create or win_create_dynamic
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
!> Lock a window, preparing for access from rank
!> =============================================================================
SUBROUTINE win_lock (mp, win, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_SHARED, rank, 0, win, mpi_err)
  call mp%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE win_lock

!> =============================================================================
!> Lock a window, preparing for exclusive access from rank
!> =============================================================================
SUBROUTINE win_lock_x (mp, win, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, rank, 0, win, mpi_err)
  call mp%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE win_lock_x

!> =============================================================================
!> Remove a lock
!> =============================================================================
SUBROUTINE win_unlock (mp, win, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Unlock the access to the window
  !----------------------------------------------------------------------------
  call MPI_Win_unlock (rank, win, mpi_err)
  call mp%assert ('MPI_Win_unlock', mpi_err)
END SUBROUTINE win_unlock

!> =============================================================================
!> Synchronize the operations on MPI open windows
!> =============================================================================
SUBROUTINE win_fence (mp, win)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call MPI_Win_fence (0, win, mpi_err)
  call mp%assert ('MPI_Win_fence', mpi_err)
END SUBROUTINE win_fence

!> =============================================================================
!> Force a sync
!> =============================================================================
SUBROUTINE win_sync (mp, win, rank)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, rank
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call MPI_Win_sync (rank, win, mpi_err)
  call mp%assert ('MPI_Win_sync', mpi_err)
END SUBROUTINE win_sync

!> =============================================================================
!> Create a window that can be changed dynamically
!> =============================================================================
SUBROUTINE win_create_dynamic (mp, win)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win
  !----------------------------------------------------------------------------
  ! Create dynamic MPI remote memory window.
  !----------------------------------------------------------------------------
#ifndef __GFORTRAN__
  call MPI_Win_create_dynamic (MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call mp%assert ('MPI_Win_create', mpi_err)
#endif
END SUBROUTINE win_create_dynamic

!> =============================================================================
!> Add a value to an integer "hosted" on a certain rank (e.g. 0 = master)
!> =============================================================================
SUBROUTINE win_add_int (mp, win, rank, add, val)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, rank, n=1, val, add, org, tmp
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  call mp%win_lock_x (win, rank)
  org = val
  call MPI_Accumulate (add, n, MPI_INTEGER, rank, offset, &
                            n, MPI_INTEGER, MPI_SUM, win, mpi_err)
  call mp%assert ('MPI_Accumulate', mpi_err)
  !call mp%win_get (win, tmp, sizeof(tmp), offset, rank)
  !call MPI_Win_flush (rank, win, mpi_err)
  !tmp = tmp+add
  !call mp%win_put (win, tmp, sizeof(tmp), offset, rank)
  call MPI_Win_flush (rank, win, mpi_err)
  !!!! FIXME call mp%win_get (win, val, sizeof(val), offset, rank)
  call mp%win_unlock (win, rank)
  !print'(a,4i5)', 'org, add, val =', org, add, val
END SUBROUTINE win_add_int

!> =============================================================================
!> Find the smallest value among windows (?)
!> =============================================================================
SUBROUTINE win_min_real (mp, win, val)
  implicit none
  class(mpi_rma_t):: mp
  real(8):: val
  integer:: win, n=1, rank=0, sz
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !----------------------------------------------------------------------------
  call mp%win_lock (win, rank)
  call MPI_Accumulate (val, n, MPI_REAL8, rank, offset, &
                            n, MPI_REAL8, MPI_MIN, win, mpi_err)
  call mp%win_unlock (win, rank)
  sz = sizeof(val)
  call mp%win_lock (win, rank)
  !!!! FIXME call mp%win_get (win, val, sz, offset, rank)
  call mp%win_unlock (win, rank)
  call mp%assert ('MPI_Accumulate', mpi_err)
END SUBROUTINE win_min_real

!> =============================================================================
!> Create an MPI window of size "nb" bytes
!> =============================================================================
SUBROUTINE win_create (mp, win, b, nb)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nb, flag, mode, mpi_err
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  ! Create MPI remote memory window, of size nb bytes.
  !----------------------------------------------------------------------------
  nbytes = nb
  call MPI_Win_create (b, nbytes, 1, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call mp%assert ('MPI_Win_create', mpi_err)
#ifndef __GFORTRAN__
  if (master .and. first_time) then
    first_time = .false.
    !call MPI_Win_get_attr (win, MPI_WIN_MODEL, mode, flag, mpi_err)
    !call mp%assert ('MPI_Win_get_attr', mpi_err)
    !if (mode==MPI_WIN_UNIFIED)  print *,'MPI_WIN_MODEL is MPI_WIN_UNIFIED'
    !if (mode==MPI_WIN_SEPARATE) print *,'MPI_WIN_MODEL is MPI_WIN_SEPARATE'
  end if
#endif
END SUBROUTINE win_create

!> =============================================================================
!> Create an MPI window of size "nb" bytes
!> =============================================================================
SUBROUTINE win_create_real (mp, win, b, nw)
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
#ifndef __GFORTRAN__
  if (master .and. first_time) then
    first_time = .false.
    !call MPI_Win_get_attr (win, MPI_WIN_MODEL, mode, flag, mpi_err)
    !call mp%assert ('MPI_Win_get_attr', mpi_err)
    !if (mode==MPI_WIN_UNIFIED)  print *,'MPI_WIN_MODEL is MPI_WIN_UNIFIED'
    !if (mode==MPI_WIN_SEPARATE) print *,'MPI_WIN_MODEL is MPI_WIN_SEPARATE'
  end if
#endif
END SUBROUTINE win_create_real

!> =============================================================================
!> Define a new range for a window
!> =============================================================================
SUBROUTINE win_attach (mp, win, b, nb)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, nb, mpi_err
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  !----------------------------------------------------------------------------
  ! Attach memory to MPI dynamic window, of size nb bytes.
  !----------------------------------------------------------------------------
  nbytes = nb
#ifndef __GFORTRAN__
  call MPI_Win_attach (win, b, nbytes, mpi_err)
#endif
  call mp%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE win_attach

!> =============================================================================
!> Cancel a previous win_attach
!> =============================================================================
SUBROUTINE win_detach (mp, win, b)
  implicit none
  class(mpi_rma_t):: mp
  integer:: win, mpi_err
  character(len=1), dimension(1):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  !----------------------------------------------------------------------------
  ! Detach memory from MPI window
  !----------------------------------------------------------------------------
#ifndef __GFORTRAN__
  call MPI_Win_detach (win, b, mpi_err)
#endif
  call mp%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE win_detach

!> =============================================================================
!> Get nb bytes from a window on rank
!> =============================================================================
SUBROUTINE win_get (mp, win, b, nb, offset, rank)
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
END SUBROUTINE win_get

!> =============================================================================
!> Put nb bytes from a window on rank
!> =============================================================================
SUBROUTINE win_put (mp, win, b, nb, offset, rank)
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
  call mp%assert ('MPI_Get', mpi_err)
END SUBROUTINE win_put

END MODULE mpi_rma
