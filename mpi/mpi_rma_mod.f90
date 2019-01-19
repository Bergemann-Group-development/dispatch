!> =============================================================================
!> $Id$
!> MPI remote-memory-access (RMA) calls, also referred to as one-sided MPI
!> =============================================================================
MODULE mpi_rma_mod
  USE mpi_mod
  implicit none
  include "mpif.h"
  integer, private:: mpi_err
  interface win_get
    module procedure get_int, get_real
  end interface
  type, extends(mpi_t):: mpi_rma_t
  contains
    procedure:: create_int
    procedure:: create_real
    procedure:: free
    procedure:: lock
    procedure:: lock_x
    procedure:: unlock
    procedure:: fence
    procedure:: sync
    procedure:: create_dynamic
    procedure:: attach
    procedure:: detach
    procedure:: add_int
    procedure:: put_int
    procedure:: put_real
    procedure:: get_int
    procedure:: get_real
  end type
PUBLIC
CONTAINS

!> =============================================================================
!> Create an MPI window of size "nb" bytes
!> =============================================================================
SUBROUTINE create_int (self, win, b, nw)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nw, flag, mode, mpi_err
  integer, dimension(nw):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  ! Create MPI remote memory window, of size nb bytes.
  !----------------------------------------------------------------------------
  nbytes = nw*4
  call MPI_Win_create (b, nbytes, 1, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call self%assert ('MPI_Win_create', mpi_err)
#ifndef __GFORTRAN__
  if (master .and. first_time) then
    first_time = .false.
    !call MPI_Win_get_attr (win, MPI_WIN_MODEL, mode, flag, mpi_err)
    !call self%assert ('MPI_Win_get_attr', mpi_err)
    !if (mode==MPI_WIN_UNIFIED)  print *,'MPI_WIN_MODEL is MPI_WIN_UNIFIED'
    !if (mode==MPI_WIN_SEPARATE) print *,'MPI_WIN_MODEL is MPI_WIN_SEPARATE'
  end if
#endif
END SUBROUTINE create_int

!> =============================================================================
!> Create an MPI window of size "nb" bytes
!> =============================================================================
SUBROUTINE create_real (self, win, b, nw)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nw, flag, mode, mpi_err
  real(4), dimension(nw):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  ! Create MPI remote memory window, of size nw words.
  !----------------------------------------------------------------------------
  nbytes = nw
  call MPI_Win_create (b, nbytes, 4, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call self%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE create_real

!> =============================================================================
!> Free a window opened with create or create_dynamic
!> =============================================================================
SUBROUTINE free (self, win)
  implicit none
  class(mpi_rma_t):: self
  integer:: win
  !----------------------------------------------------------------------------
  ! Assert there is no more window activity, and free the window
  !----------------------------------------------------------------------------
  call MPI_Win_fence (MPI_MODE_NOSUCCEED, win, mpi_err)
  call MPI_Win_free (win, mpi_err)
  call self%assert ('MPI_Win_free', mpi_err)
END SUBROUTINE free

!> =============================================================================
!> Lock a window, preparing for access from rank
!> =============================================================================
SUBROUTINE lock (self, win, rank)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_SHARED, rank, 0, win, mpi_err)
  call self%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE lock

!> =============================================================================
!> Lock a window, preparing for exclusive access from rank
!> =============================================================================
SUBROUTINE lock_x (self, win, rank)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, rank, 0, win, mpi_err)
  call self%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE lock_x

!> =============================================================================
!> Remove a lock
!> =============================================================================
SUBROUTINE unlock (self, win, rank)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Unlock the access to the window
  !----------------------------------------------------------------------------
  call MPI_Win_unlock (rank, win, mpi_err)
  call self%assert ('MPI_Win_unlock', mpi_err)
END SUBROUTINE unlock

!> =============================================================================
!> Synchronize the operations on MPI open windows
!> =============================================================================
SUBROUTINE fence (self, win)
  implicit none
  class(mpi_rma_t):: self
  integer:: win
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call MPI_Win_fence (0, win, mpi_err)
  call self%assert ('MPI_Win_fence', mpi_err)
END SUBROUTINE fence

!> =============================================================================
!> Force a sync
!> =============================================================================
SUBROUTINE sync (self, win)
  implicit none
  class(mpi_rma_t):: self
  integer:: win
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call MPI_Win_sync (win, mpi_err)
  call self%assert ('MPI_Win_sync', mpi_err)
END SUBROUTINE sync

!> =============================================================================
!> Create a window that can be changed dynamically
!> =============================================================================
SUBROUTINE create_dynamic (self, win)
  implicit none
  class(mpi_rma_t):: self
  integer:: win
  !----------------------------------------------------------------------------
  ! Create dynamic MPI remote memory window.
  !----------------------------------------------------------------------------
#ifndef __GFORTRAN__
  call MPI_Win_create_dynamic (MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call self%assert ('MPI_Win_create', mpi_err)
#endif
END SUBROUTINE create_dynamic

!> =============================================================================
!> Add a value to an integer "hosted" on a certain rank (e.g. 0 = master)
!> =============================================================================
SUBROUTINE add_int (self, win, rank, add, val)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, rank, n=1, val(:), add(:)
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  call MPI_Get_accumulate (add, n, MPI_INTEGER, &
                           val, n, MPI_INTEGER, rank, offset, &
                                n, MPI_INTEGER, MPI_SUM, win, mpi_err)
  call self%assert ('MPI_Accumulate', mpi_err)
  call MPI_Win_flush (rank, win, mpi_err)
  call self%assert ('MPI_Flush', mpi_err)
END SUBROUTINE add_int

!> =============================================================================
!> Define a new range for a window
!> =============================================================================
SUBROUTINE attach (self, win, b, nb)
  implicit none
  class(mpi_rma_t):: self
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
  call self%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE attach

!> =============================================================================
!> Cancel a previous attach
!> =============================================================================
SUBROUTINE detach (self, win, b)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, mpi_err
  character(len=1), dimension(1):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  !----------------------------------------------------------------------------
  ! Detach memory from MPI window
  !----------------------------------------------------------------------------
#ifndef __GFORTRAN__
  call MPI_Win_detach (win, b, mpi_err)
#endif
  call self%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE detach

!> =============================================================================
!> Get ni ints from a window on rank
!> =============================================================================
SUBROUTINE get_int (self, win, rank, i, ni, offset)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, ni, offset, rank, mpi_err
  integer, dimension(ni):: i
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_Get (i, ni, MPI_INT, rank, offset2, ni, MPI_INT, win, mpi_err)
  !call MPI_Win_sync (win, mpi_err)
  call self%assert ('MPI_Get (int)', mpi_err)
END SUBROUTINE get_int

!> =============================================================================
!> Get ni ints from a window on rank
!> =============================================================================
SUBROUTINE get_real (self, win, rank, r, nr, offset)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nr, offset, rank, mpi_err
  real, dimension(nr):: r
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_Get (r, nr, MPI_REAL, rank, offset2, nr, MPI_REAL, win, mpi_err)
  call MPI_Win_sync (win, mpi_err)
  call self%assert ('MPI_Get (real)', mpi_err)
END SUBROUTINE get_real

!> =============================================================================
!> Put nb bytes from a window on rank
!> =============================================================================
SUBROUTINE put_int (self, win, rank, i, nw, offset)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nw, offset, rank, mpi_err
  integer, dimension(nw):: i
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_put (i, nw, MPI_INT, rank, offset2, nw, MPI_INT, win, mpi_err)
  call self%assert ('MPI_Get', mpi_err)
END SUBROUTINE put_int

!> =============================================================================
!> Put nb bytes from a window on rank
!> =============================================================================
SUBROUTINE put_real (self, win, rank, r, nw, offset)
  implicit none
  class(mpi_rma_t):: self
  integer:: win, nw, offset, rank, mpi_err
  real(4), dimension(nw):: r
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  offset2 = offset
  call MPI_put (r, nw, MPI_REAL, rank, offset2, nw, MPI_REAL, win, mpi_err)
  call self%assert ('MPI_Get', mpi_err)
END SUBROUTINE put_real

END MODULE mpi_rma_mod
