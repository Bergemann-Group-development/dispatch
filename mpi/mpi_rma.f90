!> =============================================================================
!> $Id: 79486ae72058c4609ebd8e5f4be02792e50c8ca1 $
!> MPI remote-memory-access (RMA) calls, also referred to as one-sided MPI
!> =============================================================================
MODULE mpi_rma
  USE mpi_basic
  USE io
  implicit none
PUBLIC
CONTAINS

!> =============================================================================
!> Free a window opened with win_create_mpi or win_create_dynamic_mpi
!> =============================================================================
SUBROUTINE win_free_mpi (win)
  implicit none
  integer:: win
  !----------------------------------------------------------------------------
  ! Assert there is no more window activity, and free the window
  !----------------------------------------------------------------------------
  call trace_begin ('win_free_mpi')
  call MPI_Win_fence (MPI_MODE_NOSUCCEED, win, mpi_err)
  call MPI_Win_free (win, mpi_err)
  call assert_mpi ('MPI_Win_free', mpi_err)
  call trace_end
END SUBROUTINE win_free_mpi

!> =============================================================================
!> Lock a window, preparing for access from rank
!> =============================================================================
SUBROUTINE win_lock_mpi (win, rank)
  implicit none
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_SHARED, rank, 0, win, mpi_err)
  call assert_mpi ('MPI_Win_lock', mpi_err)
END SUBROUTINE win_lock_mpi

!> =============================================================================
!> Lock a window, preparing for exclusive access from rank
!> =============================================================================
SUBROUTINE win_lock_x_mpi (win, rank)
  implicit none
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, rank, 0, win, mpi_err)
  call assert_mpi ('MPI_Win_lock', mpi_err)
END SUBROUTINE win_lock_x_mpi

!> =============================================================================
!> Remove a lock
!> =============================================================================
SUBROUTINE win_unlock_mpi (win, rank)
  implicit none
  integer:: win, nw, offset, rank
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Unlock the access to the window
  !----------------------------------------------------------------------------
  if (do_trace.and.verbose>1) write(*,*) 'unlock waiting for rank', rank
  if (do_trace.and.verbose>0) write(output,*) 'unlock waiting for rank', rank
  call MPI_Win_unlock (rank, win, mpi_err)
  call assert_mpi ('MPI_Win_unlock', mpi_err)
END SUBROUTINE win_unlock_mpi

!> =============================================================================
!> Synchronize the operations on MPI open windows
!> =============================================================================
SUBROUTINE win_fence_mpi (win)
  implicit none
  integer:: win
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call trace_begin ('win_fence_mpi')
  call MPI_Win_fence (0, win, mpi_err)
  call assert_mpi ('MPI_Win_fence', mpi_err)
  call trace_end
END SUBROUTINE win_fence_mpi

!> =============================================================================
!> Force a sync
!> =============================================================================
SUBROUTINE win_sync_mpi (win, rank)
  implicit none
  integer:: win, rank
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call trace_begin ('win_sync_mpi')
  call MPI_Win_sync (rank, win, mpi_err)
  call assert_mpi ('MPI_Win_sync', mpi_err)
  call trace_end
END SUBROUTINE win_sync_mpi

!> =============================================================================
!> Create a window that can be changed dynamically
!> =============================================================================
SUBROUTINE win_create_dynamic_mpi (win)
  implicit none
  integer:: win
  !----------------------------------------------------------------------------
  ! Create dynamic MPI remote memory window.
  !----------------------------------------------------------------------------
  !call trace_begin ('win_create_mpi')
#ifndef __GFORTRAN__
  call MPI_Win_create_dynamic (MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call assert_mpi ('MPI_Win_create', mpi_err)
#endif
  !call trace_end
END SUBROUTINE win_create_dynamic_mpi

!> =============================================================================
!> Add a value to an integer "hosted" on a certain rank (e.g. 0 = master)
!> =============================================================================
SUBROUTINE win_add_int_mpi (win, rank, add, val)
  implicit none
  integer:: win, rank, n=1, val, add, org, tmp
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  call trace_begin ('win_add_mpi')
  call win_lock_x_mpi (win, rank)
  org = val
  call MPI_Accumulate (add, n, MPI_INTEGER, rank, offset, &
                            n, MPI_INTEGER, MPI_SUM, win, mpi_err)
  call assert_mpi ('MPI_Accumulate', mpi_err)
  !call win_get_mpi (win, tmp, sizeof(tmp), offset, rank)
  !call MPI_Win_flush (rank, win, mpi_err)
  !tmp = tmp+add
  !call win_put_mpi (win, tmp, sizeof(tmp), offset, rank)
  call MPI_Win_flush (rank, win, mpi_err)
  call win_get_mpi (win, val, sizeof(val), offset, rank)
  call win_unlock_mpi (win, rank)
  !print'(a,4i5)', 'org, add, val =', org, add, val
  call trace_end
END SUBROUTINE win_add_int_mpi

!> =============================================================================
!> Find the smallest value among windows (?)
!> =============================================================================
SUBROUTINE win_min_real_mpi (win, val)
  implicit none
  real(8):: val
  integer:: win, n=1, rank=0, sz
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !----------------------------------------------------------------------------
  call trace_begin ('win_min_real_mpi')
  call win_lock_mpi (win, rank)
  call MPI_Accumulate (val, n, MPI_REAL8, rank, offset, &
                            n, MPI_REAL8, MPI_MIN, win, mpi_err)
  call win_unlock_mpi (win, rank)
  sz = sizeof(val)
  call win_lock_mpi (win, rank)
  call win_get_mpi (win, val, sz, offset, rank)
  call win_unlock_mpi (win, rank)
  call assert_mpi ('MPI_Accumulate', mpi_err)
  call trace_end
END SUBROUTINE win_min_real_mpi

END MODULE mpi_rma

!> =============================================================================
!> Create an MPI window of size "nb" bytes
!> =============================================================================
SUBROUTINE win_create_mpi (win, b, nb)
  USE mpi_rma
  implicit none
  integer:: win, nb, flag, mode
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  ! Create MPI remote memory window, of size nb bytes.
  !----------------------------------------------------------------------------
  !call trace_begin ('win_create_mpi')
  nbytes = nb
  call MPI_Win_create (b, nbytes, 1, MPI_INFO_NULL, MPI_COMM_WORLD, win, mpi_err)
  call assert_mpi ('MPI_Win_create', mpi_err)
#ifndef __GFORTRAN__
  if (master .and. first_time) then
    first_time = .false.
    !call MPI_Win_get_attr (win, MPI_WIN_MODEL, mode, flag, mpi_err)
    !call assert_mpi ('MPI_Win_get_attr', mpi_err)
    !if (mode==MPI_WIN_UNIFIED)  print *,'MPI_WIN_MODEL is MPI_WIN_UNIFIED'
    !if (mode==MPI_WIN_SEPARATE) print *,'MPI_WIN_MODEL is MPI_WIN_SEPARATE'
  end if
#endif
  !call trace_end
END SUBROUTINE win_create_mpi

!> =============================================================================
!> Define a new range for a window
!> =============================================================================
SUBROUTINE win_attach_mpi (win, b, nb)
  USE mpi_rma
  implicit none
  integer:: win, nb
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  !----------------------------------------------------------------------------
  ! Attach memory to MPI dynamic window, of size nb bytes.
  !----------------------------------------------------------------------------
  call trace_begin ('win_create_mpi')
  nbytes = nb
#ifndef __GFORTRAN__
  call MPI_Win_attach (win, b, nbytes, mpi_err)
#endif
  call assert_mpi ('MPI_Win_create', mpi_err)
  call trace_end
END SUBROUTINE win_attach_mpi

!> =============================================================================
!> Cancel a previous win_attach
!> =============================================================================
SUBROUTINE win_detach_mpi (win, b)
  USE mpi_rma
  implicit none
  integer:: win
  character(len=1), dimension(1):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  !----------------------------------------------------------------------------
  ! Detach memory from MPI window
  !----------------------------------------------------------------------------
  call trace_begin ('win_create_mpi')
#ifndef __GFORTRAN__
  call MPI_Win_detach (win, b, mpi_err)
#endif
  call assert_mpi ('MPI_Win_create', mpi_err)
  call trace_end
END SUBROUTINE win_detach_mpi

!> =============================================================================
!> Get nb bytes from a window on rank
!> =============================================================================
SUBROUTINE win_get_mpi (win, b, nb, offset, rank)
  USE mpi_rma
  implicit none
  integer:: win, nb, offset, rank
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  if (verbose>1) call trace_begin ('win_get_mpi')
  offset2 = offset
  call MPI_Get (b, nb, MPI_CHAR, rank, offset2, nb, MPI_CHAR, win, mpi_err)
  !call MPI_Win_sync (win, mpi_err)
  call assert_mpi ('MPI_Get', mpi_err)
  if (verbose>1) call trace_end
END SUBROUTINE win_get_mpi

!> =============================================================================
!> Put nb bytes from a window on rank
!> =============================================================================
SUBROUTINE win_put_mpi (win, b, nb, offset, rank)
  USE mpi_rma
  implicit none
  integer:: win, nb, offset, rank
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: offset2
  !----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !----------------------------------------------------------------------------
  call trace_begin ('win_put_mpi')
  offset2 = offset
  call MPI_put (b, nb, MPI_CHAR, rank, offset2, nb, MPI_CHAR, win, mpi_err)
  call assert_mpi ('MPI_Get', mpi_err)
  call trace_end
END SUBROUTINE win_put_mpi
