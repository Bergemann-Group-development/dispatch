!===============================================================================
!> $Id$
!> MPI remote-memory-access (RMA) calls, also referred to as one-sided MPI
!===============================================================================
MODULE mpi_rma_mod
  USE mpi_lib
  USE mpi_mod
  implicit none
  type, extends(mpi_t):: mpi_rma_t
    integer:: window
    integer:: nr, ni
    real, dimension(:), pointer:: r
    integer, dimension(:), pointer:: i
  contains
    procedure:: create
    procedure:: free
    procedure:: lock
    procedure:: lock_x
    procedure:: unlock
    procedure:: fence
    procedure:: sync
    procedure:: flush
    procedure:: create_dynamic
    procedure:: attach
    procedure:: detach
    procedure:: accumulate
    procedure:: put
    procedure:: get
  end type
PUBLIC
CONTAINS

!===============================================================================
!> Create an MPI window of size "nb" bytes
!===============================================================================
SUBROUTINE create (self, ni, nr)
  class(mpi_rma_t):: self
  integer, optional:: ni, nr
  !.............................................................................
  integer:: nw, flag, mode, mpi_err
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  ! Create MPI remote memory window, of size nb bytes.
  !----------------------------------------------------------------------------
  if (present(nr)) then
    if (nr>0) then
      self%nr = nr
      self%ni = 0
      allocate (self%r(nr))
      nbytes = nr*4
      call MPI_Win_create (self%r, nbytes, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &
                           self%window, mpi_err)
    end if
  else if (present(ni)) then
    if (ni>0) then
      self%ni = ni
      self%nr = 0
      allocate (self%i(ni))
      nbytes = ni*4
      call MPI_Win_create (self%i, nbytes, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &
                           self%window, mpi_err)
    end if
  end if
  call self%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE create

!===============================================================================
!> Free a window opened with create or create_dynamic
!===============================================================================
SUBROUTINE free (self)
  class(mpi_rma_t):: self
  integer:: mpi_err
  !-----------------------------------------------------------------------------
  ! Assert there is no more window activity, and free the window
  !-----------------------------------------------------------------------------
  call MPI_Win_fence (MPI_MODE_NOSUCCEED, self%window, mpi_err)
  call MPI_Win_free (self%window, mpi_err)
  call self%assert ('MPI_Win_free', mpi_err)
END SUBROUTINE free

!===============================================================================
!> Lock a window, preparing for shared access from ranks
!===============================================================================
SUBROUTINE lock (self, rank)
  class(mpi_rma_t):: self
  integer:: rank
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_SHARED, rank, 0, self%window, mpi_err)
  call self%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE lock

!===============================================================================
!> Lock a window, preparing for exclusive access from rank
!===============================================================================
SUBROUTINE lock_x (self, rank)
  class(mpi_rma_t):: self
  integer:: rank
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Lock the access, to allow us to use MPI_Win_get safely
  !----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, rank, 0, self%window, mpi_err)
  call self%assert ('MPI_Win_lock', mpi_err)
END SUBROUTINE lock_x

!===============================================================================
!> Remove a lock
!===============================================================================
SUBROUTINE unlock (self, rank)
  class(mpi_rma_t):: self
  integer:: rank
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Unlock the access to the window
  !----------------------------------------------------------------------------
  call MPI_Win_unlock (rank, self%window, mpi_err)
  call self%assert ('MPI_Win_unlock', mpi_err)
END SUBROUTINE unlock

!===============================================================================
!> Synchronize the operations on MPI open windows
!===============================================================================
SUBROUTINE fence (self)
  class(mpi_rma_t):: self
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call MPI_Win_fence (0, self%window, mpi_err)
  call self%assert ('MPI_Win_fence', mpi_err)
END SUBROUTINE fence

!===============================================================================
!> Force a sync
!===============================================================================
SUBROUTINE sync (self)
  class(mpi_rma_t):: self
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Sync windows
  !----------------------------------------------------------------------------
  call MPI_Win_sync (self%window, mpi_err)
  call self%assert ('MPI_Win_sync', mpi_err)
END SUBROUTINE sync

!===============================================================================
!> Flush
!===============================================================================
SUBROUTINE flush (self)
  class(mpi_rma_t):: self
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Flush operations
  !----------------------------------------------------------------------------
  call MPI_Win_flush (self%window, mpi_err)
  call self%assert ('MPI_Win_flush', mpi_err)
END SUBROUTINE flush

!===============================================================================
!> Create a window that can be changed dynamically
!===============================================================================
SUBROUTINE create_dynamic (self)
  class(mpi_rma_t):: self
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Create dynamic MPI remote memory window.
  !----------------------------------------------------------------------------
#ifndef __GFORTRAN__
  call MPI_Win_create_dynamic (MPI_INFO_NULL, MPI_COMM_WORLD, self%window, mpi_err)
  call self%assert ('MPI_Win_create', mpi_err)
#endif
END SUBROUTINE create_dynamic

!===============================================================================
!> Add a value to an integer "hosted" on a certain rank (e.g. 0 = master)
!===============================================================================
SUBROUTINE accumulate (self, target, r, i)
  class(mpi_rma_t):: self
  integer:: target
  integer, optional:: i(:)
  real,    optional:: r(:)
  !.............................................................................
  integer:: n, mpi_err
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !-----------------------------------------------------------------------------
  ! Note that offset is counted in bytes, while numbers are storage units
  !-----------------------------------------------------------------------------
  if (self%nr > 0) then
    n = size(r)
    call MPI_Get_accumulate (r, n, MPI_REAL, &
                        self%r, n, MPI_REAL, target, offset, &
                                n, MPI_REAL, MPI_SUM, self%window, mpi_err)
  else if (self%ni > 0) then
    n = size(i)
    call MPI_Get_accumulate (i, n, MPI_INT, &
                        self%i, n, MPI_INT, target, offset, &
                                n, MPI_INT, MPI_SUM, self%window, mpi_err)
  end if
  call self%assert ('MPI_Accumulate', mpi_err)
  call MPI_Win_flush (target, self%window, mpi_err)
  call self%assert ('MPI_Flush', mpi_err)
END SUBROUTINE accumulate

!===============================================================================
!> Define a new range for a window
!===============================================================================
SUBROUTINE attach (self, b, nb)
  class(mpi_rma_t):: self
  integer:: nb
  character(len=1), dimension(nb):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  integer:: mpi_err
  !----------------------------------------------------------------------------
  ! Attach memory to MPI dynamic window, of size nb bytes.
  !----------------------------------------------------------------------------
  nbytes = nb
#ifndef __GFORTRAN__
  call MPI_Win_attach (self%window, b, nbytes, mpi_err)
#endif
  call self%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE attach

!===============================================================================
!> Cancel a previous attach
!===============================================================================
SUBROUTINE detach (self, b)
  class(mpi_rma_t):: self
  integer:: mpi_err
  character(len=1), dimension(1):: b
  integer(kind=MPI_ADDRESS_KIND):: nbytes
  !----------------------------------------------------------------------------
  ! Detach memory from MPI window
  !----------------------------------------------------------------------------
#ifndef __GFORTRAN__
  call MPI_Win_detach (self%window, b, mpi_err)
#endif
  call self%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE detach

!===============================================================================
!> Get window values from target
!===============================================================================
SUBROUTINE put (self, target)
  class(mpi_rma_t):: self
  integer:: target, mpi_err
  integer(kind=MPI_ADDRESS_KIND):: offset
  !----------------------------------------------------------------------------
  offset = 0
  if (self%ni > 0) then
    call MPI_Put (self%i, self%ni, MPI_INT,  target, offset, &
                          self%ni, MPI_INT,  self%window, mpi_err)
    call self%assert ('MPI_Put', mpi_err)
  else if (self%nr > 0) then
    call MPI_Put (self%r, self%nr, MPI_REAL, target, offset, &
                          self%nr, MPI_REAL, self%window, mpi_err)
    call self%assert ('MPI_Put', mpi_err)
  end if
END SUBROUTINE put

!===============================================================================
!> Get window values from target
!===============================================================================
SUBROUTINE get (self, target)
  class(mpi_rma_t):: self
  integer:: target, mpi_err
  integer(kind=MPI_ADDRESS_KIND):: offset
  !----------------------------------------------------------------------------
  offset = 0
  if (self%ni > 0) then
    call MPI_Get (self%i, self%ni, MPI_INT,  target, offset, &
                          self%ni, MPI_INT,  self%window, mpi_err)
    call self%assert ('MPI_Get', mpi_err)
  else if (self%nr > 0) then
    call MPI_Get (self%r, self%nr, MPI_REAL, target, offset, &
                          self%nr, MPI_REAL, self%window, mpi_err)
      call self%assert ('MPI_Get', mpi_err)
  end if
END SUBROUTINE get

END MODULE mpi_rma_mod
