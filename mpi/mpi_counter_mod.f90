!> =============================================================================
!> Simple MPI counter, which allows setting, incrementing,  and decrementing a
!> counter atomically
!> =============================================================================
MODULE mpi_counter_mod
  USE mpi_mod
  implicit none
  private
  include "mpif.h"
  type, public:: mpi_counter_t
    integer:: window
    integer(8):: i=0
    integer(8):: n=0
  contains
    procedure, private:: init4
    procedure, private:: init8
    generic:: init => init4, init8
    procedure, private:: update4
    procedure, private:: update8
    generic:: update => update4, update8
  end type
CONTAINS

!> =============================================================================
!> Create and initialize a counter to value = n. NOTE: This is COLLECTIVE call!
!> =============================================================================
SUBROUTINE init4 (self, n)
  class(mpi_counter_t):: self
  integer:: n
  !.............................................................................
  integer(kind=MPI_ADDRESS_KIND):: nbytes=8
  integer:: mpi_err
  !----------------------------------------------------------------------------
  self%i = n
  self%n = n
  call MPI_Win_create (self%i, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                       self%window, mpi_err)
  call mpi%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE init4

!> =============================================================================
!> Create and initialize a counter to value = n. NOTE: This is COLLECTIVE call!
!> =============================================================================
SUBROUTINE init8 (self, n)
  class(mpi_counter_t):: self
  integer(8):: n
  !.............................................................................
  integer(kind=MPI_ADDRESS_KIND):: nbytes=8
  integer:: mpi_err
  !----------------------------------------------------------------------------
  self%i = n
  self%n = n
  call MPI_Win_create (self%i, nbytes, 8, MPI_INFO_NULL, MPI_COMM_WORLD, &
                       self%window, mpi_err)
  call mpi%assert ('MPI_Win_create', mpi_err)
END SUBROUTINE init8

!> =============================================================================
!> Add a value i to the counter on master.  If the value reaches 0, return 0,
!> while resetting the counter to n.  Allow both 4-byte and 8-byte arguments.
!> =============================================================================
FUNCTION update4 (self, i, n) RESULT (j)
  class(mpi_counter_t):: self
  integer:: i, j
  integer, optional:: n
  !.............................................................................
  if (present(n)) then
    j= update8 (self, int(i,kind=8), int(n,kind=8))
  else
    j= update8 (self, int(i,kind=8))
  end if
END FUNCTION update4
!> =============================================================================
FUNCTION update8 (self, i, n) RESULT (j)
  class(mpi_counter_t):: self
  integer(8):: i, j
  integer(8), optional:: n
  !.............................................................................
  integer:: master=0, mpi_err
  integer(kind=MPI_ADDRESS_KIND):: offset=0
  !-----------------------------------------------------------------------------
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, master, 0, self%window, mpi_err)
  call add(i)
  call MPI_Win_flush (master, self%window, mpi_err)
  if (mpi%rank == master) then
    j = self%i - i
  else
    j = self%i
  end if
  !-----------------------------------------------------------------------------
  ! If / when the (actual) counter hits zero, add n to it, if present
  !-----------------------------------------------------------------------------
  if (j+i == 0) then
    if (present(n)) then
      call add (n)
    else
      call add (self%n)
    end if
  end if
  call MPI_Win_unlock (master, self%window, mpi_err)
contains
  !-----------------------------------------------------------------------------
  subroutine add (i)
  integer(8):: i
  call MPI_Get_accumulate (i, 1, MPI_INTEGER8, &
                      self%i, 1, MPI_INTEGER8, master, offset, &
                              1, MPI_INTEGER8, MPI_SUM, self%window, mpi_err)
  call mpi%assert ('MPI_Accumulate', mpi_err)
  end subroutine
END FUNCTION update8

END MODULE mpi_counter_mod
