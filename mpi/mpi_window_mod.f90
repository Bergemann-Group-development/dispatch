MODULE mpi_window_mod
  USE trace_mod
  USE mpi_lib
  USE mpi_mod, mp => mpi
  implicit none
  private
  type, public:: mpi_window_t
    real:: value
    integer:: window
  contains
    procedure:: init
    procedure:: get
    procedure:: put
    procedure:: test
  end type
#ifdef MPI
  integer(kind=MPI_ADDRESS_KIND):: offset=0, nbytes=4
#endif
CONTAINS

!===============================================================================
! Creat a single real word window
!===============================================================================
SUBROUTINE init (self, value)
  class(mpi_window_t):: self
  real:: value
  integer:: err
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_window_t%init')
  self%value = value
#ifdef MPI
  !$omp barrier
  !$omp single
  call MPI_Win_create (self%value, nbytes, 1, MPI_INFO_NULL, &
                       MPI_COMM_WORLD, self%window, err)
  !$omp end single
#endif
  call trace%end ()
END SUBROUTINE init

!===============================================================================
! Get the current value from another rank
!===============================================================================
FUNCTION get (self, rank) RESULT (value)
  class(mpi_window_t):: self
  integer:: rank
  real:: value
  integer:: err
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_window_t%get')
#ifdef MPI
  !$omp critical (mpi_get_cr)
  call MPI_Win_lock (MPI_LOCK_SHARED, rank, 0, self%window, err)
  call MPI_Get (value, 1, MPI_REAL, rank, offset, 1, MPI_REAL, self%window, err)
  call MPI_Win_unlock (rank, self%window, err)
  !$omp end critical (mpi_get_cr)
#endif
  call trace%end ()
END FUNCTION get

!===============================================================================
! Get the current value from another rank
!===============================================================================
SUBROUTINE put (self, rank, value)
  class(mpi_window_t):: self
  integer:: rank
  real:: value
  integer:: err
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_window_t%put')
#ifdef MPI
  !$omp critical (mpi_put_cr)
  call MPI_Win_lock (MPI_LOCK_EXCLUSIVE, rank, 0, self%window, err)
  call MPI_Put (value, 1, MPI_REAL, rank, offset, 1, MPI_REAL, self%window, err)
  call MPI_Win_unlock (rank, self%window, err)
  !$omp end critical (mpi_put_cr)
#endif
  call trace%end ()
END SUBROUTINE put

!===============================================================================
! Creat a single real word window
!===============================================================================
SUBROUTINE test (self)
  class(mpi_window_t):: self
  integer:: rank, err
  !-----------------------------------------------------------------------------
  call self%init (1.1+mp%rank)
  do rank=0,mp%size-1
    if (rank /= mp%rank) then
      print *, mp%rank, rank, self%get(rank)
    end if
  end do
END SUBROUTINE test

END MODULE mpi_window_mod
