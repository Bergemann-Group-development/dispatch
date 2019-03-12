!===============================================================================
!> Simple module for sending / receiving packages to / from other ranks
!===============================================================================
MODULE mpi_comm_mod
#ifdef MPI
  USE mpi
#endif MPI
  USE mpi_mod, mp => mpi
  implicit none
  private
  type, public:: mpi_comm_t
    integer:: req=0
    integer:: ierr
    logical:: flag
    integer, dimension(:), allocatable:: buffer
  contains
    procedure:: send
    procedure:: sent
    procedure:: recv
  end type
CONTAINS

!===============================================================================
!> Send a buffer to a rank, with a given tag, and ignore the request
!===============================================================================
SUBROUTINE send (self, rank, tag)
  class(mpi_comm_t):: self
  integer:: rank, tag
  !-----------------------------------------------------------------------------
#ifdef MPI
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_Isend (self%buffer, size(self%buffer), MPI_INTEGER, rank, tag, &
                    mp%comm, self%req, self%ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_Isend (self%buffer, size(self%buffer), MPI_INTEGER, rank, tag, &
                    mp%comm, self%req, self%ierr)
    !$omp end critical (mpi_cr)
  end if
#endif MPI
END SUBROUTINE send

!===============================================================================
!> If not send request is active, start a send.  When finished, clear request.
!===============================================================================
FUNCTION sent (self)
  logical:: sent
  class(mpi_comm_t):: self
#ifdef MPI
  !............................................................................
  integer:: stat(MPI_STATUS_SIZE)
  !-----------------------------------------------------------------------------
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    call MPI_TEST (self%req, sent, stat, self%ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_TEST (self%req, sent, stat, self%ierr)
    !$omp end critical (mpi_cr)
  end if
#else
  sent = .false.
#endif
  if (sent) self%req = 0
END FUNCTION sent

!===============================================================================
!> If no request is active, issue an MPI_Irecv.  If a request is active, clear
!> the request if it has finished.  In this way, we can receive a stream of
!> messages that we know are being sent from a specific rank.  The buffer is
!> assumed to be permanent, but we are only allowed to access the content when
!> a recv call results in a cleared request.
!===============================================================================
FUNCTION recv (self, rank, tag)
  logical:: recv
  class(mpi_comm_t):: self
  integer:: rank, tag
#ifdef MPI
  !............................................................................
  integer:: stat(MPI_STATUS_SIZE)
  !-----------------------------------------------------------------------------
  if (mp%mode == MPI_THREAD_MULTIPLE) then
    if (self%req==0) then
      call MPI_Irecv (self%buffer, size(self%buffer), MPI_INTEGER, rank, tag, &
                      mp%comm, self%req, self%ierr)
      call MPI_TEST (self%req, recv, stat, self%ierr)
    end if
  else
    !$omp critical (mpi_cr)
    if (self%req==0) then
      call MPI_Irecv (self%buffer, size(self%buffer), MPI_INTEGER, rank, tag, &
                      mp%comm, self%req, self%ierr)
      call MPI_TEST (self%req, recv, stat, self%ierr)
    end if
    !$omp end critical (mpi_cr)
  end if
#else
  recv = .false.
#endif
  if (recv) self%req = 0
END FUNCTION recv 

END MODULE mpi_comm_mod
