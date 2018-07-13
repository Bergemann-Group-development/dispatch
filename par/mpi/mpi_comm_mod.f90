!> =============================================================================
!> Simple module for sending / receiving packages to / from other ranks
!> =============================================================================
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

!> =============================================================================
!> Send a buffer to a rank, with a given tag, and ignore the request
!> =============================================================================
SUBROUTINE send (self, rank, tag)
  class(mpi_comm_t):: self
  integer:: rank, tag
  !-----------------------------------------------------------------------------
#ifdef MPI
  call MPI_Isend (self%buffer, size(self%buffer), MPI_INTEGER, rank, tag, &
                  mp%comm, self%req, self%ierr)
#endif MPI
END SUBROUTINE send

!> =============================================================================
!> If not send request is active, start a send.  When finished, clear request.
!> =============================================================================
FUNCTION sent (self)
  logical:: sent
  class(mpi_comm_t):: self
#ifdef MPI
  !............................................................................
  integer:: stat(MPI_STATUS_SIZE)
  !-----------------------------------------------------------------------------
  call MPI_TEST (self%req, sent, stat, self%ierr)
#else
  sent = .false.
#endif
  if (sent) self%req = 0
END FUNCTION sent

!> =============================================================================
!> If no request is active, issue an MPI_Irecv.  If a request is active, clear
!> the request if it has finished.  In this way, we can receive a stream of
!> messages that we know are being sent from a specific rank.  The buffer is
!> assumed to be permanent, but we are only allowed to access the content when
!> a recv call results in a cleared request.
!> =============================================================================
FUNCTION recv (self, rank, tag)
  logical:: recv
  class(mpi_comm_t):: self
  integer:: rank, tag
#ifdef MPI
  !............................................................................
  integer:: stat(MPI_STATUS_SIZE)
  !-----------------------------------------------------------------------------
  if (self%req==0) then
    call MPI_Irecv (self%buffer, size(self%buffer), MPI_INTEGER, rank, tag, &
                    mp%comm, self%req, self%ierr)
  end if
  call MPI_TEST (self%req, recv, stat, self%ierr)
#else
  recv = .false.
#endif
  if (recv) self%req = 0
END FUNCTION recv 

END MODULE mpi_comm_mod
