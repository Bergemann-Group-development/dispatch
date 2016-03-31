!> =============================================================================
!> $Id$
!> =============================================================================
MODULE mpi_mesg_mod
  USE trace_mod
  USE mpi_mod
  USE omp_lock_mod
  implicit none
  private
  include "mpif.h"
  !
  type, public:: mesg_t
    type(mesg_t), pointer:: next => null()
    integer, dimension(:), pointer:: buffer
    integer, pointer:: reqs(:)
    integer:: nbuf, nreq, req, sender, tag
  contains
    procedure:: send
    procedure:: test_all
  end type
  !
  type, public:: mesg_list_t
    type(mesg_t), pointer:: head => null()
    type(mesg_t), pointer:: tail => null()
    integer:: n=0
    type(lock_t):: lock
  contains
    procedure:: add
    procedure:: remove
  end type
  !
  type, public:: mpi_mesg_t
    logical:: initialized=.false.
    type(mesg_list_t):: sent_list, recv_list, unpack_list
  contains
    procedure:: init
    procedure:: check
    procedure:: get
  end type
  type(mpi_mesg_t), public:: mpi_mesg
CONTAINS

!> =============================================================================
!> Send a buffer, already prepared with size and tag, returning a request id
!> =============================================================================
SUBROUTINE send (self, rank)
  class(mesg_t):: self
  integer:: rank, req, ierr
  !.............................................................................
  call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, self%tag, &
    mpi%comm, req, ierr)
  self%nreq = self%nreq+1                                   ! count reqs
  self%reqs(self%nreq) = req                                ! collect reqs
END SUBROUTINE send

!> =============================================================================
!> Test if all requests are comlete
!> =============================================================================
SUBROUTINE test_all (self, flag)
  class(mesg_t):: self
  logical:: flag
  integer:: rank, req, ierr
  !.............................................................................
  call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
END SUBROUTINE test_all

!> =============================================================================
!> Append a mesg at the end of the message list
!> =============================================================================
SUBROUTINE add (self, mesg)
  class(mesg_list_t):: self
  type(mesg_t), pointer:: mesg
  !.............................................................................
  nullify (mesg%next)
  if (associated(self%tail)) self%tail%next => mesg
  self%tail => mesg
END SUBROUTINE add

!> =============================================================================
!> Remove the message at the head of the message list, and deallocate its
!> allocatable parts
!> =============================================================================
SUBROUTINE remove (self, mesg)
  class(mesg_list_t):: self
  type(mesg_t), pointer:: mesg
  !.............................................................................
  self%head => mesg%next
  deallocate (mesg%buffer)
  deallocate (mesg%reqs)
END SUBROUTINE remove

!> =============================================================================
!> Initialize the three message lists
!> =============================================================================
SUBROUTINE init (self)
  class(mpi_mesg_t):: self
  !.............................................................................
  if (self%initialized) return
  call self%sent_list%lock%init ('mpi_mesg')
  call self%recv_list%lock%init ('mpi_mesg')
  call self%unpack_list%lock%init ('mpi_mesg')
  self%initialized = .true.
END SUBROUTINE init

!> =============================================================================
!> Check for new incoming messages
!> =============================================================================
SUBROUTINE check (self)
  class(mpi_mesg_t):: self
  integer:: stat(MPI_STATUS_SIZE)
  integer:: msg ,ierr, nbuf, req
  logical:: flag
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace_begin ('mpi_mesg_mod::check')
  call self%init
  call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
  do while (flag)
    call trace_tag ('flag')
    allocate (mesg)
    call MPI_GET_COUNT (stat, MPI_INT, nbuf, ierr)
    allocate (mesg%buffer(nbuf))
    call MPI_IMRECV (mesg%buffer, nbuf, MPI_INT, msg, req, ierr)
    mesg%req = req
    mesg%tag = stat(MPI_TAG)
    mesg%sender = stat(MPI_SOURCE)
    call self%recv_list%lock%set ('check')
    self%recv_list%tail%next => mesg
    self%recv_list%tail => mesg
    self%recv_list%n = self%recv_list%n + 1
    call self%recv_list%lock%unset ('check')
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
  end do
  !-----------------------------------------------------------------------------
  ! Check the recv_list for completed receives
  !-----------------------------------------------------------------------------
  call trace_tag ('check')
  call self%recv_list%lock%set ('check')
  mesg => self%recv_list%head
  do while (associated(mesg))
    call trace_tag ('mesg')
    req = mesg%req
    call MPI_TEST (req, flag, stat, ierr)
    if (flag) then
      self%recv_list%head => mesg%next
      self%recv_list%n = self%recv_list%n - 1
      call self%unpack_list%lock%set ('check')
      self%unpack_list%tail%next => mesg
      self%unpack_list%tail => mesg
      self%unpack_list%n = self%unpack_list%n + 1
      call self%unpack_list%lock%unset ('check')
    end if
    mesg => mesg%next
  end do
  call self%recv_list%lock%unset ('check')
  call trace_end
END SUBROUTINE check

!> =============================================================================
!> Return the first buffer from the unpack_list, if non-empty
!> =============================================================================
SUBROUTINE get (self, buffer)
  class(mpi_mesg_t):: self
  integer, dimension(:), pointer:: buffer
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call self%unpack_list%lock%set ('get')
  mesg => self%unpack_list%head
  if (associated(mesg)) then
    buffer => mesg%buffer
    self%unpack_list%head => mesg%next
    self%unpack_list%n = self%unpack_list%n - 1
    call self%unpack_list%lock%unset ('get')
    deallocate (mesg)
  else
    call self%unpack_list%lock%unset ('get')
    nullify(buffer)
  end if
END SUBROUTINE get

END MODULE mpi_mesg_mod
