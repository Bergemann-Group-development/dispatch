!> =============================================================================
!> $Id$
!> =============================================================================
MODULE mpi_mesg_mod
  USE trace_mod
  USE mpi_mod
  USE omp_lock_mod
  USE io_mod
  USE io_unit_mod
  implicit none
  private
  include "mpif.h"
  !
  type, public:: mesg_t
    class(mesg_t), pointer:: next => null()
    integer, dimension(:), pointer:: buffer
    integer, pointer:: reqs(:) => null()
    integer:: nbuf, nreq, req, sender, id
  contains
    procedure:: send
    procedure:: test_all
  end type
  !
  type, public:: mesg_list_t
    class(mesg_t), pointer:: head => null()
    class(mesg_t), pointer:: tail => null()
    character(len=16):: name
    integer:: n=0
    type(lock_t):: lock
  contains
    procedure:: add
    procedure:: remove
    procedure:: print => print_list
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
  logical, save:: debug=.true.
CONTAINS

!> =============================================================================
!> Send a buffer, already prepared with size and id tag, returning a request id
!> =============================================================================
SUBROUTINE send (self, rank)
  class(mesg_t):: self
  integer:: rank, req, ierr
  !.............................................................................
  if (debug) then
    write (io_unit%log,*) ' send', self%id, ' to', rank
    flush (io_unit%log)
  end if
  call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, self%id, &
    mpi%comm, req, ierr)
  self%nreq = self%nreq+1                                   ! count reqs
  self%reqs(self%nreq) = req                                ! collect reqs
END SUBROUTINE send

!> =============================================================================
!> Test if all requests are complete
!> =============================================================================
SUBROUTINE test_all (self, flag)
  class(mesg_t):: self
  logical:: flag
  integer:: rank, req, ierr
  !.............................................................................
  if ((self%nreq<=0) .or. (.not.associated(self%reqs))) then
    print*,'WARNING: test_all', self%nreq, associated(self%reqs)
    return
  end if
  call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
  if (debug.and.io%verbose>1) then
    write (io_unit%log,*) ' test_all flag =', flag
    flush (io_unit%log)
  end if
END SUBROUTINE test_all

!> =============================================================================
!> Append a mesg at the end of the message list
!> =============================================================================
SUBROUTINE add (self, mesg)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  nullify (mesg%next)
  if (associated(self%head)) then               ! if the list is non-empty
    self%tail%next => mesg                      ! add to the tail
  else
    self%head => mesg                           ! else add at the head
  end if
  self%tail => mesg                             ! tail points to it
  self%n = self%n+1                             ! increment counter
  if (debug.and.io%verbose>1) then
    write (io_unit%log,*) trim(self%name)//' added mesg', mesg%id, self%n
    flush (io_unit%log)
  end if
END SUBROUTINE add

!> =============================================================================
!> Remove a message from a message list, and deallocate the beffer and reqs
!> The list has no prev pointers, so we need to search for the message first.
!> This routine is assumed to be called from inside a critical region.
!> =============================================================================
SUBROUTINE remove (self, this)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: this, mesg, next
  !.............................................................................
  mesg => self%head                                     ! list start
  if (associated(mesg,this)) then
    self%head => mesg%next                              ! skip this
    deallocate (this%buffer)                            ! free to dealloc
    if (associated(this%reqs)) deallocate (this%reqs)   ! send mesgs have reqs
    self%n = self%n-1                                   ! count down
    if (debug) then
      write (io_unit%log,*) trim(self%name)//' remove head OK', this%id, self%n
      flush (io_unit%log)
    end if
    nullify(self%head,self%tail)
    return
  end if
  do while (associated(mesg))                           ! as long as ..
    next => mesg%next                                   ! tentative
    if (associated(next,this)) then                     ! is it this?
      next => next%next                                 ! skip this
      mesg%next => next                                 ! bypass this
      if (associated(self%head,this)) self%head => next ! case this==head
      if (associated(self%tail,this)) self%tail => mesg ! case this==tail
      deallocate (this%buffer)                          ! free to dealloc
      if (associated(this%reqs)) deallocate (this%reqs) ! send mesgs have reqs
      self%n = self%n-1                                 ! count down
      if (debug) then
        write (io_unit%log,*) trim(self%name)//' remove mesg OK', this%id, self%n
        flush (io_unit%log)
      end if
      return                                            ! done
    end if
    mesg => next                                        ! next in list
  end do
  if (debug) then
    write (io_unit%log,*) trim(self%name)//' remove mesg failed', this%id, self%n
    flush (io_unit%log)
  end if
END SUBROUTINE remove

!> =============================================================================
!> Initialize the three message lists
!> =============================================================================
SUBROUTINE init (self)
  class(mpi_mesg_t):: self
  !.............................................................................
  if (self%initialized) return
  if (debug) write (io_unit%log,*) 'mpi_mesg_mod::init'
  call self%sent_list%lock%init ('sent')
  call self%recv_list%lock%init ('recv')
  call self%unpack_list%lock%init ('unpk')
  self%sent_list%name = 'sent_list'
  self%recv_list%name = 'recv_list'
  self%unpack_list%name = 'unpk_list'
  self%initialized = .true.
END SUBROUTINE init

!> =============================================================================
!> Check for new incoming messages
!> =============================================================================
SUBROUTINE check (self)
  class(mpi_mesg_t):: self
  integer:: stat(MPI_STATUS_SIZE)
  integer:: msg ,ierr, nbuf, req, nmsg
  logical:: flag
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1) return
  call trace_begin ('mpi_mesg_mod::check')
  if (debug.and.io%verbose>1) write (io_unit%log,*) 'mpi_mesg_mod::check'
  call self%init
  call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
  !$omp critical (recv_cr)
  nmsg = 0
  do while (flag)
    nmsg = nmsg+1
    !call trace_tag ('flag')
    allocate (mesg)
    call MPI_GET_COUNT (stat, MPI_INT, nbuf, ierr)
    allocate (mesg%buffer(nbuf))
    call MPI_IMRECV (mesg%buffer, nbuf, MPI_INT, msg, req, ierr)
    mesg%req = req
    mesg%id = stat(MPI_TAG)
    mesg%sender = stat(MPI_SOURCE)
    if (debug.and.io%verbose>1) then
      write (io_unit%log,*) ' check mesh add', mesg%id
      flush (io_unit%log)
    end if
    call self%recv_list%add (mesg)
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
  end do
  if (io%verbose>1) then
    write (io_unit%log,*) 'nmsg =', nmsg
    flush (io_unit%log)
  end if
  !-----------------------------------------------------------------------------
  ! Check the recv_list for completed receives, adding them to unpack_list
  !-----------------------------------------------------------------------------
  !call trace_tag ('check')
  mesg => self%recv_list%head
  do while (associated(mesg))
    !call trace_tag ('mesg')
    req = mesg%req
    call MPI_TEST (req, flag, stat, ierr)
    if (flag) then
      self%recv_list%head => mesg%next
      self%recv_list%n = self%recv_list%n - 1
      if (debug.and.io%verbose>1) then
        write (io_unit%log,*) ' recv complete', mesg%id
        flush (io_unit%log)
      end if
      call self%unpack_list%add (mesg)
    else
      if (debug.and.io%verbose>1) then
        write (io_unit%log,*) ' check incomplete', mesg%id
        flush (io_unit%log)
      end if
    end if
    mesg => mesg%next
  end do
  !$omp end critical (recv_cr)
  call trace_end
END SUBROUTINE check

!> =============================================================================
!> Get the head mesg from the unpack_list, if non-empty, and move the head.
!> =============================================================================
SUBROUTINE get (self, mesg)
  class(mpi_mesg_t):: self
  class(mesg_t), pointer:: mesg, next
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1) return
  call trace_begin('mpi_mesg_mod::get')
  call self%unpack_list%print
  !$omp critical (recv_cr)
  mesg => self%unpack_list%head
  if (associated(mesg)) then
    self%unpack_list%head => mesg%next
    self%unpack_list%n = self%unpack_list%n - 1
    if (debug) then
      write (io_unit%log,*) ' get ', mesg%id, ' from', mesg%sender, ' n', self%unpack_list%n
      next => self%unpack_list%head
      if (associated(next)) then
        write (io_unit%log,*) ' next ', mesg%id, ' from', mesg%sender
      else
        nullify (self%unpack_list%tail)
      end if
      flush (io_unit%log)
    end if
  end if
  call self%unpack_list%print
  !$omp end critical (recv_cr)
  call trace_end
END SUBROUTINE get

!> =============================================================================
!> =============================================================================
SUBROUTINE print_list (self)
  class (mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  mesg => self%head
  do while (associated(mesg))
    write (io_unit%log,*) 'mesg_list: '//self%name, mesg%id, mesg%sender
    mesg => mesg%next
  end do
END SUBROUTINE print_list

END MODULE mpi_mesg_mod
