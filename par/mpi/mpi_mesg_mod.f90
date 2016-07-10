!> =============================================================================
!> $Id$
!> =============================================================================
MODULE mpi_mesg_mod
  USE trace_mod
  USE mpi_mod
  USE io_mod
  USE io_unit_mod
  USE omp_timer_mod
  implicit none
  private
  include "mpif.h"
  !
  type, public:: mesg_t
    class(mesg_t), pointer:: next => null()
    class(mesg_t), pointer:: prev => null()
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
  contains
    procedure:: add
    procedure:: remove
    procedure:: check_sent
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
    procedure:: sent
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
  call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, self%id, &
    mpi%comm, req, ierr)
  self%nreq = self%nreq+1                                   ! count reqs
  self%reqs(self%nreq) = req                                ! collect reqs
END SUBROUTINE send

!> =============================================================================
!> Check messages on a sent_list for completeness
!> =============================================================================
SUBROUTINE check_sent (self)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg, next
  logical:: flag
  !.............................................................................
  mesg => self%head                                             ! top of list
  do while (associated(mesg))
    next => mesg%next                                           ! in case removed
    call mesg%test_all (flag)                                   ! done with?
    if (flag) then
      call self%remove (mesg)                                   ! yes, remove
      deallocate (mesg%buffer, mesg%reqs)
      deallocate (mesg)
    end if
    mesg => next                                                ! true next
  end do
END SUBROUTINE check_sent

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
  if (debug.and.io%verbose>2) then
    write (io_unit%log,*) wallclock(), self%id, ' test_all flag =', flag
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
    mesg%prev => self%tail                      ! prev is the previous tail
  else
    self%head => mesg                           ! else add at the head
    nullify (mesg%prev, mesg%next)
  end if
  self%tail => mesg                             ! tail points to it
  self%n = self%n+1                             ! increment counter
  if (debug.and.io%verbose>1) then
    write (io_unit%log,*) trim(self%name)//' added mesg', mesg%id, self%n
    flush (io_unit%log)
  end if
END SUBROUTINE add

!> =============================================================================
!> Remove a message from a message list. This routine is assumed to be called
!> from inside a critical region.
!> =============================================================================
SUBROUTINE remove (self, mesg)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  if (associated(mesg%prev)) then               ! associated(prev) => not head
    mesg%prev%next => mesg%next                 ! cut out forwards
  else                                          ! mesg is head
    self%head => mesg%next                      ! new head
  end if
  if (associated(mesg%next)) then               ! associated(next) => not tail
    mesg%next%prev => mesg%prev                 ! cut out backwards
  else                                          ! mesg is tail
    self%tail => mesg%prev                      ! new tail
  end if
  self%n = self%n-1                             ! decrement count
  if (io%verbose>1) then
    write (io_unit%log,*) trim(self%name)//' remove mesg OK', mesg%id, self%n
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
  self%sent_list%name = 'sent_list'
  self%recv_list%name = 'recv_list'
  self%unpack_list%name = 'unpk_list'
  self%initialized = .true.
END SUBROUTINE init

!> =============================================================================
!> Add a message to the sent_list, in a critical region
!> =============================================================================
SUBROUTINE sent (self, mesg)
  class(mpi_mesg_t):: self
  class(mesg_t), pointer:: mesg, next
  !-----------------------------------------------------------------------------
  !$omp critical (sent_cr)
  call self%sent_list%add (mesg)
  !$omp end critical (sent_cr)
END SUBROUTINE sent

!> =============================================================================
!> Check for new incoming messages
!> =============================================================================
SUBROUTINE check (self)
  class(mpi_mesg_t):: self
  integer:: stat(MPI_STATUS_SIZE)
  integer:: msg ,ierr, nbuf, req, nmsg, i
  logical:: flag
  class(mesg_t), pointer:: mesg
  integer, save:: nprint=50
  integer, save:: itimer(5)=0
  integer, save:: every=10
  real(8):: wc(12)
  real(8),save:: ws(11)
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1) return
  call trace_begin ('mpi_mesg_mod::check(1)', itimer=itimer(1))
  wc(1) = wallclock()
  if (every>0) then
    every = every-1
    wc(2:3) = wc(1)
  else
    every = 10
    !$omp critical (sent_cr)
    wc(2) = wallclock()
    call self%sent_list%check_sent
    wc(3) = wallclock()
    !$omp end critical (sent_cr)
  end if
  call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
  nmsg = 0
  wc(4) = wallclock()
  do while (flag)
    nmsg = nmsg+1
    allocate (mesg)
    call MPI_GET_COUNT (stat, MPI_INT, nbuf, ierr)
    allocate (mesg%buffer(nbuf))
    call MPI_IMRECV (mesg%buffer, nbuf, MPI_INT, msg, req, ierr)
    mesg%req = req
    mesg%id = stat(MPI_TAG)
    mesg%sender = stat(MPI_SOURCE)
    !$omp critical (recv_cr)
    call self%recv_list%add (mesg)
    !$omp end critical (recv_cr)
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
  end do
  wc(5) = wallclock()
  if (io%verbose>2) then
    if (nmsg==0) then
      if (nprint>0) then
        nprint = nprint-1
        write (io_unit%log,*) 'mpi_mesg_mod::check nmsg =', nmsg
        flush (io_unit%log)
      else if (nprint==0) then
        nprint = nprint-1
        write (io_unit%log,*) 'mpi_mesg_mod::check continuing...'
        flush (io_unit%log)
      end if
    else
      if (nprint==-1) nprint = 10
      write (io_unit%log,*) 'mpi_mesg_mod::check nmsg =', nmsg
      flush (io_unit%log)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Check the recv_list for completed receives, adding them to unpack_list
  !-----------------------------------------------------------------------------
  wc(6) = wallclock()
  if (mod(every,2)==0) then
    wc(7) = wc(6)
    !$omp critical (recv_cr)
    wc(8) = wallclock(); ws(7) = ws(7)+(wc(8)-wc(7)); wc(7) = wc(8) 
    mesg => self%recv_list%head
    i = 0
    do while (associated(mesg))
      i = i+1
      req = mesg%req
      call MPI_TEST (req, flag, stat, ierr)
      wc(8) = wallclock(); ws(8) = ws(8)+(wc(8)-wc(7)); wc(7) = wc(8) 
      if (flag) then
        call self%recv_list%remove (mesg)
        wc(8) = wallclock(); ws(9) = ws(9)+(wc(8)-wc(7)); wc(7) = wc(8) 
        !$omp critical (unpk_cr)
        call self%unpack_list%add (mesg)
        !$omp end critical (unpk_cr)
        wc(8) = wallclock(); ws(10) = ws(10)+(wc(8)-wc(7)); wc(7) = wc(8) 
        if (io%verbose>1) then
          write (io_unit%log,1) wallclock(), ' recv   complete', mesg%id, &
            i, self%recv_list%n, self%unpack_list%n
          flush (io_unit%log)
        end if
      else
        if (io%verbose>1) then
          write (io_unit%log,1) wallclock(), ' recv incomplete', mesg%id, &
            i, self%recv_list%n, self%unpack_list%n
          1 format(f15.6,2x,a,i9,3i5)
          flush (io_unit%log)
        end if
      end if
      mesg => mesg%next
    end do
    !$omp end critical (recv_cr)
  end if
  wc(7) = wallclock()
  do i=1,6
    !$omp atomic
    ws(i)=ws(i)+(wc(i+1)-wc(i))
  end do
  if (io%verbose>1) write (io_unit%log,'(a,3p,10f10.3)') 'wc:', ws(1:10)
  call trace_end (itimer(1))
END SUBROUTINE check

!> =============================================================================
!> Get the head mesg from the unpack_list, if non-empty, and move the head.
!> =============================================================================
SUBROUTINE get (self, mesg)
  class(mpi_mesg_t):: self
  class(mesg_t), pointer:: mesg, next
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1) return
  !$omp critical (unpk_cr)
  call trace_begin('mpi_mesg_mod::get', itimer=itimer)
  call self%unpack_list%print
  mesg => self%unpack_list%head
  if (associated(mesg)) then
    self%unpack_list%head => mesg%next
    if (associated(self%unpack_list%head)) &
      nullify(self%unpack_list%head%prev)
    self%unpack_list%n = self%unpack_list%n - 1
    if (debug .and. io%verbose>2) then
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
  call trace_end (itimer)
  !$omp end critical (unpk_cr)
END SUBROUTINE get

!> =============================================================================
!> =============================================================================
SUBROUTINE print_list (self)
  class (mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  if (io%verbose < 3) return
  mesg => self%head
  do while (associated(mesg))
    write (io_unit%log,*) 'mesg_list: '//self%name, mesg%id, mesg%sender
    mesg => mesg%next
  end do
END SUBROUTINE print_list

END MODULE mpi_mesg_mod
