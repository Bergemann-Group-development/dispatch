!> =============================================================================
!> $Id$
!> =============================================================================
MODULE mpi_mesg_mod
  USE trace_mod
  USE mpi_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
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
    integer:: nbuf, nreq, req, sender, id, ntry=0
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
    procedure:: count
  end type
  !
  type, public:: mpi_mesg_t
    logical:: initialized=.false.
    type(mesg_list_t):: sent_list, recv_list, unpack_list
    integer:: n_check=0
    integer:: n_ready=0
    integer:: n_update=0
    integer:: n_send=0
    integer:: n_recv=0
    integer:: n_delay=0
    real:: delay_ms=0.0
  contains
    procedure:: init
    procedure:: check
    procedure:: get
    procedure:: sent
    procedure:: delay
    procedure:: diagnostics
  end type
  type(mpi_mesg_t), public:: mpi_mesg
  real, save:: delay_ms=3.0
  integer, save:: max_probe=5                   ! max loop count in MPI_Probe
  integer, save:: every_send=5                  ! how often to check send_list
  integer, save:: every_recv=1                  ! how often to check recv_list
  logical, save:: debug=.false.
  logical, save:: wait=.false.
CONTAINS

!> =============================================================================
!> Send a buffer, already prepared with size and id tag, returning a request id
!> =============================================================================
SUBROUTINE send (self, rank)
  class(mesg_t):: self
  integer:: rank, req, ierr
  !.............................................................................
  call trace%begin ('mesg_t%send')
  if (io%verbose>1) then
    write (io_unit%log,*) &
    'mpi_mesg::send id, nbuf, to =', self%id, self%nbuf, rank
    flush (io_unit%log)
  end if
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, self%id, &
      mpi%comm, req, ierr)
    self%nreq = self%nreq+1                                 ! count reqs
    self%reqs(self%nreq) = req                              ! collect reqs
    !$omp end critical (mpi_cr)
  else
    call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, self%id, &
      mpi%comm, req, ierr)
    self%nreq = self%nreq+1                                 ! count reqs
    self%reqs(self%nreq) = req                              ! collect reqs
  end if
  !$omp atomic
  mpi_mesg%n_send = mpi_mesg%n_send+1
  call trace%end
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
    else if (io%verbose>0) then
      mesg%ntry = mesg%ntry+1
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
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
    !$omp end critical (mpi_cr)
  else
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
  end if
  if (io%verbose>2) then
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
  !$omp atomic
  self%n = self%n+1                             ! increment counter
  if (io%verbose>1) then
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
!> Print a list of messages
!> =============================================================================
SUBROUTINE print_list (self, label)
  class (mesg_list_t):: self
  character(len=*), optional:: label
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  if (present(label)) &
    write (io_unit%log,*) '------------------ '//label//' ------------------'
  mesg => self%head
  do while (associated(mesg))
    write (io_unit%log,'(a,i5,i9,2i5)') 'mesg_list: '//self%name, &
      self%n, mesg%id, mesg%sender, mesg%ntry
    mesg => mesg%next
  end do
END SUBROUTINE print_list

!> =============================================================================
!> Count messages, with warning if the %n does not agree with the actual number
!> =============================================================================
SUBROUTINE count (self, label)
  class (mesg_list_t):: self
  character(len=*):: label
  class(mesg_t), pointer:: mesg
  integer:: n
  !-----------------------------------------------------------------------------
  if (io%verbose < 1) return
  mesg => self%head
  n = 0
  do while (associated(mesg))
    n = n+1
    mesg => mesg%next
  end do
  if (n /= self%n) then
    write (io_unit%log,*) 'WARNING: inconsistent '//self%name, n, self%n
  end if
END SUBROUTINE count

!> =============================================================================
!> Initialize the three message lists
!> =============================================================================
SUBROUTINE init (self)
  class(mpi_mesg_t):: self
  namelist /mpi_mesg_params/ every_recv, every_send, delay_ms, max_probe, wait
  integer:: iostat
  !.............................................................................
  if (self%initialized) return
  rewind (io%input); read (io%input, mpi_mesg_params,iostat=iostat)
  if (io%master) write (io%output, mpi_mesg_params)
  if (io%verbose>0) write (io_unit%log,*) 'mpi_mesg_mod::init'
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
SUBROUTINE check (self, wait_req)
  class(mpi_mesg_t):: self
  logical:: wait_req
  integer:: stat(MPI_STATUS_SIZE)
  integer:: msg ,ierr, nbuf, req, nmsg, i, nr, nu
  logical:: flag
  class(mesg_t), pointer:: mesg
  real(8):: wc(12)
  integer, save:: itimer=0
  integer, save:: every=0
  integer, save:: probe
  real(8), save:: ws(11)
  !$omp threadprivate (ws, every)
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1) return
  call trace_begin ('mpi_mesg_mod::check', itimer=itimer)
  wc(1) = wallclock()
  !-----------------------------------------------------------------------------
  ! Only test for sent messages completin 'every_send' time; there is no need
  ! to rush with the deallocation of mesg buffers, and the test cost can be
  ! significant.
  !-----------------------------------------------------------------------------
  if (every>0) then
    every = every-1
    wc(2:3) = wc(1)
  else
    every = every_send
    !$omp critical (sent_cr)
    wc(2) = wallclock()
    call self%sent_list%check_sent
    wc(3) = wallclock()
    !$omp end critical (sent_cr)
  end if
  !-----------------------------------------------------------------------------
  ! Probe for incoming messages
  !-----------------------------------------------------------------------------
  wc(4) = wallclock()
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call probe_mesg
    !$omp end critical (mpi_cr)
  else
    call probe_mesg
  end if
  call trace_end (itimer)
  return
contains
  !-----------------------------------------------------------------------------
  ! Probe for incoming MPI messages, and if they are already complete add them
  ! directly to the unpack list, otherwise add them to the recv list
  !-----------------------------------------------------------------------------
  subroutine probe_mesg
    if (io%verbose>1) then
      !$omp atomic read
      nr  = self%recv_list%n
      !$omp atomic read
      nu = self%unpack_list%n
      write (io_unit%log,*) 'mpi_mesg::probe_mesg recv%n, unpk% =', nr, nu, mpi%mode
      flush (io_unit%log)
    end if
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
    nmsg = 0
    do while (flag)
      nmsg = nmsg+1
      allocate (mesg)
      call MPI_GET_COUNT (stat, MPI_INT, nbuf, ierr)
      allocate (mesg%buffer(nbuf))
      call MPI_IMRECV (mesg%buffer, nbuf, MPI_INT, msg, req, ierr)
      mesg%req = req
      mesg%id = stat(MPI_TAG)
      mesg%sender = stat(MPI_SOURCE)
      mesg%nbuf = nbuf
      call MPI_WAIT (req, stat, ierr)
      mpi_mesg%n_recv = mpi_mesg%n_recv+1
!      if (io%verbose>1) then
        write (io_unit%log,'(a,i9,2i4,i4,3x,2l1)') &
        'mpi_mesg::recv id, recv%n, unpk%n, fr =', &
        mesg%id, mpi_mesg%recv_list%n, mpi_mesg%unpack_list%n, mesg%sender, flag, wait_req
        flush (io_unit%log)
!      end if
      !$omp critical (unpk_cr)
      call self%unpack_list%add (mesg)
      !$omp end critical (unpk_cr)
      if (nmsg > max_probe) exit
      call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
    end do
  end subroutine probe_mesg
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
  !call self%unpack_list%print
  mesg => self%unpack_list%head
  if (associated(mesg)) then
    self%unpack_list%head => mesg%next
    if (associated(self%unpack_list%head)) &
      nullify(self%unpack_list%head%prev)
    !$omp atomic
    self%unpack_list%n = self%unpack_list%n - 1
    !call self%unpack_list%count ('get')
    if (io%verbose>0) then
      write (io_unit%log,'(a,i9,i6,2i4)') 'unpack msg id, sender, revc%n, unpk%n:', mesg%id, mesg%sender, &
        self%recv_list%n, self%unpack_list%n
    else if (io%verbose>2) then
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
  call trace_end (itimer)
  !$omp end critical (unpk_cr)
END SUBROUTINE get

!> =============================================================================
!> Delay for about 1 ms, to encourage incoming messages to complete
!> =============================================================================
SUBROUTINE delay (self, n)
  class(mpi_mesg_t):: self
  integer:: n
  real(8):: wc
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mpi_mesg_t%delay', itimer=itimer)
  wc = wallclock()
  write (io_unit%log,*) wc,'delay:',n
  flush (io_unit%log)
  do while ((wallclock()-wc) < 1e-3*mpi_mesg%delay_ms)
  end do
  !$omp atomic
  self%n_delay = self%n_delay + 1
  call trace%end (itimer)
END SUBROUTINE delay

!> =============================================================================
!> Diagnostic printout, called when task_list_mod::update stalls
!> =============================================================================
SUBROUTINE diagnostics (self, onoff)
  class(mpi_mesg_t):: self
  integer:: onoff
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  debug = onoff>0
  if (debug) then
    !$omp critical (abort_cr)
    call trace%begin ('mpi_mesg%diagnostics', itimer=itimer)
    call self%sent_list%print ('ABORT')
    call self%recv_list%print ('ABORT')
    flush (io_unit%log)
    call trace%end (itimer)
    !$omp end critical (abort_cr)
  end if
END SUBROUTINE diagnostics

END MODULE mpi_mesg_mod
