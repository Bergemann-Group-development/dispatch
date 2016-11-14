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
  USE timer_mod
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
    procedure:: wait_all
    procedure:: wait_for_completion
    procedure:: is_complete
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
    procedure:: delete
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
    integer:: n_unpk=0
    real:: delay_ms=0.0
    integer:: max_recv
    integer:: max_probe
    integer:: min_nq
  contains
    procedure:: init
    procedure:: get
    procedure:: sent
    procedure:: delay
    procedure:: diagnostics
  end type
  type(mpi_mesg_t), public:: mpi_mesg
  real, save:: delay_ms=3.0
  integer, save:: min_nq=10                     ! min ready q size for MPI_Test
  integer, save:: max_sent=5                    ! max outstanding patch sends
  integer, save:: max_recv=100                  ! max outstanding patch recvs
  integer, save:: max_probe=10                  ! max loop count in MPI_Probe
  integer, save:: every_send=5                  ! how often to check send_list
  integer, save:: every_recv=1                  ! how often to check recv_list
  logical, save:: wait=.false.                  ! if true, always use MPI_Wait
  logical, save:: debug=.false.
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
  integer, save:: every=0
  integer, save:: itimer=0
  real(8):: wc
  !$omp threadprivate (every)
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1 .or. .not.associated(self%head)) return
  call trace%begin ('mpi_mesg_t%check_sent', itimer=itimer)
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  mesg => self%head                                             ! top of list
  do while (self%n > max_sent)
    next => mesg%next
    call timer%end (itimer=itimer)
    call mesg%wait_all
    call self%remove (mesg)
    call timer%begin (itimer=itimer)
    mesg => next
    if (.not.associated(mesg)) then
      call trace%end (itimer)
      return
    end if
  end do
  !-----------------------------------------------------------------------------
  ! Only test for sent messages completion 'every_send' time; there is no need
  ! to rush with the deallocation of mesg buffers, and the test cost can be
  ! significant.
  !-----------------------------------------------------------------------------
  if (every>0) then
    every = every-1
    call trace%end (itimer)
    return
  else
    every = every_send
  end if
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  do while (associated(mesg))
    next => mesg%next                                           ! in case removed
    call timer%end (itimer=itimer)
    call mesg%test_all (flag)                                   ! done with?
    call timer%begin (itimer=itimer)
    if (flag) then
      call self%remove (mesg)                                   ! yes, remove
    end if
    mesg => next                                                ! true next
  end do
  call trace%end (itimer)
END SUBROUTINE check_sent

!> =============================================================================
!> Test if all requests are complete
!> =============================================================================
SUBROUTINE test_all (self, flag)
  class(mesg_t):: self
  logical:: flag
  integer:: rank, req, ierr
  integer, save:: itimer=0
  !.............................................................................
  if ((self%nreq<=0) .or. (.not.associated(self%reqs))) then
    print*,'WARNING: test_all', self%nreq, associated(self%reqs)
    return
  end if
  call trace%begin ('mesg_t%test_all', itimer=itimer)
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
  call trace_end
END SUBROUTINE test_all

!> =============================================================================
!> Wait for all requests to complete
!> =============================================================================
SUBROUTINE wait_all (self)
  class(mesg_t):: self
  integer:: rank, req, ierr
  real(8):: wc
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('mesg_t%wait_all', itimer=itimer)
  if (io%verbose>0) wc = wallclock()
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_WAITALL (self%nreq, self%reqs, MPI_STATUSES_IGNORE, ierr)
    !$omp end critical (mpi_cr)
  else
    call MPI_WAITALL (self%nreq, self%reqs, MPI_STATUSES_IGNORE, ierr)
  end if
  if (io%verbose>0) then
    wc = wallclock()-wc
    write (io_unit%log,*) wallclock(), ' wait_all:', self%id, wc
    flush (io_unit%log)
  end if
  call trace_end
END SUBROUTINE wait_all

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
  !$omp critical (mesg_list_cr)
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
  if (io%verbose>1) then
    write (io_unit%log,*) trim(self%name)//' remove mesg OK', mesg%id, self%n
    flush (io_unit%log)
  end if
  call self%delete (mesg)                       ! deallocate message
  !$omp atomic
  self%n = self%n-1                             ! decrement count
  !$omp end critical (mesg_list_cr)
END SUBROUTINE remove

!> =============================================================================
!> Deallocate the message buffer and request array, and then the message itself
!> =============================================================================
SUBROUTINE delete (self, mesg)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  deallocate (mesg%buffer)
  if (associated(mesg%reqs)) deallocate (mesg%reqs)
  deallocate (mesg)
END SUBROUTINE delete

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
  namelist /mpi_mesg_params/ min_nq, max_sent, max_probe, max_recv, every_recv, &
                             every_send, delay_ms, wait
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
  self%min_nq = min_nq
  self%max_recv = max_recv
  self%max_probe = max_probe
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
!> Check for new incoming messages.  Note that, ideally, we want to probe for
!> messages just before doing a task update, and then check afterwards, so that
!> if the delay doing the work is enough, some probed messages will have arrived.
!> This procedure returns flag=.true. as long as there are incoming messages,
!> and returns the message if comlete.   If not, the message is added to the
!> recv_list, to be check later.
!> =============================================================================
SUBROUTINE get (self, mesg)
  class(mpi_mesg_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  logical:: flag
  integer:: stat(MPI_STATUS_SIZE)
  integer:: msg ,ierr, nbuf, req
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  nullify(mesg)
  if (mpi%size <= 1) return
  call trace%begin ('mpi_mesg_t%get', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Probe for incoming messages
  !-----------------------------------------------------------------------------
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (unpk_cr)
    call probe_for
    !$omp end critical (unpk_cr)
  else
    call probe_for
  end if
  call trace%end (itimer)
  return
contains
  !-----------------------------------------------------------------------------
  ! Probe for incoming MPI messages
  !-----------------------------------------------------------------------------
  subroutine probe_for
    logical:: complete
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
    if (flag) then
      allocate (mesg)
      call MPI_GET_COUNT (stat, MPI_INT, nbuf, ierr)
      allocate (mesg%buffer(nbuf))
      call MPI_IMRECV (mesg%buffer, nbuf, MPI_INT, msg, req, ierr)
      mesg%req = req
      mesg%id = stat(MPI_TAG)
      mesg%sender = stat(MPI_SOURCE)
      mesg%nbuf = nbuf
      !$omp atomic
      mpi_mesg%n_recv = mpi_mesg%n_recv+1
    end if
  end subroutine probe_for
END SUBROUTINE get

!> =============================================================================
!> Return true or false, depending on if the mesg is complete or not.
!> =============================================================================
FUNCTION is_complete (self)
  class(mesg_t):: self
  logical:: is_complete
  !.............................................................................
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  is_complete = .false.
  if (mpi%size <= 1) return
  call trace%begin ('mesg_t%is_complete', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Check for completion
  !-----------------------------------------------------------------------------
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (complete_cr)
    call check_recv
    !$omp end critical (complete_cr)
  else
    call check_recv
  end if
  call trace%end (itimer)
  return
contains
  !-----------------------------------------------------------------------------
  ! Check for completed messages.
  !-----------------------------------------------------------------------------
  subroutine check_recv
    integer:: ierr
    integer:: stat(MPI_STATUS_SIZE)
    real(8):: wc
    !---------------------------------------------------------------------------
    if (mpi_mesg%recv_list%n > max_recv .or. wait) then
      if (io%verbose>1) then
        wc = wallclock()
        write (io_unit%log,*) wc, &
          'is_complete waiting for id', self%id, mpi_mesg%recv_list%n, wait
        flush (io_unit%log)
      end if
      call MPI_WAIT (self%req, stat, ierr)
      if (io%verbose>1) then
        wc = wallclock()-wc
        write (io_unit%log,*) wallclock(), &
          'is_complete:', self%id, wc, mpi_mesg%recv_list%n, wait
      end if
      is_complete = .true.
    else
      call MPI_TEST (self%req, is_complete, stat, ierr)
    end if
  end subroutine check_recv
END FUNCTION is_complete

!> =============================================================================
!> Return true or false, depending on if the mesg is complete or not.
!> =============================================================================
SUBROUTINE wait_for_completion (self)
  class(mesg_t):: self
  !.............................................................................
  integer:: stat(MPI_STATUS_SIZE)
  integer:: ierr
  integer, save:: itimer=0
  real(8):: wc
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1) return
  call trace%begin ('mesg_t%wait_for_completion', itimer=itimer)
  if (io%verbose>0) wc=wallclock()
  call MPI_WAIT (self%req, stat, ierr)
  if (io%verbose>0) &
    write (io_unit%log,*) 'wait_for_completion', self%id, wallclock()-wc
  call trace%end (itimer)
END SUBROUTINE wait_for_completion

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
