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
#ifdef MPI
  include "mpif.h"
#endif
  !
  type, public:: mesg_t
    class(mesg_t), pointer:: next => null()
    class(mesg_t), pointer:: prev => null()
    integer, dimension(:), pointer:: buffer
    integer, pointer:: reqs(:) => null()
    integer:: nbuf, nreq, req=0, sender, id, ntry=0
    integer:: tag=2**30
    real(8):: test_next=0d0
  contains
    procedure:: send
    procedure:: recv
    procedure:: test_all
    procedure:: wait_all
    procedure:: wait_for_completion
    procedure:: is_complete
    procedure:: completed
    procedure:: all_completed
    procedure:: irecv
    procedure:: test_recv
  end type
  integer(8) :: n_wait_all=0  , n_wait_for=0
  real(8)    :: t_wait_all=0d0, t_wait_for=0d0
  !
  type, public:: mesg_list_t
    class(mesg_t), pointer:: head => null()
    class(mesg_t), pointer:: tail => null()
    character(len=16):: name='mesg_list'
    integer:: n=0
    integer:: id=0
  contains
    procedure:: add
    procedure:: remove_completed
    procedure:: remove
    procedure:: delete
    procedure:: reset
    procedure:: check_sent
    procedure:: print => print_list
    procedure:: count
  end type
  !
  type, public:: mpi_mesg_t
    logical:: initialized=.false.
    type(mesg_list_t):: sent_list, recv_list
    integer:: n_check=0
    integer:: n_ready=0
    integer:: n_update=0
    integer:: n_send=0
    integer:: n_recv=0
    integer:: n_delay=0
    integer:: n_unpk=0
    integer:: nq_send=0
    integer:: nq_recv=0
    real:: delay_ms=0.0
    integer:: max_recv
    integer:: max_probe
    integer:: min_nq
    logical:: wait=.false.
    logical:: recv_priv=.true.
  contains
    procedure:: init
    procedure:: get
    procedure:: sent
    procedure:: delay
    procedure:: diagnostics
    procedure:: test
  end type
  type(mpi_mesg_t), public:: mpi_mesg
  type(mesg_list_t), save, public:: unpk_list, sent_list, recv_list
  !$omp threadprivate (unpk_list,sent_list,recv_list)
  real, save:: delay_ms=3.0, test_time=20d-3
  integer, save:: min_nq=10                     ! min ready q size for MPI_Test
  integer, save:: max_sent=5                    ! max outstanding patch sends
  integer, save:: max_recv=100                  ! max outstanding patch recvs
  integer, save:: max_probe=10                  ! max loop count in MPI_Probe
  integer, save:: every_send=5                  ! how often to check send_list
  integer, save:: every_recv=1                  ! how often to check recv_list
  logical, save:: wait=.false.                  ! if true, always use MPI_Wait
  logical, save:: send_wait=.false.             ! if true, wait in OMP task
  logical, save:: send_priv=.true.              ! if true, sent_list is private
  logical, save:: recv_priv=.false.             ! if true, recv_list is private
  logical, save:: debug=.false.
  integer, save:: id=0
CONTAINS

!> =============================================================================
!> Send a buffer, already prepared with size and id tag, returning a request id
!> =============================================================================
SUBROUTINE send (self, rank, tag)
  class(mesg_t):: self
  integer, optional:: tag
  integer:: rank, req, ierr,ltag
  !.............................................................................
  call trace%begin ('mesg_t%send')
  if (present(tag)) then
    ltag = tag
  else
    ltag = self%id
  end if
  if (io%verbose>2) then
    write (io_unit%log,*) &
    'mpi_mesg::send id, nbuf, to =', self%id, self%nbuf, rank
    flush (io_unit%log)
  end if
#ifdef MPI
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, ltag, &
      mpi%comm, req, ierr)
    self%nreq = self%nreq+1                                 ! count reqs
    self%reqs(self%nreq) = req                              ! collect reqs
    !$omp end critical (mpi_cr)
  else
    call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, ltag, &
      mpi%comm, req, ierr)
    self%nreq = self%nreq+1                                 ! count reqs
    self%reqs(self%nreq) = req                              ! collect reqs
  end if
#else
  !$omp critical (mpi_cr)
  self%nreq = self%nreq+1                                 ! count reqs
  self%reqs(self%nreq) = req                              ! collect reqs
  !$omp end critical (mpi_cr)
#endif
  !$omp atomic
  mpi_mesg%n_send = mpi_mesg%n_send+1
  call trace%end
END SUBROUTINE send

!> =============================================================================
!> Recv a buffer, already prepared with size and id tag, returning a request id
!> =============================================================================
SUBROUTINE recv (self, rank, tag)
  class(mesg_t):: self
  integer, optional:: tag
  integer:: rank, req, ierr,ltag
  !.............................................................................
  call trace%begin ('mesg_t%recv')
  if (present(tag)) then
    ltag = tag
  else
    ltag = self%id
  end if
  if (io%verbose>2) then
    write (io_unit%log,*) &
    'mpi_mesg%recv id, nbuf, to =', self%id, self%nbuf, rank
    flush (io_unit%log)
  end if
  !-----------------------------------------------------------------------------
  ! Start a new receive
  !-----------------------------------------------------------------------------
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_IRECV (self%buffer, self%nbuf, MPI_INTEGER, rank, ltag, &
      mpi%comm, self%req, ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_IRECV (self%buffer, self%nbuf, MPI_INTEGER, rank, ltag, &
      mpi%comm, self%req, ierr)
    !$omp end critical (mpi_cr)
  end if
#endif
  call trace%end
END SUBROUTINE recv

!> =============================================================================
!> Remove and delete completed messages in a message list
!> =============================================================================
SUBROUTINE remove_completed (self)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  call trace%begin ('mesg_list%remove_completed')
  mesg => self%head
  do while (associated(mesg))
    if (mesg%all_completed()) then
      call self%remove (mesg)
      call self%delete (mesg)
    end if
  end do
  call trace%end
END SUBROUTINE remove_completed

!> =============================================================================
!> Check messages on a sent_list for completeness
!> =============================================================================
SUBROUTINE check_sent (self, nq)
  class(mesg_list_t):: self
  integer:: nq
  class(mesg_t), pointer:: mesg, next
  type(mesg_list_t):: sent_tmp
  logical:: flag
  integer, save:: every=0
  integer, save:: itimer=0
  real(8):: wc
  integer:: n, m
  !$omp threadprivate (every)
  !-----------------------------------------------------------------------------
  if (send_priv) then
    call check_priv
    return
  end if
  if (mpi%size <= 1 .or. .not.associated(self%head)) return
  call trace%begin ('mpi_mesg_t%check_sent', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! If there are more than max_sent messages pending, bring the number down with
  ! wait_all on the oldest ones.
  !-----------------------------------------------------------------------------
  !m = merge(max_sent,0,nq >= mpi_mesg%min_nq)
  m = max_sent
  !$omp atomic read
  n = self%n
  if (n > m) then
    !$omp critical (sent_cr)
    call sent_tmp%reset
    sent_tmp%name = 'sent_tmp'
    mesg => self%head                                             ! top of list
    do while (self%n > m)
      next => mesg%next
      call self%remove (mesg)
      call sent_tmp%add (mesg)
      mesg => next
    end do
    !$omp end critical (sent_cr)
    !---------------------------------------------------------------------------
    ! Now we can wait on those messages outside of any critical region
    !---------------------------------------------------------------------------
    write (io_unit%log,*) 'WAITALL: n,m,sent_list%n,sent_tmp%n =', n, m, self%n, sent_tmp%n
    mesg => sent_tmp%head                                             ! top of list
    do while (associated(mesg))
      next => mesg%next
      call timer%end (itimer=itimer)
      call mesg%wait_all
      call timer%begin (itimer=itimer)
      call sent_tmp%remove (mesg)
      call sent_tmp%delete (mesg, .true.)
      mesg => next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Only test for sent messages completion 'every_send' time; there is no need
  ! to rush with the deallocation of mesg buffers, and the test cost can be
  ! significant.  Note that this counted 'per thread'; every is threadprivate!
  !-----------------------------------------------------------------------------
  if (every>0) then
    every = every-1
    call trace%end (itimer)
    return
  else
    every = every_send
  end if
  !-----------------------------------------------------------------------------
  ! This must be done in a critical region, since it traverses a list that may
  ! be manipulated by other threads.  Apart from the testing itself, this is a
  ! fast traversal, since the only action is deallocation / message deletion.
  !-----------------------------------------------------------------------------
  if (associated(self%head)) then
  !$omp critical (sent_cr)
    mesg => self%head                                           ! top of list
    do while (associated(mesg))
      next => mesg%next                                         ! in case removed
      call timer%end (itimer=itimer)
      call mesg%test_all (flag)                                 ! done with?
      call timer%begin (itimer=itimer)
      if (flag) then
        call self%remove (mesg)                                 ! yes, remove
        call self%delete (mesg, .true.)                         ! deallocate
      end if
      mesg => next                                              ! true next
    end do
    !$omp end critical (sent_cr)
  end if
  call trace%end (itimer)
END SUBROUTINE check_sent

!> =============================================================================
!> Check messages on a sent_list for completeness
!> =============================================================================
SUBROUTINE check_priv
  class(mesg_t), pointer:: mesg, next
  logical:: flag
  integer:: ierr
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1 .or. .not.associated(sent_list%head)) return
  call trace%begin ('mpi_mesg_t%check_priv', itimer=itimer)
  !-----------------------------------------------------------------------------
  mesg => sent_list%head                                      ! top of list
  do while (associated(mesg))
    next => mesg%next                                         ! in case removed
    !call MPI_TESTALL (mesg%nreq, mesg%reqs, flag, MPI_STATUSES_IGNORE, ierr)
    call mesg%test_all (flag)
    if (flag) then
      call sent_list%remove (mesg)                            ! yes, remove
      call sent_list%delete (mesg, .true.)                    ! deallocate
    end if
    mesg => next                                              ! true next
  end do
  call trace%end (itimer)
END SUBROUTINE check_priv

!> =============================================================================
!> Test if all requests (sending to several ranks) are complete for this mesg
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
#ifdef MPI
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
    !$omp end critical (mpi_cr)
  else
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
  end if
#else
  flag = .true.
#endif MPI
  if (io%verbose>2) then
    write (io_unit%log,*) wallclock(), self%id, ' test_all flag =', flag
    flush (io_unit%log)
  end if
END SUBROUTINE test_all

!> =============================================================================
!> Wait for all requests related to this mesg to complete
!> =============================================================================
SUBROUTINE wait_all (self)
  class(mesg_t):: self
  integer:: rank, req, ierr
  real(8):: wc
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('mesg_t%wait_all', itimer=itimer)
  if (io%verbose >= 0) wc = wallclock()
#ifdef MPI
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_WAITALL (self%nreq, self%reqs, MPI_STATUSES_IGNORE, ierr)
    !$omp end critical (mpi_cr)
  else
    call MPI_WAITALL (self%nreq, self%reqs, MPI_STATUSES_IGNORE, ierr)
  end if
#endif MPI
  if (io%verbose >= 0) then
    wc = wallclock()-wc
    !$omp atomic
    n_wait_all = n_wait_all + 1
    !$omp atomic
    t_wait_all = t_wait_all + wc
    if (io%verbose > 0) &
      write (io_unit%log,*) wallclock(), ' wait_all:', self%id, wc
  end if
  call trace_end
END SUBROUTINE wait_all

!> =============================================================================
!> Append a mesg at the end of the message list.  Since adding and removing use
!> the same critical region to protect the list handling, there is no need for
!> atomic, but let's keep it for safe measure for now.
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
  if (io%verbose>1) then
    write (io_unit%log,*) trim(self%name)//' remove mesg OK', mesg%id, self%n
    flush (io_unit%log)
  end if
  !$omp atomic
  self%n = self%n-1                             ! decrement count
END SUBROUTINE remove

!> =============================================================================
!> Deallocate the message buffer and request array, and then the message itself
!> =============================================================================
SUBROUTINE delete (self, mesg, send)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  logical, optional:: send
  !.............................................................................
  if (associated(mesg%buffer)) deallocate (mesg%buffer)
  timer%gb = timer%gb - storage_size(mesg%buffer)*mesg%nbuf/(8.*1024**3)
  if (associated(mesg%reqs)) deallocate (mesg%reqs)
  deallocate (mesg)
  if (present(send)) then
    if (send) then
      timer%nq_send_max = max(timer%nq_send_max,mpi_mesg%nq_send)
      !$omp atomic
      mpi_mesg%nq_send = mpi_mesg%nq_send-1
    else
      timer%nq_recv_max = max(timer%nq_recv_max,mpi_mesg%nq_recv)
      !$omp atomic
      mpi_mesg%nq_recv = mpi_mesg%nq_recv-1
    end if
  end if
END SUBROUTINE delete

!> =============================================================================
!> Reset the mesg_list
!> =============================================================================
SUBROUTINE reset (self)
  class(mesg_list_t):: self
  !.............................................................................
  nullify (self%head)
  nullify (self%tail)
END SUBROUTINE reset

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
    write (io_unit%log,'(a,i5,i9,2i5)') ' mesg_list: '//self%name, &
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
                             every_send, delay_ms, wait, send_wait, send_priv, &
                             recv_priv, test_time
  integer:: iostat
  !.............................................................................
  if (self%initialized) return
  rewind (io%input); read (io%input, mpi_mesg_params,iostat=iostat)
  if (io%master) write (io%output, mpi_mesg_params)
  if (io%verbose>0) write (io_unit%log,*) 'mpi_mesg_mod::init'
  self%sent_list%name = 'sent_list'
  self%recv_list%name = 'recv_list'
  !$omp parallel
  recv_list%name = 'recv_list'
  sent_list%name = 'sent_list'
  unpk_list%name = 'unpk_list'
  sent_list%id = omp_get_thread_num()
  !$omp end parallel
  self%initialized = .true.
  self%min_nq = min_nq
  self%max_recv = max_recv
  self%max_probe = max_probe
  self%recv_priv = recv_priv
  self%wait = wait
  !call self%test
END SUBROUTINE init

!> =============================================================================
!> Test the components
!> =============================================================================
SUBROUTINE test (self)
  class(mpi_mesg_t):: self
  type (mesg_list_t):: test_list
  class (mesg_t), pointer:: mesg, next
  integer:: i, n=3
  !.............................................................................
  do i=1,n
    allocate (mesg)
    allocate (mesg%buffer(10))
    mesg%id = i
    call test_list%add (mesg)
  end do
  call test_list%print ('test1')
  mesg => test_list%head
  do i=1,n-1
    next => mesg%next
    call test_list%remove (mesg)
    call test_list%delete (mesg, .true.)
    mesg => next
  end do
  call test_list%print ('test2')
  mesg => test_list%head
  call test_list%remove (mesg)
  call test_list%delete (mesg, .true.)
  if (io%master) write (io_unit%log,*) &
    associated(test_list%head), associated(test_list%tail)
  call test_list%print ('test3')
END SUBROUTINE test

!> =============================================================================
!> Add a message to the sent_list, in a critical region, or wait for all sends
!> in a background task, or use a threadprivate list
!> =============================================================================
SUBROUTINE sent (self, mesg)
  class(mpi_mesg_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  !$omp atomic
  mpi_mesg%nq_send = mpi_mesg%nq_send+1
  if (send_priv) then
    call sent_list%add (mesg)
    if (io%verbose>0) then
      write (io_unit%log,*) 'sent:', omp_get_thread_num(), mesg%id, sent_list%n
      flush (io_unit%log)
    end if
  else if (send_wait) then
    !$omp task firstprivate(mesg)
    call mesg%wait_all
    deallocate (mesg%buffer)
    deallocate (mesg)
    !$omp atomic
    mpi_mesg%nq_send = mpi_mesg%nq_send-1
    !$omp end task
  else
    !$omp critical (sent_cr)
    call self%sent_list%add (mesg)
    !$omp end critical (sent_cr)
    if (io%verbose>0) then
      write (io_unit%log,*) 'sent:', omp_get_thread_num(), mesg%id, self%sent_list%n
      flush (io_unit%log)
    end if
  end if
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
#ifdef MPI
  integer:: stat(MPI_STATUS_SIZE)
#endif
  integer:: msg ,ierr, nbuf, req
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  nullify(mesg)
  if (mpi%size <= 1) return
  call trace%begin ('mpi_mesg_t%get', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Probe for incoming messages
  !-----------------------------------------------------------------------------
#ifdef MPI
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call probe_for
    !$omp end critical (mpi_cr)
  else
    call probe_for
  end if
#endif
  call trace%end (itimer)
  return
contains
  !-----------------------------------------------------------------------------
  ! Probe for incoming MPI messages
  !-----------------------------------------------------------------------------
  subroutine probe_for
    logical:: complete
#ifdef MPI
    call MPI_IMPROBE (MPI_ANY_SOURCE, MPI_ANY_TAG, mpi%comm, flag, msg, stat, ierr)
    if (flag) then
      allocate (mesg)
      call MPI_GET_COUNT (stat, MPI_INT, nbuf, ierr)
      allocate (mesg%buffer(nbuf))
      !$omp atomic
      mpi_mesg%nq_recv = mpi_mesg%nq_recv+1
      call MPI_IMRECV (mesg%buffer, nbuf, MPI_INT, msg, req, ierr)
      mesg%req = req
      mesg%id = stat(MPI_TAG)
      mesg%sender = stat(MPI_SOURCE)
      mesg%nbuf = nbuf
      !$omp atomic
      mpi_mesg%n_recv = mpi_mesg%n_recv+1
      !$omp atomic
      timer%bytes_recv = timer%bytes_recv + 4.0_8*nbuf
      !$omp atomic
      timer%n_recv = timer%n_recv + 1_8
    end if
#endif
  end subroutine probe_for
END SUBROUTINE get

!> =============================================================================
!> Return true or false, depending on if the mesg request is complete or not
!> =============================================================================
FUNCTION completed (self) RESULT(flag)
  class(mesg_t):: self
  logical:: flag, ierr
  integer, save:: itimer=0
#ifdef MPI
  integer:: stat(MPI_STATUS_SIZE)
#endif
  !-----------------------------------------------------------------------------
  flag = .false.
  if (mpi%size <= 1) return
  call trace%begin ('mesg_t%completed', itimer=itimer)
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_TEST (self%req, flag, stat, ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_TEST (self%req, flag, stat, ierr)
    !$omp end critical (mpi_cr)
  end if
#endif
  call trace%end (itimer)
END FUNCTION completed

!> =============================================================================
!> Return true or false, depending on if the mesg request is complete or not
!> =============================================================================
FUNCTION all_completed (self) RESULT(flag)
  class(mesg_t):: self
  logical:: flag, ierr
  integer, save:: itimer=0
#ifdef MPI
  integer:: stat(MPI_STATUS_SIZE)
#endif
  !-----------------------------------------------------------------------------
  flag = .false.
  if (mpi%size <= 1) return
  call trace%begin ('mesg_t%all_completed', itimer=itimer)
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
    !$omp end critical (mpi_cr)
  end if
#endif
  call trace%end (itimer)
END FUNCTION all_completed

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
#ifdef MPI
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (complete_cr)
    call check_recv
    !$omp end critical (complete_cr)
  else
    call check_recv
  end if
#endif
  call trace%end (itimer)
  return
contains
  !-----------------------------------------------------------------------------
  ! Check for completed messages.
  !-----------------------------------------------------------------------------
  subroutine check_recv
    integer:: ierr
#ifdef MPI
    integer:: stat(MPI_STATUS_SIZE)
#endif
    real(8):: wc
    !---------------------------------------------------------------------------
    if (mpi_mesg%recv_list%n > max_recv .or. mpi_mesg%wait) then
      if (io%verbose>1) then
        wc = wallclock()
        write (io_unit%log,*) wc, &
          'is_complete waiting for id', self%id, mpi_mesg%recv_list%n, mpi_mesg%wait
        flush (io_unit%log)
      end if
#ifdef MPI
      if (mpi%mode == MPI_THREAD_SERIALIZED) then
        !$omp critical (mpi_cr)
        call MPI_WAIT (self%req, stat, ierr)
        !$omp end critical (mpi_cr)
      else
        call MPI_WAIT (self%req, stat, ierr)
      end if
#endif
      if (io%verbose>1) then
        wc = wallclock()-wc
        write (io_unit%log,*) wallclock(), &
          'is_complete:', self%id, wc, mpi_mesg%recv_list%n, mpi_mesg%wait
      end if
      is_complete = .true.
#ifdef MPI
    else
      if (mpi%mode == MPI_THREAD_SERIALIZED) then
        !$omp critical (mpi_cr)
        call MPI_TEST (self%req, is_complete, stat, ierr)
        !$omp end critical (mpi_cr)
      else
        call MPI_TEST (self%req, is_complete, stat, ierr)
      end if
#endif
    end if
  end subroutine check_recv
END FUNCTION is_complete

!> =============================================================================
!> Issue a message recieve request from a specific rank
!> =============================================================================
SUBROUTINE irecv (self, rank, tag)
  class(mesg_t):: self
  integer:: rank, tag, ierr
  !.............................................................................
  self%sender = rank
#ifdef MPI
  call MPI_IRECV (self%buffer, self%nbuf, MPI_INT, rank, tag, mpi%comm, self%req, ierr)
#endif
END SUBROUTINE irecv

!> =============================================================================
!> Issue a message recieve request from a specific rank
!> =============================================================================
FUNCTION test_recv (self) RESULT (flag)
  class(mesg_t):: self
  integer:: ierr
  logical:: flag
#ifdef MPI
  integer:: stat(MPI_STATUS_SIZE)
  real(8):: now
  !.............................................................................
  now = wallclock()
  if (now < self%test_next) return
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_TEST (self%req, flag, stat, ierr)
    !$omp end critical (mpi_cr)
  else if (wait) then
    call MPI_WAIT (self%req, stat, ierr)
    flag = .true.
  else
    call MPI_TEST (self%req, flag, stat, ierr)
  end if
  if (flag) then
    self%id = stat(MPI_TAG)
  end if
  self%test_next = now + test_time
  !$omp atomic
  timer%n_test = timer%n_test + 1
#else
  flag = .true.
#endif
END FUNCTION test_recv

!> =============================================================================
!> Return true or false, depending on if the mesg is complete or not.
!> =============================================================================
SUBROUTINE wait_for_completion (self)
  class(mesg_t):: self
  !.............................................................................
  integer:: ierr
  integer, save:: itimer=0
  real(8):: wc
#ifdef MPI
  integer:: stat(MPI_STATUS_SIZE)
#endif
  !-----------------------------------------------------------------------------
  if (mpi%size <= 1) return
  call trace%begin ('mesg_t%wait_for_completion', itimer=itimer)
  if (io%verbose>=0) wc=wallclock()
#ifdef MPI
  if (mpi%mode == MPI_THREAD_SERIALIZED) then
    !$omp critical (mpi_cr)
    call MPI_WAIT (self%req, stat, ierr)
    !$omp end critical (mpi_cr)
  else
    call MPI_WAIT (self%req, stat, ierr)
  end if
#endif
  if (io%verbose >= 0) then
    wc = wallclock()-wc
    !$omp atomic
    n_wait_for = n_wait_for + 1
    !$omp atomic
    t_wait_for = t_wait_for + wc
    if (io%verbose > 0) &
      write (io_unit%log,*) 'wait_for_completion', self%id, wc
  end if
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
  if (mpi_mesg%delay_ms==0.0) return
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
SUBROUTINE diagnostics (self, flag)
  class(mpi_mesg_t):: self
  integer:: flag
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (flag==9) then
    !$omp critical (abort_cr)
    call trace%begin ('mpi_mesg%diagnostics', itimer=itimer)
    call self%sent_list%print ('ABORT')
    call self%recv_list%print ('ABORT')
    flush (io_unit%log)
    call trace%end (itimer)
    !$omp end critical (abort_cr)
  end if
  if (flag==1) then
    !$omp single
    write (io_unit%mpi,1) 'average wait_for time:', t_wait_for/max(n_wait_for,1_8), &
      ', with', real(n_wait_for),' waits'
    write (io_unit%mpi,1) 'average wait_all time:', t_wait_all/max(n_wait_all,1_8), &
      ', with', real(n_wait_all), ' waits'
    1 format(a,1p,2(e10.2,a))
    if (io%master.and..not.io_unit%do_validate) then
      write (io_unit%output,1) 'average wait_for time:', t_wait_for/max(n_wait_for,1_8), &
        ', with', real(n_wait_for), ' waits'
      write (io_unit%output,1) 'average wait_all time:', t_wait_all/max(n_wait_all,1_8), &
        ', with', real(n_wait_all), ' waits'
    end if
    !$omp end single
  end if
END SUBROUTINE diagnostics

END MODULE mpi_mesg_mod
