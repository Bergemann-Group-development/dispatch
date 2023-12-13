!===============================================================================
!> SENDING: Each message sent is added to a thread-private sent_list, which is
!> then checked for completed messages later, to allow deallocating the send
!> buffer.  Since this is not urgent, the sent_list doesn't have to be checked
!> very often.
!>
!> RECEIVING: With the (preferred) 'isend' method, each thread (or only the
!> master thread) handles a thread-private list of 'irecv' requests, which
!> accept messages of a predetermined size from any rank.  With the 'improbe'
!> method, the size does not have to be predetermined, but can vary from case
!> to case.  With the 'irecv' method, the thread goes through the list, looking
!> for completed messages, moving them to an unpack queue, and then refills the
!> list of irecv requests.
!===============================================================================
MODULE mpi_mesg_mod
  USE trace_mod
  USE mpi_mod
  USE os_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE timer_mod
  USE mpi_lib
  USE shared_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Data type holding various attributes of MPI messages
  !-----------------------------------------------------------------------------
  type, public:: mesg_t
    class(mesg_t), pointer:: next => null()
    class(mesg_t), pointer:: prev => null()
    integer, dimension(:), pointer:: buffer => null()
    integer:: nbuf, nreq=0, req=0, sender, id, n_try=0
    integer:: n_failed=0
    integer:: i=0
    integer:: seq=0
    integer:: tag=-1
    integer(8):: offset=0
    real(8):: test_next=0d0
    logical:: anonymous=.true.
    type(lock_t), pointer:: lock => null()
    integer, pointer:: reqs(:) => null()
    integer, pointer:: fail(:) => null()
    integer, pointer:: ranks(:) => null()
  contains
    procedure:: init => init_mesg
    procedure:: dealloc => dealloc_mesg
    procedure:: send
    procedure:: test_sent
    procedure:: is_complete
    procedure:: request
  end type
  integer, dimension(:), allocatable, save:: expected
  integer(8) :: n_wait_all=0  , n_wait_for=0
  real(8)    :: t_wait_all=0d0, t_wait_for=0d0
  integer, save:: ident=0
  !-----------------------------------------------------------------------------
  ! Data type for lists of MPI messages
  !-----------------------------------------------------------------------------
  type, public:: mesg_list_t
    class(mesg_t), pointer:: head => null()
    class(mesg_t), pointer:: tail => null()
    character(len=16):: name='mesg_list'
    integer:: n=0, na=0
    integer:: id=0
    type(lock_t), pointer:: lock
  contains
    procedure:: init => init_list
    procedure:: dealloc => dealloc_list
    procedure:: add
    procedure:: remove
    procedure:: delete
    procedure:: check
    !procedure:: reset
    procedure:: print => print_list
    !procedure:: count
  end type
  !-----------------------------------------------------------------------------
  ! Data type collecting parameters and several message lists into one object
  !-----------------------------------------------------------------------------
  type, public:: mpi_mesg_t
    logical:: initialized=.false.
    type(mesg_list_t):: sent_list, recv_list, unpk_list
    integer:: n_check=0
    integer:: n_ready=0
    integer:: n_update=0
    integer:: n_send=0
    integer:: n_recv=0
    integer:: n_delay=0
    integer:: n_unpk=0
    integer:: nq_send=0
    integer:: nq_recv=0
    integer:: tag_wrap
    integer:: nbuf=0
    logical:: recv_wait=.false.
    logical:: debug=.false.
  contains
    procedure:: init
    procedure:: diagnostics
    procedure:: test
  end type
  type(mpi_mesg_t), public:: mpi_mesg
  integer, save:: every_send=2                  ! how often to check send_list
  integer, save:: id=0
  integer, save:: verbose=0
  logical, save:: detailed_timer=.false.
CONTAINS

!===============================================================================
!> Initialize mesg with allocated lock
!===============================================================================
SUBROUTINE init_mesg (self)
  class(mesg_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mesg_t%init')
  if (.not. associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init ('mesg')
  else
    write (stdlog,*) 'mesg_t%init WARNING: lock already associated', self%id
  end if
  call trace%end ()
END SUBROUTINE init_mesg

!===============================================================================
!> Recover allocated mesg memory
!===============================================================================
SUBROUTINE dealloc_mesg (self)
  class(mesg_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mesg_t%dealloc_mesg')
  if (associated(self%lock)) then
    call self%lock%destroy
    deallocate (self%lock)
  end if
  call os%deallocate(self%buffer)
  call os%deallocate(self%reqs)
  call os%deallocate(self%fail)
  call os%deallocate(self%ranks)
  call trace%end ()
END SUBROUTINE dealloc_mesg

!===============================================================================
!> Send a buffer, already prepared with size and id tag, returning a request id
!===============================================================================
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
  self%tag = ltag
  if (verbose > 0 .or. self%id==io%id_debug) then
    write (io_unit%log,'(f12.6,1x,a,3i6,i8)') wallclock(), &
      'mpi_mesg_t%send id, nbuf, to, tag =', self%id, self%nbuf, rank, ltag
  end if
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, ltag, &
      mpi%comm, req, ierr)
    call mpi%assert ('mesg_t%send: MPI_ISEND', ierr)
    self%nreq = self%nreq+1                                 ! count reqs
    self%reqs(self%nreq) = req                              ! collect reqs
  else
    !$omp critical (mpi_cr)
    call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, rank, ltag, &
      mpi%comm, req, ierr)
    call mpi%assert ('mesg_t%send: MPI_ISEND', ierr)
    self%nreq = self%nreq+1                                 ! count reqs
    self%reqs(self%nreq) = req                              ! collect reqs
    !$omp end critical (mpi_cr)
  end if
#else
  !$omp critical (mpi_cr)
  self%nreq = self%nreq+1                                 ! count reqs
  self%reqs(self%nreq) = req                              ! collect reqs
  self%ranks(self%nreq) = rank                            ! collect ranks
  !$omp end critical (mpi_cr)
#endif
  !$omp atomic
  mpi_mesg%n_send = mpi_mesg%n_send+1
  call trace%end()
END SUBROUTINE send

!===============================================================================
!> Issue a request for the next message, with the given tag
!===============================================================================
SUBROUTINE request (self, any)
  class(mesg_t):: self
  logical, optional:: any
  integer:: ierr
  !-----------------------------------------------------------------------------
  ! Make a request, optionally allowing any source
  !-----------------------------------------------------------------------------
#ifdef MPI
  call trace%begin ('mesg_t%request')
  if (present(any)) then
    self%sender = MPI_ANY_SOURCE
  end if
  !-----------------------------------------------------------------------------
  ! Make the call in a critical region, or not
  !-----------------------------------------------------------------------------
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_IRECV (self%buffer, self%nbuf, MPI_INT, self%sender, self%tag, &
                    mpi%comm, self%req, ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_IRECV (self%buffer, self%nbuf, MPI_INT, self%sender, self%tag, &
                    mpi%comm, self%req, ierr)
    !$omp end critical (mpi_cr)
  end if
  call mpi%assert ('mesg_t%request: MPI_IRECV', ierr)
  if (verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,i5,i8)') wallclock(), &
      'mpi_mesg_t%request: id, sender, tag =', &
      self%id, self%sender, self%tag
  end if
  !-----------------------------------------------------------------------------
  ! Increment counters for number of requests queued, number received,
  ! bytes received, and total number of packages received
  !-----------------------------------------------------------------------------
  !$omp atomic
  mpi_mesg%n_recv = mpi_mesg%n_recv+1
  !$omp atomic
  timer%bytes_recv = timer%bytes_recv + 4.0_8*self%nbuf
  !$omp atomic
  timer%n_recv = timer%n_recv + 1_8
  call trace%end ()
#endif
END SUBROUTINE request

!===============================================================================
!> Test if all requests (sending to several ranks) are complete for this mesg
!===============================================================================
FUNCTION test_sent (self) RESULT (flag)
  class(mesg_t):: self
  logical:: flag
  integer:: rank, req, ierr
  integer, save:: itimer=0
  !.............................................................................
  if ((self%nreq <= 0) .or. (.not.associated(self%reqs))) then
    if (verbose > 2) &
      write (stdout,*) 'mesg_t%test_sent: WARNING', &
        self%nreq, associated(self%reqs)
    return
  end if
#ifdef MPI
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
    call mpi%assert ('MPI_TESTALL', ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_TESTALL (self%nreq, self%reqs, flag, MPI_STATUSES_IGNORE, ierr)
    call mpi%assert ('MPI_TESTALL', ierr)
    !$omp end critical (mpi_cr)
  end if
#else
  flag = .true.
#endif MPI
  if (.not.flag) self%n_try = self%n_try + 1
  if (verbose > 1) then
    write (io_unit%log,'(f12.6,1x,a,2i6,l4)') wallclock(), &
      'mpi_mesg_t%test_sent: id, nreq, flag =', self%id, self%nreq, flag
    flush (io_unit%log)
  end if
END FUNCTION test_sent

!===============================================================================
!> Return true or false, depending on if the mesg is complete or not.
!===============================================================================
FUNCTION is_complete (self, parentname)
  class(mesg_t):: self
  logical:: is_complete
  character(len=*), optional:: parentname
  !.............................................................................
  integer:: tag, ierr
  real(8):: now, test_next
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  is_complete = .true.
  if (mpi%size <= 1) then
    return
  end if
  if (detailed_timer) then
    if (present(parentname) .and. detailed_timer) then
      call trace%begin ('mesg_t%is_complete('//trim(parentname)//')', 2, itimer=itimer)
    else
      call trace%begin ('mesg_t%is_complete', 2, itimer=itimer)
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Check for completion
  !-----------------------------------------------------------------------------
#ifdef MPI
  block
  integer:: stat(MPI_STATUS_SIZE)
  if (mpi%mode == MPI_THREAD_MULTIPLE) then
    !---------------------------------------------------------------------------
    call MPI_TEST (self%req, is_complete, stat, ierr)
    call mpi%assert ('is_complete', ierr)
  else
    !$omp critical (mpi_cr)
    call MPI_TEST (self%req, is_complete, stat, ierr)
    call mpi%assert ('is_complete', ierr)
    !$omp end critical (mpi_cr)
  end if
  !$omp atomic
  timer%mpi_test = timer%mpi_test + 1
  if (is_complete) then
    self%n_failed = 0
    !$omp atomic
    timer%mpi_hit = timer%mpi_hit + 1
  else
    self%n_failed = self%n_failed+1
  end if
  end block
#endif
  if (detailed_timer) then
    call trace%end (itimer)
  end if
END FUNCTION is_complete

!===============================================================================
SUBROUTINE init_list (self)
  class(mesg_list_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('msg_list_t%init_list')
  allocate (self%lock)
  call self%lock%init ('mlist')
  call trace%end ()
END SUBROUTINE init_list

!===============================================================================
!> If the list contains messages, deallocate then, then destroy and deallocate
!> the list lock
!===============================================================================
SUBROUTINE dealloc_list (self, task_mesg)
  class(mesg_list_t):: self
  class(mesg_t), pointer, optional:: task_mesg
  class(mesg_t), pointer:: mesg, next
  character(32):: label
  !-----------------------------------------------------------------------------
  call trace%begin ('msg_list_t%dealloc_list')
  mesg => self%head
  do while (associated(mesg))
    next => mesg%next
    call mesg%dealloc ()
    if (present(task_mesg)) then
      if (associated(task_mesg,mesg)) then
        nullify (task_mesg)
      end if
    end if
    if (associated(mesg)) then
      deallocate (mesg)
    end if
    mesg => next
  end do
  if (associated(self%lock)) then
    call self%lock%destroy
    deallocate (self%lock)
  end if
  call trace%end ()
END SUBROUTINE dealloc_list

!===============================================================================
!> Append a mesg at the end of the message list.  Since adding and removing use
!> the same critical region to protect the list handling, there is no need for
!> atomic, but let's keep it for safe measure for now.
!===============================================================================
SUBROUTINE add (self, mesg)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  call trace%begin ('mesg_list_t%add')
  call                                                     self%lock%set ('add')
  nullify (mesg%next)                           ! since mesg will be tail
  if (associated(self%tail)) then               ! if the list is non-empty
    self%tail%next => mesg                      ! add to the tail
    mesg%prev => self%tail                      ! prev is the previous tail
  else
    self%head => mesg                           ! else add at the head
    nullify (mesg%prev)                         ! mesg is head
  end if
  self%tail => mesg                             ! tail points to it
  !$omp atomic
  self%n = self%n + 1                           ! increment counter
#ifdef MPI
  if (mesg%tag == MPI_ANY_TAG) then
    !$omp atomic
    self%na = self%na + 1
    !...........................................................................
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,1x,a,i6)') wallclock(), &
      'mesg_t%add: na=', self%na
    end if
  end if
#endif
  call                                                   self%lock%unset ('add')
  call trace%end ()
END SUBROUTINE add

!===============================================================================
!> Remove a message from a message list. This routine is assumed to be called
!> from inside a critical region, unless executed by the master thread.
!===============================================================================
SUBROUTINE remove (self, mesg)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  !.............................................................................
  if (associated(mesg%prev)) then               ! associated(prev) => not head
    mesg%prev%next => mesg%next                 ! cut out forwards
  else                                          ! mesg is head
    self%head => mesg%next                      ! new head
    if (associated(self%head)) then             ! if non-empty
      nullify (self%head%prev)                  ! head hair-cut for safety
    end if
  end if
  if (associated(mesg%next)) then               ! associated(next) => not tail
    mesg%next%prev => mesg%prev                 ! cut out backwards
  else                                          ! mesg is tail
    self%tail => mesg%prev                      ! new tail
    if (associated(self%tail)) then             ! if non-empty
      nullify (self%tail%next)                  ! tail hair-cut, for safety
    end if
  end if
  if (verbose > 1) then
    write (io_unit%log,*) trim(self%name)//' remove mesg OK', mesg%id, self%n
    flush (io_unit%log)
  end if
  !$omp atomic
  self%n = self%n-1                             ! decrement count
END SUBROUTINE remove

!===============================================================================
!> Deallocate the allocatable message attributes, and then the message itself
!===============================================================================
SUBROUTINE delete (self, mesg, send)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  logical, optional:: send
  logical:: send_L
  !.............................................................................
  if (.not. associated(mesg)) return
  call trace%begin ('mesg_list_t%delete')
  send_L = shared%option (send)
  if (verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,i5,2i8,i4)') wallclock(), &
      'mesg_list_t%delete: id,sender,seq,tag,n_try =', &
      mesg%id, mesg%sender, mesg%seq, mesg%tag, mesg%n_try
  end if
  call mesg%dealloc ()
  if (send_L) then
    !---------------------------------------------------------------------------
    ! Do NOT deallocate mesgs used for sending here -- they will be deallocated
    ! when the task itself is deallocated.
    !---------------------------------------------------------------------------
    timer%nq_send_max = max(timer%nq_send_max,mpi_mesg%nq_send)
    !$omp atomic
    mpi_mesg%nq_send = mpi_mesg%nq_send-1
  else
    deallocate (mesg)
    timer%nq_recv_max = max(timer%nq_recv_max,mpi_mesg%nq_recv)
    !$omp atomic
    mpi_mesg%nq_recv = mpi_mesg%nq_recv-1
  end if
  call trace%end ()
END SUBROUTINE delete

!===============================================================================
!> Check list integrity
!===============================================================================
SUBROUTINE check (self, label, log)
  class(mesg_list_t):: self
  character(*):: label
  logical:: log
  class(mesg_t), pointer:: mesg
  integer:: n
  logical:: ok
  !-----------------------------------------------------------------------------
  n = 0
  ok = .true.
  mesg => self%head
  do while (associated(mesg))
    if (verbose > 1 .or. log) &
      call output (stdlog)
    if (mesg%id /= self%id) &
      ok = .false.
    n = n + 1
    mesg => mesg%next
  end do
  if (n /= self%n) then
    ok = .false.
  end if
  if (associated(self%head)) then
    if (associated(self%head%prev)) then
      write (stderr,*) 'head%prev associated'
      ok = .false.
    end if
  end if
  if (associated(self%tail)) then
    if (associated(self%tail%next)) then
      write (stderr,*) 'tail%next associated'
      ok = .false.
    end if
  end if
  if (ok) then
    if (log) then
      write (stdlog,'(a,3i7,2x,a)') &
        'mesg_list_t%check: id,s%n,n,l=', self%id, self%n, n, label
    end if
  else
    write (stdlog,'(a,3i7,2x,a)') &
      'mesg_list_t%check: WARNING, id,s%n,n=', self%id, self%n, n, label
    write (stderr,'(a,3i7,2x,a)') &
      'mesg_list_t%check: WARNING, id,s%n,n=', self%id, self%n, n, label
    mesg => self%head
    do while (associated(mesg))
      call output (stdlog)
      call output (stderr)
      mesg => mesg%next
    end do
    flush (stdlog)
    flush (stderr)
    call io%abort ('mesg_list_t%check')
  end if
contains
  !=============================================================================
  subroutine output (unit)
  integer:: unit
  write (unit,'(a,2i7,2l2,2i7,2x,a)') &
    'mesg_list_t%check: id,assoc(reqs),assoc(buffer),n,nbuf=', mesg%id, mesg%tag, &
    associated(mesg%reqs), associated(mesg%buffer), self%n, mesg%nbuf, label
  end subroutine
END SUBROUTINE check

!===============================================================================
!> Reset the mesg_list
!===============================================================================
SUBROUTINE reset (self)
  class(mesg_list_t):: self
  !.............................................................................
  nullify (self%head)
  nullify (self%tail)
  self%n = 0
END SUBROUTINE reset

!===============================================================================
!> Print a list of messages
!===============================================================================
SUBROUTINE print_list (self, label)
  class (mesg_list_t):: self
  character(len=*), optional:: label
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  if (present(label)) &
    write (io_unit%mpi,*) '------------------ '//label//' ------------------'
  mesg => self%head
  do while (associated(mesg))
    write (io_unit%mpi,'(a,i5,i9,2i5)') ' mesg_list: '//self%name, &
      self%n, mesg%id, mesg%sender, mesg%n_try
    mesg => mesg%next
  end do
END SUBROUTINE print_list

!===============================================================================
!> Count messages, with warning if the %n does not agree with the actual number
!===============================================================================
SUBROUTINE count (self, label)
  class (mesg_list_t):: self
  character(len=*):: label
  class(mesg_t), pointer:: mesg
  integer:: n
  !-----------------------------------------------------------------------------
  if (verbose < 1) return
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

!===============================================================================
!> Initialize the three message lists
!===============================================================================
SUBROUTINE init (self)
  class(mpi_mesg_t):: self
  !.............................................................................
  logical, save:: debug
  integer, save:: tag_wrap=10
  integer:: iostat
  namelist /mpi_mesg_params/ tag_wrap, debug, verbose, detailed_timer
  !-----------------------------------------------------------------------------
  ! Prevent initializing more than once
  !-----------------------------------------------------------------------------
  if (self%initialized) return
  self%initialized = .true.
  call trace%begin ('mpi_mesg_t%init')
  !-----------------------------------------------------------------------------
  ! Namelist input
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, mpi_mesg_params,iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, mpi_mesg_params)
  end if
  call io%namelist_warning ('mpi_mesg_params', iostat)
  write (io%namelist, mpi_mesg_params)
  self%debug = debug
  self%tag_wrap = tag_wrap
  if (mpi%size > 1) then
    call self%test ()
  endif
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Test the components
!===============================================================================
SUBROUTINE test (self)
  class(mpi_mesg_t):: self
  type (mesg_list_t):: test_list
  class (mesg_t), pointer:: mesg, next, save1, save2
  integer:: i, n=3
  !-----------------------------------------------------------------------------
  ! Allocate two messages, deallocate the first one, and allocate a new
  !-----------------------------------------------------------------------------
  write (io_unit%mpi,*) 'mpi_mesg_t%test:'
  call alloc_mesg (self, mesg)
  save1 => mesg
  call alloc_mesg (self, mesg)
  save2 => mesg
  write (io_unit%mpi,*) associated(save1), associated(save2), &
    associated(mesg,save1), associated(mesg,save2)
  deallocate (save1)
  write (io_unit%mpi,*) associated(save1), associated(save2), &
    associated(mesg,save1), associated(mesg,save2)
  call alloc_mesg (self, mesg)
  write (io_unit%mpi,*) associated(save1), associated(save2), &
    associated(mesg,save1), associated(mesg,save2)
  !-----------------------------------------------------------------------------
  ! Allocate n messages and add them to a list, then remove and delete n-1
  !-----------------------------------------------------------------------------
  call test_list%init ()
  do i=1,n
    call alloc_mesg (self, mesg)
    call os%allocate (mesg%buffer, 10)
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
  write (io_unit%mpi,*) &
    associated(test_list%head), associated(test_list%tail), &
    associated(test_list%head,test_list%tail)
  call test_list%print ('test3')
END SUBROUTINE test

SUBROUTINE alloc_mesg (self, mesg)
  class(mpi_mesg_t):: self
  class(mesg_t), pointer:: mesg
  !-----------------------------------------------------------------------------
  allocate (mesg)
  call mesg%init ()
END SUBROUTINE alloc_mesg

!===============================================================================
!> Diagnostic printout, called when task_list_mod::update stalls
!===============================================================================
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
    1 format(1x,a,1p,2(e10.2,a))
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
