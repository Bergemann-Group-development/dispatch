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
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE timer_mod
  USE mpi_lib
  implicit none
  private
  !
  type, public:: mesg_t
    class(mesg_t), pointer:: next => null()
    class(mesg_t), pointer:: prev => null()
    integer, dimension(:), pointer:: buffer => null()
    integer:: nbuf, nreq=0, req=0, sender, id, ntry=0
    integer:: n_failed=0
    integer:: i=0
    integer:: seq=0
    integer:: tag=2**30
    real(8):: test_next=0d0
    real:: test_time=0.0
    logical:: anonymous=.true.
    type(lock_t), pointer:: lock => null()
    integer, pointer:: reqs(:) => null()
    integer, pointer:: fail(:) => null()
    integer, pointer:: ranks(:) => null()
  contains
    procedure:: init => init_mesg
    procedure:: dealloc => dealloc_mesg
    procedure:: send
    procedure:: irecv
    procedure, private:: test_all
    procedure:: test_sent
    procedure:: is_complete
    procedure, private:: get_id
    procedure:: request
  end type
  integer, dimension(:), allocatable, save:: expected
  integer(8) :: n_wait_all=0  , n_wait_for=0
  real(8)    :: t_wait_all=0d0, t_wait_for=0d0
  integer, save:: ident=0
  !
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
    !procedure:: reset
    procedure:: check_sent
    procedure:: print => print_list
    !procedure:: count
  end type
  !-----------------------------------------------------------------------------
  ! Data type collecting parameters and message lists into one object
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
    integer:: tag_type=0
    integer:: nbuf=0
    logical:: recv_wait=.false.
    logical:: send_priv=.true.
    logical:: debug=.false.
    logical:: uniq_mesg=.true.
    real:: test_time=0.0
  contains
    procedure:: init
    procedure:: diagnostics
    procedure:: test
    procedure:: log_files
  end type
  type(mpi_mesg_t), public:: mpi_mesg
  integer, save:: max_fail=3                    ! max MPI_TEST before resend
  integer, save:: every_send=2                  ! how often to check send_list
  logical, save:: resend=.false.                ! if true, resend late mesgs
  integer, save:: id=0
  integer, save:: verbose=0
  logical, save:: detailed_timer=.false.
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init_mesg (self)
  class(mesg_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mesg_t%init')
  if (.not. associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init ('mesg')
  else
    write (io_unit%log,*) 'mesg_t%init WARNING: lock already associated', self%id
  end if
  call trace%end ()
END SUBROUTINE init_mesg

!===============================================================================
SUBROUTINE dealloc_mesg (self)
  class(mesg_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('mesg_t%dealloc_mesg')
  if (associated(self%lock)) then
    call omp_lock%destroy (self%lock)
  end if
  if (associated(self%buffer)) then
    call io%bits_mem (-storage_size(self%buffer), &
                      product(shape(self%buffer)), '-buffer')
    deallocate (self%buffer)
  end if
  if (associated(self%reqs)) then
    call io%bits_mem (-storage_size(self%reqs), &
                      product(shape(self%reqs)), '-reqs')
    deallocate (self%reqs)
  end if
  if (associated(self%fail)) then
    call io%bits_mem (-storage_size(self%fail), &
                      product(shape(self%fail)), '-fail')
    deallocate (self%fail)
  end if
  if (associated(self%ranks)) then
    call io%bits_mem (-storage_size(self%ranks), &
                      product(shape(self%ranks)), '-ranks')
    deallocate (self%ranks)
  end if
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
!> Stub
!===============================================================================
SUBROUTINE irecv (self, rank, id)
  class(mesg_t):: self
  integer:: rank, id
  stop 'irecv'
END SUBROUTINE irecv

!===============================================================================
!> Test if all requests (sending to several ranks) are complete for this mesg
!===============================================================================
SUBROUTINE test_all (self, flag)
  class(mesg_t):: self
  logical:: flag
  integer:: rank, req, ierr
  integer, save:: itimer=0
  !.............................................................................
  if ((self%nreq <= 0) .or. (.not.associated(self%reqs))) then
    if (verbose > 2) &
      write (stdout,*) 'mesg_t%test_all: WARNING', &
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
  if (verbose > 1) then
    write (io_unit%log,'(f12.6,1x,a,2i6,l4)') wallclock(), &
      'mpi_mesg_t%test_all: id, nreq, flag =', self%id, self%nreq, flag
    flush (io_unit%log)
  end if
END SUBROUTINE test_all

!===============================================================================
!> Test if all requests (sending to several ranks) are complete for this mesg,
!> cancel and resend any messages that fails more than max_fail times.
!===============================================================================
SUBROUTINE test_and_resend (self, flag)
  class(mesg_t):: self
  logical:: flag
#ifdef MPI
  logical::  completed
  integer:: rank, req, ierr, i
  integer, save:: itimer=0
  integer:: stat(MPI_STATUS_SIZE)
  !.............................................................................
  call trace%begin('mpi_mesg_t%test_and_resend')
  if (associated(self%reqs) .and. self%nreq > 0) then  
    flag = .true.
    do i=1,self%nreq
    !-----------------------------------------------------------------------------
      if (mpi%mode == MPI_THREAD_MULTIPLE) then
        call MPI_TEST (self%reqs(i), completed, stat, ierr)
      else
        !$omp critical (mpi_cr)
        call MPI_TEST (self%reqs(i), completed, stat, ierr)
        !$omp end critical (mpi_cr)
      end if
      call mpi%assert ('is_complete', ierr)
      if (completed) then
        self%fail(i) = 0
      else
        flag = .false.
        self%fail(i) = self%fail(i)+1
        if (self%fail(i) > max_fail) then
          write (stderr,*) mpi%rank, omp%thread, 'RESEND event', self%tag
          if (verbose > 0) then
            write (io_unit%log,'(f12.6,1x,a,i6,2i5,i8)') wallclock(), &
              'test_and_resend: id, rank, max_fail, tag (resend) =', &
              self%id, self%ranks(i), self%fail(i), self%tag
          end if
          if (mpi%mode == MPI_THREAD_MULTIPLE) then
            call MPI_CANCEL (self%reqs(i), ierr)
            call mpi%assert ('mesg_t%test_and_send: MPI_CANCEL', ierr)
            call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, &
              self%ranks(i), self%tag, mpi%comm, self%reqs(i), ierr)
            call mpi%assert ('mesg_t%test_and_send: MPI_ISEND', ierr)
          else
            !$omp critical (mpi_cr)
            call MPI_CANCEL (self%reqs(i), ierr)
            call mpi%assert ('mesg_t%test_and_send: MPI_CANCEL', ierr)
            call MPI_ISEND (self%buffer, self%nbuf, MPI_INTEGER, &
              self%ranks(i), self%tag, mpi%comm, self%reqs(i), ierr)
            call mpi%assert ('mesg_t%test_and_send: MPI_ISEND', ierr)
            !$omp end critical (mpi_cr)
          end if
          self%fail(i) = 0
        else if (verbose > 1) then
          write (io_unit%log,'(f12.6,1x,a,i6,2i5,i8)') wallclock(), &
            'test_and_resend: id, rank, fail, tag =', &
            self%id, self%ranks(i), self%fail(i), self%tag
        end if
      end if
    end do
    if (verbose > 1) then
      write (io_unit%log,'(f12.6,1x,a,2i6,l4)') wallclock(), &
        'mpi_mesg_t%test_and_resend: id, nreq, flag =', &
        self%id, self%nreq, flag
      flush (io_unit%log)
    end if
  else if (verbose > 2) then
    write (stdout,*) 'mesg_t%test_and_resend: WARNING', &
      self%nreq, associated(self%reqs)
  end if
  call trace%end ()
#else
  flag = .true.
#endif
END SUBROUTINE test_and_resend

!===============================================================================
!> Logical function interface to mesg_t%test_all()
!===============================================================================
FUNCTION test_sent (self) RESULT (flag)
  class(mesg_t):: self
  logical:: flag
  !.............................................................................
  if (resend) then
    call test_and_resend (self, flag)
  else
    call test_all (self, flag)
  end if
END FUNCTION test_sent

!===============================================================================
!> Get mesg%id, mesg%sender, mesg%seq.  This set self%seq to the message 
!> sequence number, whether or net the task ID part is included in the tag.
!> For messages where the tag does not encode the task ID, mesg%id is set to 0.
!===============================================================================
SUBROUTINE get_id (self, stat)
  class(mesg_t):: self
  integer:: stat(:)
#ifdef MPI
  self%tag = stat(MPI_TAG)
  self%sender = stat(MPI_SOURCE)
  if (mpi_mesg%uniq_mesg) then
    self%seq = mod(self%tag,100)
    self%id = self%tag/100
    if (verbose > 1) &
      write (io_unit%log,*) 'mpi_mesg_t%get_id: id, sender, seq =', &
        self%id, self%sender, self%seq
  else if (verbose > 1) then
    self%id = self%tag
    write (io_unit%log,*) 'mpi_mesg_t%get_id: id, sender      =', &
      self%id, self%sender
  end if
#endif
END SUBROUTINE get_id

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
  else
    if (present(parentname) .and. detailed_timer) then
      call trace%begin ('mesg_t%is_complete('//trim(parentname)//')', 2)
    else
      call trace%begin ('mesg_t%is_complete', 2)
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
    call self%get_id (stat)
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,1x,a,i6,i5,i8)') wallclock(), &
        'mpi_mesg_t%is_complete: id, sender, tag =', &
        self%id, self%sender, self%tag
    end if
    self%n_failed = 0
    !$omp atomic
    timer%mpi_hit = timer%mpi_hit + 1
  else
    if (verbose > 1) then
      write (io_unit%log,'(f12.6,1x,a,i6,i5,i8)') wallclock(), &
        'mpi_mesg_t%is_not_complete: id, sender, tag =', &
        self%id, self%sender, self%tag
    end if
    self%n_failed = self%n_failed+1
  end if
  end block
#endif
  if (detailed_timer) then
    call trace%end (itimer)
  else
    call trace%end ()
  end if
END FUNCTION is_complete

!===============================================================================
SUBROUTINE init_list (self)
  class(mesg_list_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('msg_list_t%init_list')
  !allocate (self%lock)
  !call self%lock%init ('mlist')
  call trace%end ()
END SUBROUTINE init_list

!===============================================================================
SUBROUTINE dealloc_list (self)
  class(mesg_list_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('msg_list_t%dealloc_list')
  call omp_lock%destroy (self%lock)
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
#ifdef MPI
  nullify (mesg%next)                           ! since mesg will be tail
  if (associated(self%head)) then               ! if the list is non-empty
    self%tail%next => mesg                      ! add to the tail
    mesg%prev => self%tail                      ! prev is the previous tail
  else
    self%head => mesg                           ! else add at the head
    nullify (mesg%prev)                         ! mesg is head
  end if
  self%tail => mesg                             ! tail points to it
  !$omp atomic
  self%n = self%n+1                             ! increment counter
  if (mesg%tag == MPI_ANY_TAG) then
    !$omp atomic
    self%na = self%na+1
    !...........................................................................
    if (verbose > 0) then
      write (io_unit%log,'(f12.6,1x,a,i6)') wallclock(), &
      'mesg_t%add: increased na =', self%na
    end if
    !...........................................................................
  end if
  if (verbose > 1) then
    write (io_unit%log,*) trim(self%name)//' added mesg', mesg%id, self%n
    flush (io_unit%log)
  end if
#endif
END SUBROUTINE add

!===============================================================================
!> Remove a message from a message list. This routine is assumed to be called
!> from inside a critical region.
!===============================================================================
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
  if (verbose > 1) then
    write (io_unit%log,*) trim(self%name)//' remove mesg OK', mesg%id, self%n
    flush (io_unit%log)
  end if
  !$omp atomic
  self%n = self%n-1                             ! decrement count
END SUBROUTINE remove

SUBROUTINE check_sent (self, rank, id)
  class(mesg_list_t):: self
  integer:: rank, id
  optional:: id
END SUBROUTINE check_sent

!===============================================================================
!> Deallocate the message buffer and request array, and then the message itself
!===============================================================================
SUBROUTINE delete (self, mesg, send)
  class(mesg_list_t):: self
  class(mesg_t), pointer:: mesg
  logical, optional:: send
  !.............................................................................
  if (.not. associated(mesg)) &
    return
  if (verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,i5,i8)') wallclock(), &
      'mesg_list_t%delete: id, sender, tag =', &
      mesg%id, mesg%sender, mesg%tag
  end if
  call mesg%dealloc
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
    write (io_unit%log,*) '------------------ '//label//' ------------------'
  mesg => self%head
  do while (associated(mesg))
    write (io_unit%log,'(a,i5,i9,2i5)') ' mesg_list: '//self%name, &
      self%n, mesg%id, mesg%sender, mesg%ntry
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
  logical, save:: recv_priv
  logical, save:: recv_active
  logical, save:: send_priv
  logical, save:: queue_unpack
  logical, save:: uniq_mesg
  real, save:: test_time=20e-3
  namelist /mpi_mesg_params/ every_send, resend, &
    max_fail, test_time, uniq_mesg, debug, verbose, detailed_timer
  integer:: iostat
  character(len=120):: id = &
  '$Id$ mpi_mesg_mod.f90'
  !-----------------------------------------------------------------------------
  ! Prevent initializing more than once
  !-----------------------------------------------------------------------------
  call trace%print_id (id)
  if (self%initialized) return
  self%initialized = .true.
  call trace%begin ('mpi_mesg_t%init')
  !-----------------------------------------------------------------------------
  ! Default values from dispatcher
  !-----------------------------------------------------------------------------
  uniq_mesg = self%uniq_mesg
  send_priv = self%send_priv
  !-----------------------------------------------------------------------------
  ! Namelist input
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, mpi_mesg_params,iostat=iostat)
  write (io%output, mpi_mesg_params)
  self%debug = debug
!  self%min_nq = min_nq
!  self%max_recv = max_recv
!  self%max_probe = max_probe
!  self%recv_wait = recv_wait
  self%uniq_mesg = uniq_mesg
  mpi_mesg%test_time = test_time
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Close and reopen the sent_rrrrr_tttt_x.log files once per minute, so the last
!> 1-2 minutes of logging is always available
!===============================================================================
SUBROUTINE log_files (self)
  class(mpi_mesg_t):: self
  character(len=120):: filename
  integer:: one_two
  integer, save:: previous=-1
  !-----------------------------------------------------------------------------
  if (io%log_sent > 0) then
    one_two = wallclock()/60.
    one_two = mod(one_two,2) + 1
    if (one_two /= previous) then
      if (previous > 0) &
        close (previous)
      write (filename,'(a,"/sent_",i5.5,"_",i1,".log")') &
        trim(io%outputname), mpi%rank, one_two
      open (io_unit%sent,file=filename, form='formatted', status='unknown')
      previous = one_two
    end if
  end if
END SUBROUTINE log_files

!===============================================================================
!> Test the components
!===============================================================================
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
