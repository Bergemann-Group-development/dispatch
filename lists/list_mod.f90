!> *****************************************************************************
!> Module with list handling for generic class task_t objects
!> *****************************************************************************
MODULE list_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE load_balance_mod
  USE timer_mod
  implicit none
  private
  type, public:: list_t
    character(len=64):: name=''
    class(link_t), pointer:: head => null()
    class(link_t), pointer:: tail => null()
    class(link_t), pointer:: queue => null()
    class(link_t), pointer:: active => null()
    integer::  n=0                      ! total
    integer:: nq=0                      ! ready
    integer:: nv=0                      ! virtual
    integer:: nb=0                      ! boundary
    integer:: ni=0                      ! internal
    integer:: ne=0                      ! external
    integer:: na=0                      ! active = internal+boundary
    integer:: nac=0                     ! active tasks
    real(8):: lc(3) = +huge(0d0)        ! lower corner
    real(8):: uc(3) = -huge(0d0)        ! upper corner
    real(8):: llc(3) = +huge(0d0)       ! lower left corner
    real(8):: urc(3) = -huge(0d0)       ! upper right corner
    real(8):: position(3) = 0.0_8       ! center offset
    real(8):: size(3) = 1.0_8           ! full extent = uc-lc
    integer:: n_behind=0, n_check=0, n_ready=0, n_nbor=0
    integer:: debug=0
    integer:: n_tasks=1
    integer:: dims(3)=0
    type(lock_t):: lock
    character(len=4):: kind='list'
    logical:: face_nbors=.false.
  contains
    procedure:: init
    procedure:: update
    procedure:: reset
    procedure:: append
    procedure:: append_link
    procedure:: remove_link
    procedure:: remove_head
    procedure:: check_nbors
    procedure:: check_ready
    procedure:: check_all
    procedure:: resend_bdry
    procedure:: check_oldest
    procedure:: check_queue
    procedure:: init_nbors
    procedure:: init_all_nbors
    procedure:: remove_parents
    procedure:: intpos
    procedure:: consistency
    procedure:: qconsistency
    procedure:: aconsistency
    procedure:: statistics
    procedure:: append_list
    procedure:: remove_task
    procedure:: reset_status
    procedure:: add_by_quality
    procedure:: queue_by_time
    procedure:: queue_active
    procedure:: remove_active
    procedure:: count_status
    procedure:: print_queue
    procedure:: print_queue_until
    procedure:: print_queue_times
    procedure:: print => print_list
    procedure:: print_tasks
  end type
  integer, parameter:: mnbor=100
  public link_t, task2patch
CONTAINS

!> =============================================================================
!> Initialize list
!> =============================================================================
SUBROUTINE init (self, name)
  class(list_t):: self
  character(len=*), optional:: name
  !.............................................................................
  call trace_begin ('list_t%init')
  if (present(name)) then
    self%name = name
    write (io_unit%log,*)'init_list: ', name
  end if
  call self%lock%init (self%kind)
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Update simple list
!> =============================================================================
SUBROUTINE update (self)
  class(list_t):: self
END

!> =============================================================================
!> Reset list, by deallocating memory used by tasks, and deallocating tasks and
!> links
!> =============================================================================
SUBROUTINE reset (self)
  class(list_t):: self
  class(link_t), pointer:: link, old
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin ('list_t%reset '//self%name)
  link => self%head
  do while (associated(link))
    task => link%task
    if (associated(link)) then
      call task%deallocate
      deallocate (task)
    end if
    old => link
    link => link%next
    deallocate (old)
  end do
  nullify (self%head)
  nullify (self%tail)
  call trace_end
END SUBROUTINE reset

!> =============================================================================
!> =============================================================================
SUBROUTINE append (self, task, nbor, needed, needs_me)
  class(list_t):: self
  class(link_t), pointer, optional:: nbor
  logical, optional:: needed, needs_me
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  !.............................................................................
  call trace_begin ('list_t%append '//self%name)
  allocate (link)
  link%task => task
  if (present(nbor))     link%nbor     => nbor
  if (present(needed))   link%needed   = needed
  if (present(needs_me)) link%needs_me = needs_me
  call self%append_link (link)
  call trace_end
END SUBROUTINE append

!> =============================================================================
!> Append to the task list, counting active tasks only
!> =============================================================================
SUBROUTINE append_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(*), pointer:: anon
  logical, save:: first_time=.true.
  !.............................................................................
  call trace_begin ('list_t%append_link '//self%name)
  !$omp critical (list_cr)
  nullify(link%next)
  if (associated(self%tail)) then
    link%prev => self%tail
    link%prev%next => link
    self%tail => link
  else
    self%head => link
    self%tail => link
  end if
  self%n = self%n+1
  if      (link%task%is_set (bits%virtual)) then
    self%nv = self%nv + 1
  else if (link%task%is_set (bits%boundary)) then
    self%nb = self%nb + 1
    self%na = self%na + 1
  else if (link%task%is_set (bits%internal)) then
    self%ni = self%ni + 1
    self%na = self%na + 1
  end if
  if (associated(link%task)) then
    task => link%task
    select type (task)
    class is (patch_t)
    task%link => link
    end select
  else
    error stop 'append_link: no task associated'
  end if
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE append_link

!> =============================================================================
!> Remove a link from the list, taking special care with the head and tail
!> =============================================================================
SUBROUTINE remove_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, next
  class(task_t), pointer:: task
  !.............................................................................
  call trace_begin ('list_t%remove_link '//self%name)
  next => link%next
  !-------------------------------------------------------------------------
  ! Check for tail task
  !-------------------------------------------------------------------------
  if (associated(next)) then
    next%prev => link%prev                  ! point backwards, possibly null
  else
    self%tail => link%prev                  ! new tail link
  end if
  !-------------------------------------------------------------------------
  ! Check for head task
  !-------------------------------------------------------------------------
  if (associated(link%prev)) then           ! not head task
    link%prev%next => next                  ! make prev point to next
  else                                      ! link was head
    self%head => next                       ! make head point to next
  end if
  task => link%task
  deallocate (link)                         ! free to remove link
  self%n = self%n-1
  call trace_end
END SUBROUTINE

!> =============================================================================
!> Remove the head link
!> =============================================================================
SUBROUTINE remove_head (self)
  class(list_t):: self
  class(link_t), pointer:: head, link, next
  !.............................................................................
  call trace_begin ('list_t%remove_head '//self%name)
  head => self%head
  next => head%next
  !-------------------------------------------------------------------------
  ! Check for tail task (only one task left)
  !-------------------------------------------------------------------------
  if (associated(next)) then
    next%prev => head%prev                  ! point backwards, possibly null
  else
    self%tail => head%prev                  ! new tail link
  end if
  !-------------------------------------------------------------------------
  ! Remove the head link
  !-------------------------------------------------------------------------
  self%head => next                         ! make head point to next
  self%n = self%n-1
  call trace_end
END SUBROUTINE remove_head

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready.  If the ready bit is already set it means the patch has already
!> been put in the ready queue, and should not be checked again.
!>
!> The link pointer and everything it points to are private to this task, and
!> are not at this point in time accessible from the ready queue.
!> =============================================================================
SUBROUTINE check_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('list_t%check_nbors', itimer=itimer)
  if (io%verbose > 1) &
    write (io_unit%log,*) 'checking nbors of',link%task%id,associated(link%nbor)
  !$omp critical (queue_cr)
  task => link%task                                   ! main task
  nbor => link%nbor                                   ! first nbor
  do while (associated (nbor))                        ! keep going until end
    if (nbor%task%id==io%id_debug) print *, &
      'task', task%id,' needs task ', nbor%task%id, nbor%needs_me
    if (nbor%needs_me) then                           ! needs checking?
      if (io%verbose > 1) &
        write (io_unit%log,*) 'check nbor', nbor%task%id
      call self%check_ready (nbor%link)               ! pointer back
    else
      if (io%verbose > 1) &
        write (io_unit%log,*) 'no need to', nbor%task%id
    end if
    nbor => nbor%next                                 ! next nbor
  end do
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we check if any bits are set indicating this to be the case
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%virtual+bits%external)) then
    call self%check_ready (link)                      ! finally check link task
  end if
  !$omp end critical (queue_cr)
  call trace_end (itimer)
END SUBROUTINE check_nbors

!> =============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!> =============================================================================
SUBROUTINE check_all (self, repair)
  class(list_t):: self
  logical, optional:: repair
  class(link_t), pointer:: link
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('list_t%check_all', itimer=itimer)
  write (io_unit%log,*) 'check_all: phase 1', self%nq, associated(self%queue)
  !$omp critical (list_cr)
  if (present(repair)) then
    nullify (self%queue)
  end if
  link => self%head
  do while (associated (link))                          ! keep going until end
    call link%task%clear (bits%ready)
    call self%check_ready (link)                        ! link ready?
    link => link%next
  end do
  write (io_unit%log,*) 'check_all: don    ', self%nq, associated(self%queue)
  !if (.not.associated(self%queue)) then
  !  write (io_unit%log,*) 'check_all: phase 2'
  !  link => self%head
  !  do while (associated (link))                        ! keep going until end
  !    if (link%task%is_clear (bits%virtual+bits%external)) then
  !      !if (link%task%is_set (bits%ready)) then
  !      !  write (io_unit%log,*) 'clearing ready bit on', link%task%id
  !      !  call link%task%clear (bits%ready)
  !      !end if
  !      call self%check_ready (link)                  ! link ready?
  !    end if
  !    link => link%next
  !  end do
  !end if
  !$omp end critical (list_cr)
  if (.not.associated(self%queue)) then
    call self%check_oldest
  end if
  call trace_end (itimer)
END SUBROUTINE check_all

!> =============================================================================
!> Resend all boundary patches
!> =============================================================================
SUBROUTINE resend_bdry (self)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  link => self%head
  do while (associated(link))
    if (link%task%is_set (bits%boundary)) then
      call link%send_to_vnbors
    end if
  end do
END SUBROUTINE resend_bdry

!> =============================================================================
!> Find the oldest local and virtual patches
!> =============================================================================
SUBROUTINE check_oldest (self)
  class(list_t):: self
  class(link_t), pointer:: link, oldest, oldestv
  real(8):: told, toldv
  !.............................................................................
  call trace_begin('list_t%check_oldest '//self%name)
  write (io_unit%log,*) 'check_oldest: phase 1'
  told = 1d30
  toldv = told
  nullify (oldestv)
  !$omp critical (list_cr)
  link => self%head
  do while (associated (link))                        ! keep going until end
    if (link%task%is_set(bits%virtual)) then
      if (link%task%time < toldv) then
        toldv = link%task%time
        oldestv => link
      end if
    else
      if (link%task%time < told) then
        told = link%task%time
        oldest => link
      end if
    end if
    link => link%next
  end do
  link => oldest
  write (io_unit%log,*) 'the oldest active  task is ', oldest%task%id, oldest%task%time
  flush (io_unit%log)
  !$omp end critical (list_cr)
  !
  if (associated(oldest)) then
    write (io_unit%log,*) 'check_oldest: phase 2', associated(oldest)
    write (io_unit%log,1) 'the oldest local task is', oldest%task%id, &
      'rank', oldest%task%rank, oldest%task%time, &
      oldest%task%is_set (bits%virtual), oldest%task%is_set (bits%external)
    link => oldest%nbor
    do while (associated(link))
      write (io_unit%log,1) &
        'nbor', link%task%id, 'rank', link%task%rank, link%task%time, &
        link%task%is_set(bits%boundary), link%task%is_set (bits%virtual)
      1 format (1x,a,i9,3x,a,i6,1p,e16.6,2l3)
      link => link%next
    end do
  end if
  flush (io_unit%log)
  if (associated(oldestv)) then
    write (io_unit%log,*) 'the oldest virtual task is ', oldestv%task%id, oldestv%task%time
    write (io_unit%log,*) 'check_oldest: phase 3', associated(oldestv)
    write (io_unit%log,1) 'the oldest virtual task is', oldestv%task%id, &
      'rank', oldestv%task%rank, oldestv%task%time, &
      oldest%task%is_set (bits%virtual), oldest%task%is_set (bits%external)
    link => oldestv%nbor
    do while (associated(link))
      write (io_unit%log,1) &
        'nbor', link%task%id, 'rank', link%task%rank, link%task%time, &
        link%task%is_set(bits%boundary), link%task%is_set (bits%virtual)
      link => link%next
    end do
  end if
  flush (io_unit%log)
  call trace_end
END SUBROUTINE check_oldest

!> =============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!> =============================================================================
SUBROUTINE check_ready (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we first check if any bits are set indicating this to be the case
  !-----------------------------------------------------------------------------
  task => link%task
  if (task%is_set (bits%ready+bits%virtual+bits%external)) return
  nbor => link%nbor                                     ! start on nbor list
  do while (associated (nbor))                          ! keep going until end
    if (nbor%needed) then
      if (nbor%task%is_ahead_of(task)) then             ! and is ahead in time
        nbor => nbor%next                               ! next nbor
      else                                              ! not ahead in time
        if (io%verbose > 1) &
          write (io_unit%log,'(i5,f10.6,a,i6,2f10.6,f6.2,l4)') &
            task%id, task%time, ' failed on', nbor%task%id, nbor%task%time, nbor%task%dtime, &
            nbor%task%grace, nbor%task%is_set(bits%virtual)
        return
      end if
    else                                                ! nbor not needed
      nbor => nbor%next                                 ! skip to the next
    end if
  end do
  call self%queue_by_time (link)                        ! add task to queue
END SUBROUTINE check_ready

!> =============================================================================
!> Check if the number of nbor tasks is consistent
!> =============================================================================
SUBROUTINE consistency (self, link, i)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  if (self%debug < 2) return
  call trace_begin ('list_t%consistency',1)
  !$omp critical (list_cr)
  n = 0
  task => link%task
  nbor => link%nbor                                   ! start on nbor list
  if (associated(nbor)) id1 = nbor%task%id
  do while (associated (nbor))                        ! keep going until end
    id2 = nbor%task%id
    !write (io_unit%log,*)id2
    nbor => nbor%next                                 ! next nbor
    n = n+1
  end do
  if (n /= link%task%n_nbors) then
    write (io_unit%log,*) task%id, 'ERROR: n, nbors', n, id1, id2, i
  else
!    write (io_unit%log,*) task%id, 'consistent', id1, id2, i
  end if
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE consistency

!> =============================================================================
!> Check if the number of queued tasks is consistent; i.e., of the number of
!> tasks linked at self%queue is the same as nq
!> =============================================================================
SUBROUTINE qconsistency (self, ident)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, prev
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, ident
  real(8):: time
  !.............................................................................
  !if (self%debug < 1) return
  call trace_begin ('list_t%qconsistency',4)
  n = 0
  if (io%verbose>2) then
    link => self%head
    do while (associated(link))
      print *,'task list: id, time =', link%task%id, link%task%time
      link => link%next
    end do
    print *,'nq:', self%nq
  end if
  link => self%queue
  prev => link%task
  do while (associated(link))
    n = n+1
    if (n > self%nq) exit
    task => link%task
    time = task%time
    if (io%verbose>2) then
      print *,'queue: n, id, time =', n, link%task%id, link%task%time, &
        link%task%is_set (bits%ready), link%task%is_set (bits%busy)
    end if
    if (task%time < time) then
      write (io_unit%log,'(a,i4,2(i9,1p,g15.6,3x))') &
        'ERROR: queue out of order ident, prev%id, time, task%id', ident, prev%id, time, task%id, task%time
    end if
    prev => task
    link => link%next_time
  end do
  if (n /= self%nq) then
    write (io_unit%log,*) &
      'ERROR: qconsistency ident, n, nq, ident =', ident, n, self%nq
    n = 0
    link => self%queue
    do while (associated(link))
      n = n+1
      write (io_unit%log,*) n, link%task%id, link%task%time
      link => link%next_time
      if (n > self%nq+1) exit
    end do
    !call mpi%abort ('ERROR: qconsistency size')
  end if
  call trace_end
END SUBROUTINE qconsistency

!> =============================================================================
!> Check if the number of queued tasks is consistent; i.e., of the number of
!> tasks linked at self%queue is the same as nq
!> =============================================================================
SUBROUTINE aconsistency (self, ident)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, prev
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, ident
  real(8):: time
  !.............................................................................
  !if (self%debug < 1) return
  call trace_begin ('list_t%aconsistency',4)
  n = 0
  link => self%active
  if (associated(link)) then
    time = link%task%atime
    prev => link%task
    link => link%next_active
    n = n+1
  end if
  do while (associated(link))
    task => link%task
    n = n+1
    if (n > self%nac) exit
    if (task%atime < time) then
      write (io_unit%log,'(a,i4,2(i9,1p,g15.6,3x),i5)') &
        'ERROR: active queue out of order ident, prev%id, time, task%id, task%atime =', &
        ident, prev%id, time, task%id, task%atime, self%nac
      n = 0
      link => self%active
      do while (associated(link))
        n = n+1
        write (io_unit%log,*) n, link%task%id, link%task%atime
        link => link%next_active
      end do
      exit
    end if
    time = task%atime
    prev => task
    link => link%next_active
  end do
  if (n /= self%nac) then
    write (io_unit%log,*) &
      'ERROR: aconsistency ident, n, nac =', ident, n, self%nac
    n = 0
    link => self%active
    do while (associated(link))
      n = n+1
      write (io_unit%log,*) n, link%task%id, link%task%atime
      link => link%next_active
      if (n > self%nac+2) exit
    end do
  end if
  call trace_end
END SUBROUTINE aconsistency

!> =============================================================================
!> Set non_leaf flag on non-leaf patches
!> =============================================================================
SUBROUTINE init_nonleaf (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan
  class(link_t), pointer:: nbor
  integer:: n_add
  !.............................................................................
  return
  call trace_begin('list_t%init_nonleaf ')
  !xomp critical (list_cr)
  !-----------------------------------------------------------------------------
  ! First find and mark all non-leaf patches; these are patches that the link
  ! patch is fully contained inside -- we *assume* that finer level patches
  ! cover the whole parent patch!
  !-----------------------------------------------------------------------------
  link => self%head                                     ! start checking links
  do while (associated(link))
    scan => self%head                                   ! scan all links (FIXME)
    do while (associated(scan))                         ! until end
      if (.not.associated(scan,link).and. &             ! skip same patch
        .not.scan%task%is_set(bits%not_leaf)) then      ! skip non-leaf patch
        if (scan%task%overlaps(link%task)) then         ! if patches overlap
          if (scan%task%level==link%task%level-1) then  ! and level is one below
            if (all(abs(scan%task%position-link%task%position) < scan%task%size*0.55d0)) then
              call scan%task%set (bits%not_leaf)        ! mark non-leaf
              if (io%verbose>0) &
                write (io_unit%log,*) 'patch', link%task%id, ' parent:', scan%task%id
            end if
          end if
        end if
      end if
      scan => scan%next                                 ! continue search
    end do
    link => link%next                                   ! continue list
  end do
  !xomp end critical (list_cr)
  call trace_end
END SUBROUTINE init_nonleaf

!> =============================================================================
!> The very first initialization of a neighbor list.  After this, it can be
!> maintained without searching the entire task list on the rank.
!> =============================================================================
SUBROUTINE init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, scan
  class(link_t), pointer:: nbor
  integer:: n_add
  logical:: overlaps
  !.............................................................................
  call trace_begin('list_t%init_nbors ')
  nullify (link%nbor)                                 ! initially no nbor
  n_add = 0
  scan => self%head                                   ! scan all links (FIXME)
  do while (associated(scan))                         ! until end
    !---------------------------------------------------------------------------
    ! Add ranks that differ from the process rank to the list of load balance
    ! neighbors (each rank is only added once)
    !---------------------------------------------------------------------------
    if (scan%task%rank /= mpi%rank .and. load_balance%on) then
      call load_balance%add (scan%task%rank)
    end if
    if (.not.associated(scan,link) .and. associated(scan%task)) then
      overlaps = scan%task%overlaps(link%task)
      !-------------------------------------------------------------------------
      ! Potential nbors are other patches that overlap
      !-------------------------------------------------------------------------
      if (overlaps) then  ! if its task overlaps
        if (io%verbose>2) &
          write (io_unit%log, '(3(i9,2(3f7.3,1x)),l4)') &
          scan%task%id, scan%task%position, scan%task%size, &
          link%task%id, link%task%position, link%task%size, &
          0, scan%task%distance (link%task), self%size, &
          overlaps
        !-----------------------------------------------------------------------
        ! The ifs select to add nbors to a patch if either the patch is not a
        ! virtual patch, or the nbor is not a virtual patch, thus excluding
        ! adding virtual patches as nbors to virtual patches
        !-----------------------------------------------------------------------
        if (scan%task%is_clear(bits%virtual).or.link%task%is_clear(bits%virtual)) then
          allocate (nbor)
          nbor%task => scan%task
          nbor%link => scan
          call link%add_nbor_by_rank (nbor)
          n_add = n_add+1
          if (nbor%task%id==io%id_debug) write (io_unit%log,*) io%id_debug, 'mk is nbor:', link%task%id
          if (link%task%id==io%id_debug) write (io_unit%log,*) io%id_debug, 'mk has nbor:', nbor%task%id
        end if
      end if
    end if
    scan => scan%next                                 ! continue search
  end do
  link%task%n_nbors = n_add                           ! remember n_nbors
  if (io%verbose>1) write (io_unit%log,*) &
    'added ', n_add, ' neighbors to', link%task%id, &
    link%task%is_set(bits%boundary), link%task%is_set(bits%virtual)
  call trace_end
END SUBROUTINE init_nbors

!> =============================================================================
!> The very first initialization of the neighbor lists.  After this, it can be
!> maintained without a global search.  Distribute the work over available
!> threads.
!> =============================================================================
SUBROUTINE init_all_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link
  integer:: i, i0, i1, i2, n=0, nv=0, nb=0
  integer, save:: itimer=0
  !.............................................................................
  call trace_begin('list_t%init_all_nbors', itimer=itimer)
  call init_nonleaf (self)                              ! set non_leaf bit
  i2 = self%n/omp%nthreads + 1                          ! tasks per thread
  !$omp parallel private(i,i0,i1,link) shared(self,i2) default(none)
  i0 = omp%thread*i2
  i1 = i0 + i2
  link => self%head                                     ! start checking links
  i = 0                                                 ! task counter
  do while (associated(link))
    if (i >= i0 .and. i<i1) then                        ! thread interval
      call self%init_nbors (link)
    end if
    i = i+1                                             ! increment counter
    link => link%next                                   ! continue list
  end do
  !$omp end parallel
  link => self%head                                     ! start checking links
  do while (associated(link))
    n = n+1
    if (link%task%is_set(bits%boundary)) nb = nb+1
    if (link%task%is_set(bits%virtual)) nv = nv+1
    link => link%next                                   ! continue list
  end do
  write (io_unit%mpi,*) 'tasks: n, nb, nv =', n, nb, nv
  call trace_end (itimer)
END SUBROUTINE init_all_nbors

!> =============================================================================
!> Clear and then set the virtual and boundary bits
!> =============================================================================
SUBROUTINE reset_status (self)
  class(list_t):: self
  class(link_t), pointer:: link, nbor
  integer:: n_internal, n_boundary, n_virtual, n_external
  !-----------------------------------------------------------------------------
  ! Clear virtual and boundary bits
  !-----------------------------------------------------------------------------
  write (io_unit%log,*) 'reset_status'
  n_internal = 0
  n_boundary = 0
  n_virtual  = 0
  n_external = 0
  !$omp critical (list_cr)
  !-----------------------------------------------------------------------------
  ! Set bit
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    if (link%task%rank == mpi%rank) then
      nbor => link%nbor
      call link%task%set(bits%internal)
      call link%task%clear(bits%virtual+bits%external)
      do while (associated(nbor))
        if (nbor%task%rank/=mpi%rank) then
          call link%task%set(bits%boundary)
          call link%task%clear(bits%internal)
          exit
        end if
        nbor => nbor%next
      end do
    else
      nbor => link%nbor
      call link%task%set(bits%external)
      call link%task%clear(bits%internal+bits%boundary)
      do while (associated(nbor))
        if (link%task%id==io%id_debug) write (io_unit%log,'(a,i5,2i4,3f9.4)') &
          'mk2', nbor%task%id, link%task%rank, nbor%task%rank, link%task%position-nbor%task%position
        if (nbor%task%rank==mpi%rank) then
          call link%task%set(bits%virtual)
          call link%task%clear(bits%external)
          exit
        end if
        nbor => nbor%next
      end do
    end if
    n_internal = n_internal + merge(1,0,link%task%is_set(bits%internal))
    n_boundary = n_boundary + merge(1,0,link%task%is_set(bits%boundary))
    n_virtual  = n_virtual  + merge(1,0,link%task%is_set(bits%virtual ))
    n_external = n_external + merge(1,0,link%task%is_set(bits%external))
    link => link%next
  end do
  self%nb = n_boundary
  self%nv = n_virtual
  self%ni = n_internal
  self%na = n_internal+n_boundary
  self%n  = n_internal+n_boundary+n_virtual
  if (io%verbose>0) then
    write (io_unit%log,'(4(a,i7,5x))') &
      'rank:',mpi%rank, &
      'n_internal:', n_internal, &
      'n_boundary:', n_boundary, &
      'n_virtual:' , n_virtual
    flush (io_unit%log)
  end if
  !$omp end critical (list_cr)
END SUBROUTINE reset_status

!> =============================================================================
!> Remove parent patches, which have higher level patches inside
!> =============================================================================
SUBROUTINE remove_parents (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan, parent
  !.............................................................................
  call trace_begin('list_t%remove_parents')
  !$omp critical (list_cr)
  link => self%head                                     ! start checking links
  do while (associated(link))
    scan => self%head                                   ! scan all links (FIXME)
    do while (associated(scan))                         ! until end
      !---------------------------------------------------------------------
      ! Find parents, as patches on the next lower level that have the
      ! current patch position inside
      !---------------------------------------------------------------------
      nullify (parent)
      if (scan%task%level==link%task%level-1) then
        if (all(abs(scan%task%position-link%task%position) < scan%task%size*0.55d0)) then
          parent => scan
        end if
      end if
      scan => scan%next                                 ! continue search
      if (associated(parent)) then
        write (*,'(2(1x,a,i9,5x))') &
          'removing parent task', parent%task%id, 'level', parent%task%level
        call parent%remove_from_nbors (parent)
        call self%remove_link (parent)
      end if
    end do
    link => link%next                                   ! continue list
  end do
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE remove_parents

!> =============================================================================
!> =============================================================================
FUNCTION intpos (self, task1, task2) RESULT (out)
  class(list_t):: self
  class(task_t), pointer:: task1
  class(task_t), pointer, optional:: task2
  real:: dist(3)
  integer:: out(3)
  !.............................................................................
  if (present(task2)) then
    dist = task1%position - task2%position
    dist = modulo (dist + 0.5*self%size, self%size) - 0.5*self%size
    out = floor(dist/task1%size + 0.5)
  else
    out = floor(task1%position/task1%size + 0.5)
  end if
END FUNCTION intpos

!> =============================================================================
!> =============================================================================
SUBROUTINE statistics (self)
  class(list_t):: self
  write (io_unit%log,'(a,3(5x,a,i6))') 'output_experiment:', &
    'n_ready:', self%n_ready, &
    'n_check:', self%n_check, &
    'n_nbor:' , self%n_nbor
END SUBROUTINE statistics

!> =============================================================================
!> =============================================================================
SUBROUTINE append_list (self, task_list)
  class(list_t):: self
  class(list_t), pointer:: task_list
  class(link_t), pointer:: link
  integer:: n
  !-----------------------------------------------------------------------------
  call trace_begin ('list_t%append_list')
  if (associated(self%tail)) then
    self%tail%next => task_list%head
  else
    self%head => task_list%head
  end if
  self%tail => task_list%tail
  self%n = self%n + task_list%n
  link => self%head
  n = 0
  do while (associated(link))
    link => link%next
    n = n+1
  end do
  if (n==self%n) then
    write (io_unit%log,*) trim(self%name)//': OK append of '//task_list%name
  else
    write (io_unit%log,*) trim(self%name)//': inconsistent append of '//task_list%name
  end if
  call trace_end
END SUBROUTINE append_list

!> =============================================================================
SUBROUTINE remove_task (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link, prev
  !.............................................................................
  call trace_begin('list_t%remove_task')
  !$omp critical (list_cr)
  link => self%head
  nullify (prev)
  do while (associated(link))
    if (link%task%id == task%id) then
      if (associated(prev)) then
        prev%next => link%next
      else
        self%head => link%next
      end if
      self%n = self%n - 1
      goto 9
    end if
    prev => link
    link => link%next
  end do
  print*,'ERROR: failed to remove from task list, task', task%id
9 continue
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE remove_task

!> =============================================================================
!> Add a link node into the queue chain, in ascending time order.  When a link
!> to a task with time larger than this%task%time is found, or the end of the
!> list is fount, this link is inserted before.
!>
!> -> link0 -next_time-> link1 -next_time-> link2
!>     |                  |                  |
!>     t0                 t1                 t2
!>        tnew0   tnew1         tnew2
!>
!> In the sketch above the worst case is when two tasks are being queued at
!> nearly the same time, into the same interval.  This may be handled by having
!> a lock on each link, so when the two tasks have identified the interval
!> they should be in, which means they want to modify the next_time of the
!> previous link (link0), they need to compete about who gets the lock. The one
!> that succeeds modifies the next_time of link0 and its own next_time and then
!> unlocks the link.
!> =============================================================================
SUBROUTINE queue_by_time (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  logical:: error
  integer:: verbose
  !-----------------------------------------------------------------------------
  call real_queue_by_time (self, this, error)
  if (error) then
    io%do_trace = .true.
    verbose = io%verbose
    io%verbose = 4
    call self%qconsistency (999)
    call self%check_all (repair=.true.)
    call real_queue_by_time (self, this, error)
    if (error) then
      print *,'ERROR: queue repair failed'
      error stop
    end if
    io%verbose = verbose
    io%do_trace = .false.
  end if
END SUBROUTINE queue_by_time

SUBROUTINE real_queue_by_time (self, this, error)
  class(list_t):: self
  class(link_t), pointer:: this
  logical:: error
  class(link_t), pointer:: next, prev
  class(task_t), pointer:: task
  integer, save:: itimer=0
  integer:: nit
  !-----------------------------------------------------------------------------
  ! Refuse to queue virtual tasks and finished tasks
  !-----------------------------------------------------------------------------
  error=.false.
  task => this%task
  if (task%is_set(bits%virtual) .or. task%has_finished()) return
  call trace_begin ('list_t%queue_by_time ', itimer=itimer)
  !$omp atomic
  mpi_mesg%n_ready = mpi_mesg%n_ready+1
  !-----------------------------------------------------------------------------
  ! Lock the queue while inserting the task, and mark the task as being in the
  ! ready queue (or executing -- the bit is cleared only after updating).
  ! Double check that no other thread has set the ready bit on this task.
  !-----------------------------------------------------------------------------
  !.... critical (queue_cr)
  if (task%is_clear (bits%ready)) then
    call task%set(bits%ready)
    if (io%verbose>1) &
      write (io_unit%log,'(f12.6,i5,2x,a,i9,a,f12.6,2x,a,i5)') wallclock(), omp%thread, &
        'list_t%queue_by_time: adding task', task%id, ' at time', task%time, &
        'nq', self%nq
    nullify (prev)
    next => self%queue
    nit = 0
    do while (associated(next))
      nit = nit+1
      if (nit > self%nq) then
        print *,'ERROR: hang in queue_by_time, nit =', nit
        error = .true.
        call trace_end (itimer)
        return
      end if
      if (associated(next%task, task)) then
        write (io_unit%log,*) omp_mythread, ' WARNING: task', task%id, ' is already in ready queue'
        go to 9
      else if (next%task%time > task%time) then
        this%next_time => next
        if (associated(prev)) then
          if (io%verbose > 1) then
            write (io_unit%log,'(i4,2x,a,i6,a,i6,i5,1p,3g16.6)') &
              omp%thread, 'task',task%id,' in ready queue between',prev%task%id,next%task%id, &
              prev%task%time, this%task%time, next%task%time
          end if
          prev%next_time => this
        else
          self%queue => this
          if (io%verbose > 1) then
            write (io_unit%log,*) 'task',task%id,' at ready queue head',self%nq
          end if
        end if
        !$omp atomic
        self%nq = self%nq+1
        go to 9
      end if
      prev => next
      next => next%next_time
    end do
    self%nq = self%nq+1
    if (associated(prev)) then
      if (io%verbose > 1) then
        write (io_unit%log,*) 'task',task%id,' in ready queue after',prev%task%id,self%nq
      end if
      prev%next_time => this
    else
      if (io%verbose > 1) then
        write (io_unit%log,*) 'task',task%id,' at ready queue head',self%nq
      end if
      self%queue => this
    end if
    nullify (this%next_time)
  9 continue
    call self%remove_active (this)                      ! remove from active_queue
  !-----------------------------------------------------------------------------
  ! Since we clear the ready bit in task_list%update, outside of any critical
  ! region, it is perfectly normal that several threads compete about putting
  ! a given task onto the queue, so there is no reason to print a warning here
  !-----------------------------------------------------------------------------
  else if (io%verbose > 1) then
    if (io%verbose > 2) &
      print '(i6,a,i5,a,i6,1p,e15.6)', mpi%rank,' INFO: thread',omp%thread, &
      ' found ready bit in queue_by_time for task', task%id, task%time
    write (io_unit%log,*) wallclock(),' INFO: thread',omp%thread, &
      ' found ready bit in queue_by_time for task', task%id, task%time
    flush (io_unit%log)
  end if
  !.... end critical (queue_cr)
  call trace_end (itimer)
END SUBROUTINE real_queue_by_time

!> =============================================================================
!> =============================================================================
SUBROUTINE queue_active (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  class(task_t), pointer:: task
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin ('list_t%queue_active ', itimer=itimer)
  !-----------------------------------------------------------------------------
  task => this%task
  task%atime = task%time
  if (io%verbose>1) &
    write (io_unit%log,'(1x,i4,2x,a,i9,a,f10.5,i5,l4)') omp%thread, &
      'list_t%queue_active: adding task', task%id, ' at time', task%atime, &
      self%nac, associated(self%active)
  nullify (prev)
  !$omp critical (active_cr)
  self%nac = self%nac+1
  next => self%active
  do while (associated(next))
    if (associated(next%task, task)) then
      write (io_unit%log,*) omp_mythread, ' WARNING: task', task%id, ' is already in active queue'
      exit
    else if (next%task%atime > task%atime) then
      this%next_active => next
      if (associated(prev)) then
        if (io%verbose > 1) then
          write (io_unit%log,*) 'task',task%id,' in active queue after',prev%task%id,self%nac
        end if
        prev%next_active => this
      else
        self%active => this
        if (io%verbose > 1) then
          write (io_unit%log,*) 'task',task%id,' at ready queue head',self%nac
        end if
      end if
      exit
    end if
    prev => next
    next => next%next_active
  end do
  if (.not.associated(next)) then
    if (associated(prev)) then
      if (io%verbose > 1) then
        write (io_unit%log,*) 'task',task%id,' at active queue tail',prev%task%id,self%nac
      end if
      prev%next_active => this
    else
      if (io%verbose > 1) then
        write (io_unit%log,*) 'task',task%id,' at active queue head',self%nac
      end if
      self%active => this
    end if
    nullify (this%next_active)
  end if
  call self%aconsistency (1)
  !$omp end critical (active_cr)
  call trace_end (itimer)
END SUBROUTINE queue_active

!> =============================================================================
!> Remove a task link from the active queue
!> =============================================================================
SUBROUTINE remove_active (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%remove_active')
  !$omp critical (active_cr)
  nullify (prev)
  next => self%active
  do while (associated(next))
    if (associated(next%task, this%task)) then
      if (associated(prev)) then
        prev%next_active => next%next_active
      else
        self%active => next%next_active
      end if
      self%nac = self%nac-1
      if (io%verbose > 1) &
        write (io_unit%mpi,*) 'remove_active: id, nac =', this%task%id, self%nac
      exit
    end if
    prev => next
    next => next%next_active
  end do
  !$omp end critical (active_cr)
  call trace%end
END SUBROUTINE remove_active

!> =============================================================================
!> Add a link node into the task list, in increasing quality order
!> =============================================================================
SUBROUTINE add_by_quality (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace_begin ('list_t%add_by_quality')
  !$omp critical (list_cr)
  nullify (prev)
  next => self%head
  do while (associated(next))
    if (associated(next%task,this%task)) then
      write (io_unit%log,*)'ERROR: Trying to add task', this%task%id, ' which is already there'
      call trace_end
      goto 9
    end if
    if (next%task%level >= this%task%level) then
      !-------------------------------------------------------------------------
      ! Found a task we should be ahead of, so insert here
      !-------------------------------------------------------------------------
      this%next => next
      if (associated(prev)) then
        prev%next => this
      else
        self%head => this
      end if
      self%n = self%n+1
      call trace_end
      goto 9
    end if
    prev => next
    next => next%next
  end do
  !-----------------------------------------------------------------------------
  ! This task has the highest quality, so append it
  !-----------------------------------------------------------------------------
  if (associated(prev)) then
    prev%next => this
  else
    self%head => this
  end if
  nullify (this%next)
  self%tail => this
  self%n = self%n+1
9 continue
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE add_by_quality

!> =============================================================================
!> Print a table with patch id, positions, and all neighbors
!> =============================================================================
SUBROUTINE check_queue (self)
  class(list_t):: self
  class(link_t), pointer:: link
  real:: previous
  integer:: nq
  !.............................................................................
  !$omp critical (list_cr)
  link => self%queue
  previous = -1.0
  nq = 0
  do while (associated(link))
    if (link%task%time < previous) &
      write (io_unit%log,*)'ERROR: queue out of order at', link%task%id, link%task%time, &
      previous
    previous = link%task%time
    link => link%next_time
    nq = nq+1
  end do
  if (nq /= self%nq) write (io_unit%log,*)'ERROR: inconsistent number of tasks in queue', nq, self%nq
  !$omp end critical (list_cr)
END SUBROUTINE check_queue

!> =============================================================================
!> Insert 'new' task before 'old' task
!> =============================================================================
SUBROUTINE insert (self, old, new)
  class(list_t):: self
  class(link_t), pointer:: old, new
  !.............................................................................
  call trace_begin ('list_t%insert '//self%name)
  !$omp critical (list_cr)
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    self%head => new
  end if
  new%next => old
  new%next%prev => new
  self%n = self%n + 1
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE insert

!> =============================================================================
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE count_status (self)
  class (list_t):: self
  class(link_t) , pointer:: link, nbor
  integer:: nf
  !.............................................................................
  if (io%verbose < 2) return
  !$omp critical (list_cr)
  self%nb=0; self%nv=0; self%ni=0; self%ne=0; self%n = 0
  nf = 0
  link => self%head
  do while (associated(link))
    self%n = self%n+1
    if (link%task%is_set(bits%frozen  )) nf = nf+1
    if (link%task%is_set(bits%boundary)) self%nb = self%nb+1
    if (link%task%is_set(bits%virtual )) self%nv = self%nv+1
    if (link%task%is_set(bits%internal)) self%ni = self%ni+1
    if (link%task%is_set(bits%external)) self%ne = self%ne+1
    link => link%next
  end do
  self%na = self%ni+self%nb-nf
  !$omp end critical (list_cr)
  if (io%verbose>1) write (io_unit%log,'(a,5i6)') &
    'count_status: n[ibvea] = ', self%ni, self%nb, self%nv, self%ne, self%na
END SUBROUTINE count_status

!> =============================================================================
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_queue (self)
  class (list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  if (io%verbose < 1) return
  !$omp critical (list_cr)
  link => self%queue
  write (io_unit%log,'(a,2i4,":",$)') 'ready_queue', self%nq, self%n
  do while (associated(link))
    write (io_unit%log,'(i4,$)') link%task%id
    link => link%next_time
  end do
  !$omp end critical (list_cr)
  if (io%master) write (io_unit%log,*)''
END SUBROUTINE print_queue

!> =============================================================================
!> Print the queue until and including a given task
!> =============================================================================
SUBROUTINE print_queue_until (self, task)
  class(list_t):: self
  class(task_t):: task
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !$omp critical (list_cr)
  link => self%queue
  if (io%master) write (io_unit%log,*)'print_queue_until:'
  do while (associated(link))
    write (io_unit%log,'("queue: id, t0, t, t0-t =",i7,1p,3e12.3)') &
      link%task%id, task%time, link%task%time, task%time-link%task%time
    if (link%task%id==task%id) exit
    link => link%next_time
  end do
  !$omp end critical (list_cr)
END SUBROUTINE print_queue_until

!> =============================================================================
!> Print the queue
!> =============================================================================
SUBROUTINE print_queue_times (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !$omp critical (list_cr)
  link => self%queue
  if (present(label)) then
     write (io_unit%log,*)'print_queue_times: '//trim(label)
  else
     write (io_unit%log,*)'print_queue_times:'
  end if
  do while (associated(link))
    write (io_unit%log,'("queue: id, time =",i7,1p,g14.5,l5)') link%task%id, link%task%time, &
      link%task%is_set(bits%ready)
    link => link%next_time
  end do
  !$omp end critical (list_cr)
  write (io_unit%log,*)'done'
END SUBROUTINE print_queue_times

!> =============================================================================
!> =============================================================================
SUBROUTINE print_list (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  class(link_t), pointer:: link, prev
  class(task_t), pointer:: task
  character(len=32):: type
  !.............................................................................
  call trace_begin ('list_t%print')
  if (io%verbose>0) &
    write (io_unit%log,*) 'list_t%print: ',trim(self%name), self%n
  !$omp critical (list_cr)
  link => self%head
  do while (associated(link))
    task => link%task
    type = ''
    !select type (task)
    !type is (task_t)
    !  type = 'task_t'
    !type is (particle_t)
    !  type = 'particle_t'
    !type is (planet_t)
    !  type = 'planet_t'
    !type is (star_t)
    !  type = 'star_t'
    !type is (pebble_t)
    !  type = 'pebble_t'
    !class is (patch_t)
      type = 'patch_t'
    !end select
    if (io%verbose>1) &
      write (io_unit%log,'(a,i9,f10.6,3x,3f10.6,i4,i6,l5,2x,a)') &
        'list_t%print:', task%id, task%time, task%position, task%status, &
        task%n_nbors, associated(link%nbor), type
    link => link%next
  end do
  !$omp end critical (list_cr)
  call trace_end
END SUBROUTINE print_list

!> =============================================================================
!> Print a table with task id, positions, and all neighbors
!> =============================================================================
SUBROUTINE print_tasks (self, label)
  class (list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !if (io%verbose < 2) return
  if (present(label)) then
    write (io_unit%log,'(a,i4,":")') 'task_list: '//trim(label), self%n
  else
    write (io_unit%log,'(a,i4,":")') 'task_list', self%n
  end if
  link => self%head
  do while (associated(link))
    if (io%verbose < 3) then
      write (io_unit%log,'(i4,g15.6,l5)') link%task%id, link%task%time, link%task%is_set(bits%ready)
    else
      call link%print_nbors
    end if
    link => link%next
  end do
  write (io_unit%log,*) ''
END SUBROUTINE print_tasks

END MODULE list_mod
