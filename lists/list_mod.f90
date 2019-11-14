!*******************************************************************************
!> Module with list handling for generic class task_t objects
!*******************************************************************************
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
  !USE load_balance_mod
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
    integer:: verbose=0
    type(lock_t):: lock
    character(len=4):: kind='list'
    logical:: face_nbors=.false.
    logical:: detailed_timer=.false.
  contains
    procedure:: init
    procedure:: init_bdries
    procedure:: reset
    procedure:: append_task
    procedure:: append_link
    generic:: append => append_task, append_link
    procedure:: prepend_link
    procedure:: remove_task
    procedure:: remove_link
    generic:: remove => remove_task, remove_link
    procedure:: remove_and_reset
    procedure:: update_counts
    procedure:: count_status
    procedure:: reset_status
    procedure:: check_nbors
    procedure:: check_ready
    procedure:: check_all
    procedure:: resend_bdry
    procedure:: check_oldest
    procedure:: check_queue
    procedure:: add_new_link
    procedure:: init_nbors
    procedure:: set_init_nbors
    procedure:: init_nbor_nbors
    procedure:: check_nbor_nbors
    procedure:: refresh_nbors
    procedure:: init_all_nbors
    procedure:: remove_parents
    procedure:: intpos
    procedure:: consistency
    procedure:: qconsistency
    procedure:: aconsistency
    procedure:: statistics
    procedure:: append_list
    procedure:: add_by_quality
    procedure:: queue_by_time
    procedure:: queue_active
    procedure:: remove_active
    procedure:: print_queue
    procedure:: print_queue_until
    procedure:: print_queue_times
    procedure:: print => print_list
    procedure:: print_tasks
    procedure, nopass:: check_nbor_list
    procedure:: give_to
    procedure:: update_nbor_status
    procedure:: update_status
    procedure:: send_to_vnbors
    procedure:: test_nbor_status
    procedure:: info
  end type
  integer, dimension(:), allocatable:: sequence
  integer, parameter:: mnbor=100
  type(list_t), public:: list
  public link_t, task2patch
CONTAINS

!===============================================================================
!> Initialize list
!===============================================================================
SUBROUTINE init (self, name)
  class(list_t):: self
  character(len=*), optional:: name
  !.............................................................................
  call trace%begin ('list_t%init')
  if (present(name)) then
    self%name = name
    write (io_unit%log,*)'init_list: ', name
    call self%lock%init (name(1:4))
  else
    write (io_unit%log,*)'init_list: list_t'
    call self%lock%init ('list_t')
  end if
  !$omp critical (sequence_alloc_cr)
  if (.not.allocated(sequence)) then
    allocate (sequence(0:mpi%size-1))
    sequence(:) = 0
  end if
  !$omp end critical (sequence_alloc_cr)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Reset list, by deallocating memory used by tasks, and deallocating tasks and
!> links
!===============================================================================
SUBROUTINE reset (self)
  class(list_t):: self
  class(link_t), pointer:: link, old
  class(task_t), pointer:: task
  !.............................................................................
  call trace%begin ('list_t%reset '//self%name)
  link => self%head
  do while (associated(link))
    task => link%task
    if (associated(link)) then
      call task%dealloc
      deallocate (task)
    end if
    old => link
    link => link%next
    deallocate (old)
  end do
  nullify (self%head)
  nullify (self%tail)
  self%n = 0
  call trace%end()
END SUBROUTINE reset

!===============================================================================
!> Allocate a new link and append to the list, optionally setting attributes
!===============================================================================
SUBROUTINE append_task (self, task, nbor, needed, needs_me)
  class(list_t):: self
  class(link_t), pointer, optional:: nbor
  logical, optional:: needed, needs_me
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  !.............................................................................
  call trace%begin ('list_t%append_task '//self%name)
  allocate (link)
  call link%init
  link%task => task
  if (present(nbor))     link%nbor     => nbor
  if (present(needed))   link%needed   = needed
  if (present(needs_me)) link%needs_me = needs_me
  call self%append_link (link)
  call trace%end()
END SUBROUTINE append_task

!===============================================================================
!> Append to the task list, counting active tasks only
!===============================================================================
SUBROUTINE append_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(*), pointer:: anon
  logical, save:: first_time=.true.
  !.............................................................................
  call trace%begin ('list_t%append_link ')
  call link%init
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' wait for tasklist(2)'
  call self%lock%set ('list_t%append_link')
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(2)'
  nullify(link%next)
  if (associated(self%tail)) then
    link%prev => self%tail
    link%prev%next => link
    self%tail => link
  else
    self%head => link
    self%tail => link
  end if
  if (associated(link%task)) then
    task => link%task
    select type (task)
    class is (patch_t)
    task%link => link
    end select
  else
    call io%abort('append_link: no task associated')
  end if
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  !-----------------------------------------------------------------------------
  self%n = self%n + 1
  if (self%verbose > 2) &
    print *,'list_t%append: (1) n =', trim(self%name), self%n, link%task%id
  call self%update_counts (link, +1)
  if (self%verbose > 1) &
    print *,'list _t%append:(2) n =', trim(self%name), self%n, link%task%id
  call self%lock%unset ('list_t%append_link')
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(2)'
  call trace%end()
END SUBROUTINE append_link

!===============================================================================
!> Prepend to the task list, counting active tasks only
!===============================================================================
SUBROUTINE prepend_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(*), pointer:: anon
  logical, save:: first_time=.true.
  !.............................................................................
  call trace%begin ('list_t%prepend_link '//self%name)
  call link%init
  call self%lock%set ('list_t%prepend_link')
  nullify(link%prev)
  if (associated(self%head)) then
    link%next => self%head
    link%next%prev => link
    self%head => link
  else
    self%head => link
    self%tail => link
  end if
  if (associated(link%task)) then
    task => link%task
    select type (task)
    class is (patch_t)
    task%link => link
    end select
    call self%update_counts (link, +1)
  else
    call io%abort('prepend_link: no task associated')
  end if
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  !-----------------------------------------------------------------------------
  call self%lock%unset ('list_t%prepend_link')
  call trace%end()
END SUBROUTINE prepend_link

!===============================================================================
!> Remove a task from consideration, making sure that virtual copies do the same.
!> Set bits%remove, which excludes the task from nbor lists in check_ready(),
!> and signals rank nbors to repeat this remove_and_reset().
!===============================================================================
SUBROUTINE remove_and_reset (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  integer:: count
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%remove_and_reset')
  !-----------------------------------------------------------------------------
  ! By setting bits%remove immediately, we make sure that the task is no longer
  ! considered in check_ready() calls, and also that the task is no longer
  ! counted in check_ready() calls to nbors of it.
  !-----------------------------------------------------------------------------
  call link%task%set (bits%remove)            ! cf. check_ready() and unpack()
  !-----------------------------------------------------------------------------
  ! Send a suicide note to virtual nbors, if this is a boundary patch
  !-----------------------------------------------------------------------------
  if (link%task%is_set (bits%boundary)) then
    call self%send_to_vnbors (link)
  end if
  !-----------------------------------------------------------------------------
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,2x,a,i6,2l3)') &
      wallclock(), 'list_t%remove_and_reset: id =', link%task%id, &
      link%task%is_set(bits%boundary), link%task%is_set(bits%virtual)
  end if
  !-----------------------------------------------------------------------------
  call self%lock%set ('remove_and_reset')  ! lock task list while changing
  call self%remove_link (link)             ! remove link from task list
  call self%lock%unset ('remove_and_reset')! release the task list
  call self%check_nbors (link)             ! check nbors, now link is removed
  call self%set_init_nbors (link)          ! trigger new nbor lists on nbors
  !$omp atomic update
  timer%levelpop (link%task%level) = &     ! decrement level population count
  timer%levelpop (link%task%level) - 1
  call link%garbage_collect (link)         ! add to garbage bin
  call trace%end()
END SUBROUTINE remove_and_reset

!===============================================================================
!> Remove a link from the list, taking special care with the head and tail
!===============================================================================
SUBROUTINE remove_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, next
  class(task_t), pointer:: task
  !.............................................................................
  call trace%begin ('list_t%remove_link '//self%name)
  next => link%next
  !-----------------------------------------------------------------------------
  ! Check for tail task
  !-----------------------------------------------------------------------------
  if (associated(next)) then
    !call next%qlock%set ('remove_link')
    next%prev => link%prev                  ! point backwards, possibly null
    !call next%qlock%unset ('remove_link')
  else
    !call self%tail%qlock%set ('remove_link')
    self%tail => link%prev                  ! new tail link
    !call self%tail%qlock%unset ('remove_link')
  end if
  !-----------------------------------------------------------------------------
  ! Check for head task
  !-----------------------------------------------------------------------------
  if (associated(link%prev)) then           ! not head task
    !call link%prev%qlock%set ('remove_link')
    link%prev%next => next                  ! make prev point to next
    !call link%prev%qlock%unset ('remove_link')
  else                                      ! link was head
    !call self%head%qlock%set ('remove_link')
    self%head => next                       ! make head point to next
    !call self%head%qlock%unset ('remove_link')
  end if
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  !-----------------------------------------------------------------------------
  call self%update_counts (link, -1)
  call trace_end
END SUBROUTINE remove_link

!===============================================================================
!> Updating the status counts using the delta of a single task, which is faster
!> but less robust than counting all tasks.  Should be thoroughly validated
!> before put into use, which is also not motivated unless count_status() takes
!> a noticeable time.
!===============================================================================
SUBROUTINE update_counts (self, link, delta)
  class(list_t):: self
  class(link_t), pointer:: link
  integer:: delta, n, na
  !-----------------------------------------------------------------------------
  ! The counter updates below do not need to be atomic if the list is locked,
  ! but the penalty is minimal
  !-----------------------------------------------------------------------------
  if      (link%task%is_set (bits%virtual)) then
    !$omp atomic update
    self%nv = self%nv + delta
  else if (link%task%is_set (bits%boundary)) then
    !$omp atomic update
    self%nb = self%nb + delta
    !$omp atomic update
    self%na = self%na + delta
  else if (link%task%is_set (bits%internal)) then
    !$omp atomic update
    self%ni = self%ni + delta
    !$omp atomic update
    self%na = self%na + delta
  end if
  !-----------------------------------------------------------------------------
  ! Validate using the incremental method, when task_list_t%verbose > 0
  !-----------------------------------------------------------------------------
  if (self%verbose > 0) then
    n  = self%n
    na = self%na
    call self%count_status ()
    if (n /= self%n .or. na /= self%na) then
      write (stdout,'(a,2(2i6,2x))') &
        'list_t%update_counts WARNING: inconsistent counts', &
        n, self%n, na, self%na
      self%n  = n
      self%na = na
    end if
  end if
END SUBROUTINE update_counts 

!===============================================================================
!> Count the number of internal, boundary, virtual, frozen, and external tasks
!===============================================================================
SUBROUTINE count_status (self, label)
  class (list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  integer:: n, ni, nb, nv, ne, nf
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('list_t%count_status', itimer=itimer)
  call self%lock%set ('count_status')
  nb=0; ni=0; nv=0; ne=0; nf=0; n = 0
  link => self%head
  do while (associated(link))
    n = n+1
    if (link%task%is_set(bits%frozen  )) nf = nf+1
    if (link%task%is_set(bits%boundary)) nb = nb+1
    if (link%task%is_set(bits%virtual )) nv = nv+1
    if (link%task%is_set(bits%internal)) ni = ni+1
    if (link%task%is_set(bits%external)) ne = ne+1
    link => link%next
  end do
  self%ni = ni
  self%nb = nb
  self%nv = nv
  self%n  = n
  if (self%verbose > 0) then
    if (present(label)) then
      write (io_unit%log,'(1x,a,6i6,2x,a)') &
      'list_t%count_status: n[ibveaf] = ', &
      self%ni, self%nb, self%nv, self%ne, self%na, nf, label
    else
      write (io_unit%log,'(1x,a,6i6,2x,a)') &
      'list_t%count_status: n[ibveaf] = ', &
      self%ni, self%nb, self%nv, self%ne, self%na, nf
    end if
  end if
  self%na = ni + nb - nf
  !-----------------------------------------------------------------------------
  call self%lock%unset ('count_status')
  call trace%end (itimer)
END SUBROUTINE count_status

!===============================================================================
!> Clear and then set the virtual and boundary bits
!===============================================================================
SUBROUTINE reset_status (self, check)
  class(list_t):: self
  logical, optional:: check
  logical:: flags(3,2)
  class(link_t), pointer:: link, nb
  integer:: n_internal, n_boundary, n_virtual, n_external, n_frozen
  !-----------------------------------------------------------------------------
  ! Clear virtual and boundary bits
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%reset_status')
  if (self%verbose > 1) &
    write (io_unit%output,*) 'reset_status'
  n_internal = 0
  n_boundary = 0
  n_virtual  = 0
  n_external = 0
  n_frozen   = 0
  call self%lock%set ('reset_status')
  !-----------------------------------------------------------------------------
  ! Set bit
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    if (link%task%rank == mpi%rank) then
      !-------------------------------------------------------------------------
      ! Assume internal -- this may change in set_status_and_flags()
      !-------------------------------------------------------------------------
      call link%task%set(bits%internal)
      call link%task%clear(bits%boundary+bits%virtual+bits%external)
    else
      !-------------------------------------------------------------------------
      ! Assume external -- this may change in set_status_and_flags()
      !-------------------------------------------------------------------------
      call link%task%set(bits%external)
      call link%task%clear(bits%internal+bits%boundary+bits%virtual)
    end if
    !---------------------------------------------------------------------------
    ! Run through nbors and set their status, optionally checking
    !---------------------------------------------------------------------------
    nb => link%nbor
    do while (associated(nb))
      if (present(check)) then
        flags(:,1) = [nb%needed, nb%needs_me, nb%download]
      end if
      call link%task%nbor_relations (nb%task, nb%needed, nb%needs_me, nb%download)
      if (present(check)) then
        flags(:,2) = [nb%needed, nb%needs_me, nb%download]
        if (.not.all(flags(:,2) .eqv. flags(:,1))) then
          write (stderr,*) link%task%id, flags(:,1)
          write (stderr,*)   nb%task%id, flags(:,2)
          flush (stderr)
          call io%abort ('nbor_relations() returned value inconsistent with initial values')
        end if
      end if
      nb => nb%next
    end do
    if (self%verbose > 1) &
      write(io_unit%log,'(a,i6,i4,3x,3l1)') &
        'reset_status: task, rank, IBV =', link%task%id, link%task%rank, &
        link%task%is_set(bits%internal), &
        link%task%is_set(bits%boundary), &
        link%task%is_set(bits%virtual)
    n_internal = n_internal + merge(1,0,link%task%is_set(bits%internal))
    n_boundary = n_boundary + merge(1,0,link%task%is_set(bits%boundary))
    n_virtual  = n_virtual  + merge(1,0,link%task%is_set(bits%virtual ))
    n_external = n_external + merge(1,0,link%task%is_set(bits%external))
    n_frozen   = n_frozen   + merge(1,0,link%task%is_set(bits%frozen  ))
    link => link%next
  end do
  self%nb = n_boundary
  self%nv = n_virtual
  self%ni = n_internal
  !-----------------------------------------------------------------------------
  ! Check that this does not interfere with end_time logics
  !-----------------------------------------------------------------------------
  self%na = n_internal+n_boundary-n_frozen
  self%n  = n_internal+n_boundary+n_virtual
  if (self%verbose > 1) then
    write (io_unit%log,'(5(a,i7,5x))') &
      'rank:',mpi%rank, &
      'n_internal:', n_internal, &
      'n_boundary:', n_boundary, &
      'n_virtual:' , n_virtual, &
      'n_frozen:'  , n_frozen
    flush (io_unit%log)
  end if
  !-----------------------------------------------------------------------------
  call self%lock%unset ('reset_status')
  call trace%end()
END SUBROUTINE reset_status

!===============================================================================
!> Check the nbor nbors for tasks ready to update, and the task itself
!===============================================================================
SUBROUTINE check_nbor_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%check_nbor_nbors')
  nbor => link%nbor
  do while (associated(nbor))
    call self%check_nbors (nbor%link)
    nbor => nbor%next
  end do
  call self%check_ready(link)
  call trace%end()
END SUBROUTINE check_nbor_nbors

!===============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready.  If the ready bit is already set it means the patch has already
!> been put in the ready queue, and should not be checked again.
!>
!> The link pointer and everything it points to are private to this task, and
!> are not at this point in time accessible from the ready queue.
!>
!> To prevent AMR from changing the nbor lists while they are being accessed,
!> it is necessary (and sufficient) to lock the current link
!===============================================================================
SUBROUTINE check_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor, nbors
  class(task_t), pointer:: task
  integer:: itr
  integer, save:: itimer=0
  !.............................................................................
  task => link%task                                     ! main task
  itr = 1
  call trace%begin('list_t%check_nbors', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Lock the other task link while accessing its nbor list, to prevent changes
  !-----------------------------------------------------------------------------
  nbor => link%nbor                                     ! first nbor
  do while (associated (nbor))                          ! keep going until end
    if (nbor%task%id==io%id_debug) print *, &
      'task', task%id,' needs task ', nbor%task%id, nbor%needs_me
    if (nbor%needs_me) then                             ! needs checking?
      if (self%verbose > 1) &
        write (io_unit%log,*) wallclock(), task%id, 'checks nbor', nbor%task%id
      call self%check_ready (nbor%link, lock=.true.)    ! pointer back
    else
      if (self%verbose > 1) &
        write (io_unit%log,*) wallclock(), task%id, 'no need to check', nbor%task%id
    end if
    nbor => nbor%next                                   ! next nbor
  end do
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we check if any bits are set indicating this to be the case.
  ! Note that this call to check_ready(), which concerns the task being updated,
  ! does NOT need locking of the task link.
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%virtual+bits%external)) then
    call self%check_ready (link)                      ! finally check link task
  end if
  call trace%end (itimer)
  return
END SUBROUTINE check_nbors

!===============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready, and send task data to non-locals.
!===============================================================================
SUBROUTINE check_all (self, repair)
  class(list_t):: self
  logical, optional:: repair
  !.............................................................................
  class(link_t), pointer:: link, nbor
  integer, save:: itimer=0
  integer:: nq
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%check_all', itimer=itimer)
  write (io_unit%log,*) wallclock(), 'check_all: phase 1', &
    self%nq, associated(self%queue)
  if (present(repair)) then
    nullify (self%queue)
  end if
  nq = self%nq
  link => self%head
  do while (associated (link))                          ! keep going until end
    if (link%task%is_clear (bits%virtual+bits%external)) then
      call self%check_ready (link)                      ! link ready?
      if (self%nq > nq) then
        write (io_unit%log,*) 'check_all found', link%task%id, link%task%time
        nbor => link%nbor
        do while (associated(nbor))
          write (io_unit%log,*) 'nbor, time =', nbor%task%id, nbor%task%time, &
            nbor%needed, nbor%needs_me, nbor%download
          nbor => nbor%next
        end do
        nq = self%nq
      end if
    end if
    link => link%next
  end do
  write (io_unit%log,*) wallclock(), 'check_all: done   ', &
    self%nq, associated(self%queue)
  if (.not.associated(self%queue)) then
    write (io_unit%log,*) 'check_all: phase 2, clearing ready bits'
    link => self%head
    do while (associated (link))                        ! keep going until end
      if (link%task%is_clear (bits%virtual+bits%external)) then
        if (link%task%is_set (bits%ready)) then
          write (io_unit%log,*) 'clearing ready bit on', link%task%id
          call link%task%clear (bits%ready)
        end if
        call self%check_ready (link)                    ! link ready?
      end if
      link => link%next
    end do
  end if
  call trace%end (itimer)
END SUBROUTINE check_all

!===============================================================================
!> Resend all boundary patches
!===============================================================================
SUBROUTINE resend_bdry (self)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  call trace%begin ('list_t%resend_bdry')
  call self%lock%set ('resend_bdry')
  link => self%head
  do while (associated(link))
    if (link%task%is_set (bits%boundary)) then
      call self%send_to_vnbors (link)
    end if
    link => link%next
  end do
  call self%lock%unset ('resend_bdry')
  call trace%end()
END SUBROUTINE resend_bdry

!===============================================================================
!> Find the oldest local and virtual patches
!===============================================================================
SUBROUTINE check_oldest (self)
  class(list_t):: self
  class(link_t), pointer:: link, oldest, oldestv, oldestnb, nbor
  real(8):: told, toldv
  integer, save:: phase=1
  !.............................................................................
  call trace%begin('list_t%check_oldest ')
  write (io_unit%log,*) 'check_oldest: phase 1', wallclock()
  flush (io_unit%log)
  told = 1d30
  toldv = told
  nullify (oldest)
  nullify (oldestv)
  nullify (oldestnb)
  link => self%head
  do while (associated (link))                        ! keep going until end
    if (link%task%is_clear (bits%frozen)) then
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
    end if
    link => link%next
  end do
  if (associated(oldest)) then
    write (io_unit%log,*) 'oldest active is ', oldest%task%id , told, &
      oldest%task%is_set(bits%ready), oldest%task%is_set(bits%busy) 
  else
    write (stderr,*) 'STALLED: oldest active not found '
  end if
  flush (io_unit%log)
  if (associated(oldestv)) then
    write (io_unit%log,*) 'the oldest virtual task is ', oldestv%task%id, toldv
    write (io_unit%log,*) 'it was last unpacked at', oldestv%task%unpack_time
    flush (io_unit%log)
  end if
  if (associated(oldest)) then
    nbor => oldest%nbor
    write (io_unit%log,2) oldest%task%id, oldest%task%time, oldest%task%rank, oldest%task%level, &
      oldest%needed, oldest%needs_me, oldest%task%is_set(bits%boundary), oldest%task%is_set(bits%virtual)
  2 format(i6,f12.6,2i4,5L3,f12.6)
    write (io_unit%log,*) 'nbor list:'
    do while (associated(nbor))
      write (io_unit%log,2) nbor%task%id, nbor%task%time, nbor%task%rank, nbor%task%level, &
        nbor%needed, nbor%needs_me, nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual), &
        nbor%task%is_ahead_of (oldest%task), nbor%task%unpack_time
      nbor => nbor%next
    end do
    call oldest%task%clear (bits%ready+bits%busy)
    call self%queue_by_time (oldest)                  ! add task to queue
    write (io_unit%log,*) 'task queued'
    flush (io_unit%log)
    write (io_unit%log,*) 'check_oldest: phase 1 done'
    flush (io_unit%log)
    link => oldest
  else
    call trace%end ()
    return
  end if
  !
  if (phase > 1) then
    if (associated(oldest)) then
      write (io_unit%log,*) 'check_oldest: phase 2', associated(oldest)
      write (io_unit%log,1) 'the oldest local task is', oldest%task%id, &
        'rank', oldest%task%rank, oldest%task%time, &
        oldest%task%is_set (bits%virtual), oldest%task%is_set (bits%external)
      told = 1d30
      link => oldest%nbor
      do while (associated(link))
        if (link%task%time < told) then
          told =link%task%time
          oldestnb => link
        end if
        write (io_unit%log,1) &
          'nbor', link%task%id, 'rank', link%task%rank, link%task%time, &
          link%task%is_set(bits%boundary), link%task%is_set (bits%virtual)
        1 format (1x,a,i9,3x,a,i6,1p,e16.6,2l3)
        link => link%next
      end do
    end if
    write (io_unit%log,'(f12.6,2(3x,a,i8,i6,f12.6,3x,4l1),3i5)') wallclock(), &
     'oldest:', oldest%task%id, oldest%task%rank, oldest%task%time, &
     oldest%task%is_set(bits%internal), oldest%task%is_set(bits%boundary), &
     oldest%task%is_set(bits%virtual), oldest%task%is_set(bits%external), &
     'nbor:', oldestnb%task%id, oldestnb%task%rank, oldestnb%task%time, &
     oldestnb%task%is_set(bits%internal), oldestnb%task%is_set(bits%boundary), &
     oldestnb%task%is_set(bits%virtual), oldestnb%task%is_set(bits%external), &
     self%nq, sent_list%n, recv_list%n, unpk_list%n
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
  end if
  call trace%end()
END SUBROUTINE check_oldest

!===============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!===============================================================================
SUBROUTINE check_ready (self, link, lock)
  class(list_t):: self
  class(link_t), pointer:: link
  logical, optional:: lock
  !.............................................................................
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: ignore, unit
  integer, save:: itimer=0
  logical:: ok
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we first check if any bits are set indicating this to be the case
  !-----------------------------------------------------------------------------
  task => link%task
  ignore = bits%ready + bits%busy + bits%remove + bits%virtual + bits%frozen
  if (task%is_set (ignore)) &
    return 
  if (self%detailed_timer) &
    call trace%begin ('list_t%check_ready', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Since only this thread can change the nbor list we are free to use it w/o
  ! lock in cases where task is the active task, but in other cases we need to
  ! lock the nbor list, so it isn't changed during the loop
  !-----------------------------------------------------------------------------
  if (present(lock) .and. omp_lock%links) &
    call link%lock%set ('check_ready')
  !-----------------------------------------------------------------------------
  ! Tasks that are frozen or removed should not be used in the loop
  !-----------------------------------------------------------------------------
  ignore = bits%frozen + bits%remove
  ok = .true.
  nbor => link%nbor                                     ! use original nbor list
  do while (associated (nbor))                          ! keep going until end
    !---------------------------------------------------------------------------
    ! Check nbors that are needed, not frozen, and have support
    !---------------------------------------------------------------------------
    if (nbor%needed .and. nbor%task%is_clear(ignore)) then
      if (.not. nbor%task%is_ahead_of(task)) then       ! and is ahead in time
        !$omp atomic update
        task%n_failed = task%n_failed+1
        if (self%verbose > 1 .or. task%n_failed > 10000) then
          if (self%verbose > 1) then
            unit = io_unit%log
          else
            unit = stderr
          end if
!          write (unit ,'(a,i4,1p,2(2x,a,i6,2x,a,g14.5),2(2x,a,g11.3))') &
!            'list_t%check_ready: on rank', task%rank, &
!            'task'               , task%id     , 'at t=', task%time, &
!            'failed on nbor task', nbor%task%id, 'at t=', nbor%task%time, &
!            'grace =', task%grace, ' dt=', nbor%task%dtime
        end if
        ok = .false.
        exit
      end if
    end if
    nbor => nbor%next                                   ! next nbor
  end do
  if (present(lock) .and. omp_lock%links) &
    call link%lock%unset ('check_ready')
  if (ok) then
    !$omp atomic write
    task%n_failed = 0
    call self%queue_by_time (link)                      ! add task to queue
  end if
  if (self%detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE check_ready

!===============================================================================
!> Check if the number of nbor tasks is consistent
!===============================================================================
SUBROUTINE consistency (self, link, i)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, i
  !.............................................................................
  if (self%debug < 2) return
  call trace%begin ('list_t%consistency',1)
  call self%lock%set ('list_t%consistency')
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
  call self%lock%unset ('list_t%consistency')
  call trace%end()
END SUBROUTINE consistency

!===============================================================================
!> Check if the number of queued tasks is consistent; i.e., of the number of
!> tasks linked at self%queue is the same as nq
!===============================================================================
SUBROUTINE qconsistency (self, ident)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, prev
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, ident
  real(8):: time
  !.............................................................................
  !if (self%debug < 1) return
  call trace%begin ('list_t%qconsistency',4)
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
  call trace%end()
END SUBROUTINE qconsistency

!===============================================================================
!> Check if the number of queued tasks is consistent; i.e., of the number of
!> tasks linked at self%queue is the same as nq
!===============================================================================
SUBROUTINE aconsistency (self, ident)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task, prev
  class(link_t), pointer:: nbor
  integer:: n, id1=0, id2=0, ident
  real(8):: time
  !.............................................................................
  !if (self%debug < 1) return
  call trace%begin ('list_t%aconsistency',4)
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
  call trace%end()
END SUBROUTINE aconsistency

!===============================================================================
!> Set non_leaf flag on non-leaf patches
!===============================================================================
SUBROUTINE init_nonleaf (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan
  class(link_t), pointer:: nbor
  integer:: n_add
  !.............................................................................
  return                                                ! FIXME: remove call?
  call trace%begin('list_t%init_nonleaf ')
  call self%lock%set ('list_t%init_nonleaf')
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
  call self%lock%unset ('list_t%init_nonleaf')
  call trace%end()
END SUBROUTINE init_nonleaf

!===============================================================================
!> Add a new task (e.g. originating from AMR or load balancing)
!===============================================================================
SUBROUTINE add_new_link (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%add_new_task')
  !-----------------------------------------------------------------------------
  ! Make sure class(patch_t) tasks have a pointer back to the task link
  !-----------------------------------------------------------------------------
  task => link%task
  select type (task)
  class is (patch_t)
    task%link => link
  end select
  !-----------------------------------------------------------------------------
  ! Add an nbor list, and set status bits and flags.  Add the task to the
  ! task_list, and if it turns out to be a boundary task, send it to rank nbors.
  !-----------------------------------------------------------------------------
  call self%init_nbors (link)
  call self%append (link)
  if (link%task%is_set (bits%boundary)) &
    call self%send_to_vnbors (link)
  !-----------------------------------------------------------------------------
  ! Increment task counters on the rank
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%virtual)) then
    !$omp atomic update
    timer%levelpop (task%level) = &          ! increment level population count
    timer%levelpop (task%level) + 1
  end if
  !-----------------------------------------------------------------------------
  ! Set the flag init_nbors in the nbor tasks, so the next time they are being
  ! updated, they obtain updated nbor lists
  !-----------------------------------------------------------------------------
  call self%set_init_nbors (link)
  !-----------------------------------------------------------------------------
  ! Check the ready state of the new task, which already knows its nbors.
  ! The nbors do not have the new task in their nbor lists yet, but will check
  ! correspondingly when the new task has been added.
  !-----------------------------------------------------------------------------
  call self%check_ready (link)
  call trace%end ()
END SUBROUTINE add_new_link

!===============================================================================
!> The very first initialization of a neighbor list.  After this, it can be
!> maintained without searching the entire task list on the rank.
!===============================================================================
SUBROUTINE init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor, scan, new_head, new_sort, old_head, old_sort, nbors
  integer:: n_add
  logical:: overlaps
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%init_nbors ', itimer=itimer)
  nullify (new_head, new_sort)
  n_add = 0
  call self%lock%set ('init_nbors')
  scan => self%head                                   ! scan all links (FIXME)
  do while (associated(scan))                         ! until end
    !---------------------------------------------------------------------------
    ! Add ranks that differ from the process rank to the list of load balance
    ! neighbors (each rank is only added once)
    !---------------------------------------------------------------------------
    !if (scan%task%rank /= mpi%rank .and. load_balance%on) then
      !call load_balance%add (scan%task%rank)
    !end if
    if (.not.associated(scan,link) .and. associated(scan%task)) then
      overlaps = scan%task%overlaps(link%task) &
           .and. scan%task%level <= link%task%level+1
      !-------------------------------------------------------------------------
      ! Potential nbors are other patches that overlap
      !-------------------------------------------------------------------------
      if (io%verbose>2 .or. link%task%id==io%id_debug) &
        write (io_unit%log, '(a,3(i6,2(3f7.4,1x)),2x,2l2)') &
        'init_nbors:', &
        link%task%id, link%task%position, link%task%size, &
        scan%task%id, scan%task%position, scan%task%size, &
        0, scan%task%distance (link%task), self%size, &
        overlaps, scan%task%is_set(bits%virtual)
      if (overlaps) then  ! if its task overlaps
        !-----------------------------------------------------------------------
        ! The ifs select to add nbors to a patch if either the patch is not a
        ! virtual patch, or the nbor is not a virtual patch, thus excluding
        ! adding virtual patches as nbors to virtual patches
        !-----------------------------------------------------------------------
        allocate (nbor)
        nbor%task => scan%task                      ! pointer to task
        nbor%link => scan                           ! pointer to task link
        call link%add_nbor_by_rank (new_head, nbor) ! add in rank order
        call link%task%nbor_relations (nbor%task, nbor%needed, &
          nbor%needs_me, nbor%download)             ! set flags and status bits
        n_add = n_add+1
      end if
    end if
    scan => scan%next                               ! continue search
  end do
  call self%lock%unset ('init_nbors')
  !-----------------------------------------------------------------------------
  if (self%verbose > 1 .or. link%task%id == io%id_debug) &
    write(io_unit%log,'(a,i6,i4,3x,3l1)') &
      'init_nbors: link, rank, IBV =', link%task%id, link%task%rank, &
      link%task%is_set(bits%internal), &
      link%task%is_set(bits%boundary), &
      link%task%is_set(bits%virtual)
  if (self%verbose > 2 .or. link%task%id == io%id_debug) &
    write (io_unit%log,*) &
    'added ', n_add, ' neighbors to', link%task%id, &
    link%task%is_set(bits%boundary), link%task%is_set(bits%virtual)
  !-----------------------------------------------------------------------------
  ! The links only needs to be locked the brief moment while the nbor pointers
  ! are changed.  
  !-----------------------------------------------------------------------------
  nullify (new_sort)
  call link%sort_nbors_by_level (new_head, new_sort)
  call link%lock%set ('init_nbors')
  old_head => link%nbor
  old_sort => link%nbors_by_level
  link%nbor => new_head                               ! switch nbor list
  link%nbors_by_level => new_sort
  call link%lock%unset ('init_nbors')
  !-----------------------------------------------------------------------------
  ! Remove the two previour nbor lists
  !-----------------------------------------------------------------------------
  call link%remove_nbor_list (old_head)               ! recover nbor list mem
  call link%remove_nbor_list (old_sort)
  !-----------------------------------------------------------------------------
  link%task%n_nbors = n_add                           ! remember n_nbors
  !-----------------------------------------------------------------------------
  ! Finally, clear the bits%init_nbors, in case it was set, and perform a
  ! check_nbors() on the link, in order to detect if an nbor that was waiting
  ! for an update on this task, which did not yet have it in its nbor list,
  ! now should be added to the ready queue.
  !-----------------------------------------------------------------------------
  call link%task%clear (bits%init_nbors)
  call trace%end (itimer)
END SUBROUTINE init_nbors

!===============================================================================
!> Set bits%init_nbors on nbor tasks, if the task itself is new (and hence has
!> a new nbor list, which should be matched by the nbors nbor list).  For safe
!> measure, we do this twice, even though one time should in principle be eno
!===============================================================================
SUBROUTINE set_init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%set_init_nbors')
  if (self%verbose > 0) &
    write (io_unit%log,*) &
      'set_init_nbors: id =', link%task%id
  nbor => link%nbor
  do while (associated(nbor))
    call nbor%task%set (bits%init_nbors)
    if (self%verbose > 1) then
        write (io_unit%log,*) 'bits%init_nbors set ', nbor%task%id, nbor%task%istep
    end if
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE set_init_nbors

!===============================================================================
!> Initialize the nbor lists of all nbors of link.  Note that the nbor list of
!> an nbor starts at nbor%link%nbor.   To prevent other threads from modifying
!> nbor list we need to loop over, we lock the link in the meantime, and to 
!> avoid a possible deadlock caused by having more than one lock active, we use
!> a copy of the nbor list (cf. check_nbors())  
!===============================================================================
SUBROUTINE init_nbor_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor, nbors
  !-----------------------------------------------------------------------------
  ! Re-generate the link nbor list, to make sure a correct nbor loop below
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%init_nbor_nbors')
  call link%lock%set ('init_nbnbs')           ! lock the link while modifying
  call self%init_nbors (link)                 ! initialize for loop below
  call self%check_nbors (link)
  call link%copy_nbor_list (link%nbor, nbors) ! copy to temporary nbor list
  call link%lock%unset ('init_nbnbs')         ! unlock the link
  !-----------------------------------------------------------------------------
  ! When resetting the nbors of nbors, notify rank nbors to do the same, by
  ! setting bits%init_nbors
  !-----------------------------------------------------------------------------
  nbor => nbors                               ! loop over temporary list
  do while (associated(nbor))
    if (nbor%link%task%is_set (bits%boundary)) &
      call nbor%link%task%set (bits%init_nbors)
    if (self%verbose > 0) &
      write (io_unit%log,*) 'init_nbor_nbors: id =', nbor%link%task%id
    call self%init_nbors (nbor%link)
    call self%check_nbors (nbor%link)
    nbor => nbor%next
  end do
  call link%remove_nbor_list (nbors)          ! remove the nbor list copy
  call trace%end()
END SUBROUTINE init_nbor_nbors

!===============================================================================
!> Check if the nbor list is still OK.  A first test is to check if any of the
!> nbors no longer overlaps.  One could / should also check if any nbor nbor
!> now overlaps.
!===============================================================================
SUBROUTINE refresh_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: nbor1, nbor2
  logical:: ok
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%refresh_nbors ', itimer=itimer)
  ok = .true.
  nbor1 => link%nbor
  do while (ok .and. associated(nbor1))
    ok = ok .and. nbor1%task%overlaps(link%task)
    !nbor2 => nbor1%link%nbor
    !do while (ok .and. associated(nbor2))
    !  ok = ok .and. .not. nbor2%task%overlaps(link%task)
    !  nbor2 => nbor2%next
    !end do
    nbor1 => nbor1%next
  end do
  if (.not.ok) then
    print *, link%task%id, 'list_t%refresh_nbors: init_nbprs'
    call self%init_nbors (link)
    call self%check_nbors (link)
  end if
  call trace%end (itimer)
END SUBROUTINE refresh_nbors

!===============================================================================
!> The very first initialization of the neighbor lists.  After this, it can be
!> maintained without a global search.  Distribute the work over available
!> threads.
!===============================================================================
SUBROUTINE init_all_nbors (self)
  class(list_t):: self
  class(link_t), pointer:: link
  integer:: i, i0, i1, i2, n=0, nv=0, nb=0
  real(8):: time
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin('list_t%init_all_nbors', itimer=itimer)
  call timer%tic (time)
  call init_nonleaf (self)                              ! set non_leaf bit
  call self%reset_status                                ! update task bits
  call self%count_status                                ! update task counts
  i2 = self%n/omp%nthreads + 1                          ! tasks per thread
  !$omp parallel private(i,i0,i1,link) shared(self,i2) default(none)
  i0 = omp%thread*i2
  i1 = i0 + i2
  link => self%head                                     ! start checking links
  i = 0                                                 ! task counter
  do while (associated(link))
    if (i >= i0 .and. i<i1) then                        ! thread interval
      call self%init_nbors (link)
      !write (io_unit%log,*) 'init_all_nbors:', wallclock(), i0, i, i1
    end if
    i = i+1                                             ! increment counter
    link => link%next                                   ! continue list
  end do
  !$omp end parallel
  !call self%reset_status                                ! set bits and count
  call timer%toc ('init_all_nbors', 1, time)
  call trace%end (itimer)
END SUBROUTINE init_all_nbors

!===============================================================================
!> Remove parent patches, which have higher level patches inside
!===============================================================================
SUBROUTINE remove_parents (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan, parent
  !.............................................................................
  call trace%begin('list_t%remove_parents')
  call self%lock%set ('remove_parent')
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
        io%nwrite = io%nwrite - 1
        io%ntask  = io%ntask - 1
        call parent%remove_from_nbors (parent)
        call self%remove_link (parent)
      end if
    end do
    link => link%next                                   ! continue list
  end do
  call self%lock%unset ('remove_parent')
  call trace%end()
END SUBROUTINE remove_parents

!===============================================================================
!===============================================================================
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

!===============================================================================
!===============================================================================
SUBROUTINE statistics (self)
  class(list_t):: self
  write (io_unit%log,'(a,3(5x,a,i6))') 'output_experiment:', &
    'n_ready:', self%n_ready, &
    'n_check:', self%n_check, &
    'n_nbor:' , self%n_nbor
END SUBROUTINE statistics

!===============================================================================
!===============================================================================
SUBROUTINE append_list (self, task_list)
  class(list_t):: self
  class(list_t), pointer:: task_list
  class(link_t), pointer:: link
  integer:: n
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%append_list')
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
  call trace%end()
END SUBROUTINE append_list

!===============================================================================
SUBROUTINE remove_task (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  class(link_t), pointer:: link, prev
  !.............................................................................
  call trace%begin('list_t%remove_task')
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' wait for tasklist(3)'
  call self%lock%set ('remove_task')
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(3)'
  link => self%head
  nullify (prev)
  do while (associated(link))
    if (link%task%id == task%id) then
      if (associated(prev)) then
        prev%next => link%next
      else
        self%head => link%next
      end if
      goto 9
    end if
    prev => link
    link => link%next
  end do
  print*,'ERROR: failed to remove from task list, task', task%id
9 continue
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  !-----------------------------------------------------------------------------
  call self%count_status()
  call self%lock%unset ('remove_task')
  if (self%verbose > 1) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(3)'
  call trace%end()
END SUBROUTINE remove_task

!===============================================================================
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
!===============================================================================
SUBROUTINE queue_by_time (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  logical:: error
  integer:: verbose
  !-----------------------------------------------------------------------------
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' waitfor tasklist(4)'
  call self%lock%set ('queue_by_time')
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(4)'
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
      call io%abort()
    end if
    io%verbose = verbose
    io%do_trace = .false.
  end if
  call self%lock%unset ('queue_by_time')
  if (self%verbose > 2) &
    write (io%output,*) wallclock(),' thread',omp%thread,' unlocked tasklist(4)'
END SUBROUTINE queue_by_time

SUBROUTINE real_queue_by_time (self, this, error)
  class(list_t):: self
  class(link_t), pointer:: this
  logical:: error
  class(link_t), pointer:: next, prev, link
  class(task_t), pointer:: task
  integer, save:: itimer=0
  integer:: nit, i, j
  !-----------------------------------------------------------------------------
  ! Refuse to queue virtual tasks and finished tasks
  !-----------------------------------------------------------------------------
  error=.false.
  task => this%task
  if (task%is_set (bits%virtual + bits%frozen)) then
    write(stderr,*) mpi%rank, omp%thread, &
      'ERROR: tried to queue virtual or frozen task, VF =', &
      task%is_set(bits%virtual), task%is_set(bits%frozen)
    write(io_unit%log,*) mpi%rank, omp%thread, &
      'ERROR: tried to queue virtual or frozen task, VF =', &
      task%is_set(bits%virtual), task%is_set(bits%frozen)
    return
  end if
  if (task%time > io%end_time) then
    write (io_unit%mpi,*) mpi%rank, omp%thread, &
      'ERROR: tried to queue a task beyond end_time, VF =', &
      task%id, task%time, io%end_time, &
      task%is_set(bits%virtual), task%is_set(bits%frozen)
    link => self%head
    i = 0
    j = 0
    write (io_unit%mpi,*) 'task list:'
    do while (associated(link))
      i = i+1
      task => link%task
      if (task%time > io%end_time) then
        j = j+1
        write (io_unit%mpi,*) i, j, task%id, task%time, &
          task%is_set(bits%ready), task%is_set(bits%busy), &
          task%is_set(bits%boundary), task%is_set(bits%virtual), &
          task%is_set(bits%frozen), associated(link%next_time)
      end if
      link => link%next
    end do
    j = 0
    link => self%queue
    write (io_unit%mpi,*) 'queue:'
    do while (associated(link))
      task => link%task
      j = j+1
      write (io_unit%mpi,*) j, task%id, task%time, &
        task%is_set(bits%ready), task%is_set(bits%busy), &
        task%is_set(bits%boundary), task%is_set(bits%virtual), &
        task%is_set(bits%frozen), associated(link%next_time)
      link => link%next_time
    end do
    write (io_unit%mpi,*) 'na =', self%na
    call self%reset_status()
    call self%count_status('overtime')
    flush (io_unit%mpi)
    !call mpi%delay (5e3)
    !call io%abort ('task overtime')
    return
  end if
  !-----------------------------------------------------------------------------
  if (self%detailed_timer) &
    call trace%begin ('list_t%queue_by_time ', itimer=itimer)
  !$omp atomic
  mpi_mesg%n_ready = mpi_mesg%n_ready+1
  !-----------------------------------------------------------------------------
  ! Lock the queue while inserting the task, and mark the task as being in the
  ! ready queue (or executing -- the bit is cleared only after updating).
  ! Double check that no other thread has set the ready bit on this task.
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%ready+bits%busy)) then
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
        if (self%detailed_timer) &
          call trace%end (itimer)
        return
      end if
      if (associated(next%task, task)) then
        write (io_unit%log,*) omp_mythread, ' WARNING: task', task%id, ' is already in ready queue'
        go to 9
      else if (next%task%time > task%time) then
        this%next_time => next
        if (associated(prev)) then
          if (io%verbose > 1) then
            write (io_unit%output,'(i4,2x,a,i6,a,i6,i5,1p,3g16.6)') &
              omp%thread, 'task',task%id,' in ready queue between',prev%task%id,next%task%id, &
              prev%task%time, this%task%time, next%task%time
            write (io_unit%log,'(i4,2x,a,i6,a,i6,i5,1p,3g16.6)') &
              omp%thread, 'task',task%id,' in ready queue between',prev%task%id,next%task%id, &
              prev%task%time, this%task%time, next%task%time
          end if
          !call prev%qlock%set ('queue_by_time')
          prev%next_time => this
          !call prev%qlock%unset ('queue_by_time')
        else
          !call self%queue%qlock%set ('queue_by_time')
          self%queue => this
          !call self%queue%qlock%unset ('queue_by_time')
          if (io%verbose > 1) then
            write (io_unit%output,*) 'task',task%id,' at ready queue head',self%nq
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
        write (io_unit%output,*) 'task',task%id,' in ready queue after',prev%task%id,self%nq
        write (io_unit%log,*) 'task',task%id,' in ready queue after',prev%task%id,self%nq
      end if
      prev%next_time => this
    else
      if (io%verbose > 1) then
        write (io_unit%output,*) 'task',task%id,' at ready queue head',self%nq
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
  if (self%detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE real_queue_by_time

!===============================================================================
!===============================================================================
SUBROUTINE queue_active (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  class(task_t), pointer:: task
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%detailed_timer) &
    call trace%begin ('list_t%queue_active ', itimer=itimer)
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
  if (self%detailed_timer) &
    call trace%end (itimer)
END SUBROUTINE queue_active

!===============================================================================
!> Remove a task link from the active queue
!===============================================================================
SUBROUTINE remove_active (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !-----------------------------------------------------------------------------
  if (.not.associated(self%active)) return
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
  call trace%end()
END SUBROUTINE remove_active

!===============================================================================
!> Add a link node into the task list, in increasing quality order
!===============================================================================
SUBROUTINE add_by_quality (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !.............................................................................
  call trace%begin ('list_t%add_by_quality')
  call self%lock%set ('add_by_quality')
  nullify (prev)
  next => self%head
  do while (associated(next))
    if (associated(next%task,this%task)) then
      write (io_unit%log,*)'ERROR: Trying to add task', this%task%id, ' which is already there'
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
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  !-----------------------------------------------------------------------------
  call self%count_status()
  call self%lock%unset ('add_by_quality')
  call trace%end()
END SUBROUTINE add_by_quality

!===============================================================================
!> Print a table with patch id, positions, and all neighbors
!===============================================================================
SUBROUTINE check_queue (self)
  class(list_t):: self
  class(link_t), pointer:: link
  real:: previous
  integer:: nq
  !.............................................................................
  call self%lock%set ('check_queue')
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
  call self%lock%unset ('check_queue')
END SUBROUTINE check_queue

!===============================================================================
!> Insert 'new' task before 'old' task
!===============================================================================
SUBROUTINE insert (self, old, new)
  class(list_t):: self
  class(link_t), pointer:: old, new
  !.............................................................................
  call trace%begin ('list_t%insert '//self%name)
  call self%lock%set ('insert')
  if (associated(old%prev)) then
    new%prev => old%prev
    new%prev%next => new
  else
    self%head => new
  end if
  new%next => old
  new%next%prev => new
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  !-----------------------------------------------------------------------------
  call self%count_status()
  call self%lock%unset ('insert')
  call trace%end()
END SUBROUTINE insert

!===============================================================================
!> Print a table with task id, positions, and all neighbors
!===============================================================================
SUBROUTINE print_queue (self)
  class (list_t):: self
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call self%print_queue_times ('print_queue')
END SUBROUTINE print_queue

!===============================================================================
!> Print the queue until and including a given task
!===============================================================================
SUBROUTINE print_queue_until (self, task)
  class(list_t):: self
  class(task_t):: task
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call self%lock%set ('print_queue_until')
  link => self%queue
  if (io%master) write (io_unit%log,*)'print_queue_until:'
  do while (associated(link))
    write (io_unit%log,'("queue: id, t0, t, t0-t =",i7,1p,3e12.3)') &
      link%task%id, task%time, link%task%time, task%time-link%task%time
    if (link%task%id==task%id) exit
    link => link%next_time
  end do
  call self%lock%unset ('print_queue_until')
END SUBROUTINE print_queue_until

!===============================================================================
!> Print the queue
!===============================================================================
SUBROUTINE print_queue_times (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call self%lock%set ('print_queue_times')
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
  call self%lock%unset ('print_queue_times')
  write (io_unit%log,*)'done'
END SUBROUTINE print_queue_times

!===============================================================================
!===============================================================================
SUBROUTINE print_list (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  class(link_t), pointer:: link, prev
  class(task_t), pointer:: task
  character(len=32):: type
  !.............................................................................
  call trace%begin ('list_t%print')
  if (io%verbose>0) &
    write (io_unit%log,*) 'list_t%print: ',trim(self%name), self%n
  call self%lock%set ('print_list')
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
  call self%lock%unset ('print_list')
  call trace%end()
END SUBROUTINE print_list

!===============================================================================
!> Print a table with task id, positions, and all neighbors
!===============================================================================
SUBROUTINE print_tasks (self, label)
  class (list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  !if (io%verbose < 2) return
  call trace%begin ('list_t%print_tasks')
  if (present(label)) then
    write (io_unit%log,'(a,i8,":")') 'task_list: '//trim(label), self%n
  else
    write (io_unit%log,'(a,i8,":")') 'task_list', self%n
  end if
  link => self%head
  do while (associated(link))
    if (io%verbose < 3) then
      write (io_unit%log,'(i8,g15.6,3x,3l1)') link%task%id, link%task%time, &
        link%task%is_set(bits%ready), link%task%is_set(bits%boundary), &
        link%task%is_set(bits%virtual)
    else
      call link%print_nbors
    end if
    link => link%next
  end do
  write (io_unit%log,*) ''
  call trace%end()
END SUBROUTINE print_tasks

!===============================================================================
!> Initialize boundary bits, which may (optionally) be used to apply boundary
!> conditions in experiment_t%update() procedures. 
!===============================================================================
SUBROUTINE init_bdries (self)
  implicit none
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: patch
  real(8):: position(3), limit(3)
  !-----------------------------------------------------------------------------
  call trace%begin('task_list_t%init_bdries')
  link => self%head
  do while (associated(link))
    patch => link%task
    select type (patch)
    class is (patch_t)
      call patch%init_bdries
    end select
    self%lc  = min(self%lc,  link%task%position)
    self%uc  = max(self%uc,  link%task%position)
    self%llc = min(self%llc, link%task%position-0.5000000001_8*link%task%size)
    self%urc = max(self%urc, link%task%position+0.5000000001_8*link%task%size)
    link => link%next
  end do
  call trace%end()
contains
  !-----------------------------------------------------------------------------
  function formatted (task) result (out)
  class(patch_t):: task
  character(len=6):: out
  out = '......'
  if (task%boundaries%is_set(bits%xl))  out(1:1) = 'T'
  if (task%boundaries%is_set(bits%xu))  out(2:2) = 'T'
  if (task%boundaries%is_set(bits%yl))  out(3:3) = 'T'
  if (task%boundaries%is_set(bits%yu))  out(4:4) = 'T'
  if (task%boundaries%is_set(bits%zl))  out(5:5) = 'T'
  if (task%boundaries%is_set(bits%zu))  out(6:6) = 'T'
  end function formatted
END SUBROUTINE init_bdries

!===============================================================================
!> Check nbor list of self, verifying that %task and %link are consistent
!===============================================================================
SUBROUTINE check_nbor_list (link, label, verbose)
  class(link_t):: link
  character(len=*):: label
  integer, optional:: verbose
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  !.............................................................................
  nbor => link%nbor
  do while (associated(nbor))
    if (associated(nbor%link)) then
      if (associated(nbor%task)) then
        task => nbor%task
        select type (task)
        class is (patch_t)
          if (.not.associated(task%link,nbor%link)) then
            print *, label, link%task%id, task%id, &
              'ERROR: nbor%link, task%link not associated'
          end if
        end select
      else
        print *, label, link%task%id, 'ERROR: nbor%task not associated'
      end if
    else
      print *, label, link%task%id, 'ERROR: nbor%link not associated'
    end if
    if (present(verbose)) then
      task => nbor%link%task
      select type (task)
      class is (patch_t)
      print *, link%task%id, nbor%task%id, task%id, associated(task%link,nbor%link)
      end select
    end if
    nbor => nbor%next
  end do
END SUBROUTINE check_nbor_list

!===============================================================================
!> Give a task to another rank, making sure to update the nbor relations on
!> this rank -- doing the same on the receiver rank is handled when unpacking,
!===============================================================================
SUBROUTINE give_to (self, link, rank)
  class(list_t):: self
  class(link_t), pointer:: link
  integer:: rank
  !-----------------------------------------------------------------------
  ! Set the bits and rank as they should be on the local rank, and mark
  ! the task with a swap_request bit
  !-----------------------------------------------------------------------
  call link%task%clear (bits%boundary)
  call link%task%set (bits%virtual+bits%swap_request)
  link%task%rank = rank
  !-----------------------------------------------------------------------
  ! With the new rank and status settings, update the status of the nbors
  ! and their nbor lists, which may all be affected by the change of rank.
  ! The ones that are new boundary tasks will be sent to nbor ranks by
  ! the update_nbor_status procedure.  The rank itself will be sent by
  ! task_list_mod::update, and should not be sent from here.
  !-----------------------------------------------------------------------
  call self%update_nbor_status (link)
  if (io%verbose>1) &
    call self%test_nbor_status (link)
END SUBROUTINE give_to

!===============================================================================
!> Update the nbor/vnbor status of a links nbor list, including the self link.
!> If the boundary bit gets set, send the task to its vnbors.  Before doing so,
!> the nbor list of the nbors need to be updated, to reflect the new task rank
!===============================================================================
SUBROUTINE update_nbor_status (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, nbor, this
  integer:: status
  !-----------------------------------------------------------------------------
  call trace%begin ('link_mod::update_nbor_status')
  call link%lock%set ('update_nbor_status')
  !
  nbor => link%nbor                                     ! then help others
  do while (associated(nbor))
    if (io%verbose>1) &
      write (io_unit%log,*) 'update_nbor_status:', nbor%task%id, &
      nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual)
    status = nbor%task%status                           ! current status
    call nbor%link%remove_nbor (link)                   ! remove out-of-order
    allocate (this)                                     ! new nbor link
    this%link => link                                   ! attach task link
    this%task => link%task                              ! attach task
    call nbor%link%add_nbor_by_rank (link%nbor, this)   ! add in rank order
    call self%update_status (nbor%link)                       ! new status
    if (nbor%task%is_set (bits%boundary).and.iand(status,bits%boundary)==0) then
      if (io%verbose>0) &
        write (io_unit%log,*) 'sending new bdry task',nbor%task%id,' to vnbors'
      if (.not.nbor%task%is_clear (bits%virtual)) then
        call nbor%task%clear (bits%virtual)
        write (io_unit%log,*) 'WARNING: needed to clear virtual bit(1)'
      end if
      call nbor%task%set (bits%swap_request)
      nbor%task%rank = mpi%rank                         ! since it's a boundary
      call self%send_to_vnbors (nbor%link)
      call nbor%link%task%clear (bits%swap_request + &  ! sometimes set -- make
                                 bits%init_nbors)       ! sure clear after send
    end if
    nbor => nbor%next                                   ! check next nbor
  end do
  !
  if (io%verbose>0) &
    write (io_unit%log,*) 'update_nbor_status: id =', link%task%id, associated(link)
  status = link%task%status                             ! current status
  call self%update_status (link)                        ! new status
  if (link%task%is_set (bits%boundary).and.iand(status,bits%boundary)==0) then
    if (io%verbose>0) &
      write (io_unit%log,*) 'sending new bdry task',link%task%id,' to vnbors'
    if (.not.nbor%task%is_clear (bits%virtual)) then
      write (io_unit%log,*) 'WARNING: needed to clear virtual bit(2)'
      call nbor%task%clear (bits%virtual)
    end if
    call nbor%task%set (bits%swap_request)
    nbor%task%rank = mpi%rank                           ! since it's a boundary
    call self%send_to_vnbors (link)                     ! need to send ...
    call nbor%task%clear (bits%swap_request)
  end if
  call link%lock%unset ('update_nbor_status')
  call trace%end()
END SUBROUTINE update_nbor_status

!===============================================================================
!> Update the nbor/vnbor status of a link.  If self belongs to the local rank
!> it is either an internal or boundary rank, if not it is virtual or external,
!> depending on what nbors it has.
!===============================================================================
SUBROUTINE update_status (self, link)
  class(list_t):: self
  class(link_t):: link
  class(link_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call link%lock%set ('update_status')
  nbor => link%nbor
  if      (link%task%rank == mpi%rank) then
    do while (associated(nbor))
      if (nbor%task%rank /= mpi%rank)  then
        call link%task%set   (bits%boundary)
        call link%task%clear (bits%internal+bits%external+bits%virtual)
        if (io%verbose>0) &
          write (io_unit%log,*) ' set to boundary id =', link%task%id
        return
      end if
      nbor => nbor%next
    end do
    call link%task%set   (bits%internal)
    call link%task%clear (bits%boundary+bits%external+bits%virtual)
    if (io%verbose>0) &
      write (io_unit%log,*) ' set to internal id =', link%task%id, associated(link%nbor)
  else
    do while (associated(nbor))
      if (nbor%task%rank == mpi%rank)  then
        call link%task%set   (bits%virtual)
        call link%task%clear (bits%internal+bits%external+bits%boundary)
        if (io%verbose>0) &
          write (io_unit%log,*) ' set to virtual  id =', link%task%id
        return
      end if
      nbor => nbor%next
    end do
    call link%task%set   (bits%external)
    call link%task%clear (bits%internal+bits%boundary+bits%virtual)
   !call link%remove_from_nbors
    if (io%verbose>0) &
      write (io_unit%log,*) '    set to external id =', link%task%id
  end if
  call link%lock%unset ('update_status')
END SUBROUTINE update_status

!===============================================================================
!> Send a package with updated time slice info to neighbor ranks.
!>
!> 1) assemble task data into a message
!> 3) send the message to the nbors rank (only once per rank)
!> 4) add the message to a sent_list
!> 5) check the sent list for completed sends
!> 6) remove the message from the sent_list and disassemble the message
!===============================================================================
SUBROUTINE send_to_vnbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  class(mesg_t), pointer:: mesg
  character(len=24):: label
  integer:: ierr, rank, tag, seq
  integer, save:: itimer=0
  !.............................................................................
  task => link%task
  if (task%is_clear (bits%boundary)) then
    write (stdout,*) mpi%rank, omp%thread, &
      'ERROR: trying to send non-boundary task', task%id
    return
  end if
  if (task%is_set (bits%virtual)) then
    write (stdout,*) mpi%rank, omp%thread, &
      'ERROR: trying to send virtual task', task%id
    return
  end if
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%send_to_vnbors', itimer=itimer)
  if (task%id == io%id_debug) &
    write(io%output,*) 'DBG link_t%send_to_vnbors: id, rank =', &
      task%id, mpi%rank
  !-----------------------------------------------------------------------------
  ! Pack task into mesg
  !-----------------------------------------------------------------------------
  select type (task)
  class is (patch_t)
  call task%pack (mesg)
  end select
  if (mpi_mesg%uniq_mesg .and. mpi_mesg%tag_type == 1) then
    task%seq = task%seq + 1
    tag = mod(task%seq,100) + 100*task%id
  else
    tag = task%id
  end if
  !-----------------------------------------------------------------------------
  ! Same id for all steps, and increment mesg%seq only once, for all ranks
  ! If the task is new, force the sequence number to be 1
  !-----------------------------------------------------------------------------
  mesg%id = task%id
  !-----------------------------------------------------------------------------
  ! Send to all unique ranks in the nbor list (which is sorted by rank).
  !-----------------------------------------------------------------------------
  rank = -1                                                     ! previous rank
  nbor => link%nbor                                             ! first nbor
  do while (associated(nbor))                                   ! until end
    nbtask => nbor%task
    if (nbtask%rank/=mpi%rank .and. nbtask%rank/=rank) then     ! new rank?
      rank = nbtask%rank                                        ! current rank
      if (mpi_mesg%uniq_mesg .and. mpi_mesg%tag_type == 2) then
        !$omp atomic capture
        sequence(rank) = sequence(rank)+1
        seq = sequence(rank)
        !$omp end atomic
        tag = mod(seq,100) + 100*task%id
      end if
      if (task%logging > 1) then
        write (label,'(a,i4,i8)') 'vnbor  ', rank, tag
        call task%log (label)
      end if
      call mesg%send (rank, tag=tag)                            ! send it
      !write (stdout,*) 'send: ', task%id, mesg%id, tag, mesg%tag, rank
      if (self%verbose > 0 .or. task%id == io%id_debug) then
        write (io_unit%mpi,'(f12.6,2x,a,i9,1p,e18.6,i9,2x,a,i5,2x,5l1)') &
          wallclock(), 'send_to_vnbors: sent', &
          mesg%id, task%time, tag, 'to', rank, &
          task%is_set (bits%internal), &
          task%is_set (bits%boundary), &
          task%is_set (bits%virtual), &
          task%is_set (bits%external), &
          task%is_set (bits%swap_request)
        flush (io_unit%mpi)
      end if
    end if
    nbor => nbor%next                                           ! next nbor
  end do
  !-----------------------------------------------------------------------------
  ! Add the message to mesg%sent_list, if it was sent, and clear the swap bit
  !-----------------------------------------------------------------------------
  call mpi_mesg%sent (mesg)                                     ! add to sent_list
  call task%clear (bits%swap_request)
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE send_to_vnbors

!===============================================================================
!> Check the status of a link and all nbor links
!===============================================================================
SUBROUTINE test_nbor_status (self, link)
  class(list_t):: self
  class(link_t):: link
  class(link_t), pointer:: nbor
  !.............................................................................
  write (io_unit%log,*) 'test_nbor_status: link', link%task%id
  call test_status (link)
  nbor => link%nbor
  do while (associated(nbor))
    call test_status (nbor%link)
    nbor => nbor%next
  end do
END SUBROUTINE test_nbor_status

!===============================================================================
!> Check the status of a link task  with respect to its nbors
!===============================================================================
SUBROUTINE test_status (self)
  class(link_t):: self
  class(link_t), pointer:: nbor
  integer:: status, test, n
  !.............................................................................
  write (io_unit%log,*) 'test_status: link', self%task%id
  if (self%task%rank == mpi%rank) then
    status = bits%internal
  else
    status = bits%external
  end if
  nbor => self%nbor
  do while (associated(nbor))
    if      (nbor%task%rank /= mpi%rank .and. status==bits%internal) then
      status = bits%boundary
    else if (nbor%task%rank == mpi%rank .and. status==bits%external) then
      status = bits%virtual
    end if
    nbor => nbor%next
  end do
  test = iand(self%task%status,bits%internal+bits%boundary+bits%virtual+bits%external)
  if (status /= test) then
    write (io_unit%log,'(i9,2x,a,i6,2(2x,4l1))') self%task%id, 'inconsistent:', &
      self%task%rank, &
      self%task%is_set(bits%internal), &
      self%task%is_set(bits%boundary), &
      self%task%is_set(bits%virtual) , &
      self%task%is_set(bits%external), &
      iand(status,bits%internal)/=0, &
      iand(status,bits%boundary)/=0, &
      iand(status,bits%virtual )/=0, &
      iand(status,bits%external)/=0
    n = 0
    nbor => self%nbor
    do while (associated(nbor))
      n = n+1
      write (io_unit%log,'(2i9,2x,a,2i6)') self%task%id, nbor%task%id, &
        'inconsistent ranks', self%task%rank, nbor%task%rank
      nbor => nbor%next
    end do
    if (n==0) write (io_unit%log,*) self%task%id, ' inconsistent bits: has no nbors'
  end if
END SUBROUTINE test_status

!===============================================================================
!> Dump nbor relations, flags and status bits for the list
!===============================================================================
SUBROUTINE info (self)
  class(list_t):: self
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task, nbtask
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%info')
  link => self%head
  do while (associated(link))
    call link%info
    link => link%next
  end do
  call trace%end()
END SUBROUTINE info

END MODULE list_mod
