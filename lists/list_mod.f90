!*******************************************************************************
!> Module with list handling for generic class task_t objects.
!>
!> Currently, this module depends on patch_t in the following procedures.
!> If that dependence could be removed, list_mod could lie btw extras_t and
!> patch_t, with no need for gpatch_t.
!>
!> PROCEDURE            USED IN                 POTENTIALLY MOVE TO
!> -----------------------------------------------------------------------------
!> append_link          sink_particles          not needed -- use local code
!> prepend_link         short_characteristics   not needed -- use local code
!> add_new_link         refine_mod              not needed -- use local code
!> print_list                                   task_list_mod
!> init_bdries                                  task_list_mod
!> check_nbor_lists                             task_list_mod
!> send_to_vnbors       list_mod                problem: called from remove_and_reset
!>
!> The problem is that refine_t is used by task_list_t, and refine_t needs to
!> be able to deal with removing tasks, and sending messages about that to nbors
!>
!> Important call hierarchies:
!>   add_new_link
!>     init_nbors                    [ accesses task_list ]
!>       link%task%add_nbor_by_rank  [ accesses own nbor list ]
!>       link%task%nbor_relations    [ accesses own nbor list ]
!>       link%sort_nbors_by_level    [ accesses own nbor list ]
!>     append                        [ accesses task_list ]
!>       link%init                   [ redundant, benign ]
!>       update_counts(+1)           [ accesses own nbor list ]
!>     send_to_vnbors                [ accesses own nbor list ]
!>       task%pack
!>       mesg%send
!>     set_init_nbors                [ accesses own nbor list ]
!>       nbor%task%set
!>     check_ready                   [ accesses own nbor list ]
!>       queue_by_time               [ accesses task_list ]
!>
!> update_status calls:
!> task_list%update
!>   ...
!>   task%update
!>     refine%check_current
!>       make_child_patch
!>         add_new_link
!>           init_nbors
!>             update_status
!>           update_status
!>   update_status
!>   check_mpi
!>     unpack
!>       add_new_link
!>           init_nbord
!>             update_status
!>           update_status
!*******************************************************************************
MODULE list_mod
  USE io_mod
  USE omp_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE mpi_window_mod
  USE trace_mod
  USE bits_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  !USE load_balance_mod
  USE timer_mod
  USE running_mean_mod
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
    integer:: no=0                      ! has output
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
    logical:: face_nbors=.false.
    logical:: detailed_timer=.false.
    character(len=4)              :: kind='list'
    type(mpi_window_t)            :: mpi_window
    type(running_mean_t), pointer :: running_mean
  contains
    procedure:: init
    !procedure:: init_bdries
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
    procedure:: set_io_counts
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
    procedure:: prune_nbors
    procedure:: set_init_nbors
    procedure:: set_update_status
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
    procedure:: remove_from_queue
    procedure:: print_queue
    procedure:: print_queue_until
    procedure:: print_queue_times
    procedure:: print => print_list
    procedure:: print_tasks
    procedure, nopass:: check_nbor_list
    procedure:: update_nbor_status
    procedure:: update_status
    procedure:: swap
    procedure:: send_to_vnbors
    procedure:: test_nbor_status
    procedure:: info
    procedure:: check_support
    procedure:: find_nbors
    procedure:: check_nbor_consistency
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
  class(list_t), target:: self
  character(len=*), optional:: name
  !.............................................................................
  call trace%begin ('list_t%init')
  if (present(name)) then
    self%name = name
    write (io_unit%log,*)'init_list: ', name
    call self%lock%init (name(1:4))
  else
    write (io_unit%log,*)'init_list: list_t'
    call self%lock%init ('list')
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
  class(nbor_t), pointer, optional:: nbor
  logical, optional:: needed, needs_me
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  !.............................................................................
  call trace%begin ('list_t%append_task '//self%name)
  allocate (link)
  call link%init
  link%task => task
  !-----------------------------------------------------------------------------
  ! Note that adding a task to more than one list requires care; the task%link
  ! must remain a pointer to the link that has the nbor list, and the nbor list
  ! head must be accesses as task%link%nbor, not as link%nbor, for links in the
  ! 2nd task list.
  !-----------------------------------------------------------------------------
  if (present(nbor))     link%nbor     => nbor
  if (present(needed))   nbor%needed   = needed
  if (present(needs_me)) nbor%needs_me = needs_me
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
  if (self%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' wait for tasklist(2)'
  call self%lock%set ('list_t%append_link')
  if (self%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(2)'
  !-----------------------------------------------------------------------------
  ! Check that task%link is set on patch_t class tasks
  !-----------------------------------------------------------------------------
  if (associated(link%task)) then
    task => link%task
    select type (task)
    class is (patch_t)
    if (.not.associated(task%link)) then
      task%link => link
    end if
    end select
  else
    call io%abort('append_link: no task associated')
  end if
  !-----------------------------------------------------------------------------
  ! Append in such a way that the task list does not need to be locked
  !-----------------------------------------------------------------------------
!  call link%qlock%set ('append_link1')
  nullify(link%next)
  if (associated(self%tail)) then
    link%prev => self%tail
    link%prev%next => link
!    call self%tail%qlock%set ('append_link2')
    self%tail => link
!    call self%tail%qlock%unset ('append_link2')
  else
    self%head => link
    self%tail => link
  end if
!  call link%qlock%unset ('append_link1')
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  ! Calling update_counts is robust, but may take significant time.
  ! FIXME: Maintain counts incrementally (testing & validating consistency)
  !-----------------------------------------------------------------------------
  if (self%verbose > 2) &
    write (io_unit%log,*) 'list_t%append: (1) n =', &
      trim(self%name), self%n, link%task%id
  call self%update_counts (link, +1)
  if (self%verbose > 2) &
    write (io_unit%log,*) 'list_t%append: (2) n =', &
      trim(self%name), self%n, link%task%id
  call self%lock%unset ('list_t%append_link')
  if (self%verbose > 3) &
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
    write (io_unit%log,'(f12.6,1x,a,i6,2l3)') &
      wallclock(), 'list_t%remove_and_reset: id =', link%task%id, &
      link%task%is_set(bits%boundary), link%task%is_set(bits%virtual)
    write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
      wallclock(), omp%thread, 'remove', link%task%id, link%task%time
  end if
  !-----------------------------------------------------------------------------
  call self%lock%set ('remove_and_reset')  ! lock task list while changing
  call self%remove_link (link)             ! remove link from task list
  call self%lock%unset ('remove_and_reset')! release the task list
  call self%check_nbors (link)             ! check nbors, now link is removed
  call self%set_init_nbors (link)          ! trigger new nbor lists on nbors
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
  call trace%begin ('list_t%update_counts')
  associate (self_nv=>self%nv, self_nb=>self%nb, self_ni=>self%ni, &
             self_na=>self%na, self_no=>self%no, self_n=>self%n)
  if      (link%task%is_set (bits%virtual)) then
    !$omp atomic update
    self_nv = self_nv + delta
  else if (link%task%is_set (bits%boundary)) then
    !$omp atomic update
    self_nb = self_nb + delta
  else if (link%task%is_set (bits%internal)) then
    !$omp atomic update
    self_ni = self_ni + delta
  end if
  if (link%task%is_set (bits%boundary) .or. &
      link%task%is_set (bits%internal)) then
    !$omp atomic update
    self_na = self_na + delta
  end if
  if (link%task%is_set (bits%has_output)) then
    !$omp atomic update
    self_no = self_no + delta
  end if
  n = self%na + self%nv
  !$omp atomic write
  self_n = n
  end associate
  call self%set_io_counts
  !-----------------------------------------------------------------------------
  ! Validate the incremental method when task_list_t%verbose > 0
  !-----------------------------------------------------------------------------
  if (self%verbose > 0) then
    n  = self%n
    na = self%na
    call self%count_status ()
    if (n /= self%n .or. na /= self%na) then
      write (stdout,'(a,2(2i6,2x),a,i6,i3,2x,3l1)') &
        'list_t%update_counts WARNING: inconsistent counts', &
        n, self%n, na, self%na, &
        'after updating task, delta, IBV =', link%task%id, delta, &
        link%task%is_set(bits%internal), &
        link%task%is_set(bits%boundary), &
        link%task%is_set(bits%virtual)
      write (io_unit%mpi,'(a,2(2i6,2x),a,i6,i3,2x,3l1)') &
        'list_t%update_counts WARNING: inconsistent counts', &
        n, self%n, na, self%na, &
        'after updating task, delta, IBV =', link%task%id, delta, &
        link%task%is_set(bits%internal), &
        link%task%is_set(bits%boundary), &
        link%task%is_set(bits%virtual)
    end if
  end if
  call trace%end ()
END SUBROUTINE update_counts 

SUBROUTINE set_io_counts (self)
  class (list_t):: self
  !-----------------------------------------------------------------------------
  ! Make sure task counts for I/O purposes are up-to-date
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%set_io_counts')
  associate (io_ntask=>io%ntask, io_nwrite=>io%nwrite, io_ntotal=>io%ntotal)
  !$omp atomic write
  io_ntask = self%na
  !$omp atomic write
  io_nwrite = self%no
  !$omp atomic write
  io_ntotal = self%n_tasks
  end associate
  call trace%end ()
END SUBROUTINE set_io_counts

!===============================================================================
!> Count the number of internal, boundary, virtual, frozen, and external tasks
!===============================================================================
SUBROUTINE count_status (self, label)
  class (list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link
  class(nbor_t) , pointer:: nbor
  class(task_t), pointer:: task
  integer:: n, ni, nb, nv, ne, nf, no, nq
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('list_t%count_status', itimer=itimer)
  call self%lock%set ('count_status', itimer)
  nb=0; ni=0; nv=0; ne=0; nf=0; n = 0; no = 0; nq = 0
  link => self%head
  do while (associated(link))
    task => link%task
    n = n+1
    if (task%is_set(bits%frozen  )) nf = nf+1
    if (task%is_set(bits%boundary)) nb = nb+1
    if (task%is_set(bits%virtual )) nv = nv+1
    if (task%is_set(bits%internal)) ni = ni+1
    if (task%is_set(bits%external)) ne = ne+1
    if (task%is_set(bits%ready   )) nq = nq+1
    if (task%is_set(bits%has_output)) no = no+1
    if (task%is_clear (bits%boundary + &
                       bits%virtual  + &
                       bits%internal + &
                       bits%external)) then
      write (stderr,'(a,i6,3x,"ibvef =",1x,5l1)') &
        'list_t%count_status: WARNING, status bit is missing in id, IBVEF =', &
        task%id, &
        task%is_set(bits%internal), &
        task%is_set(bits%boundary), &
        task%is_set(bits%virtual ), &
        task%is_set(bits%external), &
        task%is_set(bits%frozen  )
      write (io_unit%log,'(a,i6,3x,"ibvef =",1x,5l1)') &
        'list_t%count_status: WARNING, status bit is missing in id, IBVEF =', &
        task%id, &
        task%is_set(bits%internal), &
        task%is_set(bits%boundary), &
        task%is_set(bits%virtual ), &
        task%is_set(bits%external), &
        task%is_set(bits%frozen  )
    end if
    link => link%next
  end do
  self%ni = ni
  self%nb = nb
  self%nv = nv
  self%no = no
  self%nq = nq
  self%n  = n
  if (self%verbose > 0) then
    if (present(label)) then
      write (io_unit%log,'(1x,a,6i6,1x,a)') &
      'list_t%count_status: n[ibveaf] = ', &
      self%ni, self%nb, self%nv, self%ne, self%na, nf, label
    else
      write (io_unit%log,'(1x,a,6i6,1x,a)') &
      'list_t%count_status: n[ibveaf] = ', &
      self%ni, self%nb, self%nv, self%ne, self%na, nf
    end if
  end if
  self%na = n - nv - ne
  call self%set_io_counts
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
    call link%reset_status (link%nbor)
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
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%check_nbor_nbors')
  nbor => link%nbor
  do while (associated(nbor))
    call self%check_nbors (nbor%link)
    nbor => nbor%next
  end do
  call self%check_ready (link)
  call trace%end()
END SUBROUTINE check_nbor_nbors

!===============================================================================
!> Among a task and its neighbor tasks, move local tasks to ready_queue if they
!> are ready.  If the ready bit is already set it means the patch has already
!> been put in the ready queue, and should not be checked again.
!>
!> The link pointer and everything it points to are private to this task, and
!> are not at this point in time accessible from the ready queue.
!===============================================================================
SUBROUTINE check_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor, nbors
  class(task_t), pointer:: task
  integer:: itr, status
  integer, save:: itimer=0
  !.............................................................................
  task => link%task                                     ! main task
  itr = 1
  call trace%begin('list_t%check_nbors', itimer=itimer)
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,i4,f12.6,2x,4l1)') wallclock(), &
    'list_t%check_nbors: task', link%task%id, link%task%level, link%task%time, &
     link%task%is_set(bits%boundary), link%task%is_set(bits%virtual), &
     associated(link%nbor), associated(link%nbors_by_level)
  end if
  !-----------------------------------------------------------------------------
  ! Lock the other task link while accessing its nbor list, to prevent changes
  !-----------------------------------------------------------------------------
  nbor => link%nbor                                     ! first nbor
  do while (associated (nbor))                          ! keep going until end
    !---------------------------------------------------------------------------
    ! If the nbor needs this task, and isn't busy updating, check it.  For
    ! thread safety, the busy bit must NOT be cleared before a task has done
    ! a check_ready() on itself, lest an opportunity could be missed.
    !---------------------------------------------------------------------------
    if (self%verbose > 1) then
      write (io_unit%log,'(f12.6,1x,a,i6,i4,f12.6,2x,2l1)') wallclock(), &
        'list_t%check_nbors: nbor, level, NB =', nbor%task%id, &
        nbor%task%level, nbor%task%time, nbor%needs_me, &
        nbor%task%is_clear(bits%busy)
    end if
    if (nbor%needs_me .and. nbor%task%is_clear(bits%busy)) then
      call self%check_ready (nbor%link, lock=.true., status=status) ! pointer back
      if (self%verbose > 1) &
        write (io_unit%log,'(f12.6,i6,1x,a,i6,2i4)') wallclock(), task%id, &
        'checked nbor ', nbor%task%id, nbor%task%level, status
      if (status == 2) then
        !---------------------------------------------------------------------
        ! If the failure was due to an unsupported nbor, refresh its nbor
        ! list and try again -- it may have aquired support in the meantime.
        !---------------------------------------------------------------------
        call self%init_nbors (nbor%link)
        call self%check_ready (nbor%link, lock=.true., status=status)
        if (self%verbose > 1) &
          write (io_unit%log,'(f12.6,i6,2x,a,i6,2i4)') wallclock(), task%id, &
          'checked again', nbor%task%id, nbor%task%level, status
      end if
    else
      if (self%verbose > 1) &
        write (io_unit%log,'(f12.6,i6,2x,a,i6,i4)') wallclock(), task%id, &
        'no need to check', nbor%task%id, nbor%task%level
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
    call self%check_ready (link, status=status)                         ! RULE G
    if (self%verbose > 1) &
      write (io_unit%log,'(f12.6,i6,1x,a,i6,2i4)') wallclock(), task%id, &
      'checked self', task%id, task%level, status
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
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
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
  call self%lock%set ('check_all')
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
  call self%lock%unset ('check_all')
  if (associated(self%queue)) then
    write (io_unit%log,*) wallclock(), 'check_all: done, nq =', &
      self%nq, associated(self%queue)
  else
    write (io_unit%log,*) 'check_all: phase 2, clearing ready bits'
    call self%lock%set ('check_all(2)')
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
    call self%lock%unset ('check_all(2)')
  end if
  write (io_unit%log,*) wallclock(), 'check_all: done, nq =', &
    self%nq, associated(self%queue)
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
SUBROUTINE check_oldest (self, clear)
  class(list_t):: self
  logical, optional:: clear
  class(link_t), pointer:: link, oldest, oldestv
  class(nbor_t), pointer:: oldestnb, nbor
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
  call self%lock%set ('check_oldest')
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
  call self%lock%unset ('check_oldest')
  if (associated(oldest)) then
    write (io_unit%log,'(f12.6,2x,a,i6,f12.6," rb:",2l1)') wallclock(), &
      'the oldest active task is ', oldest%task%id , told, &
      oldest%task%is_set(bits%ready), oldest%task%is_set(bits%busy) 
  else
    write (stderr,*) 'STALLED: oldest active not found '
  end if
  flush (io_unit%log)
  if (associated(oldestv)) then
    write (io_unit%log,'(f12.6,2x,a,i6,f12.6)') wallclock(), &
      'the oldest virtual task is ', oldestv%task%id, toldv
    write (io_unit%log,'(f12.6,2x,a,f12.6)') wallclock(), &
      'it was last unpacked at', oldestv%task%unpack_time
    flush (io_unit%log)
  end if
  if (associated(oldest)) then
    nbor => oldest%nbor
    write (io_unit%log,2) oldest%task%id, oldest%task%time, oldest%task%rank, &
      oldest%task%level, &
      oldest%task%is_set(bits%boundary), oldest%task%is_set(bits%virtual)
  2 format(i6,f12.6,2i4,5L3,f12.6)
    write (io_unit%log,*) 'nbor list:'
    do while (associated(nbor))
      write (io_unit%log,2) nbor%task%id, nbor%task%time, nbor%task%rank, &
        nbor%task%level, nbor%needed, nbor%needs_me, &
        nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual), &
        nbor%task%is_ahead_of (oldest%task), nbor%task%unpack_time
      nbor => nbor%next
    end do
    if (present(clear)) then
      call oldest%task%clear (bits%ready+bits%busy)
      call self%queue_by_time (oldest)                  ! add task to queue
      write (io_unit%log,*) 'task queued'
      flush (io_unit%log)
    end if
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
      nbor => oldest%nbor
      do while (associated(nbor))
        if (nbor%task%time < told) then
          told =nbor%task%time
          oldestnb => nbor
        end if
        write (io_unit%log,1) &
          'nbor', nbor%task%id, 'rank', nbor%task%rank, nbor%task%time, &
          nbor%task%is_set(bits%boundary), nbor%task%is_set (bits%virtual)
        1 format (1x,a,i9,3x,a,i6,1p,e16.6,2l3)
        nbor => nbor%next
      end do
    end if
    write (io_unit%log,'(f12.6,2(3x,a,i8,i6,f12.6,3x,4l1),3i5)') wallclock(), &
     'oldest:', oldest%task%id, oldest%task%rank, oldest%task%time, &
     oldest%task%is_set(bits%internal), oldest%task%is_set(bits%boundary), &
     oldest%task%is_set(bits%virtual), oldest%task%is_set(bits%external), &
     'nbor:', oldestnb%task%id, oldestnb%task%rank, oldestnb%task%time, &
     oldestnb%task%is_set(bits%internal), oldestnb%task%is_set(bits%boundary), &
     oldestnb%task%is_set(bits%virtual), oldestnb%task%is_set(bits%external), &
     self%nq
    flush (io_unit%log)
    if (associated(oldestv)) then
      write (io_unit%log,*) 'the oldest virtual task is ', oldestv%task%id, &
        oldestv%task%time
      write (io_unit%log,*) 'check_oldest: phase 3', associated(oldestv)
      write (io_unit%log,1) 'the oldest virtual task is', oldestv%task%id, &
        'rank', oldestv%task%rank, oldestv%task%time, &
        oldest%task%is_set (bits%virtual), oldest%task%is_set (bits%external)
      nbor => oldestv%nbor
      do while (associated(nbor))
        write (io_unit%log,1) &
          'nbor', nbor%task%id, 'rank', nbor%task%rank, nbor%task%time, &
          nbor%task%is_set(bits%boundary), nbor%task%is_set (bits%virtual)
        nbor => nbor%next
      end do
    end if
  end if
  call trace%end()
END SUBROUTINE check_oldest

!===============================================================================
!> Check if the link task is ready, and if so, add it to self (the ready_queue)
!===============================================================================
SUBROUTINE check_ready (self, link, lock, status)
  class(list_t):: self
  class(link_t), pointer:: link
  logical, optional:: lock
  integer, optional:: status
  !.............................................................................
  logical:: ok
  class(task_t), pointer:: task
  class(link_t), pointer:: oldest
  class(nbor_t), pointer:: nbor
  real(8):: t_oldest
  integer:: ignore, unit, status_l
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we first check if any bits are set indicating this to be the case
  !-----------------------------------------------------------------------------
  if (self%detailed_timer .or. self%verbose > 0) &
    call trace%begin ('list_t%check_ready', itimer=itimer)
  status_l = -1
  task => link%task
  ignore = bits%ready + bits%remove + bits%virtual + bits%frozen
  if (task%is_clear (ignore)) then
    !...........................................................................
    if (self%verbose > 2) then
      write (io_unit%log,'(f12.6,1x,a,i6,2l3)') wallclock(), &
       'list_t%check_ready: task, lock, nbors', link%task%id, present(lock), &
       associated(link%nbor)
    end if
    !---------------------------------------------------------------------------
    ! Since only this thread can change the nbor list we are free to use it w/o
    ! lock in cases where task is the active task, but in other cases we need to
    ! lock the nbor list, so it isn't changed during the loop
    !---------------------------------------------------------------------------
    if (present(lock) .and. omp_lock%links) then
      if (self%detailed_timer) then
        call link%lock%set ('check_ready', itimer=itimer)
      else
        call link%lock%set ('check_ready')
      end if
    end if
    !---------------------------------------------------------------------------
    ! Tasks that are frozen or removed should not be used in the loop
    !---------------------------------------------------------------------------
    ignore = bits%frozen + bits%remove
    ok = .true.
    status_l = 0
    nbor => link%nbor                                   ! use original nbor list
    t_oldest = huge(1d0)
    do while (associated (nbor))                        ! keep going until end
      !-------------------------------------------------------------------------
      ! Check for any non-supported nbors and nix on them, setting also a bit
      ! that marks them to be ignored when checking nbors for readiness.  The
      ! bit is cleared when the task is again checked, by being an nbor in an
      ! nbor list used by check_nbors().  This requires that the level L-2 task
      ! that needs to be refined in order to provide support from an L-1 child
      ! is tested with a condition, list_t%check_support(), that discovers the
      ! need.
      !-------------------------------------------------------------------------
      if (nbor%task%level - link%task%level < -1) then
        status_l = 2
        ok = .false.
        call link%task%set (bits%no_support)
        !.......................................................................
        if (self%verbose > 1) &
          write (io_unit%log,'(f12.6,1x,a,i6,1x,a,i6,i4,f12.6)') wallclock(), &
            'list_t%check_ready:', link%task%id, &
            'L-2 nbor id, level, time =', nbor%task%id, nbor%task%level, &
            nbor%task%time
        !.......................................................................
      else if (nbor%task%is_set (bits%no_support)) then
        !.......................................................................
        if (self%verbose > 1) &
          write (io_unit%log,'(f12.6,1x,a,i6,1x,a,i6,i4,f12.6)') wallclock(), &
            'list_t%check_ready:', link%task%id, &
            'skip no support id, level, time =', nbor%task%id, nbor%task%level,&
            nbor%task%time         
        !.......................................................................
      else if (nbor%task%level - link%task%level > 1) then
        !.......................................................................
        if (self%verbose > 1) &
          write (io_unit%log,'(f12.6,1x,a,i6,1x,a,i6,i4,f12.6)') wallclock(), &
            'list_t%check_ready:', link%task%id, &
            'skip L+2 nbor id, level, time =', nbor%task%id, nbor%task%level
      !-------------------------------------------------------------------------
      ! Check nbors that are needed, not ignored, and have support
      !-------------------------------------------------------------------------
      else if (nbor%needed .and. nbor%task%is_clear(ignore)) then
        if (.not. nbor%task%is_ahead_of(task)) then       ! and is ahead in time
          ok = .false.
          status_l = max(1,status_l)
        end if
        !.......................................................................
        if (self%verbose > 1) &
          write (io_unit%log,'(f12.6,1x,a,i6,1x,a,i6,i4,1p,e12.3,l2)') &
            wallclock(), &
            'list_t%check_ready:', link%task%id, &
            'nbor, level, delta_t, ahead =', nbor%task%id, nbor%task%level, &
            (nbor%task%time-task%time)/max(nbor%task%dtime,tiny(1d0)), &
            nbor%task%is_ahead_of(task)
        !.......................................................................
      else
        !.......................................................................
        if (self%verbose > 1) &
          write (io_unit%log,'(f12.6,1x,a,i6,1x,a,i6,i4,2x,4l1)') &
            wallclock(), &
            'list_t%check_ready:', link%task%id, &
            'skips nbor, level, NIRF =', &
            nbor%task%id, nbor%task%level, &
            nbor%needed, nbor%task%is_set (ignore), &
            nbor%task%is_set (bits%remove ), &
            nbor%task%is_set (bits%frozen )
        !.......................................................................
      end if
      if (.not.ok) then
        if (nbor%task%time < t_oldest) then
          oldest => nbor%link; t_oldest = nbor%task%time
        end if
        !-----------------------------------------------------------------------
        ! verbose==1:  report which nbor it fails on
        ! verbose==2:  report more details
        ! verbose==3:  report all nbors it would fail on
        if (self%verbose < 4) then
          exit
        end if
      end if
      nbor => nbor%next                                   ! next nbor
    end do
    if (present(lock) .and. omp_lock%links) then
      call link%lock%unset ('check_ready')
    end if
    if (ok) then
      !$omp atomic write
      task%n_failed = 0
      call link%task%clear (bits%no_support)
      !.........................................................................
      if (self%verbose > 0) &
        write (io_unit%log,'(f12.6,1x,a,i6,f12.6)') wallclock(), &
        'list_t%check_ready: queue_by_time', task%id, task%time
      !.........................................................................
      call self%queue_by_time (link)                      ! add task to queue
    else
      !$omp atomic update
      task%n_failed = task%n_failed+1
      !-------------------------------------------------------------------------
      ! At verbose=1 we get a line in the rank log for every 10000 failures,
      ! while at verbose=2 we get a line for ever fail, and line on stderr for
      ! every 10000 failures
      !-------------------------------------------------------------------------
      if (self%verbose > 1  ) then
        write (io_unit%log,'(f12.6,i6,l3,f12.6,1x,a,i6,f12.6,i5,l3)') &
          wallclock(), task%id, task%is_set(bits%boundary), task%time, &
          'list_t%check_ready: failed with id, time, rank, V =', &
          oldest%task%id, t_oldest, oldest%task%rank, &
          oldest%task%is_set(bits%virtual)
      end if
      !.........................................................................
      if (task%n_failed > 10000) then
        if (self%verbose > 1) then
          unit = stderr
        else
          unit = io_unit%log
        end if
        write (unit ,'(a,i4,1p,2(2x,a,i6,1x,a,g14.5),2(2x,a,g11.3))') &
          'list_t%check_ready: on rank', task%rank, &
          'task'               , task%id     , 'at t=', task%time, &
          'failed on nbor task', nbor%task%id, 'at t=', nbor%task%time, &
          'grace =', task%grace, ' dt=', nbor%task%dtime
      end if
      !.........................................................................
    end if
  else if (self%verbose > 1) then
    !...........................................................................
    write (io_unit%log,'(f12.6,1x,a,i6,1x,a,5l1)') wallclock(), &
      'list_t%check_ready:', task%id, 'return because VFRr = ', &
      task%is_set (bits%virtual), &
      task%is_set (bits%frozen ), &
      task%is_set (bits%remove ), &
      task%is_set (bits%ready  )
    !...........................................................................
  end if
  if (present(status)) then
    status = status_l
  end if
  if (self%detailed_timer .or. self%verbose > 0) &
    call trace%end (itimer)
END SUBROUTINE check_ready

!===============================================================================
!> Check if the number of nbor tasks is consistent
!===============================================================================
SUBROUTINE consistency (self, link, i)
  class(list_t):: self
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
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
  class(nbor_t), pointer:: nbor
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
        'ERROR: queue out of order ident, prev%id, time, task%id', &
        ident, prev%id, time, task%id, task%time
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
  class(nbor_t), pointer:: nbor
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
        'ERROR: active queue ident, prev%id, time, task%id, task%atime =', &
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
  class(nbor_t), pointer:: nbor
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
            if (all(abs(scan%task%position-link%task%position) &
                < scan%task%size*0.55d0)) then
              call scan%task%set (bits%not_leaf)        ! mark non-leaf
              if (io%verbose > 0) &
                write (io_unit%log,*) 'patch:', link%task%id, &
                                      ' parent:', scan%task%id
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
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%add_new_link')
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,2i6)') wallclock(), &
     'list_t%add_new_link: task', link%task%id
  end if
  !-----------------------------------------------------------------------------
  ! Make sure class(patch_t) tasks have a pointer back to the task link
  !-----------------------------------------------------------------------------
  task => link%task
  select type (task)
  class is (patch_t)
    task%link => link
  end select
  !-----------------------------------------------------------------------------
  ! Increment task counters on the rank
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%virtual)) then
    !$omp atomic update
    timer%levelpop (task%level) = &          ! increment level population count
    timer%levelpop (task%level) + 1
  end if
  !-----------------------------------------------------------------------------
  ! Generate a new nbor list -- which includes running update_status() and
  ! check_ready() on it
  !-----------------------------------------------------------------------------
  call self%init_nbors (link)                                           ! RULE B
  !-----------------------------------------------------------------------------
  ! Add the task to the task_list, and if it turns out to be a boundary task,
  ! send it to rank nbors.
  !-----------------------------------------------------------------------------
  call self%append (link)
  if (link%task%is_set (bits%boundary)) &
    call self%send_to_vnbors (link)
  call self%set_init_nbors (link)                                       ! RULE C
  !-----------------------------------------------------------------------------
  ! A task that has a new nbor list needs to be checked to see if is possibly
  ! is ready to be added to the queue
  !-----------------------------------------------------------------------------
  call self%check_ready (link)
!  call self%check_nbors (link)                                          ! RULE D
  call trace%end ()
END SUBROUTINE add_new_link

!===============================================================================
!> Initialize the nbor list of a task at level L:  
!> 1) Initially, levels L+-2 are entered into the nbor list, and are 
!>    subsequently checked by prune_nbors(), which removes nbors at L+2 that
!>    have support at L+1, and nbors at L-2, if the current task has support
!>    at L-1.
!> 2) Tasks that have nbors at level L-2 lack support, and should be rejected
!>    when tested by check_ready(). The bits%no_support are set for such tasks.
!> 3) Task with bits%no_support set are not counted when check_ready loops over
!>    nbors to other tasks
!> 4) When a level L task that has a level L+2 nbor is checked by refine_mod,
!>    it will be refined, creating an L+1 nbor that provides support for the
!>    L+2 task.  When the last level L nbor of the L+2 task is refined, the
!>    level L+2 task has support, and the bits%no_support should be cleared.
!> 5) Level L task with level L+2 nbors are fine, and need no special treatment,
!>    except they should not download data from its L+2 nbors; this is enforced
!>    by giving any L+-2 nbor a setting nbor%download=.false.
!===============================================================================
SUBROUTINE init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: scan
  class(nbor_t), pointer:: nbor, new_head, new_sort, old_head, old_sort, nbors
  integer:: n_add, status
  logical:: overlaps
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%init_nbors ', itimer=itimer)
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,2i6)') wallclock(), &
     'list_t%init_nbors: task', link%task%id
  end if
  nullify (new_head, new_sort)
  n_add = 0
  call self%lock%set ('init_nbors', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! A new task, or one that changes from boundary to internal, needs to have the
  ! internal bit set initially, to be removed if finding a virtual nbor
  !-----------------------------------------------------------------------------
  call link%task%set (bits%internal)
  scan => self%head                                   ! scan all links (FIXME)
  do while (associated(scan))                         ! until end
    if (.not.associated(scan,link) .and. associated(scan%task)) then
      !-------------------------------------------------------------------------
      ! The ifs select to add nbors to a patch if either the patch is not a
      ! virtual patch, or the nbor is not a virtual patch, thus excluding
      ! adding virtual patches as nbors to virtual patches
      !-------------------------------------------------------------------------
      if (link%task%rank == mpi%rank .or. scan%task%rank == mpi%rank) then   
        !-----------------------------------------------------------------------
        ! Potential nbors are other patches that overlap
        !-----------------------------------------------------------------------
        overlaps = scan%task%overlaps(link%task) &
             .and. abs(scan%task%level - link%task%level) <= 2
        if (overlaps) then  ! if its task overlaps
          allocate (nbor)
          nbor%link => scan
          nbor%task => scan%task
          !---------------------------------------------------------------------
          ! Set status bits and flags, but note that these may have to be reset
          ! after pruning, since that may for example remove the only virtual
          ! nbor, which means the task that the list belongs to changes from a
          ! boundary task to an internal task.
          !---------------------------------------------------------------------
          call link%task%nbor_relations (nbor%task, nbor%needed, nbor%needs_me, &
            nbor%download)
          !---------------------------------------------------------------------
          ! Add in rank order and increment the count of tasks need this task
          !---------------------------------------------------------------------
          call link%add_nbor_by_rank (new_head, nbor)
          !$omp atomic update
          scan%task%n_needed = scan%task%n_needed + 1
          !.....................................................................
          n_add = n_add+1
          if (self%verbose > 1 .or. link%task%id == io%id_debug) &
            write (io_unit%log, '(2x,a,2i6,i4,3x,4l1)') &
            'init_nbors: id, nbor, level, NMDV =', link%task%id, &
            nbor%task%id, nbor%task%level, nbor%needed, nbor%needs_me, &
            nbor%download, nbor%task%is_set(bits%virtual)
          !.....................................................................
        end if
      end if
    end if
    scan => scan%next                               ! continue search
  end do
  call self%prune_nbors (link, new_head, status)
  !call self%lock%unset ('init_nbors')
  !.............................................................................
  if (self%verbose > 0 .or. link%task%id == io%id_debug) then
    write(io_unit%log,'(f12.6,1x,a,i6,i4,3x,3l1,2i4)') wallclock(), &
      'init_nbors: link, rank, IBV, n_nbor, status =', &
      link%task%id, link%task%rank, &
      link%task%is_set(bits%internal), &
      link%task%is_set(bits%boundary), &
      link%task%is_set(bits%virtual), n_add, status
    flush (io_unit%log)
  end if
  !-----------------------------------------------------------------------------
  ! The links only needs to be locked the brief moment while the nbor pointers
  ! are changed.  
  !-----------------------------------------------------------------------------
  nullify (new_sort)
  call link%sort_nbors_by_level (new_head, new_sort)
  !-----------------------------------------------------------------------------
  ! Under a brief link lock, swap the nbor lists, and update the status bits, in
  ! case for example an internal task became a boundary task
  !-----------------------------------------------------------------------------
  call link%lock%set ('init_nbors', itimer=itimer)
  old_head => link%nbor
  old_sort => link%nbors_by_level
  link%nbor => new_head                               ! switch nbor list
  link%nbors_by_level => new_sort
  call self%update_status (link)                      ! based on new list
  link%task%n_nbors = n_add                           ! remember n_nbors
  call link%lock%unset ('init_nbors')
  call self%lock%unset ('init_nbors')
  !-----------------------------------------------------------------------------
  ! Remove the two released nbor lists, to recover the memory
  !-----------------------------------------------------------------------------
  call link%remove_nbor_list (old_head)               ! recover nbor list mem
  call link%remove_nbor_list (old_sort)
  !-----------------------------------------------------------------------------
  ! Clear the init_nbors bit, unless the task is a boundary task, in which
  ! case the bit needs to be present also on the nbor ranks, where it will be
  ! cleared at this point.   This creates a potential problem, since if a task
  ! that still has the init_nbors bit set, after coming in with the swap_reuqest
  ! bit set, it will never clear its init_nbors bit.
  !-----------------------------------------------------------------------------
  if (link%task%is_clear (bits%boundary)) &
    call link%task%clear (bits%init_nbors)
  call trace%end (itimer)
END SUBROUTINE init_nbors

!===============================================================================
!> Remove nbors at level L+2 that are contained in another nbor at level L+1.
!> However, this may not be necessary and sufficient.  If an L+2 nbor that lacks
!> support has guard zones that extend into the link patch then it is the
!> patch that needs refinement.  Conversely, if an L+2 nbor lacks support in
!> guard zones that lie outside this patch, then it is an nbor level L patch
!> that needs to produce a new L+1 patch to support it.
!===============================================================================
SUBROUTINE prune_nbors (self, link, head, status)
  class(list_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: head
  integer, optional:: status
  !.............................................................................
  class(nbor_t), pointer:: nbor, next, nbor2
  class(task_t), pointer:: task
  real, pointer, contiguous:: irefine(:,:,:)
  logical:: removed
  integer:: n(2)
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%prune_nbors')
  removed = .false.
  n = 0
  !-----------------------------------------------------------------------------
  ! If link%task has support, then remove any L-2 nbors
  !-----------------------------------------------------------------------------
  task => link%task
  select type (task)
  class is (patch_t)
  allocate (irefine(task%gn(1),task%gn(2),task%gn(3)))
  end select
  if (check_support (self, link, -1, irefine) <= 0) then
    nbor => head
    do while (associated(nbor))
      next => nbor%next
      if (nbor%task%level == link%task%level-2) then
        !-------------------------------------------------------------------------
        if (self%verbose > 1) then
          write (io_unit%log,'(f12.6,1x,a,i6,i4)') wallclock(), &
            'list_t%prune_nbors: L-2 task id, level =', &
            nbor%task%id, nbor%task%level
          flush (io_unit%log)
        end if
        call link%remove_nbor (nbor%task, head)
        removed = .true.
      end if
      nbor => next
    end do
    call link%task%clear (bits%no_support)  
  else
    call link%task%set (bits%no_support)
  end if
  !-----------------------------------------------------------------------------
  ! If the task has L+2 nbors, check if those tasks have support, and remove
  ! them if they do
  !-----------------------------------------------------------------------------
  nbor => head
  do while (associated(nbor))
    next => nbor%next
    if (nbor%task%level == link%task%level+2) then
      !.........................................................................
      if (self%verbose > 1) then
        write (io_unit%log,'(f12.6,1x,a,i6,i4)') wallclock(), &
          'list_t%prune_nbors: L+2 task id, level =', &
          nbor%task%id, nbor%task%level
        flush (io_unit%log)
      end if
      !-------------------------------------------------------------------------
      ! Prune any level L+2 task from nbor list, if it has support
      !-------------------------------------------------------------------------
      n(2) = n(2)+1
      if (self%check_support (nbor%link, -1, irefine) <= 0) then
        !.......................................................................
        if (self%verbose > 1) then
          write (io_unit%log,'(f12.6,1x,a,i6,i4)') wallclock(), &
            'list_t%prune_nbors: removing L+2 id, level =', &
            nbor%task%id, nbor%task%level
          flush (io_unit%log)
        end if
        !.......................................................................
        call link%remove_nbor (nbor%task, head)
        n(2) = n(2)-1
        removed = .true.
      end if
    end if
    nbor => next
  end do
  deallocate (irefine)
  if (removed) then
    !---------------------------------------------------------------------------
    ! If any task was removed, the status bits and flags need to be recomputed,
    ! and the task list counts need to be updated.
    !---------------------------------------------------------------------------
    call self%update_status (link)
  end if
  if (present(status)) then
    status = n(1)
  end if
  !.............................................................................
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,3i4)') wallclock(), &
      'list_t%prune_nbors: id, level, n =', &
      link%task%id, link%task%level, n
    flush (io_unit%log)
  !.............................................................................
  end if
  call trace%end ()
END SUBROUTINE prune_nbors

!===============================================================================
!> Set bits%init_nbors on nbor tasks, if the task itself is new (and hence has
!> a new nbor list, which should be matched by the nbors nbor list). 
!===============================================================================
SUBROUTINE set_init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%set_init_nbors')
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,2i6)') wallclock(), &
     'list_t%set_init_nbors: task', link%task%id
  end if
  nbor => link%nbor
  do while (associated(nbor))
    !call nbor%link%lock%set ('set_init_nbors')
    call nbor%task%set (bits%init_nbors)
    if (self%verbose > 1) &
      write (io_unit%log,*) 'set_init_nbors:', link%task%id, nbor%task%id, &
        nbor%task%istep, nbor%task%is_set(bits%init_nbors)
    !call nbor%link%lock%unset ('set_init_nbors')
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE set_init_nbors

!===============================================================================
!> Trigger update_status() on nbors
!===============================================================================
SUBROUTINE set_update_status (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%set_update_status')
  nbor => link%nbor
  do while (associated(nbor))
    !call nbor%link%lock%set ('set_update_status')
    call nbor%task%set (bits%update_status)
    if (self%verbose > 1) &
      write (io_unit%log,'(1x,a,2i6)') &
        'list_t%set_update_status: task, nbtask =', link%task%id, &
        nbor%task%id, nbor%task%is_set(bits%update_status)
    !call nbor%link%lock%unset ('set_update_status')
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE set_update_status

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
  class(nbor_t), pointer:: nbor, nbors
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
  class(nbor_t), pointer:: nbor1, nbor2
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
    end if
    i = i+1                                             ! increment counter
    link => link%next                                   ! continue list
  end do
  !$omp end parallel
  !call self%reset_status                                ! set bits and count
  call timer%toc ('time to initialize nbor lists', 1, time)
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
        if (parent%task%is_set (bits%has_output)) then
          io%nwrite = io%nwrite - 1
        end if
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
    where (task1%periodic)
      dist = modulo (dist + 0.5*self%size, self%size) - 0.5*self%size
    end where
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
  call self%lock%set ('append_list')
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
  call self%lock%unset ('append_list')
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
  if (self%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' wait for tasklist(3)'
  call self%lock%set ('remove_task')
  if (self%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(3)'
  !-----------------------------------------------------------------------------
  ! First remove the task from the ready queu
  !-----------------------------------------------------------------------------
    link => self%queue
    nullify (prev)
    do while (associated(link))
      if (link%task%id == task%id) then
        if (associated(prev)) then
          prev%next_time => link%next_time
        else
          self%queue => link%next_time
        end if
        exit
      end if
      prev => link
      link => link%next_time
    end do
  !-----------------------------------------------------------------------------
  ! Now remove the task from the full task list
  !-----------------------------------------------------------------------------
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
  if (self%verbose > 3) &
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
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (self%detailed_timer) then
    call trace%begin ('list_t%queue_by_time ', itimer=itimer)
    if (self%verbose > 3) &
      write (io%output,*) wallclock(),' thread',omp%thread,' waitfor tasklist(4)'
    call self%lock%set ('queue_by_time', itimer=itimer)
  else
    if (self%verbose > 3) &
      write (io%output,*) wallclock(),' thread',omp%thread,' waitfor tasklist(4)'
    call self%lock%set ('queue_by_time')
  end if
  if (self%verbose > 3) &
    write (io%output,*) wallclock(),' thread',omp%thread,' locked tasklist(4)'
  call real_queue_by_time (self, this, error)
  if (io%task_logging > 2) then
    block
    character(len=24):: label
    write (label,'(a)') 'queue  '
    call this%task%log (label)
    end block
  end if
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
  if (self%detailed_timer) then
    call trace%end (itimer)
  end if
END SUBROUTINE queue_by_time

SUBROUTINE real_queue_by_time (self, this, error)
  class(list_t):: self
  class(link_t), pointer:: this
  logical:: error
  class(link_t), pointer:: next, prev, link
  class(task_t), pointer:: task
  integer:: nit, i, j
  !-----------------------------------------------------------------------------
  ! Refuse to queue virtual tasks
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
  !-----------------------------------------------------------------------------
  !$omp atomic
  mpi_mesg%n_ready = mpi_mesg%n_ready+1
  !-----------------------------------------------------------------------------
  ! Lock the queue while inserting the task, and mark the task as being in the
  ! ready queue (or executing -- the bit is cleared only after updating).
  ! Double check that no other thread has set the ready bit on this task.
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%ready+bits%busy)) then
    call task%set(bits%ready)
    if (self%verbose > 1) &
      write (io_unit%log,'(f12.6,i5,1x,a,i9,a,f12.6,2x,a,i5)') wallclock(), omp%thread, &
        'list_t%queue_by_time: adding task', task%id, ' at time', task%time, &
        'nq', self%nq
    !call this%qlock%set ('queue_by_time')
    nullify (prev)
    !call self%queue%qlock%set ('queue_by_time')
    next => self%queue
    !call self%queue%qlock%unset ('queue_by_time')
    nit = 0
    do while (associated(next))
      nit = nit+1
      if (nit > self%nq) then
        print *,'ERROR: hang in queue_by_time, nit =', nit
        error = .true.
        return
      end if
      if (associated(next%task, task)) then
        write (io_unit%log,*) omp_mythread, ' WARNING: task', task%id, &
          ' is already in ready queue'
        go to 9
      else if (next%task%time > task%time) then
        !call this%qlock%set ('queue_by_time')
        this%next_time => next
        !call this%qlock%unset ('queue_by_time')
        if (associated(prev)) then
          if (self%verbose > 1) then
            write (io_unit%log,'(i4,2x,a,i6,a,i6,i5,1p,3g16.6)') &
              omp%thread, 'task',task%id,' in ready queue between', &
              prev%task%id,next%task%id, prev%task%time, this%task%time, &
              next%task%time
          end if
          !call prev%qlock%set ('queue_by_time')
          prev%next_time => this
          !call prev%qlock%unset ('queue_by_time')
        else
          !call self%queue%qlock%set ('queue_by_time')
          self%queue => this
          !call self%queue%qlock%unset ('queue_by_time')
          if (self%verbose > 1) then
            write (io_unit%log,*) 'task',task%id,' at ready queue head',self%nq
          end if
        end if
        !$omp atomic
        self%nq = self%nq+1
        go to 9
      end if
      !call next%qlock%set ('queue_by_time')
      prev => next
      next => next%next_time
      !call next%qlock%unset ('queue_by_time')
    end do
    self%nq = self%nq+1
    if (associated(prev)) then
      if (self%verbose > 0) then
        1 format(1x,a,3(a,i6,5x))
        write (io_unit%log,1) 'list_t%queue_by_time: ', &
        'task',task%id, 'in queue after', prev%task%id, 'nq:', self%nq
      end if
      prev%next_time => this
    else
      if (self%verbose > 0) then
        write (io_unit%log,1) 'list_t%queue_by_time: ', &
        'task', task%id, 'at queue head    nq:', self%nq
      end if
      self%queue => this
    end if
    nullify (this%next_time)
    !...........................................................................
    if (self%verbose > 1) then
      write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,2i6)') &
      wallclock(), omp%thread, 'queued', task%id, task%time, &
      task%istep, self%nq
    end if
  9 continue
    call self%remove_active (this)                    ! remove from active_queue
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
    write (io_unit%log,'(1x,i4,1x,a,i9,a,f10.5,i5,l4)') omp%thread, &
      'list_t%queue_active: adding task', task%id, ' at time', task%atime, &
      self%nac, associated(self%active)
  nullify (prev)
  !$omp critical (active_cr)
  self%nac = self%nac+1
  next => self%active
  do while (associated(next))
    if (associated(next%task, task)) then
      write (io_unit%log,*) omp_mythread, ' WARNING: task', task%id, &
        ' is already in active queue'
      exit
    else if (next%task%atime > task%atime) then
      this%next_active => next
      if (associated(prev)) then
        if (io%verbose > 1) then
          write (io_unit%log,*) 'task',task%id,' in active queue after', &
            prev%task%id,self%nac
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
        write (io_unit%log,*) 'task',task%id,' at active queue tail', &
          prev%task%id,self%nac
      end if
      prev%next_active => this
    else
      if (io%verbose > 1) then
        write (io_unit%log,*) 'task',task%id,' at active queue head', self%nac
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
!> Remove a task link from the active queue
!===============================================================================
SUBROUTINE remove_from_queue (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  class(link_t), pointer:: next, prev
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%remove_from_queue')
  nullify (prev)
  next => self%queue
  do while (associated(next))
    if (associated(next%task, this%task)) then
      if (associated(prev)) then
        prev%next_time => next%next_time
      else
        self%queue => next%next_time
      end if
      call this%task%clear (bits%ready)
      call self%count_status
      !if (self%verbose > 0) &
        write (io_unit%mpi,*) &
         'remove_from_queue: id, nq =', this%task%id, self%nq
      call trace%end ()
      return
    end if
    prev => next
    next => next%next_time
  end do
  call io%warning ('list_t%remove_from_queue did not find task')
  call trace%end()
END SUBROUTINE remove_from_queue

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
      write (io_unit%log,*)'ERROR: Trying to add task', this%task%id, &
        ' which is already there'
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
  !.............................................................................
  call self%print_queue_times ('print_queue')
END SUBROUTINE print_queue

!===============================================================================
!> Print the queue until and including a given task
!===============================================================================
SUBROUTINE print_queue_until (self, task)
  class(list_t):: self
  class(task_t):: task
  class(link_t) , pointer:: link
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
    write (io_unit%log,'("queue: id, time =",i7,1p,g14.5,l5)') link%task%id, &
      link%task%time, link%task%is_set(bits%ready)
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
  class(link_t) , pointer:: link
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
      write (io_unit%log,'(a,1x,i8,g15.6,3x,3l1)') &
        'list_t%print_tasks: id, time, RBV =', &
        link%task%id, link%task%time, &
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
!> Check nbor list of self, verifying that %task and %link are consistent
!===============================================================================
SUBROUTINE check_nbor_list (link, label, verbose)
  class(link_t):: link
  character(len=*):: label
  integer, optional:: verbose
  class(nbor_t), pointer:: nbor
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
            write (io_unit%log,*) label, link%task%id, task%id, &
              'ERROR: nbor%link, task%link not associated'
          end if
        end select
      else
        write (io_unit%log,*) &
          label, link%task%id, 'ERROR: nbor%task not associated'
      end if
    else
      write (io_unit%log,*) &
        label, link%task%id, 'ERROR: nbor%link not associated'
    end if
    if (present(verbose)) then
      task => nbor%link%task
      select type (task)
      class is (patch_t)
        write (io_unit%log,*) link%task%id, nbor%task%id, task%id, &
        associated(task%link,nbor%link)
      end select
    end if
    nbor => nbor%next
  end do
END SUBROUTINE check_nbor_list

!===============================================================================
!> Since there are corresponding actions on the two sides, it is best to handle
!> both sides with one procedure, and it needs to be list_t procedure, since it
!> must be called from both load_balance_t and task_mesh_t.
!>
!> NBOR LISTS:
!> Tasks that change from boundary to virtual, or the other way around, need to
!> have their nbor lists reconstructed, since the nbor lists of virtual tasks
!> do not include other virtual task, but the nbor lists of boundary tasks do.
!>
!> The nbor lists of nbors also need to be updated, since the settings of the
!> status bits and download flags etc may need to be change.
!>
!> STATUS BITS:
!> The swap_request bit is set to mark the request, and needs to be cleared on
!> both the sender and receiver side.  This could be done as part of the normal
!> processing that all tasks go through.
!>
!> The virtual and boundary bits need to be swapped, and if this is done before
!> sending the task it only needs to be done once.
!>
!> PROCESS COUNTS
!> The %na  (active count): decrease on sender and increase on recipient
!> The %n    (total count): decrease on sender and increase on recipient
!> THe %nv (virtual count): increase on sender and decrease on recipient
!> The %nb   (bndry count): decrease on sender and increase on recipient
!===============================================================================
SUBROUTINE swap (self, link, sender, rank)
  class(list_t)          :: self
  class(link_t), pointer :: link
  logical                :: sender
  integer, optional      :: rank
  !.............................................................................
  integer           :: one
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%swap')
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,2i6)') wallclock(), &
     'list_t%swap: task', link%task%id
  end if
  one = merge(+1,-1,sender)
  !$omp atomic update
  self%na = self%na-one                  ! active tasks
  !$omp atomic update
  self%nb = self%nb-one                  ! boundary tasks
  !$omp atomic update
  self%nv = self%nv+one                  ! virtual tasks
  call self%set_io_counts                ! apply counts
  !$omp atomic update
  timer%levelpop (link%task%level) = &    ! level population count
  timer%levelpop (link%task%level) - one
  !-----------------------------------------------------------------------------
  ! Set the rank to the new owner on the local rank, and flip the swap_request
  ! bit. 
  !-----------------------------------------------------------------------------
  if (sender) then
    link%task%rank = rank
    call link%task%set (bits%swap_request)
    call change_running_mean (-0.1)
    if (self%verbose > 0) &
      write (io_unit%log,*) 'list_t%swap: to rank =', link%task%id, rank
  else
    link%task%rank_steps = 0
    if (self%verbose > 0) then
      write (io_unit%log,'(a,i6,3x,"ibvef =",1x,5l1)') &
        'list_t%swap: bits in', link%task%id, &
        link%task%is_set(bits%internal), &
        link%task%is_set(bits%boundary), &
        link%task%is_set(bits%virtual ), &
        link%task%is_set(bits%external), &
        link%task%is_set(bits%frozen  )
    end if
    call link%task%clear (bits%swap_request+bits%busy+bits%ready)
    call link%task%clear (bits%virtual)
    call link%task%set (bits%boundary)
    call change_running_mean (+0.1)
    if (self%verbose > 0) then
      write (io_unit%log,'(a,i6,3x,"ibvef =",1x,5l1)') &
        'list_t%swap: bits in', link%task%id, &
        link%task%is_set(bits%internal), &
        link%task%is_set(bits%boundary), &
        link%task%is_set(bits%virtual ), &
        link%task%is_set(bits%external), &
        link%task%is_set(bits%frozen  )
    end if
  end if
  !-----------------------------------------------------------------------------
  ! With the new rank and status settings, update the status bits of the link
  ! task, and also trigger updates of the status bits on nbor tasks.
  !-----------------------------------------------------------------------------
  !call self%update_status (link)
  !call self%set_update_status (link)
  !call self%count_status                    ! FIXME: redundant (and costly)?
  !call self%check_nbors (link)
  call self%init_nbors (link)               ! power hammer
  call self%set_init_nbors (link)
  call self%check_ready (link)
  call trace%end ()
contains
  !=============================================================================
  !> Facilitate immediate changing of the queue information
  !=============================================================================
  subroutine change_running_mean (increment)
    real:: increment
    real::  old_measure, old_queue, measure, queue
    integer:: n_task
    !---------------------------------------------------------------------------
    n_task = self%na
    old_measure = self%running_mean%average
    old_queue = n_task/old_measure - 1.0 
    queue = old_queue + increment
    self%mpi_window%value = queue
    measure = n_task/(queue + 1.0)
    call self%running_mean%update (measure, all=.true.)
    if (self%verbose > 0) then
      write (io_unit%log, '(1x,a,i6,4f8.1)') &
        'list_t%swap: task, old_measure, old_queue, measure, queue =', &
        link%task%id, old_measure, old_queue, measure, queue
    end if
  end subroutine change_running_mean
END SUBROUTINE swap

!===============================================================================
!> Update the nbor/vnbor status of a links nbor list, including the self link.
!> If the boundary bit gets set, send the task to its vnbors.  Before doing so,
!> the nbor list of the nbors need to be updated, to reflect the new task rank
!===============================================================================
SUBROUTINE update_nbor_status (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, this
  class(nbor_t), pointer:: nbor
  integer:: status
  !-----------------------------------------------------------------------------
  call trace%begin ('link_t%update_nbor_status')
  write (stderr,*) 'WARNING: list_t%update_nbor_status may not be up-to-date!'
  call link%lock%set ('update_nbor_status')
  !
  nbor => link%nbor                                     ! then help others
  do while (associated(nbor))
    if (io%verbose>1) &
      write (io_unit%log,*) 'update_nbor_status:', nbor%task%id, &
      nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual)
    status = nbor%task%status                           ! current status
    call io%abort ('disabled for now')
    !call link%remove_nbor (nbor%task, link%nbor)        ! remove out-of-order
    !call nbor%link%add_nbor_by_rank (link%nbor, link)   ! add in rank order
    call self%update_status (nbor%link)                 ! new status
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
  if (io%verbose> 1 ) then
    call self%test_nbor_status (link)
  end if
  call link%lock%unset ('update_nbor_status')
  call trace%end()
END SUBROUTINE update_nbor_status

!===============================================================================
!> Update the nbor/vnbor status of a link.  If self belongs to the local rank
!> it is either an internal or boundary rank, if not it is virtual or external,
!> depending on what nbors it has. 
!>
!> NOTE: This routine assumes that the initial task status bits are consistent
!> with the task list counts, before the bits were updated.  Therefore, if the
!> task is new, it should have none of the IBV bits set.
!===============================================================================
SUBROUTINE update_status (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  type (task_t):: old
  integer:: n
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%update_status')
  !call link%lock%set ('update_status')
  nbor => link%nbor
  do while (associated(nbor))
    old%status = link%task%status
    call link%task%nbor_relations (nbor%task, nbor%needed, nbor%needs_me, &
                                   nbor%download)
    call update_count (link%task, old, bits%virtual , self%nb)
    call update_count (link%task, old, bits%boundary, self%nv)
    call update_count (link%task, old, bits%internal, self%na)
    nbor => nbor%next
  end do
  n = self%na + self%nv
  !$omp atomic write
  self%n = n
  !-----------------------------------------------------------------------------
  ! Clear the update_status bit, unless the task is a boundary task, in which
  ! case the bit needs to be present also on the nbor ranks, where it will be
  ! cleared at this point.
  !-----------------------------------------------------------------------------
  if (link%task%is_clear (bits%boundary)) &
    call link%task%clear (bits%update_status)
  !call link%lock%unset ('update_status')
  call trace%end ()
contains
  subroutine update_count (old, new, bit, count)
    class(task_t):: old, new
    integer:: bit, count
    if (old%is_set(bit) .and. new%is_clear (bit)) then
      !$omp atomic update
      count = count-1
    else if (new%is_set(bit) .and. old%is_clear (bit)) then
      !$omp atomic update
      count = count+1
    end if
  end subroutine update_count
END SUBROUTINE update_status

!===============================================================================
!> Send a package with updated time slice info to neighbor ranks.
!>
!> 1) assemble task data into a message
!> 3) send the message to the nbors rank (only once per rank)
!> 4) add the message to a sent_list
!> 5) check the sent list for completed sends
!> 6) remove the message from the sent_list and disassemble the message
!>
!> To ensure scheduling consistency, the aspects below are important:
!>
!> 1) Since load balancing may have changed the nature of nbors of this task,
!> and might have transformed this task from a boundary task to an internal task
!> a list_t%update_status() must be performed, to re-establish the status of
!> the current task.
!> 
!> 2) During the brief time that this task is deciding who -- if anyone -- to
!> send to, load balancing must be forbidden, via a lock on load_balance_t,
!> which also must be applied during the handling of a swap request on the
!> receiver side.
!>
!> 3) Packages sent must also be received, since they are expected in sequential
!> order.  Therefor, when a rank takes ownership of a task, the swap message is
!> the final message, and it should be preventing the next recv request from
!> being waited upon.  Instead, that particular recv request should be canceled.
!===============================================================================
SUBROUTINE send_to_vnbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  class(mesg_t), pointer:: mesg
  character(len=24):: label
  integer:: ierr, rank, tag
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  task => link%task
  if (task%is_clear (bits%boundary)) then
    if (task%is_clear (bits%swap_request)) then
      write (stdout,*) mpi%rank, omp%thread, &
        'ERROR: trying to send non-boundary task', task%id
      return
    end if
  end if
  if (task%is_set (bits%virtual)) then
    if (task%is_clear (bits%swap_request)) then
      write (stdout,*) mpi%rank, omp%thread, &
        'ERROR: trying to send virtual task', task%id
      return
    end if
  end if
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%send_to_vnbors', itimer=itimer)
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,f12.6)') wallclock(), &
     'list_t%send_to_vnbors: task', link%task%id, link%task%time
  end if
  !-----------------------------------------------------------------------------
  ! Pack task into mesg, this also allocates task%mesg
  !-----------------------------------------------------------------------------
  select type (task)
  class is (patch_t)
  call task%pack ()
  mesg => task%mesg
  end select
  !-----------------------------------------------------------------------------
  ! The first time a task is sent it has tag = task%id*100 + 1
  !-----------------------------------------------------------------------------
  if (mpi_mesg%uniq_mesg .and. mpi_mesg%tag_type == 1) then
    !---------------------------------------------------------------------------
    ! A new child task has task%seq=-1, which becomes 0 in task_mesg_t%send_to..
    !---------------------------------------------------------------------------
    if (task%seq == 0) then
      tag = 0
    else
      tag = task%id*100 + mod(task%seq,100)
    end if
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
  call link%lock%set ('send_to_vnbors')
  nbor => link%nbor                                             ! first nbor
  if (self%verbose > 1) &
    write (io_unit%log,*) link%task%id, 'associated', associated(nbor)
  do while (associated(nbor))                                   ! until end
    nbtask => nbor%task
    if (self%verbose > 1) &
      write (io_unit%log,*) nbtask%id, 'ranks', mpi%rank, nbtask%rank, rank
    if (nbtask%rank/=mpi%rank .and. nbtask%rank/=rank) then     ! new rank?
      rank = nbtask%rank                                        ! current rank
      if (mpi_mesg%uniq_mesg .and. mpi_mesg%tag_type == 2) then
        !$omp atomic capture
        sequence(rank) = sequence(rank)+1
        mesg%seq = sequence(rank)
        !$omp end atomic
        tag = mod(mesg%seq,100)
      end if
      mesg%tag = tag
      if (io%task_logging > 1) then
        write (label,'(a,i4,i8)') 'send   ', rank, mod(tag,100)
        call task%log (label)
      end if
      call mesg%send (rank, tag=tag)                            ! send it
      if (self%verbose > 0 .or. task%id == io%id_debug) then
        write (io_unit%log,'(f12.6,1x,a,i9,1p,e18.6,i9,2x,a,i5,2x,5l1)') &
          wallclock(), 'send_to_vnbors: id, time, tag =', &
          mesg%id, task%time, tag, 'to', rank, &
          task%is_set (bits%internal), &
          task%is_set (bits%boundary), &
          task%is_set (bits%virtual), &
          task%is_set (bits%external), &
          task%is_set (bits%swap_request)
        flush (io_unit%mpi)
      end if
      if (self%verbose > 1) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,i6)') &
          wallclock(), omp%thread, 'sent', task%id, task%time, rank
      end if
    end if
    nbor => nbor%next                                           ! next nbor
  end do
  call link%lock%unset ('send_to_vnbors')
  !-----------------------------------------------------------------------------
  ! Add the message to mesg%sent_list, if it was sent, and clear the swap bit
  !-----------------------------------------------------------------------------
!  call mpi_mesg%sent (mesg)                                     ! add to sent_list
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
  class(nbor_t), pointer:: nbor
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
  class(nbor_t), pointer:: nbor
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
SUBROUTINE info (self, label)
  class(list_t):: self
  class(link_t), pointer:: link
  character(len=*), optional:: label
  !.............................................................................
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%info')
  if (present(label)) then
    write (io_unit%log,*) 'list_t%info: '//trim(label)
  else
    write (io_unit%log,*) 'list_t%info: '
  end if
  link => self%head
  do while (associated(link))
    call link%info
    link => link%next
  end do
  call trace%end()
END SUBROUTINE info

!===============================================================================
!> Check refinement support.  When sign is +1 it checks if an overlapping level
!> L+2 task requires the level L link task to split to provide support.  When
!> sign is -1 it checks if the level L link task has support at level L-1.
!>
!> The use of a bitmap to reveal if a task lacks support is general; it does not
!> require stationary patches, or exact Cartesian overlap.  It only requires
!> that the source index regions computed are correct, so the bitmap is correct.
!===============================================================================
FUNCTION check_support (self, link, sign, irefine) RESULT (refine)
  class(list_t):: self
  class(link_t), pointer:: link
  integer, intent(in):: sign
  real, pointer:: irefine(:,:,:)
  integer:: refine
  !-----------------------------------------------------------------------------
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  logical:: guards, check
  integer:: l(3), u(3), n(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%check_support', itimer=itimer)
  task => link%task
  refine = 0
  !.............................................................................
  select type (task)
  class is (patch_t)
    !---------------------------------------------------------------------------
    ! If patch memory has not been allocated yet, return with refine=0
    !---------------------------------------------------------------------------
    if (.not. task%mem_allocated) then
      call trace%end (itimer)
      return
    end if
    !---------------------------------------------------------------------------
    ! Use an external test array when called from refine_t, and a local one
    ! when called stand-alone
    !---------------------------------------------------------------------------
    irefine = 0.0
    !---------------------------------------------------------------------------
    ! Find L+-2 tasks, and check for overlap of their guard zones
    !---------------------------------------------------------------------------
    check = .false.
    !call link%lock%set ('check_support', itimer)
    nbor => link%nbor
    do while (associated(nbor))
      if (nbor%task%level == link%task%level+sign*2) then
        nbtask => nbor%task
        select type (nbtask)
        class is (patch_t)
          if (sign > 0) then
            !-------------------------------------------------------------------
            ! For L+2 nbors, find the guard zone range that the level L task
            ! would provide, and set bitmap(:,:,:)=1 there
            !-------------------------------------------------------------------
            guards = .true.
            if (task%get_overlap (nbtask, guards, l, u)) then
              if (self%verbose > 1) then
                write (io_unit%log,'(f12.6,2(1x,a,i6,i4),2(2x,3i4))') &
                  wallclock(), &
                  'list_t%check_support: task, level =', &
                  task%id, task%level, &
                  'nbtask, level, needs  l, u =', &
                  nbtask%id, nbtask%level, l, u
              end if
              irefine(l(1):u(1),l(2):u(2),l(3):u(3)) = 2.0
              check = .true.
            end if
          else
            !-------------------------------------------------------------------
            ! Find the index range in the L+2 task that would be needed from the
            ! the level L task, not counting guard zones
            !-------------------------------------------------------------------
            guards = .false.
            if (task%get_overlap (nbtask, guards, l, u)) then
              if (self%verbose > 1) then
                write (io_unit%log,'(f12.6,2(1x,a,i6,i4),2(2x,3i4))') &
                  wallclock(), &
                  'list_t%check_support: task, level =', &
                  task%id, task%level, &
                  'nbtask, level, needs  l, u =', &
                  nbtask%id, nbtask%level, l, u
              end if
              irefine(l(1):u(1),l(2):u(2),l(3):u(3)) = 2.0
              check = .true.
            end if
          end if
        end select
      end if
      nbor => nbor%next
    end do
    !---------------------------------------------------------------------------
    ! If any candidates for missing support were found, run through L+-1 nbors,
    ! check for overlap of their interiors, and clear the refine flags in those
    ! cells; this cancels the flags possibly set above, if those L+-2 patches
    ! have or provide support at L+-1
    !---------------------------------------------------------------------------
    if (check) then
      nbor => link%nbor
      do while (associated(nbor))
        if (nbor%task%level == link%task%level+sign*1) then
          nbtask => nbor%task
          select type (nbtask)
          class is (patch_t)
            guards = .false.
            if (task%get_overlap (nbtask, guards, l, u)) then
              irefine(l(1):u(1),l(2):u(2),l(3):u(3)) = 0.0
              if (self%verbose > 1) then
                write (io_unit%log,'(f12.6,2(1x,a,i6,i4),2(2x,3i4),f5.1)') &
                  wallclock(), &
                  'list_t%check_support: task, level =', &
                  task%id, task%level, &
                  'nbtask, level, covers l, u, max =', &
                  nbtask%id, nbtask%level, l, u, &
                  maxval(irefine)
              end if
            end if
          end select
        end if
        nbor => nbor%next
      end do
      l = task%mesh%li
      u = task%mesh%ui
      refine = maxval(irefine(l(1):u(1),l(2):u(2),l(3):u(3)))
    end if
    !call link%lock%unset ('check_support')
    !---------------------------------------------------------------------------
    ! Task%irefine is used additively, so make sure refinement happen if and
    ! where the local irefine is still set    
    !---------------------------------------------------------------------------
    if (associated(task%irefine)) then
      task%irefine = task%irefine + irefine
    end if
  end select
  !-----------------------------------------------------------------------------
  ! If sign is negative then the task lack support if refine > 0 and otherwise
  ! has support.  The bits%no_support is set accordingly
  !-----------------------------------------------------------------------------
  if (sign < 0) then
    if (refine > 0) then
      call task%set (bits%no_support)
    else
      call task%clear (bits%no_support)
    end if
  end if
  !.............................................................................
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,3i4)') wallclock(), &
     'list_t%check_support: task, level, sign, support', &
     task%id, task%level, sign, refine
  end if
  !.............................................................................
  call trace%end (itimer)
END FUNCTION check_support

!===============================================================================
!> Find nbors to task with the given link, and add this task and the nbor task
!> to each others nbor lists, if they have a dependence.
!===============================================================================
SUBROUTINE find_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: scan
  logical:: needed, needs_me, download, done
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%find_nbors')
  scan => self%head
  do while (associated(scan))
    if (link%needs_to_add (scan)) then
      call link%add_nbor (scan)
      call scan%add_nbor (link)
    end if
    scan => scan%next
  end do
  call trace%end ()
END SUBROUTINE find_nbors

!===============================================================================
!> Check consistency of the nbor lists of all tasks in the list
!===============================================================================
SUBROUTINE check_nbor_consistency (self, ok)
  class(list_t):: self
  logical:: ok
  !.............................................................................
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%check_nbor_consistency')
  link => self%head
  do while (associated(link))
    if (link%nbor_inconsistency()) then
      write (stderr,*) 'ERROR: nbor list is inconsistent for task', link%task%id
    end if
    link => link%next
  end do
  call trace%end ()
END SUBROUTINE check_nbor_consistency

END MODULE list_mod
