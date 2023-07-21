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
  USE os_mod
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
  USE guards_mod
  USE halo_mod
  USE timer_mod
  USE running_mean_mod
  USE region_mod
  USE hash_table_mod
  USE shared_mod
  implicit none
  private
  integer, parameter:: max_level=40
  type linkp_t
    class(link_t), pointer:: pointer => null()
  end type
  type, public:: list_t
    character(len=64):: name=''
    class(link_t), pointer:: head => null()
    class(link_t), pointer:: tail => null()
    class(link_t), pointer:: queue => null()
    class(link_t), pointer:: active => null()
    class(list_t), pointer:: list => null()
    integer::  n=0                      ! total
    integer:: nq=0                      ! ready
    integer:: nv=0                      ! virtual
    integer:: nb=0                      ! boundary
    integer:: ni=0                      ! internal
    integer:: ne=0                      ! external
    integer:: no=0                      ! has output
    integer:: na=0                      ! active = internal+boundary
    integer:: nac=0                     ! active tasks
    integer:: nfinished=0               ! finished tasks
    integer:: max_id=0                  ! max MHD idy
    integer:: levelmin=9999             ! lowest level
    integer:: levelmax=0                ! highest level
    integer:: next_free=1               ! next free ID
    real(8):: lc(3) = +huge(0d0)        ! lower corner
    real(8):: uc(3) = -huge(0d0)        ! upper corner
    real(8):: llc(3) = +huge(0d0)       ! lower left corner
    real(8):: urc(3) = -huge(0d0)       ! upper right corner
    real(8):: position(3) = 0.0_8       ! center offset
    real(8):: size(3) = 1.0_8           ! full extent = uc-lc
    integer:: n_behind=0, n_check=0, n_ready=0, n_nbor=0
    integer:: debug=0
    integer:: hash_min=100
    integer:: n_tasks=1
    integer:: dims(3)=0
    integer:: verbose=0
    integer:: hash_table_size=128000
    type(lock_t):: lock, qlock
    logical:: use_hash=.true.
    logical:: face_nbors=.false.
    logical:: detailed_timer=.false.
    logical:: hands_off_nbor_lists=.false.
    character(len=4)              :: kind='list'
    type(mpi_window_t)            :: mpi_window
    type(running_mean_t), pointer :: running_mean
    type(hash_table_t)            :: hash_table
    type(linkp_t), dimension(max_level):: first
  contains
    procedure:: init
    procedure:: dealloc
    !procedure:: init_bdries
    procedure:: reset
    procedure:: prepare
    procedure:: add_task
    procedure:: add_link
    generic:: add => add_task, add_link
    procedure:: append_task
    procedure:: append_link
    generic:: append => append_task, append_link
    procedure:: prepend
    procedure:: remove
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
    procedure:: set_init_nbors
    procedure:: set_update_status
    procedure:: init_nbor_nbors
    procedure:: check_nbor_nbors
    procedure:: refresh_nbors
    procedure:: init_all_nbors
    procedure:: init_nonleaf
    procedure:: find_parent
    procedure:: find_parents
    procedure:: check_hash_table
    procedure:: remove_parents
    procedure:: intpos
    procedure:: consistency
    procedure:: qconsistency
    procedure:: aconsistency
    procedure:: statistics
    procedure:: append_list
    procedure:: prepend_list
    procedure:: add_by_level
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
    procedure:: update_status
    procedure:: swap
    procedure:: send_to_vnbors
    procedure:: test_nbor_status
    procedure:: info
    procedure:: find_nbors
    procedure:: check_nbor_consistency
    procedure:: check_duplicates
    procedure:: make_parent_levels
    procedure:: sort_by_level
    procedure:: init_task
    procedure:: get_head
  end type
  integer, parameter:: mnbor=100
  type(list_t), public:: list
  public link_t, task2patch
CONTAINS

!===============================================================================
!> Initialize list; this is kept lightweight to allow using temporary lists;
!> hash list, mulligan, and offload support are added in task_list_t
!===============================================================================
SUBROUTINE init (self, name)
  class(list_t), target:: self
  character(len=*), optional:: name
  !.............................................................................
  call trace%begin ('list_t%init')
  if (present(name)) then
    self%name = name
    if (self%verbose > 2) &
      write (stdlog,*) 'list_t%init: ', name
  else
    if (self%verbose > 2) &
      write (stdlog,*) 'list_t%init'
  end if
  call self%lock%init ('list')
  call self%qlock%init ('list')
  call self%hash_table%init (self%hash_table_size, 3)
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Deallocate related resources
!===============================================================================
SUBROUTINE dealloc (self)
  class(list_t):: self
  !.............................................................................
  call self%lock%destroy
  call self%qlock%destroy
  call self%hash_table%dealloc
END SUBROUTINE dealloc

!===============================================================================
!> Reset list, by deallocating memory used by tasks, and deallocating tasks and
!> links
!===============================================================================
SUBROUTINE reset (self, keep)
  class(list_t):: self
  logical, optional:: keep
  class(link_t), pointer:: link, old
  class(task_t), pointer:: task
  !.............................................................................
  call trace%begin ('list_t%reset')
  link => self%head
  do while (associated(link))
    task => link%task
    if (.not.present(keep)) then
      if (associated(task)) then
        call task%dealloc
        deallocate (task)
      end if
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
!> Allocate a new link, point to it, and back again
!===============================================================================
SUBROUTINE prepare (self, task)
  class(list_t):: self
  class(task_t), pointer:: task
  !.............................................................................
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%append_patch')
  allocate (halo_t:: link)
  call os%mem_register (storage_size(link), 1, 'list_%prepare')
  call link%init
  link%task => task
  select type (task)
  class is (patch_t)
    task%link => link
  end select
  call trace%end()
END SUBROUTINE prepare

!===============================================================================
!> Allocate a new link and append to a list, optionally setting attributes
!===============================================================================
SUBROUTINE append_task (self, task, nbor, needed, needs_me, keep)
  class(list_t):: self
  class(task_t), pointer:: task
  class(nbor_t), pointer, optional:: nbor
  logical, optional:: needed, needs_me, keep
  !.............................................................................
  class(link_t), pointer:: link
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%append_task', itimer=itimer)
  allocate (halo_t:: link)
  call os%mem_register (storage_size(link), 1, 'list_%append_task')
  call link%init
  !-----------------------------------------------------------------------------
  ! Optionally keep the task%link unchanged
  !-----------------------------------------------------------------------------
  link%task => task
  if (.not.present(keep)) then
    select type (task)
    class is (patch_t)
      task%link => link
    end select
  end if
  !-----------------------------------------------------------------------------
  ! Optionally set nbor properties
  !-----------------------------------------------------------------------------
  if (present(nbor))     link%nbor     => nbor
  if (present(needed))   nbor%needed   = needed
  if (present(needs_me)) nbor%needs_me = needs_me
  !-----------------------------------------------------------------------------
  call self%append_link (link, keep)
  call trace%end (itimer)
END SUBROUTINE append_task

!===============================================================================
!> Light weight append to a list
!===============================================================================
SUBROUTINE append_link (self, link, keep)
  class(list_t):: self
  class(link_t), pointer:: link
  logical, optional:: keep
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('list_t%append_link', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Append in such a way that the list in principle does not need to be locked
  !-----------------------------------------------------------------------------
  nullify(link%next)
  call                                      self%lock%set ('list_t%append_link')
  if (associated(self%tail)) then
    link%prev => self%tail
    self%tail%next => link
    self%tail => link
  else
    self%head => link
    self%tail => link
  end if
  self%max_id = max(self%max_id,link%task%id)
  self%n = self%n + 1
  call                                    self%lock%unset ('list_t%append_link')
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE append_link

!===============================================================================
!> Light weight prepend to a list
!===============================================================================
SUBROUTINE prepend (self, link, keep)
  class(list_t):: self
  class(link_t), pointer:: link
  logical, optional:: keep
  !.............................................................................
  call trace%begin ('list_t%prepend')
  !-----------------------------------------------------------------------------
  ! Prepend in such a way that the list in principle does not need to be locked
  !-----------------------------------------------------------------------------
  nullify(link%next)
  call                                          self%lock%set ('list_t%prepend')
  if (associated(self%head)) then
    link%next => self%head
    self%head%prev => link
    self%head => link
  else
    self%head => link
    self%tail => link
  end if
  self%max_id = max(self%max_id,link%task%id)
  self%n = self%n + 1
  call                                        self%lock%unset ('list_t%prepend')
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE prepend

!===============================================================================
!> Allocate a new link and append to the list, optionally setting attributes
!===============================================================================
SUBROUTINE add_task (self, task, nbor, needed, needs_me, keep)
  class(list_t):: self
  class(nbor_t), pointer, optional:: nbor
  logical, optional:: needed, needs_me, keep
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  !.............................................................................
  call trace%begin ('list_t%add_task')
  allocate (halo_t:: link)
  call os%mem_register (storage_size(link), 1, 'list_%add_task')
  call link%init
  link%task => task
  if (.not.present(keep)) then
    select type (task)
    class is (patch_t)
      task%link => link
    end select
  end if
  !-----------------------------------------------------------------------------
  ! Note that adding a task to more than one list requires care; the task%link
  ! must remain a pointer to the link that has the nbor list, and the nbor list
  ! head must be accesses as task%link%nbor, not as link%nbor, for links in the
  ! 2nd task list.
  !-----------------------------------------------------------------------------
  if (present(nbor))     link%nbor     => nbor
  if (present(needed))   nbor%needed   = needed
  if (present(needs_me)) nbor%needs_me = needs_me
  call self%add_by_level (link, keep)
  call trace%end()
END SUBROUTINE add_task

!===============================================================================
!> Append to the task list, counting active tasks only
!===============================================================================
SUBROUTINE add_link (self, link, keep)
  class(list_t):: self
  class(link_t), pointer:: link
  logical, optional:: keep
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%add_link')
  call self%add_by_level (link, keep)
  call trace%end()
END SUBROUTINE add_link

!===============================================================================
!> Remove a task from consideration, making sure that virtual copies do the same.
!> Set bits%remove, which excludes the task from nbor lists in check_ready(),
!> and signals rank nbors to repeat this remove_and_reset().
!===============================================================================
SUBROUTINE remove_and_reset (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%remove_and_reset')
  call io%abort ('list_t%remove_and_reset: unsuccessful overload')
  call trace%end()
END SUBROUTINE remove_and_reset

!===============================================================================
!> Remove a link from the list, taking special care with the head and tail
!===============================================================================
SUBROUTINE remove (self, link, count)
  class(list_t):: self
  class(link_t), pointer:: link
  logical, optional:: count
  !.............................................................................
  logical:: done
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%remove')
  !-----------------------------------------------------------------------------
  ! Check for head task
  !-----------------------------------------------------------------------------
  if (associated(link,self%head)) then
    self%head => link%next
    if (associated(self%head)) then
      nullify (self%head%prev)
    end if
    done = .true.
  end if
  !-----------------------------------------------------------------------------
  ! Check for tail task
  !-----------------------------------------------------------------------------
  if (associated(link,self%tail)) then
    self%tail => link%prev
    if (associated(self%tail)) then
      nullify (self%tail%next)
    end if
    done = .true.
  end if
  !-----------------------------------------------------------------------------
  ! If none of the above applied, the link can be removed by cross-over
  !-----------------------------------------------------------------------------
  if (.not.done) then
    link%prev%next => link%next
    link%next%prev => link%prev
  end if
  !-----------------------------------------------------------------------------
  ! If this was the first task on the level, reassign
  !-----------------------------------------------------------------------------
  if (associated (link, self%first(link%task%level)%pointer)) then
    self%first(link%task%level)%pointer => link%next
    if (self%verbose > 0) &
      write (stdlog,*) 'self%first(link%task%level)%pointer => next'
  end if
  call                                    self%lock%unset ('list_t%remove_link')
  !-----------------------------------------------------------------------------
  ! Remove the task key from the hash table
  !-----------------------------------------------------------------------------
  call self%hash_table%remove (link%task%key)
  !-----------------------------------------------------------------------------
  ! IMPORTANT: Some procedures depend on knowing the size of the list!
  !-----------------------------------------------------------------------------
  !call self%update_counts (link, -1)
  if (shared%option(count)) then
    call self%count_status ()
  end if
  call trace%end ()
END SUBROUTINE remove

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
                                call self%lock%set_read ('list_t%update_counts')
  if (link%task%is_set (bits%ready)) then
    self%nq = self%nq + delta
  end if
  if      (link%task%is_set (bits%virtual)) then
    self%nv = self%nv + delta
  else if (link%task%is_set (bits%boundary)) then
    self%nb = self%nb + delta
    self%na = self%na + delta
  else if (link%task%is_set (bits%internal)) then
    self%ni = self%ni + delta
    self%na = self%na + delta
  end if
  if (link%task%is_set (bits%has_output) .and. &
      link%task%is_clear(bits%virtual)) then
    self%no = self%no + delta
  end if
  !$omp atomic update
  self%n = self%n + delta
                              call self%lock%unset_read ('list_t%update_counts')
  call trace%end ()
END SUBROUTINE update_counts

!===============================================================================
!> Set the io% attributes needed to keep track of I/O progress
!===============================================================================
SUBROUTINE set_io_counts (self)
  class (list_t):: self
  !-----------------------------------------------------------------------------
  ! Make sure task counts for I/O purposes are up-to-date
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%set_io_counts')
  call self%count_status ('set_io_counts')
  associate (io_ntask=>io%ntask, io_nwrite=>io%nwrite)
  !$omp atomic write
  io_ntask = self%na
  !$omp atomic write
  io_nwrite = self%no
  if (self%verbose > 0) &
    write (stdlog,*) 'set_io_counts: io_ntask,io_nwrite=',io_ntask,io_nwrite
  end associate
  call trace%end ()
END SUBROUTINE set_io_counts

!===============================================================================
!> Count the number of internal, boundary, virtual, frozen, external, and
!> output tasks
!===============================================================================
SUBROUTINE count_status (self, label)
  class (list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link
  class(nbor_t) , pointer:: nbor
  class(task_t), pointer:: task
  integer:: n, ni, nb, nv, ne, nf, no, nq, i, level
  logical:: debug
  integer, save:: itimer=0
  !.............................................................................
  call trace%begin ('list_t%count_status', itimer=itimer)
  call                              self%lock%set_read  ('count_status', itimer)
  nb=0; ni=0; nv=0; ne=0; nf=0; n=0; no=0; nq=0
  link => self%head
  debug = .false.
  do i=1,2
  if (debug) then
    write (stdout,*) 'first(level):', self%n
    do level=shared%levelmin,shared%levelmax
      if (associated(self%first(level)%pointer)) then
        write (stdout,*) 'level:', level, self%first(level)%pointer%task%id
        flush (stdout)
      end if
    end do
  end if
  n = 0
  do while (associated(link))
    task => link%task
    if (debug) then
      write (stdlog,*) 'id:', task%id, task%level
      flush (stdlog)
    end if
    n = n+1
    if (n > self%n+1) then
      if (debug) then
        call task%log ('inconsistent(2)')
        write (stdout,*) 'n:', n, self%n
        call io%abort ('list_t%count_status: inconsistent list')
      else
        call task%log ('inconsistent(1)')
        write (stdout,*) 'n:', n, self%n
        call io%warning ('list_t%count_status: inconsistent list')
      end if
      exit
    end if
    if (task%is_set(bits%frozen  )) nf = nf+1
    if (task%is_set(bits%boundary)) nb = nb+1
    if (task%is_set(bits%virtual )) nv = nv+1
    if (task%is_set(bits%internal)) ni = ni+1
    if (task%is_set(bits%external)) ne = ne+1
    if (task%is_set(bits%ready   )) nq = nq+1
    if (task%is_set(bits%has_output) .and. &
        task%is_clear(bits%virtual)) no = no+1
    if (task%is_clear (bits%boundary + &
                       bits%virtual  + &
                       bits%internal + &
                       bits%external)) then
      write (stderr,'(a,i6,3x,"ibvef =",1x,5l1,4x,z8.8)') &
        'list_t%count_status: WARNING, status bit is missing in id, IBVEF =', &
        task%id, &
        task%is_set(bits%internal), &
        task%is_set(bits%boundary), &
        task%is_set(bits%virtual ), &
        task%is_set(bits%external), &
        task%is_set(bits%frozen  ), task%status
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
  if (n <= self%n+1) then
    exit
  else
    debug = .true.
  end if
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
  call                                     self%lock%unset_read ('count_status')
  !-----------------------------------------------------------------------------
  ! Set the counters that affect I/O decisions
  !-----------------------------------------------------------------------------
  !$omp atomic write
  io%nwrite = no
  !$omp atomic write
  io%ntask = self%na
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE count_status

!===============================================================================
!> Clear, set, and count the status bits
!===============================================================================
SUBROUTINE reset_status (self, check)
  class(list_t):: self
  logical, optional:: check
  class(link_t), pointer:: link
  integer:: n_internal, n_boundary, n_virtual, n_external, n_frozen, n, id
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
  !-----------------------------------------------------------------------------
  ! Clear and set status bits
  !-----------------------------------------------------------------------------
  call                                            self%lock%set ('reset_status')
  n = 0
  link => self%head
  do while (associated(link))
    id = link%task%id
    call link%reset_status (link%nbor)
    n_internal = n_internal + merge(1,0,link%task%is_set(bits%internal))
    n_boundary = n_boundary + merge(1,0,link%task%is_set(bits%boundary))
    n_virtual  = n_virtual  + merge(1,0,link%task%is_set(bits%virtual ))
    n_external = n_external + merge(1,0,link%task%is_set(bits%external))
    n_frozen   = n_frozen   + merge(1,0,link%task%is_set(bits%frozen  ))
    n = n + 1
    link => link%next
  end do
  call                                          self%lock%unset ('reset_status')
  !-----------------------------------------------------------------------------
  ! Optionally, check consistency
  !-----------------------------------------------------------------------------
  !$omp atomic write
  self%na = n_internal+n_boundary-n_frozen
  !$omp atomic write
  self%n  = n
  call self%count_status ()
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
  class(task_t), pointer:: task, nbtask
  integer:: itr, status, ignore
  character(24):: label
  logical:: needs_me
  integer, save:: itimer=0
  !.............................................................................
  task => link%task                                     ! main task
  itr = 1
  call trace%begin('list_t%check_nbors', itimer=itimer)
  call link%task%log ('check_nbors', 4)
  if (self%verbose > 0) then
    write (io_unit%log,'(f12.6,1x,a,i6,i4,f12.6,2x,4l1)') wallclock(), &
    'list_t%check_nbors: task', link%task%id, link%task%level, link%task%time, &
     link%task%is_set(bits%boundary), link%task%is_set(bits%virtual), &
     associated(link%nbor)
  end if
  !-----------------------------------------------------------------------------
  ! Lock the other task link while accessing its nbor list, to prevent changes
  !-----------------------------------------------------------------------------
  ignore = bits%virtual + bits%external + bits%ready + bits%busy
  nbor => link%nbor                                     ! first nbor
  do while (associated (nbor))                          ! keep going until end
    nbtask => nbor%task
    !---------------------------------------------------------------------------
    ! If the nbor needs this task, and isn't busy updating, check it.  For
    ! thread safety, the busy bit must NOT be cleared before a task has done
    ! a check_ready() on itself, lest an opportunity could be missed.
    !---------------------------------------------------------------------------
    !$omp atomic read
    needs_me = nbor%needs_me
    if (nbtask%is_clear(ignore) .and. needs_me) then
      ! call                                   self%qlock%set_read ('check_nbors')
      call self%check_ready (nbor%link, lock=.true.)
      !.........................................................................
      if (self%verbose > 1) &
        write (io_unit%log,'(f12.6,i6,1x,a,i6,i4)') wallclock(), task%id, &
        'checked nbor ', nbtask%id, nbtask%level
      !.........................................................................
      ! call                                 self%qlock%unset_read ('check_nbors')
    else
      !.........................................................................
      if (io%task_logging > 3 .or. self%verbose > 1) then
        write (label,'(a,6x,5a1)') 'skip', &
          merge('N','n',needs_me), &
          merge('R','r',nbtask%is_set(bits%ready  )), &
          merge('B','b',nbtask%is_set(bits%busy   )), &
          merge('V','v',nbtask%is_set(bits%virtual))
        call nbtask%log (label, 4)
      end if
      if (self%verbose > 1) &
        write (io_unit%log,'(f12.6,i6,2x,a,i6)') wallclock(), task%id, &
        label, nbtask%id
      !.........................................................................
    end if
    nbor => nbor%next                                   ! next nbor
  end do
  call trace%end (itimer)
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
  integer:: nq
  logical:: is_ahead
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%check_all', itimer=itimer)
  write (io_unit%log,*) wallclock(), 'check_all: phase 1', &
    self%nq, associated(self%queue)
  if (present(repair)) then
    nullify (self%queue)
  end if
  !-----------------------------------------------------------------------------
  ! Phase 1: check all tasks with a normal check_ready()
  !-----------------------------------------------------------------------------
  nq = self%nq
  call                                              self%qlock%set ('check_all')
  link => self%head
  call link%task%log ('check_all(1)')
  write (stderr,'(a,i5,1x,a,i5)') &
    'WARNING: check_all phase 1, rank', mpi%rank, 'thread', omp%thread
  write (stdlog,*) 'check_all: phase 1, check_ready()'
  do while (associated (link))                          ! keep going until end
    if (link%task%is_clear (bits%virtual+bits%external)) then
      call self%check_ready (link, lock=.true.)         ! link ready?
      if (self%nq > nq) then
        write (io_unit%log,*) 'check_all found', link%task%id, link%task%time
        nbor => link%nbor
        do while (associated(nbor))
          is_ahead = nbor%task%is_ahead_of (link%task)
          write (io_unit%log,'(a,i6,3x,4l1,2f8.3)') &
            'nbor,NND,ahead,delta/dtime,time=', nbor%task%id, &
            nbor%needed, nbor%needs_me, nbor%download, is_ahead, &
            (nbor%task%time - link%task%time)/max(tiny(1d0),nbor%task%dtime), &
            nbor%task%grace + shared%extra_grace
          nbor => nbor%next
        end do
        nq = self%nq
      end if
    end if
    link => link%next
  end do
  call                                            self%qlock%unset ('check_all')
  !-----------------------------------------------------------------------------
  ! If the queue has tasks, return
  !-----------------------------------------------------------------------------
  if (associated(self%queue)) then
    write (stdlog,*) wallclock(), 'check_all: done, nq =', &
      self%nq, associated(self%queue)
  !-----------------------------------------------------------------------------
  ! Phase 2: clearing ready bits, if set and preventing queing
  !-----------------------------------------------------------------------------
  else if (shared%extra_grace == 0.0) then
    write (stderr,'(a,i5,1x,a,i5)') &
      'WARNING: check_all phase 2, rank', mpi%rank, 'thread', omp%thread
    write (stdlog,*) 'check_all: phase 2, clearing ready bits'
    call                                          self%lock%set ('check_all(2)')
    link => self%head
    call link%task%log ('check_all(2)')
    do while (associated (link))                        ! keep going until end
      if (link%task%is_clear (bits%virtual+bits%external)) then
        if (link%task%is_set (bits%ready)) then
          write (io_unit%log,*) 'clearing ready bit on', link%task%id
          call link%task%clear (bits%ready)
        end if
        call self%check_ready (link, lock=.true.)
      end if
      link => link%next
    end do
    call                                        self%lock%unset ('check_all(2)')
  end if
  write (io_unit%log,*) wallclock(), 'check_all: done, nq =', &
    self%nq, associated(self%queue), shared%extra_grace
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
  call                                        self%lock%set_read ('resend_bdry')
  link => self%head
  do while (associated(link))
    if (link%task%is_set (bits%boundary)) then
      call self%send_to_vnbors (link)
    end if
    link => link%next
  end do
  call                                      self%lock%unset_read ('resend_bdry')
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
  integer:: it, istep
  logical:: is_ahead
  integer, save:: phase=1
  !.............................................................................
  call trace%begin('list_t%check_oldest')
  write (io_unit%log,*) 'check_oldest: phase 1', wallclock()
  flush (io_unit%log)
  told = 1d30
  toldv = told
  nullify (oldest)
  nullify (oldestv)
  nullify (oldestnb)
  call                                       self%lock%set_read ('check_oldest')
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
  call                                     self%lock%unset_read ('check_oldest')
  if (associated(oldest)) then
    call oldest%task%log ('oldest')
    write (io_unit%log,'(f12.6,2x,a,i6,f12.6," rb:",2l1)') wallclock(), &
      'the oldest active task is ', oldest%task%id , told, &
      oldest%task%is_set(bits%ready), oldest%task%is_set(bits%busy)
  else
    write (stderr,*) 'STALLED: oldest active not found '
  end if
  flush (io_unit%log)
  if (associated(oldestv)) then
    call oldestv%task%log ('oldestv')
    write (io_unit%log,'(f12.6,2x,a,i6,f12.6)') wallclock(), &
      'the oldest virtual task is ', oldestv%task%id, toldv
    write (io_unit%log,'(f12.6,2x,a,f12.6)') wallclock(), &
      'it was last unpacked at', oldestv%task%unpack_time
    flush (io_unit%log)
  end if
  if (associated(oldest)) then
    nbor => oldest%nbor
    it = max(1,oldest%task%nt-1)
    it = oldest%task%iit(it)
    write (io_unit%log,2) oldest%task%id, oldest%task%istep, &
      oldest%task%time, oldest%task%t(it), &
      oldest%task%rank, oldest%task%level, &
      oldest%task%is_set(bits%boundary), oldest%task%is_set(bits%virtual)
  2 format(2i6,2f12.6,2i4,5L3,f12.6)
    write (io_unit%log,*) 'nbor list:'
    do while (associated(nbor))
      it = max(1, nbor%task%nt-1)
      it = nbor%task%iit(it)
      is_ahead = nbor%task%is_ahead_of (oldest%task) 
      write (io_unit%log,2) nbor%task%id, nbor%task%istep, &
        nbor%task%time, nbor%task%t(it), &
        nbor%task%rank, nbor%task%level, nbor%needed, nbor%needs_me, &
        nbor%task%is_set(bits%boundary), nbor%task%is_set(bits%virtual), &
        is_ahead, nbor%task%unpack_time
      nbor => nbor%next
    end do
    if (present(clear)) then
      call oldest%task%clear (bits%ready+bits%busy)
      !$omp atomic read
      istep = oldest%task%istep
      call self%queue_by_time (oldest, istep)
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
SUBROUTINE check_ready (self, link, lock)
  class(list_t):: self
  class(link_t), pointer:: link
  logical, optional:: lock
  !.............................................................................
  logical:: ok, is_ahead
  class(task_t), pointer:: task, nbtask
  class(link_t), pointer:: oldest
  class(nbor_t), pointer:: nbor
  real(8):: t_oldest
  character(24):: label
  logical:: lock_l
  integer:: ignore, unit, istep
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Since we may be called upon to check nbors of a patch that should not itself
  ! be checked, we first check if any bits are set indicating this to be the case
  !-----------------------------------------------------------------------------
  if (self%detailed_timer .or. self%verbose > 0) &
    call trace%begin ('list_t%check_ready', itimer=itimer)
  call link%task%log ('check_ready', 4)
  lock_l = shared%option (lock)
  task => link%task
  !$omp atomic read
  istep = task%istep
  !-----------------------------------------------------------------------------
  ! Ignore tasks with bits set that indicate they are not candidates for
  ! scheduling.  Tasks with bits%ready are in the queue, and tasks with
  ! bits%busy are updating, with bits%busy cleared by the current thread
  ! before calling this procedure, so no chance to miss a ready state.
  !-----------------------------------------------------------------------------
  ignore = bits%ready + bits%busy + bits%remove + bits%virtual + bits%frozen &
         + bits%external
  ok = task%is_clear (ignore)
  if (ok) then
    !...........................................................................
    if (self%verbose > 2 .or. task%id == io%id_debug) then
      write (stdlog,'(f12.6,1x,a,i6,2l3)') wallclock(), &
       'list_t%check_ready: task, lock, nbors', link%task%id, present(lock), &
       associated(link%nbor)
    end if
    !---------------------------------------------------------------------------
    ! Since only this thread can change the nbor list we are free to use it w/o
    ! lock in cases where task is the active task, but in other cases we need to
    ! lock the nbor list, so it isn't changed during the loop. Even when the
    ! current thread owns the task we need to prevent other tasks to start
    ! other checks, which could even start before the lock and end after it,
    ! unless the lock is a write-lock.
    !---------------------------------------------------------------------------
    if (lock_l) call                               link%lock%set ('check_ready')
    !---------------------------------------------------------------------------
    ! Tasks that are frozen or removed should not be used in the loop
    !---------------------------------------------------------------------------
    ignore = bits%frozen + bits%remove + bits%external
    ok = .true.
    t_oldest = huge(1d0)
    nbor => link%nbor                                   ! use original nbor list
    do while (associated (nbor))                        ! keep going until end
      nbtask => nbor%task
      nbor%check_step = nbtask%istep                    ! diagnostics
      !-------------------------------------------------------------------------
      ! Check nbors that are needed, not ignored, and have support
      !-------------------------------------------------------------------------
      if (nbor%needed .and. nbor%task%is_clear(ignore)) then
         !-----------------------------------------------------------------------
        ! This extra test is only needed if task_t%nbor_relations does not
        ! set nbor%needed to .false.
        !-----------------------------------------------------------------------
        if (abs(task%level-nbtask%level) <= 1) then
          if (.not. nbor%task%is_ahead_of (task)) then
            ok = .false.
            call nbtask%log ('not ahead', 4)
          end if
        end if
        !.......................................................................
        if (self%verbose > 1 .or. task%id == io%id_debug) then
          is_ahead = nbor%task%is_ahead_of(task)
          write (stdlog,'(f12.6,1x,2(a,i6,":",i2.2,f10.6,2x),1pg12.3,2x,3l1)') &
            wallclock(), &
            'list_t%check_ready: task,t=', task%id, task%level, task%time, &
            'nbtask,t,delta,ahead =', nbtask%id, nbtask%level, nbtask%time, &
            (nbtask%time-task%time)/max(nbtask%dt(nbtask%it),tiny(1d0)), &
            is_ahead, nbor%needed, nbor%needs_me
        end if
        !.......................................................................
      else
        !.......................................................................
        if (io%task_logging > 3) then
          write (label,'(a,4x,3a1)') 'nbignore', &
            merge('N','n',nbor%needed), &
            merge('D','d',nbtask%is_set(bits%remove )), &
            merge('F','f',nbtask%is_set(bits%frozen ))
          call nbtask%log (label, 4)
        end if
        if (self%verbose > 1 .or. task%id == io%id_debug) &
          write (stdlog,'(f12.6,1x,a,i6,1x,a,i6,i4,2x,4l1)') &
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
    if (lock_l) call                             link%lock%unset ('check_ready')
    if (ok) then
      !$omp atomic write
      task%n_failed = 0
      associate (check_step => task%check_step)
      !$omp atomic write
      check_step = task%istep
      end associate
      call task%log ('ready', 4)
      !.........................................................................
      if (self%verbose > 0 .or. task%id == io%id_debug) &
        write (stdlog,'(f12.6,1x,a,i6,f12.6,1x,a)') wallclock(), &
        'list_t%check_ready: queue_by_time', task%id, task%time, trim(task%kind)
      !.........................................................................
      call self%queue_by_time (link, istep)
    else
      !$omp atomic update
      task%n_failed = task%n_failed+1
      !-------------------------------------------------------------------------
      ! At verbose=1 we get a line in the rank log for every 10000 failures,
      ! while at verbose=2 we get a line for ever fail, and line on stderr for
      ! every 10000 failures
      !-------------------------------------------------------------------------
      if (self%verbose > 1) then
        if (associated(oldest)) then
          if (io%task_logging > 2) then
            write (label,'(a,2i8)') &
              'fail', oldest%task%id, oldest%task%level
            call task%log (label, 3)
          end if
          write (stdlog,'(f12.6,i6,l3,f12.6,1x,a,i6,f12.6,i5,l3)') &
            wallclock(), task%id, task%is_set(bits%boundary), task%time, &
            'list_t%check_ready: failed with id,time,rank,V=', &
            oldest%task%id, t_oldest, oldest%task%rank, &
            oldest%task%is_set(bits%virtual)
        else
          if (io%task_logging > 2) then
            write (label,'(a,i8,f12.6)') &
              'fail', nbtask%id, nbtask%time
            call task%log (label, 3)
          end if
          write (stdlog,'(f12.6,i6,l3,f12.6,1x,a,i6,f12.6,i5,l3)') &
            wallclock(), task%id, task%is_set(bits%boundary), task%time, &
            'list_t%check_ready: failed with id,time,rank,V=', &
            nbtask%id, nbtask%time, nbtask%rank, &
            nbtask%is_set(bits%virtual)
        end if
      end if
      !.........................................................................
      if (task%n_failed > 10000) then
        if (self%verbose > 1) then
          unit = stderr
        else
          unit = stdlog
        end if
        write (unit ,'(a,i4,1p,2(2x,a,i6,1x,a,g14.5),2(2x,a,g11.3))') &
          'list_t%check_ready: on rank', task%rank, &
          'task'               , task%id     , 'at t=', task%time, &
          'failed on nbor task', nbor%task%id, 'at t=', nbor%task%time, &
          'grace =', task%grace, ' dt=', nbor%task%dtime
      end if
      !.........................................................................
    end if
  else 
    !...........................................................................
    if (io%task_logging > 3) then
      write (label,'(a,4x,5a1)') 'ignore', &
        merge('R','r',task%is_set(bits%ready  )), &
        merge('B','b',task%is_set(bits%busy   )), &
        merge('D','d',task%is_set(bits%remove )), &
        merge('V','v',task%is_set(bits%virtual)), &
        merge('F','f',task%is_set(bits%frozen ))
      call task%log (label, 4)
    end if
    if (self%verbose > 1 .or. task%id == io%id_debug) then
      write (stdlog,'(f12.6,1x,a,i6,1x,a,5l1)') wallclock(), &
        'list_t%check_ready:', task%id, 'return because VFRr = ', &
        task%is_set (bits%virtual), &
        task%is_set (bits%frozen ), &
        task%is_set (bits%remove ), &
        task%is_set (bits%ready  )
    end if
    !...........................................................................
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
  call                                 self%lock%set_read ('list_t%consistency')
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
  call                               self%lock%unset_read ('list_t%consistency')
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
  call trace%begin ('list_t%qconsistency',4)
  n = 0
  if (io%verbose > 2) then
    link => self%head
    do while (associated(link))
      print *,'task list: id,time=', link%task%id, link%task%time
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
    if (io%verbose > 2) then
      print *,'queue: n,id,time,rb=', n, link%task%id, link%task%time, &
        link%task%is_set (bits%ready), link%task%is_set (bits%busy)
    end if
    if (task%time < time) then
      write (io_unit%log,'(a,i4,2(i9,1p,g15.6,3x))') &
        'ERROR: queue out of order ident,prev%id,time,task%id=', &
        ident, prev%id, time, task%id, task%time
    end if
    prev => task
    link => link%next_time
  end do
  if (n /= self%nq) then
    write (io_unit%log,*) &
      'ERROR: qconsistency ident,n,nq,ident=', ident, n, self%nq
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
  if (self%debug < 2) return
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
  call trace%begin('list_t%init_nonleaf')
  call                                     self%lock%set ('list_t%init_nonleaf')
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
  call                                   self%lock%unset ('list_t%init_nonleaf')
  call trace%end()
END SUBROUTINE init_nonleaf

!===============================================================================
!> Add a new task (e.g. originating from AMR or load balancing). This procedure
!> contains the actions that need to be performed on both the send and receive
!> side in a multi-process setting.
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
  ! Initialize task attributes from the task list class; this is important eg
  ! to ensure correctness of hash tags
  !-----------------------------------------------------------------------------
  call self%init_task (task)
  !-----------------------------------------------------------------------------
  ! Increment task counters on the rank, unless the task is virtual
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%virtual)) then
    !$omp atomic update
    timer%levelpop (task%level) = &         ! increment level population count
    timer%levelpop (task%level) + 1
    !$omp atomic update
    self%n = self%n + 1
    !$omp atomic update
    self%na = self%na + 1
    call self%set_io_counts ()
  end if
  !-----------------------------------------------------------------------------
  ! Generate a new nbor list
  !-----------------------------------------------------------------------------
  call self%init_nbors (link)               ! add nbor list
  call link%check_status ('add_new_link')   ! detect changes
  call self%find_parent (link)              ! connect to parent
  call self%set_init_nbors (link)           ! trigger new nbor lists
  call self%check_nbors (link)              ! check new nbors
  if (io_unit%do_validate) &
    call link%print_nbors()                 ! compact nbor list to stdout
  !.............................................................................
  if (self%verbose > 0) &
    call link%info ('add_new_link')
  !-----------------------------------------------------------------------------
  ! If it is a boundary task, send it to rank nbors
  !-----------------------------------------------------------------------------
  if (task%is_set (bits%boundary)) &
    call self%send_to_vnbors (link)
  call trace%end ()
END SUBROUTINE add_new_link

!===============================================================================
!> Use the fast hash table method to find +-1 level nbors, and expand to +-2
!> level nbors by searching for overlap among nbors of nbors. Then connect the
!> new nbor list under a brief lock and deallocate the previous one.
!===============================================================================
SUBROUTINE init_nbors (self, link)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(nbor_t), pointer:: new_nbors, old_nbors, new_nbors2
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%init_nbors', itimer=itimer)
  call link%task%log ('init_nbors', 2)
  nullify (new_nbors)
  !-----------------------------------------------------------------------------
  ! Optionally, check that the hash table method gives consistent results
  !-----------------------------------------------------------------------------
  if (self%debug == 1) then
    call init_nbors_by_level (self, link, new_nbors)
    call init_nbors_hash     (self, link, new_nbors2)
    call compare_nbors (new_nbors, new_nbors2, 'hash nbors')
    call compare_nbors (new_nbors2, new_nbors, 'by_level nbors')
    call link%remove_nbor_list (new_nbors2)
  !-----------------------------------------------------------------------------
  ! If the task list is long, use the hash table to speed up the search,
  ! otherwise use a direct search of the task list
  !-----------------------------------------------------------------------------
  else if (self%n > self%hash_min .and. self%use_hash) then
    call init_nbors_hash   (self, link, new_nbors)
  else
    call init_nbors_by_level (self, link, new_nbors)
  end if
  !-----------------------------------------------------------------------------
  ! Swap the link%nbor pointer under lock
  !-----------------------------------------------------------------------------
  call                                              link%lock%set ('init_nbors')
  old_nbors => link%nbor
  link%nbor => new_nbors
  !-----------------------------------------------------------------------------
  ! Update the status bits based on the new nbor list and adjust counts
  !-----------------------------------------------------------------------------
  call link%reset_status (link%nbor)
  call self%update_status (link)
  call                                            link%lock%unset ('init_nbors')
  !-----------------------------------------------------------------------------
  ! Clear the init_nbors bit and remove the old nbor list
  !-----------------------------------------------------------------------------
  call link%task%clear (bits%init_nbors)
  call link%remove_nbor_list (old_nbors)
  !.............................................................................
  if (self%verbose > 0) &
    call link%info ('init_nbors')
  call trace%end (itimer)
contains
  !=============================================================================
  !> Check if all nbors1 tasks are present in nbors2
  !=============================================================================
  subroutine compare_nbors (nbors1, nbors2, label)
    class(nbor_t), pointer:: nbors1, nbors2, nbor1, nbor2
    logical:: ok
    character(*):: label
    !---------------------------------------------------------------------------
    nbor1 => nbors1
    do while (associated(nbor1))
      ok = .false.
      nbor2 => nbors2
      do while (associated(nbor2))
        if (nbor2%task%id == nbor1%task%id) then
          ok = .true.
          exit
        end if
        nbor2 => nbor2%next
      end do
      if (.not. ok) then
        write (stderr,'(2(i7,":",i2.2,3i4,2x,a))') &
          link%task%id, link%task%key, 'compare_nbors: nbor ', &
          nbor1%task%id, nbor1%task%key, ' missing in '//label
        write (io_unit%mpi,'(2(i7,":",i2.2,3i4,2x,a))') &
          link%task%id, link%task%key, 'compare_nbors: nbor ', &
          nbor1%task%id, nbor1%task%key, ' missing in '//label
        call init_nbors_hash (self, link, new_nbors2)
      end if
      nbor1 => nbor1%next
    end do
  end subroutine compare_nbors
END SUBROUTINE init_nbors


!===============================================================================
!> Update the nbor list of task, by searching among all tasks within +-1 of
!> link%task%level
!>
!> With gfortran, this is about twice as fast as init_nbors, which searches
!> nbors of nbors, using a hash table to avoid trying the same repeatedly.
!> The balance of cost could change if there are many more tasks per level
!> per core.
!===============================================================================
SUBROUTINE init_nbors_by_level (self, link, nbors)
  class(list_t)          :: self
  class(link_t)          :: link
  class(nbor_t), pointer :: nbors
  !.............................................................................
  class(link_t), pointer :: scan
  class(task_t), pointer :: task
  class(nbor_t), pointer :: nbor, nbor2
  integer                :: n, level, levelmin, levelmax, tid, sid
  real(8)                :: start
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  ! Save a pointer to the current nbor list, to remove after adding a new
  ! (to ensure that task%n_needed does not touch zero prematurely)
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%init_nbors_by_level', itimer=itimer)
  nullify (nbors)
  !-----------------------------------------------------------------------------
  ! Run through (only) tasks at the relevant levels
  !-----------------------------------------------------------------------------
  start = wallclock()
  n = 0
  task => link%task
  tid = task%id
  !.............................................................................
  if (self%verbose == -1) then
    select type(task)
    class is (patch_t)
    write (stdlog,'(f12.6,a,i6,":",i2.2,2x,4i4)') wallclock(), &
      pad_trace(1:indent)//'init_nbors: id,key=', &
      task%id, task%level, task%hash_key()
    end select
  end if
  !.............................................................................
  if (task%level > 0) then
    levelmin = max(task%level-1, shared%levelroot)
    levelmax = min(task%level+1, shared%levelmax)
  else
    levelmin = shared%levelroot
    levelmax = shared%levelmax
  end if
  do level=levelmin,levelmax
    if (associated(self%first(level)%pointer)) then
      call                            self%lock%set_read ('init_nbors_by_level')
      scan => self%first(level)%pointer
      do while (associated(scan))
        sid = scan%task%id
        if (scan%task%level /= level) then
          exit
        end if
        if (scan%task%id /= task%id) then
          select type (link)
          class is (halo_t)
            call link%conditionally_add_nbor (scan, nbors)
          end select
          n = n+1
        end if
        scan => scan%next
      end do
      call                          self%lock%unset_read ('init_nbors_by_level')
    end if
  end do
  call trace%end (itimer)
END SUBROUTINE init_nbors_by_level

!===============================================================================
SUBROUTINE init_nbors_direct (self, link, nbors)
  class(list_t):: self
  class(link_t), pointer:: link
  !.............................................................................
  class(link_t), pointer:: scan
  class(nbor_t), pointer:: nbor, nbors
  logical:: overlaps, rejected
  integer:: check, count
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%init_nbors_direct', itimer=itimer)
  check = 0
  nullify (nbors)
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
             .and. abs(scan%task%level - link%task%level) <= 1
        if (overlaps) then
          call link%add_nbor (scan, nbors, nbor, rejected=rejected)
          if (rejected) then
            call self%check_duplicates (scan)
          end if
        end if
      end if
    end if
    scan => scan%next                               ! continue search
  end do
  call trace%end (itimer)
END SUBROUTINE init_nbors_direct

!===============================================================================
!> Fast method to generate +-1 nbor lists in cases with very long task lists
!===============================================================================
SUBROUTINE init_nbors_hash (self, link, nbors)
  class(list_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer:: i1, i2, i3, n
  real(8):: size(3), position(3)
  real(8):: times(4)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%init_nbors_hash', itimer=itimer)
  task => link%task
  nullify (nbors)
  !$omp atomic write
  task%n_nbors = 0
  !-----------------------------------------------------------------------------
  ! For a level above levelmin, check all L-1 positions
  !-----------------------------------------------------------------------------
  if (task%level > shared%levelmin) then
    size = task%size*shared%refine_ratio
    do i3=-1,1
    do i2=-1,1
    do i1=-1,1
      !-------------------------------------------------------------------------
      ! Position of level L nbors, and key for a L-1 patches there
      !-------------------------------------------------------------------------
      position = task%position + [i1,i2,i3]*size
      call check_nbor (task%level-1, position, size, link)
    end do
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! At the same level, check all level L positions except self
  !-----------------------------------------------------------------------------
  times(1) = wallclock()
  size = task%size
  do i3=-1,1
  do i2=-1,1
  do i1=-1,1
    if (i1==0 .and. i2==0 .and. i3==0) cycle
    position = task%position + [i1,i2,i3]*size
    call check_nbor (task%level, position, size, link)
  end do
  end do
  end do
  times(1) = wallclock() - times(1)
  !-----------------------------------------------------------------------------
  ! For a level below levelmax, check all L+1 positions
  !-----------------------------------------------------------------------------
  if (task%level < shared%levelmax) then
    times(2) = wallclock()
    size = task%size/shared%refine_ratio
    n = nint(shared%refine_ratio) + 1
    do i3=0,n
    do i2=0,n
    do i1=0,n
      position = task%position + ([i1,i2,i3]-n/2.0)*size
      call check_nbor (task%level+1, position, size, link)
    end do
    end do
    end do
    times(2) = wallclock() - times(2)
    times(3) = wallclock()
  end if
  if (self%verbose > 2 .or. task%id==io%id_debug) &
    write (stdout,*) 'init_nbors_hash:', link%task%id, ' added', &
      task%n_nbors, 'nbors'
  call trace%end (itimer)
contains
!===============================================================================
!> Check if nbor task exists, and if it overlaps -- if so add to nbor list
!===============================================================================
subroutine check_nbor (level, position, size, link)
  integer:: level
  real(8):: position(3), size(3)
  class(link_t), pointer:: link
  integer:: key(4)
  class(link_t), pointer:: nblink
  class(task_t), pointer:: nbtask
  logical:: added
  !-----------------------------------------------------------------------------
  ! Level and positions key, with periodic wrapping
  !-----------------------------------------------------------------------------
  key = task%hash_key (level=level, position=position, size=size)
  call self%hash_table%get (key, nblink)
  if (associated(nblink)) then
    select type (link)
    class is (halo_t)
      call link%conditionally_add_nbor (nblink, nbors, added)
      if (added) then
        !$omp atomic update
        task%n_nbors = task%n_nbors + 1
        nbtask => nblink%task
        if (self%verbose > 2 .or. task%id==io%id_debug) &
          write (stderr,'(a,2(i6,":",i2.2))') &
          'init_nbors_hash:', task%id, task%level, nbtask%id, nbtask%level
      end if
    end select
  end if
end subroutine check_nbor
END SUBROUTINE init_nbors_hash

!===============================================================================
!> Update the nbor list of task, by searching among the nbors of nbors.  The
!> search may be limited by excluding levels larger that self, since their nbors
!> must also be nbors of their parent task.
!===============================================================================
SUBROUTINE init_nbors_expand (self, link, nbors)
  class(list_t)          :: self
  class(link_t)          :: link
  class(nbor_t), pointer :: nbors
  !.............................................................................
  class(task_t), pointer :: task, nbtask
  class(nbor_t), pointer :: nbor, nbor2
  type(hash_table_t)     :: hash_table
  logical                :: added
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  ! Save a pointer to the current nbor list, to remove after adding a new
  ! (to ensure that task%n_needed does not touch zero prematurely)
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%init_nbors_expand', itimer=itimer)
  call hash_table%init (5000, 1)
  !-----------------------------------------------------------------------------
  ! Run through (only) tasks at the relevant levels
  !-----------------------------------------------------------------------------
  task => link%task
  !.............................................................................
  if (self%verbose > 0) then
    select type(task)
    class is (patch_t)
    write (stdlog,'(f12.6,a,i6,":",i2.2,2x,4i4)') wallclock(), &
      pad_trace(1:indent)//'init_nbors: id,key=', &
      task%id, task%level, task%hash_key()
    end select
  end if
  !---------------------------------------------------------------------------
  ! Add the self task to the hash table, to avoid trying to actually add it
  !---------------------------------------------------------------------------
  call hash_table%set (link%task%id, link%task%id)
  select type (link)
  class is (halo_t)
    nbor => link%nbor
    do while (associated(nbor))
      if (abs(nbor%task%level - task%level) == 1) then
        !-----------------------------------------------------------------------
        ! Run through the nbor list of +-1 nbors
        !-----------------------------------------------------------------------
        nbor2 => nbor%link%nbor
        do while (associated(nbor2))
          nbtask => nbor2%task
          if (abs(nbtask%level - task%level) == 2) then
            if (.not. hash_table%exists (nbtask%id)) then
              !-----------------------------------------------------------------
              ! Mark +-2 nbors as tested and add if they overlap
              !-----------------------------------------------------------------
              call hash_table%set (nbtask%id, nbtask%id)
              call link%conditionally_add_nbor (nbor2%link, nbors, added)
              if (added) then
                write (stderr,'(a,2(i6,":",i2.2))') 'init_nbors_expand:', &
                  task%id, task%level, nbtask%id, nbtask%level
              end if
            end if
          end if
          nbor2 => nbor2%next
        end do
      end if
      nbor => nbor%next
    end do
  end select
  call hash_table%dealloc
  call trace%end (itimer)
END SUBROUTINE init_nbors_expand

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
  call link%task%log ('set_init_nbors')
  nbor => link%nbor
  do while (associated(nbor))
    call nbor%task%log ('set init_nbors')
    call nbor%task%set (bits%init_nbors)
    if (self%verbose > 1) &
      write (io_unit%log,*) 'set_init_nbors:', link%task%id, nbor%task%id, &
        nbor%task%istep, nbor%task%is_set(bits%init_nbors)
    !call nbor%link%lock%unset_read ('set_init_nbors')
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
    !call nbor%link%lock%set_read ('set_update_status')
    call nbor%task%set (bits%update_status)
    if (self%verbose > 1) &
      write (io_unit%log,'(1x,a,2i6)') &
        'list_t%set_update_status: task, nbtask =', link%task%id, &
        nbor%task%id, nbor%task%is_set(bits%update_status)
    !call nbor%link%lock%unset_read ('set_update_status')
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
!print *, link%task%id, 'set init_nbor_nbors'
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
    print *, link%task%id, 'list_t%refresh_nbors: init_nbors'
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
  integer:: i, i0, i1, i2, m, n=0, nv=0, nb=0
  real(8):: time
  integer, save:: itimer=0
  !.............................................................................
  if (self%hands_off_nbor_lists) return
  call trace%begin('list_t%init_all_nbors', itimer=itimer)
  if (self%verbose > 0) &
    call timer%tic (time)
  call init_nonleaf (self)                              ! set non_leaf bit
  i2 = self%n/omp%nthreads + 1                          ! tasks per thread
  n = 0
  !$omp parallel private(i,m,i0,i1,link) shared(self,i2,n,stdout) default(none)
  i0 = omp%thread*i2
  i1 = i0 + i2
  if (omp%master) &
    write (stdout,*) self%n, 'task to generate nbor lists for'
  link => self%head                                     ! start checking links
  i = 0                                                 ! task counter
  do while (associated(link))
    if (i >= i0 .and. i < i1) then                      ! thread interval
      call init_nbors (self, link)
      !$omp atomic capture
      n = n+1
      m = n
      !$omp end atomic
      if (mod(m,5000)==0) write (stdout,*) m, 'nbor lists'
    end if
    i = i+1                                             ! increment counter
    link => link%next                                   ! continue list
  end do
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! With new nbor lists it makes sense to reset status bits and count status
  !-----------------------------------------------------------------------------
  call self%reset_status                                ! update task bits
  if (self%verbose > 0) &
    call timer%toc ('time to initialize nbor lists', 1, time)
  call trace%end (itimer)
END SUBROUTINE init_all_nbors

!===============================================================================
!> Find parent patches, which have higher level patches inside
!===============================================================================
SUBROUTINE find_parent (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, parent, item
  integer:: key(4), key1(4)
  class(*), pointer:: ptr
  !-----------------------------------------------------------------------------
  ! Check that a parent is actually needed
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%find_parent')
  if (link%task%level == max(shared%levelmin,shared%levelroot)) then
    if (self%verbose > 2) &
      write (stderr,*) link%task%id, 'needs no parent, levels=', &
      link%task%level, shared%levelmin, shared%levelroot
    call trace%end ()
    return
  end if
  !.............................................................................
  if (self%verbose > 0) &
    write (stdlog,*) link%task%id, 'searching for parent'
  !-----------------------------------------------------------------------------
  ! Call to hash_key fully overloadable here
  !-----------------------------------------------------------------------------
  key = link%task%hash_key (level=link%task%level-1, &
                            size=link%task%size*shared%refine_ratio)  
  call self%hash_table%get (key, parent)
  if (associated(parent)) then
    link%parent => parent
    link%task%parent => parent%task
    if (self%verbose > 0) then
      write (stdlog,'(a,2(i6,":",i2.2,3i3))') 'list_t%find_parent: id,parent=', &
        link%task%id, link%task%key, parent%task%id, parent%task%key
    end if
  else
    write (stderr,'(a,i6,2(":",i2.2,3i3,4x))') &
      'list_t%find_parent: id, no parent=', &
      link%task%id, link%task%key, key
    call check_hash_table (self)
    call io%abort ('list_t%find_parent: missing')
  end if
  !.............................................................................
  call trace%end ()
END SUBROUTINE find_parent

!===============================================================================
!> Find parent patches, which have higher level patches inside
!===============================================================================
SUBROUTINE find_parents (self)
  class(list_t):: self
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin('list_t%find_parents')
  link => self%head
  do while (associated(link))
    if (link%task%level > shared%levelmin) then
      call self%find_parent (link)
    end if
    link => link%next
  end do
  call trace%end ()
END SUBROUTINE find_parents

!===============================================================================
!> Check task list hash table
!===============================================================================
SUBROUTINE check_hash_table (self, label)
  class(list_t):: self
  character(*), optional:: label
  !.............................................................................
  class(link_t), pointer:: link
  integer:: key(4)
  logical:: assoc
  class(*), pointer:: ptr
  !-----------------------------------------------------------------------------
  if (self%verbose > 1) then
    if (present(label)) then
      write (stdlog,*) 'check_hash_table: '//label
    else
      write (stdlog,*) 'check_hash_table: '
    end if
    call trace%begin ('list_t%check_hash_table')
    link => self%head
    do while (associated(link))
      key = link%task%key
      call self%hash_table%get_ptr (key, ptr)
      assoc = associated(ptr)
      select type (ptr)
      class is (link_t)
        write (stdlog,'(a,i6,4i3,l3,i6,2x,a)') &
          'check_hash_table:', link%task%id, key, &
          link%task%is_set (bits%virtual), ptr%task%id, 'link'
      class is (patch_t)
        write (stdlog,'(a,i6,4i3,l3 ,i6,2x,a)') &
          'check_hash_table:', link%task%id, key, &
          link%task%is_set (bits%virtual), ptr%id, 'patch'
      class default
        write (stdlog,'(a,i6,4i3,l3,l1,2x,a)') &
          'check_hash_table:', link%task%id, key, &
          link%task%is_set (bits%virtual), assoc, 'other'
      end select
      link => link%next
    end do
    call trace%end ()
  end if
END SUBROUTINE check_hash_table

!===============================================================================
!> Remove parent patches, which have higher level patches inside
!===============================================================================
SUBROUTINE remove_parents (self)
  class(list_t):: self
  class(link_t), pointer:: link, scan, parent
  !.............................................................................
  call trace%begin('list_t%remove_parents')
  call                                          self%lock%set ('remove_parents')
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
        if (all(abs(scan%task%position-link%task%position) &
            < scan%task%size*0.55d0)) then
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
        call self%remove (parent)
      end if
    end do
    link => link%next                                   ! continue list
  end do
  call                                        self%lock%unset ('remove_parents')
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
  !$omp atomic update
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
    write (io_unit%log,*) &
      trim(self%name)//': inconsistent append of '//task_list%name
  end if
                                            call self%lock%unset ('append_list')
  call trace%end()
END SUBROUTINE append_list

!===============================================================================
!===============================================================================
SUBROUTINE prepend_list (self, task_list)
  class(list_t):: self
  class(list_t), pointer:: task_list
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: n, min_level, max_level
  integer, save:: nprint=100
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%append_list')
              call io%header ('-list_t%prepend_list', left=25)
              call timer%tic ()
  call self%lock%set ('append_list')
  write (stdout,*) 'append_list', self%n, task_list%n, associated(self%tail)
  task_list%tail%next => self%head
  self%head => task_list%head
  self%n = self%n + task_list%n
  write (stdout,*) 'append_list', self%n, associated(self%tail)
  link => self%head
  n = 0
  min_level = 2**30
  max_level = 0
  do while (associated(link))
    task => link%task
    n = n+1
    min_level = min (task%level, min_level)
    max_level = max (task%level, max_level)
    if (task%level < 0) then
      if (nprint > 0) then
        nprint = nprint - 1
        write (stdout,*) n, task%id, task%level, task%key
        select type (task)
        class is (patch_t)
          write(stdout,*) associated(task%erot)
        end select
      end if
      task%level = task%key(1)
    end if
    link => link%next
  end do
  call self%lock%unset ('append_list')
  if (n==self%n) then
    write (stdout,*) trim(self%name)//': OK append of '//task_list%name
    write (stdout,*) 'min/max_level=', min_level, max_level
  else
    write (stderr,*) trim(self%name)//': inconsistent append of '//task_list%name
  end if
              call timer%toc ()
  call trace%end()
END SUBROUTINE prepend_list

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
SUBROUTINE queue_by_time (self, this, istep)
  class(list_t):: self
  class(link_t), pointer:: this
  logical:: error
  integer:: verbose, istep
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%queue_by_time ', itimer=itimer)
  if (self%verbose > 3) then
    write (io%output,*) wallclock(),' thread',omp%thread,' waitfor tasklist(4)'
  end if
  call real_queue_by_time (self, this, istep, error)
  if (error) then
    io%do_trace = .true.
    verbose = io%verbose
    io%verbose = 4
                                             call self%qlock%set ('queue error')
    call self%qconsistency (999)
    call self%check_all (repair=.true.)
    call real_queue_by_time (self, this, istep, error)
                                           call self%qlock%unset ('queue error')
    if (error) then
      print *,'ERROR: queue repair failed'
      call io%abort('queue_by_time: repair failed')
    end if
    io%verbose = verbose
    io%do_trace = .false.
  end if
  call trace%end (itimer)
END SUBROUTINE queue_by_time

!===============================================================================
!> Actually queue the task, unless the task status bits have changed in the
!> meantime.  This can happen in the time window between when the task was
!> first accepted for a call to check_ready() .. even with additional checks
!> in check_ready() cleared -- during the time it takes to loop through the
!> ready-queue to find the proper place to add the task.
!===============================================================================
SUBROUTINE real_queue_by_time (self, this, istep, error)
  class(list_t):: self
  class(link_t), pointer:: this
  integer:: istep
  logical:: error, ok
  class(link_t), pointer:: next, prev, link
  class(task_t), pointer:: task
  character(16):: label
  integer:: n, i, j, nq, ignore
  !-----------------------------------------------------------------------------
  ! Refuse to queue virtual tasks
  !-----------------------------------------------------------------------------
  error=.false.
  ok = .true.
  task => this%task
  ignore = bits%virtual + bits%external + bits%frozen + bits%remove
  if (task%is_set (ignore) .or. task%rank /= mpi%rank) then
    write(stderr,'(3i6,4x,a,2x,4l1)') task%rank, mpi%rank, omp%thread, &
      'ERROR: tried to queue virtual, frozen or removed task, VFRE=', &
      task%is_set(bits%virtual), task%is_set(bits%frozen), &
      task%is_set(bits%remove), task%is_set(bits%external)
    write(io_unit%log,'(3i6,4x,a,2x,4l1)') task%rank, mpi%rank, omp%thread, &
      'ERROR: tried to queue virtual, frozen or removed task, VFRE=', &
      task%is_set(bits%virtual), task%is_set(bits%frozen), &
      task%is_set(bits%remove), task%is_set(bits%external)
    return
  end if
  !-----------------------------------------------------------------------------
  ! Lock the queue while inserting the task, and mark the task as being in the
  ! ready queue (or executing -- the bit is cleared only after updating).
  ! Double check that no other thread has set the ready bit on this task.
  !-----------------------------------------------------------------------------
  if (task%is_clear (bits%ready)) then
    if (self%verbose > 1 .and. task%id == omp%thread_task) then
      write (stdlog,'(a,i6,":",i2.2)') 'queue_by_time: self-schedule task', &
        task%id, task%level
    end if
    !---------------------------------------------------------------------------
    ! Briefly lock the queuw, in case two threads arrive here at the same time.
    ! Setting bits%ready on a task must be done before the task is put back
    ! onto the ready-queue, to ensure a consistent state there. Checking must
    ! be on both %ready and %busy, since the task could both have been queued
    ! and started running in the meantime.
    !---------------------------------------------------------------------------
    if (task%is_set (bits%ready+bits%busy)) then
      call task%log ('reject1', 2)
      return
    end if
    !---------------------------------------------------------------------------
    ! The actual manipulation of the ready-queue must be done under %qlock, and
    ! must not overlap with the brief task locking above
    !---------------------------------------------------------------------------
    call                                        self%qlock%set ('queue_by_time')
    nullify (prev)
    next => self%queue
    n = 0
    do while (associated(next))
      if (next%task%is_clear (bits%remove)) &
        n = n+1
      if (associated(next%task, task)) then
        !.......................................................................
        if (self%verbose > 0) &
          write (io_unit%log,*) omp%thread, ' WARNING: task', task%id, &
            ' is already in ready queue'
        !.......................................................................
        call task%log ('reject2', 2)
        ok = .false.
        exit
      else if (next%task%time > task%time) then
        exit
      end if
      prev => next
      next => next%next_time
    end do
    !...........................................................................
    if (self%verbose > 0) then
      !$omp atomic read
      nq = self%nq
      if (n > nq) then
        write (stdlog,*) 'WARNING:  queue_by_time, nq,n=', self%nq, n
      end if
    end if
    !---------------------------------------------------------------------------
    ! Check one more time if the task busy bit has been set in the meantime.
    ! or if possibly the task time has been updated (can happen if this thread
    ! has been delayed or suspended)
    !---------------------------------------------------------------------------
    if (task%is_set(bits%busy+bits%ready)) then
      call task%log ('reject3', 2)
      ok = .false.
    end if
    if (task%istep /= istep) then
      call task%log ('reject3', 2)
      ok = .false.
    end if
    !---------------------------------------------------------------------------
    ! If OK, add the task to the queue, AFTER setting the %ready bit
    !---------------------------------------------------------------------------
    if (ok) then
      !-------------------------------------------------------------------------
      ! Keep track of queue size
      !-------------------------------------------------------------------------
      !$omp atomic capture
      self%nq = self%nq+1
      nq = self%nq
      !$omp end atomic
      !-------------------------------------------------------------------------
      ! Set the ready bit before actually queuing the task.  The ready bit is
      ! used only as advice that this task is either in the queue or else
      ! is being updated, and that it should therefore not be tested or queued
      !-------------------------------------------------------------------------
      call task%set (bits%ready)
      call task%log ('set ready', 2)
      if (associated(next)) then
        this%next_time => next
        this%prev_time => prev
        if (associated(prev)) then
          prev%next_time => this
          call task%log ('queue1', 2)
          !.....................................................................
          if (self%verbose > 1) &
            write (stdlog,1) 'queue_by_time: id=', task%id, 'nq=', nq, &
            'between', prev%task%id, 'and', next%task%id
          1 format(1x,4(a,i6,5x))
        else
          self%queue => this
          call task%log ('queue2', 2)
          !.....................................................................
          if (self%verbose > 1) &
            write (stdlog,1) 'queue_by_time: id=', task%id, 'nq=', nq, &
            'at head'
        end if
      !-------------------------------------------------------------------------
      ! No place found, so add after prev if exists, or else at head
      !-------------------------------------------------------------------------
      else
        this%prev_time => prev
        nullify (this%next_time)
        if (associated(prev)) then
          prev%next_time => this
          call task%log ('queue3', 2)
          !.....................................................................
          if (self%verbose > 0) then
            write (stdlog,1) 'queue_by_time: id=', task%id, 'nq=', nq, &
            'after', prev%task%id
          end if
        else
          self%queue => this
          call task%log ('queue4', 2)
          !.....................................................................
          if (self%verbose > 0) then
            write (stdlog,1) 'queue_by_time: id=', task%id, 'nq=', nq, &
            'is queue'
          end if
        end if
      end if
      !-------------------------------------------------------------------------
      !task%ostep = task%istep
      !$omp atomic
      mpi_mesg%n_ready = mpi_mesg%n_ready+1
      call self%remove_active (this)                 ! remove from active_queue
    end if
    call                                      self%qlock%unset ('queue_by_time')
  !-----------------------------------------------------------------------------
  ! Since we clear the ready bit in task_list%update, outside of any critical
  ! region, it is perfectly normal that several threads compete about putting
  ! a given task onto the queue, so there is no reason to print a warning here
  !-----------------------------------------------------------------------------
  else
    ok = .false.
    call task%log ('reject5', 2)
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
  if (io%verbose > 1) &
    write (io_unit%log,'(1x,i4,1x,a,i9,a,f10.5,i5,l4)') omp%thread, &
      'list_t%queue_active: adding task', task%id, ' at time', task%atime, &
      self%nac, associated(self%active)
  nullify (prev)
  !$omp critical (active_cr)
  self%nac = self%nac+1
  next => self%active
  do while (associated(next))
    if (associated(next%task, task)) then
      write (io_unit%log,*) omp%thread, ' WARNING: task', task%id, &
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
!> Add a link node into the task list, in increasing level order.  This
!> procedure relies on the assumption that levels are introduced in increasing
!> order, so if the first level pointer does not exit, the insert-after point
!> is tail.  Since this is much faster than a linear search, one should make
!> that the assumption is correct, rather then use a slower procedure
!===============================================================================
SUBROUTINE add_by_level (self, this, count, keep)
  class(list_t)          :: self
  class(link_t), pointer :: this
  logical, optional      :: count, keep
  !.............................................................................
  class(link_t), pointer :: next, scan
  class(task_t), pointer :: task
  logical                :: found
  integer                :: level, l
  integer, save          :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%add_by_level', itimer=itimer)
  call                                            self%lock%set ('add_by_level')
  level = this%task%level
  self%levelmin = min(level,shared%levelmin)
  self%levelmax = max(level,shared%levelmax)
  if (level < 1 .or. level > max_level) then
    write (stderr,*) 'list_t%add_by_level: level=', level
    call io%abort ('list_t%add_by_level: level outside range')
  end if
  !-----------------------------------------------------------------------------
  ! Check presence of task%link, unless keep is true
  !-----------------------------------------------------------------------------
  if (.not. shared%option(keep)) then
    task => this%task
    select type (task)
    class is (patch_t)
        task%link => this
    end select
  end if
  !-----------------------------------------------------------------------------
  ! If head is null this is the first list element
  !-----------------------------------------------------------------------------
  found = .false.
  if (.not.associated(self%head)) then
    nullify (this%prev, this%next)
    self%head => this
    self%tail => this
    self%first(level)%pointer => this
    self%levelmin = level
    self%levelmax = level
    found = .true.
  !-----------------------------------------------------------------------------
  ! If the level pointer exists, prepend to it and move it
  !-----------------------------------------------------------------------------
  else if (associated(self%first(level)%pointer)) then
    next => self%first(level)%pointer
    call insert_before (next, this, level, self)
    if (associated(next,self%head)) then
      self%head => this
    end if
    self%first(level)%pointer => this
    found = .true.
  !-----------------------------------------------------------------------------
  ! If not, searc for an insert point and set it
  !-----------------------------------------------------------------------------
  else
    !---------------------------------------------------------------------------
    ! Search for a valid pointer to start from
    !---------------------------------------------------------------------------
    scan => self%head
    do l=level-1,self%levelmin,-1
      if (associated(self%first(l)%pointer)) then
        scan => self%first(l)%pointer
        exit
      end if
    end do
    !---------------------------------------------------------------------------
    ! Search for insert point
    !---------------------------------------------------------------------------
    do while (associated(scan))
      if (scan%task%level > level) then
        call insert_before (scan, this, level, self)
        found = .true.
        exit
      end if
      scan => scan%next
    end do
    !---------------------------------------------------------------------------
    ! If no task found with higher level, insert as tail
    !---------------------------------------------------------------------------
    if (.not.found) then
      this%prev => self%tail
      nullify (this%next)
      self%tail%next => this
      self%tail => this
    end if
    self%first(level)%pointer => this
  end if
  self%levelmin = min(level,self%levelmin)
  self%levelmax = max(level,self%levelmax)
  self%n = self%n+1
  self%max_id = max(self%max_id,this%task%id)
                                           call self%lock%unset ('add_by_level')
  !-----------------------------------------------------------------------------
  ! Count and unlock
  !-----------------------------------------------------------------------------
  if (present(count)) then
    call self%count_status ()
  end if
  call self%hash_table%add (this%task%key, this)
  call this%task%log ('add_by_level', 4)
  call trace%end (itimer)
contains
!===============================================================================
subroutine insert_before (next, this, level, self)
  class(list_t):: self
  class(link_t), pointer:: this, next
  integer:: level
  !-----------------------------------------------------------------------------
  this%prev => next%prev          ! might be null
  this%next => next               ! exists
  if (associated(this%prev)) then ! might be null
    this%prev%next => this        ! wrap backwards
  end if
  this%next%prev => this          ! wrap forwards
end subroutine
!===============================================================================
END SUBROUTINE add_by_level

!===============================================================================
!> Add a link node into the task list, in increasing quality order
!===============================================================================
SUBROUTINE add_by_quality (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  !.............................................................................
  call trace%begin ('list_t%add_by_quality')
  this%task%level = floor(this%task%quality)
  call self%add_by_level (this)
  call trace%end ()
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
  call                                        self%lock%set_read ('check_queue')
  link => self%queue
  previous = -1.0
  nq = 0
  do while (associated(link))
    if (link%task%time < previous) &
      write (io_unit%log,*) &
        'ERROR: queue out of order at', link%task%id, link%task%time, &
      previous
    previous = link%task%time
    link => link%next_time
    nq = nq+1
  end do
  if (nq /= self%nq) &
    write (io_unit%log,*) &
      'ERROR: inconsistent number of tasks in queue', nq, self%nq
  call                                      self%lock%unset_read ('check_queue')
END SUBROUTINE check_queue

!===============================================================================
!> Insert 'new' task before 'old' task
!===============================================================================
SUBROUTINE insert (self, old, new)
  class(list_t):: self
  class(link_t), pointer:: old, new
  !.............................................................................
  call trace%begin ('list_t%insert')
  call                                                  self%lock%set ('insert')
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
  call                                                self%lock%unset ('insert')
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
  call                                  self%lock%set_read ('print_queue_until')
  link => self%queue
  if (io%master) write (io_unit%log,*)'print_queue_until:'
  do while (associated(link))
    write (io_unit%log,'("queue: id, t0, t, t0-t =",i7,1p,3e12.3)') &
      link%task%id, task%time, link%task%time, task%time-link%task%time
    if (link%task%id==task%id) exit
    link => link%next_time
  end do
  call                                self%lock%unset_read ('print_queue_until')
END SUBROUTINE print_queue_until

!===============================================================================
!> Print the queue
!===============================================================================
SUBROUTINE print_queue_times (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call                                  self%lock%set_read ('print_queue_times')
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
  call                                self%lock%unset_read ('print_queue_times')
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
  if (io%verbose > 0) &
    write (io_unit%log,*) 'list_t%print: ',trim(self%name), self%n
  call                                         self%lock%set_read ('print_list')
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
    if (io%verbose > 1) &
      write (io_unit%log,'(a,i9,f10.6,3x,3f10.6,i4,i6,l5,2x,a)') &
        'list_t%print:', task%id, task%time, task%position, task%status, &
        task%n_nbors, associated(link%nbor), type
    link => link%next
  end do
  call                                       self%lock%unset_read ('print_list')
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
    write (stdlog,'(f12.6,1x,a,i6,l3,1x,3l1)') wallclock(), &
      'list_t%swap: task,sender,BVS=', link%task%id, sender, &
      link%task%is_set(bits%boundary), link%task%is_set(bits%virtual), &
      link%task%is_set(bits%swap_request)
  end if
  one = merge(+1,-1,sender)
  !$omp atomic update
  self%na = self%na-one                  ! active tasks
  !$omp atomic update
  self%nb = self%nb-one                  ! boundary tasks
  !$omp atomic update
  self%nv = self%nv+one                  ! virtual tasks
  !$omp atomic update
  timer%levelpop (link%task%level) = &    ! level population count
  timer%levelpop (link%task%level) - one
  !-----------------------------------------------------------------------------
  ! Set rank to the new owner on the local rank, and flip the swap_request bit
  !-----------------------------------------------------------------------------
  if (sender) then
    call link%task%log ('swap (send)')
    link%task%rank = rank
    call link%task%set (bits%swap_request)
    call change_running_mean (-0.1)
    if (self%verbose > 0) &
      write (io_unit%log,*) 'list_t%swap: to rank =', link%task%id, rank
  !-----------------------------------------------------------------------------
  ! On recv side, clear VS, set B, and set_init_nbors
  !-----------------------------------------------------------------------------
  else
    call link%task%log ('swap (recv)')
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
    call link%task%clear (bits%swap_request+bits%ready)
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
  !call                                           link%lock%set ('update_status')
  nbor => link%nbor
  do while (associated(nbor))
    old%status = link%task%status
    call link%check_nbor_relations (nbor%link%task, nbor%needed,  &
                                    nbor%needs_me, nbor%download)
    call update_count (link%task, old, bits%virtual , self%nb)
    call update_count (link%task, old, bits%boundary, self%nv)
    call update_count (link%task, old, bits%internal, self%na)
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! Clear the update_status bit, unless the task is a boundary task, in which
  ! case the bit needs to be present also on the nbor ranks, where it will be
  ! cleared at this point.
  !-----------------------------------------------------------------------------
  if (link%task%is_clear (bits%boundary)) &
    call link%task%clear (bits%update_status)
  !call                                         link%lock%unset ('update_status')
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
  character(28):: label
  integer:: ierr, rank, i
  logical, allocatable:: active(:)
  integer, save, allocatable:: sequence(:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Error traps
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
    if (associated(task%parent)) then
      write (io_unit%log,'(f12.6,1x,a,2(i6,":",i2.2,3i3),f12.6)') &
        wallclock(), &
        'list_t%send_to_vnbors: task,time,parent=', link%task%id, &
        link%task%key, link%parent%task%id, link%parent%task%key, link%task%time
    else
      write (io_unit%log,'(f12.6,1x,a,i6,":",i2.2,3i3,f12.6)') wallclock(), &
        'list_t%send_to_vnbors: task,time=', link%task%id, &
        link%task%key, link%task%time
    end if
  end if
  !-----------------------------------------------------------------------------
  ! Make sure the task%mesg tag is unique (needs to be here to allow calls from
  ! within list_mod.f90)
  !-----------------------------------------------------------------------------
  mesg => task%mesg
  mesg%id = task%id
  !-----------------------------------------------------------------------------
  ! Collect active ranks in a table containing only the ranks to send to
  !-----------------------------------------------------------------------------
  allocate (active(0:mpi%size-1))
  call os%mem_register (storage_size(active), mpi%size, 'active')
  rank = -1                                                     ! previous rank
  active = .false.
  !-----------------------------------------------------------------------------
  ! No need to lock the nbor list since thris threads owns it
  !-----------------------------------------------------------------------------
  !call link%lock%set ('send_to_vnbors')
  nbor => link%nbor                                             ! first nbor
  do while (associated(nbor))
    if (nbor%needs_me) then
      rank = nbor%task%rank
      if (rank /= mpi%rank) then
        active(rank) = .true.
      end if
      if (io%task_logging > 3) then
        write (label,'(a,i6,i4,l3)') 'nbor   ', nbor%task%id, rank, active(rank)
        call task%log (label, 4)
      end if
    endif
    nbor => nbor%next
  end do
  !call link%lock%unset ('send_to_vnbors')
  !-----------------------------------------------------------------------------
  ! Send the task to the active ranks
  !-----------------------------------------------------------------------------
  do rank=0,mpi%size-1
    if (active(rank)) then
      !-------------------------------------------------------------------------
      ! Use the given task%tag as the message tag
      !-------------------------------------------------------------------------
      mesg%id = task%id
      mesg%seq = task%seq
      !.........................................................................
      if (io%task_logging > 1) then
        write (label,'(a,3i7)') 'send', rank, mesg%tag, task%tag
        call task%log (label, 2)
      end if
      call mesg%send (rank, tag=mesg%tag)
      if (task%is_set (bits%send_new)) then
        call task%clear (bits%send_new)
        call task%log ('clear send_new', 3)
      end if
      !.........................................................................
      ! Ranks file tracking
      !.........................................................................
      if (self%verbose > 0 .or. task%id == io%id_debug) then
        write (stdlog,'(f12.6,1x,a,i9,1p,e18.6,i9,2x,a,i5,2x,5l1)') &
          wallclock(), 'send_to_vnbors: id,time,tag =', &
          mesg%id, task%time, mesg%tag, 'to', rank, &
          task%is_set (bits%internal), &
          task%is_set (bits%boundary), &
          task%is_set (bits%virtual), &
          task%is_set (bits%external), &
          task%is_set (bits%swap_request)
        flush (stdlog)
      end if
      !.........................................................................
      ! Task event file tracking
      !.........................................................................
      if (self%verbose > 1) then
        write (io_unit%queue,'(f12.6,i4,1x,a12,i6,f12.6,2i6)') &
          wallclock(), omp%thread, 'sent', task%id, task%time, rank, mesg%tag
      end if
    end if
  end do
  deallocate (active)
  call task%clear (bits%swap_request)
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
contains
  subroutine log_tag (label)
    character(*):: label
    character(32):: s
    write (s,'(a,2i6)') label, link%task%seq, mesg%tag
    call link%task%log (trim(s), 3)
  end subroutine
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
  integer:: status, test, n, ignore
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
  ignore = bits%internal + bits%boundary + bits%virtual+bits%external
  test = iand(self%task%status,ignore)
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
    if (n==0) &
      write (io_unit%log,*) self%task%id, ' inconsistent bits: has no nbors'
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
  integer:: level, n, n1
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%info', itimer=itimer)
  if (present(label)) then
    write (stdlog,*) 'list_t%info: '//trim(label)
  else
    write (stdlog,*) 'list_t%info: '
  end if
  link => self%head
  do while (associated(link))
    call link%info
    link => link%next
  end do
  !write (stdlog,'(a,2x,2i3)') 'list_t%info first(level)%pointer:', &
  !  shared%levelmin, shared%levelmax
  !n = 0
  !do level=shared%levelmin,shared%levelmax
  !  link => self%first(level)%pointer
  !  if (associated(link)) then
  !    task => link%task
  !    if (self%verbose > 0) &
  !      write (stdlog,1) wallclock(), &
  !        pad_trace(1:indent)//'list_t%info: level,task%id,task%level=', &
  !        level, task%id, task%level
  !      1 format(f12.6,a,i3,i7,":",i2.2)
  !    n1 = 0
  !    do while (link%task%level == level)
  !      n1 = n1+1
  !      link => link%next
  !      if (.not.associated(link)) &
  !        exit
  !    end do
  !    write (stdlog,*) 'level, n =', level, n1
  !    n = n+n1
  !  else
  !    write (stdlog,1) wallclock(), &
  !      pad_trace(1:indent)//'list_t%info: empty level=', level
  !  end if
  !end do
  !if (n /= self%n) then
  !  write (stderr,*) 'list_t%info: number of tasks in task list:', self%n
  !  write (stderr,*) 'list_t%info: number of tasks in level lists:', n
  !  call io%abort ('list_t%info: inconsistent task list level pointers')
  !end if
  call trace%end(itimer)
END SUBROUTINE info

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
SUBROUTINE check_nbor_consistency (self, label)
  class(list_t):: self
  character(len=*), optional:: label
  !.............................................................................
  logical:: ok, found
  class(link_t), pointer:: link1, link2
  class(task_t), pointer:: task1, task2
  class(nbor_t), pointer:: nbor1, nbor2
  logical:: test
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%check_nbor_consistency')
  if (present(label) .and. self%verbose > 0) then
    write (stdout,*) 'list_t%check_nbor_consistency: ', label
  end if
  call timer%tic ()
  if (io%id_debug > 0 .and. first_time) then
    !---------------------------------------------------------------------------
    ! Find the debug link and check all tasks for overlap -- to be compared
    !---------------------------------------------------------------------------
    link1 => self%head
    task1 => link1%task
    select type (task1)
    class is (patch_t)
      do while (associated(link1))
        if (task1%id == io%id_debug) then
          link2 => self%head
          do while (associated(link2))
            task2 => link2%task
            select type (task2)
            class is (patch_t)
              test = task1%overlaps2 (task2)
            end select
            link2 => link2%next
          end do
          exit
        end if
        link1 => link1%next
      end do
    end select
    first_time = .false.
  end if
  ok = .true.
  found = .false.
  link1 => self%head
  do while (associated(link1))
    nbor1 => link1%nbor
    do while (associated(nbor1))
      !-------------------------------------------------------------------------
      ! Find this nbor on the task list, and find the link1%task nbor
      !-------------------------------------------------------------------------
      nbor2 => nbor1%link%nbor
      do while (associated(nbor2))
        if (nbor2%task%id == link1%task%id) then
          found = .true.
          if ((nbor2%needed   .neqv. nbor1%needs_me) .or. &
              (nbor2%needs_me .neqv. nbor1%needed)) then
            write (stdlog,'(a,2(i6,2x,2l1,3x))') &
              'list_t%check_nbor_consistency: inconsistency btw tasks', &
              nbor1%task%id, nbor1%needed, nbor1%needs_me, &
              nbor2%task%id, nbor2%needed, nbor2%needs_me
            ok = .false.
          end if
        end if
        nbor2 => nbor2%next
      end do
      if (.not. found) then
        write (stderr,'(a,i6,a,i6)') &
          ' list_t%check_nbor_consistency: adding missing nbor', &
          link1%task%id, ' to task', nbor1%task%id
        call nbor1%link%add_nbor (link1)
      end if
      nbor1 => nbor1%next
    end do
    link1 => link1%next
  end do
  call timer%toc ('list_t: check_nbor_consistency')
  call io%assert (ok, 'list_t%check_nbor_consistency: not consistent')
  call trace%end ()
END SUBROUTINE check_nbor_consistency

!===============================================================================
!> Check for duplicated task
!===============================================================================
SUBROUTINE check_duplicates (self, link)
  class(list_t):: self
  class(link_t), pointer:: link, scan
  integer:: count
  !-----------------------------------------------------------------------------
  count = 0
  scan => self%head
  do while (associated(scan))
    if (scan%task%id == link%task%id) then
      count = count+1
    end if
    scan => scan%next
  end do
  if (count > 1) then
    write (stderr,'(a,i5,a,i5,a)') &
      'task',link%task%id,' occurs',count,' times in task list'
  else if (count==1) then
    write (stderr,'(a,i5,a,i5,a)') &
      'task',link%task%id,' occurs only once in task list'
  else
    write (stderr,'(a,i5,a,i5,a)') &
      'task',link%task%id,' is not in the task list'
  end if
END SUBROUTINE check_duplicates

!===============================================================================
!> Allow extras level procedures to access task_list_t%make_parent_levels()
!===============================================================================
SUBROUTINE make_parent_levels (self)
  class(list_t):: self
  !-----------------------------------------------------------------------------
  call io%abort ('list_t%make_parent_levels should never be called')
END SUBROUTINE make_parent_levels

!===============================================================================
!> Sort the current task list by level, by copying to a new task list and
!> add back by level
!===============================================================================
SUBROUTINE sort_by_level (self)
  class(list_t):: self
  class(list_t), pointer:: list
  class(link_t), pointer:: link, next
  integer:: n
  !.----------------------------------------------------------------------------
              call io%header ('-list_t%sort_by_level', left=25)
              call timer%tic ()
  call                                    self%lock%set ('list_t%sort_by_level')
  allocate (list)
  call os%mem_register (storage_size(list), 1, 'sort_by_level')
  call list%init
  link => self%head
  do while (associated(link))
    next => link%next
    call list%append (link)
    link => next
  end do
  n = 0
  nullify (self%head, self%tail)
  link => list%head
  do while (associated(link))
    next => link%next
    call self%add_by_level (link)
    link => next
    n = n+1
    if (mod(n,5000)==0) &
      write (stdout,*) 'sort_by_level', n
  end do
  call self%count_status ()
  call list%dealloc ()
  deallocate (list)
  call                                  self%lock%unset ('list_t%sort_by_level')
              call timer%toc ()
END SUBROUTINE sort_by_level

!===============================================================================
!> Initialize attributes of a task that are only available at the task list
!> level.  To be overloaded in more complex task list setups.
!===============================================================================
SUBROUTINE init_task (self, task)
  class(list_t):: self
  class(task_t):: task
  !-----------------------------------------------------------------------------
  call trace%begin ('list_t%init_task')
  task%box = self%size
  task%origin = self%position - self%size/2d0
  call trace%end ()
END SUBROUTINE init_task

!===============================================================================
!> Stub, for overloading
!===============================================================================
SUBROUTINE get_head (self, this)
  class(list_t):: self
  class(link_t), pointer:: this
  !-----------------------------------------------------------------------------
END SUBROUTINE get_head

END MODULE list_mod
