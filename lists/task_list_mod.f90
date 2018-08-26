!> =============================================================================
!> Task list data type, with methods for startup and updates.  Message handling
!> is inherited from the task_mesg_t data type.  At this level, task handling
!> is aware of experiment_t data types, and all data types it has inheritet from.
!> =============================================================================
MODULE task_list_mod
  USE io_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  USE mpi_mod
  USE mpi_mesg_mod
  USE list_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE bits_mod
  USE refine_mod
  USE load_balance_mod
  USE task_mesg_mod
  implicit none
  private
  type, public, extends(task_mesg_t):: task_list_t
    logical:: syncing=.false.
    real(8):: sync_time=-1.0, sync_next=-1.0
  contains
    procedure:: init
    procedure:: startup
    procedure:: execute
    procedure:: update
    procedure:: average
    procedure:: append_task_list
    procedure:: print => print1
  end type
  public:: task2patch
CONTAINS

!> =============================================================================
!> Initialize the task list, by first initializing the list tasks, then making
!> neighbor lists, and finally checking if they are ready to execute
!> =============================================================================
SUBROUTINE init (self, name)
  class(task_list_t):: self
  character(len=*), optional:: name
  class(link_t), pointer:: link
  integer:: iostat
  real(8), save:: job_seconds=1d30
  real(8), save:: sync_time=0.0
  namelist /task_list_params/ sync_time, job_seconds
  !-----------------------------------------------------------------------------
  ! An optional namelist can be used to turn debugging on
  !-----------------------------------------------------------------------------
  call trace%begin('task_list_t%init')
  call self%list_t%init (name)
  call mpi_mesg%init
  rewind (io%input); read(io%input, task_list_params, iostat=iostat)
  if (io%master) write (*,task_list_params)
  if (sync_time > 0.0) &
    self%sync_next = sync_time
  self%sync_time = sync_time
  io%job_seconds = job_seconds
  write (io_unit%mpi,*) 'task_list%init: next_sync =', self%sync_next, self%sync_time
  call load_balance%init
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!> =============================================================================
SUBROUTINE startup (self)
  class(task_list_t):: self
  !.............................................................................
  class(link_t), pointer:: link, rem
  class(task_t), pointer:: task
  real(8):: sec
  integer:: dims(4)
  integer, save:: itimer=0
  logical:: was_refined
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%startup', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Initialize nbor lists and setting virtual and boundary bit should NOT be
  ! done here, since at this point the task list may be a splice between the
  ! task lists of several components, with nbor relations that should not be
  ! overwritten.  Each component should return a task list with appropriate
  ! nbor lists, which may then be combined into a more comprehensive setup, with
  ! special handling of nbor relations
  !-----------------------------------------------------------------------------
  !call self%init_bdries
  !call self%init_all_nbors
  !call self%reset_status
  !-----------------------------------------------------------------------------
  ! Check the ready status of all tasks
  !-----------------------------------------------------------------------------
  link => self%head
  do while (associated(link))
    call link%task%clear(bits%ready)
    if (link%task%is_set(bits%virtual)) link%task%wc_last = wallclock()
    call self%check_ready (link)
    if (self%sync_time > 0.0) &
      link%task%sync_time = link%task%time + self%sync_time
    link => link%next
  end do
  if (mpi_mesg%recv_active .and. mpi_mesg%recv_priv) &
    call self%init_virtual
  !-----------------------------------------------------------------------------
  ! Initialize the data I/O, and start timer
  !-----------------------------------------------------------------------------
  if (mpi%master) &
    print '(a,f8.3,a)', ' Memory per process:', io%gb, ' GB'
  call mpi%barrier ('self%execute')
  call tic (time=sec)
  call timer%print
  self%n_tasks = self%na
  call trace%end (itimer)
END SUBROUTINE startup

!> =============================================================================
!> Execute the task list, updating it until it is empty.  With !$omp parallel here,
!> everything local to self%update is thread private.
!> =============================================================================
SUBROUTINE execute (self)
  USE omp_timer_mod
  class(task_list_t):: self
  real(8):: sec
  integer:: dims(4)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%execute', itimer=itimer)
  call self%startup
  call tic (time=sec)
  call timer%print
  !-----------------------------------------------------------------------------
  ! Loop over tasks
  !-----------------------------------------------------------------------------
  !$omp parallel
  do while (self%na > 0 .and. .not.io%do_stop)
    call self%update (self%head, .false.)
  end do
  call timer%print
  call mpi_mesg%diagnostics(1)
  write (io_unit%log,*) 'thread',omp%thread,' arrived'
  flush (io_unit%log)
  !$omp end parallel
  call toc ('wall time', n_update, time=sec)
  call mpi%barrier ('end')
  !if (io%master) &
  !  write (io_unit%log,*) 'download_differ % =', 100.*download%n_differ/ &
  !    max(1,(download%n_differ+download%n_same))
  call trace%end (itimer)
END SUBROUTINE execute

!> =============================================================================
!> Common part of task list update, used by dispatchers
!> =============================================================================
SUBROUTINE update (self, head, test, was_refined)
  class(task_list_t):: self
  class(link_t), pointer:: head
  logical:: test
  logical, optional:: was_refined
  !.............................................................................
  class(task_t), pointer:: task
  real(8):: wc
  logical:: refined
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  ! Check for flag files.  If io%out_time is set, set a new out_next.
  ! If task%time is slightly smaller than a multiple of out_time, then
  ! out_next will becoe that multiple.
  !----------------------------------------------------------------------------
  call trace%begin('task_list_t%update', itimer=itimer)
  task => head%task
  call io%check_flags
  if (task%id==io%id_track) then
    task%track =  .not.task%track
    io%id_track = 0
  end if
  !-----------------------------------------------------------------------------
  ! Check if refinement is needed on the task; if so this will push new tasks
  ! onto the queue, with the same task time; i.e., to the head of the queue.
  ! If the task is virtual, it will not be checked by this rank.
  !-----------------------------------------------------------------------------
  call refine%check_current(self, head, refined)
  if (present(was_refined)) then
    was_refined = refined
  end if
  !---------------------------------------------------------------------------
  ! If the task is frozen (e.g. because it has finished), return
  !---------------------------------------------------------------------------
  if (task%is_set (bits%frozen+bits%virtual)) then
    !if (task%iout == 0) call task%output(self%name)
    ! if a frozen task is set to ready, make it not ready
    if (task%is_set(bits%ready)) then
      !$omp atomic
      self%na = self%na - 1
    end if
    call trace%end (itimer)
    return
  end if
  !----------------------------------------------------------------------------
  ! Download nbor info -- this may be an experiment_t procedure, but if not,
  ! the call should be answered by a solver-specific procedure (which may or
  ! may not choose to call the generic patch guard zone handler in tasks/)
  !----------------------------------------------------------------------------
  if (.not.test) &
    call task%dnload                                  ! download nbor info
  !----------------------------------------------------------------------------
  ! Update the task, whatever that means (may include call to task%output)
  !----------------------------------------------------------------------------
  if (mpi_mesg%delay_ms==0.0) then
    wc = wallclock()
  end if
  task%rotated = .false.
  if (test) then
    call task%test_update
  else
    call task%update                                  ! update the task
  end if
  !----------------------------------------------------------------------------
  ! Default delay time = first update time
  !----------------------------------------------------------------------------
  if (mpi_mesg%delay_ms==0.0) then
    !$omp atomic write
    mpi_mesg%delay_ms = 1d3*(wallclock()-wc)
  end if
  !$omp atomic
  mpi_mesg%n_update = mpi_mesg%n_update+1
  if (io%verbose>1) &
    write (io_unit%log,'(a,i4,2x,a,i7,2x,a,2x,a,i7,2x,a,1p,2g14.6,2x,a,2i5,l3)') &
    'thread', omp_mythread, 'task', task%id, trim(task%type), &
    'step', task%istep, 'dt, time:', task%dtime, task%time, &
    'n, nq', self%n, self%nq, associated(head%nbor)
  !-----------------------------------------------------------------------------
  ! The rotate procedure in the generic patch_t module is responsible for
  ! updating the time, and rotating the memory slots where information about
  ! previous time steps are saved.  Task may signal that they have already done
  ! the rotate internally, or else are not yet ready to do so, by setting the
  ! task%rotated flag
  !-----------------------------------------------------------------------------
  if (.not.task%rotated) then
    call task%rotate
    select type (task)
    class is (patch_t)
      !$omp atomic
      n_update = n_update + product(task%n)
    end select
  end if
  call task%info (self%nq, self%na)                   ! print info on stdout
  !-----------------------------------------------------------------------------
  ! Tasks that have finished are subtracted from the task list count but are
  ! not removed, since their data may be needed by other tasks (including on
  ! other ranks).  They are, however, not added back onto the ready queue.
  !-----------------------------------------------------------------------------
  if (task%has_finished()) then                       ! finished:
    call task%set (bits%frozen)
    call load_balance%active (.false.)                ! turn off load balancing
    !$omp atomic
    self%na = self%na-1
  end if
  !-----------------------------------------------------------------------------
  ! If the task is a boundary patch, first check if it should be given to
  ! another rank.  A patch that has been given to another rank becomes
  ! a virtual patch, and the load balance procedure has already sent those
  ! of its nbors that became new boundary patches over to relevant ranks,
  ! while the patch itself is sent here, with bits%virtual+swap_reqest set.
  !-----------------------------------------------------------------------------
  if (task%is_set(bits%boundary)) then                ! boundary patch?
    head%task%nq = self%nq                            ! make sure to pass on
    if (load_balance%check_load (head)) then          ! sell?
      !$omp atomic
       self%na = self%na-1
      !$omp atomic
       self%nb = self%nb-1
      !$omp atomic
       self%nv = self%nv-1
       call self%count_status                         ! redundant? (FIXME)
    end if
    call head%send_to_vnbors                          ! send to virtual nbors
  end if
  !----------------------------------------------------------------------------
  ! Periodic task syncronization:  The first time a task%time is exactly equal
  ! to load_balance%sync_next, we wait for all ranks to arrive at the same
  ! time. All tasks within a rank are then also at this time, since as long as
  ! one has not arrived there yet, the head%ask%time is smaller.
  !----------------------------------------------------------------------------
  if (self%syncing) then
    write (io_unit%mpi,*) task%id,omp%thread, &
      ' is waiting on a sync at t =', self%sync_next
    !call trace%end; call trace%begin ('mpi%barrier', itimer=itimer)
    call mpi%barrier ('sync')
    !call trace%end (itimer); call trace%begin('dispatcher0_t%update')
    write (io_unit%mpi,*) task%id,omp%thread, &
      ' finished wating on a sync at t =', self%sync_next
    self%sync_next = self%sync_next + self%sync_time
    self%syncing = .false.
  end if
  task%sync_time = self%sync_next
  call trace%end (itimer)
END SUBROUTINE update

!> =============================================================================
!> Average over a variable with index idx
!> =============================================================================
FUNCTION average (self, idx, time)
  class (task_list_t):: self
  integer:: idx
  real(8):: average, time
  !.............................................................................
  class(link_t) , pointer:: link
  class(task_t) , pointer:: task
  real(8), dimension(:,:,:), allocatable:: buffer
  integer:: n, m(3), jt(2), l(3), u(3)
  real:: pt(2)
  !-----------------------------------------------------------------------------
  ! Compute smallest time on the list
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%average')
  link => self%head
  time = link%task%time
  do while (associated(link))
    time = min(time, link%task%time)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Interpolate to that time and sum up averages
  !-----------------------------------------------------------------------------
  link => self%head
  n = 0
  average = 0d0
  do while (associated(link))
    task => link%task
    select type (task)
    class is (patch_t)
      m = task%ncell
      l = task%li
      u = l + m - 1
      call task%time_interval (time, jt, pt)
      if (.not.allocated(buffer)) &
        allocate (buffer(m(1),m(2),m(3)))
      buffer = task%mem(l(1):u(1),l(2):u(2),l(3):u(3),idx,jt(1),1)*pt(1) &
             + task%mem(l(1):u(1),l(2):u(2),l(3):u(3),idx,jt(2),1)*pt(2)
      average = average + sum(buffer)
      n = n+1
      !print '(5i5,2f8.5,3f12.6,1p,e12.4)', idx, n, task%id, jt, pt, task%t(jt(1)), time, task%t(jt(2)), average
    end select
    link => link%next
  end do
  average = average/n/product(shape(buffer))
  call trace%end
END FUNCTION average

!> =============================================================================
!> Append a task list to self
!> =============================================================================
SUBROUTINE append_task_list (self, task_list)
  class (task_list_t):: self
  class (task_list_t), pointer:: task_list
  class (list_t), pointer:: list
  !-----------------------------------------------------------------------------
  call trace%begin ('task_list_t%append_task_list')
  list => task_list
  call self%append_list (list)
  call trace%end
END SUBROUTINE append_task_list

!> =============================================================================
!> Print a table with task id and all neighbors
!> =============================================================================
SUBROUTINE print1 (self, label)
  class (task_list_t):: self
  character(len=*), optional:: label
  class(link_t) , pointer:: link, nbor
  !.............................................................................
  call trace%begin ('task_list_t%print_tasklist')
  if (io%master.and.present(label)) write (io_unit%output,*) label
  link => self%head
  do while (associated(link))
    call link%print_nbors ('')
    link => link%next
  end do
  call self%print_queue
  call trace%end
END SUBROUTINE print1

END MODULE task_list_mod
