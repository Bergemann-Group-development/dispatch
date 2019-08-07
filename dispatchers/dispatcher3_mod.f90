!===============================================================================
!> Execute a task list.  Since only the master thread is calling check_mpi(), we
!> need to make sure that it handles all send/recv, hence send_priv=.false. and
!> recv_priv=.false., and we can use the simplest receive mechanism (recv_active
!> =.true.), with no buffering needed (queue_unpack=.false.).  These values are
!> imposed after the input namelists are read, to make sure the choices are
!> correct.
!>
!> Strategy: Keep track of, separately, the number of ready tasks, the number
!> of currently spawned tasks, and the number of currently busy tasks.  To this
!> end, define the following states:
!>
!> state=0: task is inactive or new
!> state=1: the task is ready to update, but not yet spawned
!> state=2: the task has been spawned
!> state=3: the task is busy updating
!> state=4: the task has finished updating
!>
!> The simplest strategy for counting is to keep track of the number of tasks
!> in each state, by incrementing the count for the new state and decrementing
!> the count for the old state, whenever the state changes. These counts should
!> be part of the dispatcher data type.
!===============================================================================
MODULE dispatcher3_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE mpi_io_mod
  USE timer_mod
  USE link_mod
  USE mpi_mesg_mod
  USE task_list_mod
  USE bits_mod
  USE global_mod
  USE task_mod
  implicit none
  private
  type, public:: dispatcher3_t
    integer:: verbose=0
    integer:: n_spawn=0
    !type(task_list_t), pointer:: task_list => null()
  contains
    procedure:: init
    procedure:: execute
  end type
  integer:: n_fail=0
  integer:: n_state(0:4)=0
  integer:: d_state(0:4)=0
  integer:: verbose=0
  integer:: n_write=2000
  integer(8):: n_check=0, n_ready=0
  real(8):: wc_failed=0d0
  real(8):: start
  logical:: use_critical=.false.
  type(dispatcher3_t), public:: dispatcher3
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, task_list)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  integer:: iostat
  namelist /dispatcher3_params/ verbose, use_critical, n_write
  !.............................................................................
  rewind (io%input)
  read (io%input, dispatcher3_params, iostat=iostat)
  write (io%output, dispatcher3_params)
END SUBROUTINE init

!===============================================================================
!> Simple dispatcher, which spawns other threads that handle tasks that are
!> ready. If / when the master thread has spawned a number of tasks that depend
!> on the implementation, it will participate in execution of the tasks.
!===============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher3_t):: self
  type(task_list_t), pointer:: task_list
  logical:: test, fell_through
  !.............................................................................
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  logical:: mytask
  integer:: n1, n2
  real(8):: start_fail, start_iter
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dispatcher3_t%execute', itimer=itimer)
  call io%header('begin dispatcher3_t%execute: Any opportunity dispatcher')
  call self%init (task_list)
  n_state = 0; n_state(0) = io%ntask
  call tic (time=start)
  call timer%print()
  !$omp parallel private(link,nbor,task,mytask) shared(n_state) default(shared)
  !$omp master
  do while (task_list%na>0 .and. wallclock()-start < io%job_seconds)
    start_iter = wallclock()
    !---------------------------------------------------------------------------
    ! Check if there is MPI work
    !---------------------------------------------------------------------------
    call task_list%check_mpi()
    !---------------------------------------------------------------------------
    ! Run through the task list, looking for tasks that are ready
    !---------------------------------------------------------------------------
    fell_through = .true.
    link => task_list%head
    do while (associated(link))
      task => link%task
      if (task%is_clear (bits%virtual)) then
        !-----------------------------------------------------------------------
        ! If task is updated, set it to idle; this should ONLY be done by the
        ! master thread
        !-----------------------------------------------------------------------
        if (link%task%state==4) then
          call set_state (link%task, 0)
        end if
        !-----------------------------------------------------------------------
        ! If task is idle, check if it is ready to update
        !-----------------------------------------------------------------------
        if (task%state==0) then
          if (is_ready (link)) &
            call set_state (link%task, 1)
        end if
        !-----------------------------------------------------------------------
        ! If task is ready to update, spawn a background task
        !-----------------------------------------------------------------------
        if (task%state==1) then
          fell_through = .false.
          !$omp atomic read
          task_list%nq = n_state(1)
          !$omp atomic
          task_list%nq = task_list%nq  + n_state(2)
          call set_state (link%task, 2)
          !$omp task firstprivate(task_list,link,test) default(shared)
          call set_state (link%task, 3)
          call task_list%update (link, test)
          call set_state (link%task, 4)
          !$omp end task
        end if
      end if
      link => link%next
    end do
    !---------------------------------------------------------------------------
    ! Start a thread doing I/O, if buffers from all tasks pending
    !---------------------------------------------------------------------------
    if (mpi_io%iwrite_list%n == io%nwrite) then
      !$omp task default(shared)
      call mpi_io%iwrite_list%check()
      !$omp end task
    end if
    call write_state (io_unit%queue, '1')
    if (verbose>1) call print_state ('1')
    !---------------------------------------------------------------------------
    ! Diagnostics
    !---------------------------------------------------------------------------
    if (fell_through) then
      if (verbose>3) print *,'dispatcher3: fell through'
      wc_failed = wc_failed + (wallclock()-start_fail - start_iter)
      n_fail = n_fail+1
    else
      if (n_fail>0 .and. verbose>0) &
        call print_state ('2')
      n_fail = 0
    end if
  end do
  !$omp end master
  !$omp end parallel
  call timer%print()
  if (omp%master) &
    write(io%output,'(a,i6,4(a,f7.4))') &
      ' dispatcher3_t%finalize: rank =',mpi%rank, &
      ', fraction of ready tasks =', real(n_ready)/real(n_check), &
      ', fraction of dispatcher idle time =', wc_failed/max(1d-10,wallclock()-start)
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=start)
  call trace%end (itimer)
END SUBROUTINE execute

!===============================================================================
!===============================================================================
SUBROUTINE print_state (label)
  character(len=*):: label
  !-----------------------------------------------------------------------------
  call write_state (io_unit%output, label)
END SUBROUTINE

!===============================================================================
!===============================================================================
SUBROUTINE write_state (unit, label)
  integer:: unit
  character(len=*):: label
  !-----------------------------------------------------------------------------
  if (n_write <= 0) return
  !$omp critical (set_state_cr)
  n_write = n_write-1
  write(unit,'("DISPATCHER3: ",a,f12.6,f9.3,2x,3(2x,a,i6),9(2x,a,i5))') &
    label, wallclock(), &
    wc_failed/max(1d-10,wallclock()-start), &
    'fail:',n_fail    , &
    'idl:', n_state(0), &
    'rdy:', n_state(1), &
    'spw:', n_state(2), &
    'act:', n_state(3), &
    'upd:', n_state(4), &
    'idl:', d_state(0), &
    'rdy:', d_state(1), &
    'spw:', d_state(2), &
    'act:', d_state(3), &
    'upd:', d_state(4)
  d_state = 0
  !$omp end critical (set_state_cr)
END SUBROUTINE write_state

!===============================================================================
!> Set a new task state, while keeping track of the number of tasks in each state
!===============================================================================
SUBROUTINE set_state (task, in)
  class(task_t), pointer:: task
  integer:: in, ip
  !...........................................................................
  if (use_critical) then
    !$omp critical (set_state_cr)
    ip = task%state
    task%state = in
    n_state(ip) = n_state(ip) - 1
    n_state(in) = n_state(in) + 1
    d_state(in) = d_state(in) + 1
    if (verbose>1) &
      print '(a,i4,f12.6,2x,5i4,2x,2i3)', ' dispatcher3:', omp%thread, wallclock(), n_state, ip, in
    !$omp end critical (set_state_cr)
  else
   !$omp atomic read
   ip = task%state
   !$omp atomic
   n_state(ip) = n_state(ip) - 1
   !$omp atomic
   n_state(in) = n_state(in) + 1
   !$omp atomic
   d_state(in) = d_state(in) + 1
   !$omp atomic write
   task%state = in
   if (verbose>1) &
     write (io%output,'(a,i4,f12.6,2x,5i4,2x,2i3)') &
       ' dispatcher3:', omp%thread, wallclock(), n_state, ip, in
  end if
END SUBROUTINE set_state

!===============================================================================
!> Check if a task is ready to update, by changing if all nbor task are "ahead
!> of" the task in question
!===============================================================================
LOGICAL FUNCTION is_ready (link)
  class(link_t), pointer:: link, nbor
  integer:: state
  !...........................................................................
  !$omp atomic
  n_check = n_check+1
  call io%assert (associated(link), 'is_ready: link missing')
  call io%assert (associated(link%task), 'is_ready: link%task missing')
  !$omp atomic read
  state = link%task%state
  if (link%task%is_set (bits%virtual)) then
    is_ready = .false.
  else
    is_ready = .true.
    nbor => link%nbor
    do while (associated(nbor))
      if (nbor%needed) then
        !$omp atomic read
        state = nbor%task%state
        if (state==3) then
          is_ready = .false.
        else
          is_ready = is_ready .and. nbor%task%is_ahead_of (link%task)
        end if
        if (.not.is_ready) then
          if (verbose > 2) &
            print *, link%task%id, 'is_ready: failed  on', nbor%task%id, &
              nbor%task%time, link%task%time
          exit
        end if
      end if
      nbor => nbor%next
    end do
  end if
  if (is_ready) then
    !$omp atomic
    n_ready = n_ready+1
  end if
END FUNCTION is_ready

END MODULE dispatcher3_mod
