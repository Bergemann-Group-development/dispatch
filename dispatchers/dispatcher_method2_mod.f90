!> =============================================================================
!> A simple task dispatcher, which handles a list of tasks that may be started
!> as background OpenMP tasks.  A task can be in any of 4 states:
!>
!> 1: ready for updating
!> 2: undergoing updating by some thread
!> 3: has been updated, but nbors not yet checked
!> 4: finished
!>
!> When starting, all tasks are normally in state 1, and the dispatcher starts
!> a number (omp%nthreads) of background tasks that update them, first changing
!> the state to 2.  When finished, the updating thread changes the state from 2
!> to 3. When the dispatcher finds a task in state 3, it runs a check on the 
!> nbor tasks
!>
!> dispatcher%execute
!>   dispatcher%update
!>     link%check_task
!>       dispatcher%activate
!>         link%check_mpi
!>         link%updater
!> =============================================================================
MODULE dispatcher_method2_mod
  USE mpi_mesg_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE link_mod
  USE list_mod
  USE task_list_mod
  USE task_mod
  USE bits_mod
  USE patch_mod
  USE experiment_mod
  implicit none
  private
  type, public :: dispatcher_method2_t
    integer   :: running
    type(task_list_t), pointer:: task_list => null()
  contains
    procedure :: execute
    procedure :: activate
    procedure :: update
  end type
  type(task_list_t), save, pointer:: task_list => null()
  type(dispatcher_method2_t), save, public:: dispatcher_method2
CONTAINS

!> =============================================================================
!> Initialize the list, by populating it from a task_list
!> =============================================================================
SUBROUTINE execute (self, list)
  class(dispatcher_method2_t)    :: self
  type(task_list_t), pointer :: list
  class(link_t), pointer :: link
  !-----------------------------------------------------------------------------
  call trace_begin ('dispatcher_method2_t%execute')
  call list%startup
  self%running = 0
  task_list => list
  list%nq = list%n
  list%na = list%n
  link => list%head
  do while (associated(link))
    call link%task%set (bits%ready)
    link => link%next
  end do
  !$omp parallel
  !$omp master
  do while (task_list%na > 0 .and. wallclock() < task_list%job_seconds)
    call self%update
  end do
  !$omp end master
  !$omp end parallel
  call trace_end
  call timer%print
  call mpi_mesg%diagnostics(1)
END SUBROUTINE execute

!> =============================================================================
!> Activate a task pointed to by a link, by starting a task%update as a
!> background OpenMP task.  Set the state of task to 1 when dormant, 2 when
!> ready, and 3 when active.
!> =============================================================================
SUBROUTINE activate (self, link)
  class(dispatcher_method2_t)    :: self
  class(link_t), pointer :: link
  !-----------------------------------------------------------------------------
  if (omp%nthreads==1) then
    link%task%state = 3
    call task_list%updater (link)
    link%task%state = 1
  else
    !$omp task default(shared) firstprivate(link)
    !$omp atomic write
    link%task%state = 3
    call task_list%check_mpi ()
    call task_list%updater (link)
    !$omp atomic write
    link%task%state = 1
    !$omp atomic
    task_list%nq = task_list%nq - 1
    !$omp end task
  end if
END SUBROUTINE activate

!> =============================================================================
!> Update the list of running tasks.  For each task on the running list, if it
!> is no longer busy, check if any of its passive nbors has become executable.
!> Finally, check the task itself.
!> =============================================================================
SUBROUTINE update (self)
  class(dispatcher_method2_t)     :: self
  class(link_t), pointer  :: link, nbor
  integer                 :: state
  logical                 :: ok
  integer, save           :: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin ('dispatcher_method2_t%update', itimer=itimer)
  link => task_list%head
  do while (associated(link))
    if (io%debug(2)) then
      print *, 'checking link', link%task%id, link%task%state
    end if
    !$omp atomic read
    state = link%task%state
    if (state==4) then
      link => link%next
      cycle
    end if
    if (link%task%time > io%end_time) then
      !$omp atomic
      task_list%na = task_list%na - 1
      !$omp atomic write
      link%task%state = 4
      cycle
    end if
    if (state<=1) then
      nbor => link%nbor
      do while (associated(nbor))
        !$omp atomic read
        state = nbor%link%task%state
        if (state==1) then
          if (io%debug(3)) then
            print '(i3,f12.6,a,i6,i3)', omp%thread, wallclock(), '       checking nbor', nbor%link%task%id, state
          end if
          call check_task (nbor%link, ok)
        end if
        nbor => nbor%next
      end do
      call check_task (link, ok)
    end if
    link => link%next
  end do
  call trace_end (itimer)
contains
  !-----------------------------------------------------------------------------
  ! If all nbors of a task are ahead in time activate the task.  Only the master
  ! thread is checking, so there should be no chance of starting multiple 
  ! threads on the same update. The only state change performed by the many
  ! production threads is to go from state 2 to state 3, which is of no concern
  ! to the master thread, and then from state 3 to state 1, which will
  ! eventually be detected (typically much later) by the master thread. The
  ! delay before checking is not a disadvantage, since it actually increases
  ! the chances that the patch is ready for updating when it checked. At the end
  ! a change to state 4 indicates a finished task.
  !-----------------------------------------------------------------------------
  subroutine check_task (link, ok)
    class(link_t), pointer:: link, nbor
    class(task_t), pointer:: task
    logical:: ok
    integer:: state
    !---------------------------------------------------------------------------
    ok = .true.
    task => link%task
    nbor => link%nbor
    do while (ok.and.associated (nbor))
      ok = ok .and. nbor%task%is_ahead_of(task)
      nbor => nbor%next
    end do
    if (ok) then
      !$omp atomic read
      state = link%task%state
      if (state<=1) then
        if (state==1) then
          !$omp atomic
          task_list%nq = task_list%nq + 1
        end if
        !$omp atomic write
        link%task%state = 2
        call self%activate (link)
      else
         print '(i3,f12.6,a,i6,i3)', omp%thread, wallclock(), ' found unexpected state', link%task%id, link%task%state
      end if
    end if
  end subroutine check_task
END SUBROUTINE update

END MODULE dispatcher_method2_mod
