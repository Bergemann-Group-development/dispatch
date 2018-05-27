MODULE dispatcher_mod
  USE io_unit_mod
  USE trace_mod
  USE omp_mod
  USE link_mod
  USE omp_timer_mod
  USE task_list_mod
  implicit none
  private
  type, public:: dispatcher_t
  contains
    procedure, nopass:: execute
  end type
  integer, save:: verbose=2
  type(dispatcher_t), public:: dispatcher
CONTAINS

!> =============================================================================
!> Simple dispatcher, which in each loop through the task list takes on the
!> first task that is ready, and then spawns other threads that handle the
!> remaining tasks that are ready.
!> =============================================================================
SUBROUTINE execute (task_list)
  type(task_list_t):: task_list
  class(link_t), pointer:: link, nbor
  logical:: ready, mytask
  real(8):: start
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  start = wallclock()
  write (io_unit%dispatcher,*) wallclock(), 'dispatcher_t%execute'
  !$omp parallel
  !$omp master
  mytask = .false.
  do while (task_list%na>0 .and. wallclock()-start < task_list%job_seconds)
    call trace%begin ('dispatcher_t%execute', itimer=itimer)
    if (verbose>2) &
      write (io_unit%dispatcher,*) wallclock(), 'check_mpi'
    call task_list%check_mpi()
    if (verbose>2) &
      write (io_unit%dispatcher,*) wallclock(), 'check task_list'
    link => task_list%head
    do while (associated(link))
      if (link%task%state==2) then
        if (verbose>3) &
          write (io_unit%dispatcher,*) wallclock(), link%task%id, 'busy'
      else
        ready = .true.
        nbor => link%nbor
        do while (associated(nbor))
          ready = ready .and. nbor%task%is_ahead_of (link%task)
          if (.not.ready) then
            if (verbose>2) &
              write (io_unit%dispatcher,*) wallclock(), link%task%id, 'failed', nbor%task%id, nbor%task%time, link%task%time
            exit
          end if
          nbor => nbor%next
        end do
        if (ready) then
          if (mytask) then
            if (verbose>1) &
              write (io_unit%dispatcher,*) wallclock(), link%task%id, 'my task'
            call task_list%updater (link)
            mytask = .false.
          else if (omp%nthreads>1) then
            if (verbose>1) &
              write (io_unit%dispatcher,*) wallclock(), link%task%id, 'spawning'
            link%task%state = 2
            !$omp task firstprivate(link)
            call task_list%updater (link)
            link%task%state = 1
            !$omp end task
          else
            if (verbose>1) &
              write (io_unit%dispatcher,*) wallclock(), link%task%id, 'updating'
            call task_list%updater (link)
          end if
        end if
      end if
      if (verbose>0) call flush(io_unit%dispatcher)
      link => link%next
    end do
    mytask = .true.
    call trace%end (itimer)
  end do
  if (verbose>0) then
    write (io_unit%dispatcher,*) wallclock(), 'end do while', task_list%na, wallclock()-start
    call flush(io_unit%dispatcher)
  end if
  !$omp end master
  !$omp end parallel
END SUBROUTINE execute

END MODULE dispatcher_mod
