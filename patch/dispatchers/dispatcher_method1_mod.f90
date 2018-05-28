MODULE dispatcher_method1_mod
  USE io_unit_mod
  USE trace_mod
  USE omp_mod
  USE link_mod
  USE omp_timer_mod
  USE task_list_mod
  implicit none
  private
  type, public:: dispatcher_method1_t
    integer:: verbose=0
    integer:: n_spawn=0
  contains
    procedure:: execute
  end type
  type(dispatcher_method1_t), public:: dispatcher_method1
CONTAINS

!> =============================================================================
!> Simple dispatcher, which in each loop through the task list takes on the
!> first task that is ready, and then spawns other threads that handle the
!> remaining tasks that are ready.
!> =============================================================================
SUBROUTINE execute (self, task_list)
  class(dispatcher_method1_t):: self
  class(task_list_t):: task_list
  class(link_t), pointer:: link, nbor
  logical:: ready, mytask
  real(8):: start
  integer, save:: itimer=0
  integer:: n
  !-----------------------------------------------------------------------------
  call task_list%init_all_nbors
  start = wallclock()
  !$omp parallel private(link,nbor,ready,mytask) default(shared)
  !$omp master
  mytask = .false.
  do while (task_list%na>0 .and. wallclock()-start < task_list%job_seconds)
    call trace%begin ('dispatcher_t%execute', itimer=itimer)
    call task_list%check_mpi()
    link => task_list%head
    do while (associated(link))
      ready = .false.
      if (link%task%state>2) then
      else
        ready = .true.
        nbor => link%nbor
        do while (associated(nbor))
          ready = ready .and. nbor%task%is_ahead_of (link%task)
          if (.not.ready) then
            exit
          end if
          nbor => nbor%next
        end do
        if (ready) then
          if (link%task%state<2) then
            !$omp atomic
            task_list%nq = task_list%nq+1
            link%task%state = 2
          end if
          if (mytask) then
            call task_list%updater (link)
            mytask = .false.
          else if (omp%nthreads>1) then
            if (self%n_spawn < omp%nthreads+3) then
              !$omp atomic
              self%n_spawn = self%n_spawn + 1
              link%task%state = 3
              !$omp task default(shared) firstprivate(link)
              call task_list%updater (link)
              !$omp atomic
              self%n_spawn = self%n_spawn - 1
              link%task%state = 1
              !$omp atomic
              task_list%nq = task_list%nq-1
              !$omp end task
            else
              !$omp taskyield
            end if
          else
            call task_list%updater (link)
          end if
        end if
      end if
      link => link%next
    end do
    mytask = .true.
    call trace%end (itimer)
  end do
  !$omp end master
  !$omp end parallel
END SUBROUTINE execute

END MODULE dispatcher_method1_mod
