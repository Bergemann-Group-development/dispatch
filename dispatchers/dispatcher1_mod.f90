MODULE dispatcher1_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE omp_mod
  USE omp_timer_mod
  USE mpi_mod
  USE timer_mod
  USE link_mod
  USE mpi_mesg_mod
  USE task_list_mod
  USE bits_mod
  USE process_mod
  USE global_mod
  USE index_mod
  implicit none
  private
  type, public:: dispatcher1_t
    integer:: verbose=0
    integer:: n_spawn=0
    type(task_list_t), pointer:: task_list => null()
    type(process_t):: process
  contains
    procedure:: execute
  end type
  type(dispatcher1_t), public:: dispatcher1
  type(global_t):: average_density
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self, task_list)
  class(dispatcher1_t):: self
  type(task_list_t), pointer:: task_list
  !.............................................................................
END SUBROUTINE init

!===============================================================================
!> Simple dispatcher, which in each loop through the task list takes on the
!> first task that is ready, and then spawns other threads that handle the
!> remaining tasks that are ready.
!===============================================================================
SUBROUTINE execute (self, task_list, test)
  class(dispatcher1_t):: self
  type(task_list_t), pointer:: task_list
  logical:: test
  !.............................................................................
  class(link_t), pointer:: link, nbor
  logical:: ready, mytask
  real(8):: start, v, t
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  self%task_list => task_list
  call tic (time=start)
  !call self%process%update (task_list)
  !call average_density%init ('density', mpi%size, 5, cadence=0.01)
  call timer%print()
  !$omp parallel private(link,nbor,ready,mytask) default(shared)
  !$omp master
  mytask = .false.
  task_list%nq = 0
  do while (task_list%na>0 .and. wallclock()-start < io%job_seconds)
    call trace%begin ('dispatcher1_t%execute', itimer=itimer)
                                        if (self%verbose>2) write (io_unit%dispatcher,*) wallclock(), 'check_mpi'
    call task_list%check_mpi()
    call self%process%update (task_list)
    !v = task_list%average (idx%d, t)
    !print *, 'AVER', mpi%rank, v, t
    !average_density%order = min (3, task_list%head%task%istep)
    !call average_density%update (v, 1d0, t)
    link => task_list%head
    do while (associated(link))
      ready = .false.
      if (link%task%state>2) then
                                        if (self%verbose>3) write (io_unit%dispatcher,*) wallclock(), link%task%id, 'busy'
      else if (link%task%is_set (bits%internal) .or. link%task%is_set (bits%boundary)) then
        ready = .true.
        nbor => link%nbor
        do while (associated(nbor))
          if (io%verbose > 2) &
            print *, 'dispatcher1: task, nbor, is_ahead =', link%task%id, nbor%task%id, nbor%task%is_ahead_of (link%task)
          ready = ready .and. nbor%task%is_ahead_of (link%task)
          if (.not.ready) then
                                        if (self%verbose>2) write (io_unit%dispatcher,*) wallclock(), link%task%id, 'failed', nbor%task%id, nbor%task%time, link%task%time
            exit
          end if
          nbor => nbor%next
        end do
        if (io%verbose > 2) &
          print *, 'dispatcher1: task, ready =', link%task%id, ready
        !-----------------------------------------------------------------------
        ! If a task is ready, change its state to 2, so we can keep track of nq
        !-----------------------------------------------------------------------
        if (ready) then
          mpi_mesg%n_ready = mpi_mesg%n_ready+1
          if (link%task%state<2) then
            !$omp atomic write
            link%task%state = 2
            !$omp atomic
            task_list%nq = task_list%nq+1
          end if
          !---------------------------------------------------------------------
          ! The master thread takes on at least one task update per loop
          !---------------------------------------------------------------------
          if (mytask) then
                                        if (self%verbose>2) write (io_unit%dispatcher,*) wallclock(), link%task%id, 'my task'
            call task_list%check_mpi()
                                        if (self%verbose>2) write (io_unit%dispatcher,*) wallclock(), link%task%id, 'check task_list'
                                        if (self%verbose>1) write (io_unit%log,*) wallclock(), link%task%id, 'update start'
            if (io%verbose > 1) &
              print *, omp%thread, 'dispatche1: updating  task', link%task%id
            call task_list%update (link, test)
            !$omp atomic
            task_list%nq = task_list%nq-1
            !$omp atomic write
            link%task%state = 1
                                        if (self%verbose>1) write (io_unit%log,*) wallclock(), link%task%id, 'update end'
            mytask = .false.
          !---------------------------------------------------------------------
          ! If multi-threaded, spawn background tasks up to OMP_NUM_THREADS + 1
          !---------------------------------------------------------------------
          else if (omp%nthreads>2) then
            if (self%n_spawn <= omp%nthreads+1) then
                                        if (self%verbose>1) write (io_unit%dispatcher,*) wallclock(), link%task%id, 'spawning'
                                        if (self%verbose>2) write (io_unit%log,*) wallclock(), link%task%id, 'spawn begin'
              !$omp atomic write
              link%task%state = 3
              !$omp atomic
              self%n_spawn = self%n_spawn + 1
              if (io%verbose > 2) &
                print *, omp%thread, 'dispatche1: spawning task', link%task%id
              !$omp task default(shared) firstprivate(link)
              call self%task_list%check_mpi ()
                                        if (self%verbose>2) write (io_unit%log,*) wallclock(), link%task%id, 'check task_list'
                                        if (self%verbose>1) write (io_unit%log,*) wallclock(), link%task%id, 'update start', self%n_spawn
              if (io%verbose > 1) &
                print *, omp%thread, 'dispatche1: updating  task', link%task%id
              call task_list%update (link, test)
              !$omp atomic write
              link%task%state = 1
              !$omp atomic
              self%n_spawn = self%n_spawn - 1
                                        if (self%verbose>1) write (io_unit%log,*) wallclock(), link%task%id, 'update end', self%n_spawn
              !$omp end task
              !$omp atomic
              task_list%nq = task_list%nq-1
                                        if (self%verbose>2) write (io_unit%log,*) wallclock(), link%task%id, 'spawn end', self%n_spawn
            else
                                        if (self%verbose>2) write (io_unit%log,*) wallclock(), link%task%id, 'taskyield'
              !$omp taskyield
            end if
          !---------------------------------------------------------------------
          ! Single-threaded: just update
          !---------------------------------------------------------------------
          else
                                        if (self%verbose>2) write (io_unit%dispatcher,*) wallclock(), link%task%id, 'updating'
                                        if (self%verbose>2) write (io_unit%log,*) wallclock(), link%task%id, 'update start'
            call task_list%update (link, test)
            link%task%state = 1
            task_list%nq = task_list%nq-1
                                        if (self%verbose>2) write (io_unit%log,*) wallclock(), link%task%id, 'update end'
          end if
        end if
      end if
      link => link%next
    end do
                                        if (self%verbose>0) call flush(io_unit%dispatcher)
    mytask = .true.
    call trace%end (itimer)
  end do
  !$omp end master
  !$omp end parallel
  call timer%print()
  call mpi_mesg%diagnostics(1)
  call toc ('wall time', timer%n_update, time=start)
END SUBROUTINE execute

END MODULE dispatcher1_mod
