MODULE mulligan_mod
!===============================================================================
!> Check if any nbor of the link task should do a mulligan.  The criterion is
!> triggered if an nbor task is leading the current task with more than a factor
!> (1.0+grace) times the maximum time step it would currently take, based on
!> the Courant condition computed from the u_phase(:,:,:) values in its guard
!> cells, owned by the link task.
!>
!===============================================================================
  USE trace_mod
  USE io_mod
  USE omp_mod
  USE task_mod
  USE bits_mod
  USE mesh_mod
  USE patch_mod
  USE link_mod
  USE list_mod
  implicit none
  private
  type, public:: mulligan_t
  contains
    procedure:: init
    procedure:: check
    procedure, nopass:: trigger
    procedure, nopass:: update
  end type
  integer, save:: verbose=0, level=0
  type(mulligan_t), public:: mulligan
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(mulligan_t):: self
  integer:: iostat
  namelist /mulligan_params/ verbose, level
  !-----------------------------------------------------------------------------
  call trace%begin ('mulligan_t%init')
  rewind (io%input)
  read(io%input, mulligan_params, iostat=iostat)
  call io%namelist_warning ('mulligan_params', iostat)
  write (io%namelist, mulligan_params)
  if (level > 1 .and. omp%nthreads > 1) then
    call io%warning ('The current Mulligan mechanism may not be entirely multi-thread-safe')
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Check if any nbor needs to do a Mulligan
!===============================================================================
SUBROUTINE check (self, task_list, link)
  class(mulligan_t):: self
  class(list_t):: task_list
  class(link_t), pointer:: link
  !.............................................................................
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task, nbtask
  real(8):: gsize(3), dist(3), dt, ds
  real:: u_max
  integer:: l(3), u(3)
  logical:: trigger
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (level <= 0) return
  call trace%begin ('mulligan_t%check', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Create a border mask
  !-----------------------------------------------------------------------------
  task => link%task
  select type (task)
  class is (patch_t)
    !$omp atomic write
    task%dt_fixed = 0d0
    ds = minval(task%ds,mask=task%mesh%n > 1)
    !---------------------------------------------------------------------------
    ! Check if the nbor%u_max calls for a do-over of the nbor task
    !---------------------------------------------------------------------------
    nbor => link%nbor
    do while (associated(nbor))
      nbtask => nbor%task
      select type (nbtask)
      class is (patch_t)
        !-----------------------------------------------------------------------
        ! Compute the index range where the nbtask guard zones overlap
        !-----------------------------------------------------------------------
        dist = nbtask%task_t%distance (task)
        gsize = nbtask%mesh%s/2d0 + nbtask%mesh%d*nbtask%mesh%ng
        where (task%n == 1)
          l = 1
          u = 1
        else where
          l =   floor((dist - gsize)/task%mesh%d + task%mesh%o)
          u = ceiling((dist + gsize)/task%mesh%d + task%mesh%o)
          l = max (l, task%mesh%li)
          u = min (u, task%mesh%ui)
        end where
        u_max = maxval(task%u_phase(l(1):u(1),l(2):u(2),l(3):u(3))) 
        dt = task%courant*ds/u_max
        !.......................................................................
        if (verbose > 2) then
          write (io_unit%log,'(1x,a,2i6,3f6.2,2(2x,3i4))') &
            'mulligan_t%nbor_check: id, nbor, dist, l, u =', &
            task%id, nbtask%id, dist/task%size, l, u
        end if
        !-----------------------------------------------------------------------
        ! The dt above is the Courant time step of the region from which the
        ! nbor task picks up its guard zones.  If the nbor task is way ahead,
        ! relative to that time step, it is likely that it advanced too far,
        ! and it should retrace its steps, with a smaller time step.  If the
        ! nbor task is actually advancing with the same time step as this task,
        ! it could possibly be up to (1+grace)*dt ahead, so we use that as the
        ! limit for "way ahead".
        !-----------------------------------------------------------------------
        trigger = (nbtask%time-task%time) > dt*(1d0+task%grace)
        if (verbose > 1) then
          write (io_unit%log,'(1x,a,1p,2(i6,e12.3),1x,3i3,3g11.2,l3)') &
            'mulligan_t%nbor_check: id, time, nbor, time, u_max, lead/dt =', &
            task%id, task%time, nbtask%id, nbtask%time, nint(dist/task%size), &
            nbtask%u_max, u_max, (nbtask%time-task%time)/dt, trigger
        end if
        if (trigger) then
          !---------------------------------------------------------------------
          ! Write diagnostics if level is 1
          !---------------------------------------------------------------------
          if (level == 1) then
            if (verbose > 0) then
              write (io_unit%log,'(1x,a,2i6,1p,4e12.2)') &
              'list_t%mulligan_check: id, nbor, time, lead/dt, dt, u_max =', &
              task%id, nbtask%id, task%time, &
              (nbtask%time-task%time)/dt, dt, nbor%u_max
            end if
          !---------------------------------------------------------------------
          ! Trigger a Mulligan if level is 2
          !---------------------------------------------------------------------
          else if (level > 1) then
            call self%trigger (task_list, task, nbtask, dt)
          end if
        end if
      end select
      nbor => nbor%next
    end do
  end select
  call trace%end (itimer)
END SUBROUTINE check

!===============================================================================
!> Setup a task to do a Mulligan, under the control of a thread that owns it,
!> just like any other task, thus avoiding all problem with task locking,
!> associated with trying to manipulate the status of a task not owned by
!> the thread.
!>
!> The nbtask could be in an idle state, a ready state, or a busy state, as
!> shown by bits%ready and bits%busy.  In practice it is extremely unlikely
!> to be in a busy state, since it was selected exactly by being too far ahead
!> in time.  The only aspect that needs to be protected against is then that
!> the task time will become updated after it is modified here.  This can be
!> prevented by changing the "rotated" state from false to true, and by 
!> testing (under lock), in task%rotate that it is not set.
!>
!> This procedure requires only single locks, first on the task meta data,
!> then on the task list, so there is no change that this will cause a deadlock.
!===============================================================================
SUBROUTINE trigger (task_list, task, nbtask, dt)
  class(list_t):: task_list
  class(patch_t):: task, nbtask
  real(8):: dt
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('mulligan_t%trigger')
  !-----------------------------------------------------------------------------
  ! Lock the task meta-data, set the mulligan_time and modify the task time,
  ! in order to put the task near the head of the queue
  !-----------------------------------------------------------------------------
  call nbtask%lock%set ('mulligan')
  nbtask%mulligan_time =  task%time + dt
  nbtask%time = task%time
  if (nbtask%is_set (bits%busy)) then
    nbtask%rotated = .true.
  end if
  !.............................................................................
  if (verbose > 1) then
    write (io_unit%log,'(a,l4)') &
    'mulligan_t%trigger: busy =', nbtask%is_set (bits%busy)
  end if
  !.............................................................................
  call nbtask%lock%unset ('mulligan')
  !-----------------------------------------------------------------------------
  ! Make sure the task is near the head of the ready queue, by locking the task
  ! list, removing the task from the queue if it is there, and adding it back at
  ! the new, older time.  If the task should happen to be in the busy state, all
  ! that will happen is that when the task is checked in check_ready, it will
  ! be found to already have the ready bit set, and it will remain in the queue.
  !-----------------------------------------------------------------------------
  call task_list%lock%set ('mulligan')
  if (nbtask%is_set (bits%ready)) then
    call task_list%remove_from_queue (nbtask%link)
    !...........................................................................
    if (verbose > 1) then
      write (io_unit%log,'(a,l4)') &
      'mulligan_t%trigger: remove_from_queue nbtask =', nbtask%id
    end if
    !...........................................................................
  end if
  call task_list%queue_by_time (nbtask%link)
  !-----------------------------------------------------------------------------
  call task_list%lock%unset ('mulligan')
  call trace%end()
END SUBROUTINE trigger

!===============================================================================
!> Do-over on a task; this is signalled by task%mulligan_time being non-zero,
!> and results in a modified task%update, implemented here.
!===============================================================================
SUBROUTINE update (task_list, link)
  class(list_t):: task_list
  class(link_t), pointer:: link
  !.............................................................................
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer:: n
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('mulligan_t%update', itimer=itimer)
  !-----------------------------------------------------------------------------
  task => link%task
  if (level > 1 .and. task%mulligan_time > 0d0) then
    select type (task)
    class is (patch_t)
      !-------------------------------------------------------------------------
      ! De-rotate the task, effectively backing up one or more timesteps, until
      ! a time that is behind the current task%mulligan_time
      !-------------------------------------------------------------------------
      call task%lock%set ('mulligan')
      task%time = task%t(task%it)
      n = 0
      do while (task%time > task%mulligan_time)
        call task%derotate ()
        n = n+1
        !.......................................................................
        if (verbose > 0) then
          write (io_unit%log,1) 'mulligan_t%update: id,  time =', &
          task%id, n, task%time, task%mulligan_time
          1 format(a,i6,i3,1p,3e12.3)
        end if
      end do
      !-------------------------------------------------------------------------
      ! Set a fixed timestep, which brings the task to the same time as the patch
      !-------------------------------------------------------------------------
      task%dt_fixed = task%mulligan_time - task%time
      task%mulligan_time = 0d0
      !.........................................................................
      if (verbose > 0) then
        write (io_unit%log,1) &
        'mulligan_t%update: id, n, time, dtime =', &
         task%id, n, task%time, task%dt_fixed
      end if
      if (verbose > 1) then
        write (stdout,1) &
        'mulligan_t%update: id, n, time, dtime =', &
         task%id, n, task%time, task%dt_fixed
      end if
      call task%lock%unset ('mulligan')
      !-------------------------------------------------------------------------
      ! Cancel any nbors in the ready queue; they will be checked again
      !-------------------------------------------------------------------------
      call task_list%lock%set ('mulligan_t%update')
      nbor => link%nbor
      do while (associated(nbor))
        if (nbor%task%is_set (bits%ready)) then
          !call task_list%remove_from_queue (nbor%link)
          write (io_unit%log,*) 'should be removed:', nbor%task%id, nbor%task%time
        end if
        nbor => nbor%next
      end do
      call task_list%lock%unset ('mulligan_t%update')
    end select
  end if
  call trace%end (itimer)
END SUBROUTINE update 

END MODULE mulligan_mod
