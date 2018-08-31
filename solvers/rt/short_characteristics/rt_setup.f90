!> =============================================================================
!> Add RT subtasks into the general task_list.  Each subtask is added to the
!> nbor list of the parent task, and the parent task should be added to the nbor
!> list of the subtask.
!>
!> Reminder about links and nbors:   A task should be represented on the task_list
!> by a link, which has %next and %prev pointers that are not used for any other
!> purpose.   The link must contain an %nbor pointer, which is the head of a chain
!> of links whose tasks are to provide guard zones.  The main and nbor links have
!> the relations
!>
!> link%nbor%task => nbtask
!> nbtask%link%task => nbtask
!>
!> The nbor lists have special behavior, in that even though the parent RT task
!> needs to check that the sub-tasks are up-to-date, it does not need to, and
!> should not, use a standard download_link call visavis the RT sub-tasks
!> =============================================================================
SUBROUTINE rt_setup (task_list)
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  USE task_list_mod
  USE list_mod
  USE link_mod
  USE task_mod
  USE patch_mod
  USE trace_mod
  USE rt_mod
  USE bits_mod
  USE data_io_mod
  USE download_mod
  implicit none
  type(task_list_t):: task_list
  !.............................................................................
  class(link_t), pointer:: link, nbor
  class(task_t), pointer:: task
  class(rt_t), pointer:: omega=>null(), mhd
  class(patch_t), pointer:: patch
  class(task_list_t), pointer:: rt_task_list
  integer:: i_omega, ntask, iter, iostat, id_save, n_mhd, n_rt
  integer, save:: verbose=0, n_warmup=5, n_omega=0
  real(8):: boundary_distance
  logical, save:: unique_ids=.false.
  logical, save:: only_vertical=.false.
  namelist /rt_setup_params/ n_warmup, verbose, unique_ids, only_vertical
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_setup')
  call io%header ('begin rt_setup')
  !-----------------------------------------------------------------------------
  rewind (io_unit%input)
  read (io_unit%input, rt_setup_params, iostat=iostat)
  if (mpi%size > 1 .and. .not.unique_ids) then
    print *, 'WARNING: unique_ids set to .true. when using MPI size =', mpi%size
    unique_ids = .true.
  end if
  if (io%master) write (io_unit%output, rt_setup_params)
  !-----------------------------------------------------------------------------
  allocate (rt_task_list)
  call rt_task_list%init
  if (verbose > 1) &
    call task_list%print ('task_list')
  n_rt = 0
  n_mhd = 0
  link => task_list%head
  do while (associated(link))
    call io%assert (associated(link%task), 'link%task missing')
    task => link%task
    if (verbose > 0) & 
      print *,'processing MHD task', task%id
    select type (task)
    class is (rt_t)
      mhd => task
      call io%assert (associated(mhd%link%task), 'mhd%link%task missing')
      !-------------------------------------------------------------------------
      ! Measure distance from lower or upper boundary, in patch size units
      !-------------------------------------------------------------------------
      boundary_distance = min(mhd%position(3)-task_list%lc(3), &
                              task_list%uc(3)-mhd%position(3))/mhd%size(3)
      !-------------------------------------------------------------------------
      ! Add nbor lists to the sub-tasks
      !-------------------------------------------------------------------------
      call mhd%init_nbors
      !call mhd%rt_eos
      nbor => link%nbor
      do while (associated(nbor))
        call io%assert (associated(nbor%task), 'nbor%task missing')
        call io%assert (associated(nbor%link), 'nbor%link missing')
        ! -- download MHD nbor values to the MHD parent tasks
        nbor%download = .true.
        ! -- call check_ready() on MHD nbor task when MHD task updates
        nbor%needs_me = .true.
        ! -- include MHD nbor in check_ready() of MHD task updates
        nbor%needed = .true.
        nbor => nbor%next
      end do
      if (mhd%is_clear(bits%no_rt)) then
        !-----------------------------------------------------------------------
        ! Prepend subtasks to the MHD task nbor list, prepend the MHD task to the
        ! subtask nbor list, and finally append the sub-task link to the task list
        !-----------------------------------------------------------------------
        n_omega = mhd%n_omega
        do i_omega=1,mhd%n_omega
          omega => mhd%omega(i_omega)
          !---------------------------------------------------------------------
          ! Reassign unique IDs to RT tasks
          !---------------------------------------------------------------------
          if (unique_ids) then
            id_save = omega%id
            omega%id = 0
            call omega%task_t%init
          end if
          omega%parent => mhd
          if (verbose > 1) &
            print *,'appending task', omega%id, omega%nw, omega%it, omega%new
          call io%assert (associated(omega%link), 'omega%link missing')
          call prepend_task (link%nbor, omega)
          call io%assert (associated(mhd%link), 'mhd%link missing')
          call io%assert (associated(mhd%link%task), 'mhd%link%task missing')
          call prepend_task (omega%link%nbor, mhd)
          ! -- do not download RT values to MHD task
          link%nbor%download = .false.
          ! -- do not download MHD values to RT task
          omega%link%nbor%download = .false.
          ! -- call check_ready() on RT task when MHD task updates
          omega%link%nbor%needs_me = .true.
          ! -- include RT sub-task in check_ready() of MHD task
          omega%link%nbor%needed = .true.
          !omega%rt_next = omega%time + (1.0-boundary_distance*0.5)*mhd%rt_time
          omega%rt_next = omega%time + mhd%rt_time
          !---------------------------------------------------------------------
          call omega%clear (bits%ready)
          call rt_task_list%append_link (omega%link)
        end do
        n_rt = n_rt+1
      else
        n_mhd = n_mhd+1
      end if
    end select
    if (verbose > 1) &
      call link%print_nbors
    link => link%next
  end do
  print *, 'rt_setup: n_rt, n_mhd =', n_rt, n_mhd
  if (associated(omega)) then
    patch => omega
    call data_io%register_kind (patch)
  end if
  !-----------------------------------------------------------------------------
  ! Run through all nbor tasks of all RT tasks
  !-----------------------------------------------------------------------------
  if (only_vertical) then
    link => rt_task_list%head
    do while (associated(link))
      nbor => link%nbor
      do while (associated(nbor))
        !-------------------------------------------------------------------------
        ! If the nbor task has a different horizontal position, do not check it's
        ! time -- accepting a possible delay of information sideways.
        !-------------------------------------------------------------------------
        if (any(link%task%position(1:2) /= nbor%task%position(1:2))) then
          nbor%needed = .false.
        end if
        nbor => nbor%next
      end do
      link => link%next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Reassign unique IDs to RT tasks
  !-----------------------------------------------------------------------------
  if (unique_ids) then
    link => task_list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (rt_t)
        mhd => task
        if (mhd%is_clear(bits%no_rt)) then
          do i_omega=1,n_omega
            omega => mhd%omega(i_omega)
            omega%id = 0
            call omega%task_t%init
          end do
        end if
      end select
      link => link%next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! EOS for warm up
  !-----------------------------------------------------------------------------
  print *, 'EOS for warm-up'
  !$omp parallel default (shared) private(link, task)
  !$omp single
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (rt_t)
      mhd => task
      if (task%is_clear(bits%no_rt)) then
        !-----------------------------------------------------------------------
        ! Spawn OMP tasks to compute EOS
        !-----------------------------------------------------------------------
        !$omp task default(shared) firstprivate(mhd)
          call mhd%rt_eos
        !$omp end task
      end if
    end select
    link => link%next
  end do
  !$omp end single
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! Merge task lists, if RT is present
  !-----------------------------------------------------------------------------
  if (n_omega>0) then
    if (verbose > 1) &
      call rt_task_list%print ('rt_task_list')
    call task_list%append_task_list (rt_task_list)
    !---------------------------------------------------------------------------
    ! Warm up the RT
    !---------------------------------------------------------------------------
    do iter=1,n_warmup
      print *, 'warm-up, iteration', iter
      link => task_list%head
      do while (associated(link))
        task => link%task
        select type (task)
        class is (rt_t)
        if (task%i_omega > 0) then
          call download%download_link (link)
          call task%solve
        end if
        end select
        link => link%next
      end do
    end do
  end if
  !-----------------------------------------------------------------------------
  if (verbose > 0) &
    print *, 'rt_setup: io%ntask =', io%ntask
  if (verbose > 0) &
    call task_list%print ('new task_list')
  call trace%end
contains

!> =============================================================================
!> =============================================================================
subroutine prepend_task (nbor, mhd)
  class(link_t), pointer:: nbor
  class(rt_t), pointer:: mhd
  class(patch_t), pointer:: patch
  class(link_t), pointer:: new
  !.............................................................................
  allocate (new)
  patch => mhd
  new%task => patch
  new%next => nbor
  new%link => patch%link
  nbor => new
end subroutine prepend_task

END SUBROUTINE
