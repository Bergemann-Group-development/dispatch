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
  class(task_t), pointer:: mhd, task
  class(rt_t), pointer:: omega
  class(patch_t), pointer:: patch
  class(task_list_t), pointer:: rt_task_list
  integer:: i_omega, ntask, iter, iostat
  integer, save:: verbose=0, n_warmup=5
  namelist /rt_setup_params/ n_warmup
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_setup')
  !-----------------------------------------------------------------------------
  rewind (io_unit%input)
  read (io_unit%input, rt_setup_params, iostat=iostat)
  if (io%master) write (io_unit%output, rt_setup_params)
  !-----------------------------------------------------------------------------
  allocate (rt_task_list)
  call rt_task_list%init
  if (verbose > 1) &
    call task_list%print ('task_list')
  link => task_list%head
  do while (associated(link))
    call io%assert (associated(link%task), 'link%task missing')
    mhd => link%task
    if (verbose > 0) & 
      print *,'processing MHD task', mhd%id
    select type (mhd)
    class is (rt_t)
      call io%assert (associated(mhd%link%task), 'mhd%link%task missing')
      !-------------------------------------------------------------------------
      ! Add nbor lists to the sub-tasks
      !-------------------------------------------------------------------------
      call mhd%init_nbors
      call mhd%rt_eos
      nbor => link%nbor
      do while (associated(nbor))
        call io%assert (associated(nbor%task), 'nbor%task missing')
        call io%assert (associated(nbor%link), 'nbor%link missing')
        ! -- do not download RT values to MHD parent tasks
        nbor%download = .true.
        ! -- call check_ready() on MHD task when RT task updates
        nbor%needs_me = .true.
        ! -- include MHD in check_ready() of RT task updates
        nbor%needed = .true.
        nbor => nbor%next
      end do
      !-------------------------------------------------------------------------
      ! Prepend subtasks to the MHD task nbor list, prepend the MHD task to the
      ! subtask nbor list, and finally append the sub-task link to the task list
      !-------------------------------------------------------------------------
      do i_omega=1,mhd%n_omega
        omega => mhd%omega(i_omega)
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
        !-----------------------------------------------------------------------
        call omega%clear (bits%ready)
        call rt_task_list%append_link (omega%link)
      end do
    end select
    if (verbose > 1) &
      call link%print_nbors
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! EOS for warm up
  !-----------------------------------------------------------------------------
   print *, 'EOS for warm-up'
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (rt_t)
      call task%rt_eos
    end select
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Merge task lists
  !-----------------------------------------------------------------------------
  patch => omega
  call data_io%register_kind (patch)
  if (verbose > 1) &
    call rt_task_list%print ('rt_task_list')
  call task_list%append_task_list (rt_task_list)
  !-----------------------------------------------------------------------------
  ! Warm up the RT
  !-----------------------------------------------------------------------------
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
