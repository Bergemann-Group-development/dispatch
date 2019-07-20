!===============================================================================
!> Module handling setup of radiative transfer tasks, adding them to the
!> already existing task_list.
!===============================================================================
MODULE rt_setup_mod
  USE io_mod
  USE io_unit_mod
  USE omp_mod
  USE mpi_mod
  USE task_list_mod
  USE link_mod
  USE task_mod
  USE patch_mod
  USE solver_mod
  USE trace_mod
  USE rt_mod
  USE bits_mod
  USE data_io_mod
  implicit none
  private
  type, public:: rt_setup_t
  contains
    procedure:: init
    procedure:: init_rt_nbors_sc
    procedure:: init_rt_nbors_hybrid
  end type
  type(rt_setup_t), public:: rt_setup
CONTAINS

!===============================================================================
!> Add radiative transfer (RT) tasks to a new rt_task_list, and then append it
!> to an incoming task_list.  Each MHD task that should do RT adds 1+n_omega
!> tasks.  The n_omega sub-tasks each handle radiation in one direction (given
!> by mu=cos(theta) and phi), while the main task sums up the contributions from
!> the n_omega sub-tasks.
!>
!> The order of execution and interdependence of tasks along ray directions is
!> handled by having the main RT task wait until all sub-tasks are finished,
!> and then 1) sum up the cooling/heating, and 2) compute new EOS data.  To be
!> able to do that it also needs to wait for the MHD task to pass by in time,
!> so the density and internal energy values can be interpolated to the new RT
!> time.
!>
!> The new cooling/heating rate is thus used by the MHD tasks until the next
!> dedicated RT time, which (for now) is a contant time rt_timestep ahead
!> of the previous RT time.  The MHD tasks are thus using a constant cooling/
!> heating rate, for a time interval = rt_timestep.  By saving time slices of
!> the cooling/heating one could improve on this, by extrapolating forward in
!> time.
!===============================================================================
SUBROUTINE init (self, task_list)
  class (rt_setup_t):: self
  type(task_list_t):: task_list
  !.............................................................................
  class(link_t), pointer:: link
  class(task_t), pointer:: task, rt_task
  class(patch_t), pointer:: patch
  class(rt_t), pointer:: omega, rt
  class(task_list_t), pointer:: rt_task_list
  integer:: i_omega, iostat, n_mhd, n_rt, id
  integer, save:: verbose=0
  logical:: on = .false.
  namelist /rt_setup_params/ verbose, on
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_setup')
  call io%header ('begin rt_setup')
  !-----------------------------------------------------------------------------
  rewind (io_unit%input)
  read (io_unit%input, rt_setup_params, iostat=iostat)
  if (io%master) write (io_unit%output, rt_setup_params)
  if (.not.on) then
    call trace%end
    return
  end if
  !-----------------------------------------------------------------------------
  ! Allocate and initialize an RT task list
  !-----------------------------------------------------------------------------
  allocate (rt_task_list)
  call rt_task_list%init ('RT_list')
  !-----------------------------------------------------------------------------
  ! Add the (1+n_omega) RT tasks that belong to each MHD task with RT
  !-----------------------------------------------------------------------------
  n_rt = 0
  n_mhd = 0
  !-----------------------------------------------------------------------------
  ! Add the (1+n_omega) RT tasks that belong to each MHD task with RT
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    call io%assert (associated(link%task), 'link%task missing')
    task => link%task
    if (verbose > 0) & 
      write(io%output,*) 'processing MHD task', task%id
    select type (task)
    class is (solver_t)
      if (task%is_clear(bits%no_rt)) then
        call io%assert (associated(task%link%task), 'mhd%link%task missing')
        !-----------------------------------------------------------------------
        ! task%rt was allocated via solver_t%init, which calls extras_t%pre_init
        !-----------------------------------------------------------------------
        call task%rt%init
        !-----------------------------------------------------------------------
        ! Add the main RT task to the rt_task_list
        !-----------------------------------------------------------------------
        call rt_task_list%append_link (task%rt%link)
        n_mhd = n_mhd + 1
        n_rt = n_rt + 1
        !-----------------------------------------------------------------------
        ! Add the n_omega RT sub-tasks to the rt_task_list
        !-----------------------------------------------------------------------
        do i_omega=1,task%rt%n_omega
          call rt_task_list%append_link (task%rt%omega(i_omega)%link)
          n_rt = n_rt + 1
        end do
        omega => task%rt%omega(1)
        rt => task%rt
      end if
    end select
    link => link%next
  end do
  write (stdout,*) 'rt_setup: n_mhd, n_rt =', n_mhd, n_rt
  if (associated(omega)) then
    patch => omega
    call data_io%register_kind(patch)
  else if (associated(rt)) then
    patch => rt
    call data_io%register_kind(patch)
  end if
  !-----------------------------------------------------------------------------
  ! Add nbor lists (in a separate loop so all rt tasks are initialized first)
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (solver_t)
      if (task%rt%on) then
        if (trim(task%rt%kind)=='rt_sc') then
          call self%init_rt_nbors_sc (task)
        else if (trim(task%rt%kind)=='rt_sc') then
          call self%init_rt_nbors_hybrid (task)
        else
          call io%abort("rt_setup_mod:: unknown RT solver kind")
        end if
      end if
    end select
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Merge task lists
  !-----------------------------------------------------------------------------
  if (verbose > 1) &
    call rt_task_list%print ('rt_task_list')
  call task_list%append_task_list (rt_task_list)
  if (verbose > 0) &
    call task_list%print ('new task_list')
  !-----------------------------------------------------------------------------
  ! Reset task_list status to re-establish nbor relations, checking consistency
  !-----------------------------------------------------------------------------
  call task_list%reset_status (check=.true.)
  !-----------------------------------------------------------------------------
  ! Add nbor lists sorted by level
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    call link%sort_nbors_by_level(link%nbors_by_level,link%nbors_by_level)
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  call trace%end
END SUBROUTINE init

!===============================================================================
!> In the new, synced RT/MHD arrangement the steps and relations needed are:
!>
!> 1) The omega0 task needs the MHD nbors to dnload and call EOS for RT
!> 2) Each omega task needs the upstream omega task
!> 3) Each omega task neds the omega0 task
!> 4) The MHD task needs all omega task, while not vice versa.
!>
!> The upstream nbors should have the settings needed=T, needs_me=F, download=T,
!> while the downstream  tasks should have needed=F, needs_me=T, download=F.
!===============================================================================
SUBROUTINE init_rt_nbors_sc (self, mhd)
  class (rt_setup_t):: self
  class (solver_t):: mhd
  !.............................................................................
  integer:: i_omega
  class(rt_t), pointer:: omega0, omega, upstream
  class(link_t), pointer:: nbor
  class(task_t), pointer:: nb_task
  logical:: needs, needs_me, download
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_setup_t%init_nbors')
  omega0 => mhd%rt
  call omega0%clear(bits%init_nbors)
  !-----------------------------------------------------------------------------
  ! 1) The omega0 task needs the MHD nbors, and needs download
  !-----------------------------------------------------------------------------
  needs=.true.; needs_me=.false.; download=.true.
  call init_nbor_pair (omega0, needs, mhd, needs_me, download)
  !-----------------------------------------------------------------------------
  ! Loop over all ray directions, and add nbors and nbor relations. The only
  ! downloads that are needed are from upstream omega tasks
  !-----------------------------------------------------------------------------
  do i_omega=1,omega0%n_omega
    omega => omega0%omega(i_omega)
    !---------------------------------------------------------------------------
    ! Loop over all MHD nbors
    !---------------------------------------------------------------------------
    nbor => mhd%link%nbor
    do while (associated(nbor))
      !-------------------------------------------------------------------------
      ! Cast the nbors to solver_t type
      !-------------------------------------------------------------------------
      nb_task => nbor%task
      select type (nb_task)
      class is (solver_t)
      !-------------------------------------------------------------------------
      ! Select only omega tasks
      !-------------------------------------------------------------------------
      if (nb_task%rt%n_omega > 0) then
        !-----------------------------------------------------------------------
        ! The omega0%needs() function tests for upstream nbor omega
        !-----------------------------------------------------------------------
        upstream => nb_task%rt%omega(i_omega)
        if (omega%needs (omega, upstream) .and. upstream%on) then
          !---------------------------------------------------------------------
          ! 2) The omega task needs the upstream task, and should download from 
          ! it, but the upstream task does not need the downstream omega task.
          !---------------------------------------------------------------------
          needs=.true.; needs_me=.false.; download=.true.
          call init_nbor_pair (omega, needs, upstream, needs_me, download)
        end if
      end if
      end select
      nbor => nbor%next
    end do
    !---------------------------------------------------------------------------
    ! 3) The omega task needs the RT task, which does not need the omega tasks
    !---------------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.false.
    call init_nbor_pair (omega, needs, omega0, needs_me, download)
    !---------------------------------------------------------------------------
    ! 4) The MHD task needs the omega tasks, which do not need the MHD task
    !---------------------------------------------------------------------------
    needs=.true.; needs_me=.false.; download=.false.
    call init_nbor_pair (mhd, needs, omega, needs_me, download)
  end do
  call trace%end
END SUBROUTINE init_rt_nbors_sc

!===============================================================================
!> Collect RT nbors. Each MHD depends on the RT task as well.
!===============================================================================
SUBROUTINE init_rt_nbors_hybrid (self, mhd)
  class (rt_setup_t):: self
  class (solver_t):: mhd
  !.............................................................................
  integer:: i_omega
  class(rt_t), pointer:: rt, omega, upstream
  class(link_t), pointer:: nbor, new
  class(task_t), pointer:: nb_task
  logical:: needed, needs_me, download
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_setup_t%init_nbors')
  rt => mhd%rt
  call rt%clear(bits%init_nbors)
  !---------------------------------------------------------------------------
  ! Loop over all MHD nbors, looking at their RT tasks
  !---------------------------------------------------------------------------
  nbor => mhd%link%nbor
  do while (associated(nbor))
    nb_task => nbor%task
    select type (nb_task)
    class is (solver_t)
      if (nb_task%rt%on) then
        new%link => nb_task%rt%link
        new%task => new%link%task
        call rt%link%add_nbor_by_rank(rt%link%nbor, new)
      end if
    end select
    nbor => nbor%next
  end do
  !-----------------------------------------------------------------------------
  ! The MHD task needs the main RT task (rt), which also needs the MHD task
  !-----------------------------------------------------------------------------
  call init_nbor_pair (mhd, needed, rt, needs_me, download)
  call trace%end
END SUBROUTINE init_rt_nbors_hybrid

!===============================================================================
!> Add mutual nbor relations to two nbor lists, making sure that the needed and
!> needs_me flags are mutually consistent.
!===============================================================================
SUBROUTINE init_nbor_pair (task1, needed, task2, needs_me, download)
  class(task_t):: task1, task2
  logical:: needed, needs_me, download
  class(link_t), pointer:: link1, link2, new
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%init_nbor_pairs')
  select type (task1)
  class is (patch_t)
  link1 => task1%link
  end select
  select type (task2)
  class is (patch_t)
  link2 => task2%link
  end select
  !-----------------------------------------------------------------------------
  ! Add the 2nd linked task as an nbor to the 1st linked task
  !-----------------------------------------------------------------------------
  allocate (new)
  new%link => link2
  new%task => new%link%task
  new%needed = needed
  new%needs_me = needs_me
  new%download = download
  call link1%add_nbor_by_rank (link1%nbor, new)
  !-----------------------------------------------------------------------------
  ! Add the 1st linked task as an nbor to the 2nd linked task
  !-----------------------------------------------------------------------------
  allocate (new)
  new%link => link1
  new%task => new%link%task
  new%needed = needs_me
  new%needs_me = needed
  new%download = .false.
  call link2%add_nbor_by_rank (link2%nbor, new)
  !-----------------------------------------------------------------------------
  call trace%end
END SUBROUTINE init_nbor_pair

END MODULE rt_setup_mod
