!===============================================================================
!> Module handling setup of radiative transfer tasks, adding them to the
!> already existing task_list.
!===============================================================================
MODULE rt_setup_mod
  USE nbors_mod
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
          call nbors%init (task)
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

END MODULE rt_setup_mod
