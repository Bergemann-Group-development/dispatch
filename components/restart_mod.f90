MODULE restart_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE data_io_mod
  USE kinds_mod
  USE amr_io_mod
  USE task_factory_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE link_mod
  USE list_mod
  USE task_list_mod
  USE read_bifrost_mod
  USE read_stagger_mod
  USE load_balance_mod
  USE bits_mod
  implicit none
  private
  type, public:: restart_t
  contains
    procedure:: init
    procedure:: read
  end type
  character(120) :: data        = 'data/'
  character(120) :: run         = ''
  character(16)  :: method      = 'amr'
  integer        :: snapshot    = -1
  integer        :: verbose     = 0
  integer        :: nml_version = 4
  type(restart_t), public:: restart
CONTAINS

!===============================================================================
!> Select I/O method
!===============================================================================
SUBROUTINE init (self, task_list)
  class(restart_t)             :: self
  class(task_list_t),  pointer :: task_list
  integer                      :: iostat
  character(120)               :: filename
  namelist /restart_params/ snapshot, run, data, method, verbose
  !.............................................................................
  call trace%begin ('restart_t%init')
  !-----------------------------------------------------------------------------
  ! Namelist
  !-----------------------------------------------------------------------------
  run = io%run_name
  data = io%data_dir
  rewind (io%input)
  read (io%input, restart_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, restart_params)
  end if
  call io%namelist_warning ('restart_params', iostat)
  call io%ensure_dirname (data)
  write (io%namelist, restart_params)
  !-----------------------------------------------------------------------------
  ! Only allocate a task list if actuall restarting (snapshot non-negative)
  !-----------------------------------------------------------------------------
  io%restart = snapshot
  io%inputdir = trim(data)//trim(run)
  call io%ensure_dirname (io%inputdir)
  if (snapshot >= 0) then
    !-----------------------------------------------------------------------------
    ! Select method
    !-----------------------------------------------------------------------------
    allocate (task_list)
    call task_list%init
    select case (method)
    case ('bifrost')
      call read_bifrost%init (task_list)
    case ('stagger')
      call read_stagger%init (task_list)
    case default
      call self%read (task_list)
    end select
    !-----------------------------------------------------------------------------
    ! Prepare task list for execution
    !-----------------------------------------------------------------------------
    call task_list%init_all_nbors
    call task_list%queue_all
  end if
  !.............................................................................
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Generic method
!===============================================================================
SUBROUTINE read (self, task_list)
  class(restart_t)             :: self
  class(task_list_t),  pointer :: task_list
  !.............................................................................
  class(experiment_t), pointer :: experiment
  class(patch_t),      pointer :: patch
  class(task_t),       pointer :: task
  class(link_t),       pointer :: link, next
  class(list_t),       pointer :: list
  class(task_list_t),  pointer :: tlist
  integer                      :: n, unit
  logical                      :: ok
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read')
  !-----------------------------------------------------------------------------
  ! Populate task list, so far with only light-weight tasks, based on metadata
  !-----------------------------------------------------------------------------
  n = 0
  unit = 0
  do while (.true.)
    allocate (experiment)
    call experiment%pre_init ()
    patch => experiment
    call data_io%input_nml (patch, unit)
    io%start_time = max(io%start_time, task%time)
    if (data_io%eof (unit)) then
      call experiment%dealloc
      exit
    end if
    task => experiment
    call task_list%append (task)
    n = n+1
  end do
  close (unit)
  !-----------------------------------------------------------------------------
  ! Decide which tasks to keep, and append them to a new task list
  !-----------------------------------------------------------------------------
  call task_list%init_all_nbors
  call task_list%reset_status
  list => task_list
  call load_balance%initial (list)
  !-----------------------------------------------------------------------------
  ! If the task rank matches, append the task to a new task list
  !-----------------------------------------------------------------------------
  allocate (tlist)
  call tlist%init
  link => task_list%head
  do while (associated(link))
    task => link%task
    if (task%rank == mpi%rank .or. task%is_set(bits%virtual)) then
      call tlist%append (task)
    else
      select type (task)
      class is (experiment_t)
        call task%dealloc ()
      end select
      deallocate (task)
    end if
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Deallocate the initial task list and its links, and point it to the new list
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    next => link%next
    deallocate (link)
    link => next
  end do
  deallocate (task_list)
  task_list => tlist
  !-----------------------------------------------------------------------------
  ! Prepare nbors, status bits, status counts, I/O coounts
  !-----------------------------------------------------------------------------
  call task_list%init_all_nbors
  call task_list%reset_status
  call task_list%count_status
  call task_list%set_io_counts
  !-----------------------------------------------------------------------------
  ! Allocate and read in remaining tasks
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (experiment_t)
      call task%init ()
      task%t = task%time
      task%dt = 0d0
      call task%input (ok)
    end select
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  ! Prepare task list for execution
  !-----------------------------------------------------------------------------
  call task_list%queue_all
  call trace%end ()
END SUBROUTINE read

END MODULE restart_mod
