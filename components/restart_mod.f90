MODULE restart_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE data_io_mod
  USE amr_io_mod
  USE task_factory_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE link_mod
  USE task_list_mod
  implicit none
  private
  type, public:: restart_t
  contains
    procedure:: init
    procedure:: read_amr
    procedure:: read_legacy
  end type
  character(120):: filename, from_dir='data/'
  integer:: id, rank, unit, iostat
  integer:: snapshot=-1, nml_version=4
  type(restart_t), public:: restart
CONTAINS

!===============================================================================
!> Select I/O method
!===============================================================================
SUBROUTINE init (self, task_list)
  class(restart_t)            :: self
  type(task_list_t),  pointer :: task_list
  namelist /restart_params/ from_dir, snapshot
  !.............................................................................
  call trace%begin ('restart_t%init')
  !-----------------------------------------------------------------------------
  ! Namelist
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, restart_params, iostat=iostat)
  call io%namelist_warning ('restart_params', iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, restart_params)
  end if
  write (io%namelist, restart_params)
  !-----------------------------------------------------------------------------
  ! Only allocate a task list if actuall restarting (snapshot non-negative)
  !-----------------------------------------------------------------------------
  io%restart = snapshot
  if (snapshot >= 0) then
    !-----------------------------------------------------------------------------
    ! Select method
    !-----------------------------------------------------------------------------
    allocate (task_list)
    select case (io%method)
    case ('amr')
      call self%read_amr (task_list)
    case ('legacy')
      call self%read_legacy (task_list)
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
!> MPI parallel I/O based method
!===============================================================================
SUBROUTINE read_amr (self, task_list)
  class(restart_t)        :: self
  type(task_list_t),  pointer :: task_list
  !.............................................................................
  class(task_t),      pointer :: task
  class(link_t),      pointer :: link
  character(16)               :: type
  integer(8)                  :: offset
  namelist /amr_nml/ id, type, rank, offset, nml_version
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_amr')
  !-----------------------------------------------------------------------------
  ! Populate task lists, so far with only light-weight tasks, based on metadata
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/amr.nml")') trim(from_dir), snapshot
  open (newunit=unit, file=filename, form='formatted', status='old')
  do while (.true.)
    read (unit, amr_nml, iostat=iostat)
    if (iostat /= 0) then
      exit
    end if
    if (rank == mpi%rank) then
      call task_factory%allocate (type, task)
      task%id = id
      call task%pre_init
      task%amr_offset = offset
      io%nml_version = nml_version
      if (first_time) then
        call task%input_metadata (unit, snapshot=.true.)
        first_time = .false.
      end if
      call task%input_metadata (unit)
      call task_list%append (task)
    end if
  end do
  close (unit)
  !-----------------------------------------------------------------------------
  ! Read binary data for all tasks in the task list; this is kept separate here
  ! to allow for adding initial load balancing
  !-----------------------------------------------------------------------------
  call amr_io%open_file (from_dir, snapshot)
  link => task_list%head
  do while (associated(link))
    task => link%task
    call task%init
    call task%allocate_mesg
    call amr_io%input (task)
    link => link%next
  end do
  call amr_io%close_file
  !-----------------------------------------------------------------------------
  ! Prepare task list for execution
  !-----------------------------------------------------------------------------
  call task_list%init_all_nbors
  call task_list%queue_all
  call trace%end ()
END SUBROUTINE read_amr

!===============================================================================
!> Process the data/run/SSSSS/amr.nml with task metadata, allocate tasks, and
!> add them to the task list
!===============================================================================
SUBROUTINE read_legacy (self, task_list)
  class(restart_t)        :: self
  type(task_list_t),  pointer :: task_list
  !.............................................................................
  class(task_t),      pointer :: task
  class(link_t),      pointer :: link
  character(16)               :: type
  integer(8)                  :: offset
  namelist /restart_params/ from_dir, snapshot
  namelist /amr_nml/ id, type, rank, offset, nml_version
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_legacy')
  !-----------------------------------------------------------------------------
  ! Populate task lists, so far with only light-weight tasks, based on metadata
  !-----------------------------------------------------------------------------
  do rank=1,mpi%size
    write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
      trim(from_dir), snapshot, rank
    open (newunit=unit, file=filename, form='formatted', status='old')
    do while (.true.)
      if (rank == mpi%rank) then
        allocate (experiment_t:: task)
        select type (task)
        class is (patch_t)
          call data_io%input_nml (task)
          call task%pre_init
        end select
        call task_list%append (task)
      end if
    end do
    close (unit)
  end do
  !-----------------------------------------------------------------------------
  ! Read binary data for all tasks in the task list
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => link%task
    call task%init
    select type (task)
    class is (patch_t)
      write (filename,'(a,i5.5,"/rank_",i5.5,".dat")') &
        trim(from_dir), snapshot, rank
      open (newunit=unit, file=filename, form='unformatted', status='old', &
        access='direct', recl=product(shape(task%mem))/task%nw)
      read (unit, rec=1) task%mem
    end select
    close (unit)
    link => link%next
  end do
  call trace%end ()
END SUBROUTINE read_legacy

END MODULE restart_mod
