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
  USE task_list_mod
  USE read_bifrost_mod
  USE read_stagger_mod
  implicit none
  private
  type, public:: restart_t
  contains
    procedure:: init
    procedure:: read_amr
    procedure:: read_legacy
  end type
  character(120) :: from_dir    = 'data/'
  character(16)  :: io_format   = 'amr'
  integer        :: snapshot    = -1
  integer        :: verbose     = 0
  integer        :: nml_version = 4
  type(restart_t), public:: restart
CONTAINS

!===============================================================================
!> Select I/O method
!===============================================================================
SUBROUTINE init (self, task_list)
  class(restart_t)            :: self
  type(task_list_t),  pointer :: task_list
  integer                     :: iostat
  character(120)              :: filename
  namelist /restart_params/ from_dir, snapshot, io_format, verbose
  !.............................................................................
  call trace%begin ('restart_t%init')
  !-----------------------------------------------------------------------------
  ! Namelist
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, restart_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, restart_params)
  end if
  call io%namelist_warning ('restart_params', iostat)
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
    call task_list%init
    select case (io_format)
    case ('amr')
      call self%read_amr (task_list)
    case ('legacy')
      call self%read_legacy (task_list)
    case ('bifrost')
      call read_bifrost%init (task_list)
    case ('stagger')
      call read_stagger%init (task_list)
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
  character(120)              :: filename
  character(16)               :: type
  integer(8)                  :: offset
  integer                     :: unit, iostat, id, rank
  namelist /amr_nml/ id, type, rank, offset, nml_version
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
      task%amr_offset = offset
      io%nml_version = 4
      call task%input_metadata (unit)   ! default: gpatch_t%input_metadata
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
  class(restart_t)             :: self
  type(task_list_t),   pointer :: task_list
  !.............................................................................
  class(experiment_t), pointer :: experiment
  class(patch_t),      pointer :: patch
  class(task_t),       pointer :: task
  class(link_t),       pointer :: link
  character(120)               :: filename
  character(16)                :: type
  integer(8)                   :: offset
  integer                      :: unit, iostat, recl, id, rank, n, l(3), u(3)
  real(var_k)        , pointer :: d(:,:,:), p(:,:,:), mem(:,:,:,:)
  namelist /restart_params/ from_dir, snapshot
  namelist /amr_nml/ id, type, rank, offset, nml_version
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_legacy')
  !-----------------------------------------------------------------------------
  ! Populate task lists, so far with only light-weight tasks, based on metadata.
  ! Load-balancing is a separate issue; this step only assumes that there is
  ! enough memory to read.
  !-----------------------------------------------------------------------------
  write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
    trim(from_dir), snapshot, mpi%rank
  open (newunit=unit, file=filename, form='formatted', status='old', &
    iostat=iostat)
  if (iostat /= 0) then
    write (stderr,*) 'file does not exist: '//trim(filename)
  else
    !---------------------------------------------------------------------------
    ! Loop until end-of-file on the rank file
    !---------------------------------------------------------------------------
    n = 0
    do while (.true.)
      allocate (experiment)
      call experiment%pre_init
      patch => experiment
      call data_io%input_nml (patch, unit)
      if (data_io%eof (unit)) then
        call experiment%dealloc
        exit
      end if
      n = n+1
      task => experiment
      call task_list%append (task)
    end do
    close (unit)
    if (verbose > 0) then
      write (stdout,'(a,i4,a)') 'legacy metadata for',n,' snapshot'
      flush (stdout)
    end if
    !---------------------------------------------------------------------------
    ! Read binary data for all tasks in the task list
    !---------------------------------------------------------------------------
    link => task_list%head
    do while (associated(link))
      task => link%task
      select type (task)
      class is (patch_t)
        call task%init
        task%iout     = task%iout + 1
        task%out_next = task%time + io%out_time
        l = merge (task%mesh%lb, task%mesh%li, io%guard_zones)
        u = merge (task%mesh%ub, task%mesh%ui, io%guard_zones)
        write (filename,'(a,i5.5,"/",i5.5,".dat")') &
          trim(from_dir), snapshot, task%id

        mem => task%mem(l(1):u(1),l(2):u(2),l(3):u(3),:,1,1)
        recl=(storage_size(1.0_var_k)/8)*product(shape(mem))
        open (newunit=unit, file=filename, form='unformatted', status='old', &
          access='direct', recl=recl)
        !-----------------------------------------------------------------------
        ! Initialize the task, which allocates memory, and read binary data
        !-----------------------------------------------------------------------
        read (unit, rec=1) mem
        close (unit)
        if (verbose > 2) then
          call task%stats ('restart')
        else if (verbose > 1) then
          call task%stats (task%mem(:,:,:,task%idx%d,1,1), 'restart')
        end if
      end select
      link => link%next
    end do
  end if
  call trace%end ()
END SUBROUTINE read_legacy

END MODULE restart_mod
