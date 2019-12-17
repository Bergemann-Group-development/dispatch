!===============================================================================
!> Generic restart module, which attemps to recreate the setup from a given
!> run directory and snapshot.
!===============================================================================
MODULE restart_mod
  USE io_mod
  USE mpi_coords_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE experiment_mod
  USE link_mod
  USE task_list_mod
  implicit none
  private
  type, public:: restart_t
    integer:: verbose
    type(task_list_t), pointer:: task_list => null()
    type(mpi_coords_t):: mpi
    real(8):: size(3), origin(3)
    integer:: dims(3)
  contains
    procedure:: init
    procedure:: dealloc
    procedure:: read_metadata
    procedure:: read_cartesian_metadata
    procedure:: read_amr_metadata
    procedure:: read_patch_metadata
  end type
  integer:: verbose=1
  integer:: snapshot=-1
  character(len=32):: run
  character(len=64):: file
  character(len=64):: fromdir
  type(restart_t), public:: restart
CONTAINS

!===============================================================================
!> Initialize restart parameters
!===============================================================================
SUBROUTINE init (self, task_list)
  class(restart_t):: self
  type(task_list_t), optional, pointer:: task_list
  integer:: iostat
  namelist /restart_params/ run, snapshot, verbose
  !-----------------------------------------------------------------------------
  if (present(task_list)) then
    if (associated(task_list)) then
      write (stdout,*) 'restart_t%init: return -- task_list already allocated'
      return
    end if
  end if
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%init')
  rewind (io%input)
  read (io%input, restart_params, iostat=iostat)
  !-----------------------------------------------------------------------------
  if (iostat >= 0 .and. snapshot > 0) then
    write (io%output, restart_params)
    self%verbose = verbose
    call self%read_metadata (task_list)
    !call self%select_tasks
    !call self%read_tasks
    !call self%dealloc (task_list)
  end if
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Deallocate task_list, if populated
!===============================================================================
SUBROUTINE dealloc (self, task_list)
  class(restart_t):: self
  type(task_list_t), pointer:: task_list
  class(link_t), pointer:: link, next
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%dealloc')
  if (associated(task_list)) then
    link => task_list%head
    do while (associated(link))
      next => link%next
      task => link%task
      select type (task)
      class is (experiment_t)
        call task%dealloc
      class is (patch_t)
        call task%dealloc
      class default
        call task%dealloc
      end select
      deallocate (link%task)
      deallocate (link)
      link => next
    end do
    deallocate (task_list)
    nullify (task_list)
  end if
  call trace%end ()
END SUBROUTINE dealloc

!===============================================================================
!> Read metadata for all tasks, collecting into a task_list data type
!===============================================================================
SUBROUTINE read_metadata (self, task_list)
  class(restart_t):: self
  type(task_list_t), pointer:: task_list
  type(experiment_t), pointer:: patch
  class(link_t), pointer:: link
  integer:: iostat, rank
  logical:: do_debug, do_trace, do_output, do_flags
  namelist /io_nml/ verbose, do_debug, do_trace, do_output, do_flags
  namelist /snapshot_nml/ verbose, do_debug, do_trace, do_output, do_flags
  namelist /io_params/ verbose, do_debug, do_trace, do_output, do_flags
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_metadata')
  write (fromdir,'(a,"/",a,"/")') trim(io%datadir), trim(run)
  file = trim(fromdir)//'params.nml'
  if (self%verbose > 0) &
    write(stdout,*) 'reading '//trim(file)
  open (io_unit%data,file=file,form='formatted',status='old',iostat=iostat)
  if (iostat /= 0) then
    write (stdout,*) 'restart_t%init: return -- did not find '//trim(file)
    return
  end if
  read (io_unit%data, io_params, iostat=iostat)
  call namelist_assert ('io_params', file, iostat)
  !-----------------------------------------------------------------------------
  ! Try a cartesian setup
  !-----------------------------------------------------------------------------
  if (self%read_cartesian_metadata()) then
    write (stdout, *) ' found cartesian patch geometry'
  else if (self%read_amr_metadata()) then
    write (stdout, *) ' found AMR patch geometry'
  end if
  !-----------------------------------------------------------------------------
  !> Look for a snapshot.nml file
  !-----------------------------------------------------------------------------
  write (fromdir,'(a,i5.5,"/")') trim(fromdir), snapshot
  file = trim(fromdir)//'/snapshot.nml'
  if (self%verbose > 0) &
    write (stdout,*) 'reading '//trim(file)
  open (io_unit%data,file=file,form='formatted',status='old')
  !-----------------------------------------------------------------------------
  !> Read the io_nml and snapshot_nml namelists
  !-----------------------------------------------------------------------------
  read (io_unit%data, io_nml, iostat=iostat)
  call namelist_assert ('io_nml', file, iostat)
  read (io_unit%data, snapshot_nml, iostat=iostat)
  call namelist_assert ('snapshot_nml', file, iostat)
  !-----------------------------------------------------------------------------
  ! Read metadata for all tasks
  !-----------------------------------------------------------------------------
  allocate (task_list)
  self%task_list => task_list
  do rank=0,self%mpi%size-1
    write (file,'(a,"rank_",i5.5,"_patches.nml")') trim(fromdir), rank
    if (self%verbose > 0) &
      write (stdout,*) 'reading '//trim(file)
    open (io_unit%data,file=file,form='formatted',status='old')
    do while (.true.)
      allocate (link)
      allocate (patch)
      link%task => patch
      patch%link => link
      if (self%read_patch_metadata (patch)) then
        patch%rank = rank
        call task_list%append_link (link)
      else
        exit
      end if
    end do
    close (io_unit%data)
  end do
  close (io_unit%data)
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call self%task_list%init_all_nbors
  call self%task_list%reset_status
  call self%task_list%count_status
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  io%ntask = self%task_list%na
  io%nwrite = self%task_list%na
  io%ntotal = product(self%dims)
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE read_metadata

!===============================================================================
!> Read the patch metadata and impose on a newly created patch
!===============================================================================
FUNCTION read_patch_metadata (self, patch) RESULT (ok)
  class(restart_t):: self
  type(experiment_t), pointer:: patch
  logical:: ok
  integer:: iostat
  !.............................................................................
  integer:: id, level
  real(8):: size(3), position(3), time, dtime
  namelist /patch_nml/ id, level, size, position, time, dtime
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_patch_metadata')
  read (io_unit%data, patch_nml, iostat=iostat)
  if (self%verbose > 2) &
    write (stdout,*) iostat
  call namelist_assert ('patch_nml', file, 0)
  ok = iostat >= 0
  if (ok) then
    patch%id = id
    patch%level = level
    patch%time = time
    patch%dtime = dtime
    patch%size = size
    patch%position = position
    if (self%verbose > 2) &
      print 1, patch%id, patch%size, patch%ds, 'before'
      1 format(i6,2(2x,1p,3e12.3),4x,a)
    call patch%init
    if (self%verbose > 2) &
      print 2, patch%id, patch%size, patch%ds, patch%mesh%d, 'after'
      2 format(i6,3(2x,1p,3e12.3),4x,a)
    patch%id = id
    patch%level = level
    patch%time = time
    patch%dtime = dtime
    patch%size = size
    print 2, patch%id, patch%size, patch%ds, patch%mesh%d, 'reset'
    patch%position = position
    if (self%verbose > 1) &
      write (stdout,patch_nml)
  end if
  call trace%end ()
  !-----------------------------------------------------------------------------
END FUNCTION read_patch_metadata

!===============================================================================
!> Try reading a cartesian_params namelist, and the corresponding tasks
!===============================================================================
FUNCTION read_cartesian_metadata (self) RESULT (ok)
  class(restart_t):: self
  type(patch_t), pointer:: patch
  logical:: ok
  integer:: iostat, rank
  real(8), dimension(3):: size, origin
  integer, dimension(3):: mpi_dims, dims
  namelist /cartesian_params/ size, origin, mpi_dims, dims, rank
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_cartesian_metadata')
  read (io_unit%data, cartesian_params, iostat=iostat)
  call namelist_assert ('cartesian_params', file, iostat)
  ok = (iostat >= 0)
  if (ok) then
    self%size = size
    self%origin = origin
    self%dims = dims
    self%mpi%dims = mpi_dims
    self%mpi%size = product(mpi_dims)
    self%mpi%rank = rank
    close (io_unit%data)
  end if
  call trace%end ()
END FUNCTION read_cartesian_metadata

!===============================================================================
!> Try reading an AMR run namelist, and the corresponding tasks
!===============================================================================
FUNCTION read_amr_metadata (self) RESULT (ok)
  class(restart_t):: self
  type(patch_t), pointer:: patch
  logical:: ok
  integer:: iostat, rank
  real(8), dimension(3):: size
  integer, dimension(3):: mpi_dims, dims
  namelist /cartesian_params/ size, mpi_dims, dims, rank
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_cartesian_metadata')
  read (io_unit%data, cartesian_params, iostat=iostat)
  call namelist_assert ('amr_params', file, iostat)
  ok = iostat==0
  if (ok) then
    close (io_unit%data)
  end if
  call trace%end ()
END FUNCTION read_amr_metadata

!===============================================================================
!> Issue messages about namelists -- missing or with errors
!===============================================================================
SUBROUTINE namelist_assert (label, file, iostat)
  character(len=*):: label
  character(len=64):: file
  integer:: iostat
  if (iostat < 0) then
    call io%abort (trim(label)//' missing in '//trim(file))
  else if (iostat > 0) then
    write (stdout,*) ' some read error in '//trim(label)//' in '//trim(file)
  end if
END SUBROUTINE namelist_assert 

END MODULE restart_mod
