!===============================================================================
!> Generic restart module, which attemps to recreate the setup from a given
!> run directory and snapshot.
!===============================================================================
MODULE restart_mod
  USE io_mod
  USE mpi_mod
  USE mpi_coords_mod
  USE data_io_mod
  USE trace_mod
  USE task_mod
  USE bits_mod
  USE patch_mod
  USE experiment_mod
  USE link_mod
  USE task_list_mod
  implicit none
  private
  type, public:: restart_t
    integer:: verbose
    integer:: istep
    type(task_list_t), pointer:: task_list => null()
    type(mpi_coords_t):: mpi
    real(8):: size(3), origin(3), velocity(3)
    integer:: dims(3)
    logical:: periodic(3)
  contains
    procedure:: init
    procedure:: read_metadata
    procedure:: read_cartesian_metadata
    procedure:: read_amr_metadata
    procedure:: read_patch_metadata
    procedure:: read_task_data
    procedure:: select_tasks
    procedure:: dealloc
  end type
  !-- locals -------------------------------------------------------------------
  integer:: verbose=1
  integer:: snapshot=-1
  integer:: nml_version
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
  if (present (task_list)) then
    if (associated(task_list)) then
      write (stdout,*) 'restart_t%init: return -- task_list already allocated'
      return
    end if
  end if
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%init')
  call io%header ('restart_t%init.')
  rewind (io%input)
  read (io%input, restart_params, iostat=iostat)
  !-----------------------------------------------------------------------------
  if (iostat >= 0 .and. snapshot >= 0) then
    allocate (self%task_list)
    task_list => self%task_list
    call self%task_list%init
    write (io%output, restart_params)
    self%verbose = verbose
    call self%read_metadata                     ! get size & position
    call self%task_list%initialize              ! build nbor lists
    call self%select_tasks                      ! select tasks by MPI rank
    !---------------------------------------------------------------------------
    ! Initialize the I/O counters
    !---------------------------------------------------------------------------
    io%ntask  = self%task_list%na
    io%nwrite = self%task_list%na
    io%ntotal = product(self%dims)*product(self%mpi%dims)
    if (self%verbose > 0) &
      write (stdout,*) 'restart_t%init: ntask, nwrite, ntotal =', &
      io%ntask, io%nwrite, io%ntotal
    call self%read_task_data                    ! read task data
  end if
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Read metadata for all tasks, collecting into a task_list data type
!===============================================================================
SUBROUTINE read_metadata (self)
  class(restart_t):: self
  type(experiment_t), pointer:: patch
  class(link_t), pointer:: link
  integer:: iostat, rank, n_patch
  !-- io_nml -------------------------------------------------------------------
  integer:: format, time_derivs, verbose, nml_version
  logical:: do_debug, do_trace, do_output, do_flags
  character(len=16):: method
  namelist /io_nml/ format, ntotal, out_time, guard_zones, time_derivs, method, &
    nml_version
  !-- snapshot_nml -------------------------------------------------------------
  integer:: ioformat, restart, iout, ntotal, istep, mesh_type, level, ncell(3), &
    li(3), ui(3), n(3), ng(3), gn(3), nv, nw, nt, omp_nthreads, mpi_size
  real(8):: position(3), size(3), ds(3), box(3), velocity(3), gamma, time, &
    out_time
  real(4):: quality
  character(len=16):: kind, eos, opacity
  logical:: periodic(3), guard_zones, no_mans_land
  namelist /snapshot_nml/ ioformat, iout, time, ntotal, istep, mesh_type, &
    position, size, ds, box, velocity, level, quality, gamma, ncell, li, ui, n, &
    ng, gn, nv, nt, nw, kind, eos, opacity, periodic, guard_zones, time_derivs, &
    omp_nthreads, mpi_size, no_mans_land, top
  !-- io_params ----------------------------------------------------------------
  logical:: do_validate, do_legacy, do_direct, omp_trace, namelist_errors
  integer:: id_debug, levelmax, log_sent, task_logging
  real(8):: log_interval
  character(len=64):: top, datadir, inputdir
  namelist /io_params/ verbose, do_debug, do_trace, do_output, do_flags, &
    do_validate, do_legacy, do_direct, levelmax, omp_trace, id_debug, &
    method, restart, format, guard_zones, time_derivs, log_sent, task_logging, &
    namelist_errors, nml_version, log_interval, top, datadir, inputdir
  !-----------------------------------------------------------------------------
  ! Read namelists from params.nml and nnnnnn/snapshot.nml
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_metadata')
  write (fromdir,'(a,"/",a,"/")') trim(io%datadir), trim(run)
  io%inputdir = trim(fromdir)
  file = trim(fromdir)//'restart.nml'
  if (self%verbose > 0) &
    write(stdout,*) 'restart_t: reading '//trim(file)
  open (io_unit%data,file=trim(file),form='formatted',status='old',iostat=iostat)
  if (iostat /= 0) then
    write (stdout,*) 'restart_t%init: return -- did not find '//trim(file)
    call io%abort ('missing '//trim(file))
  end if
  !-----------------------------------------------------------------------------
  ! io_params
  !-----------------------------------------------------------------------------
  read (io_unit%data, io_params, iostat=iostat)
  !...................=========................
  if (self%verbose > 2) &
    write (stdout, io_params)
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
    write (stdout,*) 'restart_t: reading '//trim(file)
  open (io_unit%data,file=file,form='formatted',status='old')
  !-----------------------------------------------------------------------------
  !> Read the io_nml and snapshot_nml namelists
  !-----------------------------------------------------------------------------
  read (io_unit%data, io_nml, iostat=iostat)
  !...................======................
  if (self%verbose > 2) &
    write (stdout, io_nml)
  call namelist_assert ('io_nml', file, iostat)
  io%method = method
  read (io_unit%data, snapshot_nml, iostat=iostat)
  !...................============................
  if (self%verbose > 2) &
    write (stdout, snapshot_nml)
  call namelist_assert ('snapshot_nml', file, iostat)
  !-----------------------------------------------------------------------------
  ! Parameters used while reading task data, later reset for writing
  !-----------------------------------------------------------------------------
  self%periodic  = periodic
  self%velocity  = velocity
  self%istep     = istep
  io%format      = ioformat
  io%ntotal      = ntotal
  io%time_derivs = time_derivs
  io%guard_zones = guard_zones
  !-----------------------------------------------------------------------------
  ! Read metadata for all tasks
  !-----------------------------------------------------------------------------
  do rank=0,self%mpi%size-1
    write (file,'(a,"rank_",i5.5,"_patches.nml")') trim(fromdir), rank
    if (self%verbose > 0) &
      write (stdout,*) 'restart_t: reading metadata from '//trim(file)
    open (io_unit%data,file=file,form='formatted',status='old')
    n_patch = 0
    do while (.true.)
      allocate (link)
      allocate (patch)
      link%task => patch
      patch%link => link
      if (self%read_patch_metadata (patch)) then
        patch%rank = rank
        call self%task_list%append_link (link)
        n_patch = n_patch+1
      else
        exit
      end if
    end do
    if (self%verbose > 0) &
      write (stdout,'(1x,a,i6,1x,a)') 'metadata for', n_patch, 'patches read'
    close (io_unit%data)
  end do
  close (io_unit%data)
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
  integer:: id, level, istep
  real(8):: size(3), position(3), time, dtime
  real(8), save:: velocity(3)=0d0
  namelist /patch_nml/ id, level, size, position, time, dtime, velocity, istep
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_patch_metadata')
  !-----------------------------------------------------------------------------
  ! Read, after setting defaults values that may not exist in patch_nml
  !-----------------------------------------------------------------------------
  velocity = self%velocity
  istep = self%istep
  read (io_unit%data, patch_nml, iostat=iostat)
  !...................=========................
  if (self%verbose > 3) &
    write (stdout, patch_nml)
  call namelist_assert ('patch_nml', file, 0)
  !-----------------------------------------------------------------------------
  ! Pick up values and use them to initialize the patch task
  !-----------------------------------------------------------------------------
  ok = iostat >= 0
  if (ok) then
    patch%id       = id
    patch%time     = time
    patch%size     = size
    patch%position = position
    patch%origin   = position - size/2d0
    !...........................................................................
    if (self%verbose > 2) &
      print 1, patch%id, patch%position, patch%size, patch%ds, 'before'
      1 format(i6,3(2x,1p,3e10.2),3x,a)
    !...........................................................................
    call patch%init
    !...........................................................................
    if (self%verbose > 2) &
      print 2, patch%id, patch%position, patch%mesh%p, patch%size, &
        patch%mesh%s, patch%ds, patch%mesh%d, 'after'
      2 format(i6,6(2x,1p,3e10.2),3x,a)
    !...........................................................................
    patch%id       = id
    patch%level    = level
    patch%time     = time
    patch%t        = time
    patch%dtime    = dtime
    patch%dt       = dtime
    patch%istep    = istep
    patch%size     = size
    patch%position = position
    patch%origin   = position - size/2d0
    patch%velocity = velocity
    patch%iit      = 1
    patch%it       = 1
    patch%new      = 2
    patch%iit(patch%nt) = 2
    patch%out_next = (nint(patch%time/io%out_time+1.0))*io%out_time
print *, 'restart_t%read_patch_metadata: id, out_next =', &
patch%id, patch%out_next, patch%unsigned
    !---------------------------------------------------------------------------
    ! Derived
    !---------------------------------------------------------------------------
    patch%guard_zones = io%guard_zones
    patch%time_derivs = io%time_derivs
    patch%llc_cart    = patch%origin
    patch%centre_nat  = patch%position
    patch%llc_nat     = patch%llc_cart
    patch%periodic    = self%periodic
    call patch%set(bits%static)
    !...........................................................................
    if (self%verbose > 1) &
      print 2, patch%id, patch%position, patch%mesh%p, patch%size, &
        patch%mesh%s, patch%ds, patch%mesh%d, 'reset'
  end if
  if (self%verbose > 1) then
  end if
  call trace%end ()
END FUNCTION read_patch_metadata

!===============================================================================
!> Try reading a cartesian_params namelist, and the corresponding tasks, saving
!> useful information in self
!===============================================================================
FUNCTION read_cartesian_metadata (self) RESULT (ok)
  class(restart_t):: self
  type(patch_t), pointer:: patch
  logical:: ok
  integer:: iostat
  real(8), dimension(3):: size, origin, position
  integer, dimension(3):: mpi_dims, dims, per_rank(3)
  logical:: omp_init, periodic(3), face_nbors(3)
  namelist /cartesian_params/ size, position, dims, mpi_dims, per_rank, origin, &
    face_nbors, omp_init, periodic, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_cartesian_metadata')
  read (io_unit%data, cartesian_params, iostat=iostat)
  !...................================................
  if (self%verbose > 2) &
    write (stdout, cartesian_params)
  call namelist_assert ('cartesian_params', file, iostat)
  ok = (iostat >= 0)
  if (ok) then
    self%size     = size                        ! size of the box
    self%origin   = origin                      ! origin of the box
    self%dims     = dims                        ! patch arrangement
    self%mpi%dims = mpi_dims                    ! MPI split
    self%mpi%size = product(mpi_dims)           ! number of MPI ranks
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
  real(8), dimension(3):: size, origin
  integer, dimension(3):: mpi_dims, dims, per_rank(3)
  logical:: omp_init, periodic(3), face_nbors(3)
  namelist /cartesian_params/ size, dims, mpi_dims, per_rank, origin, &
    face_nbors, omp_init, periodic
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_cartesian_metadata')
  read (io_unit%data, cartesian_params, iostat=iostat)
  if (self%verbose > 2) &
    write (stdout, cartesian_params)
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
SUBROUTINE select_tasks (self)
  class(restart_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  logical:: keep
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%select_tasks')
  !-----------------------------------------------------------------------------
  ! Keep real and virtual tasks on this rank
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  do while (associated(link))
    keep = .false.
    if (link%task%rank == mpi%rank) then
      keep = .true.
      exit
    end if
    nbor => link%nbor
    do while (associated (nbor))
      if (nbor%task%rank == mpi%rank) then
        keep = .true.
        exit
      end if
      nbor => nbor%next
    end do
    !---------------------------------------------------------------------------
    ! Deallocate and remove other tasks, with their nbor_lists
    !---------------------------------------------------------------------------
    if (.not. keep) then
      task => link%task
      if (self%verbose > 1) &
        write (stdout,*) 'restart_t: remove task', task%id
      select type (task)
      class is (experiment_t)
      call task%dealloc
      end select
      call link%remove_nbor_list2 (link%nbor)
      call self%task_list%remove (link)
    end if
    link => link%next
  end do
  call trace%end ()
END SUBROUTINE select_tasks

!===============================================================================
!> Issue messages about namelists -- missing or with errors
!===============================================================================
SUBROUTINE read_task_data (self)
  class(restart_t):: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  logical:: ok
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%read_task_data')
  if (self%verbose > 0) &
    write (stdout,'(1x,a,i5.5,a,i6,a)') &
      'restart_t: reading task data from '//trim(io%inputdir), snapshot, &
      '/ for', self%task_list%n,' tasks'
  io%restart = snapshot
  io%nml_version = nml_version
  i = 0
  link => self%task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (experiment_t)
      if (self%verbose > 2) &
        write (stdout,*) 'reading patch data for task', task%id
      !-------------------------------------------------------------------------
      ! The experiment input procedure should always work, for any task, but
      ! may rely on task%iout or task%restart being set correctly
      !-------------------------------------------------------------------------
      task%restart = snapshot
      task%iout    = snapshot
      call task%input (ok)
      call task%gpatch_t%pre_update ()
      i = i+1
    class default
      if (self%verbose > 1) &
        write (stdout,*) 'no data for task', task%id
    end select
    link => link%next
  end do
  io%restart = -9
  call trace%end ()
END SUBROUTINE read_task_data

!===============================================================================
!> Issue messages about namelists -- missing or with errors
!===============================================================================
SUBROUTINE namelist_assert (label, file, iostat)
  character(len=*):: label
  character(len=64):: file
  integer:: iostat
  !-----------------------------------------------------------------------------
  if (iostat < 0) then
    call io%abort (trim(label)//' missing in '//trim(file))
  else if (iostat > 0) then
    write (stdout,*) ' some read error in '//trim(label)//' in '//trim(file)
  end if
END SUBROUTINE namelist_assert 

!===============================================================================
!> Deallocate task_list, if populated
!===============================================================================
SUBROUTINE dealloc (self)
  class(restart_t):: self
  class(link_t), pointer:: link, next
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('restart_t%dealloc')
  if (associated(self%task_list)) then
    link => self%task_list%head
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
    deallocate (self%task_list)
    nullify (self%task_list)
  end if
  call trace%end ()
END SUBROUTINE dealloc

END MODULE restart_mod
