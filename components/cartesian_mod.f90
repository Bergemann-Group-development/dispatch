MODULE cartesian_mod
  USE task_list_mod
  USE experiment_mod
  USE gpatch_mod
  USE patch_mod
  USE link_mod
  USE list_mod
  USE task_mod
  USE trace_mod
  USE io_mod
  USE timer_mod
  USE mpi_coords_mod
  USE mpi_mod
  USE omp_mod
  USE omp_timer_mod
  USE bits_mod
  implicit none
  private
  type, public, extends(task_t):: rank_t
  end type
  type, public:: cartesian_t
    real(8):: position(3)
    real(8):: size(3)
    integer:: dims(3)
    type(task_list_t), pointer:: task_list => null()
  contains
    procedure init
    procedure, nopass:: diagnostics
    procedure:: initial_refine
  end type
  integer:: verbose=0
  logical:: omp_init=.true.
  type(rank_t):: rank
  integer, save:: np=0                  ! patches per process
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!===============================================================================
!> Distribute patches in a Cartesian arrangement in a box
!===============================================================================
SUBROUTINE init (self, task_list, label)
  class(cartesian_t):: self
  type(task_list_t), optional, pointer:: task_list
  character(len=*), optional:: label
  !.............................................................................
  class(list_t), pointer:: tlist
  class(link_t), pointer:: link, next
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i, j, k, id, ip, iostat
  integer:: dims(3) = [4,4,4]
  integer:: patches_per_mpi(3), refine
  integer:: mpi_dims(3) = [1,1,1]
  integer:: per_rank(3) = [0,0,0]
  real(8):: size(3) = [1.0,1.0,1.0]
  real(8):: origin(3) = [-0.5,-0.5,-0.5]
  real(8):: position(3) = tiny(1d0)
  real(8):: wc
  logical:: face_nbors=.false.
  logical:: periodic(3)=.true.
  integer, save:: itimer=0
  namelist /cartesian_params/ size, position, dims, mpi_dims, per_rank, origin, &
    face_nbors, omp_init, periodic, verbose
  character(len=120):: ids = &
    '$Id$'
  !-----------------------------------------------------------------------------
  if (present (task_list)) then
    if (associated(task_list)) then
      write (stdout,*) 'cartesian_t%init: return -- task_list already allocated'
      return
    else
      allocate(task_list)
      self%task_list => task_list
    end if
  else
    allocate(self%task_list)
  end if
  !-----------------------------------------------------------------------------
  call trace%begin ('cartesian_t%init')
  call io%header('cartesian_t%init: Cartesian patch arrangement:')
  call trace%print_id (ids)
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, cartesian_params, iostat=iostat)
  call io%namelist_warning ('cartesian_params', iostat)
  if (any(position /= tiny(1d0))) then
    origin = position - size/2d0
    call io%warning ('cartesian_t%init: position used, origin ignored')
  else
    position = origin + size/2d0
  end if
  if (any(per_rank/=0)) then
    where (per_rank/=0) dims=per_rank*mpi_dims
  end if
  per_rank = dims/mpi_dims
  if (any(dims == 2 .and. periodic)) then
    call io%abort ('two patches per dimension is not allowed')
  end if
  self%position = position
  self%dims = dims
  self%size = size
  io%dims     = dims
  if (mpi%size==1) &
    mpi_dims = 1
  !-----------------------------------------------------------------------------
  ! Generate a (periodic) MPI geometry, with dims specified in the input file,
  ! or assigned automatically by mpi_coords%init.
  !-----------------------------------------------------------------------------
  if (mpi%size == 1) then
     mpi_coords%dims = 1
  else
    if (product(mpi_dims)==1) then
      call mpi_coords%init (dims=dims)
    else
      call mpi_coords%init (mpi_dims=mpi_dims)
    end if
  end if
  mpi_dims = mpi_coords%dims
  patches_per_mpi = self%dims/mpi_dims
  mpi_coords%npatch = patches_per_mpi
  io%mpi_dims = mpi_dims
  write (io%namelist, cartesian_params)
  !-----------------------------------------------------------------------------
  ! Set attritbutes of the rank_t data type -- just enough to use distance().
  !-----------------------------------------------------------------------------
  rank%periodic = .true.
  rank%box = self%size
  rank%size = rank%box/mpi_dims
  rank%ds = rank%size/patches_per_mpi
  rank%position = (mpi_coords%rank_to_coords(mpi%rank) + 0.5_8)*rank%size + origin
  !-----------------------------------------------------------------------------
  ! Begin on task list
  !-----------------------------------------------------------------------------
  call self%task_list%init ('tlist')
  tlist => self%task_list
  self%task_list%size = self%size
  self%task_list%dims = self%dims
  self%task_list%position = self%position
  io%dims = self%dims
  self%task_list%n_tasks = product(self%dims)
  self%task_list%face_nbors = face_nbors
  !-----------------------------------------------------------------------------
  ! Allocate an array of links, pointing to patches, and give the patches each
  ! id, rank, size, position, integer position, and status bits.  Discard
  ! external patches and add the rest to the task_list.
  !-----------------------------------------------------------------------------
  ip = 0
  do k=0,self%dims(3)-1
  do j=0,self%dims(2)-1
  do i=0,self%dims(1)-1
    allocate (link, patch); link%task => patch; patch%link => link
    task => link%task
    ! --- patch 3D and 1D coordinates ---
    patch%ipos = [i,j,k]
    patch%rank = mpi_coords%coords_to_rank ([i,j,k]/patches_per_mpi)
    if (patch%rank == mpi%rank) then
      ip = ip + 1
      patch%ip = ip
    end if
    !---------------------------------------------------------------------------
    ! Generate the same IDs on all ranks, pruning the external tasks below
    !---------------------------------------------------------------------------
    call patch%task_t%init_unique (same=.true.)
    patch%task_list => self%task_list
    patch%box = self%size
    patch%size = self%size / self%dims
    patch%position = ([i,j,k]+0.5d0)*patch%size + origin
    !if (any(abs(rank%distance (patch)) > 0.5*(rank%size + 1.5*rank%ds))) then
    !  deallocate (link, patch)
    !  cycle
    !end if
    call patch%set (bits%virtual)
    patch%ds = patch%size*0.001
    patch%origin = origin
    patch%llc_cart = patch%position - 0.5 * patch%size
    patch%centre_nat = patch%position
    patch%llc_nat = patch%llc_cart
    patch%periodic = periodic
    call patch%set(bits%static)
    if (all(self%dims==1)) call task%set(bits%root_grid)
    ! -- initially append all tasks to the task list, to build nbor lists ------
    call self%initial_refine(tlist, patch,patch%level,ip)
    !call self%task_list%append_link (link)
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Generate nbor lists and set status bits for all tasks
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(stdout,*) 'cartesian_t%init: init_all_nbors'
  call self%task_list%init_all_nbors
  if (verbose > 0) write(stdout,*) 'cartesian_t%init: reset_status'
  call self%task_list%reset_status
  !-----------------------------------------------------------------------------
  ! Remove tasks that are external, and make sure to also remove them from the
  ! nbor lists of virtual tasks
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  do while (associated(link))
    next => link%next
    if (link%task%is_set (bits%external)) then
      call link%remove_nbor_list2 (link%nbor)
      call self%task_list%remove (link)
      call link%task%dealloc()
      deallocate (link%task)
      deallocate (link)
    end if
    link => next
  end do
  !-----------------------------------------------------------------------------
  ! This counts are not fully correct until the patch meshes have been allocated
  !-----------------------------------------------------------------------------
  if (verbose > 0) write(stdout,*) 'cartesian_t%init: count status'
  call self%task_list%count_status
  !-----------------------------------------------------------------------------
  ! The I/O conunters must be set before initializing patches
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  write (stdout,'(a,i8)') ' number of tasks in task list:', self%task_list%n
  write (stdlog,*) 'cartesian_mod: ntask,nwrite,ntotal=', &
    io%ntask, io%nwrite, io%ntotal
  !-----------------------------------------------------------------------------
  ! Initialize tasks
  !-----------------------------------------------------------------------------
  !$omp parallel default(shared)
  call init_exp (self, origin)
  !$omp end parallel
  !-----------------------------------------------------------------------------
  ! Refresh the counts, now with bits%has_output set on tasks
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  write (stdlog,*) 'cartesian_mod: ntask,nwrite,ntotal=', &
    io%ntask, io%nwrite, io%ntotal
  !-----------------------------------------------------------------------------
  ! This call is not valid until the patch mem and meshes have been allocated
  !-----------------------------------------------------------------------------
  if (verbose > 0) &
    write(stdout,*) 'cartesian_t%init: init boundaries'
  call self%task_list%init_bdries
  !-----------------------------------------------------------------------------
  if (io_unit%master) then
    write (io%namelist, cartesian_params)
    flush (io%namelist)
  end if
  if (verbose > 0) &
    write(stdout,*) 'cartesian_t%init: diagnostics'
  call self%diagnostics()
  if (verbose > 0) &
    write(stdout,*) 'cartesian_t%init: print tasks'
  call self%task_list%print_tasks
  if (verbose > 0) &
    write(stdout,*) 'cartesian_t%init: done'
  flush (stdout)
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

recursive SUBROUTINE initial_refine(self,tlist, task,level_min, ip)
  class(cartesian_t):: self
  class(list_t), pointer:: tlist
  class(experiment_t), pointer:: task
  integer:: level_min, ip
  ! ----------------------------------------------------------------------------
  class(link_t), pointer:: link
  class(experiment_t), pointer:: patch
  logical, save:: on = .false.
  real(8), save:: center(3)
  real(8), save:: extent(3), ratio=2.0
  integer, save:: n_refine(3)
  logical, save:: first_time=.true.
  logical, save:: keep_parent=.true.
  real(8), dimension(3):: size, pos
  integer:: refine=-1, iostat
  integer:: level_max, i, j, k
  namelist /initial_refine_params/ on, n_refine, center, extent, ratio, keep_parent
  ! ............................................................................
  if (first_time) then
    first_time = .false.
    rewind (io%input); read (io%input, initial_refine_params, iostat=iostat)
    write (io%output, initial_refine_params)
    if (on) call io%header('begin initial mesh refinement')
  end if
  if (.not.on) then
    !call task%pre_init()
    call self%task_list%append_link (task%link)
    !call task%add_to_task_list (tlist)
    write (io_unit%mpi,'(3i4,3f8.3,i4,a)') i,j,k,task%position,task%rank,'i,j,k,pos,rank'
    return
  end if
  level_max = level_min+maxval(n_refine)
  ! check X
  refine=-1
  refine = max(refine,nint(level_min+(level_max-level_min)* &
           exp(-abs(task%position(1)-center(1))/(extent(1)+1e-20))) - task%level)
  ! check Y
  refine = max(refine,nint(level_min+(level_max-level_min)* &
           exp(-abs(task%position(2)-center(2))/(extent(2)+1e-20))) - task%level)
  ! check Z
  refine = max(refine,nint(level_min+(level_max-level_min)* &
           exp(-abs(task%position(3)-center(3))/(extent(3)+1e-20))) - task%level)
  if (refine>0 .and. task%level<level_max) then
    do k=1,ratio
      do j=1,ratio
        do i=1,ratio
          !---------------------------------------------------------------------
          ! Size and position
          !---------------------------------------------------------------------
          size = task%size/ratio
          pos = task%position + ([i,j,k]-(ratio+1d0)*0.5d0)*size
          !---------------------------------------------------------------------
          ! Create a refined patch here
          !---------------------------------------------------------------------
          allocate (link, patch); link%task => patch; patch%link => link
          ! --- patch 3D and 1D coordinates ---
          call patch%task_t%init_unique (same=.true.)
          patch%ipos = task%ipos
          patch%rank = task%rank
          patch%task_list => self%task_list
          if (patch%rank == mpi%rank) then
            ip = ip + 1
            patch%ip = ip
          end if
          !---------------------------------------------------------------------
          ! Generate the same IDs on all ranks, pruning the external tasks below
          !---------------------------------------------------------------------
          patch%box = self%size
          patch%size = size
          patch%position = pos
          !if (any(abs(rank%distance (patch)) > 0.5*(rank%size + 1.5*rank%ds))) then
          !  deallocate (link, patch)
          !  cycle
          !end if
          call patch%set (bits%virtual)
          patch%ds = patch%size*0.001
          patch%origin = task%origin
          patch%llc_cart = patch%position - 0.5 * patch%size
          patch%centre_nat = patch%position
          patch%llc_nat = patch%llc_cart
          patch%periodic = task%periodic
          patch%level = task%level + 1
        !  if (any(abs(rank%distance (patch)) > 0.5*(rank%size + 1.5*rank%ds))) then
        !    !call patch%dealloc()
        !    deallocate (link, patch)
        !    cycle
        !  end if 
          call patch%set(bits%static)
          if (keep_parent) then
            patch%link%parent => task%link
            patch%parentid    = task%id
          end if
          !---------------------------------------------------------------------
          call self%initial_refine(tlist, patch,level_min,ip)
        end do
      end do
    end do
    if (keep_parent) then
      !call task%pre_init()
      call self%task_list%append_link (task%link)
      !call task%add_to_task_list (tlist)
      write (io_unit%mpi,'(3i4,3f8.3,i4,a)') i,j,k,task%position,task%rank,'i,j,k,pos,rank'
    else
      link => task%link
      call link%task%dealloc()
      deallocate(task,link)
    end if
  else
    call task%pre_init()
    call self%task_list%append_link (task%link)
    !call task%add_to_task_list (tlist)
    write (io_unit%mpi,'(3i4,3f8.3,i4,a)') i,j,k,task%position,task%rank,'i,j,k,pos,rank'
  end if
END SUBROUTINE initial_refine


SUBROUTINE init_exp (self, origin)
  class(cartesian_t):: self
  real(8):: origin(3)
  !.............................................................................
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  class(experiment_t), pointer:: patch
  logical:: periodic(3)
  integer:: id
  !-----------------------------------------------------------------------------
  id = 0
  link => self%task_list%head
  do while (associated(link))
    id = id+1
    if (mod(id,omp%nthreads)==omp%thread) then
      task => link%task
      select type (task)
      class is (experiment_t)
        patch => task
        periodic = patch%periodic
        call patch%init
        if (any(patch%periodic .neqv. periodic)) then
          call io%abort ('cartesian_t%init: experiment init tried to change %periodic')
        end if
        if (omp%nthreads>1) then
          task%mem_thread = omp%thread
        else
          task%mem_thread = -1
        end if
        patch%mesh%origin = origin
        !$omp atomic
        np = np + 1
        ! Set physical boundary bits if required.
        if (.not.patch%periodic(1).and.patch%ipos(1)==0)              call patch%boundaries%set(bits%xl)
        if (.not.patch%periodic(2).and.patch%ipos(2)==0)              call patch%boundaries%set(bits%yl)
        if (.not.patch%periodic(3).and.patch%ipos(3)==0)              call patch%boundaries%set(bits%zl)
        if (.not.patch%periodic(1).and.patch%ipos(1)==self%dims(1)-1) call patch%boundaries%set(bits%xu)
        if (.not.patch%periodic(2).and.patch%ipos(2)==self%dims(2)-1) call patch%boundaries%set(bits%yu)
        if (.not.patch%periodic(3).and.patch%ipos(3)==self%dims(3)-1) call patch%boundaries%set(bits%zu)
        if (verbose>0) then
          write (io_unit%log,1) &
            'cartesian_t%init: id, rank, size, pos =', &
            patch%id, patch%rank, patch%size, patch%position, &
            patch%is_set(bits%boundary), patch%is_set(bits%virtual)
          1 format(a,i8,i6,1p,2(2x,3e14.6),2x,2l1,2x,a)
          flush (io_unit%log)
        end if
        if (patch%id == io%id_debug) then
          write (io_unit%output,1) &
            'cartesian_t%init: id, rank, size, pos =', &
            patch%id, patch%rank, patch%size, patch%position, &
            patch%is_set(bits%boundary), patch%is_set(bits%virtual), 'DBG'
          flush (io_unit%log)
        else if (io_unit%verbose>1 .and. io_unit%master) then
          write (io_unit%output,1) &
            'cartesian_t%init: id, rank, size, pos =', &
            patch%id, patch%rank, patch%size, patch%position, &
            patch%is_set(bits%boundary), patch%is_set(bits%virtual)
          flush (io_unit%log)
        end if
      end select
    end if
    link => link%next
  end do
END SUBROUTINE init_exp

!===============================================================================
!> Diagnostics: Number of patches and GB of mem
!===============================================================================
SUBROUTINE diagnostics
  if (mpi%master) then
    print '(a,i8)', ' patches per process:', np
  end if
END SUBROUTINE diagnostics

END MODULE cartesian_mod
