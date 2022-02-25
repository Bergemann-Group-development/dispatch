MODULE cartesian_vb_mod
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
  USE volleyball_layers_mod
  USE sphere_points_mod
  USE dll_mod
  implicit none
  private
  type, public, extends(task_t):: rank_t
  end type
  type, public:: cartesian_vb_t
    real(8):: position(3)
    real(8):: size(3)
    integer:: dims(3)
    type(task_list_t), pointer:: task_list => null()
  contains
    procedure init
    procedure, nopass:: diagnostics
  end type
  integer:: verbose=0
  logical:: omp_init=.false.
  type(rank_t):: rank
  integer, save:: np=0                  ! patches per process
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!===============================================================================
!> Distribute patches in a volley ball arrangement
!===============================================================================
SUBROUTINE init (self, task_list, label)
  class(cartesian_vb_t):: self
  class(task_list_t), optional, pointer:: task_list
  character(len=*), optional:: label
  !.............................................................................
  class(task_list_t), pointer:: patch_list
  class(list_t), pointer:: tlist
  class(link_t), pointer:: link, next
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  type(volleyball_layers_t):: vb_layers
  type(sphere_points_t), pointer:: info
  class(dll_node_t), pointer:: node
  class(*), pointer:: ptr
  integer:: i, j, k, id, ip, iostat
  integer:: dims(3) = [4,4,4]
  integer:: patches_per_mpi(3)
  integer:: mpi_dims(3) = [1,1,1]
  integer:: per_rank(3) = [0,0,0]
  real(8):: size(3) = [1.0,1.0,1.0]
  real(8):: origin(3) = [-0.5,-0.5,-0.5]
  real(8):: position(3) = tiny(1d0)
  real(8):: wc
  logical:: periodic(3)=.true.
  integer, save:: itimer=0
  namelist /cartesian_vb_params/ size, position, dims, mpi_dims, per_rank, origin, &
    omp_init, periodic
  character(len=120):: ids = &
    '$Id$ components/cartesian_mod.f90'
  !-----------------------------------------------------------------------------
  if (present (task_list)) then
    if (associated(task_list)) then
      write (stdlog,*) 'cartesian_t%init: return -- task_list already allocated'
      return
    else
      allocate(task_list)
      self%task_list => task_list
    end if
  else
    allocate(self%task_list)
  end if
tlist => self%task_list
  !-----------------------------------------------------------------------------
  call trace%begin ('cartesian_vb_t%init')
  call io%header('begin cartesian_vb_t%init: Cartesian patch arrangement')
  call trace%print_id (ids)
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, cartesian_vb_params, iostat=iostat)
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
  write (io%namelist, cartesian_vb_params)
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
  call self%task_list%init ('volleyball')
  self%task_list%size = self%size
  self%task_list%dims = self%dims
  self%task_list%position = self%position
  self%task_list%n_tasks = product(self%dims)
  !-----------------------------------------------------------------------------
  ! Generate centers and orientations of all patches, stored in vb_layers_t%list
  !-----------------------------------------------------------------------------
  call vb_layers%init
print *, 'layers', vb_layers%list%n, 'list%n'
  node => vb_layers%list%head
  do while (associated(node))
    !---------------------------------------------------------------------------
    ! Initially append all tasks to the task list, to build nbor lists
    !---------------------------------------------------------------------------
    allocate (patch)
    task => patch
    call self%task_list%prepare (task)
    call patch%task_t%init_unique (same=.true.)
    patch%task_list => self%task_list
    !---------------------------------------------------------------------------
    ! Pick up location and orientation
    !---------------------------------------------------------------------------
    ptr => node%car
    select type (ptr)
    class is (sphere_points_t)
      info => ptr
    end select
    !---------------------------------------------------------------------------
    ! Set standard patch attributes
    !---------------------------------------------------------------------------
    patch%periodic     = .false.
    patch%origin       = origin
    patch%box          = self%size
    patch%size         = info%size
    patch%ds           = info%size/16                ! FIXME
    patch%position     = [info%x, info%y, info%z]
!print *, patch%id, patch%position, 'pos'
    patch%orient%theta = info%theta
    patch%orient%phi   = info%phi
    patch%orient%p     = patch%position
    patch%orient%on    = .true.
    patch%llc_cart     = patch%position - 0.5 * patch%size
    patch%centre_nat   = patch%position
    patch%llc_nat      = patch%llc_cart
    call patch%set(bits%static)
    !patch%rank = ...
    !---------------------------------------------------------------------------
    ! Setup local coordinate system, with unit vectors
    !---------------------------------------------------------------------------
    ! Up
    patch%orient%e(:,3) = normalize(patch%position)
    if (patch%position(3) > 0d0) then
      patch%orient%e(:,2) = [0d0,0d0,1d0] - patch%position
    else
      patch%orient%e(:,2) = patch%position - [0d0,0d0,-1d0]
    end if
    ! East
    patch%orient%e(:,1) = normalize(cross(patch%orient%e(:,2),patch%orient%e(:,3)))
    ! North
    patch%orient%e(:,2) = cross(patch%orient%e(:,3),patch%orient%e(:,1))
! -- finalize patch
    call task%pre_init()
    call self%task_list%add (patch%link)
    call patch%init_task_list (tlist)
    node => node%next
  end do
  !-----------------------------------------------------------------------------
  ! Generate nbor lists and set status bits for all tasks
  !-----------------------------------------------------------------------------
  if (io%verbose > 0) write(stdout,*) 'cartesian_vb_t%init: init_all_nbors'
  call self%task_list%init_all_nbors
  if (io%verbose > 0) write(stdout,*) 'cartesian_vb_t%init: reset_status'
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
  if (verbose > 0) write(stdout,*) 'cartesian_vb_t%init: count status'
  call self%task_list%count_status
  !-----------------------------------------------------------------------------
  ! The I/O conunters must be set before initializing patches
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  write (stdout,'(a,i8)') ' number of tasks in task list:', self%task_list%n
  write (stdlog,*) 'cartesian_vb_mod: ntask,nwrite=', &
    io%ntask, io%nwrite
  !-----------------------------------------------------------------------------
  ! Initialize tasks
  !-----------------------------------------------------------------------------
  if (omp_init) then
    !$omp parallel default(shared)
    call init_exp (self, origin)
    !$omp end parallel
  else
    call init_exp (self, origin)
  end if
  !-----------------------------------------------------------------------------
  ! Refresh the counts, now with bits%has_output set on tasks
  !-----------------------------------------------------------------------------
  call self%task_list%set_io_counts
  write (stdlog,*) 'cartesian_mod: ntask,nwrite=', &
    io%ntask, io%nwrite
  !-----------------------------------------------------------------------------
  ! This call is not valid until the patch mem and meshes have been allocated
  !-----------------------------------------------------------------------------
  if (verbose > 0) &
    write(stdout,*) 'cartesian_vb_t%init: init boundaries'
  call self%task_list%init_bdries
  !-----------------------------------------------------------------------------
  if (io_unit%master) then
    write (io%namelist, cartesian_vb_params)
    flush (io%namelist)
  end if
  if (verbose > 0) &
    write(stdout,*) 'cartesian_vb_t%init: diagnostics'
  call self%diagnostics()
  if (verbose > 0) &
    write(stdout,*) 'cartesian_vb_t%init: print tasks'
  call self%task_list%print_tasks
  call self%task_list%check_hash_table
  if (verbose > 0) &
    write(stdout,*) 'cartesian_vb_t%init: done'
  flush (stdout)
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

SUBROUTINE init_exp (self, origin)
  class(cartesian_vb_t):: self
  real(8):: origin(3)
  !.............................................................................
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  class(experiment_t), pointer:: patch
  integer:: id
  !-----------------------------------------------------------------------------
  call trace%begin('cartesian_vb_t%init_exp')
  id = 0
  link => self%task_list%head
  do while (associated(link))
    id = id+1
    if (.not.omp_init .or. mod(id,omp%nthreads)==omp%thread) then
      task => link%task
      select type (task)
      class is (experiment_t)
        patch => task
        call patch%init
        !-----------------------------------------------------------------------
        ! Redefine the mesh position, to a locally vertical coordinate system
        !-----------------------------------------------------------------------
        patch%mesh(1)%p = 0d0
        patch%mesh(2)%p = 0d0
        patch%mesh(3)%p = sqrt(sum(patch%position**2))
        patch%mesh%origin = origin
        !$omp atomic
        np = np + 1
        if (io%verbose>0) then
          write (io_unit%log,1) &
            'cartesian_vb_t%init: id, rank, size, pos =', &
            patch%id, patch%rank, patch%size, patch%position, &
            patch%is_set(bits%boundary), patch%is_set(bits%virtual)
          1 format(a,i8,i6,1p,2(2x,3e14.6),2x,2l1,2x,a)
          flush (io_unit%log)
        end if
        if (patch%id == io%id_debug) then
          write (io_unit%output,1) &
            'cartesian_vb_t%init: id, rank, size, pos =', &
            patch%id, patch%rank, patch%size, patch%position, &
            patch%is_set(bits%boundary), patch%is_set(bits%virtual), 'DBG'
          flush (io_unit%log)
        else if (io_unit%verbose>1 .and. io_unit%master) then
          write (io_unit%output,1) &
            'cartesian_vb_t%init: id, rank, size, pos =', &
            patch%id, patch%rank, patch%size, patch%position, &
            patch%is_set(bits%boundary), patch%is_set(bits%virtual)
          flush (io_unit%log)
        end if
      end select
    end if
    link => link%next
  end do
  call trace%end
END SUBROUTINE init_exp

!===============================================================================
!> Diagnostics: Number of patches and GB of mem
!===============================================================================
SUBROUTINE diagnostics
  if (mpi%master) then
    print '(a,i8)',     ' Patches per process:', np
  end if
END SUBROUTINE diagnostics

!===============================================================================
!> Utility functions for single 3-D vector operations
!===============================================================================
FUNCTION dot (a, b) RESULT (c)
  real(8):: a(3), b(3), c
  c = sum(a*b)
END FUNCTION dot
!===============================================================================
FUNCTION normalize (a) RESULT (b)
  real(8):: a(3), b(3)
  b = a/sqrt(sum(a**2))
END FUNCTION normalize
!===============================================================================
FUNCTION cross (a, b) RESULT (c)
  real(8):: a(3), b(3), c(3)
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
END FUNCTION cross

END MODULE cartesian_vb_mod
