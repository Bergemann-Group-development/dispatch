MODULE volleyball_mod
  USE task_list_mod
  USE experiment_mod
  USE patch_mod
  USE gpatch_mod
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
  USE sphere_point_mod
  USE dll_mod
  USE os_mod
  implicit none
  private
  type, public, extends(task_t):: rank_t
  end type
  type, public:: volleyball_t
    real(8):: size(3)
    integer:: dims(3)
    class(task_list_t), pointer:: task_list
  contains
    procedure:: init
    procedure, nopass:: diagnostics
  end type
  integer:: verbose=0, dbglog=0
  logical:: omp_init=.false.
  type(rank_t):: rank
  integer, save:: np=0                  ! patches per process
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!===============================================================================
!> Distribute patches in a volleyball arrangement
!===============================================================================
SUBROUTINE init (self, label)
  class(volleyball_t):: self
  character(len=*), optional:: label
  !.............................................................................
  type(volleyball_layers_t)         :: vb_layers
  type(volleyball_layer_t), pointer :: layer
  class(link_t),            pointer :: link, next
  class(experiment_t),      pointer :: patch
  class(task_t),            pointer :: task
  type(sphere_point_t),     pointer :: point
  class(dll_node_t),        pointer :: node1, node2
  class(*),                 pointer :: ptr
  integer :: i, j, k, id, ip, n_freeze
  integer :: dims(3) = [4,4,4]
  integer :: patches_per_mpi(3)
  integer :: mpi_dims(3) = [1,1,1]
  integer :: per_rank(3) = [0,0,0]
  real(8) :: size(3) = [1.0,1.0,1.0]
  real(8) :: origin(3) = [-0.5,-0.5,-0.5]
  real(8) :: wc, r_min, radius, e(3,3)
  integer :: verbose = 1
  integer :: itimer = 0
  namelist /volleyball_params/ size, dims, mpi_dims, per_rank, origin, &
                               omp_init, verbose
  character(len=120):: ids = &
    '$Id$ components/cartesian_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('volleyball_t%init')
  call io%header('begin volleyball_t%init: Cartesian patch arrangement.')
  call trace%print_id (ids)
  call timer%tic ()
  !-----------------------------------------------------------------------------
  rewind (io%input); read (io%input, volleyball_params)
  dbglog = merge (io%output, io_unit%log, verbose > 1)
  !-----------------------------------------------------------------------------
  ! MPI setup
  !-----------------------------------------------------------------------------
  if (any(per_rank /= 0)) then
    where (per_rank /= 0) dims=per_rank*mpi_dims
  end if
  per_rank = dims/mpi_dims
  self%dims = dims
  self%size = size
  if (mpi%size==1) mpi_dims = 1
  !-----------------------------------------------------------------------------
  ! Generate an MPI geometry, with dims specified in the input file,
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
  !-----------------------------------------------------------------------------
  ! Set attritbutes of the rank_t data type -- just enough to use distance().
  !-----------------------------------------------------------------------------
  rank%periodic = .true.
  rank%box = self%size
  rank%size = rank%box/mpi_dims
  rank%ds = rank%size/patches_per_mpi
  rank%position = (mpi_coords%rank_to_coords(mpi%rank) + 0.5_8)*rank%size
  !-----------------------------------------------------------------------------
  ! Begin on task list
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  call os%mem_register (storage_size(self%task_list), 1, 'task_list')
  call self%task_list%init ('tlist')
  self%task_list%size = self%size
  self%task_list%dims = self%dims
  self%task_list%n_tasks = product(self%dims)
  !-----------------------------------------------------------------------------
  ! Generate centers and orientations of all patches, stored in vb_layers_t%list
  !-----------------------------------------------------------------------------
  call vb_layers%init ()
  write (dbglog,'(a,i4)') ' volleyball_t%init: layers=', vb_layers%list%n
  !-----------------------------------------------------------------------------
  ! Loop over layers -- starting at the surface
  !-----------------------------------------------------------------------------
  node1 => vb_layers%list%head
  do while (associated(node1))
    layer => to_layer (node1)
    !---------------------------------------------------------------------------
    ! Loop over patches in layer
    !---------------------------------------------------------------------------
    nullify (layer%link)
    node2 => layer%list%head
    do while (associated(node2))
      point => to_point (node2)
      !-------------------------------------------------------------------------
      ! Initially append all tasks to the task list, to build nbor lists
      !-------------------------------------------------------------------------
      allocate (patch)
      call os%mem_register (storage_size(patch), 1, 'patch')
      task => patch
      call self%task_list%prepare (task)
      patch%layer        = layer%number
      patch%box          = 1d0
      patch%size         = point%size
      radius             = sum(point%p**2)**0.5
      patch%llc_cart     = [0d0,0d0,radius] - 0.5*point%size
      patch%llc_nat      = patch%llc_cart
      patch%radius       = radius
      patch%position     = point%p
      call patch%pre_init()
      call patch%task_t%init_unique (same=.true.)
      call self%task_list%add(patch%link)
      if (.not.associated(layer%link)) then
        layer%link => patch%link
      end if
      !-------------------------------------------------------------------------
      ! Set standard patch attributes, so far in terms of solar radius units
      !-------------------------------------------------------------------------
      patch%periodic     = .false.
      patch%origin       = origin
      patch%box          = 1d0
      patch%size         = point%size
      patch%ds           = point%size/5     ! No need to match the final value
      patch%position     = point%p
      patch%centre_nat   = point%p          ! Necessary for patch%distance()!t
      patch%point        => point
      radius             = sum(point%p**2)**0.5
      patch%llc_cart     = [0d0,0d0,radius] - 0.5*point%size
      patch%llc_nat      = patch%llc_cart
      patch%radius       = radius
      patch%r_min = vb_layers%r_min
      patch%r_max = vb_layers%r_max
      call patch%set(bits%static)
      if (.not.associated(layer%prev)) call patch%boundaries%set (bits%zu)
      if (.not.associated(layer%next)) call patch%boundaries%set (bits%zl)
      patch%boundaries%set_bc_bits=.false.
      !-------------------------------------------------------------------------
      ! Initialize the three unit vector: e(:,1)=East, e(:,2)=North, e(:,3)=Up
      !-------------------------------------------------------------------------
      call patch%point%unit_vectors ()
      call os%allocate (patch%erot, [3,3])
      patch%erot = patch%point%e
      !-------------------------------------------------------------------------
      ! Make the node%car point to the task list instead, to use in loops
      !-------------------------------------------------------------------------
      node2%car => patch%link
      node2 => node2%next
    end do
    node1 => node1%next
  end do
  self%task_list%hands_off_nbor_lists = .true.
  call timer%toc ('volleyball_t%init: tasks allocate')
  !-----------------------------------------------------------------------------
  ! Generate nbor lists and set status bits for all tasks
  !-----------------------------------------------------------------------------
  call vb_layers%init_all_nbors (self%task_list)
  call timer%toc ('volleyball_t%init: init_all_nbors')
  call self%task_list%reset_status
  call timer%toc ('volleyball_t%init: reset_status')
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
      deallocate (link%task)
      deallocate (link)
    end if
    link => next
  end do
  call timer%toc ('volleyball_t%init: remove nbor lists')
  !-----------------------------------------------------------------------------
  ! Count task statistics
  !-----------------------------------------------------------------------------
  call self%task_list%count_status
  call timer%toc ('volleyball_t%init: count status')
  !-----------------------------------------------------------------------------
  ! These must be set before (possibly) reading input snapshots
  !-----------------------------------------------------------------------------
  io%ntask = self%task_list%na
  io%nwrite = self%task_list%na
  io%ntotal = product(self%dims)
  write (stdout,*) 'Number of tasks in task list:', self%task_list%n
  write (stdout,*) 'ntask, nwrite, ntotal =', io%ntask, io%nwrite, io%ntotal
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
  write (stdout,*) np, ' patches'
  call timer%toc ('volleyball_t%init: init_exp')
  !-----------------------------------------------------------------------------
  ! Save namelist
  !-----------------------------------------------------------------------------
  write (io%namelist, volleyball_params)
  !-----------------------------------------------------------------------------
  ! These calls are not valid until the patch mem and meshes have been allocated
  !-----------------------------------------------------------------------------
  call self%task_list%init_bdries
  call timer%toc ('volleyball_t%init: init boundaries')
  call self%diagnostics()
  call timer%toc ('volleyball_t%init: diagnostics')
  call self%task_list%print_tasks
  call timer%toc ('volleyball_t%init: print tasks')
  !-----------------------------------------------------------------------------
  flush (stdout)
  link => self%task_list%head
  r_min = huge(1.)
  do while (associated(link))
    radius = sqrt(sum(link%task%position**2))
    r_min = min (r_min, radius)
    link => link%next
  end do
  call timer%toc ('volleyball_t%init: r_min')
  flush (stdout)
  n_freeze = 0
  link => self%task_list%head
  do while (associated(link))
    radius = sqrt(sum(link%task%position**2))
    if (radius < 1.05*r_min) then
      !call link%task%set (bits%frozen)
      !n_freeze = n_freeze+1
    end if
    link => link%next
  end do
  call timer%toc ('volleyball_t%init: freeze bottom (off)')
  write (stdout,*) 'r_min, n_freeze =', r_min, n_freeze
  call trace%end ()
contains
!===============================================================================
function to_layer (node) result (layer)
  class(dll_node_t),         pointer:: node
  class(volleyball_layer_t), pointer:: layer
  class(*),                  pointer:: ptr
  !-----------------------------------------------------------------------------
  ptr => node%car
  select type (ptr)
  class is (volleyball_layer_t)
    layer => ptr
  end select
end function to_layer
!===============================================================================
function to_point (node) result (point)
  class(dll_node_t),     pointer:: node
  class(sphere_point_t), pointer:: point
  class(*),              pointer:: ptr
  !-----------------------------------------------------------------------------
  ptr => node%car
  select type (ptr)
  class is (sphere_point_t)
    point => ptr
  end select
end function to_point
END SUBROUTINE init

!===============================================================================
! Set initial values from a 1-D model, eg Model S (defined in experiment)
!===============================================================================
SUBROUTINE init_exp (self, origin)
  class(volleyball_t):: self
  real(8):: origin(3)
  !.............................................................................
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  integer:: id
  !-----------------------------------------------------------------------------
  call trace%begin('volleyball_t%init_exp')
  call timer%tic ()
  id = 0
  link => self%task_list%head
  do while (associated(link))
    id = id+1
    if (.not.omp_init .or. mod(id,omp%nthreads)==omp%thread) then
      task => link%task
      select type (task)
      class is (gpatch_t)
        !-----------------------------------------------------------------------
        ! Initialize, which includes interpolating from a 1-D model
        !-----------------------------------------------------------------------
        call task%init ()
        !$omp atomic
        np = np + 1
        if (mod(np,5000) == 0) &
          write (stdout,*) np, ' patches'
      end select
    end if
    link => link%next
  end do
  call timer%toc ('time to initialize patches')
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

END MODULE volleyball_mod
