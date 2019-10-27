MODULE cartesian_mod
  USE task_list_mod
  USE experiment_mod
  USE patch_mod
  USE link_mod
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
    real(8):: size(3)
    integer:: dims(3)
    type(task_list_t), pointer:: task_list
  contains
    procedure init
    procedure, nopass:: diagnostics
  end type
  logical:: omp_init=.false.
  type(rank_t):: rank
  integer, save:: np=0                  ! patches per process
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!===============================================================================
!> Distribute patches in a Cartesian arrangement in a box
!===============================================================================
SUBROUTINE init (self, label)
  class(cartesian_t):: self
  character(len=*), optional:: label
  !.............................................................................
  type(task_list_t), pointer:: patch_list
  class(link_t), pointer:: link, next
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i, j, k, id, ip
  integer:: dims(3) = [4,4,4]
  integer:: patches_per_mpi(3)
  integer:: mpi_dims(3) = [1,1,1]
  integer:: per_rank(3) = [0,0,0]
  real(8):: size(3) = [1.0,1.0,1.0]
  real(8):: origin(3) = [-0.5,-0.5,-0.5]
  real(8):: wc
  logical:: face_nbors=.false.
  logical:: periodic(3)=.true.
  integer, save:: itimer=0
  namelist /cartesian_params/ size, dims, mpi_dims, per_rank, origin, face_nbors, &
    omp_init, periodic
  character(len=120):: ids = &
    '$Id$ components/cartesian_mod.f90'
  !-----------------------------------------------------------------------------
  call trace%begin ('cartesian_t%init')
  call io%header('begin cartesian_t%init: Cartesian patch arrangement')
  call trace%print_id (ids)
  !-----------------------------------------------------------------------------
  rewind (io%input); read (io%input, cartesian_params)
  write (io%output, cartesian_params)
  if (any(per_rank/=0)) then
    where (per_rank/=0) dims=per_rank*mpi_dims
  end if
  per_rank = dims/mpi_dims
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
  call self%task_list%init ('tlist')
  self%task_list%size = self%size
  self%task_list%dims = self%dims
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
    patch%box = self%size
    patch%size = self%size / self%dims
    patch%position = ([i,j,k]+0.5d0)*patch%size + origin
    if (any(abs(rank%distance (patch)) > 0.5*(rank%size + 1.5*rank%ds))) then
      deallocate (link, patch)
      cycle
    end if
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
    call self%task_list%append_link (link)
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Generate nbor lists and set status bits for all tasks
  !-----------------------------------------------------------------------------
  if (io%verbose > 0) write(stdout,*) 'cartesian_t%init: init_all_nbors'
  call self%task_list%init_all_nbors
  if (io%verbose > 0) write(stdout,*) 'cartesian_t%init: reset_status'
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
      deallocate (link%task)
      deallocate (link)
    end if
    link => next
  end do
  !-----------------------------------------------------------------------------
  ! This call is not valid until the patch meshes have been allocated
  !-----------------------------------------------------------------------------
  if (io%verbose > 0) write(stdout,*) 'cartesian_t%init: count status'
  call self%task_list%count_status
  !-----------------------------------------------------------------------------
  ! These must be set before reading input snapshots
  !-----------------------------------------------------------------------------
  io%ntask = self%task_list%na
  io%nwrite = self%task_list%na
  io%ntotal = product(self%dims)
  write (stdout,*) 'Number of tasks in task list:', self%task_list%n
  write(stdout,*) 'ntask, nwrite, ntotal =', io%ntask, io%nwrite, io%ntotal
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
  ! This call is not valid until the patch mem and meshes have been allocated
  !-----------------------------------------------------------------------------
  if (io%verbose > 0) &
    write(stdout,*) 'cartesian_t%init: init boundaries'
  call self%task_list%init_bdries
  !-----------------------------------------------------------------------------
  if (io_unit%master) then
    write (io_unit%nml, cartesian_params)
    flush (io_unit%nml)
  end if
  if (io%verbose > 0) &
    write(stdout,*) 'cartesian_t%init: diagnostics'
  call self%diagnostics()
  if (io%verbose > 0) &
    write(stdout,*) 'cartesian_t%init: print tasks'
  call self%task_list%print_tasks
  if (io%verbose > 0) &
    write(stdout,*) 'cartesian_t%init: done'
  flush (stdout)
  !-----------------------------------------------------------------------------
  call trace%end ()
END SUBROUTINE init

SUBROUTINE init_exp (self, origin)
  class(cartesian_t):: self
  real(8):: origin(3)
  !.............................................................................
  type(link_t), pointer:: link
  class(task_t), pointer:: task
  class(experiment_t), pointer:: patch
  integer:: id
  !-----------------------------------------------------------------------------
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
        if (omp_init .and. omp%nthreads>1) then
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
        if (io%verbose>0) then
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
    print '(a,i8)',     ' Patches per process:', np
  end if
END SUBROUTINE diagnostics

END MODULE cartesian_mod
