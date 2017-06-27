MODULE cartesian_mod
  USE omp_mod
  USE mpi_mod
  USE mpi_coords_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE task_list_mod
  USE experiment_mod
  USE timer_mod
  implicit none
  private
  type, public:: cartesian_t
    real(8):: size(3)
    integer:: dims(3)
    type(task_list_t), pointer:: task_list
  contains
    procedure init
    procedure, nopass:: diagnostics
  end type
  integer, save:: np=0                  ! patches per process
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!> =============================================================================
!> Distribute patches in a Cartesian arrangement in a box
!> =============================================================================
SUBROUTINE init (self, label)
  class(cartesian_t):: self
  character(len=*), optional:: label
  type(task_list_t):: task_list
  type(task_list_t), pointer:: patch_list
  class(link_t), pointer:: link, next, nbor
  integer:: patches_per_mpi(3)
  real:: position(3)
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i, j, k, i1, j1, k1, rank, nbrank, dim, nn(2,3)
  !
  type links_t
    type(link_t), pointer:: p
  end type
  type(links_t), dimension(:,:,:), pointer:: links
  real(8):: size(3) = [1.0,1.0,1.0]
  integer:: dims(3) = [4,4,4]
  integer:: mpi_dims(3) = [1,1,1]
  integer:: per_rank(3) = [0,0,0]
  integer:: npatch, nvpatch, nbpatch, nlpatch
  real(8):: origin(3) = [-0.5,-0.5,-0.5]
  logical:: face_nbors=.false.
  integer:: thread, icoords(3),ipos(3)
  integer, save:: itimer=0
  namelist /box_params/ size, dims, mpi_dims, per_rank, origin, face_nbors
  !.............................................................................
  call trace_begin ('cartesian_mod::init',itimer=itimer)
  rewind (io%input); read (io%input, box_params)
  if (any(per_rank/=0)) then
    where (per_rank/=0) dims=per_rank*mpi_dims
  end if
  per_rank = dims/mpi_dims
  if (io%master) write (*, box_params)
  self%dims = dims
  self%size = size
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
  !-----------------------------------------------------------------------------
  ! Begin on task list
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  call self%task_list%init
  self%task_list%size = self%size
  self%task_list%n_tasks = product(self%dims)
  self%task_list%face_nbors = face_nbors
  do dim=1,3
    nn(:,dim) = merge([-1,1],[0,0],self%dims(dim)>1)
  end do
  npatch=0; nvpatch=0; nbpatch=0; nlpatch=0
  !-----------------------------------------------------------------------------
  ! Allocate an array of links, pointing to patches, and give the patches each
  ! id, rank, size, position, integer position, and status bits.  Discard
  ! external patches and add the rest to the task_list.
  !-----------------------------------------------------------------------------
  allocate (links(0:self%dims(1)-1,0:self%dims(2)-1,0:self%dims(3)-1))
  do k=0,self%dims(3)-1
  do j=0,self%dims(2)-1
  do i=0,self%dims(1)-1
    allocate (link, patch); link%task => patch; patch%link => link
    task => link%task
    links(i,j,k)%p => null()
    patch%ipos = [i,j,k]
    patch%rank = mpi_coords%coords_to_rank ([i,j,k]/patches_per_mpi)
    patch%id = 1+i+self%dims(1)*(j+k*self%dims(2))
    patch%box = self%size
    patch%size = self%size / self%dims
    patch%position = ([i,j,k]+0.5d0)*patch%size + origin
    patch%llc_cart = patch%position - 0.5 * patch%size
    patch%llc_nat = patch%llc_cart
    call patch%set(bits%static)
    if (all(self%dims==1)) call task%set(bits%root_grid)
    rank = patch%rank
    if (rank==mpi%rank) then
      call patch%set(bits%internal)
    else
      call patch%set(bits%external)
    end if
    !---------------------------------------------------------------------------
    ! Set the patch status bits, by checking the ranks of neighbor patches
    !---------------------------------------------------------------------------
    do k1=nn(1,3),nn(2,3)
    do j1=nn(1,2),nn(2,2)
    do i1=nn(1,1),nn(2,1)
      if (i1==0.and.j1==0.and.k1==0) cycle
      ipos = modulo([i+i1,j+j1,k+k1],self%dims)
      nbrank = mpi_coords%coords_to_rank(ipos/patches_per_mpi)
      if (rank==mpi%rank .and. nbrank/=mpi%rank) then
        call patch%set(bits%boundary)
        call patch%clear(bits%internal)
      end if
      if (rank/=mpi%rank .and. nbrank==mpi%rank) then
        call patch%set(bits%virtual)
        call patch%clear(bits%external)
      end if
    end do
    end do
    end do
    if (rank == mpi%rank) nlpatch=nlpatch+1
    if (patch%is_set(bits%external)) then
      deallocate (patch, link)
    else
      links(i,j,k)%p => link
      task => patch
      call self%task_list%append_link (link)
      if (patch%is_set(bits%boundary))  nbpatch = nbpatch+1
      if (patch%is_set(bits%virtual ))  nvpatch = nvpatch+1
      npatch = npatch+1
    end if
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Add an nbor list to each patch.  The modulo below assumes a period mesh in
  ! all directions.  If this is not the case on one or more directions, the
  ! nbor lists must be adjusted accordingly
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    class is (experiment_t)
      patch => task
    class default
      error stop 'xx'
    end select
    do k1=nn(1,3),nn(2,3)
    do j1=nn(1,2),nn(2,2)
    do i1=nn(1,1),nn(2,1)
      if (i1==0.and.j1==0.and.k1==0) cycle
      ipos = modulo(patch%ipos+[i1,j1,k1],self%dims)
      if (associated(links(ipos(1),ipos(2),ipos(3))%p)) then
        allocate (nbor)
        nbor%task => links(ipos(1),ipos(2),ipos(3))%p%task
        nbor%link => links(ipos(1),ipos(2),ipos(3))%p
        nbor%nbor => nbor%link
        call link%add_nbor_by_rank (nbor)
      end if
    end do
    end do
    end do
    link => link%next
  end do
  deallocate (links)
  !-----------------------------------------------------------------------------
  ! Initialize patches on our rank and virtual patches, in parallel OMP tasks
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  !$omp parallel
  !$omp master
  do while (associated(link))
    !xomp task firstprivate(link) private(task,patch) shared(np)
    task => link%task
    select type (task)
    class is (experiment_t)
      patch => task
      call timer%end (itimer=itimer)
      call patch%init
      call timer%begin (itimer=itimer)
      !$omp atomic
      np = np + 1
      if (io%verbose>0) then
        write (io_unit%log,1) &
          'init_cartesian: id, rank, size, pos =', &
          patch%id, patch%rank, patch%size, patch%position, &
          patch%is_set(bits%boundary), patch%is_set(bits%virtual)
        1 format(a,i8,i6,1p,2(2x,3e14.6),2x,2l1)
        flush (io_unit%log)
      end if
      if (io_unit%verbose>1 .and. io_unit%master) then
        write (io_unit%output,1) &
          'init_cartesian: id, rank, size, pos =', &
          patch%id, patch%rank, patch%size, patch%position, &
          patch%is_set(bits%boundary), patch%is_set(bits%virtual)
        flush (io_unit%log)
      end if
    end select
    !xomp end task
    link => link%next
  end do
  !$omp end master
  !$omp end parallel
  !call task_list%reset_status
  call self%diagnostics
  write (io_unit%mpi,*) npatch,' patches',nvpatch,' virtual',nbpatch,' boundary'
  write (io_unit%mpi,*) nlpatch,' local patches'
  flush (io%log_unit)
  call self%task_list%print_tasks
  call trace_end (itimer)
END SUBROUTINE init

!> =============================================================================
!> Diagnostics: Number of patches and GB of mem
!> =============================================================================
SUBROUTINE diagnostics
  if (mpi%master) then
    print '(a,i8)',     ' Patches per process:', np
  end if
END SUBROUTINE

END MODULE cartesian_mod
