MODULE cartesian_mod
  USE omp_mod
  USE mpi_mod
  USE mpi_coords_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE task_mod
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
  real, save:: gb=0.0                   ! GB per process
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
  class(link_t), pointer:: link, next
  integer:: patches_per_mpi(3)
  real:: position(3)
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i, j, k, i1, j1, k1, i2, j2, k2, rank, dim, nn(2,3)
  !
  real(8):: size(3) = [1.0,1.0,1.0]
  integer:: dims(3) = [4,4,4]
  integer:: mpi_dims(3) = [1,1,1]
  integer:: npatch, nvpatch, nbpatch
  real(8):: origin(3) = [-0.5,-0.5,-0.5]
  logical:: face_nbors=.false.
  integer:: thread, icoords(3)
  integer, save:: itimer=0
  namelist /box_params/ size, dims, mpi_dims, origin, face_nbors
  !.............................................................................
  call trace_begin ('cartesian_mod::init',itimer=itimer)
  rewind (io%input); read (io%input, box_params)
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
  allocate (self%task_list)
  call self%task_list%init
  self%task_list%size = self%size
  self%task_list%n_tasks = product(self%dims)
  self%task_list%face_nbors = face_nbors
  do dim=1,3
    nn(:,dim) = merge([-1,1],[0,0],self%dims(dim)>1)
  end do
  npatch=0; nvpatch=0; nbpatch=0
  !-----------------------------------------------------------------------------
  ! Find the location of all patches, based on the rank_to_coords procedure
  !-----------------------------------------------------------------------------
  allocate (patch_list)
  call patch_list%init ('patch_list')
  do rank=0,mpi_coords%size-1
    do k=0,patches_per_mpi(3)-1
    do j=0,patches_per_mpi(2)-1
    do i=0,patches_per_mpi(1)-1
      !------------------------------------------------------------------------
      ! Allocate all patches, give them a unique id, and find their position,
      ! counted in integer multiples of their size by icoords, and centered
      ! relative to the full box by origin.
      !---------------------------------------------------------------------------
      allocate (patch)
      patch%rank = rank
      icoords = mpi_coords%rank_to_coords(rank)*patches_per_mpi + [i,j,k]
      patch%id = 1 + icoords(1) + self%dims(1)*(icoords(2) + icoords(3)*self%dims(2))
      patch%i = i; patch%j = j; patch%k = k
      patch%box = self%size
      patch%size = self%size / self%dims
      patch%position = (icoords+0.5d0)*patch%size + origin
      patch%llc = patch%position - 0.5 * patch%size
      patch%origin = patch%llc
      call patch%set(bits%static)
      task=> patch
      call patch_list%append (task)
    end do
    end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Set the boundary and virtual bits, based on overlap with patches on our rank
  !-----------------------------------------------------------------------------
  call patch_list%reset_status
  !-----------------------------------------------------------------------------
  ! Initialize patches on our rank and virtual patches. This needs to be done
  ! without changing the link chain in the patch_list for other threads.
  !-----------------------------------------------------------------------------
  !xomp parallel private(link,thread,task)
  link => patch_list%head
  thread = 0
  do while (associated(link))
    task => link%task
    if ((task%rank == mpi%rank .or. task%is_set(bits%virtual)) &
        .and. thread==omp%thread) then
      !print'(a,5i6)','task%id, thread=', task%id, thread, omp%nthreads
      select type (task)
      type is (experiment_t)
        call task%init
        call trace_begin ('cartesian_mod::init',itimer=itimer)
        !$omp atomic
        gb = gb + product(task%gn)*task%nv*task%nt*task%nw*4.0/1024.**3
        !$omp atomic
        np = np + 1
        if (all(self%dims==1)) call task%set(bits%root_grid)
        if (io%verbose>0) then
          write (io_unit%log,1) &
            'init_cartesian: id, rank, size, pos =', &
            task%id, task%rank, task%size, task%position, &
            task%is_set(bits%boundary), task%is_set(bits%virtual)
          1 format(a,i8,i6,1p,2(2x,3e14.6),2x,2l1)
          flush (io_unit%log)
        end if
        if (io_unit%verbose>1 .and. io_unit%master) then
          write (io_unit%output,1) &
            'init_cartesian: id, rank, size, pos =', &
            task%id, task%rank, task%size, task%position, &
            task%is_set(bits%boundary), task%is_set(bits%virtual)
          flush (io_unit%log)
        end if
      end select
    end if
    if (omp_in_parallel()) &
      thread = mod(thread+1,omp%nthreads)
    link => link%next
  end do
  call self%diagnostics
  !xomp barrier
  !xomp end parallel
  !-----------------------------------------------------------------------------
  ! Append initialized patches to the task_list, and deallocate unused patches
  ! and links, and the patch_list.
  !-----------------------------------------------------------------------------
  link => patch_list%head
  do while (associated(link))
    task => link%task
    next => link%next
    if (task%rank == mpi%rank .or. task%is_set(bits%virtual)) then
      call self%task_list%append_link (link)
    else
      deallocate (patch)
      deallocate (link)
    end if
    link => next
  end do
  deallocate (patch_list)
  write (io_unit%log,*) npatch,' patches',nvpatch,' virtual',nbpatch,' boundary'
  flush (io%log_unit)
  !-----------------------------------------------------------------------------
  ! Construct neighbor lists used to perform checks for ready status
  !-----------------------------------------------------------------------------
  call self%task_list%init_all_nbors
  call trace_end (itimer)
END SUBROUTINE init

!> =============================================================================
!> Diagnostics: Number of patches and GB of mem
!> =============================================================================
SUBROUTINE diagnostics
  if (mpi%master) then
    print '(a,i6)',     'Patches per process:', np
    print '(a,f6.2,a)', ' Memory per process:', gb, ' GB'
  end if
END SUBROUTINE

END MODULE cartesian_mod
