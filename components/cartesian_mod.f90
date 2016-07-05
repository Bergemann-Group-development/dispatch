MODULE cartesian_mod
  USE omp_mod
  USE mpi_mod
  USE mpi_coords_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE task_mod
  USE task_list_mod
  USE experiment_mod
  implicit none
  private
  type, public:: cartesian_t
    real(8):: size(3)
    integer:: dims(3)
    type(task_list_t), pointer:: task_list
  contains
    procedure init
  end type
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!> =============================================================================
!> Distribute patches in a Cartesian arrangement in a box
!> =============================================================================
SUBROUTINE init (self, label)
  class(cartesian_t):: self
  character(len=*), optional:: label
  type(task_list_t):: task_list
  integer:: patches_per_mpi(3)
  real:: position(3)
  type(experiment_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link
  class(task_t), pointer:: task
  integer:: i, j, k, i1, j1, k1, i2, j2, k2, rank, dim, nn(2,3)
  type(experiment_t), pointer, dimension(:,:,:):: patches
  !
  real(8):: size(3) = [1.0,1.0,1.0]
  integer:: dims(3) = [4,4,4]
  integer:: mpi_dims(3) = [1,1,1]
  real(8):: origin(3) = [-0.5,-0.5,-0.5]
  logical:: face_nbors=.false.
  namelist /box_params/ size, dims, mpi_dims, origin, face_nbors
  !.............................................................................
  call trace_begin ('cartesian_mod::init')
  rewind (io%input); read (io%input, box_params)
  if (io%master) write (*, box_params)
  self%dims = dims
  self%size = size
  !-----------------------------------------------------------------------------
  ! Generate a (periodic) MPI geometry, with dims specified in the input file
  !-----------------------------------------------------------------------------
  if (mpi%size > 1) then
    call mpi_coords%init (mpi_dims=mpi_dims)
  else
    mpi_dims = 1
  end if
  patches_per_mpi = self%dims/mpi_dims
  allocate (self%task_list)
  call self%task_list%init
  self%task_list%size = self%size
  self%task_list%n_tasks = product(self%dims)
  self%task_list%face_nbors = face_nbors
  do dim=1,3
    nn(:,dim) = merge([-1,1],[0,0],self%dims(dim)>1)
  end do
  !xomp parallel do private(i,j,k,i1,j1,k1,patch,task,rank)
  do k=1,self%dims(3)
  do j=1,self%dims(2)
  do i=1,self%dims(1)
    !---------------------------------------------------------------------------
    ! Initialize each patch, via patch_t:init
    !---------------------------------------------------------------------------
    allocate (patch)
    patch%box = self%size
    patch%position = ([i,j,k]-0.5_8)*self%size/self%dims + origin
    patch%size = self%size/self%dims
    call patch%set(bits%static)
    patch%rank = pos2rank (i, j, k)
    !---------------------------------------------------------------------------
    ! Decide if it is a virtual or boundary patch, by checking the 26 neighbors
    !---------------------------------------------------------------------------
    do k1=nn(1,3),nn(2,3)
    do j1=nn(1,2),nn(2,2)
    do i1=nn(1,1),nn(2,1)
      if (i1==0.and.j1==0.and.k1==0) cycle
      rank = pos2rank (i+i1, j+j1, k+k1)
      if (patch%rank/=mpi%rank .and. rank==mpi%rank) call patch%set (bits%virtual)
      if (patch%rank==mpi%rank .and. rank/=mpi%rank) call patch%set (bits%boundary)
      if (omp_master .and. io%verbose > 2) then
        write (io%log_unit,'(3(3i4,2x),2l6)') &
          i,j,k, i1,j1,k1, mpi%rank,patch%rank,rank, &
          patch%is_set(bits%virtual), patch%is_set(bits%boundary)
        flush (io%log_unit)
      end if
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! If belongs to a different rank and is not a virtual patch, skip it
    !---------------------------------------------------------------------------
    if (patch%rank /= mpi%rank .and. patch%is_clear(bits%virtual)) then
      if (omp_master) then
        write (io%log_unit,*) 'skipping', patch%id
      end if
      deallocate (patch)
    !---------------------------------------------------------------------------
    ! Otherwise, choose a unique patch id, based on position and rank
    !---------------------------------------------------------------------------
    else
      call patch%init
      patch%id = 1 + mpi%rank + &
        mpi%size*((i-1) + self%dims(1)*((j-1) + (k-1)*self%dims(2)))
      patch%id = patch%id - mpi%rank + patch%rank             ! adjust id of vpatch
      patch%i = i; patch%j = j; patch%k = k
      if (all(self%dims==1)) call patch%set(bits%root_grid)
      if (omp_master) then
        if (io%verbose > 2) &
          write (io%log_unit,'(i5,3(3f9.4,3x))') &
            patch%id, patch%position, patch%size
        if (io%verbose>0) then
          write (io%log_unit,'(a,2i6,2(2x,3f10.6),2l6)') &
            'init_cartesian: id, rank, size, pos =', &
            patch%id, patch%rank, patch%size, patch%position, &
            patch%is_set(bits%virtual), patch%is_set(bits%boundary)
          flush (io%log_unit)
        end if
      end if
      task => patch
      call self%task_list%append (task)
    end if
  end do
  end do
  end do
  call trace_end
contains
  integer function pos2rank (i, j, k)
  integer:: i, j, k
    pos2rank = modulo((i-1),self%dims(1))/patches_per_mpi(1) + &
               modulo((j-1),self%dims(2))/patches_per_mpi(2)*mpi_dims(1) + &
               modulo((k-1),self%dims(3))/patches_per_mpi(3)*mpi_dims(1)*mpi_dims(2)
  end function
END SUBROUTINE init

END MODULE cartesian_mod
