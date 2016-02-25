MODULE cartesian_mod
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
  integer:: dims(3) = [4,4,1]
  real:: origin(3) = [-0.5,-0.5,-0.5]
  real:: size(3) = [1.0,1.0,1.0]
  real:: position(3)
  type(experiment_t), pointer:: patch
  type(link_t) , pointer:: prev=>null(), link
  class(task_t), pointer:: task
  integer:: i, j, k, i1, j1, k1, i2, j2, k2
  type(experiment_t), pointer, dimension(:,:,:):: patches
  namelist /box_params/ size, dims, origin
  !.............................................................................
  call trace_begin ('task_list_t%init_cartesian')
  rewind (io%input); read (io%input, box_params)
  if (io%master) write (*, box_params)
  self%dims = dims
  self%size = size
  allocate (self%task_list)
  allocate (patches(self%dims(1),self%dims(2),self%dims(3)))
  do k=1,self%dims(3)
  do j=1,self%dims(2)
  do i=1,self%dims(1)
    !---------------------------------------------------------------------------
    ! Initialize each patch, via patch_t:init
    !---------------------------------------------------------------------------
    patch => patches(i,j,k)
    patch%box = self%size
    patch%position = ([i,j,k]-0.5_8)*self%size/self%dims + origin
    patch%size = self%size/self%dims
    call patch%init
    if (all(self%dims==1)) call patch%set(bits%root_grid)
    if (io%verbose>1) &
      print'(a,i6,2(2x,3f10.6))', 'init_cartesian: id, size, pos =', &
        patch%id, patch%size, patch%position
    patch%i = i; patch%j = j; patch%k = k
    if (io%verbose > 2) &
      print'(i5,3(3f9.4,3x))', patch%id, patch%position, patch%size
    !---------------------------------------------------------------------------
    ! Add to a list of links, with pointers from patch to patch
    !---------------------------------------------------------------------------
    task => patches(i,j,k)
    call self%task_list%append (task)
  end do
  end do
  end do
  call self%task_list%init (label)
  call trace_end
END SUBROUTINE init

END MODULE cartesian_mod
