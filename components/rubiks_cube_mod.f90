!> =============================================================================
!> Create a nested set of patches, in a hierarchy of Rubik's cubes, where the
!> innermost patch has a given size, and there are n_levels-1 sets of larger
!> Rubik's cubes, each three times larger than the previous, with the previous
!> sets of patches inside the center patch.
!> =============================================================================
MODULE rubiks_cube_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE task_list_mod
  USE experiment_mod
  implicit none
  private
  type, public:: rubiks_cube_t
    real(8):: size(3) = [1.0, 1.0, 1.0]
    real(8):: position(3) = [0.0, 0.0, 0.0]
    integer:: n_levels = 5
    type(task_list_t), pointer:: task_list
  contains
    procedure init
  end type
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!> =============================================================================
!> Distribute patches in a Rubik's cube arrangement around a box
!> =============================================================================
SUBROUTINE init (self, label)
  class(rubiks_cube_t):: self
  character(len=*), optional:: label
  !.............................................................................
  type(task_list_t):: task_list
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: n_levels, iostat
  real(8):: position(3)
  real(8):: size(3)
  integer:: i, j, k, i1, j1, k1, i2, j2, k2, level
  namelist /rubiks_params/ size, position, n_levels
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from 
  ! existing data type values
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_cube_t%init')
  size = self%size
  position = self%position
  n_levels = self%n_levels
  rewind (io%input); read (io%input, rubiks_params, iostat=iostat)
  if (io%master) write (*, rubiks_params)
  self%n_levels = n_levels
  self%position = position
  self%size = size
  allocate (self%task_list)
  do level=1,n_levels
  do k=1,3
  do j=1,3
  do i=1,3
    !---------------------------------------------------------------------------
    ! Initialize each patch, via patch_t:init
    !---------------------------------------------------------------------------
    allocate (experiment_t:: task)
    select type (task)
    class is (experiment_t)
    patch => task
    end select
    patch%box = self%size*3d0**(n_levels+1)
    patch%size = self%size*3d0**(level-1)
    patch%position = ([i,j,k]-0.5_8)*patch%size + position
    call patch%init
    if (io%verbose>1) &
      print'(a,i6,i4,2(2x,3f10.6))', 'rubiks_cube%init: id, level, size, pos =', &
        patch%id, level, patch%size, patch%position
    !---------------------------------------------------------------------------
    ! Add to a list of links, with pointers from patch to patch
    !---------------------------------------------------------------------------
    call self%task_list%append (task)
  end do
  end do
  end do
  end do
  call self%task_list%init (label)
  call trace_end
END SUBROUTINE init

END MODULE rubiks_cube_mod
