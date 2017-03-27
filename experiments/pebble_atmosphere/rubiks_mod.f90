!> =============================================================================
!> Create a nested set of patches, in a hierarchy of Rubik's cubes, where the
!> innermost patch has a given size, and there are n_levels-1 sets of larger
!> Rubik's cubes, each three times larger than the previous, with the previous
!> sets of patches inside the center patch.
!> =============================================================================
MODULE rubiks_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE task_list_mod
  USE experiment_mod
  USE units_mod
  USE scaling_mod
  USE scene_mod
  USE link_mod
  implicit none
  private
  type, public:: rubiks_t
    real(8):: size(3) = [1.0, 1.0, 1.0]
    real(8):: position(3) = [0.0, 0.0, 0.0]
    integer:: n_levels, n_repeat
    integer:: n_patch, m_zoom
    type(task_list_t), pointer:: task_list
  contains
    procedure init_hierarchy
    procedure init_level
  end type
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!> =============================================================================
!> Distribute patches in a Rubik's cube arrangement around a box
!> =============================================================================
SUBROUTINE init_hierarchy (self, label)
  class(rubiks_t):: self
  character(len=*), optional:: label
  !.............................................................................
  type(task_list_t):: task_list
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  class(link_t), pointer:: link
  integer, save:: n_patch=3, m_zoom=3, n_levels=5
  integer, save:: n_repeat=1
  integer:: iostat
  real(8), save:: position(3)
  real(8), save:: size(3)
  integer:: i, j, k, ir, level, levels
  namelist /rubiks_params/ size, position, n_levels, n_patch, m_zoom, n_repeat
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from
  ! existing data type values
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_t%init')
  size = self%size
  position = self%position
  n_levels = self%n_levels
  rewind (io%input); read (io%input, rubiks_params, iostat=iostat)
  if (io%master) write (*, rubiks_params)
  self%n_levels = n_levels
  self%position = position
  self%size = size
  self%n_patch = n_patch
  self%m_zoom = m_zoom
  !-----------------------------------------------------------------------------
  ! Save the outer size of the hierarchy as the task list size, for distance
  ! measurements.  The self%task_list%size is the full size of the objects in
  ! the task list, and the size of the period in the periodic wrapping.  The
  ! self%size, on the other hand, is the size of the smallest, central patch.
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  self%task_list%size = self%size*self%n_patch*self%m_zoom**(n_levels-1)*n_repeat
  n_repeat = (n_repeat/2)*2 + 1
  scene%n_repeat = n_repeat
  !xomp parallel do private (ir,levels,level)
  do ir=0,n_repeat-1
    position(2) = ((ir-n_repeat/2)*self%task_list%size(2))/n_repeat
    if (ir==n_repeat/2) then
      levels = n_levels
      if (io%verbose>1) &
        print*,'generating a hiearchy of ',n_patch,' patches',' at position', position
    else
      levels = 1
      if (io%verbose>1) &
        print*,'generating a cube of ',n_patch,' patches',' at position', position
    end if
    if (io%verbose>0) &
      print'("repeated Rubiks Cubes: ir, levels, position =",i4,i3,1p,3g14.4)', &
      ir, levels, position
    do level=1,levels
      call self%init_level (level, position)
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Initialize the task list
  !-----------------------------------------------------------------------------
  call self%task_list%init (label)
  call self%task_list%init_all_nbors
  !-----------------------------------------------------------------------------
  ! Initialize the boundary bits
  !-----------------------------------------------------------------------------
  link => self%task_list%head
  print *,'task_list%lc =', self%task_list%lc
  print *,'task_list%uc =', self%task_list%uc
  do while (associated(link))
    task => link%task
    select type (task)
    class is (experiment_t)
    patch => task
    end select
    !---------------------------------------------------------------------------
    ! If at edge, set boundary bits
    !---------------------------------------------------------------------------
    if (patch%position(1) <= self%task_list%lc(1)) call patch%boundaries%set (bits%xl)
    if (patch%position(2) <= self%task_list%lc(2)) call patch%boundaries%set (bits%yl)
    if (patch%position(3) <= self%task_list%lc(3)) call patch%boundaries%set (bits%zl)
    if (patch%position(1) >= self%task_list%uc(1)) call patch%boundaries%set (bits%xu)
    if (patch%position(2) >= self%task_list%uc(2)) call patch%boundaries%set (bits%yu)
    if (patch%position(3) >= self%task_list%uc(3)) call patch%boundaries%set (bits%zu)
    if (patch%boundaries%status /= 0) print '(2i6,2(2x,a6),o9,2(2x,3f10.0))', patch%id, patch%level, 'BCs: ', &
      formatted(patch%boundaries%status), patch%boundaries%status, patch%size, patch%position
    link => link%next
  end do
  !-----------------------------------------------------------------------------
  call trace_end
contains
  function formatted (status) result (out)
  integer:: status
  character(len=6):: out
  out = '......'
  if (patch%boundaries%is_set(bits%xl))  out(1:1) = 'T'
  if (patch%boundaries%is_set(bits%xu))  out(2:2) = 'T'
  if (patch%boundaries%is_set(bits%yl))  out(3:3) = 'T'
  if (patch%boundaries%is_set(bits%yu))  out(4:4) = 'T'
  if (patch%boundaries%is_set(bits%zl))  out(5:5) = 'T'
  if (patch%boundaries%is_set(bits%zu))  out(6:6) = 'T'
  end function formatted
END SUBROUTINE init_hierarchy

!> =============================================================================
!> Distribute patches in a Rubik's cube arrangement around a box
!> =============================================================================
SUBROUTINE init_level (self, level, position)
  class(rubiks_t):: self
  integer:: level
  real(8):: position(3)
  !.............................................................................
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i, j, k
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from
  ! existing data type values
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_revenge_t%init_level')
  do k=1,self%n_patch
  do j=1,self%n_patch
  do i=1,self%n_patch
    !---------------------------------------------------------------------------
    ! Initialize each patch, via patch_t:init
    !---------------------------------------------------------------------------
    allocate (experiment_t:: task)
    select type (task)
    class is (experiment_t)
    patch => task
    end select
    !---------------------------------------------------------------------------
    ! The patch size should equal the self%size when level=n_levels, since that
    ! is the level of the smallest, central patch.  The size of the largest
    ! patches is m_zoom*(n_levels-1) larger.
    !---------------------------------------------------------------------------
    patch%size = self%size*real(self%m_zoom,kind=8)**(self%n_levels-level)
    patch%position = ([i,j,k]-(1.+self%n_patch)/2.0)*patch%size + position
    !print '(a,i4,1p,2(3x,3g15.5))', 'XX: level, size, pos', level, patch%size, &
    !  patch%position
    !---------------------------------------------------------------------------
    ! The patch%box value is used in periodic wrap-around, so this implements
    ! periodic conditions along the x-direction (= phi-direction in cylindrical
    ! coordinates).
    !---------------------------------------------------------------------------
    patch%box = self%task_list%size
    call patch%init
    patch%level = level
    self%task_list%lc = min(self%task_list%lc, patch%position)
    self%task_list%uc = max(self%task_list%uc, patch%position)
    call patch%set (bits%static)
    if (j==1.or.k==self%n_patch) call patch%set (bits%top_bottom)
    if (io%verbose>0) &
      print'(a,i6,i4,2(2x,3g16.6))', 'rubiks_mod%init_level: id, level, size, pos =', &
        patch%id, level, patch%size, patch%position - position
    !---------------------------------------------------------------------------
    ! Add to a list of links, with pointers from patch to patch
    !---------------------------------------------------------------------------
    call self%task_list%append (task)
  end do
  end do
  end do
  call trace_end
END SUBROUTINE init_level

END MODULE rubiks_mod
