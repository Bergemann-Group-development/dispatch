!===============================================================================
!> Create a nested set of patches, in a hierarchy of Rubik's cubes, where the
!> innermost patch has a given size, and there are n_levels-1 sets of larger
!> Rubik's cubes, each three times larger than the previous, with the previous
!> sets of patches inside the center patch.
!===============================================================================
MODULE rubiks_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE mesh_mod
  USE task_list_mod
  USE experiment_mod
  USE units_mod
  USE scaling_mod
  USE scene_mod
  USE parameters_mod
  implicit none
  private
  type, public:: rubiks_t
    real(8) :: size(3) = [1.0, 1.0, 1.0]     ! the min. size of a patch in the hierarchy
    real(8) :: position(3) = [0.0, 0.0, 0.0] ! centre of the hierarchy
    integer :: n_levels                      ! the number of nesting levels
    integer :: n_repeat                      ! ?
    integer :: repeat_dir
    integer :: n_patch(3)                    ! the number of patches per level per direction
    integer :: m_zoom                        ! the zoom factor applied to succesive levels
    integer :: n_per_patch(3)                ! number of zones per patch
    logical :: rotating = .true.
    type(task_list_t), pointer:: task_list => null()
  contains
    procedure init_hierarchy
    procedure init_level
    procedure init_ring
  end type
  public task_list_t                    ! export, to avoid excessive USE
CONTAINS

!===============================================================================
!> Distribute patches in a Rubik's cube arrangement around a box
!===============================================================================
SUBROUTINE init_hierarchy (self, num, rootlvl, label)
  class(rubiks_t):: self
  character(len=*), optional:: label
  integer, intent(in):: num
  integer, optional:: rootlvl
  !.............................................................................
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer, save:: n_patch(3)=3, m_zoom=3
  integer :: n_levels(num), rootlevel
  integer, save:: n_repeat=1
  integer:: iostat
  integer :: n_per_patch(3) = 32, np
  real(8) :: position(3,num), size(3,num), planet_cylpos(3,num)
  integer:: i, j, k, ir, level, levels
  namelist /rubiks_params/ size, position, n_levels, n_patch, m_zoom, n_repeat, &
                           n_per_patch, planet_cylpos
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from
  ! existing data type values.
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_t%init_hierarchy')

  rootlevel = 0
  if (present(rootlvl)) rootlevel = rootlvl

  n_levels = 5

  rewind (io%input); read (io%input, rubiks_params, iostat=iostat)
  if (io%master) write (*, rubiks_params)

  self%n_patch = n_patch
  self%m_zoom = m_zoom
  self%n_per_patch = n_per_patch

  !-----------------------------------------------------------------------------
  ! Save the outer size of the hierarchy as the task list size, for distance
  ! measurements.  The self%task_list%size is the full size of the objects in
  ! the task list, and the size of the period in the periodic wrapping.  The
  ! size, on the other hand, is the size of the smallest, central patch.
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  self%task_list%size = size(:,1)*self%n_patch*self%m_zoom**(n_levels(1)-1)*n_repeat
  n_repeat = (n_repeat/2)*2 + 1
  scene%n_repeat = n_repeat

  do np=1,num ! loop over sets of hierarchies
    do ir=0,n_repeat-1
      position(1,np) = ((ir-n_repeat/2)*self%task_list%size(1))/n_repeat
      levels = n_levels(np)

      if (io%verbose>0) &
        print'("repeated Rubiks Cubes: ir, levels, position =",i4,i3,1p,3es12.4)', &
        ir, levels, position(:,np)

      ! fool-proofing
      if (levels <= 0) then
        print*,"WARNING: The number of nesting levels is <= 0;"
        print*,"No Rubik's cubes will be created (and the task list will be empty!)."
      end if

      do level=rootlevel+1,rootlevel+levels
        call self%init_level (level, n_levels(np)+rootlevel, position(:,np), size(:,np), planet_cylpos(:,np))
      end do
    end do
  end do

  call trace_end
END SUBROUTINE init_hierarchy

!===============================================================================
!> Distribute patches in a Rubik's cube arrangement around a box
!===============================================================================
SUBROUTINE init_level (self, level, n_levels, position, size, planet_cylpos)
  class(rubiks_t):: self
  integer, intent(in) :: level, n_levels
  real(8):: position(3), offset(3), size(3), planet_cylpos(3)
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: i, j, k
  !.............................................................................
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from
  ! existing data type values
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_t%init_level')

  do k=1,self%n_patch(3)
    do j=1,self%n_patch(2)
      do i=1,self%n_patch(1)
        !---------------------------------------------------------------------------
        ! Initialize each patch, via patch_t:init
        !---------------------------------------------------------------------------
        allocate (experiment_t:: task)
        select type (task)
        class is (experiment_t)
          patch => task
        end select

        !---------------------------------------------------------------------------
        ! The patch size should equal the size when level=n_levels, since that
        ! is the level of the smallest, central patch.  The size of the largest
        ! patches is m_zoom*(n_levels-1) larger.
        !---------------------------------------------------------------------------
        patch%n = self%n_per_patch
        patch%size = size*real(self%m_zoom,kind=8)**(n_levels-level)
        where (patch%n == 1)
          patch%size = size
        end where
        patch%position = 0.0 ! in Cartesian coords.
        offset = ([i,j,k]-(1.+self%n_patch)/2.0)*patch%size + patch%position

        patch%llc_nat = offset + planet_cylpos - 0.5 * patch%size
        where (patch%n == 1)
          patch%llc_nat = 0.0
        end where
        patch%mesh_type = mesh_types%cylindrical
        patch%llc_cart = patch%position + CurrentToCartesian(patch%mesh_type, offset(1), offset(2), offset(3))

        !---------------------------------------------------------------------------
        ! The patch%box value is used in periodic wrap-around, so this implements
        ! periodic conditions along the x-direction (= phi-direction in cylindrical
        ! coordinates).
        !---------------------------------------------------------------------------
        patch%box = self%task_list%size

        if (self%rotating) then
            patch%grid_vel = [0.0d0, 0.0d0, sqrt(scaling%grav * experiment_params%mstar &
                                                / planet_cylpos(2)**3)]
          else
            patch%grid_vel = [0.0d0, 0.0d0, 0.0d0]
        end if
        call patch%init()
        patch%level = level
        call patch%set (bits%static)
        if (k==1.or.k==self%n_patch(3)) call patch%set (bits%top_bottom)
        if (io%verbose>0) &
        print'(a,i6.1,i4.1,2(1x,3es14.6))', 'rubiks_t%init_level: id, level, size, pos =', &
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

SUBROUTINE init_ring (self, rootlvl, label)
  class(rubiks_t):: self
  character(len=*), optional:: label
  integer, optional:: rootlvl
  !.............................................................................
  class(experiment_t), pointer:: patch
  class(task_t), pointer:: task
  integer, save:: n_patch(3)=3, m_zoom=3
  integer :: n_levels=0, rootlevel
  integer, save:: n_repeat=1, repeat_dir=0
  integer:: iostat
  integer :: n_per_patch(3) = 32
  real(8) :: position(3)=0.0, size(3)=0.0, planet_cylpos(3)=0.0
  integer:: i, j, k, ir, level, levels
  namelist /ring_params/ size, position, n_levels, n_patch, m_zoom, n_repeat, &
                           n_per_patch, planet_cylpos, repeat_dir
  !-----------------------------------------------------------------------------
  ! Allow size and position to be set by an overloading data type, and allow
  ! an empty namelist input here. If present, take the initial values from
  ! existing data type values.
  !-----------------------------------------------------------------------------
  call trace_begin ('rubiks_t%init_ring')

  rootlevel = 0
  if (present(rootlvl)) rootlevel = rootlvl

  rewind (io%input); read (io%input, ring_params, iostat=iostat)
  if (io%master) write (*, ring_params)

  if (any(n_patch == 0)) return

  self%n_patch = n_patch
  self%m_zoom = m_zoom
  self%n_per_patch = n_per_patch

  !-----------------------------------------------------------------------------
  ! Save the outer size of the hierarchy as the task list size, for distance
  ! measurements.  The self%task_list%size is the full size of the objects in
  ! the task list, and the size of the period in the periodic wrapping.  The
  ! size, on the other hand, is the size of the smallest, central patch.
  !-----------------------------------------------------------------------------
  allocate (self%task_list)
  self%task_list%size = size(:)*self%n_patch*self%m_zoom**(n_levels-1)*n_repeat
  n_repeat = (n_repeat/2)*2 + 1
  scene%n_repeat = n_repeat
  self%n_repeat = n_repeat
  self%repeat_dir = repeat_dir

    do ir=0,n_repeat-1
      position(1) = ((ir-n_repeat/2)*self%task_list%size(1))/n_repeat
      levels = n_levels

      if (io%verbose>0) &
        print'("repeated Rubiks Cubes: ir, levels, position =",i4,i3,1p,3es12.4)', &
        ir, levels, position(:)

      ! fool-proofing
      if (levels <= 0) then
        print*,"WARNING: The number of nesting levels is <= 0;"
        print*,"No Rubik's cubes will be created (and the task list will be empty!)."
      end if

      do level=rootlevel+1,rootlevel+levels
        call self%init_level (level, n_levels+rootlevel, position(:), size(:), planet_cylpos(:))
      end do
    end do

  call trace_end
END SUBROUTINE init_ring

END MODULE rubiks_mod
