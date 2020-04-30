!===============================================================================
!> Experiment-specific set up
!===============================================================================
MODULE experiment_setup_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE task_list_mod
  USE experiment_mod
  USE scaling_mod
  USE math_mod
  USE mpi_mod
  USE bits_mod
  USE mesh_mod, only: mesh_types, CurrentToCartesian
  USE units_mod
  implicit none
  private

  type, public:: experiment_setup_t
    type(task_list_t), pointer:: task_list
    character(len=64) :: experiment_name = 'jet launching disk'
    integer, dimension(3):: n_patches    ! number of patches in each direction
    integer, dimension(3):: n_per_patch  ! number of cells per patch in each direction
    real(8), dimension(3):: extent       ! total extent of the domain in each direction
    real(8):: r_offset                   ! radial coordinate offset; should be non-zero when
                                         ! using spherical coords.
    real(8):: theta_offset               ! theta coordinate offset given as a factor of pi;
                                         ! can only be used with spherical coords.
    character(len=16):: coord_type       ! 'spherical' or 'cylindrical'
  contains
    procedure:: init
  end type
  type(experiment_setup_t), public:: experiment_setup

CONTAINS

!===============================================================================
!> Initialise the experiment and initial set of patches
!===============================================================================
SUBROUTINE init (self)
  class(experiment_setup_t):: self
  logical, save:: first_time = .true.
  integer, dimension(3), save:: n_patches = [1,1,1]
  integer, dimension(3), save:: n_per_patch = [32,32,1]
  real(8), dimension(3), save:: extent = [9.98,1.0,2.0]
  real(8), save:: r_offset = 0.0
  real(8), save:: theta_offset = 0.0
  character(len=16), save:: coord_type = 'spherical'
  namelist /experiment_params/ n_patches, n_per_patch, extent, r_offset, &
                               coord_type

  integer:: i1, i2, i3, ip
  integer:: iostat
  real(8), dimension(3):: pextent, curr_offset
  class(task_list_t), pointer:: task_list
  class(task_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_setup_t%init')

  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    if (iostat/=0) call io%namelist_warning ('experiment_params')
    if (io%master) write (*, experiment_params)
  end if
  call io%lock%unset

  if (coord_type == 'spherical') then
    if (extent(2) == 1.0) then
      ! even in "no-no-mans-land", we want symmetry axes to be cospatial with
      ! cell faces, not cell centres (and thus the corrections below)!
      extent(2) = extent(2) * (1.0d0 - 0.9999d0 / (n_per_patch(2)-1) / n_patches(2))
      theta_offset = 0.5d0 * extent(2) / (n_per_patch(2)-1) / n_patches(2)
      print*,'A:',extent(2),theta_offset,extent(2)*math%pi,theta_offset*math%pi
    end if
    ! the theta and phi extents are given in units of pi
    extent(2) = extent(2) * math%pi
    extent(3) = extent(3) * math%pi
    ! the radial extent is given in units of au
    extent(1) = extent(1) * (1.0d0 - 0.5d0 / (n_per_patch(1)-1) / n_patches(1))
    r_offset = r_offset + 0.5d0 * extent(1) / (n_per_patch(1)-1) / n_patches(1)
    extent(1) = extent(1) * cgs%au
  else if (coord_type == 'cylindrical') then
    ! the phi extent is given in units of pi
    extent(3) = extent(3) * math%pi
    ! the radial and vertical extents are given in units of au
    extent(1) = extent(1) * cgs%au
    extent(2) = extent(2) * cgs%au
  end if

  ! fool-proofing
  if (trim(coord_type) /= 'spherical' .and. trim(coord_type) /= 'cylindrical') then
    print*,trim(coord_type)
    call mpi%abort('Unsupported coordinate type for this experiment!')
  end if
  if (trim(coord_type) == 'spherical' .and. r_offset <= 0.0) then
    print*,r_offset
    call mpi%abort('Spherical coords. need a non-zero radial offset!')
  end if
  if (trim(coord_type) == 'cylindrical' .and. r_offset /= 0.0) then
    print*,r_offset
    print*,'WARNING! Radial offset is not zero and cylindrical coords. are in use!'
  end if
  if (trim(coord_type) == 'cylindrical' .and. theta_offset /= 0.0) then
    print*,theta_offset
    print*,'WARNING! theta_offset is ignored in cylindrical coords.!'
  end if
  if (any(n_patches <= 0)) then
    print*,n_patches
    call mpi%abort('n_patches <= 0 in one or more directions!')
  end if
  if (any(n_per_patch <= 0)) then
    print*,n_per_patch
    call mpi%abort('n_per_patch <= 0 in one or more directions!')
  end if

  self%n_patches = n_patches
  self%n_per_patch = n_per_patch
  self%extent = extent
  self%r_offset = r_offset * cgs%au ! convert to cgs
  self%theta_offset = theta_offset * math%pi ! given as a factor of pi
  self%coord_type = trim(coord_type)
  pextent = extent / n_patches ! extent of each patch; domain is assumed to be
                               ! initially evenly divided

  allocate (self%task_list)
  call self%task_list%init ('jet-launching-disk')

  curr_offset(3) = 0.0
  if (n_per_patch(3) == 1) curr_offset(3) = -math%pi
  if (trim(coord_type) == 'spherical') then
    curr_offset(1) = self%r_offset
  else if (trim(coord_type) == 'cylindrical') then
    curr_offset(1) = 0.0
  end if
  if (trim(coord_type) == 'spherical') then
    curr_offset(2) = self%theta_offset
  else if (trim(coord_type) == 'cylindrical') then
    curr_offset(2) = self%r_offset
  end if

  ip = 0
  do i3=1,n_patches(3)
    do i2=1,n_patches(2)
      do i1=1,n_patches(1)
        allocate (experiment_t:: patch)
        select type (patch)
        type is (experiment_t)
          patch%rank = mpi%rank
          if (patch%rank == mpi%rank) then
            ip = ip + 1
            patch%ip = ip
          end if
          patch%size = pextent
          patch%position = curr_offset + ([i1,i2,i3] - 0.5) * pextent
          patch%n = n_per_patch
          if (trim(coord_type) == 'spherical') then
            patch%mesh_type = mesh_types%spherical
          else if (trim(coord_type) == 'cylindrical') then
            patch%mesh_type = mesh_types%cylindrical
          end if
          patch%box = extent
          patch%llc_nat = curr_offset
          patch%llc_cart = patch%position &
                         + CurrentToCartesian(patch%mesh_type, curr_offset(1), curr_offset(2), curr_offset(3))
          ! remainder of patch initialization
          call patch%init()
          if (trim(coord_type) == 'spherical') then
            patch%mesh(:)%origin = [self%r_offset, self%theta_offset, 0.0_8]
          else if (trim(coord_type) == 'cylindrical') then
            patch%mesh(:)%origin = [0.0_8, self%r_offset, 0.0_8]
          end if
          patch%origin = patch%mesh%origin
        end select
        ! only "periodic" in the collapsed dimension (x3)
        patch%periodic(1) = .false.
        patch%periodic(2) = .false.
        call patch%clear(bits%root_grid)
        ! The `list_t%append` procedure only accepts objects of type `task_t`.
        call self%task_list%append(patch)
      end do
    end do
  end do
  self%task_list%n_tasks = product(n_patches)

  !-----------------------------------------------------------------------------
  ! Done. Perform final list initialisations.
  !-----------------------------------------------------------------------------
  call self%task_list%reset_status
  call self%task_list%count_status
  call self%task_list%init_bdries
  call self%task_list%init_all_nbors
  call self%task_list%print_tasks

  call trace%end()
END SUBROUTINE init

END MODULE experiment_setup_mod
