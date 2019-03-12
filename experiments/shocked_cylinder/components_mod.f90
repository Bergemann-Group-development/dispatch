!===============================================================================
!===============================================================================
MODULE components_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE task_list_mod
  USE mhd_mod
  USE experiment_mod
  USE mesh_mod
  USE bits_mod
  USE immersed_boundary_mod
  implicit none
  private
  
  ! the computational domain and Cartesian box where everything is contained
  type, public :: box_t
    integer, dimension(3) :: n_per_patch, n_patch
    real(8), dimension(3) :: extent
    type(task_list_t), pointer :: task_list
  contains
    procedure :: init => init_box
  end type

  ! an immersed boundary that acts like a solid cylinder
  type, public :: solid_cylinder_t
    real(4) :: radius, height, mass
    real(4), dimension(3) :: position
    type(task_list_t), pointer :: task_list
  contains
    procedure :: init => init_solid_cylinder
  end type

  ! cylindrical annulus surrounding the solid cylinder
  type, public :: envelope_annulus_t
    integer, dimension(3) :: n_patch, n_per_patch
    real(4) :: iradius, oradius, height
    type(task_list_t), pointer :: task_list
  contains
    procedure :: init => init_envelope
  end type

  real(8), parameter :: pi = 3.14159265358979323846
  real(8), save :: wholebox(3)

CONTAINS

!===============================================================================
!> Initialize the Cartesian domain ("box").
!===============================================================================
SUBROUTINE init_box (self)
  class(box_t) :: self
  class(task_t), pointer :: patch
  integer :: ix, iy, iz
  integer, dimension(3) :: n_patch, n_per_patch
  real(8), dimension(3) :: extent, centre, llc, urc
  logical, save :: first_time = .true.
  namelist /cartesian_params/ n_per_patch, n_patch, extent, centre
  !-----------------------------------------------------------------------------
  call trace_begin ('box%init')
  
  if (first_time) then
    first_time = .false.
    ! default values
    extent = [4.0, 4.0, 1.0]
    centre = [0.0, 0.0, 0.0]
    n_per_patch = [32, 32, 1]
    n_patch = [1, 1, 1]
    rewind (io%input); read (io%input, cartesian_params)
    wholebox = extent
  end if

  self%n_patch = n_patch
  self%extent = extent
  self%n_per_patch = n_per_patch
  where (self%n_per_patch /= 1)
    llc = centre - extent * 0.5
    urc = centre + extent * 0.5
  elsewhere
    llc = centre
    urc = centre
  end where

  allocate(self%task_list)
  call self%task_list%list_t%init ('box task list')
  do iz=1,self%n_patch(3)
    do iy=1,self%n_patch(2)
      do ix=1,self%n_patch(1)
        allocate(experiment_t :: patch)
        select type (patch)
        class is (mhd_t)
          patch%box = extent
          patch%size = self%extent/self%n_patch
          where (self%n_per_patch /= 1)
            patch%position = ([ix,iy,iz] - 0.5) * patch%size + llc
          elsewhere
            patch%position = llc
          end where
          patch%name => self%task_list%name
          patch%n = n_per_patch
          patch%llc_cart = llc
          patch%llc_nat  = llc
          patch%mhd = .false.
          patch%staggered = .true.
          patch%mesh_type = mesh_types%Cartesian
          ! search `gp boundary` in the ZEUScore source for a description of
          ! boundary types and their indices.
          patch%btypes%iib = 10 ! inflow
          patch%btypes%oib = 9  ! outflow
          patch%btypes%ijb = 1  ! reflection
          patch%btypes%ojb = 1
          patch%btypes%ikb = 5  ! periodic
          patch%btypes%okb = 5
        end select
        call patch%init()
        call patch%set(bits%static) ! for the time being
        !!!call patch%clear(bits%periodic)
        ! the Cartesian domain is assumed to have the lowest quality of the
        ! three components in this experiment.
        patch%quality = 1.0
        patch%level = 6 ! FIXME
        call self%task_list%append(patch)
      end do
    end do
  end do

  call trace_end
END SUBROUTINE init_box

!===============================================================================
!> Initialize a solid cylinder with a given radius and position
!===============================================================================
SUBROUTINE init_solid_cylinder (self)
  class(solid_cylinder_t) :: self
  class(task_t), pointer :: task
  !-----------------------------------------------------------------------------
  ! Namelist parameters for the component
  !-----------------------------------------------------------------------------
  real(8) :: radius, height, mass
  real, dimension(3) :: position
  namelist /cylinder_params/ position, radius, height, mass
  !-----------------------------------------------------------------------------
  call trace_begin ('solid_cylinder%init')

  ! default values
  radius = 0.5
  height = 1.0
  position = [-0.5,0.0,0.0] ! in Cartesian coords.
  mass = 0.5

  rewind (io%input); read (io%input, cylinder_params)
  self%radius = radius
  self%height = height
  self%position = position
  self%mass = mass

  allocate(self%task_list)
  call self%task_list%list_t%init ('cylinder task list')

  allocate(immersed_cylinder :: task)
  select type (task)
  type is (immersed_cylinder)
    task%size = [height, radius, 2.0*pi]
    task%gsize = task%size
    task%position = position
    task%name => self%task_list%name
    task%radius = self%radius
    task%height = self%height
    task%mass = self%mass
    task%cyl_density = self%mass / (pi * task%radius**2 * task%height)
    call task%init()
    task%mesh_type = mesh_types%cylindrical
    task%box = wholebox
  end select
  call task%set(bits%frozen)
  call task%set(bits%static) ! for the time being
  ! the immersed boundary is assumed to have the highest quality of the components
  ! in this experiment so that it will be preferentially used in `download_mod`;
  ! using `level` will not accomplish this unless you short circuit its normal
  ! assignment.
  task%level = 6
  task%quality = 3.0
  call self%task_list%append(task)
  ! because the solid cylinder is set to "frozen", we need to decrement the number
  ! of patches in the task list.
  self%task_list%n = self%task_list%n - 1

  call trace_end
END SUBROUTINE init_solid_cylinder

!===============================================================================
!> Initialize a cylindrical annulus which envelops a solid cylinder.
!===============================================================================
SUBROUTINE init_envelope (self, cylinder)
  class(envelope_annulus_t) :: self
  class(solid_cylinder_t), intent(in) :: cylinder
  class(task_t), pointer :: patch
  integer :: ix, iy, iz
  real(8) :: zextent, rextent, phiextent, phioffset, zoffset
  !-----------------------------------------------------------------------------
  ! Namelist parameters for the component
  !-----------------------------------------------------------------------------
  integer, dimension(3) :: n_patch, n_per_patch
  real(8) :: radius
  namelist /envelope_params/ n_patch, n_per_patch, radius
  !-----------------------------------------------------------------------------
  call trace_begin ('envelope%init')

  ! default values
  n_patch = 1
  n_per_patch = [1, 50, 180]
  radius = cylinder%radius + 0.6
  
  rewind (io%input); read (io%input, envelope_params)
  self%oradius = radius
  self%n_patch = n_patch
  self%n_per_patch = n_per_patch
  ! inherit properties of the cylinder lying at the centre of the envelope
  self%iradius = cylinder%radius
  self%height  = cylinder%height

  allocate(self%task_list)
  call self%task_list%list_t%init ('envelope task list')

  zextent = self%height / n_patch(1)
  rextent = (self%oradius - self%iradius) / n_patch(2)
  phiextent = 2.0*pi / n_patch(3)
  phioffset = 0.0
  if (n_per_patch(1) == 1) then
    zoffset = 0.0
  else
    zoffset = -0.5 * zextent
  end if

  do iz=1,self%n_patch(3)
    do iy=1,self%n_patch(2)
      do ix=1,self%n_patch(1)
        allocate(experiment_t :: patch)
        select type (patch)
        class is (mhd_t)
          patch%size = [zextent, rextent, phiextent]
          patch%position = cylinder%position ! still in Cartesian coords.; absolute pos.
          patch%name => self%task_list%name
          patch%n = n_per_patch
          ! because absolute positions are sometimes required in curvilinear
          ! coords., set the cylindrical offset (lower left corner).
          patch%llc_nat  = [zoffset, real(self%iradius,kind=8), phioffset]
          patch%llc_cart = patch%position + CylindricalToCartesian(zoffset,real(self%iradius,kind=8), phioffset)
          patch%mhd = .false.
          patch%staggered = .true.
          patch%mesh_type = mesh_types%cylindrical
          if (phiextent .ge. 2.0*pi) then
            patch%btypes%ikb = 5 ! periodic
            patch%btypes%okb = 5
          end if
          patch%box = wholebox
        end select
        call patch%init()
        call patch%set(bits%static) ! for the time being
        patch%level = 6
        ! intermediate quality; better than `box`, but worse than `cylinder`.
        patch%quality = 2.0
        call self%task_list%append(patch)
      end do
    end do
    ! update phioffset for next patch
    phioffset = phioffset + phiextent
  end do

  call trace_end
END SUBROUTINE init_envelope

END MODULE components_mod
