!> =============================================================================
!> Set up the de Val-Borro protoplanetary disk benchmarks (2006, MNRAS, 370, 529).
!> =============================================================================
MODULE components_mod
  USE io_mod
  USE trace_mod
  USE scaling_mod
  USE parameters_mod
  USE task_mod
  USE task_list_mod
  USE mhd_mod
  USE experiment_mod
  USE mesh_mod
  USE bits_mod
  implicit none
  private

  ! the computational domain and cylindrical annulus
  type, public :: annulus_t
    integer, dimension(3) :: n_per_patch, n_patch
    real(8), dimension(3) :: extent
    logical:: rotating = .true.
    type(task_list_t), pointer :: task_list => null()
  contains
    procedure :: init => init_annulus
  end type

  real(8), parameter :: pi = asin(1.0) * 2.0

CONTAINS

!> =============================================================================
!> Initialize the cylindrical domain ("annulus").
!> =============================================================================
SUBROUTINE init_annulus (self, params)
  class(annulus_t) :: self
  type(experiment_parameters_t) :: params
  class(task_t), pointer :: patch
  integer :: ix, iy, iz
  integer, dimension(3) :: n_patch, n_per_patch
  real(8), dimension(3) :: extent, llc, urc
  real(8) :: zextent, rextent, phiextent, phioffset, zoffset, roffset, &
             curr_zoff, curr_roff, curr_phioff
  logical, save :: first_time = .true., rotating

  namelist /annulus_params/ n_per_patch, n_patch, extent, roffset, rotating
  !-----------------------------------------------------------------------------
  call trace_begin ('annulus_t%init')

  if (first_time) then
    first_time = .false.
    ! default values
    extent = [1.0, 2.1, 2.0]
    roffset = 0.4
    ! axisymmetry in the z-direction is assumed
    ! (z, r , phi)
    n_per_patch = [1, 128, 384]
    n_patch = [1, 1, 1]
    rotating = .true.
    rewind (io%input); read (io%input, annulus_params)
    if (io%master) write (*, annulus_params)
    ! the phi extent is given in units of pi
    extent(3) = extent(3) * pi
  end if

  self%n_patch = n_patch
  self%extent = extent
  self%n_per_patch = n_per_patch
  experiment_params%roffset = roffset

  allocate(self%task_list)
  call self%task_list%init ('annulus task list')

  zextent   = extent(1) / n_patch(1)
  rextent   = extent(2) / n_patch(2)
  phiextent = extent(3) / n_patch(3)
  phioffset = -pi
  if (n_per_patch(1) == 1) then
    zoffset = 0.0
  else
    zoffset = -0.5 * zextent
  end if

  curr_phioff = phioffset
  do iz=1,self%n_patch(3)
    curr_roff = roffset
    do iy=1,self%n_patch(2)
      curr_zoff = zoffset
      do ix=1,self%n_patch(1)
        allocate(experiment_t :: patch)
        select type (patch)
        class is (mhd_t)
          patch%size = [zextent, rextent, phiextent]
          patch%position = 0.0 ! still in Cartesian coords.; absolute pos.
          patch%name => self%task_list%name
          patch%n = n_per_patch
          ! because absolute positions are sometimes required in curvilinear
          ! coords., set the cylindrical offset (lower left corner).
          patch%llc_nat = [curr_zoff, curr_roff, curr_phioff]
          patch%llc_cart = patch%position &
                         + CurrentToCartesian(patch%mesh_type, curr_zoff, curr_roff, curr_phioff)
          patch%mhd = .false.
          patch%periodic(2) = .false.
          patch%mesh_type = mesh_types%cylindrical
          patch%box = extent
          ! set angular (phi) velocity of the annulus (in rad/s), if desired.
          if (rotating) then
            patch%grid_vel = [0.0d0, 0.0d0, &
                              sqrt(scaling%grav * params%mstar / params%aplanet**3)]
          else
            patch%grid_vel = [0.0d0, 0.0d0, 0.0d0]
          end if

          call patch%init()
          !if (patch%kind == 'zeus_mhd_patch') then
          !  patch%staggered = .true.
          !  if (phiextent .ge. 2.0*pi) then
          !    patch%btypes%ikb = 5 ! periodic
          !    patch%btypes%okb = 5
          !  else
          !    patch%btypes%ikb = 6 ! self-computing (but see `boundary%conditions`)
          !    patch%btypes%okb = 6
          !  end if
          !  patch%btypes%ijb = 6
          !  if (curr_roff .le. roffset) then
          !    patch%btypes%ijb = 2
          !  end if
          !  patch%btypes%ojb = 6
          !  if (curr_roff+rextent .ge. extent(2)) then
          !    patch%btypes%ojb = 2
          !  end if
          !end if
        end select
        call patch%set(bits%static) ! for the time being (should be okay)
        call patch%clear(bits%root_grid)
        call self%task_list%append(patch)
        ! update offset for next patch
        if (n_patch(1) > 1) curr_zoff = curr_zoff + zextent
      end do
      if (n_patch(2) > 1) curr_roff = curr_roff + rextent
    end do
    if (n_patch(3) > 1) curr_phioff = curr_phioff + phiextent
  end do

  call trace_end
END SUBROUTINE init_annulus

END MODULE components_mod
