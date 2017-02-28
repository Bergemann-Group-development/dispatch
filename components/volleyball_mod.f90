!> =============================================================================
!> Set up potentially nested set of spherical patches centred on a sphere.
!> The minimum number of patches is six, which is designed to cover 4*pi with
!> uniform-mesh patches and a minimum of overlap.
!> This set up is inspired by the patches of a volleyball, but the patches of a
!> volleyball cannot be described by curvilinear spherical coordinates.
!>
!> nsplit: Starting with the six base patches, split each patch into nsplit**3
!>         additional patches.
!>
!> =============================================================================
MODULE volleyball_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE mesh_mod
  USE task_mod
  USE task_list_mod
  USE experiment_mod
  implicit none
  private

  type, public:: volleyball_t
    integer:: nsplit
    integer, dimension(3):: n_per_patch
    real(8):: rmin, rmax
    type(task_list_t), pointer:: task_list
  contains
    procedure:: init
    procedure:: make_patch
  end type

  real(8), parameter :: pi = asin(1.0) * 2.0

CONTAINS
!> =============================================================================
SUBROUTINE init (self, label)
  class(volleyball_t):: self
  character(len=*), optional:: label
  integer:: nsplit=1, ibase, ispl
  integer, dimension(3):: n_per_patch=[32,32,32]
  real(8):: rmin=1.0, rmax=1.1, bt, bp, cr
  real(8), dimension(3):: extent(3)
  real(8), dimension(6):: base_theta, base_phi
  namelist /volleyball_params/ rmin, rmax, nsplit, n_per_patch
  !-----------------------------------------------------------------------------
  call trace_begin ('volleyball_t%init')

  rewind (io%input); read (io%input, volleyball_params)
  if (io%master) write (*, volleyball_params)
  self%rmin = rmin
  self%rmax = rmax
  self%nsplit = nsplit
  self%n_per_patch = n_per_patch

  base_theta = [0.5d0*pi, 0.5d0*pi, 0.5d0*pi, 0.5d0*pi, 0.0d0,  pi]
  base_phi   = [   0.0d0, 0.5d0*pi,       pi, 1.5d0*pi, 0.0d0, 0.0d0]
  extent     = [rmax-rmin, 0.5d0*pi, 0.5d0*pi] / real(nsplit,kind=8)

  allocate (self%task_list)
  self%task_list%size = [rmax-rmin,pi,2.0d0*pi]

  do ibase=1,6
    bt = base_theta(ibase)
    bp = base_phi(ibase)
    do ispl=1,nsplit
      cr = rmin + extent(0) * real(ispl - 1,kind=8)
      call self%make_patch (ibase, cr, bt, bp, extent, ispl, n_per_patch)
    end do
  end do

  call self%task_list%init (label)
  call self%task_list%init_all_nbors

  call trace_end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE make_patch (self, ibase, cr, ctheta, cphi, extent, ispl, npp)
  class(volleyball_t):: self
  integer, intent(in):: ispl, ibase
  integer, dimension(3), intent(in):: npp
  real(8), intent(in):: cr, ctheta, cphi
  real(8), dimension(3), intent(in):: extent
  class(task_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace_begin ('volleyball_t::make_patch')

  allocate (experiment_t:: patch)
  select type (patch)
  class is (experiment_t)
    patch%size = extent
    patch%position = 0.0d0
    patch%n = npp
    patch%mesh_type = mesh_types%spherical
    patch%llc_nat = [cr, ctheta, cphi] - 0.5d0 * extent
    patch%llc_cart = patch%position &
                   + CurrentToCartesian(patch%mesh_type, patch%llc_nat(1), patch%llc_nat(2), patch%llc_nat(3))
    patch%box = self%task_list%size
    patch%grid_vel = 0.0d0
    call patch%init()
  end select
  call patch%set(bits%static)
  call patch%clear(bits%root_grid)
  call self%task_list%append(patch)

  call trace_end
END SUBROUTINE make_patch

END MODULE volleyball_mod
