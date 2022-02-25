!===============================================================================
!> Initial conditions for Whole Sun convection benchmark
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE index_mod
  USE math_mod
  USE eos_mod
  USE gpatch_mod
  USE random_mod
  USE math_mod
  USE benchmark_mod
  USE patch_mod
  implicit none
  private
  type, public :: initial_conditions_t
    real:: theta        = 10.0
    real:: mpoly        = 1.0
    real:: z0           = 0.2
    real:: amplitude(2) = [1e-3, 0.0]
    integer:: seed      = 11
    logical:: mhd       = .false.
    character(len=64):: solver
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
  end type
  type(random_t):: random
CONTAINS

!===============================================================================
!> Pre-initialize the initial conditions.
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_conditions_t%pre_init')
  self%patch => patch
  self%solver = patch%kind
  call random%init (self%seed)
  call trace%end()
END subroutine pre_init

!===============================================================================
!> Initial conditions, callable also from refine_mod.f90
!===============================================================================
SUBROUTINE init (self)
  class(initial_conditions_t) :: self
  !.............................................................................
  class(patch_t), pointer:: patch
  real(8), dimension(:), pointer:: x, y, z
  real, dimension(:,:,:), pointer:: d, e, u, v, w, s
  !-----------------------------------------------------------------------------
  call trace_begin ('initial_t%initial_conditions')
  patch => self%patch
  d => patch%mem(:,:,:,patch%idx%d ,patch%it,1)
  u => patch%mem(:,:,:,patch%idx%px,patch%it,1)
  v => patch%mem(:,:,:,patch%idx%py,patch%it,1)
  w => patch%mem(:,:,:,patch%idx%pz,patch%it,1)
  allocate (x(patch%gn(1)),y(patch%gn(2)),z(patch%gn(3)))
  x = patch%mesh(1)%r + patch%mesh(1)%p
  y = patch%mesh(2)%r + patch%mesh(2)%p
  z = patch%mesh(3)%r + patch%mesh(3)%p
  if (patch%idx%e > 0) then
    e => patch%mem(:,:,:,patch%idx%e,patch%it,1)
    call benchmark%initial (z, d, e, u, v, w)
    call benchmark%perturbations (x, y, z, d)
    e = e + 0.5*(u**2+v**2+w**2)*d
  else if (patch%idx%s > 0) then
    allocate (e(patch%gn(1),patch%gn(2),patch%gn(3)))
    call benchmark%initial (z, d, e, u, v, w)
    call benchmark%perturbations (x, y, z, d)
    s => patch%mem(:,:,:,patch%idx%s ,patch%it,1)
    s = d*log(e*(eos%gamma-1d0)/d**eos%gamma)/(eos%gamma-1d0)
    deallocate (e)
  end if
  u = d*u
  v = d*v
  w = d*w
  deallocate (x,y,z)
  !-----------------------------------------------------------------------------
  call trace_end
END SUBROUTINE init

END MODULE initial_conditions_mod
