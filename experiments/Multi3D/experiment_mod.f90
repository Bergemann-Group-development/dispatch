!===============================================================================
!> Top level MHD task definition -- can either rely entirely on lower levels,
!> including extras_t%initial_conditions, or can optionally do such things here
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE multi3d_mod
  use units_mod, only: cgs
  use scaling_mod, only: scaling
  USE kinds_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure init
    procedure update
    procedure init_level
    procedure bndries
  end type
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  call self%solver_t%init
  call trace_end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(experiment_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%bndries')
  mem => self%mem(:,:,:,:,self%it,:)
  call self%boundaries%boundary(mem)
  call trace%end()
END SUBROUTINE bndries

!===============================================================================
!> This update only concerns the MHD task; the RT tasks are called directly
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  integer :: it_l, in_l, in_lte_l
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  !call self%solver_t%update

  !self%rt%temp8 = dble(self%mem(:,:,:,self%idx%e,  self%it, 1))
  !self%rt%rho8  = dble(self%mem(:,:,:,self%idx%d,  self%it, 1)) * scaling%d
  !self%rt%ne8   = dble(self%mem(:,:,:,6         ,  self%it, 1)) * scaling%p / (cgs%k_b * self%rt%temp8)
  !self%rt%vx8   = dble(self%mem(:,:,:,self%idx%px, self%it, 1)) / self%rt%rho8 *scaling%d*scaling%u
  !self%rt%vy8   = dble(self%mem(:,:,:,self%idx%py, self%it, 1)) / self%rt%rho8 *scaling%d*scaling%u
  !self%rt%vz8   = dble(self%mem(:,:,:,self%idx%pz, self%it, 1)) / self%rt%rho8 *scaling%d*scaling%u

  call trace_end
END SUBROUTINE update

!===============================================================================
SUBROUTINE init_level (self)
  class(experiment_t):: self
  !-----------------------------------------------------------------------------
  if (any(self%n > 1)) then
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
  end if
END SUBROUTINE init_level


END MODULE
