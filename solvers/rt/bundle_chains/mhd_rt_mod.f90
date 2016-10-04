!> A wrapper module for different mhd solvers. It seamlesly incorporates
!> radiative transfer into the solver/
MODULE mhd_rt_mod
  USE mhd_mod
  USE trace_mod
  USE microphysics_mod
  USE rt_data_mod
  USE rt_internal_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: mhd_rt_t
    type(rt_data_t), pointer:: rt_data => null()
    integer:: BC_direction=0
  contains
    procedure :: init_rt
    procedure :: update => update_mhd_rt
  end type
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> This is currently a place holder for setting up the RT for one patch
!> =====================================================================
SUBROUTINE init_rt (self)
  class (mhd_rt_t) :: self
  integer :: nt = 5
  ! --------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%init_mhd_rt')

  call trace_end
END SUBROUTINE init_rt

!> =====================================================================
!> This is currently a place holder for solving the RT for one patch --
!> it actually does not do the RT part yet (FIXME)
!> =====================================================================
SUBROUTINE update_mhd_rt(self)
  class (mhd_rt_t) :: self
  real, dimension(:,:,:), pointer:: q
  ! --------------------------------------------------------------------
  call trace_begin ('mhd_rt_t_%update_mhd_rt')
  ! --- Update MHD, whatever it means ---
  call self%rt_data%update
  call self%mhd_t%update
  ! --- Update dsdt by addind Q ---
  q => self%rt_data%q(:,:,:,self%it)  ! FIXME: should be interpolation in time
  !call self%apply_heating (q)
  call trace_end
END SUBROUTINE update_mhd_rt

END MODULE mhd_rt_mod
