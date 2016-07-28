!> =============================================================================
!> This module describes the boundary conditions for RT
!> module.
!> =============================================================================
MODULE rt_boundary_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  implicit none
  private
  real(8), parameter:: pi = asin(1.0d0) * 2.0d0, radian = 180./pi

  type, public:: rt_boundary_t
  !-------------- Default Parameters -----------------------------------
  integer :: ibc, obc, n_wl_bins=1
  real    :: lbox(3), ubox(3)
  logical :: first_time=.true.
  !-------------- Custom Parameters ------------------------------------
  real    :: hot_temp, cool_temp
  contains
    procedure :: init
    procedure :: set_ibc
    procedure :: set_obc
  end type
  type(rt_boundary_t), public :: rt_boundary
CONTAINS
!> =============================================================================
!> Read namelist
!> =============================================================================
SUBROUTINE init (self)
  class(rt_boundary_t):: self
  character(len=16)   :: obc_treatment ='custom' 
  character(len=16)   :: ibc_treatment ='ignore' 
  real :: hot_temp = 1000, cool_temp = 10
  namelist /rt_bc_params/ obc_treatment, ibc_treatment
  namelist /custom_bc_params/ hot_temp, cool_temp
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  !$omp critical (input_cr)
  if (self%first_time) then
    self%first_time = .false.
    rewind (io%input)
    read (io%input, rt_bc_params)
    if (io%master) write (*, rt_bc_params)
    call scaling%init()
    select CASE(obc_treatment)
      case('custom')
!> =====================================================================      
!> Modify only this part and appropriate namelist      
!> ---------------------------------------------------------------------
        self%obc = 0
        rewind (io%input)
        read (io%input, custom_bc_params)
        self%hot_temp = hot_temp/scaling%temp
        self%cool_temp = cool_temp/scaling%temp
!> =====================================================================
      case('reflective')
        self%obc = 1
      case('periodic')
        self%obc = 2
      case default  ! custom is the default one
        self%obc = 0      
    end SELECT
    select CASE(ibc_treatment)
      case('custom')
        self%ibc = 0
      case('zero')
        self%ibc = 1
      case('reflective')
        self%ibc = 2
      case default ! ignore is the default one.
        self%ibc = 3
    end SELECT
  end if
  !$omp end critical (input_cr)

  call trace_end
END SUBROUTINE init

SUBROUTINE set_ibc (self)
class (rt_boundary_t) :: self

END SUBROUTINE set_ibc

SUBROUTINE set_obc (self, ppos, ds, i_in)
  class (rt_boundary_t) :: self
  real(kind=8) :: ppos(3), ds(3), c
  real(kind=8), dimension(self%n_wl_bins) :: i_in
!> =====================================================================
!> Should be modified only from here on  
!> =====================================================================
  c = cgs%stefan/cgs%pi*scaling%temp**4/scaling%f
  if (ppos(1) <= -0.5*self%lbox(1)) then
    i_in(:) = c*self%hot_temp**4
  else if (ppos(1) >= 0.5*self%ubox(1)) then
    i_in(:) = c*self%cool_temp**4
  end if
  
  if (ppos(2) <= -0.5*self%lbox(2)) then
    i_in(:) = c*self%cool_temp**4
  else if (ppos(2) >= 0.5*self%ubox(2)) then
    i_in(:) = c*self%cool_temp**4
  end if
  
  if (ppos(3) <= -0.5*self%lbox(3)) then
    i_in(:) = c*self%cool_temp**4
  else if (ppos(3) >= 0.5*self%ubox(3)) then
    i_in(:) = c*self%cool_temp**4
  end if
END SUBROUTINE set_obc

END MODULE
