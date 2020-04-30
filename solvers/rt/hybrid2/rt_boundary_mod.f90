!===============================================================================
!> This module describes the boundary conditions for RT
!> module.
!===============================================================================
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
    !-------------- Custom Parameters ------------------------------------
    real    :: hot_temp, cool_temp
  contains
    procedure :: init
    procedure :: set_ibc
    procedure :: set_obc
  end type
  type(rt_boundary_t), public :: rt_boundary
CONTAINS
!===============================================================================
!> Read namelist
!===============================================================================
SUBROUTINE init (self)
  class(rt_boundary_t):: self
  character(len=16), save :: obc_treatment ='custom'
  character(len=16), save :: ibc_treatment ='ignore'
  real, save :: hot_temp = 1000, cool_temp = 10
  logical, save:: first_time=.true.
  integer, save:: obc, ibc
  namelist /rt_bc_params/ obc_treatment, ibc_treatment
  namelist /custom_bc_params/ hot_temp, cool_temp
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_bc_params)
    if (io%master) write (*, rt_bc_params)
    call scaling%init()
    select CASE(obc_treatment)
      case('custom')
!=======================================================================
!> Modify only this part and appropriate namelist
!> ---------------------------------------------------------------------
        obc = 0
        rewind (io%input)
        read (io%input, custom_bc_params)
        hot_temp = hot_temp/scaling%temp
        cool_temp = cool_temp/scaling%temp
!=======================================================================
      case('reflective')
        obc = 1
      case('periodic')
        obc = 2
      case default  ! custom is the default one
        obc = 0
    end SELECT
    select CASE(ibc_treatment)
      case('custom')
        ibc = 0
      case('zero')
        ibc = 1
      case('reflective')
        ibc = 2
      case default ! ignore is the default one.
        ibc = 3
    end SELECT
  end if
  call io%lock%unset
  self%ibc = ibc
  self%obc = obc
  if (obc==0) then
    self%hot_temp = hot_temp
    self%cool_temp = cool_temp
  end if

  call trace_end
END SUBROUTINE init

SUBROUTINE set_ibc (self)
class (rt_boundary_t) :: self

END SUBROUTINE set_ibc

SUBROUTINE set_obc (self, ppos, ds, i_in)
  class (rt_boundary_t) :: self
  real(kind=8) :: ppos(3), ds(3), c
  real(kind=4), dimension(self%n_wl_bins) :: i_in
!=======================================================================
!> Should be modified only from here on
!=======================================================================
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
