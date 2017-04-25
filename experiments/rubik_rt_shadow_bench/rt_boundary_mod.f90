!> =============================================================================
!> This module describes the boundary conditions for RT
!> module.
!> =============================================================================
MODULE rt_boundary_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE eos_mod
  implicit none
  private
  real(8), parameter:: pi = asin(1.0d0) * 2.0d0, radian = 180./pi

  type, public:: rt_boundary_t
    !-------------- Default Parameters -----------------------------------
    integer :: ibc, obc, n_wl_bins=1
    real(8)    :: lbox(3), ubox(3)
    real(8)    :: ibnd(3), ubnd(3)
    real(8)    :: c
    real:: gamma = 1.4                  ! in case eos%init is never called
    !-------------- Custom Parameters ------------------------------------
    real, dimension(:), pointer    :: hot_temp, cool_temp
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
SUBROUTINE init (self, lbox, ubox)
  class(rt_boundary_t):: self
  character(len=16), save :: obc_treatment ='custom'
  character(len=16), save :: ibc_treatment ='ignore'
  real, save :: hot_temp = 1000, cool_temp = 10
  logical, save:: first_time=.true.
  integer, save:: obc, ibc
  real(kind=8):: lbox(3), ubox(3)
  real, dimension(:,:,:), pointer :: lnd, ee
  real, dimension(:,:,:,:), pointer :: src
  namelist /rt_bc_params/ obc_treatment, ibc_treatment
  namelist /custom_bc_params/ hot_temp, cool_temp
  !----------------------------------------------------------------------------
  call trace_begin('initial_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_bc_params)
    if (io%master) write (*, rt_bc_params)
    call scaling%init()
    select CASE(obc_treatment)
      case('custom')
!> =====================================================================
!> Modify only this part and appropriate namelist
!> ---------------------------------------------------------------------
        obc = 0
        rewind (io%input)
        read (io%input, custom_bc_params)
        allocate(self%hot_temp(eos%n_lambda))
        allocate(self%cool_temp(eos%n_lambda))
     self%c = cgs%stefan/cgs%pi*scaling%temp**4/scaling%f
      self%hot_temp = self%c*(hot_temp/scaling%temp)**4
    self%cool_temp = self%c*(cool_temp/scaling%temp)**4
    self%cool_temp = 2.31530362E-06
!> =====================================================================
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
  !allocate(lnd(1,1,2))
  !allocate(ee(1,1,2))
  !allocate(src(1,1,2, eos%n_lambda))
  !lnd = 0.
  print *, 'before',self%cool_temp, self%hot_temp
  !ee(:,:,1) = self%cool_temp(1)/(self%gamma - 1.)
  !ee(:,:,2) = self%hot_temp(1)/(self%gamma - 1.)
  
  !call eos%lookup (shape(ee), ee=ee, lnd=lnd, &       ! call the table lookup, ...
  !                  src=src)  
  !self%cool_temp = self%cool_temp*src(1,1,1,:)
  !self%hot_temp  = self%hot_temp*src(1,1,2,:)
  !deallocate(lnd,ee,src)
  end if
  !$omp end critical (input_cr)
  self%ibc = ibc
  self%obc = obc

  
print *, 'cool', self%cool_temp
print *, 'hot', self%hot_temp
  call trace_end
END SUBROUTINE init

SUBROUTINE set_ibc (self, ppos, i_in)
class (rt_boundary_t) :: self
real(kind=8) :: ppos(3)
real(kind=4), dimension(self%n_wl_bins) :: i_in


END SUBROUTINE set_ibc

SUBROUTINE set_obc (self, ppos, i_in)
  class (rt_boundary_t) :: self
  real(kind=8) :: ppos(3)
  real, dimension(:) :: i_in
!> =====================================================================
!> Should be modified only from here on
!> =====================================================================

  if (ppos(1) <= self%lbox(1)) then
    i_in = self%hot_temp
  else if (ppos(1) >= self%ubox(1)) then
    i_in = self%cool_temp
  end if
  
  if (ppos(2) <= self%lbox(2)) then
    i_in = self%cool_temp
  else if (ppos(2) >= self%ubox(2)) then
    i_in = self%cool_temp
  end if

  if (ppos(3) <= self%lbox(3)) then
    i_in = self%cool_temp
  else if (ppos(3) >= self%ubox(3)) then
    i_in = self%cool_temp
  end if
  
 

END SUBROUTINE set_obc

END MODULE
