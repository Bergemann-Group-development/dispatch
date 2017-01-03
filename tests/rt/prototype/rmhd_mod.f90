!> =====================================================================
!> This module combines the MHD and RT solver solutions, to produce
!> time derivatives that include the net radiative cooling
!>
!> Radiative transfer module structure:
!>
!> experiment_mod                       ! generic time advance modul
!>   rmhd_mod                           ! wrapper for MHD and RT
!>     mhd_mod                          ! MHD solver
!>       eos_mod                        ! computing gas pressure
!>     rt_mod                           ! RT organizer
!>       eos_mod                        ! computing opacity and source fn
!>       rt_bundle_mod                  ! handles bundles of rays
!>         rt_solver_mod                ! solver the RT equation
!> =====================================================================
MODULE rmhd_mod
  USE io_mod
  USE trace_mod
!  USE eos_mod
  USE mhd_mod
  USE rt_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: rmhd_t
    type(rt_t)                        :: rt
  contains
    procedure :: init
    procedure :: update
  end type
! ----------------------------------------------------------------------
CONTAINS

!> =====================================================================
!> Initialization of radiative transfer, preparing for two main tasks:
!> 1) most of the time, just inter-/extrapolate in the exising RT solutions
!> 2) on the rt%dtime intervals, solve the RT and add the new net
!>    heating / cooling in a time slot
!> =====================================================================
SUBROUTINE init (self)
  class (rmhd_t) :: self
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t_%init')
  !---------------------------------------------------------------------
  ! Initialize the MHD data type, and allocate arrays for RT
  !---------------------------------------------------------------------
  call self%mhd_t%init
  !---------------------------------------------------------------------
  ! Initialize the rt_mod object (one could create more than one)
  !---------------------------------------------------------------------
  self%rt%dim      = self%mesh%gn                                       ! RT dimensions
  self%rt%position = self%position
  call self%rt%init                                                     ! allocate weights
  self%rt%position = self%position
  self%mesh%p      = self%position
  call trace%end
END SUBROUTINE init

!> =====================================================================
!> If the RT solution depends on other patches, then adding separate RT
!> tasks into the task list is a better method, so to be general we do
!> this here.   In this case, this update procedure isn't called until
!> the RT solution is already done.
!> =====================================================================
SUBROUTINE update (self)
  class (rmhd_t) :: self
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%update', itimer=itimer)
  call self%mhd_t%update                                                ! get df/dt from MHD
  call self%rt%interpolate (self%time)                                  ! interpolate in time and space
  self%dedt = self%dedt + self%rt%dedt                                  ! add the RT contribution to de/dt
  call trace%end (itimer)
END SUBROUTINE update

END MODULE rmhd_mod
