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
  USE eos_mod
  USE mhd_mod
  USE rt_mod
  implicit none
  private
  ! --------------------------------------------------------------------
  type, public, extends(mhd_t):: rmhd_t
    type(rt_t)                        :: rt
    real, dimension(:,:,:),   pointer :: lnd, lne
    real, dimension(:,:,:,:), pointer :: rk, src
  contains
    procedure :: init
    procedure :: update
  end type
  real, parameter:: pi4 = 16.*atan(1.)
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
  associate (mx=>self%mesh(1)%gn, my=>self%mesh(2)%gn, mz=>self%mesh(3)%gn)
  allocate (self%lnd (mx,my,mz), &                                      ! log density
            self%lne (mx,my,mz), &                                      ! log internal energy
            self%rk  (mx,my,mz,eos%n_lambda), &                         ! opacity
            self%src (mx,my,mz,eos%n_lambda))                           ! source function
  end associate
  !---------------------------------------------------------------------
  ! Initialize the rt_mod object (one could create more than one)
  !---------------------------------------------------------------------
  self%rt%dim      = self%mesh%gn                                       ! RT dimensions
  self%rt%n_lambda = eos%n_lambda                                       ! number of wavelengths
  self%rt%n_bundle = 1                          ! FIXME                 ! one direction
  call self%rt%init                                                     ! allocate weights
  self%rt%mesh(2)%dsup => self%mesh(2)%dsup                             ! mesh spacing
  self%rt%w_lambda = 1.0/self%rt%n_lambda       ! FIXME                 ! wavelength weights
  self%rt%w_bundle = pi4/2./self%rt%n_bundle    ! FIXME                 ! space angle weights
  !
  call trace%end
END SUBROUTINE init

!> =====================================================================
!> Update the RMHD time derivatives, by calling the MHD module, and
!> adjusting by time-interpolation of the RT memory values.
!>
!> The logic below shows the alternaive where RT tasks are never added
!> to the task queue, but are executed as a subprocess of the RMHD
!> update step.  This works on only if the RT solution does not depend
!> on any other patch.
!> =====================================================================
SUBROUTINE update_both (self)
  class (rmhd_t) :: self
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%update')
  call self%mhd_t%update                                                ! get df/dt from MHD
  !if (.not.self%rt%is_ahead_of (self)) then                             ! if RT is not up-to-date
  !  call self%solve                                                     ! solve the RT for a new slot
  !end if
  call self%rt%interpolate (self%time)                                  ! interpolate in time and space
  print *, shape(self%rt%dedt)
  self%dedt = self%dedt + self%rt%dedt                                  ! add the RT contribution to de/dt
  call trace%end
END SUBROUTINE update_both

!> =====================================================================
!> If the RT solution depends on other patches, then adding separate RT
!> tasks into the task list is a better method, so to be general we do
!> this here.   In this case, this update procedure isn't called until
!> the RT solution is already done.
!> =====================================================================
SUBROUTINE update (self)
  class (rmhd_t) :: self
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%update')
  call self%mhd_t%update                                                ! get df/dt from MHD
  call self%rt%interpolate (self%time)                                  ! interpolate in time and space
  self%dedt = self%dedt + self%rt%dedt                                  ! add the RT contribution to de/dt
  call trace%end
END SUBROUTINE update

END MODULE rmhd_mod
