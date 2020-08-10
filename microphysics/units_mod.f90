
!===============================================================================
!> Fundamental constants in CGS and SI units
!===============================================================================
MODULE units_mod
  USE io_unit_mod
  USE math_mod
  USE omp_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Data type holding constants of nature for a system
  !-----------------------------------------------------------------------------
  type, public:: const_t
    character(len=16):: name='not set'
    real(kind=8):: m_u, m_h, m_he, m_e ! proton, hydrogen, helium, electron mass
    real(kind=8):: amu                 ! Atomic mass unit
    real(kind=8):: k_b                 ! Boltzmann's constant
    real(kind=8):: h_p                 ! Planck's constant
    real(kind=8):: c                   ! speed of light
    real(kind=8):: e                   ! electron charge
    real(kind=8):: stefan              ! Stefan-Boltzmann constant
    real(kind=8):: grav                ! constant of gravity
    real(kind=8):: AU                  ! astronomical unit
    real(kind=8):: pc                  ! parsec
    real(kind=8):: yr                  ! year
    real(kind=8):: kms                 ! kilometres per second
    real(kind=8):: mum                 ! micrometer
    real(kind=8):: Myr                 ! Megayear
    real(kind=8):: m_sun               ! Solar mass
    real(kind=8):: r_sun               ! Solar radius
    real(kind=8):: m_earth             ! Earth mass
    real(kind=8):: r_earth             ! Earth radius
    real(kind=8):: pi
    real(kind=8):: mu_0                ! vacuum permeability
    real(kind=8):: eps_0               ! vacuum permittivity
  contains
    procedure:: init
  end type
  type(const_t), public:: cgs, si
  !-----------------------------------------------------------------------------
  ! Data type for transmitting system info as meta-data
  !-----------------------------------------------------------------------------
  type, public:: units_t
    character(len=4):: system='code'
    real(kind=8):: l=1d0
    real(kind=8):: t=1d0
    real(kind=8):: m=1d0
  contains
    procedure:: output
  end type
CONTAINS

SUBROUTINE init (self)
  class(const_t):: self
  type(units_t):: units
  !-----------------------------------------------------------------------------
  ! CGS units
  !-----------------------------------------------------------------------------
  call trace%begin ('const_t%init')
  cgs%name    = 'cgs'
  cgs%m_u     = 1.6726219d-24
  cgs%m_e     = 9.10938356d-28
  cgs%m_h     = cgs%m_u
  cgs%m_he    = 6.65e-24
  cgs%amu     = 1.66054d-24
  cgs%k_b     = 1.380658d-16           ! [erg/K]
  cgs%h_p     = 6.6260755d-27
  cgs%c       = 2.999792d10
  cgs%e       = 4.80320451d-10  !1.602177e-12
  cgs%stefan  = 5.67d-5
  cgs%grav    = 6.6743d-8
  cgs%AU      = 1.496d13
  cgs%pc      = 3.086d18
  cgs%yr      = 3.15542d7
  cgs%kms     = 1d5
  cgs%mum     = 1d-4
  cgs%Myr     = 3.15542d13
  cgs%m_sun   = 1.989d33
  cgs%r_sun   = 6.9598d10
  cgs%m_earth = 5.972d27
  cgs%r_earth = 6.371d8
  cgs%mu_0    = 1.0
  cgs%eps_0   = 1.0
  cgs%pi      = math%pi
  !-----------------------------------------------------------------------------
  ! SI units
  !-----------------------------------------------------------------------------
  si%name    = 'SI'
  si%m_u     = 1.6726d-27
  si%m_e     = 9.10938356d-31
  si%mu_0    = 1.2566370614359173e-6   ! vacuum permeability
  si%amu     = 1.66054d-27
  si%k_b     = 1.3807d-23
  si%h_p     = 6.6260755d-34
  si%c       = 2.999792d8
  si%e       = 1.602177e-19            ! C
  si%stefan  = 5.67d-8
  si%grav    = 6.6743d-11
  si%AU      = 1.496d11
  si%pc      = 3.086d16
  si%yr      = 3.15542d7
  si%kms     = 1d3
  si%mum     = 1d-6
  si%Myr     = 3.15542d13
  si%m_sun   = 1.989d30
  si%r_sun   = 6.9598d8
  si%m_earth = 5.972d25
  si%r_earth = 6.371d6
  si%mu_0    = math%pi4*1d-7
  si%eps_0   = 8.8542d-12
  si%pi      = math%pi
  call units%output
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Write out unit system info
!===============================================================================
SUBROUTINE output (self)
  class(units_t):: self
  character(len=4):: system
  real(kind=8):: l, t, m
  real:: m_sun, r_sun, pc, au, grav, stefan, yr, myr
  namelist /units_nml/ system, l, t, m
  namelist /cgs_nml/ m_sun, r_sun, pc, au, grav, stefan, yr, myr
  namelist /si_nml/ m_sun, r_sun, pc, au, grav, stefan, yr, myr
  !----------------------------------------------------------------------------
  call trace%begin ('units_t%output')
  system = self%system
  l      = self%l
  t      = self%t
  m      = self%m
  if (omp%master) then
    !$omp critical (write_nml_cr)
    write (io_unit%nml, units_nml)
    m_sun  = cgs%m_sun
    r_sun  = cgs%r_sun
    pc     = cgs%pc
    au     = cgs%au
    grav   = cgs%grav
    stefan = cgs%stefan
    yr     = cgs%yr
    myr    = cgs%myr
    write (io_unit%nml, cgs_nml)
    m_sun  = si%m_sun
    r_sun  = si%r_sun
    pc     = si%pc
    au     = si%au
    grav   = si%grav
    stefan = si%stefan
    yr     = si%yr
    myr    = si%myr
    write (io_unit%nml, si_nml)
    !$omp end critical (write_nml_cr)
  end if
  call trace%end ()
END SUBROUTINE output

END MODULE units_mod
