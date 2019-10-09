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
    real(kind=8):: m_u, m_h, m_he
    real(kind=8):: amu
    real(kind=8):: k_b
    real(kind=8):: h_p
    real(kind=8):: c
    real(kind=8):: e
    real(kind=8):: stefan
    real(kind=8):: grav
    real(kind=8):: AU
    real(kind=8):: pc
    real(kind=8):: yr
    real(kind=8):: kms
    real(kind=8):: mum
    real(kind=8):: Myr
    real(kind=8):: m_sun
    real(kind=8):: r_sun
    real(kind=8):: m_earth
    real(kind=8):: r_earth
    real(kind=8):: pi
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
  type(units_t), public:: units
CONTAINS

SUBROUTINE init (self)
  class(const_t):: self
  !-----------------------------------------------------------------------------
  ! CGS units
  !-----------------------------------------------------------------------------
  cgs%name    = 'cgs'
  cgs%m_u     = 1.6726219d-24
  cgs%m_h     = cgs%m_u
  cgs%m_he    = 6.65e-24
  cgs%amu     = 1.66054d-24
  cgs%k_b     = 1.380658d-16                                            ! [erg/K]
  cgs%h_p     = 6.6260755d-27
  cgs%c       = 2.999792d10
  cgs%e       = 1.602177e-12
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
  cgs%pi      = math%pi
  !-----------------------------------------------------------------------------
  ! SI units
  !-----------------------------------------------------------------------------
  si%name    = 'SI'
  si%m_u     = 1.6726d-27
  si%amu     = 1.66054d-27
  si%k_b     = 1.3807d-23
  si%h_p     = 6.6260755d-34
  si%c       = 2.999792d8
  si%e       = 1.602177e-19
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
  si%pi      = math%pi
END SUBROUTINE init

!===============================================================================
!> Write out unit system info
!===============================================================================
SUBROUTINE output (self)
  class(units_t):: self
  character(len=4):: system
  real(kind=8):: l, t, m
  namelist /units_nml/ system, l, t, m
  !----------------------------------------------------------------------------
  call trace%begin ('units_t%init')
  system = self%system
  l      = self%l
  t      = self%t
  m      = self%m
  if (omp%master) then
    !$omp critical (write_nml_cr)
    write (io_unit%nml, units_nml)
    !$omp end critical (write_nml_cr)
  end if
  call trace%end()
END SUBROUTINE output

END MODULE units_mod
