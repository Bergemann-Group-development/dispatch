!> -----------------------------------------------------------------------------
!> Fundamental constants in CGS and SI units
!> -----------------------------------------------------------------------------
MODULE units_mod
  implicit none
  private
  type, public:: const_t
    character(len=16):: name='not set'
    real(kind=8):: m_u
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
CONTAINS

SUBROUTINE init (self)
  class(const_t):: self
  !-----------------------------------------------------------------------------
  ! CGS units
  !-----------------------------------------------------------------------------
  cgs%name    = 'cgs'
  cgs%m_u     = 1.6726d-24
  cgs%amu     = 1.66054d-24
  cgs%k_b     = 1.3807d-16
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
  cgs%pi      = asin(1.0) * 2.0
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
  si%pi      = asin(1.0) * 2.0
END SUBROUTINE init

END MODULE units_mod
