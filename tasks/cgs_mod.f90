!> -----------------------------------------------------------------------------
!> Fundamental constants in CGS units
!> -----------------------------------------------------------------------------
MODULE units_mod
  implicit none
  private
  type, public:: const_t
    character(len=16):: name='cgs'
    real(kind=8):: m_u = 1.6726d-24
    real(kind=8):: amu = 1.66054d-24
    real(kind=8):: k_b = 1.3807d-16
    real(kind=8):: h_p = 6.6260755d-27
    real(kind=8):: c = 2.999792d10
    real(kind=8):: stefan = 5.67d-5
    real(kind=8):: grav = 6.6743d-8
    real(kind=8):: AU = 1.496d13
    real(kind=8):: pc = 3.086d18
    real(kind=8):: yr = 3.15542d7
    real(kind=8):: kms = 1d5
    real(kind=8):: mum = 1d-4
    real(kind=8):: Myr = 3.15542d13
    real(kind=8):: m_sun = 1.989d33
    real(kind=8):: r_sun = 6.9598d10
    real(kind=8):: m_earth = 5.972d27
    real(kind=8):: r_earth = 6.371d8
    real(kind=8):: pi = asin(1.0) * 2.0
  end type
  type(const_t), public:: cgs
END MODULE units_mod
