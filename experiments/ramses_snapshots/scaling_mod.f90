!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: scaling_t
    real(kind=8):: l, t, m, d, e, p, u, temp, grav, mu
    real(kind=8):: b
  contains
    procedure:: init
  end type
  type(scaling_t), public :: scaling
  type(const_t), public:: code
CONTAINS

!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  logical, save:: first_time=.true.
  real:: l_pc=40.0, d_cgs=4e-22, v_kms=1.0
  namelist /scaling_params/ l_pc, d_cgs, v_kms
  !----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace_begin('scaling_t%init')
  first_time = .false.
  call cgs%init
  self%l = l_pc*cgs%pc                                  ! length in pc
  self%d = d_cgs                                        ! mass density in cgs
  self%u = v_kms*cgs%kms                                ! speed in km/s
  self%m = self%d*self%l**3                             ! mass
  self%t = self%l/self%u                                ! time ~ 5 years
  self%p = self%d*self%u**2                             ! pressure = energy density
  self%e = self%m*self%u**2                             ! energy
  self%b = sqrt(self%d)*self%u                          ! magnetic flux
  self%mu = 2d0
  self%temp = self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e) ! temperature for given mu
  code%grav = cgs%grav*(self%m/self%l**3)*self%t**2     ! constant of gravity
  code%r_earth = cgs%r_earth/self%l
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'     length: ',self%l
    print 1,'   velocity: ',self%u
    print 1,'       mass: ',self%m
    print 1,'       time: ',self%t
    print 1,'    density: ',self%d
    print 1,'   pressure: ',self%p
    print 1,'     energy: ',self%e
    print 1,'    gravity: ',code%grav
    print 1,'temperature: ',self%temp
  1 format(1x,a,1p,e11.3)
  end if
  call trace_end
END SUBROUTINE init

END MODULE scaling_mod
