!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE math_mod
  implicit none
  private
  type, public:: scaling_t
    real(kind=8):: l, t, m, d, e, p, u, temp, grav, mu, stefan
    real(kind=8):: b
  contains
    procedure:: init
  end type
  type(scaling_t), public :: scaling
  !type(const_t), public:: code
CONTAINS

!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  logical, save:: first_time=.true.
  if (.not.first_time) return
  first_time = .false.
  call cgs%init
  self%l = cgs%r_earth                                  ! length = 1 R_Earth
  self%m = cgs%m_earth                                  ! mass
  self%d = self%m/self%l**3                             ! mass density
  self%u = cgs%kms                                      ! speed = 1 km/s
  self%t = self%l/self%u                                ! time ~ 6700 sec < 2 hours
  self%p = self%d*self%u**2                             ! pressure = energy density
  self%e = self%m*self%u**2                             ! energy
  self%b = sqrt(self%d)*self%u                          ! magnetic flux
  self%mu = 2.0                                         ! molecular weight
  self%temp = self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e) ! temperature for given mu
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2     ! constant of gravity
  self%stefan = cgs%stefan/(scaling%p*scaling%u)*self%temp**4
  !code%r_earth = cgs%r_earth/self%l
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'     length: ',self%l
    print 1,'   velocity: ',self%u
    print 1,'       mass: ',self%m
    print 1,'       time: ',self%t
    print 1,'    density: ',self%d
    print 1,'   pressure: ',self%p
    print 1,'     energy: ',self%e
    print 1,'    gravity: ',self%grav
    print 1,'temperature: ',self%temp
    print 1,'     stefan: ',self%stefan
  1 format(1x,a,1p,e11.3)
    print *, 'omega consistency:', &
      sqrt(cgs%grav*cgs%m_sun/cgs%au**3)/ &
      (2.*math%pi/cgs%yr) - 1.d0
  end if
END SUBROUTINE init

END MODULE scaling_mod
