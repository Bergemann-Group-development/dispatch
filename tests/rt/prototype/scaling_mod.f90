!===============================================================================
!> Define code units for stellar atmospheres, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  implicit none
  private
  type, public:: scaling_t
    real(kind=8):: l, t, m, d, e, p, u, temp, grav, mu
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
  !-----------------------------------------------------------------------------
  if (.not.first_time) return
  if (io%master) print *, &
 '-------------------------------- scaling_mod ---------------------------------'
  first_time = .false.
  call cgs%init
  !--------------------------- defining quantities -----------------------------
  self%l  = 1d8                                         ! length = 1Mm = 1e8 cm
  self%d  = 1d-7                                        ! photospheric density
  self%t  = 1d2                                         ! time = 100 sec
  !--------------------------- derived quantities ------------------------------
  self%u  = self%l/self%t                               ! speed
  self%m  = self%d*self%l**3                            ! mass
  self%p  = self%d*self%u**2                            ! pressure = energy density
  self%e  = self%m*self%u**2                            ! energy
  self%b = sqrt(self%d)*self%u                          ! magnetic flux
  self%mu = 2.
  self%temp = self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e) ! temperature for given mu
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2     ! constant of gravity
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
  1 format(1x,a,1p,e11.3)
  end if
END SUBROUTINE init

END MODULE scaling_mod
