!===============================================================================
!> Define code units for stellar atmospheres, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Extend the units_t data type, which contains only the three fundamental
  ! units: length, mass, and time (cf. CGS and MKS)
  !-----------------------------------------------------------------------------
  type, public, extends(units_t):: scaling_t
    real(kind=8):: d, e, p, u, temp, grav, mu, ee, kr, b
  contains
    procedure:: init
    procedure:: output
  end type
  type(scaling_t), public :: scaling
CONTAINS

!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace%begin ('scaling_t%init')
  first_time = .false.
  call cgs%init
  !--------------------------- defining quantities -----------------------------
  self%system = 'cgs'                                   ! pass on to Python
  self%l  = 1d8                                         ! length = 1Mm = 1e8 cm
  self%d  = 1d-7                                        ! photospheric density
  self%t  = 1d2                                         ! time = 100 sec
  !--------------------------- derived quantities ------------------------------
  self%u  = self%l/self%t                               ! speed
  self%m  = self%d*self%l**3                            ! mass
  self%p  = self%d*self%u**2                            ! pressure = energy density
  self%e  = self%d*self%u**2                            ! energy
  self%ee = self%u**2                                                   ! specific energy
  self%kr = 1.0/(self%d*self%l)                                         ! rosseland opacity [cm2/g]
  self%mu = 1.3
  self%b = self%u*sqrt(cgs%pi*4.0*self%d)
  !self%b = self%l**(-0.5)*self%m**0.5*self%t**(-1)                   ! magnetic flux density - this one is a factor of sqrt(4pi) smaller
  self%temp = 1.0!self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e) ! temperature for given mu
  !self%temp = self%mu*(cgs%m_u/cgs%k_b)*self%ee                      ! temperature for given mu
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2     ! constant of gravity
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'     length: ',self%l
    print 1,'       mass: ',self%m
    print 1,'       time: ',self%t
    print 1,'   velocity: ',self%u
    print 1,'    density: ',self%d
    print 1,'   pressure: ',self%p
    print 1,'     energy: ',self%e
    print 1,'    gravity: ',self%grav
    print 1,'temperature: ',self%temp
    print 1,'magnetic flux density:', self%b
  1 format(1x,a,1p,e11.3)
  end if
  call self%output
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Output scaling units. Do NOT call units%output from here!
!===============================================================================
SUBROUTINE output (self)
  class(scaling_t):: self
  real:: l, t, d
  namelist /scaling_nml/ l, t, d
  !-----------------------------------------------------------------------------
  l = self%l
  t = self%t
  d = self%d
  write (io%namelist,scaling_nml)
END SUBROUTINE output

END MODULE scaling_mod
