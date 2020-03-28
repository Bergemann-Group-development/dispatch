!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  USE math_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Extend the units_t data type, which contains only the three fundamental
  ! units: length, mass, and time (cf. CGS and MKS)
  !-----------------------------------------------------------------------------
  type, public, extends(units_t):: scaling_t
    real(kind=8):: d, e, ee, p, u, temp, grav, mu, kr, b
  contains
    procedure:: init
  end type
  type(scaling_t), public:: scaling
  public:: cgs
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
  if (io%master) print '(a)', &
 '-------------------------------- scaling_mod --------------------------------'
  first_time = .false.
  call cgs%init
  !--------------------------- defining quantities -----------------------------
  self%system = 'cgs'                               ! pass on to Python
  self%l  = cgs%au                                  !  length unit: 1 AU
  self%d  = 1e-9                                    ! density unit: 1e-9 g/cm**3
  self%t  = cgs%yr                                  !    time unit: 1 yr
  !--------------------------- derived quantities -----------------------------
  self%u  = self%l/self%t                           ! speed
  self%m  = self%d*self%l**3                        ! mass
  self%p  = self%d*self%u**2                        ! pressure [dyne / cm2]
  self%kr = 1.0/(self%d*self%l)                     ! opacity [cm2/g]
  self%ee = self%u**2                               ! specific energy
  self%e  = self%m*self%u**2                        ! energy
  self%b  = self%l**(-0.5)*self%m**0.5*self%t**(-1) ! magnetic flux density
  self%mu = 2.4                                     ! mean molecular mas
  self%temp = self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e) ! temperature, given mu
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2 ! constant of gravity
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'               length: ',self%l
    print 1,'                 mass: ',self%m
    print 1,'                 time: ',self%t
    print 1,'             velocity: ',self%u
    print 1,'              density: ',self%d
    print 1,'             pressure: ',self%p
    print 1,'magnetic flux density: ',self%b
    print 1,'               energy: ',self%e
    print 1,'              gravity: ',self%grav
    print 1,'          temperature: ',self%temp
  1 format(1x,a,1p,e11.3)
  end if
  call self%output
  call trace%end()
END SUBROUTINE init
END MODULE scaling_mod
