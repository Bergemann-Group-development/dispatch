!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  USE parameters_mod
  implicit none
  private
  type, public:: scaling_t
    real(kind=8):: l = 1.0
    real(kind=8):: t = 1.0
    real(kind=8):: m = 1.0
    real(kind=8):: d = 1.0
    real(kind=8):: e = 1.0
    real(kind=8):: p = 1.0
    real(kind=8):: f = 1.0
    real(kind=8):: u = 1.0
    real(kind=8):: grav = 1.0
    real(kind=8):: temp = 1.0
    real(kind=8):: mu = 2.34                             ! mean molecular weight
    real(kind=8):: b
  contains
    procedure:: init
  end type
  type(scaling_t), public:: scaling
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
  self%l = cgs%au                                       ! length
  self%m = cgs%m_sun                                    ! mass
  self%d = self%m/self%l**3                             ! mass density
  self%t = sqrt((experiment_params%R*self%l)**3 / &
           (cgs%grav*experiment_params%mstar*self%m))   ! time
  self%u = self%l/self%t                                ! speed
  self%p = self%d*self%u**2                             ! pressure = energy density
  self%e = self%m*self%u**2                             ! energy
  self%b = sqrt(self%d)*self%u                          ! magnetic flux
  self%f = self%p*self%u                                ! energy flux = P u
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2     ! constant of gravity
  self%mu   = experiment_params%mu
  self%temp = self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e) ! temperature for give mu
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'     length: ',self%l
    print 1,'   velocity: ',self%u
    print 1,'       mass: ',self%m
    print 1,'       time: ',self%t
    print 1,'    density: ',self%d
    print 1,'   pressure: ',self%p
    print 1,'     energy: ',self%e
    print 1,'energy flux: ',self%f
    print 1,'    gravity: ',self%grav
    print 1,'temperature: ',self%temp
  1 format(1x,a,1p,e11.3)
  end if
END SUBROUTINE init
END MODULE scaling_mod
