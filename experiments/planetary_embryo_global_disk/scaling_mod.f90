!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE io_mod
  USE trace_mod
  USE units_mod
  implicit none
  private
  type, public:: scaling_t
    real(kind=8):: l, t, m, d, e, p, u, temp, grav, mu
    real(kind=8):: b
  contains
    procedure:: init
  end type
  type(scaling_t), public :: scaling
CONTAINS
!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  logical, save:: first_time=.true.
  character(len=15) :: fstring
  !----------------------------------------------------------------------------
  call trace_begin('scaling_t%init')
  if (.not.first_time) return
  first_time = .false.
  call cgs%init
  self%mu = 2.3                         ! mean molecular weight
  self%l = cgs%au                       ! length
  self%m = cgs%m_sun                    ! mass
  self%t = cgs%yr                       ! time
  self%d = self%m/self%l**2             ! mass density
  self%u = self%l / self%t              ! speed
  self%p = self%d*self%u**2             ! pressure/energy density
  self%e = self%m*self%u**2             ! energy
  self%b = sqrt(self%d)*self%u          ! magnetic flux
  self%temp = self%mu * (cgs%m_u / self%m) / (cgs%k_b / self%e) ! temperature
  self%grav = self%l**3 / (self%m * self%t**2)
  self%grav = cgs%grav / self%grav  ! gravitational constant (scaled)
  if (io%master) then
    fstring = '(1x,a,1p,e11.3)'
    print fstring,' CODE UNITS:   (cgs):'
    print fstring,'     length: ',self%l
    print fstring,'   velocity: ',self%u
    print fstring,'       mass: ',self%m
    print fstring,'       time: ',self%t
    print fstring,'    density: ',self%d
    print fstring,'   pressure: ',self%p
    print fstring,'     energy: ',self%e
    print fstring,'temperature: ',self%temp
  end if

  call trace_end
END SUBROUTINE init

END MODULE scaling_mod
