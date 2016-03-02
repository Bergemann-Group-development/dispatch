!> -----------------------------------------------------------------------------
!> Define code units, in terms of (here) CGS units
!> -----------------------------------------------------------------------------
MODULE scaling_mod
  USE units_mod, only: cgs, const
  USE io
  implicit none
  type:: units_t
    real(kind=8):: l = 1.0
    real(kind=8):: t = 1.0
    real(kind=8):: m = 1.0
    real(kind=8):: d = 1.0
    real(kind=8):: e = 1.0
    real(kind=8):: p = 1.0
    real(kind=8):: u = 1.0
    real(kind=8):: grav = 1.0
    real(kind=8):: temp = 1.0
    real(kind=8):: mu = 2.4
  end type
  type(units_t):: code, real
  logical, save:: first_time=.true.
CONTAINS

!> -----------------------------------------------------------------------------
!> Initialize comprehensive code units:
!> -----------------------------------------------------------------------------
SUBROUTINE init_scaling
  implicit none
  real(kind=4):: scale_l, scale_t, scale_d

  !.............................................................................
  if (.not.first_time) return
  first_time = .false.

  ! NOTE: code%l, code%d, code%t are set in `read_ramses`
  code%u = code%l / code%t
  code%m = code%d * code%l**3
  code%p = code%d * code%u**2
  code%e = code%m * code%u**2
  code%grav = cgs%grav*(code%m/code%l**3)*code%t**2
  code%temp = code%mu*(cgs%m_u/code%m)/(cgs%k_b/code%e)
  if (master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'     length: ',code%l
    print 1,'   velocity: ',code%u
    print 1,'       mass: ',code%m
    print 1,'       time: ',code%t
    print 1,'    density: ',code%d
    print 1,'   pressure: ',code%p
    print 1,'     energy: ',code%e
    print 1,'    gravity: ',code%grav
    print 1,'temperature: ',code%temp
  1 format(1x,a,1p,e11.3)
  end if

  ! scale physical constants using the current unit system
  const%name = 'code'

END SUBROUTINE init_scaling

END MODULE scaling_mod
