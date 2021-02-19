!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
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
CONTAINS
!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  logical, save:: first_time=.true.
  !----------------------------------------------------------------------------
  call trace_begin('scaling_t%init')
  if (.not.first_time) return
  first_time = .false.
  self%l = 1.0
  self%m = 1.0
  self%d = 1.0
  self%u = 1.0
  self%t = self%l/self%u
  self%p = self%d*self%u**2
  self%e = self%m*self%u**2
  self%b = sqrt(self%d)*self%u        ! magnetic flux
  self%mu = 1.0
  self%temp = 1.0
  self%grav = 1.0
  if (io%master) then
    print 1,' CODE UNITS:   (dimensionless):'
    print 1,'     length: ',self%l
    print 1,'   velocity: ',self%u
    print 1,'       mass: ',self%m
    print 1,'       time: ',self%t
    print 1,'    density: ',self%d
    print 1,'   pressure: ',self%p
    print 1,'     energy: ',self%e
    print 1,'temperature: ',self%temp
  1 format(1x,a,1p,e11.3)
  end if

  call trace_end
END SUBROUTINE init

END MODULE scaling_mod
