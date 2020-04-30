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
    real(kind=8):: l, t, m, d, e, ee, p, u, temp, grav, mu
  contains
    procedure:: init
    procedure:: m_sun
    procedure:: myr
  end type
  type(scaling_t), public :: scaling
CONTAINS

!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  real:: l_pc=4.0, d_cgs=3200e-22, v_kms=0.18, mu=2.4
  namelist /scaling_params/ l_pc, d_cgs, v_kms, mu
  !----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace_begin('scaling_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, scaling_params, iostat=iostat)
    write (io%output, scaling_params)
  end if
  call io%lock%unset
  call cgs%init
  self%l = l_pc*cgs%pc                                  ! length in pc
  self%d = d_cgs                                        ! mass density in cgs
  self%u = v_kms*cgs%kms                                ! speed in km/s
  self%m = self%d*self%l**3                             ! mass
  self%t = self%l/self%u                                ! time ~ 5 years
  self%p = self%d*self%u**2                             ! pressure = energy density
  self%e = self%m*self%u**2                             ! energy
  self%ee = self%u**2                                   ! energy per unit mass
  self%mu = mu                                          ! menan molecular weight
  self%temp = self%mu*(cgs%m_u*self%u**2)/cgs%k_b       ! temperature for given mu
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
  call trace_end
END SUBROUTINE init

REAL FUNCTION m_sun (self, mass)
  class(scaling_t):: self
  real(8):: mass
  m_sun = mass*self%m/cgs%m_sun
END FUNCTION m_sun

REAL FUNCTION myr (self, time)
  class(scaling_t):: self
  real(8):: time
  myr = time*self%t/cgs%myr
END FUNCTION myr

END MODULE scaling_mod
