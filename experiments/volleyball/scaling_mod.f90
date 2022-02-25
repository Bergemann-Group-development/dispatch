!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
  implicit none
  private
  type, public:: scaling_t
    real(kind=8):: l, t, m, d, e, p, u, temp, grav, mu, b
  contains
    procedure:: init
    procedure:: m_sun
    procedure:: myr
    procedure:: output
  end type
  real(8):: l_cgs = 1e8
  real(8):: t_cgs = 1e3
  real(8):: d_cgs = 1e-7
  type(scaling_t), public :: scaling
  !type(const_t), public:: code
CONTAINS

!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  integer:: iostat
  namelist /scaling_params/l_cgs, t_cgs, d_cgs
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call io%lock%set ('scaling')
  if (first_time) then
    rewind (io%input)
    read (io%input, scaling_params, iostat=iostat)
    call io%namelist_warning ('scaling_params', iostat)
    write (io%namelist, scaling_params) 
    if (io%master) print *, &
  '-------------------------------- scaling_mod ---------------------------------'
    call cgs%init
    first_time = .false.
  end if
  call io%lock%unset ('scaling')
  !--------------------------- defining quantities -----------------------------
  self%l  = l_cgs                                       ! 1 Mm
  self%t  = t_cgs                                       ! 1000s
  self%d  = d_cgs                                       ! photospheric
  !--------------------------- derived quantities ------------------------------
  self%u  = self%l/self%t                               ! km/s
  self%m  = self%d*self%l**3                            ! mass
  self%p  = self%d*self%u**2                            ! pressure = energy density
  self%e  = self%d*self%u**2                            ! energy
  self%mu = 0.6                                         ! typical?
  self%temp = 1.0                                       ! using K
  self%temp = self%mu*(cgs%m_u/cgs%k_b)*self%e/self%d   ! temperature for given mu
  self%b =1.0
  self%grav = cgs%grav*(self%m/self%l**3)*self%t**2     ! constant of gravity
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'     length: ',self%l
    print 1,'   velocity: ',self%u
    print 1,'       mass: ',self%m
    print 1,'       time: ',self%t
    print 1,'    density: ',self%d
    print 1,'   pressure: ',self%p
    print 1,'energy/mass: ',self%u**2
    print 1,'temperature: ',self%temp
  1 format(1x,a,1p,e11.3)
  end if
  call self%output
  !call units%output
END SUBROUTINE init

REAL FUNCTION m_sun (self, mass)
  class(scaling_t):: self
  real(8):: mass
  m_sun = mass*self%m/cgs%m_sun
END FUNCTION m_sun

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

REAL FUNCTION myr (self, time)
  class(scaling_t):: self
  real(8):: time
  myr = time*self%t/cgs%myr
END FUNCTION myr

END MODULE scaling_mod
