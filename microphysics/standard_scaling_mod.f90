!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE standard_scaling_mod
  USE units_mod
  USE io_mod
  USE trace_mod
  USE shared_mod
  USE math_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Extend the units_t data type, which handles scaling metadata output
  !-----------------------------------------------------------------------------
  type, public, extends(units_t):: standard_scaling_t
    logical:: in_K=.false.
    real(kind=8):: v=1d0, u=1d0, d=1d0, e=1d0, em=1d0, ev=1d0, ee=1d0, p=1d0, &
      te=1d0, temp=1d0, grav=1d0, kr=1d0, b=1d0, s=1d0, ss=1d0, sv=1d0, sm=1d0
  contains
    procedure:: set
    procedure:: output
  end type
CONTAINS

!===============================================================================
!> Set code units.  Any local clone of this should define (only) the
!> fundamental units of length, time, and mass
!===============================================================================
SUBROUTINE set (self, system, l, t, m, in_K)
  class(standard_scaling_t):: self
  character(*), optional:: system
  real(8), optional:: l, t, m
  logical, optional:: in_K
  !-----------------------------------------------------------------------------
  call trace%begin ('standard_scaling_t%init')
  call io%lock%set ('standard_scaling')
  self%system = shared%option (system, 'code')
  if (self%system=='cgs' .or. self%system=='CGS') then
    self%system = 'CGS'
  else
    self%system = 'SI'
  end if
  !------------------------- defining quantities -------------------------------
  if (present(l))      self%l  = l
  if (present(t))      self%t  = t
  if (present(m) )     self%m  = m
  if (present(in_K)) self%in_K = in_K
  call self%init ()
  !--------------------------- derived quantities ------------------------------
  self%v  = self%l/self%t                             ! velocity
  self%u  = self%v                                    ! alias
  self%l  = self%v*self%t                             ! time
  self%d  = self%m/self%l**3                          ! density
  self%p  = self%d*self%u**2                          ! pressure
  self%kr = 1.0/(self%d*self%l)                       ! opacity
  self%em = self%u**2                                 ! energy per mass
  self%ee = self%em                                   ! alias
  self%ev = self%d*self%em                            ! energy per vol
  self%e  = self%m*self%em                            ! energy 
  !.............................................................................
  if (self%system == 'CGS') then
    self%b  = sqrt(4d0*math%pi*self%d)*self%u         ! magnetic flux density
    if (self%in_K) then
      self%temp = 1d0                                 ! for clarity
    else
      self%temp = (cgs%m_u/self%m)/(cgs%k_b/self%e)   ! T/mu
    end if
    self%grav = cgs%grav*(self%m/self%l**3)*self%t**2 ! constant of gravity
  else
    self%b  = sqrt(si%mu_0*self%d)*self%u             ! magnetic flux density
    if (self%in_K) then
      self%temp = 1d0                                 ! for clarity
    else
      self%temp = (si%m_u/self%m)/(si%k_b/self%e)     ! T/mu
    end if
    self%grav = si%grav*(self%m/self%l**3)*self%t**2  ! constant of gravity
  end if
  self%sm      = self%em/self%temp                    ! entropy per unit mass
  self%ss      = self%sm                              ! alias
  self%sv      = self%ev/self%temp                    ! entropy per unit volume
  call self%output ()
  call io%lock%unset ('standard_scaling')
  call trace%end ()
END SUBROUTINE set

!===============================================================================
!> Write out unit system info
!===============================================================================
SUBROUTINE output (self)
  class(standard_scaling_t):: self
  character(len=4):: system
  real(8):: l, t, m, d, v, u, b, e, ee, em, ev, temp, grav, sm
  real:: m_sun, r_sun, pc, au, stefan, yr, myr
  namelist /scaling_nml/ system, l, t, d, v, m, e, em, ev, b, temp, grav, sm
  namelist /cgs_nml/ m_sun, r_sun, pc, au, grav, stefan, yr, myr
  namelist /si_nml/ m_sun, r_sun, pc, au, grav, stefan, yr, myr
  !----------------------------------------------------------------------------
  call trace%begin ('units_t%output')
  system = self%system
  if (io%master) then
    print 1,' CODE UNITS:              ('//trim(system)//')             (ASTRO)'
    if (self%l > 3e27) then
      print 1,'               length: ',self%l, self%l/(1e9*cgs%pc), 'Gpc'
    else if (self%l > 3e24) then
      print 1,'               length: ',self%l, self%l/(1e6*cgs%pc), 'Mpc'
    else if (self%l > 3e21) then
      print 1,'               length: ',self%l, self%l/(1e3*cgs%pc), 'kpc'
    else if (self%l > 3e18) then
      print 1,'               length: ',self%l, self%l/(1e0*cgs%pc), 'pc'
    else if (self%l > 1.5e13) then
      print 1,'               length: ',self%l, self%l/(1e0*cgs%AU), 'AU'
    else
      print 1,'               length: ',self%l
    end if
    if (self%t > 3e13) then
      print 1,'                 time: ',self%t, self%t/(1e6*cgs%yr), 'Myr'
    else if (self%t > 3e10) then
      print 1,'                 time: ',self%t, self%t/(1e3*cgs%yr), 'kyr'
    else if (self%t > 3e7) then
      print 1,'                 time: ',self%t, self%t/(1e0*cgs%yr), 'yr'
    else
      print 1,'                 time: ',self%t
    end if
    if (self%m > 2e32) then
      print 1,'                 mass: ',self%m, self%m/cgs%M_sun, 'M_Sun'
    else if (self%m > 2e27) then
      print 1,'                 mass: ',self%m, self%m/cgs%M_Earth, 'M_Earth'
    else
      print 1,'                 mass: ',self%m
    end if
    print 1,'             velocity: ',self%u
    print 1,'              density: ',self%d
    print 1,'             pressure: ',self%p
    print 1,' energy per unit mass: ',self%em
    print 1,'entropy per unit mass: ',self%sm
    print 1,'magnetic flux density: ',self%b
    print 1,'              gravity: ',self%grav
    print 1,'          temperature: ',self%temp
  1 format(1x,a,1p,e11.3,0p,f16.2,1x,a)
  end if
  !-----------------------------------------------------------------------------
  t  = self%t
  v  = self%v
  d  = self%d
  l  = self%l
  m  = self%m
  b  = self%b
  e  = self%e
  em = self%em
  sm = self%sm
  ee = self%em
  ev = self%ev
  temp = self%temp
  grav = self%grav
  system = trim(self%system)
  write (io%namelist, scaling_nml)
  call self%units_t%output
END SUBROUTINE output

END MODULE standard_scaling_mod
  