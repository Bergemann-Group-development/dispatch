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
    real(kind=8):: b
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
  real, save:: l_pc=250.0, n_h2=4, v_kms=0.18, mu=2.4
  real(kind=8), save:: l, t, m, d, e, ee, p, u, temp, grav
  namelist /scaling_params/ l_pc, n_h2, v_kms, mu
  namelist /scaling_nml/ l_pc, n_h2, v_kms, &
    l, t, m, d, e, ee, p, u, temp, grav, mu
  !----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace_begin('scaling_t%init')
  call io%header ('scaling_t%init:')
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, scaling_params, iostat=iostat)
    call io%namelist_warning ('scaling_params', iostat)
    write (io%output, scaling_params)
    call cgs%init
    u = v_kms*cgs%kms
    l = l_pc*cgs%pc                                  ! length in pc
    d = n_h2*mu*cgs%m_u                              ! mass density in cgs
    u = v_kms*cgs%kms                                ! speed in km/s
    m = d*l**3                                       ! mass
    t = l/u                                          ! time
    p = d*u**2                                       ! pressure = energy density
    e = m*u**2                                       ! energy
    ee = u**2                                        ! energy per unit mass
    mu = mu                                          ! mean molecular weight
    temp = (cgs%m_u*u**2)/cgs%k_b                    ! temperature for mu=1
    grav = cgs%grav*(m/l**3)*t**2                    ! constant of gravity
    write (io%namelist, scaling_nml)
    first_time = .false.
  end if
  call io%lock%unset
  self%l = l
  self%d = d
  self%u = u
  self%m = m
  self%t = t
  self%p = p
  self%e = e
  self%b = sqrt(self%d)*self%u        ! magnetic flux
  self%ee = ee
  self%mu = mu
  self%temp = temp
  self%grav = grav
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

!===============================================================================
!> Function returning mass in solar masses
!===============================================================================
REAL FUNCTION m_sun (self, mass)
  class(scaling_t):: self
  real(8):: mass
  m_sun = mass*self%m/cgs%m_sun
END FUNCTION m_sun

!===============================================================================
!> Function returning time in Myr
!===============================================================================
REAL FUNCTION myr (self, time)
  class(scaling_t):: self
  real(8):: time
  myr = time*self%t/cgs%myr
END FUNCTION myr

END MODULE scaling_mod
