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
  type, public:: scaling_t
    real(kind=8):: l, t, m, d, e, ee, p, u, temp, grav, mu, gamma
  contains
    procedure:: init
    procedure:: m_sun
    procedure:: myr
  end type
  real:: mass_msun = 1.0
  real:: orbit_au = 1.0
  real:: surface_density_cgs = 1.0
  real:: temperature_K = 300.
  real:: gamma = 1.0
  real:: mu = 2.4
  type(scaling_t), public :: scaling
CONTAINS

!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self, l_o, d_o, u_o, mu_o)
  class(scaling_t):: self
  real, optional:: l_o, d_o, u_o, mu_o
  real, save:: l=1, d=1, u=1, mu=1
  real:: t, m, e, ee, p
  real:: l_au, v_kms, temp, grav, n_H2
  real:: d_midplane_cgs, orbit_cgs, orbit_speed_cgs, sound_speed_cgs, scale_height_cgs
  integer:: iostat
  namelist /scaling_params/ mass_msun, orbit_au, surface_density_cgs, temperature_K, &
    mu, gamma
  namelist /scaling_nml/ l_au, v_kms, &
    l, t, m, d, e, ee, p, u, temp, grav, mu
  logical:: first_time = .true.
  !----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace_begin('scaling_t%init')
  call io%lock%set
  if (first_time) then
    !-----------------------------------------------------------------------------
    ! Compute CGS values of orbit speed, sound speed, and disk scale height
    !-----------------------------------------------------------------------------
    rewind (io%input)
    read (io%input, scaling_params, iostat=iostat)
    call io%namelist_warning ('scaling_params', iostat)
    write (io%output, scaling_params)
    call cgs%init
    orbit_cgs = orbit_au*cgs%au
    mass_msun = mass_msun*cgs%m_sun
    orbit_speed_cgs = sqrt(cgs%grav*mass_msun/orbit_cgs)
    sound_speed_cgs = sqrt(gamma*cgs%k_B*temperature_K/(mu*cgs%m_u))
    scale_height_cgs = orbit_cgs*sound_speed_cgs/orbit_speed_cgs
    d_midplane_cgs = surface_density_cgs/(scale_height_cgs*2.0*math%pi)  
    !-----------------------------------------------------------------------------
    ! Choose code units = scale height, midplane density, and sound speed
    !-----------------------------------------------------------------------------
    l = scale_height_cgs
    d = d_midplane_cgs
    u = sound_speed_cgs
    !-----------------------------------------------------------------------------
    ! Derived values
    !-----------------------------------------------------------------------------
    n_H2 = d/(mu*cgs%m_u)                            ! H2 number density
    m = d*l**3                                       ! mass unit
    t = l/u                                          ! time unit
    p = d*u**2                                       ! pressure = energy density
    e = m*u**2                                       ! energy
    ee = u**2                                        ! energy per unit mass
    temp = (cgs%m_u*u**2)/cgs%k_b                    ! temperature for mu=1
    grav = cgs%grav*(m/l**3)*t**2                    ! constant of gravity
    l_au = l/cgs%au                                  ! length unit in AU
    n_H2 = d/(mu*cgs%m_u)                            ! H2 density
    write (io_unit%nml, scaling_nml)                 ! pass on to Python
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
  self%ee = ee
  self%mu = mu                               ! assumed mean molecular weight
  self%gamma = gamma                         ! assumed ideal gas gamma
  self%temp = temp                           ! temerature unit in K
  self%grav = grav                           ! constant of gravity in code units
  if (io%master) then
    print 1,' CODE UNITS:   (cgs):'
    print 1,'     length: ',self%l, ' R_Earth:', self%l/cgs%r_earth
    print 1,'   velocity: ',self%u, '    km/s:', self%u/cgs%kms
    print 1,'       mass: ',self%m, ' M_Earth:', self%m/cgs%m_earth
    print 1,'       time: ',self%t, '      yr:', self%t/cgs%yr
    print 1,'   pressure: ',self%p, '     bar:', self%p/1e6
    print 1,'    density: ',self%d
    print 1,'     energy: ',self%e
    print 1,'    gravity: ',self%grav
    print 1,'temperature: ',self%temp
  1 format(1x,2(a,1p,e11.3,4x))
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
