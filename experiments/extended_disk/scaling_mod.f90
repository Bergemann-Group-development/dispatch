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
    real(kind=8):: l, t, v, m, d, e, ee, p, u, temp, grav, mu
    real(kind=8):: b
  contains
    procedure:: init
  end type
  real(8):: l_au=1d0, v_kms=1d0, d_cgs=1d9
  type(scaling_t), public :: scaling
CONTAINS

!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  class(scaling_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  real, save:: mu=2.33
  real(kind=8), save:: l, t, v, m, d, e, ee, p, u, temp, grav
  namelist /scaling_params/ l_au, v_kms, d_cgs
  namelist /scaling_nml/ &
    l, t, m, d, e, ee, p, u, temp, grav, mu
  !----------------------------------------------------------------------------
  if (.not.first_time) return
  call trace_begin('scaling_t%init')
  !----------------------------------------------------------------------------
  ! Note that this is called from setup_mod.f90, while the code is still
  ! single-threaded, so there should be no critical region or io%lock here
  !----------------------------------------------------------------------------
  if (first_time) then
    rewind (io%input)
    read (io%input, scaling_params, iostat=iostat)
    call io%namelist_warning ('scaling_params', iostat)
    write (io%namelist, scaling_params)
    call cgs%init
    l = l_au*cgs%au                                  ! length in cm
    v = v_kms*1e5                                    ! velocity in km/s
    t = l/v                                          ! time in s
    d = d_cgs                                        ! density in g/cm^3
    u = l/t                                          ! speed
    m = d*l**3                                       ! mass
    p = d*u**2                                       ! pressure = energy density
    e = m*u**2                                       ! energy
    ee = u**2                                        ! energy per unit mass
    temp = (cgs%m_u*u**2)/cgs%k_b                    ! temperature for mu=1
    grav = cgs%grav*(m/l**3)*t**2                    ! constant of gravity
    write (io_unit%nml, scaling_nml)
    first_time = .false.
  end if
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
    print 1,'    density: ',self%d
    print 1,'       mass: ',self%m
    print 1,'       time: ',self%t
    print 1,'   pressure: ',self%p
    print 1,'     energy: ',self%e
    print 1,'    gravity: ',self%grav
    print 1,'temperature: ',self%temp
  1 format(3(1x,a,1p,e11.3,5x))
  end if
  call trace_end
END SUBROUTINE init

END MODULE scaling_mod
