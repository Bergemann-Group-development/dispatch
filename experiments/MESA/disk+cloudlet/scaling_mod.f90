!===============================================================================
!> Define code units, in terms of (here) CGS units
!===============================================================================
MODULE scaling_mod
  USE units_mod
  USE io_mod
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
    real(kind=8):: mu = 2.4                             ! mean molecular weight
    real(kind=8):: m_h
    real(kind=8):: k_B
  contains
    procedure:: init
  end type
  type(scaling_t), public:: scaling
CONTAINS
!===============================================================================
!> Initialize code units:
!===============================================================================
SUBROUTINE init (self)
  USE math_mod
  class(scaling_t):: self
  real(kind=8):: l0=3.12e16
  integer:: iostat
  logical, save:: first_time=.true.
  namelist / scaling_params / l0
  !----------------------------------------------------------------------------
  if (.not.first_time) return
  call io%lock%set ('scaling_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, scaling_params, iostat=iostat)
    call io%namelist_warning ('scaling_params', iostat)
    write (io%namelist, scaling_params)
  end if
  call io%lock%unset ('scaling_t%init')
  first_time = .false.

  call cgs%init
  self%l = l0                                           ! length; cm
  self%m = cgs%m_sun                                    ! mass; g
  self%d = self%m/self%l**3                             ! mass density; g cm**-3
  self%grav = 0.25 / math%pi                            ! gravitational constant
  self%t = sqrt(self%grav / cgs%grav / self%d)          ! time; s
  self%u = self%l/self%t                                ! speed
  self%p = self%d*self%u**2                             ! pressure = energy density
  self%e = self%m*self%u**2                             ! energy
  self%f = self%p*self%u                                ! energy flux = P u
  self%temp = self%mu*(cgs%m_u/self%m)/(cgs%k_b/self%e) ! temperature for give mu
  self%k_B  = cgs%k_B                                   ! Boltzmann constant in cgs
  self%m_h  = cgs%m_h                                   ! proton mass in g
  if (io%master) then
    print 1,' CODE UNITS:   (to cgs):'
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
