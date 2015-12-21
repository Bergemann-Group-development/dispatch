!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE microphysics_mod
  USE scaling_mod
  USE units_mod
  USE io_mod
  implicit none
  private
  type, public:: microphysics_t
    procedure(source_planck), pointer :: source => null()
    procedure(opacity_grey),  pointer :: opacity => null()
    real :: kappa_grey
    real :: lambda
  contains
    procedure          :: init
    procedure, private :: source_planck
    procedure, private :: source_teff
    procedure          :: opacity_grey
    procedure          :: dtau
  end type
  type(microphysics_t), public:: microphysics
CONTAINS

!> =============================================================================
!> Initialize a microphysics object, which includes a choice of source function
!> and opacity, and one or more methods to compute optical depth increments
!> =============================================================================
SUBROUTINE init (self)
  implicit none
  class(microphysics_t) :: self
  logical :: first_time=.true.                  ! namelist on first invocation
  character(len=12) :: source='planck'          ! choice of source function
  character(len=12) :: opacity='grey'           ! choice of opacity
  real(kind=8)      :: kappa_grey_cgs=1.0       ! grey absorption in [cm^2 / g]
  real(kind=8)      :: wavelength_mu=150        ! wavelength in micro-meter
  namelist /microphysics_params/ source, opacity, kappa_grey_cgs, wavelength_mu
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, microphysics_params)
    if (io%master) write (*,microphysics_params)
  end if
  self%kappa_grey = kappa_grey_cgs / (scaling%l**2 / scaling%m) 
  !-----------------------------------------------------------------------------
  ! For multi frequency RT lambda may be redefined by calling function
  !-----------------------------------------------------------------------------
  self%lambda     = wavelength_mu * cgs%mum / scaling%l
  !.............................................................................
  select case(source)
  case('planck')
    self%source => source_planck
  case('teff')
    self%source => source_teff
  end select
  !
  select case(opacity)
  case('grey')
    self%opacity => opacity_grey
  end select
END SUBROUTINE

!> =============================================================================
!> compute dtau given density, temperature, and step ds along the ray
!> =============================================================================
FUNCTION dtau (self, ds, opacity)
  implicit none
  class(microphysics_t) :: self
  real                  :: ds
  real, dimension(:,:)  :: opacity
  real, dimension(size(opacity,1),size(opacity,2)):: dtau
  !-----------------------------------------------------------------------------
  ! dtau is the increment in tau from one ray-point to the next. Therefore dtau
  ! is downstaggered and has only a half value on the lower edge of the patch
  !-----------------------------------------------------------------------------
associate (n=>size(opacity,2))
  dtau = opacity
  dtau(:,2:) = dtau(:,1:n-1) + opacity(:,1:n-1)
  dtau = 0.5*ds*dtau 
end associate
END FUNCTION

!> =============================================================================
!> Set the opacity to be constant 
!> =============================================================================
FUNCTION opacity_grey (self, density, temperature) result(opacity)
  implicit none
  class(microphysics_t) :: self
  real, dimension(:,:):: density, temperature
  real, dimension(size(density,1),size(density,2)):: opacity
  !.............................................................................
  opacity = density * self%kappa_grey
END FUNCTION

!> =============================================================================
!> Set the source function equal to the planck function
!> =============================================================================
FUNCTION source_planck (self, density, temperature) result(source)
  implicit none
  class(microphysics_t) :: self
  real, dimension(:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2)):: source
  !.............................................................................
  integer, dimension(2) :: n
  n = shape(temperature)
  call planck(n, self%lambda, temperature, source)
END FUNCTION

!> =============================================================================
!> The planck function, normalized to give radiation per unit \ln\lambda,  for
!> a given temperature and wavelength. Input and output are in code units, but
!> we use lambda in micrometers, to avoid underflow and divide fault in the
!> expression when computing lambda**4.
!> =============================================================================
SUBROUTINE planck (n, lambda, temperature, source)
  implicit none
  integer:: n(2)
  real:: lambda, lambda_mum, lambda_4
  real, dimension(n(1),n(2)):: temperature, source
  integer :: i
  real(8) :: c1, c2
  real, dimension(n(1)):: e
!...............................................................................
  c1 = 2d0*cgs%h_p*cgs%c**2/(cgs%mum**4)                ! 2*h*c**2/lambda**4
  c2 = cgs%h_p*cgs%c/(cgs%k_b*cgs%mum*scaling%temp)     ! h*c/(lambda*T)
  lambda_mum = lambda*scaling%l/cgs%mum                 ! lambda in micrometers
  lambda_4 = lambda_mum**4
  do i=1,n(2)
    e(:) = exp(-c2/(temperature(:,i)*lambda_mum))
    !e(:) = c2/(temperature(:,i)*lambda_mum)
    source(:,i) = c1*e(:)/((1.0-e(:))*lambda_4)
  end do
END SUBROUTINE
!> =============================================================================

!> =============================================================================
!> Set the source function equal to the planck function
!> =============================================================================
FUNCTION source_teff (self, density, temperature) result(source)
  implicit none
  class(microphysics_t) :: self
  real, dimension(:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2)):: source
  !.............................................................................
  integer, dimension(2) :: n
  n = shape(temperature)
  call teff (n, self%lambda, temperature, source)
END FUNCTION

!> =============================================================================
!> The integrated Planck function is proportional to Stefan's constant times the
!> fourth power of the temperature.  This is the relevant source function when
!> the opacity is constant as a function of wavelength, so one can integrate 
!> over wavelengthe before solving the radiative transfer equation
!> =============================================================================
SUBROUTINE teff (n, lambda, temperature, source)
  implicit none
  integer:: n(2)
  real, dimension(n(1),n(2)):: temperature, source
  real:: lambda, c
!...............................................................................
  c = cgs%stefan/cgs%pi*scaling%temp**4/(scaling%e/(scaling%l**2*scaling%t))
  source = c*temperature**4
END SUBROUTINE
!> =============================================================================
END MODULE microphysics_mod
