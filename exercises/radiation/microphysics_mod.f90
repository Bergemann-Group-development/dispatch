!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE microphysics_mod
  USE scaling_mod
  USE units_mod
  USE io_mod
  USE planck_mod
  implicit none
  private
  type, public:: microphysics_t
    procedure(opacity_grey),  pointer :: opacity => null()
    real :: kappa_grey
    real :: lambda
  contains
    procedure          :: init
    procedure          :: source
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
  class(microphysics_t) :: self
  logical :: first_time=.true.                  ! namelist on first invocation
  character(len=12) :: opacity='grey'           ! choice of opacity
  real(kind=8)      :: kappa_grey_cgs=1.0       ! grey absorption in [cm^2 / g]
  real(kind=8)      :: wavelength_mu=150        ! wavelength in micro-meter
  namelist /microphysics_params/ opacity, kappa_grey_cgs, wavelength_mu
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
!> Set the opacity to be constant per unit mass
!> =============================================================================
FUNCTION opacity_grey (self, density, temperature) result(opacity)
  class(microphysics_t) :: self
  real, dimension(:,:):: density, temperature
  real, dimension(size(density,1),size(density,2)):: opacity
  !.............................................................................
  opacity = density * self%kappa_grey
END FUNCTION

!> =============================================================================
!> Set the source function equal to a Planck function lookup table
!> =============================================================================
FUNCTION source (self, density, temperature) result(out)
  class(microphysics_t) :: self
  real, dimension(:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2)):: out
  logical, save:: first_time=.true.
  integer, dimension(2) :: n
  !.............................................................................
  if (first_time) then
    first_time = .false.
    call planck%init
  end if
  n = shape(temperature)
  call planck%source (n, self%lambda, temperature, out)
END FUNCTION source

!> =============================================================================
END MODULE microphysics_mod
