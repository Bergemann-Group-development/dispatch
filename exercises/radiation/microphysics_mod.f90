MODULE microphysics_mod
  USE patches_mod
  USE scaling_mod
  USE io_mod
  implicit none
  private
  type, public:: microphysics_t
    procedure(source_planck), pointer :: source
    procedure(opacity_grey),  pointer :: opacity
    real :: kappa_grey
    real :: lambda
  contains
    procedure, private :: source_planck
    procedure          :: opacity_grey
    procedure          :: dtau
  end type
CONTAINS
!> =============================================================================
SUBROUTINE init (self)
  implicit none
  class(microphysics_t) :: self
  logical :: first_time=.true.                                                  ! only read namelist on first invocation
  character(len=12) :: source='planck', opacity='grey'                          ! choice of source and opacity functions
  real(kind=8)      :: kappa_grey_cgs=1.0                                       ! grey absorption in [cm^2 / g]
  real(kind=8)      :: wavelength_mu=150                                        ! wavelength in micro-meter
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
  self%lambda     = wavelength_mu * 1e-4 / scaling%l 
  !.............................................................................
  select case(source)
  case('planck')
    self%source => source_planck
  end select
  !
  select case(opacity)
  case('grey')
    self%opacity => opacity_grey
  end select
END SUBROUTINE
!> compute dtau given a patch, density, temperature, and direction
!> direction is 1,2,3 for x,y,z
!> =============================================================================
FUNCTION dtau (self, direction, patch, density, temperature)
  implicit none
  class(microphysics_t) :: self
  integer               :: direction
  class(patch_t)        :: patch
  real, dimension(patch%gn(1),patch%gn(2),patch%gn(3)):: density, temperature, dtau
  !
  real, dimension(patch%gn(1),patch%gn(2),patch%gn(3)):: opacity
  !
  opacity = self%opacity(density, temperature)
  ! dtau is the increment in tau from one ray-point to the next.
  ! Therefore dtau is downstaggered and has only a half value on the lower boundary
  dtau = density * opacity
  select case(direction)
  case(1)
    dtau(2:,:,:) = dtau(1:,:,:) + density(1:,:,:) * opacity(1:,:,:)
  case(2)
    dtau(:,2:,:) = dtau(:,1:,:) + density(:,1:,:) * opacity(:,1:,:)
  case(3)
    dtau(:,:,2:) = dtau(:,:,1:) + density(:,:,1:) * opacity(:,:,1:)
  end select
  dtau = 0.5 * patch%ds(direction) * dtau 
END FUNCTION
!> Set the opacity to be constant 
!> =============================================================================
FUNCTION opacity_grey (self, density, temperature) result(opacity)
  implicit none
  class(microphysics_t) :: self
  real, dimension(:,:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2),size(temperature,3)):: opacity
  !
  opacity = self%kappa_grey
END FUNCTION
!> Set the source function equal to the planck function
!> =============================================================================
FUNCTION source_planck (self, density, temperature) result(source)
  implicit none
  class(microphysics_t) :: self
  real, dimension(:,:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2),size(temperature,3)):: source
  !
  real, dimension(size(temperature,1),size(temperature,2),size(temperature,3)):: temperature_K
  real                  :: wavelength_mu
  integer, dimension(3) :: n
  n = shape(temperature)
  wavelength_mu = self%lambda * scaling%l * 1e4
  temperature_K = temperature * scaling%t
  call planck(n, wavelength_mu, temperature_K, source)
  ! source function is returned in micro-meter^-4 Kelvin and has to be converted to code units
  source = real(1.0_8 / (scaling%t / (1e-4_8*scaling%l)**4),kind=kind(source)) * source
END FUNCTION
!> The planck function for a given temperature and wavelength
!> =============================================================================
SUBROUTINE planck (n, lambda, temperature, source)
  implicit none
  integer:: n(3)
  real, dimension(n(1),n(2),n(3)):: temperature, source
  real:: lambda
  !
  integer :: iy, iz
  real, parameter :: cp=1.191042870e7                                           ! first radiation constant: 1.191042870e-16 W m^2 sr^-1
  real, parameter :: c2=1.438777e4                                              ! second radiation constant [micro-meter Kelvin]
  real, dimension(n(1)):: e
!...............................................................................
  do iz=1,n(3)
  do iy=1,n(2)
    e(:) = exp(-c2/(temperature(:,iy,iz)*lambda))
    source(:,iy,iz) = cp*e(:)/((1.0-e(:))*lambda**5)
   end do
   end do
END SUBROUTINE
!> =============================================================================
END MODULE microphysics_mod
