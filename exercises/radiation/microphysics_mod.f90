!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE microphysics_mod
  USE patches_mod
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
    procedure          :: opacity_grey
    procedure          :: dtau
  end type
  type(microphysics_t), public:: microphysics
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
  self%lambda     = wavelength_mu * cgs%mum / scaling%l
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
!> =============================================================================
!> compute dtau given a patch, density, temperature, and direction
!> direction is 1,2,3 for x,y,z
!> =============================================================================
FUNCTION dtau (self, direction, patch, density, temperature)
  implicit none
  class(microphysics_t) :: self
  integer               :: direction
  class(patch_t)        :: patch
  real, dimension(patch%n(1),patch%n(2),patch%n(3)):: density, temperature, dtau
  real, dimension(patch%n(1),patch%n(2),patch%n(3)):: opacity
  !.............................................................................
associate (n=>patch%n)
  opacity = self%opacity(density, temperature)
  !-----------------------------------------------------------------------------
  ! dtau is the increment in tau from one ray-point to the next. Therefore dtau
  ! is downstaggered and has only a half value on the lower edge of the patch
  !-----------------------------------------------------------------------------
  dtau = density * opacity
  select case(direction)
  case(1)
    dtau(2:,:,:) = dtau(1:n(1)-1,:,:) + density(1:n(1)-1,:,:) * opacity(1:n(1)-1,:,:)
  case(2)
    dtau(:,2:,:) = dtau(:,1:n(2)-1,:) + density(:,1:n(2)-1,:) * opacity(:,1:n(2)-1,:)
  case(3)
    dtau(:,:,2:) = dtau(:,:,1:n(3)-1) + density(:,:,1:n(3)-1) * opacity(:,:,1:n(3)-1)
  end select
end associate
  dtau = 0.5 * patch%ds(direction) * dtau 
END FUNCTION
!> =============================================================================
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
!> =============================================================================
!> Set the source function equal to the planck function
!> =============================================================================
FUNCTION source_planck (self, density, temperature) result(source)
  implicit none
  class(microphysics_t) :: self
  real, dimension(:,:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2),size(temperature,3)):: source
  !
  integer, dimension(3) :: n
  n = shape(temperature)
  call planck(n, self%lambda, temperature, source)
END FUNCTION
!> =============================================================================
!> The planck function for a given temperature and wavelength input in code units
!> output is also in code units.  We use lambda in micrometers, to avoid underflow
!> and divide fault in the expression when computing lambda**5
!> =============================================================================
SUBROUTINE planck (n, lambda, temperature, source)
  implicit none
  integer:: n(3)
  real:: lambda, lambda_mum
  real, dimension(n(1),n(2),n(3)):: temperature, source
  !
  integer :: iy, iz
  real(8) :: c1, c2
  real, dimension(n(1)):: e
!...............................................................................
  c1 = 2d0*cgs%h_p*cgs%c**2/(cgs%mum**5)                ! 2*h*c**2/lambda**5
  c2 = cgs%h_p*cgs%c/(cgs%k_b*cgs%mum*scaling%temp)     ! h*c/(lambda*T)
  lambda_mum = lambda*scaling%l/cgs%mum                 ! lambda in micrometers
  do iz=1,n(3)
  do iy=1,n(2)
    e(:) = exp(-c2/(temperature(:,iy,iz)*lambda_mum))
    source(:,iy,iz) = c1*e(:)/((1.0-e(:))*lambda_mum**5)
  end do
  end do
  print*,source(1,1,1)
END SUBROUTINE
!> =============================================================================
END MODULE microphysics_mod
