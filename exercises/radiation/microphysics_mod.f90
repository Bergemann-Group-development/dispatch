MODULE microphysics_mod
  implicit none
  private
  type, public:: microphysics_t
    procedure(source_planck), pointer :: source
    procedure(opacity_grey),  pointer :: opacity
    real :: kappa_grey
  contains
    procedure, private :: source_planck
    procedure          :: opacity_grey
  end type
CONTAINS
!> =============================================================================
SUBROUTINE init (self, source_name, opacity_name)
  implicit none
  class(microphysics_t) :: self
  character(len=12) :: source_name, opacity_name
  !.............................................................................
  select case(source_name)
  case('planck')
    self%source => source_planck
  end select
  !
  select case(opacity_name)
  case('grey')
    self%opacity => opacity_grey
  end select
END SUBROUTINE
!> Set the opacity to be constant 
!> =============================================================================
FUNCTION opacity_grey (self, lambda, density, temperature) result(opacity)
  implicit none
  class(microphysics_t) :: self
  real                  :: lambda                                               ! wavelength
  real, dimension(:,:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2),size(temperature,3)):: opacity
  !
  opacity = self%kappa_grey
END FUNCTION
!> Set the source function equal to the planck function
!> =============================================================================
FUNCTION source_planck (self, lambda, density, temperature) result(source)
  implicit none
  class(microphysics_t) :: self
  real                  :: lambda                                               ! wavelength
  real, dimension(:,:,:):: density, temperature
  real, dimension(size(temperature,1),size(temperature,2),size(temperature,3)):: source
  !
  integer, dimension(3) :: n
  n = shape(temperature)
  call planck(n, lambda, temperature, source)
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
