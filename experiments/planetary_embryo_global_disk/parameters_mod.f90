!> =============================================================================
MODULE parameters_mod
  USE units_mod
  USE scaling_mod
  USE planetesimal_mod
  implicit none
  private

  type, public :: experiment_parameters_t
    real(8) :: mstar = 1.0  ! mass of the star (in solar masses)
    real(4) :: kill_rmin = 0.5   ! inner and outer radii where one starts to
    real(4) :: kill_rmax = 2.1   ! apply the "wave killing" boundary conditions
    real(8) :: aspect_rat = 0.05 ! disk aspect ratio (= H/R)
    type(spherical_planetesimal_t), pointer :: planet => null() ! pointer to the planetesimal
    type(spherical_planetesimal_t), pointer :: planetesimal => null() ! pointer to the planetesimal
  end type
  type(experiment_parameters_t),  public :: experiment_params

PUBLIC Keplerian_speed, MMSN_surface_density, MMSN_temperature
CONTAINS

!> =============================================================================
!> Calculate the Keplerian velocity for a given cylindrical radius.
!> `cylr` should be in units of AU.
!> The product GM should have units of `scaling%l**3 / scaling%t**2`.
!> =============================================================================
FUNCTION Keplerian_speed (cylr) RESULT (vkep)
  real(8), intent(in) :: cylr
  real(8) :: vkep
  !----------------------------------------------------------------------------
  vkep = sqrt(experiment_params%mstar * cgs%m_sun / scaling%m * scaling%grav / cylr)
END FUNCTION Keplerian_speed

!> =============================================================================
!> Calculate the surface density for the MMSN model, given a cylindrical radius.
!> Result will have units of `scaling%m / scaling%l**2`.
!> `cylr` should be in units of AU.
!> =============================================================================
FUNCTION MMSN_surface_density (cylr) RESULT (sigma)
  real(8), intent(in) :: cylr
  real(8) :: sigma
  !----------------------------------------------------------------------------
  ! the coefficient has units of g cm**-2, but is then converted to code units.
  sigma = 1700.0d0 * cylr**(-1.5d0) / scaling%d
END FUNCTION MMSN_surface_density

!> =============================================================================
!> Calculate the locally (i.e., vertically) isothermal temperature for the MMSN
!> model, given a cylindrical radius.
!> `cylr` should be in units of AU.
!> Result will have units of `scaling%temp`.
!> =============================================================================
FUNCTION MMSN_temperature (cylr) RESULT (temp)
  real(8), intent(in) :: cylr
  real(8) :: temp
  !----------------------------------------------------------------------------
  ! the coefficient has units of K, but is then converted to code units.
  temp = 300.0d0 * cylr**(-0.5d0) / scaling%temp
END FUNCTION MMSN_temperature

END MODULE parameters_mod
