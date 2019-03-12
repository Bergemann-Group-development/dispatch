!===============================================================================
MODULE parameters_mod
  USE io_mod
  USE trace_mod
  implicit none
  private

  type, public :: experiment_parameters_t
    real(4) :: R          = 30.0          ! Radius of corotation in AU
    real(8) :: mstar      = 1.7           ! mass of the star in solar mass
    real(8) :: mu         = 2.34          ! mean molecular weight in units of amu
    real(8) :: sigma0     = 0.51          ! characteristic density for gas in cgs units (g/cm²)
    real(4) :: Rc         = 80.2          ! characteristic radius for gas and dust in AU
    real(8) :: gamma_c    =-0.2           ! exponent for the exponential cut-off model of dust density
    real(4) :: kill_rmin  = 15.0          ! inner and outer radii where one starts to
    real(4) :: kill_rmax  = 60.0          ! apply the "wave killing" boundary conditions in AU
    real(4) :: roffset    = 20.0          ! offset of the inner edge of the domain from the (cylindrical) origin in AU
    real(8) :: dust2gas   = 0.01          ! dust-to-gas mass ratio if pressureless dust is included in the experiment
    real(8) :: T0         = 30.0          ! characteristic temperature for gas in Kelvin
    real(4) :: RT         = 20.0          ! characteristic radius for gas temperature in AU
    real(8) :: gamma_T    = 0.65          ! exponent for the exponential cut-off model of temperature
    real(8) :: nu         = 1.0e-5        ! intensity of viscosity
    real(8) :: a          = 1.0           ! size of dust grains in cgs units (cm)
    real(8) :: d_int      = 3.0           ! internal density of dust grains in cgs units (g/cm³)
    integer :: nv         = 14            ! number of variables.
    integer :: nc         = 2             ! number of fluid components
    logical :: wave_kill  = .False.       ! activate or deactivate wave killing zones
    real(8) :: amp        = 1.0e-3        ! amplitude of initial perturbation
    real(8) :: lambda     = 0.5           ! wavelength of initial perturbation
  end type
  type(experiment_parameters_t),  public :: experiment_params

END MODULE parameters_mod
