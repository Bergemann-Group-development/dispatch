!===============================================================================
!> Set up initial conditions for the Truelove grav. collapse problem.
!===============================================================================
MODULE initial_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE index_mod
  USE patch_mod
  USE scaling_mod
  implicit none
  private

  type, public:: initial_t
    character(len=64):: solver
    logical:: mhd
    real:: gamma
  contains
    procedure:: init
    procedure:: condition
  end type

  real(8), parameter:: pi = asin(1.0) * 2.0
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self, solver, gamma)
  class(initial_t):: self
  character(len=64):: solver
  real, optional:: gamma
  !----------------------------------------------------------------------------
  self%solver = solver
  if (present(gamma)) then
    self%gamma = max(gamma, 1.0001)
  else
    self%gamma = 1.0001
  end if
END SUBROUTINE init

!===============================================================================
!> Set up a uniform cloud configurations as described in Truelove et al. (1998,
!> ApJ, 495, 821) for self-gravitational collapse.
!> This set up assumes pure or near isothermality (eg., gamma = 1.001). Lastly,
!> magnetic fields have not been implemented in any form.
!>
!> This set up was ported from AZEuS; writing credits go to Jean-Pierre De Villiers
!> and Jon Ramsey.
!>
!>  INPUT/NAMELIST VARIABLES:
!>    rho0         density of uniform cloud           (default =  1.00)
!>    r0           cloud radius                       (default =  0.25)
!>    truealph     energy ratio: thermal    to grav.  (default =  0.25)
!>    truebeta     energy ratio: rotational to grav.  (default =  0.00)
!>    chipres      pressure ratio: cloud to ambient   (default =  1.00)
!>    chirho       density ratio:  cloud to ambient   (default =  1.00)
!>    perta        amplitude of initial perturbation  (default =  0.00)
!>    pertm        order (m) of initial perturbation  (default =  2   )
!>    x10          1-coordinate of cloud center       (default =  0.50)
!>    x20          2-coordinate of cloud center       (default =  0.50)
!>    x30          3-coordinate of cloud center       (default =  0.50)
!>    gamma        ratio of specific heats            (default = 1.0001) !FIXME
!>
!>  LOCAL VARIABLES:
!>  i) Cloud structural properties:
!>    vol          volume
!>    mass         mass
!>    momi         moment of inertia (about 3-axis)
!>    omega        initial angular velocity (rigidly rotating)
!>    pcloud       the cloud pressure
!>  ii) Cloud energy variables:
!>    egrav        gravitational energy
!>    etherm       thermal energy
!>    erot         rotational energy
!>  iii) Ambient properties:
!>    pamb         pressure
!>    rhoamb       density
!>  iv) Jeans-related variables:
!>    tff          free-fall time (for entire cloud)
!>    cjlen        initial Jeans length   (for entire cloud)
!>    cjnum        initial Jeans number   (for entire cloud)
!>  v) Miscellaneous variables:
!>    rz           distance of zone from cloud center
!>    kappa        constant of proportionality in the polytropic eqn.
!>    tff          free-fall time
!>    vff          free-fall speed
!>    csqiso       square of the isothermal sound speed
!>    phi0         azimuthal angle for current zone (rel. to cloud ctr.)
!>    gammaad      ratio of specific heats for an ideal gas (gamma = 5/3)
!>    c1           matching constant for cloud potential
!>    c2           matching constant for ambient potential
!>    c3           matching constant for ambient potential
!>
!===============================================================================
SUBROUTINE condition (self, m, f, idx)
  class(initial_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  class(mesh_t), pointer:: m(:), m1, m2, m3
  type(index_t):: idx
  integer:: i, j, k
  integer:: pertm=2
  real(8):: rho0=15.27026, r0=0.25, truealpha=0.0475, truebeta=0.0, perta=0.0, &
            chipres=1.0, chirho=100.0, x0=0.0, y0=0.0, z0=0.0
  real(8):: phi0    , q1      , kappa   , c1      , c2      , c3       &
          , pcloud  , vol     , mass    , momi    , omega   , gammad   &
          , egrav   , etherm  , erot    , pamb    , rhoamb  , rz       &
          , rzsq    , tff     , cjlen   , cjnum   , csqiso  , vff      &
          , nrg     , x       , y       , z       , rzxy    , rzxm1    &
          , rzym1
  real(8), parameter:: gammaad = 5.0 / 3.0
  real(8), parameter:: r2 = 1.0
  logical, save:: first_time=.true.
  namelist / experiment_params / rho0, r0, truealpha, truebeta, chipres &
                               , chirho, x0, y0, z0, perta, pertm
  !----------------------------------------------------------------------------
  call trace%begin('initial_t%condition')

  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*, experiment_params)
  end if
  !$omp end critical (input_cr)

!-----------------------------------------------------------------------
!      Calculate uniform cloud parameters.
!-----------------------------------------------------------------------
  rhoamb = rho0 / chirho
  mass   = 4.0 * pi * r0**3 * rho0 / 3.0
  momi   = 0.4 * mass * r0**2
  vol    = 4.0 * pi * r0**3 / 3.0
  omega  = sqrt ( 4.0 * pi * scaling%grav * rho0 * truebeta )
  erot   = 0.5 * momi * omega**2
  csqiso = 4.0 * pi * rho0 * mass**2 / 3.0
  csqiso = csqiso**(1.0/3.0) * 0.4 * truealpha * scaling%grav
  vff    = sqrt ( 2.0 * scaling%grav * mass / r0 )
!-----------------------------------------------------------------------
!      Calculate ambient parameters.
!-----------------------------------------------------------------------
  kappa  = csqiso / self%gamma / rho0**(self%gamma - 1.0)
  pcloud = kappa * rho0**self%gamma
  pamb   = pcloud / chipres
  rhoamb = rho0   / chirho
!-----------------------------------------------------------------------
!      Calculate Jeans-related parameters
!-----------------------------------------------------------------------
  tff    = sqrt( 3.0 * pi / ( 32.0 * scaling%grav * rho0) )
  cjlen  = tff * sqrt( 32.0 * csqiso / 3.0 )
  cjnum  = maxval(m%d) / cjlen
!-----------------------------------------------------------------------
!      Calculate analytical gravitational potential variables
!-----------------------------------------------------------------------
  c1     = -2.0 * scaling%grav * pi * rhoamb * (r2**2 - r0**2 * (1.0 - chirho))
  c2     = -2.0 * scaling%grav * pi * r2**2 * rhoamb
  c3     = (4.0/3.0) * scaling%grav * pi * r0**3 * (rho0 - rhoamb)
!-----------------------------------------------------------------------
!      Print out diagnostic information
!-----------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    print "('Gravitational constant =',1pe10.3)", scaling%grav
    print "('Cloud mass =',1pe10.3)", mass
    print "('Cloud free-fall time =',1pe10.3)", tff
    print "('Cloud free-fall speed =',1pe10.3)", vff
    print "('Initial isothermal sound speed =',1pe10.3)", sqrt(csqiso)
    print "('Initial Jeans length = ',1pe10.3)", cjlen
    print "('Initial Jeans number = ',1pe10.3)", cjnum
    print "('Initial angular speed =',1pe10.3)", omega
  end if

  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  ! initialise field variables
  do k=m3%lb,m3%ub
    z = m3%centre_nat + m3%r(k)
    do j=m2%lb,m2%ub
      y = m2%centre_nat + m2%r(j)
      do i=m1%lb,m1%ub
        x = m1%centre_nat + m1%r(i)
        rzsq = (x - x0)**2 + (y - y0)**2 + (z - z0)**2 ! zone-centred
        rz = sqrt(rzsq)
        rzxy = sqrt((x - x0)**2 + (y - y0)**2) ! zone-centred

        if (rz <= r0) then
          phi0             = datan2(y - y0, x - x0)
          q1               = 1.0 + perta * cos(real(pertm,kind=8) * phi0)
          nrg              = pcloud / (self%gamma - 1.0)
          f(i,j,k,idx%d)   = q1 * rho0
          !!!f(i,j,k,idx%phi) = c1 + 2.0 * pi * scaling%grav * rho0 * rzsq / 3.0
        else
          nrg              = pamb / (self%gamma - 1.0)
          f(i,j,k,idx%d)   = rhoamb
          !!!f(i,j,k,idx%phi) = c2 - c3 / rz + 2.0 * scaling%grav * pi * rhoamb * rzsq / 3.0
        endif
        f(i,j,k,idx%phi) = 0.0
        ! set momenta for staggered solvers appropriately
        f(i,j,k,idx%px)  = 0.0
        f(i,j,k,idx%py)  = 0.0
        f(i,j,k,idx%pz)  = 0.0
        if (self%solver(1:7) == 'stagger' .or. trim(self%solver) == 'zeus_mhd_patch') then
          rzxm1 = sqrt((x - m1%d - x0)**2 + (y - y0)**2 + (z - z0)**2)
          rzym1 = sqrt((x - x0)**2 + (y - m2%d - y0)**2 + (z - z0)**2)
          if ((rz <= r0 .or. rzxm1 <= r0) .and. i > 1) then
            f(i,j,k,idx%px)  =-omega * (y - y0) * exp(0.5 * (log(f(i,j,k,idx%d)) + log(f(i-1,j,k,idx%d))))
          end if
          if ((rz <= r0 .or. rzym1 <= r0) .and. j > 1) then
            f(i,j,k,idx%py)  = omega * (x - x0) * exp(0.5 * (log(f(i,j,k,idx%d)) + log(f(i,j-1,k,idx%d))))
          end if
        else
          if (rz < r0) then
            f(i,j,k,idx%px)  =-omega * rzxy * sin(phi0) * f(i,j,k,idx%d)
            f(i,j,k,idx%py)  = omega * rzxy * cos(phi0) * f(i,j,k,idx%d)
          end if
        end if

        if (self%solver(1:9) == 'stagger2_') then
          f(i,j,k,idx%s) = (log(nrg * (self%gamma - 1.0)) - log(f(i,j,k,idx%d)) * self%gamma) &
                         * f(i,j,k,idx%d) / (self%gamma - 1.0)
        else if (self%solver(1:10) == 'stagger2e_') then
          f(i,j,k,idx%e) = nrg
        else if (self%solver(1:7) == 'ramses_') then
          f(i,j,k,idx%e) = nrg + 0.5 * (f(i,j,k,idx%px)**2 + f(i,j,k,idx%py)**2) / f(i,j,k,idx%d)
        else if (trim(self%solver) == 'zeus_mhd_patch') then
          f(i,j,k,idx%e) = nrg
        end if

        ! this experiment is hydro and gravity only.
        if (self%mhd) then
          f(i,j,k,idx%bx) = 0.0
          f(i,j,k,idx%by) = 0.0
          f(i,j,k,idx%bz) = 0.0
        end if
      end do
    end do
  end do

  call trace%end
END SUBROUTINE condition

END MODULE initial_mod
