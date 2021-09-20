!===============================================================================
!> Set up initial conditions for the Truelove grav. collapse problem.
!===============================================================================
MODULE initial_conditions_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE index_mod
  USE patch_mod
  USE gpatch_mod
  USE scaling_mod
  USE eos_mod
  implicit none
  private
  type, public:: initial_conditions_t
    character(len=64):: solver
    logical:: mhd
    real(8):: gamma, g1
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure:: pre_init
    procedure:: init
  end type
  real(8), parameter:: pi = asin(1.0) * 2.0
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(initial_conditions_t):: self
  class(gpatch_t), pointer:: patch
  !----------------------------------------------------------------------------
  self%patch => patch
  if (eos%gamma == 1d0) then
    self%g1 = 1d0
  else
    self%g1 = eos%gamma-1d0
  end if
END SUBROUTINE pre_init

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
SUBROUTINE init (self)
  class(initial_conditions_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  class(mesh_t), pointer:: m(:), m1, m2, m3
  type(index_t), pointer:: idx
  integer:: i, j, k, iostat
  integer:: pertm=2, verbose=0
  real(8):: rho0=15.27026, omega_rot0, r0=0.25, truealpha=0.0475, truebeta=0.0, perta=0.0, &
            chipres=1.0, chirho=100.0, x0=0.0, y0=0.0, z0=0.0, B0
  real(8):: phi0    , q1      , kappa   , c1      , c2      , c3       &
          , pcloud  , vol     , mass    , momi    , omega   , gammad   &
          , egrav   , etherm  , erot    , pamb    , rhoamb  , rz       &
          , rzsq    , tff     , cjlen   , cjnum   , ciso  , vff      &
          , nrg     , x       , y       , z       , rzxy    , rzxm1    &
          , rzym1   , d
  real(8), parameter:: gammaad = 5.0 / 3.0
  real(8), parameter:: r2 = 1.0
  logical, save:: first_time=.true.
  namelist / experiment_params / rho0, r0, truealpha, truebeta, chipres &
                               , chirho, x0, y0, z0, perta, pertm, ciso, B0, verbose
  !----------------------------------------------------------------------------
  call trace%begin ('initial_conditions_t%init')
  call io%lock%set ('initial_conditions_t%initt')
  m => self%patch%mesh
  f => self%patch%mem(:,:,:,:,1,1)
  idx => self%patch%idx
  !.............................................................................
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
  end if
  call io%lock%unset ('initial_conditions_t%initt')

!-----------------------------------------------------------------------
!      Set up sphere with uniform density following Zhao et al. 2021 
!      https://ui.adsabs.harvard.edu/abs/2021MNRAS.505.5142Z/abstract  
!-----------------------------------------------------------------------
  rho0 = 4.77e-19/scaling%d !0.239819  g cm^-3
  r0 = 5d-1 ! 
  omega_rot0 =  1e-13 * scaling%t ! 0.27446s^-1
  !vphi0 = 2e4 ! cm s^-1 ! Li et al. 2014
  !vphi = vphi0 * (exp(2*R/R_rot0) - 1d0) / (exp(2*R/R_rot0) + 1d0) ! vphi=vphi0*tanh(R/R0) ! Li et al. 2014
  B0 = 10.6e-6 / scaling%B  ! 0.20628526774 microG
!-----------------------------------------------------------------------
!      Calculate ambient parameters.
!-----------------------------------------------------------------------
  pcloud = ciso**2 * rho0*eos%gamma
  pamb   = pcloud / chipres
  rhoamb = rho0   / chirho
!-----------------------------------------------------------------------
!      Calculate Jeans-related parameters
!-----------------------------------------------------------------------
  tff    = sqrt( 3.0 * pi / ( 32.0 * scaling%grav * rho0) )
  cjlen  = tff * sqrt( 32.0 * ciso / 3.0 )
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
    call io%header ('initial_conditions:')
    print "(' Gravitational constant =',1pe10.3)", scaling%grav
    !print "('             Cloud mass =',1pe10.3)", mass
    print "('   Cloud free-fall time =',1pe10.3)", tff
    !print "('  Cloud free-fall speed =',1pe10.3)", vff
    print "(' Isothermal sound speed =',1pe10.3)", ciso
    print "('   Initial Jeans length =',1pe10.3)", cjlen
    print "('   Initial Jeans number =',1pe10.3)", cjnum
    !print "('  Initial angular speed =',1pe10.3)", omega
    print "('              Gas gamma =',1pe10.3)", eos%gamma
  end if
  !
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  print * ,'x0,y0,z0', x0,y0,z0
  ! initialise field variables
  do k=m3%lb,m3%ub
    z = (m3%centre_nat + m3%r(k))
    print * ,'z', z
    do j=m2%lb,m2%ub
      y = (m2%centre_nat + m2%r(j))
      do i=m1%lb,m1%ub
        x = (m1%centre_nat + m1%r(i))
        rzsq = (x - x0)**2 + (y - y0)**2 + (z - z0)**2 ! zone-centred
        rz = sqrt(rzsq)
        rzxy = sqrt((x - x0)**2 + (y - y0)**2) ! zone-centred
        !
        if (rz <= r0) then
          phi0             = datan2(y - y0, x - x0)
          q1               = 1.0 + perta * cos(real(pertm,kind=8) * phi0)
          nrg              = pcloud / self%g1
          d                = rho0
        else
          nrg              = pamb /self%g1
          d                = rhoamb
        endif
        f(i,j,k,idx%d)   = d
        !print*,'f(i,j,k,idx%d)',f(i,j,k,idx%d)
        if (idx%phi > 0) f(i,j,k,idx%phi) = 0.0
        ! set momenta for staggered solvers appropriately
        f(i,j,k,idx%px)  = 0.0
        f(i,j,k,idx%py)  = 0.0
        f(i,j,k,idx%pz)  = 0.0
        if (rz < r0) then
          f(i,j,k,idx%px)  = -omega_rot0 * rzxy * sin(phi0) * f(i,j,k,idx%d)
          f(i,j,k,idx%py)  =  omega_rot0 * rzxy * cos(phi0) * f(i,j,k,idx%d)
        end if
        !
        nrg = nrg + 0.5 * (f(i,j,k,idx%px)**2 + f(i,j,k,idx%py)**2)/d
        ! B-field
        if (idx%bx > 0) then
          f(i,j,k,idx%bx) = 0.0
          f(i,j,k,idx%by) = 0.0
          f(i,j,k,idx%bz) = B0
          nrg = nrg + 0.5 * (f(i,j,k,idx%bx)**2 + f(i,j,k,idx%by)**2 + f(i,j,k,idx%bz)**2)
        end if
        ! internal energy or entropy
        if (idx%e > 0) then
          f(i,j,k,idx%e) = nrg
        end if
        if (idx%s > 0) then
          f(i,j,k,idx%s) = (log(nrg*self%g1) - log(d)*eos%gamma)*d/self%g1
        end if
        !.......................................................................
        if (verbose > 0) then
          if (i==j .and. j==k) write (stdlog,'("IC:",i4,6f12.6)') &
            i, d, f(i,j,k,idx%px:idx%pz), f(i,j,k,idx%e), nrg
        end if
        !.......................................................................
      end do
    end do
  end do

  call trace%end()
END SUBROUTINE init

END MODULE initial_conditions_mod
