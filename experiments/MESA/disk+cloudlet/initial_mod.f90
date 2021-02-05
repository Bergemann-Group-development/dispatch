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
    real(8):: gamma, g1
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
  real(8), optional:: gamma
  !----------------------------------------------------------------------------
  self%solver = solver
  if (present(gamma)) then
    self%gamma = gamma
  else
    self%gamma = 1d0
  end if
  if (self%gamma == 1d0) then
    self%g1 = 1d0
  else
    self%g1 = self%gamma-1d0
  end if
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE condition (self, m, f, idx)
  class(initial_t):: self
  real(kind=KindScalarVar), dimension(:,:,:,:), pointer:: f
  class(mesh_t), pointer:: m(:), m1, m2, m3
  type(index_t):: idx
  integer:: i, j, k, iostat
  integer:: verbose=0
  real(8):: rhobg, r0, cs, cs0, Sigma_AU, pamb, zau  &
        ,   x, y, z, rzsq, rz, rzxy, vKep, Omeg, temp, Tdisk, udisk  &
        ,   theta, vx, vy, vz, nrg, H, Sigma0, rhoeq, rho_disk, d  &
        ,   kconst, ddec, Sigma0in, rhoeq_in, rho_disk_in
  real(8):: Msun_in_g=1.991e33
  real(8), parameter:: gammaad = 5.0 / 3.0
  real(8), parameter:: r2 = 1.0
  real(8), parameter:: G_grav = 6.6743e-8
  real(8), parameter:: Mstar = 2.5
  real(8), parameter:: au_in_cm = 1.496e+13
  real(8), parameter:: mmsn_fac = 1d0
  real(8), parameter:: alpha_Sigma = 15d-1
  real(8), parameter:: rlim_in = 8d1 * au_in_cm
  logical, save:: first_time=.true.
  namelist / experiment_params / r0, verbose
  !----------------------------------------------------------------------------
  call trace%begin ('initial_t%condition')
  call io%lock%set ('initial_t%condition')
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
  end if
  call io%lock%unset ('initial_t%condition')
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!      Print out diagnostic information
!-----------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    print "('Gravitational constant =',1pe10.3)", scaling%grav
    print "('             Gas gamma =',1pe10.3)", self%gamma
  end if
  !
  m1 => m(1)
  m2 => m(2)
  m3 => m(3)
  rhobg = 1e-22
  r0 = 300.*au_in_cm
  cs = ( scaling%k_B * 10d0 / (2.3 * scaling%m_h) )**0.5 ! 2.3 is mean molecular weight
  pamb = rhobg*cs**2


  Sigma_AU = mmsn_fac*17d2 ! surface density at 1 au; units of g/cm^2

  cs0 = 0.05 * 3e6 ! sound speed in cm for aspect ratio of 0.05 at 1 AU

  ! initialise field variables
  do k=m3%lb,m3%ub
    z = m3%centre_nat + m3%r(k)
    do j=m2%lb,m2%ub
      y = m2%centre_nat + m2%r(j)
      do i=m1%lb,m1%ub
        x = m1%centre_nat + m1%r(i)

        zau  = z*scaling%l
        rzsq = x**2 + y**2 + z**2 ! zone-centred
        rz = sqrt(rzsq) * scaling%l
        rzxy = sqrt(x**2 + y**2) * scaling%l ! zone-centred
        !
        vKep = (G_grav * Mstar*scaling%m / rzxy)**0.5
        Omeg = (G_grav * Mstar*scaling%m / rzxy**3)**0.5
        !cs  = cs0 * (rzxy/au_in_cm)**(-3./8) ! in cm/s (scaling see lecture notes by Armitage p.16)
        Tdisk = 1d1
        temp = Tdisk * (rzxy/au_in_cm)**(-3d0/4d0)
        udisk = 1d0 / (self%gamma-1d0) * scaling%k_B / (2.3 * scaling%m_h) * temp

        theta = atan(y,x)
        vx = vKep * sin(theta)
        vy = -vKep* cos(theta)
        vz = 0d0
        nrg              = pamb / self%g1 / scaling%e
        H = cs / Omeg
        if (rzxy <= r0) then

          ! assign density
          Sigma0 = Sigma_AU * (rzxy / au_in_cm)**(-alpha_Sigma)
          !Sigma0end = Sigma_AU * (Rdisk/AU)**(-alpha_Sigma)
          rhoeq = Sigma0/( (2d0*pi)**0.5 * H )
          rho_disk = rhoeq * exp(-z**2/ (2d0*H**2d0) )

          d  = rho_disk
          !!!f(i,j,k,idx%phi) = c1 + 2.0 * pi * scaling%grav * rho0 * rzsq / 3.0

          if (rzxy .lt. rlim_in) then
          ! assign density at inner radius
             Sigma0in = Sigma_AU * (rlim_in/au_in_cm)**(-alpha_Sigma)
             rhoeq_in = Sigma0in/( sqrt(2d0*pi)*H )
             rho_disk_in = rhoeq_in * exp(-z**2/ (2d0*H**2) )

             kconst = 1d0
             ddec = 2d0 ! decay interval: at rlim_in-ddec. au the disk falls off by a factor of 50%
             rho_disk = rho_disk_in + (rhobg - rho_disk_in) / (1d0 + exp(-kconst * ( (-rzxy + rlim_in)/au_in_cm - ddec)  ) )
             d = rho_disk
          endif
        else
          d  = rhobg
          !!!f(i,j,k,idx%phi) = c2 - c3 / rz + 2.0 * scaling%grav * pi * rhoamb * rzsq / 3.0
        endif
        f(i,j,k,idx%d)   = d / scaling%d
        ! set momenta for staggered solvers appropriately
        f(i,j,k,idx%px)  = d / scaling%d * (vx / scaling%u)
        f(i,j,k,idx%py)  = d / scaling%d * (vy / scaling%u)
        f(i,j,k,idx%pz)  = d / scaling%d * (vz / scaling%u)
        !
        if (self%solver(1:7) == 'ramses_') then
          nrg = nrg + 0.5 * (f(i,j,k,idx%px)**2 + f(i,j,k,idx%py)**2)/ (d/ scaling%d)
        end if
        ! internal energy or entropy
        if (idx%e > 0) then
          f(i,j,k,idx%e) = nrg
        end if
        if (idx%s > 0) then
          f(i,j,k,idx%s) = (log(nrg*self%g1) - log(d)*self%gamma)*(d/scaling%d) /self%g1
        end if
        ! this experiment is hydro and gravity only.
        if (idx%bx > 0) then
          f(i,j,k,idx%bx) = 0.0
          f(i,j,k,idx%by) = 0.0
          f(i,j,k,idx%bz) = 0.0
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
END SUBROUTINE condition

END MODULE initial_mod
