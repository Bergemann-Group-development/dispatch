!===============================================================================
!> Port and test dust evolution from Podolak_1.py
!===============================================================================
MODULE dust_evolve_mod
  USE omp_lib
  USE units_mod
  USE math_mod
  USE trace_mod
  USE io_mod
  USE scaling_mod
  USE particles_mod
  USE omp_timer_mod
  implicit none
  private
  type, public:: dust_evolve_t
    integer:: nbin
    real(KIND=8):: start_time=0d0                           ! start time
    logical:: on=.false.
  contains
    procedure:: init
    procedure:: update
    procedure:: test
  end type

  integer, parameter :: precision_e = selected_real_kind(15, 100)
  integer:: verbose=0
  integer:: n_sub = 10
  integer:: nbin = 1
  logical:: on=.false.
  logical:: do_implicit=.true.
  logical:: do_brownian=.false.
  logical:: do_in_cell_turbulence=.false.
  logical:: do_set_frag_vel =.false.
  logical:: detailed_timer=.false.
  real(KIND=8):: start_time=0d0                             ! start time
  real:: max_fraction1 = 1.0
  real:: max_fraction2 = 0.1
  real:: size_monomer  = 1e-7                               ! Size of monomer used for determining when we have monomer-monomer collision in context of fragmentation probability
  real:: monomer_frag  = 1e5                                ! When monomer-monomer collision causes fragmentation set very high 1km/s
  real:: size_aggregate = 1e4                               ! Aggregate size used for determining aggregate-aggregate collision in context of fragmentation
  real:: aggregate_frag = 1e3                               ! Aggregate-Aggregate collision fragmentation velocity set low 10 m/s
  real(KIND=precision_e):: tt_dust = 10                     ! Dust temperature
  real(KIND=precision_e):: A_br = 2.8e3                     ! Experimental parameter for breaking energy of particles
  real(KIND=precision_e):: Eps_star = 3.7e10                ! reduced elasticity modulus in dyn cm^-2
  real(KIND=precision_e):: zeta_crit = 2e-7                 ! Critical roll in cm
  real(KIND=precision_e):: gamma_sil = 25                   ! surface energy density of silicates, in erg/(cm^-2)
  real(KIND=precision_e):: alpha = 1e-3                     ! Turbulence parameter, used when do_in_cell_turbulence is true for calculation of turbulent velocities
  real(KIND=precision_e):: sigma_d = 0.18                   ! Column density of the dust, used when do_in_cell_turbulence is true for calculation of turbulent velocities
  real(KIND=precision_e):: xi = 1.83                        ! The exponent used for the distribution in the case of cratering
  real(KIND=precision_e):: trans_width = 0.1                ! Transition width for fragmentation probability in m/s
  real(KIND=precision_e):: set_frag_vel = 1e2               ! If do_set_frag_vel then this is the set fragmentation velocity for all particles.
  real(KIND=precision_e):: dt = 1e-3 * 60*60*24*365         ! Time step in dust_evolve in years
  real(KIND=precision_e):: save_time = 1 * 60*60*24*365
  real(KIND=precision_e):: eps_time = 1.10 * 60*60*24*365
  real(KIND=precision_e), ALLOCATABLE:: C(:,:,:)
  real(KIND=precision_e), ALLOCATABLE:: SS(:,:,:)
  real(KIND=precision_e), ALLOCATABLE:: ms(:,:)
  real(KIND=precision_e), ALLOCATABLE:: ep(:,:)
  real(KIND=precision_e), ALLOCATABLE:: sigma_col(:,:)
  real(KIND=precision_e), ALLOCATABLE:: mn(:,:)
  real(KIND=precision_e), ALLOCATABLE:: mm(:,:)
  real(KIND=precision_e), ALLOCATABLE:: m_cgs(:)
  real(KIND=precision_e), ALLOCATABLE:: a(:)
  real, ALLOCATABLE:: st(:)
  type(dust_evolve_t), public:: dust_evolve
CONTAINS

!===============================================================================
!> Read parameters and initialize -- add more parameters, if they are known
!> in the particle solver
!===============================================================================
SUBROUTINE init (self, ln_bin, lm_cgs, rho_s)
  class(dust_evolve_t) :: self
  integer              :: ln_bin               ! number of bins
  real                 :: lm_cgs(:)            ! particle masses
  real                 :: rho_s               ! particle solid density
  !.............................................................................
  integer              :: iostat              ! missing(<0), syntax errors (>0)
  logical, save        :: first_time=.true.   ! shared memory


  integer:: i, j, k, l, kn, minij, maxij
  real(precision_e):: ratio_r, min_a, max_a, m_crat, eps, Stot
  real(precision_e):: m_big, m_small
  real(precision_e), ALLOCATABLE, DIMENSION(:):: m
 
  !real(precision_e), ALLOCATABLE, DIMENSION(:,:,:):: SS, C
  !--- add parameters here -----------------------------------------------------
  namelist /dust_evolve_params/ on, verbose, do_brownian, do_in_cell_turbulence, &
    tt_dust, sigma_d, alpha, xi, trans_width, set_frag_vel, Eps_star, do_set_frag_vel, &
    dt, do_implicit, start_time, gamma_sil, size_monomer, monomer_frag,size_aggregate, aggregate_frag,    &
    max_fraction1, max_fraction2, n_sub, detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%init')
  !-----------------------------------------------------------------------------
  ! Make sure that only one thread reads the namelist
  !-----------------------------------------------------------------------------
  call io%lock%set ('dust_evolve_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_evolve_params, iostat=iostat)
    if (iostat > 0) then                       ! if syntax errors
      rewind (io%input)                        ! rewind to top of file
      read (io%input, dust_evolve_params)      ! read againg to get more info
    end if
    call io%namelist_warning ('dust_evolve_params', iostat)
    write (io%namelist, dust_evolve_params)    ! write to data/run/params.nml

    first_time = .false.
    if (.not.on) then
      call io%lock%unset ('dust_evolve_t%init')
      call trace%end ()
      return
    end if

    call cgs%init
    nbin = ln_bin
    ALLOCATE(SS(nbin,nbin,nbin), C(nbin,nbin,nbin), ms(nbin,nbin), ep(nbin,nbin), a(nbin), st(nbin), m(nbin))
    ALLOCATE(sigma_col(nbin,nbin), mn(nbin,nbin), mm(nbin,nbin))

    allocate(m_cgs(nbin))
    m_cgs = lm_cgs
    !call self%test                             ! run a self%test once

    m = m_cgs
    a = (3.0/4.0 * m/(rho_s * math%pi) )**(1.0/3.0)       ! Initial setup of bin distribution
    st = rho_s * a * math%pi/(2*sigma_d*100)


    C = 0
    SS = 0
    ms = 0
    ep = 0
    sigma_col = 0
    mn = 0
    mm = 0
    do i=1,nbin
      do j=1,nbin                                       ! Setup of matrices
        minij = min(i,j)
        maxij = max(i,j)

        ms(i,j) = m(i) + m(j)

        do kn = nbin,1,-1
          if ( m(kn)<ms(i,j) ) exit
        end do

        ! ubm(i,j) = sqrt(8 * cgs%k_b* tt * (m(i) + m(j) )/(math%pi*m(i)*m(j) ))
        ! utm(i,j) = max(utm_ini(i), utm_ini(j))
        ! du(i,j) = ubm(i,j) + utm(i,j)
        sigma_col(i,j) = math%pi * (a(i) + a(j))**2

        if (m(i)/m(j) > 10 .or. m(i)/m(j) < 0.1 ) then        !!!
          !Cratering
          SS(1:minij,i,j) = a(1:minij) ** (-3*xi)
          Stot = SUM(m * SS(:,i,j) ) / (m(minij) *2)
          SS(1:minij,i,j) = SS(1:minij,i,j) / Stot
          m_crat =  abs(m(i)-m(j))
          m_big = max(m(i),m(j))
          m_small = min(m(i),m(j))

          do l = nbin,1,-1
            if (m(l)-m_big<-m_small) exit
          end do


          eps = (m_crat - m(l))/(m(l+1) - m(l))
          SS(l,i,j) = 1-eps
          SS(l+1,i,j) = eps

        else
          !Fragmentation
          SS(1:maxij,i,j) = a(1:maxij) ** (-3*xi)
          Stot = SUM(m * SS(:,i,j) ) / (m(i) + m(j))
          SS(1:maxij,i,j) = SS(1:maxij,i,j) / Stot
        end if


        if (kn+1 <= nbin) then
          mm(i,j) = m(kn)
          mn(i,j) = m(kn+1)
        else
          sigma_col(i,j) = 0
        end if

      end do
    end do

    
    !pc = 1-pf
    where (sigma_col == 0) ep = 0
    where (sigma_col /= 0) ep = (mn - ms) / (mn - mm)

    do k=1,nbin
      where (mm == m(k) ) C(k,:,:) = ep
      where (mn == m(k) ) C(k,:,:) = 1 - ep
    end do

  end if
  call io%lock%unset ('dust_evolve_t%init')

  self%nbin = nbin
  self%start_time = start_time
  self%on = on                                 ! task private copy

  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the dust state for one cell, with _in and _out arguments
!> holding particle mass density dd(n_bin) and velocity vd(n_bin,3)
!===============================================================================
SUBROUTINE update (self, dt_dispatch, vrms_in, dd_in, vd_in, dd_out, vd_out, ds)
  class(dust_evolve_t):: self
  real(kind=8), intent(in)  :: dt_dispatch, ds, vrms_in(:)
  real,         intent(in)  :: dd_in(:), vd_in(:,:)
  real,         intent(out) :: dd_out(:), vd_out(:,:)

  real:: tt, mu, cs, r, msolar, Omega_K, rho_s, t, t_end

  integer:: i, j, k, l, kn, itime
  real(precision_e):: ratio_r, min_a, max_a, m_crat, eps, Stot, a_mu, m_mu, u_f
  real(precision_e):: m_big, m_small
  real(precision_e), ALLOCATABLE, DIMENSION(:):: n, n_old, dn, dpx, dpy, dpz, dd, &
    mv, px, py, pz, v, Sppx, Sppy, Sppz, frac_dn, vrms
  real(precision_e), ALLOCATABLE, DIMENSION(:,:):: nn, kk, &
    ll, du, ubm, utm, vd, pc, pf, identity
  real, ALLOCATABLE, DIMENSION(:):: utm_ini

  ! CGS units below
  real(precision_e) :: tend
  real(8):: start, used, time(5)
  integer, save:: itimer=0, ntimer=10
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%update', itimer=itimer)
  start = wallclock()

  nbin = self%nbin
  if (verbose > 0) &
    print *, 'nbin', nbin
  ALLOCATE(n(nbin), dn(nbin), n_old(nbin), frac_dn(nbin))
  ALLOCATE(vd(nbin,3), vrms(nbin))
  ALLOCATE(utm_ini(nbin))
  ALLOCATE(du(nbin,nbin), pc(nbin,nbin), pf(nbin,nbin), identity(nbin,nbin))
  ALLOCATE(ubm(nbin,nbin), utm(nbin,nbin))
  ALLOCATE(kk(nbin,nbin))
  ALLOCATE(ll(nbin,nbin), nn(nbin,nbin))
  ALLOCATE(dpx(nbin), dpy(nbin), dpz(nbin))
  ALLOCATE(mv(nbin), px(nbin), py(nbin), pz(nbin), v(nbin))
  ALLOCATE(Sppx(nbin), Sppy(nbin), Sppz(nbin))
  time(1) = wallclock()-start

  ! convert from code units to cgs
  tend = scaling%t * dt_dispatch
  n = scaling%d * dd_in / m_cgs
  vd = scaling%u * vd_in * ds
  vrms = scaling%u * vrms_in * ds

  n_old = n

  ! define an identity matrix
  ForAll(i = 1:nbin, j = 1:nbin) identity(i,j) = (i/j)*(j/i)

  px = n * vd(:,1)
  py = n * vd(:,2)
  pz = n * vd(:,3)

  if (verbose > 1) then
    print *, "Starting density"
    print '(5ES13.3)', n
    print "(/)"

    print *, "Starting velocity x:"
    print '(5ES13.3)', vd(:,1)
    print "(/)"

    print *, "Starting velocity y:"
    print '(5ES13.3)', vd(:,2)
    print "(/)"

    print *, "Starting velocity z:"
    print '(5ES13.3)', vd(:,3)
    print "(/)"

    print *, "Starting vrms:"
    print '(5ES13.3)', vrms
    print "(/)"
  end if



  ! compute velocity difference between different size bins
  ! 1) use vd for bulk difference:
  !        vdiff(i,j) = sqrt((vd(i,1) - vd(j,1))**2 + (vd(i,2) - vd(j,2))**2 + (vd(i,3) - vd(j,3))**2)
  ! 2) add vrms in square: vdiff(i,j) += sqrt(vrms(i)**2 + vrms(j)**2)
  ! 3) possibly add brownian motion
  ! 4) possibly add a sub-cell turbulent motion contribution

  ! move forward one timestep
  ! call dust_evolve_step(tend,nd,vd,vdiff)

  ! du = vdiff + brown + turb bliver regnet i init

  if (do_in_cell_turbulence) then
    cs = sqrt(tt_dust*cgs%k_B/(cgs%m_u*scaling%mu))
    utm_ini = cs * sqrt(2*alpha*st)                     ! Handle velocities every update
    where (st > 1.0) utm_ini = cs * sqrt(2*alpha/ st)
  end if
  time(2) = wallclock()-start

  do i=1,nbin
    do j=1,nbin
      ! Velocity differences for all possible collision
      du(i,j) = sqrt((vd(i,1) - vd(j,1))**2 + (vd(i,2) - vd(j,2))**2 + (vd(i,3) - vd(j,3))**2)
      du(i,j) = du(i,j) + sqrt(vrms(i)**2 + vrms(j)**2)


      if (do_brownian) then
        ubm(i,j) = sqrt(8 * cgs%k_b* tt_dust * (m_cgs(i) + m_cgs(j) )/(math%pi*m_cgs(i)*m_cgs(j) ))
        du(i,j) = du(i,j) + ubm(i,j)
      end if

      if (do_in_cell_turbulence) then
        utm(i,j) = max(utm_ini(i), utm_ini(j))
        du(i,j) = du(i,j) + utm(i,j)
      end if


      ! Fragmentation velocity u_f is calculated based on given parameters
      ! the parameters are set in input file to decide fragmentation velocities for different sizes
      a_mu = a(i)*a(j)/(a(i)+a(j))
      m_mu = m_cgs(i)*m_cgs(j)/(m_cgs(i)+m_cgs(j))

      if (do_set_frag_vel) then
        u_f = set_frag_vel
      else
        if ((a(i) < size_monomer ) .and. (a(j) < size_monomer )) then
          u_f = monomer_frag
        else if ((a(i) > size_aggregate ) .and. (a(j) > size_aggregate )) then
          u_f = aggregate_frag
        else
          !u_f  = sqrt(2*A_br) * (gamma_sil**(5.0/6.0) * a_mu**(2.0/3.0)) / (m_mu**(1.0/2.0) * Eps_star**(1.0/3.0))             ! Monomer
          u_f  = sqrt(60 * math%pi**2 * zeta_crit * gamma_sil * a_mu /m_mu )                                                    ! aggregates 5E_roll

        end if
      end if


      ! Decide by the fragmentation velocity is the probability of fragmentation p_f
      ! and by extension the probability for coagulation p_c
      if (du(i,j) < u_f - trans_width) then
        pf(i,j) = 0
      else if (du(i,j) > u_f) then
        pf(i,j) = 1
      else
        pf(i,j) = 1-(u_f-du(i,j))/trans_width
      end if

    end do
  end do
  time(3) = wallclock()-start

  if (verbose>2) then
    print *, "pf"
    print '(5ES13.1)', pf
    print "(/)"
  end if
  pc = 1-pf

  ! Coagulation kernel and fragmentation kernel for Smoluchowski equation
  kk = du * sigma_col * pc
  ll = du * sigma_col * pf

 
  ! Do a explicit step and compare the change in n if too large a change is detected
  ! or the mass is not conserved then do a implicit step instead
  ! if the change is only above the second fraction set then do several explicit steps in smaller time steps.
  call explicit_step(nbin, kk, ll, C, SS, n, px, py, pz, tend, dn, dpx, dpy, dpz)
  frac_dn = 0
  where (n /= 0) frac_dn = ABS(dn/n)
  if (do_implicit .or. ANY( frac_dn > max_fraction1 ) .or. SUM(n*m_cgs)/SUM(n_old*m_cgs) - 1 > 1e-6 .or. ANY(n < 0.0) )  then
    if (verbose > 0) &
      print '(a,(1p,10e10.1))', 'implicit step:', frac_dn
    call implicit_step(nbin, identity, kk, ll, C, SS, n, px, py, pz, tend, dn, dpx, dpy, dpz)
    n = n + dn
    px = px + dpx
    py = py + dpy
    pz = pz + dpz
  else if ( ANY( frac_dn > max_fraction2) ) then
    n_sub = maxval(frac_dn)/max_fraction2 + 1
    if (verbose > 0) &
      print '(a,i3,(1p,10e10.1))', 'substeps:', n_sub, frac_dn
    do i=1,n_sub
      call explicit_step(nbin, kk, ll, C, SS, n, px, py, pz, tend/n_sub, dn, dpx, dpy, dpz)
      n = n + dn
      px = px + dpx
      py = py + dpy
      pz = pz + dpz
    end do
  else
    if (verbose > 0) &
      print '(a,(1p,10e10.1))', 'explicit step:', frac_dn
    ! Keep explicit step


    ! Update the number density and the momentum in each direction
    n = n + dn
    px = px + dpx
    py = py + dpy
    pz = pz + dpz
  end if


  time(4) = wallclock()-start

  if (verbose > 0) then
    print *,"Final n, number density:"
    print '(5ES13.3)', n
    print "(/)"

    print *,"Final dn, number density:"
    print '(5ES13.3)', dn
    print "(/)"
  end if

  if (ANY(n_old>0)) then
    if (verbose > 0) then
      print *,"Mass conservation:", SUM(n*m_cgs)/SUM(n_old*m_cgs) - 1
    else if (abs(SUM(n*m_cgs)/SUM(n_old*m_cgs) - 1) > 1e-6) then
      print *,"Mass conservation:", SUM(n*m_cgs)/SUM(n_old*m_cgs) - 1
    end if
  end if

  if (verbose>1) then
    print "(/)"
    print *,"Final velocity dx:"
    print '(5ES13.3)', dpx/(n+1e-30)
    print "(/)"

    print *,"Final velocity x:"
    print '(5ES13.3)', px/(n+1e-30)
    print "(/)"
  end if

  where (n > 0.0)
    vd(:,1) = px/n
    vd(:,2) = py/n
    vd(:,3) = pz/n
  else where
    vd(:,1) = 0.0
    vd(:,2) = 0.0
    vd(:,3) = 0.0
  end where

  ! convert back to code units
  dd_out = m_cgs * n / scaling%d
  vd_out = vd / (scaling%u * ds)


  if (ntimer > 0) then
    ntimer = ntimer-1
    time(5) = wallclock()-start
    print '(a,6p,5f10.3)', 'times[mus]:', time/nbin
  end if
  used = wallclock() - start

  call trace%end (itimer)
contains

subroutine explicit_step (nbin, kk, ll, C, SS, n, px, py, pz, dt, dn, &
                          dpx, dpy, dpz)
  integer, INTENT(IN)                                      :: nbin
  real(precision_e), INTENT(IN)                            :: dt
  real(precision_e), INTENT(IN), DIMENSION(nbin)           :: n, px, py, pz
  real(precision_e), INTENT(IN), DIMENSION(nbin,nbin)      :: kk, ll
  real(precision_e), INTENT(IN), DIMENSION(nbin,nbin,nbin) :: C, SS

  real(precision_e), intent(out), DIMENSION(:)             :: dn, dpx, dpy, dpz
  
  real(precision_e), ALLOCATABLE, DIMENSION(:)             :: S, Spx, Spy, Spz, Kn_k, Ln_k, mpx, mpy, mpz
  real(precision_e), ALLOCATABLE, DIMENSION(:,:)           :: KC_k, LS_k, Q, Q_frag, pnx, pny, pnz
  integer, save:: itimer=0


  ALLOCATE(mpx(nbin), mpy(nbin), mpz(nbin), pnx(nbin,nbin), pny(nbin,nbin), pnz(nbin,nbin))
  ALLOCATE(S(nbin), Spx(nbin), Spy(nbin), Spz(nbin), Kn_k(nbin), Ln_k(nbin) )
  ALLOCATE(KC_k(nbin,nbin), LS_k(nbin,nbin), Q(nbin,nbin), Q_frag(nbin,nbin))

  if (detailed_timer) call trace%begin ('explicit_step', itimer=itimer)
  mpx = m_cgs * px
  mpy = m_cgs * py
  mpz = m_cgs * pz                                                  ! xyz
  do i=1,nbin
    do j=1,nbin
      nn(i,j) = n(i)*n(j)
      pnx(i,j) = (mpx(i) * n(j) + n(i) * mpx(j))/ms(i,j)
      pny(i,j) = (mpy(i) * n(j) + n(i) * mpy(j))/ms(i,j)
      pnz(i,j) = (mpz(i) * n(j) + n(i) * mpz(j))/ms(i,j)
    end do
  end do

  dn = 0
  dpx = 0
  dpy = 0
  dpz = 0
  S = 0
  Spx = 0
  Spy = 0
  Spz = 0

  Q = kk * nn
  Q_frag = ll * nn

  do k=1,nbin
    KC_k = kk * C(k,:,:)
    LS_k = ll * SS(k,:,:)

    Kn_k = kk(k,:) * n
    Ln_k = ll(k,:) * n

    S(k) = ( 0.5 * SUM(Q*C(k,:,:)) - SUM(Q(k,:)) ) + ( 0.5 * SUM( Q_frag*SS(k,:,:)) - SUM(Q_frag(k,:)) )
    Spx(k) = 0.5 * SUM(KC_k * pnx) - SUM(Kn_k * px(k)) + 0.5 * SUM(LS_k * pnx) - SUM(Ln_k * px(k))
    Spy(k) = 0.5 * SUM(KC_k * pny) - SUM(Kn_k * py(k)) + 0.5 * SUM(LS_k * pny) - SUM(Ln_k * py(k))
    Spz(k) = 0.5 * SUM(KC_k * pnz) - SUM(Kn_k * pz(k)) + 0.5 * SUM(LS_k * pnz) - SUM(Ln_k * pz(k))
  end do

  dn = S * dt
  dpx = Spx * dt
  dpy = Spy * dt
  dpz = Spz * dt
  if (detailed_timer) call trace%end (itimer)
end subroutine explicit_step

subroutine implicit_step (nbin, identity, kk,ll,C,SS,n,px,py,pz, dt, dn, &
                          dpx, dpy, dpz)
  implicit none
  integer, INTENT(IN)                                      :: nbin
  real(precision_e), INTENT(IN)                            :: dt
  real(precision_e), INTENT(IN), DIMENSION(nbin)           :: n, px, py, pz
  real(precision_e), INTENT(IN), DIMENSION(nbin,nbin)      :: kk, ll, identity
  real(precision_e), INTENT(IN), DIMENSION(nbin,nbin,nbin) :: C, SS

  real(precision_e), intent(out), DIMENSION(nbin)          :: dn, dpx, dpy, dpz

  INTEGER                                                  :: js_n, js_nrhs, js_lda, js_ldb
  real(precision_e), ALLOCATABLE, DIMENSION(:)             :: S, Spx, Spy, Spz, Kn_k, Ln_k, js_ipiv, js_ipivx, js_ipivy, js_ipivz, js_info, mpx, mpy, mpz
  real(precision_e), ALLOCATABLE, DIMENSION(:,:)           :: KC_k, LS_k, s1, s2, Jp, Jmx, Jmy, Jmz , I_dt
  real(precision_e), ALLOCATABLE, DIMENSION(:,:)           :: JJ, IJJ, IJpx, IJpy, IJpz, Q, Q_frag, pnx, pny, pnz
  integer, save:: itimer=0

  ALLOCATE(mpx(nbin), mpy(nbin), mpz(nbin), pnx(nbin,nbin), pny(nbin,nbin), pnz(nbin,nbin), I_dt(nbin,nbin))
  ALLOCATE(S(nbin), Spx(nbin), Spy(nbin), Spz(nbin), Kn_k(nbin), Ln_k(nbin) )
  ALLOCATE(KC_k(nbin,nbin), LS_k(nbin,nbin), Q(nbin,nbin), Q_frag(nbin,nbin))

  ALLOCATE(s1(nbin,nbin), s2(nbin,nbin), Jp(nbin,nbin), Jmx(nbin,nbin), Jmy(nbin,nbin), Jmz(nbin,nbin) )
  ALLOCATE(JJ(nbin,nbin), IJJ(nbin,nbin), IJpx(nbin,nbin), IJpy(nbin,nbin), IJpz(nbin,nbin))
  ALLOCATE(js_ipiv(nbin), js_ipivx(nbin), js_ipivy(nbin), js_ipivz(nbin), js_info(nbin))
  if (detailed_timer) call trace%begin ('implicit_step', itimer=itimer)

  mpx = m_cgs * px
  mpy = m_cgs * py
  mpz = m_cgs * pz                                                  ! xyz
  do i=1,nbin
    do j=1,nbin
      nn(i,j) = n(i)*n(j)
      pnx(i,j) = (mpx(i) * n(j) + n(i) * mpx(j))/ms(i,j)
      pny(i,j) = (mpy(i) * n(j) + n(i) * mpy(j))/ms(i,j)
      pnz(i,j) = (mpz(i) * n(j) + n(i) * mpz(j))/ms(i,j)
    end do
  end do

  dn = 0
  dpx = 0
  dpy = 0
  dpz = 0
  S = 0
  Spx = 0
  Spy = 0
  Spz = 0
  Jp = 0
  Jmx = 0
  Jmy = 0
  Jmz = 0

  JJ = 0
  IJJ = 0
  IJpx = 0
  IJpy = 0
  IJpz = 0
  Q = kk * nn
  Q_frag = ll * nn

  do k=1,nbin
    KC_k = kk * C(k,:,:)
    LS_k = ll * SS(k,:,:)

    Kn_k = kk(k,:) * n
    Ln_k = ll(k,:) * n

    S(k) = ( 0.5 * SUM(Q*C(k,:,:)) - SUM(Q(k,:)) ) + ( 0.5 * SUM( Q_frag*SS(k,:,:)) - SUM(Q_frag(k,:)) )
    Spx(k) = 0.5 * SUM(KC_k * pnx) - SUM(Kn_k * px(k)) + 0.5 * SUM(LS_k * pnx) - SUM(Ln_k * px(k))
    Spy(k) = 0.5 * SUM(KC_k * pny) - SUM(Kn_k * py(k)) + 0.5 * SUM(LS_k * pny) - SUM(Ln_k * py(k))
    Spz(k) = 0.5 * SUM(KC_k * pnz) - SUM(Kn_k * pz(k)) + 0.5 * SUM(LS_k * pnz) - SUM(Ln_k * pz(k))

    s1 = KC_k * 1/ms
    s2 = LS_k * 1/ms
    Jp(k,:) = (matmul(s1, n) + matmul(s2, n)) * m_cgs
    Jp(k,k) = Jp(k,k) - (SUM(Kn_k) + SUM(Ln_k) )
    Jmx(k,:) = (matmul(s1, px*m_cgs) + matmul(s2, px*m_cgs)) - px(k) * kk(k,:) - px(k) * ll(k,:)
    Jmy(k,:) = (matmul(s1, py*m_cgs) + matmul(s2, py*m_cgs)) - py(k) * kk(k,:) - py(k) * ll(k,:)
    Jmz(k,:) = (matmul(s1, pz*m_cgs) + matmul(s2, pz*m_cgs)) - pz(k) * kk(k,:) - pz(k) * ll(k,:)
    JJ(k,:) = (matmul(KC_k , n) - n(k) * kk(k,:) ) + (matmul(LS_k,  n) - n(k) * ll(k,:) )
    JJ(k,k) = JJ(k,k) - (SUM(Kn_k) + SUM(Ln_k) )
  end do

  do i=1,nbin
    JJ(i,i)=0.0
    JJ(i,i) = JJ(i,i) - SUM(m_cgs*JJ(:,i))/m_cgs(i)
  end do

  I_dt = identity / dt

  IJJ = I_dt - JJ
  IJpx = I_dt - Jp
  IJpy = I_dt - Jp
  IJpz = I_dt - Jp


  js_n = nbin
  js_nrhs = 1
  js_lda = nbin
  js_ldb = nbin

  ! Solve the inverse matrix
  CALL DGESV( js_n, js_nrhs, IJJ, js_lda, js_ipiv, S, js_ldb, js_info )
  Sppx = Spx - matmul(Jmx , S)
  Sppy = Spy - matmul(Jmy , S)
  Sppz = Spz - matmul(Jmz , S)
  CALL DGESV( js_n, js_nrhs, IJpx, js_lda, js_ipivx, Sppx, js_ldb, js_info )
  CALL DGESV( js_n, js_nrhs, IJpy, js_lda, js_ipivy, Sppy, js_ldb, js_info )
  CALL DGESV( js_n, js_nrhs, IJpz, js_lda, js_ipivz, Sppz, js_ldb, js_info )

  ! Update the parameters
  dn = S
  dpx = Sppx
  dpy = Sppy
  dpz = Sppz
  if (detailed_timer) call trace%end (itimer)
end subroutine implicit_step

END SUBROUTINE update

!===============================================================================
!> Self test
!===============================================================================
SUBROUTINE test (self)
  class(dust_evolve_t):: self
  real:: tt, mu, cs, r, msolar, Omega_K, sigma_d, rho_d, rho_s, u_f, t
  real:: trans_width, alpha, Cdt, xi, dt, t_end, save_time, eps_time, Hp
  real:: micron, mps

  integer:: i, j, k, l,  kn, minij, maxij, itime, js_n, js_nrhs, js_lda, &
    js_ldb, js_info, nbin
  integer, ALLOCATABLE, DIMENSION(:):: js_ipiv
  real(precision_e):: ratio_r, min_a, max_a, m_crat, eps, Stot
  real(precision_e):: m_big, m_small, m_tot
  real(precision_e), ALLOCATABLE, DIMENSION(:):: m, n, S, dn, dp, dd, a, Sp, &
    mv, p, v, mp, Ln_k, Kn_k, Spp
  real(precision_e), ALLOCATABLE, DIMENSION(:,:):: mn, mm, ms, nn, Q, Q_frag, jj, kk, &
    ll, KC_k, LS_k, IJJ, sigma_col, ep, identity, pf, pc, du, ubm, utm, &
    Jp, Jm, IJp, pn, s1, s2
  real(precision_e), ALLOCATABLE, DIMENSION(:):: st, utm_ini
  real(precision_e), ALLOCATABLE, DIMENSION(:,:,:):: SS, C                                                      ! Necesarry?
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%test')
  nbin = 5

  ALLOCATE(n(nbin), st(nbin), dn(nbin), S(nbin), a(nbin), mn(nbin,nbin), mm(nbin,nbin), ms(nbin,nbin))
  ALLOCATE(utm_ini(nbin), sigma_col(nbin,nbin), ep(nbin,nbin))
  ALLOCATE(du(nbin,nbin), pf(nbin,nbin), pc(nbin,nbin))
  ALLOCATE(ubm(nbin,nbin), utm(nbin,nbin), Q(nbin,nbin), Q_frag(nbin,nbin))
  ALLOCATE(jj(nbin,nbin), kk(nbin,nbin))
  ALLOCATE(ll(nbin,nbin), nn(nbin,nbin), identity(nbin,nbin), js_ipiv(nbin))
  ALLOCATE(KC_k(nbin,nbin), LS_k(nbin,nbin) , IJJ(nbin,nbin))
  ALLOCATE(Jp(nbin,nbin), Jm(nbin,nbin), Sp(nbin), dp(nbin))
  ALLOCATE(IJp(nbin,nbin), Ln_k(nbin), Kn_k(nbin), Spp(nbin))
  ALLOCATE(mv(nbin), p(nbin), v(nbin), mp(nbin), pn(nbin,nbin), s1(nbin,nbin),s2(nbin,nbin))
  ALLOCATE(SS(nbin,nbin,nbin), C(nbin,nbin,nbin))

  !print *,"Initializing starting conditions"
  call cgs%init
  !print *,cgs%grav,cgs%au,cgs%k_b,cgs%m_u,cgs%m_sun
  tt = 196
  mu = 2.37
  !cs = sqrt(tt*cgs%k_b/(mu*1.008*cgs%m_u))
  cs = 82592.515
  r = 1.0*cgs%au
  msolar = 0.5*cgs%m_sun
  Omega_K = 1.4078381e-7
  Hp = cs / Omega_K
  sigma_d = 0.18                                      ! g / cm**2
  rho_d = sigma_d/(sqrt(2*math%pi)*Hp)
  rho_s = 1.6                                         ! g / cm**3
  micron = 1e-4

  mps = 1e2
  u_f = 1*mps
  trans_width = 0.001*mps
  alpha = 1e-3
  Cdt = 0.1
  xi = 1.83

  dt = 1e10 * 60*60*24*365
  t_end = 1e13 * 60*60*24*365

  !dt=1e10
  !t_end=dt

  save_time = 1*cgs%yr
  eps_time = 1.10

                                                      ! Define binning
  min_a = 1e-4
  max_a = 1e2

  ratio_r = (max_a/min_a)**(1.0/(nbin-1))
  a(1) = min_a
  n = 0
  C = 0

  do i=2,nbin
    a(i) = a(i-1)*ratio_r                             ! Setup of bins
  end do

  st = rho_s * a * math%pi/(2*sigma_d*100)
  utm_ini = cs * sqrt(2*alpha*st)
  where (st > 1.0) utm_ini = cs * sqrt(2*alpha/ st)
  m = 4.0/3.0 * a**3 * rho_s * math%pi
  n(1) = rho_d / m(1)
  v = 1
  p = n * v
  mv = m * v

  SS = 0

  print *,"Starting n, number density:"
  print *, n
  print "(/)"

  print *,"Initial velocity:"
  print *, p
  print "(/)"

  ForAll(i = 1:nbin, j = 1:nbin) identity(i,j) = (i/j)*(j/i)
  identity = identity / dt

  do i=1,nbin
    do j=1,nbin                                       ! Setup of matrices
      minij = min(i,j)
      maxij = max(i,j)

      ms(i,j) = m(i) + m(j)

      do kn = nbin,1,-1
        if ( m(kn)<ms(i,j) ) exit
      end do

      ubm(i,j) = sqrt(8 * cgs%k_b* tt * (m(i) + m(j) )/(math%pi*m(i)*m(j) ))
      utm(i,j) = max(utm_ini(i), utm_ini(j))
      du(i,j) = ubm(i,j) + utm(i,j)
      sigma_col(i,j) = math%pi * (a(i) + a(j))**2

      if (m(i)/m(j) > 10 .or. m(i)/m(j) < 0.1 ) then        !!!
        !Cratering
        SS(1:minij,i,j) = a(1:minij) ** (-3*xi)
        Stot = SUM(m * SS(:,i,j) ) / (m(minij) *2)
        SS(1:minij,i,j) = SS(1:minij,i,j) / Stot
        m_crat =  abs(m(i)-m(j))
        m_big = max(m(i),m(j))
        m_small = min(m(i),m(j))

        do l = nbin,1,-1
          if (m(l)-m_big<-m_small) exit
        end do

        eps = (m_crat - m(l))/(m(l+1) - m(l))
        SS(l,i,j) = 1-eps
        SS(l+1,i,j) = eps

      else
        !Fragmentation
        SS(1:maxij,i,j) = a(1:maxij) ** (-3*xi)
        Stot = SUM(m * SS(:,i,j) ) / (m(i) + m(j))
        SS(1:maxij,i,j) = SS(1:maxij,i,j) / Stot
      end if


      if (kn+1 <= nbin) then
        mm(i,j) = m(kn)
        mn(i,j) = m(kn+1)
      else
        sigma_col(i,j) = 0
      end if

      if (du(i,j) < u_f - trans_width) then
        pf(i,j) = 0
      else if (du(i,j) > u_f) then
        pf(i,j) = 1
      else
        pf(i,j) = 1-(u_f-du(i,j))/trans_width
      end if

    end do
  end do

  !pf = 0
  pc = 1-pf

  where (sigma_col == 0) ep = 0
  where (sigma_col /= 0) ep = (mn - ms) / (mn - mm)

  do k=1,nbin
    where (mm == m(k) ) C(k,:,:) = ep
    where (mn == m(k) ) C(k,:,:) = 1 - ep
  end do


  kk = du * sigma_col * pc
  ll = du * sigma_col * pf                              ! Checked

  t = 0
  itime = 0

  !=============================================================================
  ! A single step, without updating, to diagnose the coagulation factors
  !=============================================================================
  !dt = t_end
  ForAll(i = 1:nbin, j = 1:nbin) identity(i,j) = (i/j)*(j/i)
  !-----------------------------------------------------------------------------
  ! Matrix notation, eliminating intermediate variables
  !-----------------------------------------------------------------------------
  do i=1,nbin
    nn(i,:) = n(i)*n
  end do
  do k=1,nbin
    S(k)    = 0.5*SUM(kk*C(k,:,:)*nn) -  SUM(kk(k,:)*nn(k,:))
    JJ(k,:) = matmul (kk*C(k,:,:), n) - n(k)*kk(k,:)
    JJ(k,k) = JJ(k,k) - SUM(kk(k,:)*n)
  end do
  !call print_2 (JJ, 'JJ')
  IJJ  = identity/dt - JJ
  !-----------------------------------------------------------------------------
  ! Solve
  !-----------------------------------------------------------------------------
  js_n = nbin
  js_nrhs = 1
  js_lda = nbin
  js_ldb = nbin
  m_tot = SUM(S*m)
  CALL DGESV( js_n, js_nrhs, IJJ, js_lda, js_ipiv, S, js_ldb, js_info )
  dn = S
  !-----------------------------------------------------------------------------
  ! Index notation:
  !           n(k):  number density in bin k
  !       C(k,i,j):  rate of transfer to bin k from i and j together
  !        kk(i,j):  velocity difference times cross-section
  !-----------------------------------------------------------------------------
  do k=1,nbin
    S(k) = 0.0
    JJ(k,:) = 0.0
    do i=1,nbin
      !-------------------------------------------------------------------------
      ! Terms corresponding to addition of particles to bin k from bins i and j
      !-------------------------------------------------------------------------
      do j=1,nbin
        S(k)    = S(k) + 0.5*kk(i,j)*C(k,i,j)*n(j)*n(i)
        JJ(k,i) = JJ(k,i) +  kk(i,j)*C(k,i,j)*n(j)
      end do
      !-------------------------------------------------------------------------
      ! Terms corresponding to removal of particles from bin k to bin i and
      ! from bin i to bin k
      !-------------------------------------------------------------------------
      S(k) = S(k)       - kk(k,i)*n(k)*n(i)
      JJ(k,i) = JJ(k,i) - kk(k,i)*n(k)
      JJ(k,k) = JJ(k,k) - kk(k,i)*n(i)
    end do
  end do
  !call print_2 (JJ, 'JJ')
  IJJ  = identity/dt - JJ
  !-----------------------------------------------------------------------------
  ! Solve
  !-----------------------------------------------------------------------------
  js_n = nbin
  js_nrhs = 1
  js_lda = nbin
  js_ldb = nbin
  m_tot = SUM(S*m)
  CALL DGESV( js_n, js_nrhs, IJJ, js_lda, js_ipiv, S, js_ldb, js_info )
  dn = S
  !=============================================================================
  !=============================================================================

  do while (t < t_end)

    mp = m * p
    do i=1,nbin
      do j=1,nbin
        nn(i,j) = n(i)*n(j)
        pn(i,j) = (mp(i) * n(j) + n(i) * mp(j))/ms(i,j)
      end do
    end do

    dn = 0
    dp = 0
    S = 0
    Sp = 0
    Jp = 0
    Jm = 0

    JJ = 0
    IJJ = 0
    IJp = 0
    Q = kk * nn
    Q_frag = ll * nn


    if (.not. do_implicit) then

      do k=1,nbin
      KC_k = kk * C(k,:,:)
      LS_k = ll * SS(k,:,:)

      Kn_k = kk(k,:) * n
      Ln_k = ll(k,:) * n

      S(k) = ( 0.5* SUM(Q*C(k,:,:)) - SUM(Q(k,:)) ) + ( 0.5 * SUM( Q_frag*SS(k,:,:)) - SUM(Q_frag(k,:)) )
      Sp(k) = 0.5 * SUM(KC_k * pn) - SUM(Kn_k * p(k)) + 0.5 * SUM(LS_k * pn) - SUM(Ln_k * p(k))
      end do

      dn = S * dt
      dp = Sp * dt
    else
      do k=1,nbin
        KC_k = kk * C(k,:,:)
        LS_k = ll * SS(k,:,:)

        Kn_k = kk(k,:) * n
        Ln_k = ll(k,:) * n

        S(k) = ( 0.5* SUM(Q*C(k,:,:)) - SUM(Q(k,:)) ) + ( 0.5 * SUM( Q_frag*SS(k,:,:)) - SUM(Q_frag(k,:)) )
        Sp(k) = 0.5 * SUM(KC_k * pn) - SUM(Kn_k * p(k)) + 0.5 * SUM(LS_k * pn) - SUM(Ln_k * p(k))
      ! end do

      ! do k=1,nbin
        s1 = KC_k * 1/ms
        s2 = LS_k * 1/ms
        Jp(k,:) = (matmul(s1, n) + matmul(s2, n)) * m
        Jp(k,k) = Jp(k,k) - (SUM(Kn_k) + SUM(Ln_k) )
        Jm(k,:) = (matmul(s1, p*m) + matmul(s2, p*m)) - p(k) * kk(k,:) - p(k) * ll(k,:)
        JJ(k,:) = (matmul(KC_k , n) - n(k) * kk(k,:) ) + (matmul(LS_k,  n) - n(k) * ll(k,:) )
        !print *, "matmul(LS_k,  n)"
        !print *, matmul(LS_k,  n)
        JJ(k,k) = JJ(k,k) - (SUM(Kn_k) + SUM(Ln_k) )                                ! Checked
      end do
      !call print_2 (JJ, 'JJ')

      ! print *, S
      ! print "(/)"
      !JJJ(:nbin,:nbin) = identity - JJ
      !JJJ(nbin:, nbin:) = identity - Jp
      !JJJ(:nbin, nbin:) = Jm
      !SSS(:,:nbin) = S
      !SSS(:,nbin:) = Sp
      IJJ = identity - JJ
      IJp = identity - Jp

      !call print_2 (IJJ, "IJJ")

      js_n = nbin
      js_nrhs = 1
      js_lda = nbin
      js_ldb = nbin

      CALL DGESV( js_n, js_nrhs, IJJ, js_lda, js_ipiv, S, js_ldb, js_info )
      Spp = Sp - matmul(Jm , S)
      CALL DGESV( js_n, js_nrhs, IJp, js_lda, js_ipiv, Spp, js_ldb, js_info )

      dn = S
      dp = Spp

      ! print *,"dp"
      ! do i = 1, nbin
      !   print '(5ES13.3)', dp(i)
      ! end do
      ! print "(/)"

    end if

    n = n + dn
    p = p + dp
    t = t + dt

    if (t > save_time) then
      save_time = save_time * eps_time
      itime = itime + 1
      if (mod(itime,1)==0) then
        print *,"Time left: ", (t_end-t)/(60*60*24*365)
        print *,"Stepping: ", dt/(60*60*24*365)
        print *, "n", n
        print *, "Implicit", do_implicit
      end if
    end if

  end do

  print *, "Lost mass to total ratio : "
  a = dn*m
  a = n*m
  print *, (SUM(a)-rho_d)/rho_d

  print *,"Final n, number density:"
  print *, n
  print "(/)"

  print *,"Final velocity:"
  print *, p/(n+1e-30)
  print "(/)"


  ! print *, "cs"
  ! print '(10ES13.3)', cs
  ! print "(/)"


  ! print *, "st"
  ! open(1, file = 'st2.txt')
  ! write(*,*) "hello"
  ! do i = 1, nbin
  !   write(*,'(5ES13.3)') st(i)
  ! end do
  ! close(1)
  ! print "(/)"

  ! print *, "utm_ini"
  ! do i = 1, nbin
  !   print '(5ES13.3)', utm_ini(i)
  ! end do
  ! print "(/)"

  ! print *, "utm"
  ! do i = 1, nbin
  !   print '(5ES13.3)', utm(i,:)
  ! end do
  ! print "(/)"

  ! print *, "ubm"
  ! do i = 1, nbin
  !   print '(5ES13.3)', ubm(i,:)
  ! end do
  ! print "(/)"

  ! print *, "du"
  ! do i = 1, nbin
  !   print '(5ES13.3)', du(i,:)
  ! end do
  ! print "(/)"

  ! print *, "L"
  ! do i = 1, nbin
  !   print '(5ES13.3)', ll(i,:)
  ! end do
  ! print "(/)"

  ! print *, "LS_k"
  ! do i = 1, nbin
  !   print '(5ES13.3)', LS_k(i,:)
  ! end do
  ! print "(/)"

  ! print *, "K"
  ! do i = 1, nbin
  !   print '(5ES13.3)', kk(i,:)
  ! end do
  ! print "(/)"

  !  print *, "KC_k"
  ! do i = 1, nbin
  !   print '(5ES13.3)', KC_k(i,:)
  ! end do
  ! print "(/)"

  ! print *, "SS"
  ! do k = 1, nbin
  !   do i = 1, nbin
  !     print '(5ES13.3)', SS(k,i,:)
  !   end do
  !   print "(/)"
  ! end do
  ! print "(/)"

  ! print *, "C"
  ! do k = 1, nbin
  !   do i = 1, nbin
  !     print '(5ES13.3)', C(k,i,:)
  !   end do
  !   print "(/)"
  ! end do
  ! print "(/)"

  ! ! print *, "Sp"
  ! ! do i = 1, nbin
  ! !   print '(5ES13.3)', Sp(i)
  ! ! end do
  ! ! print "(/)"

  ! print *, "LS_k"
  ! do i = 1, nbin
  !   print '(5ES13.3)', LS_k(i,:)
  ! end do
  ! print "(/)"

  ! print *, "S"
  ! do i = 1, nbin
  !   print '(5ES13.3)', S(i)
  ! end do
  ! print "(/)"

  ! ! print *, "Sp"
  ! ! do i = 1, nbin
  ! !   print '(5ES13.3)', Sp(i)
  ! ! end do
  ! ! print "(/)"

  ! print *, "JJ"
  ! do i = 1, nbin
  !   print '(5ES13.3)', JJ(i,:)
  ! end do
  ! print "(/)"


  ! ! print *, "Jp"
  ! ! do i = 1, nbin
  ! !   print '(5ES13.3)', Jp(i,:)
  ! ! end do
  ! ! print "(/)"

  ! ! print *, "Jm"
  ! ! do i = 1, nbin
  ! !   print '(5ES13.3)', Jm(i,:)
  ! ! end do
  ! ! print "(/)"

  ! print *, "Q"
  ! do i = 1, nbin
  !   print '(5ES13.3)', Q(i,:)
  ! end do
  ! print "(/)"

  ! print *, "Q_frag"
  ! do i = 1, nbin
  !   print '(5ES13.3)', Q_frag(i,:)
  ! end do
  ! print "(/)"

  call trace%end ()

contains

subroutine print_2 (a, label)
  real(8):: a(:,:)
  character(*):: label
  print *, label
  do i = 1, nbin
    print '(5ES13.3)', a(i,:)
  end do
  print "(/)"
end subroutine
END SUBROUTINE test

END MODULE dust_evolve_mod
