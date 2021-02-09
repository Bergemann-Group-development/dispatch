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
  implicit none
  private
  type, public:: dust_evolve_t
    logical:: on=.false.
    integer:: n_bin
    real, allocatable :: m_cgs(:)
  contains
    procedure:: init
    procedure:: update
    procedure:: test
  end type
  integer, parameter :: precision_e = selected_real_kind(15, 100)
  integer:: verbose=0
  logical:: on=.false.
  logical:: do_brownian=.true.
  logical:: do_in_cell_turbulence=.false.
  real(KIND=precision_e):: tt_dust = 10             ! Dust temperature
  real(KIND=precision_e):: alpha = 1e-3             ! Turbulence parameter, used when do_in_cell_turbulence is true for calculation of turbulent velocities
  real(KIND=precision_e):: sigma_d = 0.18           ! Column density of the dust, used when do_in_cell_turbulence is true for calculation of turbulent velocities
  real(KIND=precision_e):: xi = 1.83                ! The exponent used for the distribution in the case of cratering
  real(KIND=precision_e):: trans_width = 0.1        ! Transition width for fragmentation probability in m/s
  real(KIND=precision_e):: u_f = 1e2                ! Fragmentation speed, free parameter, assumed value in m/s
  real(KIND=precision_e):: dt = 1e-1 * 60*60*24*365 ! Time step in dust_evolve in years
  type(dust_evolve_t), public:: dust_evolve
CONTAINS

!===============================================================================
!> Read parameters and initialize -- add more parameters, if they are known
!> in the particle solver
!===============================================================================
SUBROUTINE init (self, n_bin, m_cgs, rho_s)
  class(dust_evolve_t) :: self
  integer              :: n_bin               ! number of bins
  real                 :: m_cgs(:)            ! particle masses
  real                 :: rho_s               ! particle solid density
  !.............................................................................
  integer              :: iostat              ! missing(<0), syntax errors (>0)
  logical, save        :: first_time=.true.   ! shared memory

  real:: cs, r, msolar, Omega_K, t
  real:: save_time, eps_time, Hp
  real:: micron, mps

  integer:: i, j, k, l, nbin, kn, minij, maxij, itime, js_n, js_nrhs, js_lda, &
    js_ldb, js_info
  integer, ALLOCATABLE, DIMENSION(:):: js_ipiv
  real(precision_e):: ratio_r, min_a, max_a, m_crat, eps, Stot
  real(precision_e):: m_big, m_small
  real(precision_e), ALLOCATABLE, DIMENSION(:):: m, n, S, dn, dpx, dpy, dpz, dd, a, Spx, Spy, Spz, &
    mv, px, py, pz, v, mpx, mpy, mpz, Ln_k, Kn_k, Sppx, Sppy, Sppz
  real(precision_e), ALLOCATABLE, DIMENSION(:,:):: mn, mm, ms, nn, Q, Q_frag, jj, kk, &
    ll, KC_k, LS_k, IJJ, sigma_col, ep, identity, pf, pc, du, ubm, utm, &
    Jp, Jmx, Jmy, Jmz, IJp, pnx, pny, pnz, s1, s2
  real, ALLOCATABLE, DIMENSION(:):: st, utm_ini
  real(precision_e), ALLOCATABLE, DIMENSION(:,:,:):: SS, C                                                      
  logical :: do_implicit=.true.
  !--- add parameters here -----------------------------------------------------
  namelist /dust_evolve_params/ on, verbose, do_brownian, do_in_cell_turbulence, tt_dust, sigma_d, alpha, xi, trans_width, u_f, dt
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%init')
  !-----------------------------------------------------------------------------
  ! Make sure that only one thread reads the namelist
  !-----------------------------------------------------------------------------
  call io%lock%set ('dust_evolve_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, dust_evolve_params, iostat=iostat)
    if (iostat < 0) then                       ! if syntax errors
      rewind (io%input)                        ! rewind to top of file
      read (io%input, dust_evolve_params)      ! read againg to get more info
    end if
    call io%namelist_warning ('dust_evolve_params', iostat)
    write (io%namelist, dust_evolve_params)    ! write to data/run/params.nml
    self%n_bin = n_bin
    allocate(self%m_cgs(n_bin))
    self%m_cgs = m_cgs 
    call self%test                             ! run a self%test once
    
    first_time = .false.
  end if
  call io%lock%set ('dust_evolve_t%init')
  self%on = on                                 ! task private copy
  self%n_bin = n_bin
  if (.not. allocated(self%m_cgs)) allocate(self%m_cgs(n_bin))
  self%m_cgs = m_cgs


  !class(dust_evolve_t):: self

  nbin = self%n_bin
  ALLOCATE(a(nbin), n(nbin), m(nbin), st(nbin), dn(nbin), S(nbin), ms(nbin,nbin))
  ALLOCATE(utm_ini(nbin), sigma_col(nbin,nbin), mn(nbin,nbin), mm(nbin,nbin))
  ALLOCATE(pf(nbin,nbin), pc(nbin,nbin), du(nbin,nbin), ep(nbin,nbin))
  ALLOCATE(ubm(nbin,nbin), utm(nbin,nbin), Q(nbin,nbin), Q_frag(nbin,nbin))
  ALLOCATE(jj(nbin,nbin), SS(nbin,nbin,nbin), C(nbin,nbin,nbin), kk(nbin,nbin))
  ALLOCATE(ll(nbin,nbin), nn(nbin,nbin), identity(nbin,nbin), js_ipiv(nbin))
  ALLOCATE(KC_k(nbin,nbin), LS_k(nbin,nbin) , IJJ(nbin,nbin))
  ALLOCATE(Jp(nbin,nbin), Jmx(nbin,nbin), Jmy(nbin,nbin), Jmz(nbin,nbin), Spx(nbin), Spy(nbin), Spz(nbin), dpx(nbin), dpy(nbin), dpz(nbin))
  ALLOCATE(IJp(nbin,nbin), Ln_k(nbin), Kn_k(nbin), Sppx(nbin), Sppy(nbin), Sppz(nbin))
  ALLOCATE(mv(nbin), px(nbin), py(nbin), pz(nbin), v(nbin), mpx(nbin), mpy(nbin), mpz(nbin), pnx(nbin,nbin), pny(nbin,nbin), pnz(nbin,nbin), s1(nbin,nbin),s2(nbin,nbin))

  !print *,"Initializing starting conditions"
  call cgs%init
  !print *,cgs%grav,cgs%au,cgs%k_b,cgs%m_u,cgs%m_sun
  !mu = 2.37
  ! cs = sqrt(tt_dust*cgs%k_b/(mu*1.008*cgs%m_u))
  ! r = 1.0*cgs%au
  ! msolar = 0.5*cgs%m_sun
  ! Omega_K = sqrt(cgs%grav*msolar)/r**1.5
  ! Hp = cs / Omega_K
  ! sigma_d = 0.18                                      ! g / cm**2
  ! rho_d = sigma_d/(sqrt(2*math%pi)*Hp)
  !micron = 1e-4

  !mps = 1e2 
  !u_f = 1*mps
  !trans_width = 0.001*mps
  !alpha = 1e-3
  !Cdt = 0.1
  !xi = 1.83

  !dt = 1e-1 * 60*60*24*365
  !t_end = 1e2 * 60*60*24*365
  save_time = 1*cgs%yr
  eps_time = 1.10

                                                      ! Define binning
  !min_a = 1e-4
  !max_a = 1e2

  !ratio_r = (max_a/min_a)**(1.0/(nbin-1))
  !a(1) = min_a
  !n = 0

  ! do i=2,nbin
  !   a(i) = a(i-1)*ratio_r                             ! Setup of bins
  ! end do

  !m = 4.0/3.0 * a**3 * rho_s * math%pi
  !n(1) = rho_d / m(1)
  !v = 1
  !p = n * v
  !mv = m * v 

  a = (3.0/4.0 * m_cgs/(rho_s * math%pi) )**(1/3)       ! Initial setup of bin distribution
  st = rho_s * a * math%pi/(2*sigma_d*100)              ! 


  C = 0
  SS = 0

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

      ! if (du(i,j) < u_f - trans_width) then
      !   pf(i,j) = 0
      ! else if (du(i,j) > u_f) then
      !   pf(i,j) = 1
      ! else
      !   pf(i,j) = 1-(u_f-du(i,j))/trans_width
      ! end if


    end do
  end do

  pc = 1-pf
  where (sigma_col == 0) ep = 0
  where (sigma_col /= 0) ep = (mn - ms) / (mn - mm)

  do k=1,nbin
    where (mm == m(k) ) C(k,:,:) = ep
    where (mn == m(k) ) C(k,:,:) = 1 - ep
  end do

  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the dust state for one cell, with _in and _out arguments
!> holding particle mass density dd(n_bin) and velocity vd(n_bin,3)
!===============================================================================
SUBROUTINE update (self, dt_dispatch, vrms_in, dd_in, vd_in, dd_out, vd_out)
  class(dust_evolve_t):: self
  real(kind=8), intent(in)  :: dt_dispatch
  real,         intent(in)  :: vrms_in(:), dd_in(:), vd_in(:,:)
  real,         intent(out) :: dd_out(:), vd_out(:,:)

  real:: tt, mu, cs, r, msolar, Omega_K, sigma_d, rho_d, rho_s, u_f, t
  real:: trans_width, alpha, Cdt, xi, dt, t_end, save_time, eps_time, Hp
  real:: micron, mps

  integer:: i, j, k, l, nbin, kn, minij, maxij, itime, js_n, js_nrhs, js_lda, &
    js_ldb, js_info
  integer, ALLOCATABLE, DIMENSION(:):: js_ipiv
  real(precision_e):: ratio_r, min_a, max_a, m_crat, eps, Stot
  real(precision_e):: m_big, m_small
  real(precision_e), ALLOCATABLE, DIMENSION(:):: m, n, S, dn, dpx, dpy, dpz, dd, a, Spx, Spy, Spz, &
    mv, px, py, pz, v, mpx, mpy, mpz, Ln_k, Kn_k, Sppx, Sppy, Sppz
  real(precision_e), ALLOCATABLE, DIMENSION(:,:):: mn, mm, ms, nn, Q, Q_frag, jj, kk, &
    ll, KC_k, LS_k, IJJ, sigma_col, ep, identity, pf, pc, du, ubm, utm, &
    Jp, Jmx, Jmy, Jmz, IJp, pnx, pny, pnz, s1, s2
  real, ALLOCATABLE, DIMENSION(:):: st, utm_ini
  real(precision_e), ALLOCATABLE, DIMENSION(:,:,:):: SS, C                                                      ! Necesarry?
  logical :: do_implicit=.true.

  ! CGS units below
  real(precision_e) :: tend, vd(self%n_bin,3), vrms(self%n_bin)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%update', itimer=itimer)
  ! convert from code units to cgs
  tend = scaling%t * dt_dispatch
  n = scaling%d * dd_in / self%m_cgs
  vd = scaling%u * vd_in
  vrms = scaling%u * vrms_in

  px = n * vd(:,1)
  py = n * vd(:,2)
  pz = n * vd(:,3)

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
    !cs = sqrt(tt_dust*cgs%k_b/(mu*1.008*cgs%m_u))
    cs = sqrt(tt_dust*cgs%k_B/(cgs%m_u*scaling%mu))/scaling%u
    utm_ini = cs * sqrt(2*alpha*st)                     ! Handle velocities every update
    where (st > 1.0) utm_ini = cs * sqrt(2*alpha/ st)
  end if

  do i=1,nbin
    do j=1,nbin
      du(i,j) = sqrt((vd(i,1) - vd(j,1))**2 + (vd(i,2) - vd(j,2))**2 + (vd(i,3) - vd(j,3))**2)
      du(i,j) = du(i,j) + sqrt(vrms(i)**2 + vrms(j)**2)
      
      if (do_brownian) then
      ubm(i,j) = sqrt(8 * cgs%k_b* tt * (m(i) + m(j) )/(math%pi*m(i)*m(j) ))
      du(i,j) = du(i,j) + ubm(i,j)
      end if

      if (do_in_cell_turbulence) then
      utm(i,j) = max(utm_ini(i), utm_ini(j))
      du(i,j) = du(i,j) + utm(i,j)
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
  

  kk = du * sigma_col * pc
  ll = du * sigma_col * pf                              

  t = 0
  itime = 0
  do while (t < tend)
    
    mpx = m * px
    mpy = m * py
    mpz = m * pz                                                  ! xyz
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
      Spx(k) = 0.5 * SUM(KC_k * pnx) - SUM(Kn_k * px(k)) + 0.5 * SUM(LS_k * pnx) - SUM(Ln_k * px(k))
      Spy(k) = 0.5 * SUM(KC_k * pny) - SUM(Kn_k * py(k)) + 0.5 * SUM(LS_k * pny) - SUM(Ln_k * py(k))
      Spz(k) = 0.5 * SUM(KC_k * pnz) - SUM(Kn_k * pz(k)) + 0.5 * SUM(LS_k * pnz) - SUM(Ln_k * pz(k))
      end do  
      
      dn = S * dt
      dpx = Spx * dt
      dpy = Spy * dt
      dpz = Spz * dt
    else
      do k=1,nbin
        KC_k = kk * C(k,:,:)
        LS_k = ll * SS(k,:,:)

        Kn_k = kk(k,:) * n
        Ln_k = ll(k,:) * n

        S(k) = ( 0.5* SUM(Q*C(k,:,:)) - SUM(Q(k,:)) ) + ( 0.5 * SUM( Q_frag*SS(k,:,:)) - SUM(Q_frag(k,:)) )
        Spx(k) = 0.5 * SUM(KC_k * pnx) - SUM(Kn_k * px(k)) + 0.5 * SUM(LS_k * pnx) - SUM(Ln_k * px(k))
        Spy(k) = 0.5 * SUM(KC_k * pny) - SUM(Kn_k * py(k)) + 0.5 * SUM(LS_k * pny) - SUM(Ln_k * py(k))
        Spz(k) = 0.5 * SUM(KC_k * pnz) - SUM(Kn_k * pz(k)) + 0.5 * SUM(LS_k * pnz) - SUM(Ln_k * pz(k))

        s1 = KC_k * 1/ms
        s2 = LS_k * 1/ms
        Jp(k,:) = (matmul(s1, n) + matmul(s2, n)) * m
        Jp(k,k) = Jp(k,k) - (SUM(Kn_k) + SUM(Ln_k) )                            
        Jmx(k,:) = (matmul(s1, px*m) + matmul(s2, px*m)) - px(k) * kk(k,:) - px(k) * ll(k,:)
        Jmy(k,:) = (matmul(s1, py*m) + matmul(s2, py*m)) - py(k) * kk(k,:) - py(k) * ll(k,:)
        Jmz(k,:) = (matmul(s1, pz*m) + matmul(s2, pz*m)) - pz(k) * kk(k,:) - pz(k) * ll(k,:)
        JJ(k,:) = (matmul(KC_k , n) - n(k) * kk(k,:) ) + (matmul(LS_k,  n) - n(k) * ll(k,:) )
        !print *, "matmul(LS_k,  n)"
        !print *, matmul(LS_k,  n)
        JJ(k,k) = JJ(k,k) - (SUM(Kn_k) + SUM(Ln_k) )                                ! Checked
      end do

      IJJ = identity - JJ
      IJp = identity - Jp
      
      
      js_n = nbin
      js_nrhs = 1
      js_lda = nbin
      js_ldb = nbin
      
      CALL DGESV( js_n, js_nrhs, IJJ, js_lda, js_ipiv, S, js_ldb, js_info )
      Sppx = Spx - matmul(Jmx , S)
      Sppy = Spy - matmul(Jmy , S)
      Sppz = Spz - matmul(Jmz , S)
      CALL DGESV( js_n, js_nrhs, IJp, js_lda, js_ipiv, Sppx, js_ldb, js_info )
      CALL DGESV( js_n, js_nrhs, IJp, js_lda, js_ipiv, Sppy, js_ldb, js_info )
      CALL DGESV( js_n, js_nrhs, IJp, js_lda, js_ipiv, Sppz, js_ldb, js_info )

      dn = S
      dpx = Sppx
      dpy = Sppy
      dpz = Sppz

    end if

    n = n + dn
    px = px + dpx
    py = py + dpy
    pz = pz + dpz
    t = t + dt

    if (t > save_time) then
      save_time = save_time * eps_time
      itime = itime + 1
      if (mod(itime,10)==0) then
        print *,"Time left: ", (t_end-t)/(60*60*24*365)
        ! print *, "n", n(1:10)
        ! print *, "p", p(1:10)
      end if
    end if

  end do

  print *,"Final n, number density:"
  print *, n
  print "(/)"

  print *,"Final velocity x:"
  print *, px/(n+1e-30)
  print "(/)"

  vd(:,1) = px/n
  vd(:,2) = py/n
  vd(:,3) = pz/n

  ! convert back to code units
  dd_out = self%m_cgs * n / scaling%d
  vd_out = vd / scaling%u

  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Self test
!===============================================================================
SUBROUTINE test (self)
  class(dust_evolve_t):: self
  real:: tt, mu, cs, r, msolar, Omega_K, sigma_d, rho_d, rho_s, u_f, t
  real:: trans_width, alpha, Cdt, xi, dt, t_end, save_time, eps_time, Hp
  real:: micron, mps

  integer:: i, j, k, l, nbin, kn, minij, maxij, itime, js_n, js_nrhs, js_lda, &
    js_ldb, js_info
  integer, ALLOCATABLE, DIMENSION(:):: js_ipiv
  real(precision_e):: ratio_r, min_a, max_a, m_crat, eps, Stot
  real(precision_e):: m_big, m_small
  real(precision_e), ALLOCATABLE, DIMENSION(:):: m, n, S, dn, dp, dd, a, Sp, &
    mv, p, v, mp, Ln_k, Kn_k, Spp
  real(precision_e), ALLOCATABLE, DIMENSION(:,:):: mn, mm, ms, nn, Q, Q_frag, jj, kk, &
    ll, KC_k, LS_k, IJJ, sigma_col, ep, identity, pf, pc, du, ubm, utm, &
    Jp, Jm, IJp, pn, s1, s2
  real, ALLOCATABLE, DIMENSION(:):: st, utm_ini
  real(precision_e), ALLOCATABLE, DIMENSION(:,:,:):: SS, C                                                      ! Necesarry?
  logical :: do_implicit=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%test')
  nbin = self%n_bin
  !nbin = 5

  ALLOCATE(a(nbin), n(nbin), m(nbin), st(nbin), dn(nbin), S(nbin), ms(nbin,nbin))
  ALLOCATE(utm_ini(nbin), sigma_col(nbin,nbin), mn(nbin,nbin), mm(nbin,nbin))
  ALLOCATE(pf(nbin,nbin), pc(nbin,nbin), du(nbin,nbin), ep(nbin,nbin))
  ALLOCATE(ubm(nbin,nbin), utm(nbin,nbin), Q(nbin,nbin), Q_frag(nbin,nbin))
  ALLOCATE(jj(nbin,nbin), SS(nbin,nbin,nbin), C(nbin,nbin,nbin), kk(nbin,nbin))
  ALLOCATE(ll(nbin,nbin), nn(nbin,nbin), identity(nbin,nbin), js_ipiv(nbin))
  ALLOCATE(KC_k(nbin,nbin), LS_k(nbin,nbin) , IJJ(nbin,nbin))
  ALLOCATE(Jp(nbin,nbin), Jm(nbin,nbin), Sp(nbin), dp(nbin))
  ALLOCATE(IJp(nbin,nbin), Ln_k(nbin), Kn_k(nbin), Spp(nbin))
  ALLOCATE(mv(nbin), p(nbin), v(nbin), mp(nbin), pn(nbin,nbin), s1(nbin,nbin),s2(nbin,nbin))

  !print *,"Initializing starting conditions"
  call cgs%init
  !print *,cgs%grav,cgs%au,cgs%k_b,cgs%m_u,cgs%m_sun
  tt = 196
  mu = 2.37
  cs = sqrt(tt*cgs%k_b/(mu*1.008*cgs%m_u))
  r = 1.0*cgs%au
  msolar = 0.5*cgs%m_sun
  Omega_K = sqrt(cgs%grav*msolar)/r**1.5
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

  dt = 1e-1 * 60*60*24*365
  t_end = 1e2 * 60*60*24*365
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

  pc = 1-pf
  where (sigma_col == 0) ep = 0
  where (sigma_col /= 0) ep = (mn - ms) / (mn - mm)

  do k=1,nbin
    where (mm == m(k) ) C(k,:,:) = ep
    where (mn == m(k) ) C(k,:,:) = 1 - ep
  end do


  kk = du * sigma_col * pc
  ll = du * sigma_col * pf                              ! Checked

  ! print *, "cs"
  ! print '(10ES13.3)', cs
  ! print "(/)"

  ! print *, "st"
  ! do i = 1, nbin
  !   print '(5ES13.3)', st(i)
  ! end do
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

  !  print *, "ubm"
  ! do i = 1, nbin
  !   print '(5ES13.3)', ubm(i,:)
  ! end do
  ! print "(/)"

  ! print *, "L"
  ! do i = 1, nbin
  !   print '(5ES13.3)', ll(i,:)
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

  t = 0
  itime = 0

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
    
      ! print *, "Sp"
      ! do i = 1, nbin
      !   print '(5ES13.3)', Sp(i)
      ! end do
      ! print "(/)"
      
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

      ! print *, "Sp"
      ! do i = 1, nbin
      !   print '(5ES13.3)', Sp(i)
      ! end do
      ! print "(/)"

      ! print *, "JJ"
      ! do i = 1, nbin
      !   print '(5ES13.3)', JJ(i,:)
      ! end do
      ! print "(/)"

      ! print *, "Jp"
      ! do i = 1, nbin
      !   print '(5ES13.3)', Jp(i,:)
      ! end do
      ! print "(/)"

      ! print *, "Jm"
      ! do i = 1, nbin
      !   print '(5ES13.3)', Jm(i,:)
      ! end do
      ! print "(/)"


      ! print *, S
      ! print "(/)"
      !JJJ(:nbin,:nbin) = identity - JJ
      !JJJ(nbin:, nbin:) = identity - Jp
      !JJJ(:nbin, nbin:) = Jm
      !SSS(:,:nbin) = S
      !SSS(:,nbin:) = Sp
      IJJ = identity - JJ
      IJp = identity - Jp
      
      ! print *,"JJ"
      ! do i = 1, nbin
      !   print '(5ES13.3)', JJ(i,:)
      ! end do
      ! print "(/)"
      
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
      if (mod(itime,10)==0) then
        print *,"Time left: ", (t_end-t)/(60*60*24*365)
        print *, "n", n(1:10)
        print *, "p", p(1:10)
      end if
    end if

  end do

  print *,"Final n, number density:"
  print *, n
  print "(/)"

  print *,"Final velocity:"
  print *, p/(n+1e-30)
  print "(/)"

  call trace%end ()
END SUBROUTINE test

END MODULE dust_evolve_mod
