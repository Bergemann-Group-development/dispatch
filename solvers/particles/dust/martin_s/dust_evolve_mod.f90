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
  contains
    procedure:: init
    procedure:: update
    procedure:: test
  end type
  integer:: verbose=0
  logical:: on=.false.
  type(dust_evolve_t), public:: dust_evolve
CONTAINS

!===============================================================================
!> Read parameters and initialize -- add more parameters, if they are known
!> in the particle solver
!===============================================================================
SUBROUTINE init (self, n_bin, m_cgs)
  class(dust_evolve_t) :: self
  integer              :: n_bin               ! number of bins
  real                 :: m_cgs(:)            ! particle masses
  !.............................................................................
  integer              :: iostat              ! missing(<0), syntax errors (>0)
  logical, save        :: first_time=.true.   ! shared memory
  !--- add parameters here -----------------------------------------------------
  namelist /dust_evolve_params/ on, verbose
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
    call self%test                             ! run a self%test once
    first_time = .false.
  end if
  call io%lock%set ('dust_evolve_t%init')
  self%on = on                                 ! task private copy
  self%n_bin = n_bin
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update the dust state for one cell, with _in and _out arguments
!> holding particle mass density dd(n_bin) and velocity vd(n_bin,3)
!===============================================================================
SUBROUTINE update (self, vrms, dd_in, vd_in, dd_out, vd_out)
  class(dust_evolve_t):: self
  real:: vrms(:), dd_in(:), dd_out(:), vd_in(:,:), vd_out(:,:)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%update', itimer=itimer)
  !....= ... dd_in
  !....= vd_in(i) ...
  !... dd_out = ...
  !... vd_out(i) = ...
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Self test
!===============================================================================
SUBROUTINE test (self)
  class(dust_evolve_t):: self
  real:: tt, mu, cs, r, msolar, Omega_K, sigma_d, rho_d, rho_s, u_f, t
  real:: trans_width, alpha, Cdt, xi, dt, t_end, save_time, eps_time, Hp
  real:: micron, mps

  integer, parameter :: dp = selected_real_kind(15, 100)
  integer:: i, j, k, l, nbin, kn, minij, maxij, itime, js_n, js_nrhs, js_lda, &
    js_ldb, js_info
  integer, ALLOCATABLE, DIMENSION(:):: js_ipiv
  real(dp):: ratio_r, min_a, max_a, m_crat, eps, Stot
  real(dp):: m_big, m_small
  real(dp), ALLOCATABLE, DIMENSION(:):: m, n, S, dn, a
  real(dp), ALLOCATABLE, DIMENSION(:,:):: mn, mm, ms, nn, Q, Q_frag, jj, kk, &
    ll, KC_k, LS_k, IJJ, sigma_col, ep, identity, pf, pc, du, ubm, utm
  real, ALLOCATABLE, DIMENSION(:):: st, utm_ini
  real(dp), ALLOCATABLE, DIMENSION(:,:,:):: SS, C                                                      ! Necesarry?
  logical :: do_implicit=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_evolve_t%test')
  nbin = self%n_bin

  ALLOCATE(a(nbin), n(nbin), m(nbin), st(nbin), dn(nbin), S(nbin), ms(nbin,nbin))
  ALLOCATE(utm_ini(nbin), sigma_col(nbin,nbin), mn(nbin,nbin), mm(nbin,nbin))
  ALLOCATE(pf(nbin,nbin), pc(nbin,nbin), du(nbin,nbin), ep(nbin,nbin))
  ALLOCATE(ubm(nbin,nbin), utm(nbin,nbin), Q(nbin,nbin), Q_frag(nbin,nbin))
  ALLOCATE(jj(nbin,nbin), SS(nbin,nbin,nbin), C(nbin,nbin,nbin), kk(nbin,nbin))
  ALLOCATE(ll(nbin,nbin), nn(nbin,nbin), identity(nbin,nbin), js_ipiv(nbin))
  ALLOCATE(KC_k(nbin,nbin), LS_k(nbin,nbin) , IJJ(nbin,nbin))

  !print *,"Initializing starting conditions"
  call cgs%init
  !print *,cgs%grav,cgs%au,cgs%k_b,cgs%m_u,cgs%m_sun
  tt = 196
  mu = 2.37
  cs = sqrt(tt*cgs%k_b/(mu*cgs%m_u))
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

  dt = 1 * 60*60*24*365
  t_end = 1e3 * 60*60*24*365
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

  SS = 0

  print *,"Starting n, number density:"
  print *, n
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

  ! print *, "ms"
  ! do i = 1, nbin
  !   print '(5ES20.6)', ms(i,:)
  ! end do
  ! print "(/)"

  ! print *, "mn"
  ! do i = 1, nbin
  !   print '(5ES20.6)', mn(i,:)
  ! end do
  ! print "(/)"

  ! print *, "sigma_col"
  ! do i = 1, nbin
  !   print '(5ES20.6)', sigma_col(i,:)
  ! end do
  ! print "(/)"

  do k=1,nbin
    where (mm == m(k) ) C(k,:,:) = ep
    where (mn == m(k) ) C(k,:,:) = 1 - ep
  end do


  kk = du * sigma_col * pc
  ll = du * sigma_col * pf                              ! Checked

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
    dn = 0

    do i=1,nbin
      do j=1,nbin
        nn(i,j) = n(i)*n(j)
      end do
    end do

    S = 0
    JJ = 0
    IJJ = 0
    Q = kk * nn
    Q_frag = ll * nn

    do k=1,nbin
      S(k) = ( 0.5* SUM(Q*C(k,:,:)) - SUM(Q(k,:)) ) + ( 0.5 * SUM( Q_frag*SS(k,:,:)) - SUM(Q_frag(k,:)) )
    end do

    if (.not. do_implicit) then
      dn = S * dt
    else
      do k=1,nbin

        KC_k = kk * C(k,:,:)
        LS_k = ll * SS(k,:,:)
        JJ(k,:) = (matmul(KC_k , n) - n(k) * kk(k,:) ) + (matmul(LS_k,  n) - n(k) * ll(k,:) )
        JJ(k,k) = JJ(k,k) - (SUM(kk(k,:) * n) + SUM(ll(k,:) * n) )                                ! Checked
      end do

      ! print *, S
      ! print "(/)"

      IJJ = identity - JJ

      ! print *,"JJ"
      ! do i = 1, nbin
      !   print '(5ES13.3)', JJ(i,:)
      ! end do
      ! print "(/)"

      ! print *, "I-JJ"
      ! do i = 1, nbin
      !   print '(5ES13.3)', IJJ(i,:)
      ! end do
      ! print "(/)"

      js_n = nbin
      js_nrhs = 1
      js_lda = nbin
      js_ldb = nbin

      CALL DGESV( js_n, js_nrhs, IJJ, js_lda, js_ipiv, S, js_ldb, js_info )

      ! print *, S
      ! print "(/)"
      ! print *, js_ipiv
      ! print "(/)"

      dn = S
    end if

    n = n + dn
    t = t + dt

    if (t > save_time) then
      save_time = save_time * eps_time
      itime = itime + 1
      if (mod(itime,10)==0) then
        print *,"Time left: ", (t_end-t)/(60*60*24*365), n(1)
      end if
    end if

  end do

  print *,"Final n, number density:"
  print *, n
  print "(/)"

  call trace%end ()
END SUBROUTINE test

END MODULE dust_evolve_mod
