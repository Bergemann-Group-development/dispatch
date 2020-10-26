!===============================================================================
!> Port and test dust evolution from Podolak_1.py
!===============================================================================
PROGRAM test
  USE omp_lib
  USE units_mod
  USE math_mod
  implicit none
  real:: T, mu, cs, r, msolar, Omega_K, sigma_d, rho_d, rho_s, a_s, m_s, u_f
  real:: trans_width, alpha, Cdt, xi, dt, t_end, save_time, eps_time, Hp
  real:: micron, mps


  integer:: nbin
  integer:: i, j, k ! input?
  real:: a(200), n(200), m(200), ratio_r, min_a, max_a, st(200), ubm(200), utm(200)
  real:: ms(200,200)
  ! real:: a(nbin), n(nbin), st(nbin), ms(nbin,nbin), m(nbin), utm(nbin), ubm(nbin)
  ! real:: sigma_col(nbin,nbin), mn(nbin,nbin), mm(nbin,nbin), pf(nbin,nbin)
  ! real:: SS(nbin,nbin,nbin), ep(nbin,nbin), C(nbin,nbin,nbin)

  ! real:: dn(nbin), S(nbin), J(nbin,nbin), nn(nbin,nbin), Q(nbin,nbin), Q_frag(nbin,nbin)
  !-----------------------------------------------------------------------------
  print *,"Initializing starting conditions"
  call cgs%init
  print *,cgs%grav,cgs%au,cgs%k_b,cgs%m_u,cgs%m_sun
  T = 196
  mu = 2.37
  cs = sqrt(T*cgs%k_b/(mu*cgs%m_u))
  print *,'cs', cs
  r = 1.0*cgs%au
  msolar = 0.5*cgs%m_sun
  Omega_K = sqrt(cgs%grav*msolar)/r**1.5
  print *,'Omega_K', Omega_K
  Hp = cs / Omega_K
  sigma_d = 0.18
  rho_d = sigma_d/(2*math%pi)*Hp
  rho_s = 1.6
  micron = 1e-4

  mps = 1e2
  u_f = 1*mps
  trans_width = 0.001*mps
  alpha = 1e-3
  Cdt = 0.1
  xi = 1.83

  dt = 0.1 * 60*60*24*365
  t_end = 1e3 * 60*60*24*365
  save_time = 1*cgs%yr
  eps_time = 1.10

  nbin = 200
  min_a = 1e-4
  max_a = 1e2
  
  ratio_r = (max_a/min_a)**(1.0/(nbin-1))
  a(1) = min_a
  !a(2:nbin) = (/ (a(i-1)*ratio_r, i=2,nbin) /)

  n = (/ (0, i=1,nbin) /)
  n(1) = rho_d / m(1)

  do i=2,nbin
    a(i) = a(i-1)*ratio_r
  end do

  st = rho_s * a * math%pi/(2*sigma_d*100)
  utm = cs * sqrt(2*alpha*st)
  where (st > 1.0) utm = cs * sqrt(2*alpha/ st)
  m = 4.0/3.0 * a**3 * rho_s * math%pi

  do i=1,nbin
    do j=1,nbin
      ms(i,j) = m(i) + m(j)
      
    end do
  end do
  
END PROGRAM
