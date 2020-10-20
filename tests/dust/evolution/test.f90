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
  !-----------------------------------------------------------------------------
  print *,"Initializing starting conditions"
  T = 196
  mu = 2.37
  cs = sqrt(T*cgs%k_b/mu*cgs%m_u)
  r = 1.0*cgs%au
  msolar = 0.5*cgs%m_sun
  Omega_K = sqrt(cgs%grav*msolar/r**3)
  Hp = cs / Omega_K
  sigma_d = 0.18
  !rho_d = sigma_d/(2*math%pi)*Hp
  rho_s = 1.6
  micron = 1e-4
  a_s = 1*micron
  m_s = 4./3.*math%pi*a_s**3 * rho_s
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

END PROGRAM
