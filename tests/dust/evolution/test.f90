!===============================================================================
!> Port and test dust evolution from Podolak_1.py
!===============================================================================
PROGRAM test
  USE omp_lib
  USE units_mod
  USE math_mod
  implicit none
  real:: T, mu, cs, r, msolar, Omega_K, sigma_d, rho_d, rho_s, u_f
  real:: trans_width, alpha, Cdt, xi, dt, t_end, save_time, eps_time, Hp
  real:: micron, mps

  integer:: i, j, k, l, nbin, kn, minij, maxij
  real:: ratio_r, min_a, max_a, m_crat, eps, Stot
  real, ALLOCATABLE, DIMENSION(:):: a, n, m, st, dn, S, utm_ini
  real, ALLOCATABLE, DIMENSION(:,:):: ms, sigma_col, mn, mm, pf, pc, ep, ubm, utm, du, Q, Q_frag, jj, kk, ll
  real, ALLOCATABLE, DIMENSION(:,:,:):: SS, C
  !-----------------------------------------------------------------------------
  nbin = 5
  ALLOCATE(a(nbin), n(nbin), m(nbin), st(nbin), dn(nbin), S(nbin), ms(nbin,nbin), utm_ini(nbin))
  ALLOCATE(sigma_col(nbin,nbin), mn(nbin,nbin), mm(nbin,nbin), pf(nbin,nbin), pc(nbin,nbin), du(nbin,nbin))
  ALLOCATE(ep(nbin,nbin), ubm(nbin,nbin), utm(nbin,nbin), Q(nbin,nbin))
  ALLOCATE(Q_frag(nbin,nbin), jj(nbin,nbin), SS(nbin,nbin,nbin), C(nbin,nbin,nbin))
  ALLOCATE(kk(nbin,nbin),ll(nbin,nbin))


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

  do i=1,nbin
    do j=1,nbin                                       ! Setup of matrices
      minij = min(i,j)
      maxij = max(i,j)

      ms(i,j) = m(i) + m(j)

      do kn = nbin,1,-1
        if ( m(kn)<ms(i,j) ) exit
      end do 

      ubm(i,j) = sqrt(8 * cgs%k_b* T * ((m(i) + m(j)))/(math%pi*( m(i)*m(j) )))
      utm(i,j) = max(utm_ini(i), utm_ini(j))
      du(i,j) = ubm(i,j) + utm(i,j)
      sigma_col(i,j) = math%pi * (a(i) + a(j))**2

      if (m(i)/m(j) > 10 .or. m(i)/m(j) < 0.1 ) then
        !Cratering 
        SS(:minij,i,j) = a(:minij) ** (-3*xi)
        Stot = SUM(m * SS(:,i,j) ) / (m(minij) *2)
        SS(:minij,i,j) = SS(:minij,i,j) / Stot
        m_crat =  abs(m(i)-m(j))

        do l = nbin,1,-1
          if (m(l)<ms(i,j) ) exit
        end do 

        eps = (m_crat - m(l-1))/(m(l) - m(l-1))
        SS(l-1,i,j) = 1-eps
        SS(l,i,j) = eps   

      else
        !Fragmentation
        SS(:maxij,i,j) = a(:maxij) ** (-3*xi)
        Stot = SUM(m * SS(:,i,j) ) / (m(i) + m(j)) 
        SS(:maxij,i,j) = SS(:maxij,i,j) / Stot
      end if


      if (kn+1 < nbin) then
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
    where (m(k) == mm) C(k,:,:) = ep
    where (m(k) == mn) C(k,:,:) = 1 - ep
  end do


  kk = du * sigma_col * pc
  ll = du * sigma_col * pf
  
END PROGRAM
