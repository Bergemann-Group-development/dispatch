!> =============================================================================
!> =============================================================================
MODULE mhd_mod
  USE io_mod
  USE bits_mod
  USE omp_mod
  USE trace_mod
  USE units_mod
  USE scaling_mod
  USE scalar_mod
  USE vector_mod
  USE stagger_mod
  USE gpatch_mod
  USE force_mod
  USE index_mod
  USE timestep_mod
  USE planet_mod
  USE eos_mod
  implicit none
  private
  type, extends(gpatch_t), public:: mhd_t
    real, dimension(:,:,:),   pointer:: d, dddt, s, dsdt
    real, dimension(:,:,:,:), pointer:: p, dpdt, B, dBdt
    real:: nu(6)
    real:: csound
    real:: cdtd
    logical:: mhd=.true.
    real:: gamma_AD, rho_ions
    real(8):: dt_Ohm, dt_AD, nonideal_MHD_courant
    logical:: first_time=.true., mhd_AD=.false., mhd_Ohm=.false.
    type(force_t):: force
  contains
    procedure:: init
    procedure:: pde
    procedure:: nonideal_MHD_courant_condition
    procedure:: update
    procedure:: output
    procedure:: p2U
    procedure:: U2p
    procedure:: e2E_th
    procedure:: E_th2e
    procedure:: e2S
    procedure:: S2e
    procedure:: log_density => void
    procedure:: log_pressure => void
    procedure:: gas_pressure
    procedure:: gas_pressure_ngz
    procedure:: cell_gas_pressure
    procedure:: gas_temperature
    procedure:: gas_temperature_ngz
    procedure:: cell_gas_temperature
    procedure:: velocity_magnitude => void
    procedure:: magnetic_field_magnitude => void
    procedure:: grav_potential => void
    procedure:: apply_heating
  end type
CONTAINS

!> =============================================================================
!> The construction with a double test on first_time avoids entering a critical
!> region once the variable has been set.
!> =============================================================================
SUBROUTINE init (self)
  class(mhd_t):: self
  type(index_t):: idx
  integer:: i
  real, save:: nu(6)=[0.1,1.0,0.0,0.5,0.5,0.5], gamma=1.0, csound=1.0, cdtd=0.5
  real, save:: gamma_AD=75., rho_ions=1.
  logical, save:: first_time=.true., mhd_AD=.false., mhd_Ohm=.false.
  namelist /mhd_params/ nu, gamma, csound, cdtd, mhd_AD, mhd_Ohm, gamma_AD, &
    rho_ions, hardwire
  !----------------------------------------------------------------------------
  call trace%begin('mhd_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, mhd_params)
    if (io%master) write (*, mhd_params)
    call timestep%init
  end if
  !$omp end critical (input_cr)
  self%nu = nu
  self%gamma = gamma
  self%csound = csound
  self%cdtd = cdtd
  if (self%mhd) then
    self%kind = 'stagger2_mhd_patch'
  else
    self%kind = 'stagger2_hd_patch'
  end if
  if (mhd_Ohm) then
    self%mhd_Ohm = mhd_Ohm
    self%kind = 'stagger2_mhd_Ohm_patch'
  end if
  if (mhd_AD) then
    self%mhd_AD = mhd_AD
    self%gamma_AD  = gamma_AD
    self%rho_ions  = rho_ions
    self%kind = 'stagger2_mhd_AD_patch'
  endif
  call self%gpatch_t%init
  self%i_d  = idx%d
  self%i_e  = idx%e
  self%i_s  = idx%s
  self%i_px = idx%px
  self%i_py = idx%py
  self%i_pz = idx%pz
  self%i_bx = idx%bx
  self%i_by = idx%by
  self%i_bz = idx%bz
  self%type = 'mhd_t'
  do i=1,3
    self%mesh(i)%h(2+i) = -0.5
    if (self%mhd) self%mesh(i)%h(5+i) = -0.5
  end do
  !self%unsigned(idx%d) = .true.
  call self%force%init(self%kind)
  call trace%end
END SUBROUTINE init

!> =============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(mhd_t):: self
  !.............................................................................
  real, dimension(:,:,:,:), pointer:: U, E, J, p, B, dpdt, dBdt
  real, dimension(:,:,:,:), pointer:: Ud, E_AD
  real, dimension(:,:,:,:), pointer:: fx, fy, fz, fd, ld, dd, fm
  real, dimension(:,:,:),   pointer:: d, s, dddt, dsdt
  real, dimension(:,:,:),   pointer:: du, lnd, ss, cs, pg, pa, pp, pb, Q
  real, dimension(:,:,:),   pointer:: dpdtx, dpdty, dpdtz, dBdtx, dBdty, dBdtz
  real, dimension(:,:,:),   pointer:: Bx, By, Bz, px, py, pz
  real, dimension(:,:,:),   pointer:: Ux, Uy, Uz, Ex, Ey, Ez, Jx, Jy, Jz
  real, dimension(:,:,:),   pointer:: fdx, fdy, fdz, ldx, ldy, ldz, ddx, ddy, ddz
  real, dimension(:,:,:),   pointer:: fxx, fxy
  real, dimension(:,:,:),   pointer::      fyy, fyz
  real, dimension(:,:,:),   pointer:: fzx,      fzz
  real, dimension(:,:,:),   pointer:: Ud_x, Ud_y, Ud_z, eta_non, eta_th, eta_Ohm, ambi, H_AD
  real, dimension(:,:,:),   pointer:: Txx, Txy
  real, dimension(:,:,:),   pointer::      Tyy, Tyz
  real, dimension(:,:,:),   pointer:: Tzx,      Tzz
  integer :: OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
  logical, save:: first_time=.true., first_AD=.true.
  !.............................................................................
  real(8):: ds(3), dx
  integer, save:: itimer=0
  type(index_t):: i
  !-----------------------------------------------------------------------------
  call trace%begin('mhd_t%update',itimer=itimer)
  self%d => self%mem(:,:,:,i%d      ,self%it,1); self%dddt => self%mem(:,:,:,i%d      ,self%it,2)
  self%s => self%mem(:,:,:,i%s      ,self%it,1); self%dsdt => self%mem(:,:,:,i%s      ,self%it,2)
  self%p => self%mem(:,:,:,i%px:i%pz,self%it,1); self%dpdt => self%mem(:,:,:,i%px:i%pz,self%it,2)
  d=>self%d; dddt=>self%dddt; s=>self%s; dsdt=>self%dsdt;
  p=>self%p; dpdt=>self%dpdt
  px=>p(:,:,:,1); py=>p(:,:,:,2); pz=>p(:,:,:,3);
  dpdtx=>dpdt(:,:,:,1); dpdty=>dpdt(:,:,:,2); dpdtz=>dpdt(:,:,:,3);
  if (self%mhd) then
    self%B => self%mem(:,:,:,i%bx:i%bz,self%it,1); self%dBdt => self%mem(:,:,:,i%bx:i%bz,self%it,2)
    B=>self%B; dBdt=>self%dBdt
    Bx=>B(:,:,:,1); By=>B(:,:,:,2); Bz=>B(:,:,:,3);
    dBdtx=>dBdt(:,:,:,1); dBdty=>dBdt(:,:,:,2); dBdtz=>dBdt(:,:,:,3);
  end if
  call allocate_vectors (self%gn, ld, dd, E, J, fd, fx, fy, fz, fm)
  call allocate_scalars (self%gn, lnd, pa, du, ss, cs, pp, pb, Q)
  if (self%mhd_Ohm) call allocate_scalars (self%gn, eta_Ohm)
  if (self%do_pebbles) then
    U => self%vgas(:,:,:,:,self%it)
    pg => self%pressure(:,:,:,self%it)
    self%density => self%mem(:,:,:,1,:,1)
  else
    call allocate_scalars (self%gn, pg)
    call allocate_vectors (self%gn, U)
  end if
  Ux=>U(:,:,:,1); Uy=>U(:,:,:,2); Uz=>U(:,:,:,3)
  if (self%mhd_AD) then
    call allocate_vectors (self%gn, Ud, E_AD)
    Ud_x=>Ud(:,:,:,1)
    Ud_y=>Ud(:,:,:,2)
    Ud_z=>Ud(:,:,:,3)
    call allocate_scalars (self%gn, ambi, H_AD)
  end if
  Ex=>E(:,:,:,1); Ey=>E(:,:,:,2); Ez=>E(:,:,:,3)
  Jx=>J(:,:,:,1); Jy=>J(:,:,:,2); Jz=>J(:,:,:,3)
  fxx=>fx(:,:,:,1); fxy=>fx(:,:,:,2)
                    fyy=>fy(:,:,:,2); fyz=>fy(:,:,:,3)
  fzx=>fz(:,:,:,1);                   fzz=>fz(:,:,:,3)
  fdx=>fd(:,:,:,1); fdy=>fd(:,:,:,2); fdz=>fd(:,:,:,3)
  ldx=>ld(:,:,:,1); ldy=>ld(:,:,:,2); ldz=>ld(:,:,:,3)
  ddx=>dd(:,:,:,1); ddy=>dd(:,:,:,2); ddz=>dd(:,:,:,3)
  !-----------------------------------------------------------------------------
  ! Compute velocity U=momentum/density, leaving down shifted density and log
  ! density for reuse below.  U is valid in (2:gn)
  !-----------------------------------------------------------------------------
  lnd = log(d)
  ld = down(lnd)
  dd = exp(ld)
  U = p/dd
  !-----------------------------------------------------------------------------
  ! Gas pressure. The 2nd table variable is entropy per unit mass, effectively
  ! a log variable, from the table point of view.
  !-----------------------------------------------------------------------------
  ss = s/d
  !!!call eos%lookup (shape(d), lnx=lnd, lny=ss, pg=pg)
  pg = self%gas_pressure()
  !-----------------------------------------------------------------------------
  ! Divergence of velocity, and artificial pressure.  The divergence is valid
  ! in (2,gn-1)
  !-----------------------------------------------------------------------------
  ds = self%ds
  dx = minval(ds)                                       ! smallest {dx,dy,dz}
  du = -dx*div(ds,U)                                    ! convergence rate
  du = max(du,0.0)                                      ! positive part
  !-----------------------------------------------------------------------------
  ! Conservation of momentum -- first the diagonal part of the stress tensor
  ! The diagonal fxx is valid in (2:gn-1)
  !-----------------------------------------------------------------------------
  call allocate_scalars (self%gn, Txx, Tyy, Tzz)
  Txx = ddxup(ds,Ux)
  Tyy = ddyup(ds,Uy)
  Tzz = ddzup(ds,Uz)
  if (self%mhd) then
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Diagonal (negative) Maxwell stress (borrowing E for scratch)
    Ex = xup(Bx)**2
    Ey = yup(By)**2
    Ez = zup(Bz)**2
    ! Magnetic pressure
    pb = 0.5*(Ex+Ey+Ez)
    ! Fast mode speed + velocity
    pa = self%nu(2)*d*du**2
    cs = sqrt((self%gamma*pg+2.0*pb + pa)/d)
    if (io%verbose>3) then
      self%u_max = self%fmaxval(cs)
      print*,'u_max',self%u_max
      self%u_max = self%fmaxval(cs + sqrt(pp))
      print*,'u_max',self%u_max
    end if
    self%u_max = self%fmaxval(cs + sqrt(pp) + du)
    if (io%verbose>3) print*,'u_max',self%u_max
    ! Isotropic pressure part
    pp = pg + pa + pb
    ! Total diagonal stress
    du = dx*self%nu(1)*cs
    fxx = pp + d*(fxx - du*Txx) - Ex
    fyy = pp + d*(fyy - du*Tyy) - Ey
    fzz = pp + d*(fzz - du*Tzz) - Ex
  else
    ! Diagonal Reynolds stress
    fxx = xup(Ux)**2
    fyy = yup(Uy)**2
    fzz = zup(Uz)**2
    pp = (fxx+fyy+fzz)
    ! Isotropic pressure part
    pa = self%nu(2)*d*du**2                               ! dimension: pressure
    ! Sound speed + velocity
    cs = sqrt((self%gamma*pg+pa)/d)
    if (io%verbose>3) then
      self%u_max = self%fmaxval(cs)
      print*,'u_max',self%u_max
      self%u_max = self%fmaxval(cs + sqrt(pp))
      print*,'u_max',self%u_max
    end if
    self%u_max = self%fmaxval(cs + sqrt(pp) + du)
    if (io%verbose>3) print*,'u_max',self%u_max
    pp = pg + pa
    ! Total diagondal stress, with density factor
    du = dx*self%nu(1)*cs
    !---------------------------------------------------------------------------
    ! Incoming fxx is valid in (2:gn-1), as is du and ddxup(Ux)
    !---------------------------------------------------------------------------
    fxx = pp + d*(fxx - du*Txx)
    fyy = pp + d*(fyy - du*Tyy)
    fzz = pp + d*(fzz - du*Tzz)
  end if
  !-----------------------------------------------------------------------------
    ! du has units area per unit time, which gives velocity squared times mass
    ! density per unit time = energy per unit volume and time
  !-----------------------------------------------------------------------------
  Q = d*du*(Txx**2 + Tyy**2 + Tzz**2)
  call deallocate_scalars (Txx, Tyy, Tzz)
  !-----------------------------------------------------------------------------
  ! Conservation of mass
  !-----------------------------------------------------------------------------
  fd = -(self%nu(4)*dd*down(du))*ddown(ds,lnd)          ! diffusive mass flux
  fm = p+fd                                             ! mass flux
  dddt = -div(ds,fm)                                    ! density derivative
  self%u_max = max(self%u_max,self%cdtd*dx*self%fmaxval(abs(dddt/d)))
  if (io%verbose>3) print*,'u_max',self%u_max
  !---------------------------------------------------------------------------
  ! Off-diagonal part of Reynolds stress, for now without the density factor.
  ! Note that, because of the symmmetry, we define only three components
  !---------------------------------------------------------------------------
  fxy = ydn(Ux)*xdn(Uy)
  fyz = zdn(Uy)*ydn(Uz)
  fzx = xdn(Uz)*zdn(Ux)
  !---------------------------------------------------------------------------
  ! Add viscous part -- edge-centered
  !---------------------------------------------------------------------------
  call allocate_scalars (self%gn, Txy, Tyz, Tzx)
  Txy = ddxdn(ds,Uy)+ddydn(ds,Ux)
  Tyz = ddydn(ds,Uz)+ddzdn(ds,Uy)
  Tzx = ddzdn(ds,Ux)+ddxdn(ds,Uz)
  fxy = fxy - self%nu(3)*0.5*xdn1(ydn1(du))*Txy
  fyz = fyz - self%nu(3)*0.5*ydn1(zdn1(du))*Tyz
  fzx = fzx - self%nu(3)*0.5*zdn1(xdn1(du))*Tzx
  Q = Q + self%nu(3)*0.5*d*du*(xup(yup(Txy))**2 + &
                               yup(zup(Tyz))**2 + &
                               zup(xup(Tzx))**2)
  call deallocate_scalars (Txy, Tyz, Tzx)
  !---------------------------------------------------------------------------
  ! Multiply with edge-centered density
  !---------------------------------------------------------------------------
  fxy = exp(ydn(ldx))*fxy
  fyz = exp(zdn(ldy))*fyz
  fzx = exp(xdn(ldz))*fzx
  !---------------------------------------------------------------------------
  ! Off-diagonal Maxwell stress
  !---------------------------------------------------------------------------
  if (self%mhd) then
    fxy = fxy - ydn(Bx)*xdn(By)
    fyz = fyz - zdn(By)*ydn(Bz)
    fzx = fzx - xdn(Bz)*zdn(Bx)
  end if
  !-----------------------------------------------------------------------------
  ! Consistent addition of mass flux related momentum flux, symmetrized. The
  ! diagonal part is analogous, a centered version of \nu \rho U_i \grad\ln\rho,
  ! while the off-diagonal part needs to be properly symmetrized.
  !-----------------------------------------------------------------------------
  if (self%nu(4) > 0.0) then
    pa = d*du*self%nu(4)
    !---------------------------------------------------------------------------
    ! fxx, pa, ddxup(ldx), and xup(Ux) are all valid in (2:gn-1)
    !---------------------------------------------------------------------------
    fxx = fxx - pa*ddxup(ds,ldx)*xup(Ux)
    fyy = fyy - pa*ddyup(ds,ldy)*yup(Uy)
    fzz = fzz - pa*ddzup(ds,ldz)*zup(Uz)
    fxy = fxy + (xdn1(fdy)*ydn(Ux)+ydn1(fdx)*xdn(Uy))
    fyz = fyz + (ydn1(fdz)*zdn(Uy)+zdn1(fdy)*ydn(Uz))
    fzx = fzx + (zdn1(fdx)*xdn(Uz)+xdn1(fdz)*zdn(Ux))
  end if
  !-----------------------------------------------------------------------------
  ! Momentum equations.  The ddxdn(fxx) part is valid in (3:gn-1), and hence px
  ! only needs guard zones in (1:2) and (gn:gn)
  !-----------------------------------------------------------------------------
  dpdtx = - ddxdn1(ds,fxx) - ddyup1(ds,fxy) - ddzup1(ds,fzx)
  dpdty = - ddxup1(ds,fxy) - ddydn1(ds,fyy) - ddzup1(ds,fyz)
  dpdtz = - ddxup1(ds,fzx) - ddyup1(ds,fyz) - ddzdn1(ds,fzz)
  !-----------------------------------------------------------------------------
  ! External force; e.g. point source gravity.  Note that the density factor is
  ! applied here, and should not be included in the function definition
  !-----------------------------------------------------------------------------
  if (associated(self%force%selected)) &
    dpdt = dpdt + dd*self%force%selected (self%time, d, p, Ux, Uy, Uz, self%mesh)
  !-----------------------------------------------------------------------------
  ! Conservation of entropy.  The formulation that works best so far is to make
  ! the diffusive entropy flux directly proportional to the diffusive mass flux,
  ! plus a contribution proportional to the gradient of specific entropy ss.
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    dsdt = 0.0
  else
    fm = fm*down(ss) - (self%nu(5)*dd*down(du))*ddown(ds,ss)
    !---------------------------------------------------------------------------
    ! Q is heating per unit volume (energy per unit volume and time), which when
    ! divided by presssure (energy per unit volume) becomes entropy (unitless)
    ! per unit mass.  Multiplied by mass density it becomes mass-entropy per
    ! unit volume, which is our units for s.
    !---------------------------------------------------------------------------
    dsdt = -div(ds,fm) + Q*d/pg
  end if
  !-----------------------------------------------------------------------------
  ! Induction equation
  !-----------------------------------------------------------------------------
  if (self%mhd) then
    !---------------------------------------------------------------------------
    ! Electric current, resistive electric field and magnetic dissipation
    !---------------------------------------------------------------------------
    J = curl_dn(ds,B)
    if (self%mhd_Ohm) then
      eta_Ohm = self%nu(6)
      self%dt_Ohm = 0.1*minval(dx**2/eta_Ohm)
      E = eta_Ohm*J
    else
      self%dt_Ohm = 100.
      E = self%nu(6)*edge(du)*J
    end if
    Q = yup(zup(Jx*Ex)) + zup(xup(Jy*Ey)) + xup(yup(Jz*Ez))
    if (self%gamma/=1.0) &
      dsdt = dsdt + Q*d/pg
    !---------------------------------------------------------------------------
    ! Advective electric field on cell edges
    !---------------------------------------------------------------------------
    Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)
    Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)
    Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)
    dBdt = -curl_up(ds,E)
    if (self%mhd_AD) then
      !---------------------------------------------------------------------------
      ! Below is the actual ambipolar drift parameter as given in equation 20 in
      ! Padoan, Zweibel and Nordlund (2000).
      ! To only write ambi, the definition in Padoan et al. is already multiplied
      ! by density**(-3/2) as it appears in the calculation of the drift velocity.
      ! The term d**1.5 is a cosnequence of the single fluid approximation and the
      ! relation n_i = K(n_n/10^5cm^-3)^k + K' (n_n/10^3 cm^-3)^-2
      ! K = 3*10^-3cm^-3, k=0.5 and K'=4.64*10^-4cm^-3.
      ! based on assuming that dust grains do NOT affect ionization balance significantly
      ! => second term falls off quickly => n_i~n_n^0.5
      ! In order to compute this parameter correctly, we need to know TEMPERATURE,
      ! RESOLUTION (not directly applicable here because we use AMR), BOX LENGTH,
      ! AVERAGE DENSITY, COSMIC-RAY IONISATION (K_i parameter:
      ! K_i = n_i*n_H**(-0.5) = x_i*n_H**(-0.5) ) and AVERAGE VOLUMETRIC FLOW RATE av_sigv:
      !---------------------------------------------------------------------------
      ! ambi = 0.28 * (0.1*T)**0.5 * (N/128d0) * (L/10d0)**(-1) * (av_n/200d0)**(-0.5)
      !      * (K_i/1e-5)**(-1) * (av_sigv/2e-9)**(-1) * d**(-1.5)
      !
      ! Alternative definition in Duffin and Pudritz (2008) assumes constant parameter
      ! for gamma_AD
      ! gamma_AD = av_sigv / (m_i + m_n) = 3.28e13 g^-1cm^3s^-1
      ! ambi = 1.4 / (gamma_AD*d**1.5)
      ! rho_i = 1.
      ambi = 1./(self%gamma_AD*d*self%rho_ions)
      Ud_x = xdn(ambi) * (zup(Jy)*xdn(zup(Bz)) - yup(Jz)*xdn(yup(By)))
      Ud_y = ydn(ambi) * (xup(Jz)*ydn(xup(Bx)) - zup(Jx)*ydn(zup(Bz)))
      Ud_z = zdn(ambi) * (yup(Jx)*zdn(yup(By)) - xup(Jy)*zdn(xup(Bx)))
      E_AD(:,:,:,1) = (Ud_y*Bz - Ud_z*By)
      E_AD(:,:,:,2) = (Ud_z*Bx - Ud_x*Bz)
      E_AD(:,:,:,3) = (Ud_x*By - Ud_y*Bx)
      dBdt(:,:,:,1) = dBdt(:,:,:,1) + (ddyup(ds, (Ud_x*By - Ud_y*Bx)) - ddzup(ds, (Ud_z*Bx - Ud_x*Bz)) )
      dBdt(:,:,:,2) = dBdt(:,:,:,2) + (ddzup(ds, (Ud_y*Bz - Ud_z*By)) - ddxup(ds, (Ud_x*By - Ud_y*Bx)) )
      dBdt(:,:,:,3) = dBdt(:,:,:,3) + (ddxup(ds, (Ud_z*Bx - Ud_x*Bz)) - ddyup(ds, (Ud_y*Bz - Ud_z*By)) )
      self%dt_AD = 0.1*minval(self%gamma_AD*self%rho_ions*d/(Bx**2+By**2+Bz**2)*dx**2)
      if (self%gamma==1.0) then
        dsdt = 0.
      else
        ! ambipolar heating H_AD = rho_i*rho_n*nu_in*v_D^2. nu_in is the collisional frequency between ions and neutrals
        ! and so far this is all hidden in the constant gamma_AD.
        H_AD = self%rho_ions * d * self%gamma_AD * (xup(Ud_x)**2 + yup(Ud_y)**2 + zup(Ud_z)**2)
        dsdt = dsdt + H_AD*d/pg
      end if
    end if
  end if
  !-----------------------------------------------------------------------------
  call deallocate_vectors (ld, dd, E, J, fd, fx, fy, fz, fm)
  call deallocate_scalars (lnd, pa, du, ss, cs, pp, pb, Q)
  if (self%mhd_Ohm) then
    call deallocate_scalars (eta_Ohm)
  end if
  if (self%mhd_AD) then
    call deallocate_vectors (Ud, E_AD)
    call deallocate_scalars (ambi,H_AD)
  end if
  if (.not.self%do_pebbles) then
    call deallocate_vectors (U)
    call deallocate_scalars (pg)
  end if
  call trace%end (itimer)
END SUBROUTINE pde

!> =============================================================================
!> Courant condition, one way or another
!> =============================================================================
SUBROUTINE nonideal_MHD_courant_condition (self,dt_Ohm_inp,dt_AD_inp)
  class(mhd_t):: self
  real(8):: dt_AD_inp
  real(8):: dt_Ohm_inp
  !.............................................................................
  call trace%begin('mhd_t%nonideal_MHD_courant')
  call self%lock%set ('nonideal_MHD_courant')
  if (io%verbose>2) &
    print*, 'nonideal_MHD_courant: u_max, ds', self%u_max, minval(self%ds)
  if (self%dt_fixed > 0.0) then
    self%dtime = self%dt_fixed
    self%nonideal_MHD_courant = self%dtime*self%u_max/minval(self%ds)
  else if (self%u_fixed > 0.0) then
    self%dtime = self%nonideal_MHD_courant*minval(self%ds)/self%u_fixed
  else
    self%dtime = self%nonideal_MHD_courant*minval(self%ds)/self%u_max
  endif
  !-----------------------------------------------------------------------------
  ! Prevent large (cheap) patches from taking time steps longer than the out_time
  !-----------------------------------------------------------------------------
  if (self%out_time > 0.0) &
    self%dtime = min (self%dtime, self%out_time)
  if (self%mhd_AD) then
    self%dtime = min (self%dtime, dt_AD_inp)
  end if
  if (self%mhd_Ohm) then
    self%dtime = min (self%dtime, dt_Ohm_inp)
  end if
  self%dt(self%it) = self%dtime
  !self%dt(:) = self%dtime
  !print*,'dt', self%dt
  if (io%verbose>3) &
    print'(a,i7,f12.6,f8.4,1p,e10.2)', 'nonideal_MHD_courant:', &
      self%id, self%dtime, self%nonideal_MHD_courant, self%u_max
  call self%lock%unset ('nonideal_MHD_courant')
  call trace%end
END SUBROUTINE nonideal_MHD_courant_condition

!> =============================================================================
!> Take a complete update step.  This can be overloaded in higher levels, e.g.
!> in experiment_mod.
!> =============================================================================
SUBROUTINE update (self)
  class(mhd_t):: self
  !-----------------------------------------------------------------------------
  ! The call to output below catches the state after both guard zone downloads
  ! and boundary conditions.  The call here is a fall back; if another call is
  ! placed at an earlier point in the task handling the call here will have no
  ! effect, since out_next has already been updated.
  !-----------------------------------------------------------------------------
  call self%output
  call self%pde
  call self%courant_condition
  if (io%verbose>2) &
    call self%check_density('before timestep')
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
  call self%check_density(' after timestep')
END SUBROUTINE update

!> =============================================================================
!> Add extra, MHD-related information to the patch .txt output file
!> =============================================================================
SUBROUTINE output (self, experiment_name)
  class(mhd_t):: self
  character(len=64):: experiment_name, filename
  optional:: experiment_name
  logical:: exists
  !-----------------------------------------------------------------------------
  if (self%time < self%out_next .or. .not.io%do_output .or. .not.io%do_legacy) return
  call self%gpatch_t%output (experiment_name)
  !$omp critical (out_cr)
  write (filename,'(a,i4.4,"_",i4.4,".txt")') trim(io%outputname), &
    self%id, self%iout
  inquire(file=filename, exist=exists)
  if (exists) then
    open (io%data_unit, file=trim(filename), form='formatted', status='old', &
      position='append')
    write (io%data_unit,'(8f10.4)') self%nu
    close (io%data_unit)
  end if
  !$omp end critical (out_cr)
END SUBROUTINE output

!> =====================================================================
!> Get velocities from momenta
!> =====================================================================
SUBROUTINE p2U(self, U, it)
  class(mhd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,           it,1), &
            p => self%mem(:,:,:,self%i_px:self%i_pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  U = p/exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE p2U

!> =====================================================================
!> Put velocities to momenta
!> =====================================================================
SUBROUTINE U2p(self, U, it)
  class(mhd_t) :: self
  real, dimension(:,:,:,:), pointer:: U
  integer:: it
  !.....................................................................
  real, dimension(:,:,:,:), pointer:: dd
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d           ,it,1), &
            p => self%mem(:,:,:,self%i_px:self%i_pz,it,1))
  allocate(dd(size(d,1),size(d,2),size(d,3),3))
  dd=down(log(d))
  p = U*exp(dd)
  deallocate(dd)
  end associate
END SUBROUTINE U2p

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE e2E_th(self, E_th, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            s => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  E_th = d**g1*exp(s*g1/d)/g1
  end associate
END SUBROUTINE e2E_th

!> =====================================================================
!> Put thermal energy to entropy per unit volume
!> =====================================================================
SUBROUTINE E_th2e(self, E_th, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: E_th
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            s => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  s = d*log(E_th*g1/d**g1)/g1
  end associate
END SUBROUTINE E_th2e

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE e2S(self, s, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  s = log(e*g1/d**self%gamma)/g1
  end associate
END SUBROUTINE e2S

!> =====================================================================
!> Get thermal energy per unit mass from entropy
!> =====================================================================
SUBROUTINE S2e(self, s, it)
  class(mhd_t) :: self
  integer:: it
  !.....................................................................
  real, dimension(:,:,:), pointer:: s
  real:: g1
  !---------------------------------------------------------------------
  associate(d => self%mem(:,:,:,self%i_d,it,1), &
            e => self%mem(:,:,:,self%i_s,it,1))
  g1 = self%gamma-1.0
  e = d**self%gamma*exp(s*g1)/g1
  end associate
END SUBROUTINE S2e

!> =============================================================================
FUNCTION gas_pressure (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: pg
  !-----------------------------------------------------------------------------
  if (self%gamma==1.0) then
    pg = self%d*self%csound**2
  else
    pg = self%d**self%gamma*exp(self%s/self%d*(self%gamma-1))
  end if
END FUNCTION gas_pressure

FUNCTION gas_pressure_ngz (self) RESULT (pg)
  class(mhd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: pg
  integer :: l(3), u(3)
  !-----------------------------------------------------------------------------
  l = self%li; u = self%ui
  associate (d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), 1,self%it,1), &
             s => self%mem(l(1):u(1),l(2):u(2),l(3):u(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d
  else
    pg = d**self%gamma*exp(s/d*(self%gamma-1))
  end if
  end associate
END FUNCTION gas_pressure_ngz

FUNCTION cell_gas_pressure(self, c) RESULT (pg)
  class(mhd_t) :: self
  real :: pg
  integer :: c(3)
  !---------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3), 1,self%it,1), &
             s => self%mem(c(1),c(2),c(3), 2,self%it,1))
  if (self%gamma==1.0) then
    pg = d
  else
    pg = d**self%gamma*exp(s/d*(self%gamma-1))
  end if
  end associate
END FUNCTION cell_gas_pressure

!> =============================================================================
!> temperature from ideal gas law, taking into account that the temperature is
!> normalised by k_B and mu (so that specific gas constant = 1)
!> =============================================================================
FUNCTION gas_temperature (self) RESULT (tmp)
  class(mhd_t):: self
  real, dimension(self%gn(1),self%gn(2),self%gn(3)):: tmp
  !-----------------------------------------------------------------------------
  associate (d => self%mem(:,:,:,  1,self%it,1))
  tmp = gas_pressure(self)/d

  if (io%verbose > 2) then
    if (any(tmp<0)) then
      error stop "mhd_t%gas_temperature: T<0"
    end if
  end if
  end associate
END FUNCTION gas_temperature

FUNCTION gas_temperature_ngz (self) RESULT (tmp)
  class(mhd_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: tmp
  integer :: l(3), u(3)
  !---------------------------------------------------------------------
  l = self%li; u = self%ui
  associate (d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),  1,self%it,1))
  tmp = self%gas_pressure_ngz()/d

  if (io%verbose > 2) then
    if (any(tmp<0)) then
      error stop "mhd_t%gas_temperature: T<0"
    end if
  end if
  end associate
END FUNCTION gas_temperature_ngz

FUNCTION cell_gas_temperature (self, c) RESULT (tmp)
  class(mhd_t):: self
  real :: tmp
  integer :: c(3)
  !-----------------------------------------------------------------------------
  associate (d => self%mem(c(1),c(2),c(3),  1,self%it,1))
  tmp = cell_gas_pressure(self,c)/d
  end associate
END FUNCTION cell_gas_temperature

!> =============================================================================
SUBROUTINE log_density (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(self%d)
END SUBROUTINE

!> =============================================================================
SUBROUTINE log_pressure (self, v)
  class(mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = log(gas_pressure(self))
END SUBROUTINE

!> =============================================================================
SUBROUTINE velocity_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%p/exp(down(log(self%d))))
END SUBROUTINE

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = norm(self%B)
END SUBROUTINE

!> =============================================================================
SUBROUTINE void (self, v)
  class (mhd_t):: self
  real, dimension(:,:,:), pointer:: v
  !-----------------------------------------------------------------------------
  v = 0.0
END SUBROUTINE

SUBROUTINE apply_heating (self, q)
  class (mhd_t):: self
  real, dimension(:,:,:):: q
END SUBROUTINE

END MODULE mhd_mod
