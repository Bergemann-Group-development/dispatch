! *******************************************************************
MODULE atom_mod
    USE m3io_mod
    USE mconst
    USE coll_mod
    USE trans_mod
    USE mmath
    USE op_math_mod, only: searchsorted
    USE os_mod
    USE trace_mod
    USE spectrum_mod, only: spectrum_type
    ! *******************************************************************
    !
    !     Contains variables connected to the atom(s) 
    !     included in the calculations
    !
    ! *******************************************************************
    implicit none

    !-------------------------------------------------------------------------------
    type, extends(spectrum_type) :: atomic_spectrum_type
        ! maximum number of active lines or continua at any given frequency
        integer :: maxal,maxac
        ! which continua (ac) and lines (al) are active at each frequency
        integer,dimension(:,:),  allocatable :: ac,al
        ! how many continua (nac) and lines (nal) are active at each frequency
        integer,dimension(:),  allocatable :: nac,nal
        ! has this frequency an active prd line
        !logical,dimension(:),allocatable :: hasPRDline
    end type atomic_spectrum_type
    !---------------------------------------------------------------------
  
    type atom_type
  
       integer :: nrad   ! total number of radiative transitions
       integer :: nrfix  ! number of fixed transitions
       integer :: ncont  ! number of bound-free transitions
       integer :: nline  ! number of bound-bound transitions
       integer :: nlevel ! number of levels
       integer :: nprd=0 ! number of PRD transitions
       
       character(len=20) :: id    ! atom id string 
       character(len=20) :: crout ! which collisional SUBROUTINE to use 
  
       character(len=20),dimension(:),allocatable :: label ! level label
  
       integer,dimension(:),allocatable :: ion ! ionization stage of each level
  
       integer,dimension(:,:), allocatable :: ilin,icon ! ilin(l1,l2) is the index of the line
                                                          ! between l1 and l2, if any, 
                                                          ! likewise for icon
  
       real(8) :: abnd ! abundance ( in the [hydrogen]=12 scale)
       real(8) :: awgt ! atomic weight [g]
  
       real(8),dimension(:),allocatable :: ev ! energy of level above ground state [eV]
       real(8),dimension(:),allocatable :: g ! multiplicity of level
  
        type(line_type), dimension(:), allocatable :: line
        type(cont_type), dimension(:), allocatable :: cont
        type(coll_type), dimension(:), allocatable :: coll
        type(atomic_spectrum_type), pointer :: sp
  
        ! rather constant doppler width factor, i.e. dnud = line(i)%nu0 * doppler_factor [Hz]
        real(4), dimension(:,:,:), pointer, contiguous :: doppler_factor
  
        ! LTE, NLTE and total population of active atom
        real(4), dimension(:,:,:,:), pointer, contiguous :: nstar=>null()
        real(4), dimension(:,:,:,:), pointer, contiguous :: n    =>null()
        real(4), dimension(:,:,:), pointer, contiguous :: totn
  
        real(8), dimension(:)  , allocatable  :: e
        real(8), dimension(:,:), allocatable  :: C, W
        
        real(8) :: qnorm
        integer :: fcheck, pcheck(3), linecheck
     contains
        procedure :: add_line_opacities
        procedure :: add_continuum_opacities
        procedure :: add_all_gammas
        procedure :: add_all_gammas2
        procedure :: compute_line_damping
        procedure :: compute_doppler_factor
        procedure :: partition_function
        procedure :: lte_pop
        procedure :: compute_lte_populations
        procedure :: compute_cont_factors
        procedure :: compute_wphi
        procedure :: get_trans
        procedure :: solve_stat_equilibrium
        procedure :: create_spectrum
        procedure :: set_active_transitions
        procedure :: read
     
    end type atom_type
  
    ! *******************************************************************
  
    !type(atom_type),target :: atom
    !type(line_type),dimension(:),allocatable,target :: line
    !type(cont_type),dimension(:),allocatable,target :: cont
  
  contains
  
     !====================================================================================
     SUBROUTINE add_continuum_opacities(self, fs_patch, fe_patch, nu, temp, chi, eta)
        class(atom_type), intent(in), target :: self
        class(cont_type), pointer :: cnt
        integer, intent(in) :: fs_patch, fe_patch
        real(8), dimension(:), intent(in) :: nu
        real(4), dimension(:,:,:), intent(in) :: temp
        real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
        real(8), dimension(:), allocatable :: twohnu3_c2, hnu_k
        real(4), dimension(:,:,:), allocatable :: t_inv
        real(4) :: blz
        integer :: i, j, k, f, fcs, fce, fs, fe, nx, ny, nz, ffs, ffe, ff, p(3)
        logical :: skip
  
        ! all_chi and all_eta only range from fs:fe
        ! nu however contains all frequencies from atom
  
        fcs = max(minval(self%cont%ired),  fs_patch)
        fce = min(maxval(self%cont%iblue), fe_patch)
  
        if (fcs >= fe_patch) return
        if (fce <= fs_patch) return
  
        nx = size(temp,1)
        ny = size(temp,2)
        nz = size(temp,3)
  
        allocate(twohnu3_c2(fcs:fce), hnu_k(fcs:fce))
        allocate(t_inv(nx,ny,nz))
  
        twohnu3_c2 = twoh_c2 * nu(fcs:fce)**3
        hnu_k = h_k * nu(fcs:fce)
        t_inv = 1 / temp
  
        do i=1, size(self%cont)
           cnt => self%cont(i)
           call cnt%check_range(fcs,fce,fs,fe,skip)
           if (skip) cycle
        
           ffs = fs-fs_patch+1
           ffe = fe-fs_patch+1
  
           call cnt%opacities(cnt%alpha(fs:fe), self%n, self%nstar, eta(:,:,:,ffs:ffe), chi(:,:,:,ffs:ffe))
  
           !if (self%fcheck >= fs .and. self%fcheck <= fe) then
           !   p(:) = self%pcheck
           !   f = self%fcheck
           !   print *, i, cnt%i, cnt%j, 'chi cont'
           !   print '(8(ES15.7))', chi(p(1),p(2),p(3),f-fs_patch+1) &
           !         - eta(p(1),p(2),p(3),f-fs_patch+1) &
           !         * exp(- hnu_k(f) * t_inv(p(1),p(2),p(3))), &
           !            cnt%alpha(f) * self%n(p(1),p(2),p(3), cnt%i) &
           !         * (1 - exp(- hnu_k(f) * t_inv(p(1),p(2),p(3)))), &
           !         self%n(p(1),p(2),p(3),cnt%i), &
           !         self%n(p(1),p(2),p(3),cnt%j)
           !endif
        end do
  
        do f=fcs,fce
           ff = f-fs_patch+1
           do k=1,nz
           do j=1,ny
           !$omp simd
           do i=1,nx
              blz = - hnu_k(f) * t_inv(i,j,k)
              if (blz < -80) then
                 eta(i,j,k,ff) = 0
              else
                 if (blz < -1e-8) then
                    eta(i,j,k,ff) = eta(i,j,k,ff) * exp(blz)
                 endif
                 chi(i,j,k,ff) = chi(i,j,k,ff) - eta(i,j,k,ff)
                 eta(i,j,k,ff) = eta(i,j,k,ff) * twohnu3_c2(f)
              endif
           end do
           end do
           end do
        end do
     END SUBROUTINE add_continuum_opacities
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE add_line_opacities(self, fs_patch, fe_patch, nu, temp, vproj, chi, eta)
        class(atom_type), intent(in), target :: self
        integer, intent(in) :: fs_patch, fe_patch
        real(8), dimension(:), intent(in) :: nu
        real(4), dimension(:,:,:), intent(in) :: temp, vproj
        real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
        real(4), dimension(:,:,:,:), allocatable, target :: scratch_phi
        real(4), dimension(:,:,:,:), pointer :: phi
        real(4), dimension(:,:,:), allocatable :: vel_dop
        real(8), dimension(:), allocatable :: hnu_4pi
        class(line_type), pointer :: line
        integer :: i, fcs, fce, fs, fe, nx, ny, nz, ffs, ffe, f, p(3)
        logical :: skip
  
        ! all_chi and all_eta only range from fs:fe
        ! nu however contains all frequencies from atom
  
        fcs = max(minval(self%line%ired),  fs_patch)
        fce = min(maxval(self%line%iblue), fe_patch)
  
        if (fcs >= fe_patch) return
        if (fce <= fs_patch) return
  
        nx = size(temp,1)
        ny = size(temp,2)
        nz = size(temp,3)
  
        allocate(hnu_4pi(fcs:fce))
        allocate(scratch_phi(nx,ny,nz,fcs:fce))
  
        hnu_4pi = h_4pi * nu(fcs:fce)
        vel_dop = (1 + vproj) / self%doppler_factor
  
        do i=1, self%nline
           line => self%line(i)
           call line%check_range(fcs,fce,fs,fe,skip)
           if (skip) cycle
        
           ffs = fs-fs_patch+1
           ffe = fe-fs_patch+1
  
           phi => scratch_phi(:,:,:,fs:fe)
  
           call line%opacities(nu(fs:fe), hnu_4pi(fs:fe), self%n, vel_dop, self%doppler_factor, phi, eta(:,:,:,ffs:ffe), chi(:,:,:,ffs:ffe))
  
           !if (self%fcheck >= fs .and. self%fcheck <= fe) then
           !   p(:) = self%pcheck
           !   f = self%fcheck
           !   print *, i, line%i, line%j, 'chi line'
           !   print '(8(ES15.7))', chi(p(1),p(2),p(3),f-fs_patch+1), &
           !             scratch_phi(p(1),p(2),p(3), f) &
           !          * (line%Bij * self%n(p(1),p(2),p(3),line%i) &
           !          -  line%Bji * self%n(p(1),p(2),p(3),line%j)), &
           !          scratch_phi(p(1),p(2),p(3), f), &
           !          self%n(p(1),p(2),p(3),line%i), &
           !          self%n(p(1),p(2),p(3),line%j), nu(f)
           !endif
        end do
     END SUBROUTINE add_line_opacities
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE compute_wphi (self, fs_patch, fe_patch, nu, temp, vproj, wphi)
        class(atom_type), intent(in), target :: self
        integer, intent(in) :: fs_patch, fe_patch
        real(8), dimension(:), intent(in) :: nu
        real(4), dimension(:,:,:), intent(in) :: temp, vproj
        real(4), dimension(:,:,:), allocatable :: vel_dop
        real(4), dimension(:,:,:,:), allocatable, target :: scratch_phi
        real(4), dimension(:,:,:,:), pointer :: phi
        real(4), dimension(:,:,:,:), pointer, contiguous, intent(out) :: wphi
        class(line_type), pointer :: line
        integer :: i, fcs, fce, fs, fe, nx, ny, nz, f
        logical :: skip
  
        nx = size(temp,1)
        ny = size(temp,2)
        nz = size(temp,3)
  
        call os%allocate(wphi, [nx,ny,nz,self%nline])
        wphi = 0
        
        fcs = max(minval(self%line%ired),  fs_patch)
        fce = min(maxval(self%line%iblue), fe_patch)
  
        if (fcs >= fe_patch) return
        if (fce <= fs_patch) return
  
        allocate(scratch_phi(nx,ny,nz,fcs:fce))
  
        vel_dop = (1 + vproj) / self%doppler_factor
  
        do i=1, self%nline
           line => self%line(i)
           call line%check_range(fcs,fce,fs,fe,skip)
           if (skip) cycle
           
           phi => scratch_phi(:,:,:,fs:fe)
  
           call line%profile(nu(fs:fe), vel_dop, self%doppler_factor, line%adamp, phi)
  
           do f=fs,fe
              scratch_phi(:,:,:,f) = scratch_phi(:,:,:,f) * line%wnu(f)
           end do
  
           wphi(:,:,:,i) = sum(phi,4) ! angle weight is multiplied later
        end do
  
        deallocate(scratch_phi)
     END SUBROUTINE
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE add_all_gammas(self, ff, temp, vel_dop, Inu, Psi, chi_l, chi_c, eta_l, eta_c, gamma_ij, gamma_ji)
        class(atom_type), target, intent(in) :: self
        integer, intent(in) :: ff
        real(4), dimension(:,:,:), intent(in) :: temp, vel_dop
        real(4), dimension(:,:,:), intent(in) :: Inu, Psi, chi_l, chi_c, eta_l, eta_c
        real(4), dimension(:,:,:,:), intent(inout), target :: gamma_ij, gamma_ji
        real(4), dimension(:,:,:), pointer :: Gij, Gji
        real(4), dimension(:,:,:), allocatable :: blzfac, vij, vji, eta, chi, Ieff
        real(4), dimension(:,:,:,:), allocatable :: phi
        real(4), dimension(:,:,:,:), allocatable :: SumU, SumChiDn, SumChiUp
        real(8), dimension(:), allocatable :: wnu, fac
        integer, dimension(:), allocatable :: lower, upper, trans
        real(8) :: twohnu3_c2, hnu_inv
        real(8) :: hnu_4pi
        integer :: i, j, f, kr, kr2, kc, kl, nx, ny, nz, n_active, p(3)
        logical :: recompute_opac, cross_coupling=.false.
        !------------------------------------------------------------------------------------
        nx = size(temp,1)
        ny = size(temp,2)
        nz = size(temp,3)
        n_active = self%sp%nac(ff) + self%sp%nal(ff)
  
        allocate(phi(nx,ny,nz,1))
        allocate(SumU(nx,ny,nz,self%nlevel))
        allocate(SumChiDn(nx,ny,nz,self%nlevel))
        allocate(SumChiUp(nx,ny,nz,self%nlevel))
        allocate(blzfac, vij, vji, eta, chi, Ieff, mold=temp)
        allocate(wnu(n_active), fac(n_active))
        allocate(lower(n_active), upper(n_active), trans(n_active))
  
        SumU     = 0
        SumChiUp = 0
        SumChiDn = 0
  
        hnu_inv = 1 / (hh*self%sp%nu(ff))
        hnu_4pi = h_4pi * self%sp%nu(ff)
        
        Ieff = Inu - Psi * (eta_c + eta_l)
  
        if (self%sp%nac(ff) > 0) then
           twohnu3_c2 = twoh_c2 * self%sp%nu(ff)**3
           blzfac = exp(- h_k * self%sp%nu(ff) / temp)
        endif
        
        do kr=1, n_active
           if (kr <= self%sp%nac(ff)) then
              kc = self%sp%ac(kr,ff)
              lower(kr) = self%cont(kc)%i
              upper(kr) = self%cont(kc)%j
              trans(kr) = self%cont(kc)%ntrans
              fac  (kr) = twohnu3_c2
              wnu  (kr) = self%cont(kc)%wnu(ff) * hnu_inv ! wm needs to be multiplied later
           else
              kl = self%sp%al(kr - self%sp%nac(ff),ff)
              lower(kr) = self%line(kl)%i
              upper(kr) = self%line(kl)%j
              trans(kr) = self%line(kl)%ntrans
              fac  (kr) = self%line(kl)%Aji / self%line(kl)%Bji
              wnu  (kr) = self%line(kl)%wnu(ff) * hnu_inv ! wm needs to be multiplied later
           endif
        end do
  
        do kr=1, n_active
           do kr2=1, n_active
              if (lower(kr) == upper(kr2)) then
                 cross_coupling=.true.
              endif
           end do
        end do
  
        recompute_opac = self%sp%nac(ff) > 1 .or. self%sp%nal(ff) > 1
  
        do kr=1, n_active
           i = lower(kr)
           j = upper(kr)
  
           if (kr <= self%sp%nac(ff)) then
              if (recompute_opac) then
                 kc = self%sp%ac(kr,ff)
                 vij = self%cont(kc)%alpha(ff)
                 vji = self%nstar(:,:,:,i) / self%nstar(:,:,:,j) * blzfac * vij
              else
                 chi = chi_c
                 eta = eta_c
              endif
           else
              if (recompute_opac) then
                 kl = self%sp%al(kr - self%sp%nac(ff),ff)
                 call self%line(kl)%profile(self%sp%nu(ff:ff), vel_dop, self%doppler_factor, self%line(kl)%adamp, phi)
                 phi(:,:,:,1) = phi(:,:,:,1) * hnu_4pi
                 vij = self%line(kl)%Bij * phi(:,:,:,1)
                 vji = self%line(kl)%Bji * phi(:,:,:,1)
              else
                 chi = chi_l
                 eta = eta_l
              endif
           endif
  
           if (recompute_opac) then
              chi = vij * self%n(:,:,:,i) - vji * self%n(:,:,:,j)
           else
              vji = eta / (fac(kr) * self%n(:,:,:,j))
              vij = (chi + vji * self%n(:,:,:,j)) / self%n(:,:,:,i)
           endif
  
           SumU    (:,:,:,j) = SumU    (:,:,:,j) + fac(kr) * vji
           SumChiUp(:,:,:,i) = SumChiUp(:,:,:,i) + chi
           
           if (cross_coupling) then
              SumChiDn(:,:,:,j) = SumChiDn(:,:,:,j) + chi
           endif
  
           Gij => gamma_ij(:,:,:,trans(kr))
           Gji => gamma_ji(:,:,:,trans(kr))
  
           ! add to upward rate coefficient Gamma_ij   
           Gij = Gij + wnu(kr) * Ieff * vij
  
           ! add to downward rate coefficient Cont(Kc)%gamma_ji 
           Gji = Gji + wnu(kr) * (fac(kr) + Ieff) * vji
  
           !if (ff==self%fcheck) then
           !   p(:) = self%pcheck
           !   f = self%fcheck
           !   print *, trans(kr), i, j, 'rate'
           !   print '(8(ES15.7))', wnu(kr), vij(p(1),p(2),p(3)), vji(p(1),p(2),p(3))
           !endif
        end do
  
        do kr=1, n_active
           i = lower(kr)
           j = upper(kr)
           Gij => gamma_ij(:,:,:,trans(kr))
           Gji => gamma_ji(:,:,:,trans(kr))
  
           Gji = Gji - wnu(kr) * SumChiUp(:,:,:,i) * Psi * SumU(:,:,:,j)
  
           if (cross_coupling) then
              Gij = Gij + wnu(kr) * SumChiDn(:,:,:,j) * Psi * SumU(:,:,:,i)
           endif
  
           !if (ff==self%fcheck) then
           !   p(:) = self%pcheck
           !   f = self%fcheck
           !   print *, trans(kr), i, j, 'cross coupling'
           !   print '(8(ES15.7))', SumChiUp(p(1),p(2),p(3),i), SumChiDn(p(1),p(2),p(3),j), SumU(p(1),p(2),p(3),j), &
           !      Ieff(p(1),p(2),p(3)), eta_c(p(1),p(2),p(3)) + eta_l(p(1),p(2),p(3))
           !endif
        end do
     END SUBROUTINE add_all_gammas
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE add_all_gammas2(self, fs, temp, Inu, all_Psi, gamma_ij, gamma_ji, lines, idx_map)
        class(atom_type), target, intent(in) :: self
        integer, intent(in) :: fs
        real(4), dimension(:,:,:), intent(in) :: temp
        real(4), dimension(:,:,:,:), intent(in) :: Inu
        real(4), dimension(:,:,:,:), intent(in), target :: all_Psi
        real(4), dimension(:,:,:,:), intent(inout), target :: gamma_ij, gamma_ji
        class(line_type), dimension(:), intent(in) :: lines
        integer, dimension(:), intent(in) :: idx_map
        real(4), dimension(:,:,:), allocatable :: blzfac, chi, Ieff, SumEta, hk_over_temp
        real(4), dimension(:,:,:), pointer :: vij, vji, Psi
        real(4), dimension(:,:,:,:), allocatable :: SumU, SumChiDn, SumChiUp, phi
        real(4), dimension(:,:,:,:), allocatable, target :: Gij, Gji
        real(8), dimension(:), allocatable :: wnu, fac
        integer, dimension(:), allocatable :: lower, upper, trans
        logical, dimension(:), allocatable :: cross, up_done, lo_done
        real(8) :: twohnu3_c2, hnu_inv, ufac, Bij, Bji, alpha
        real(4) :: U
        integer :: i, j, k, kr, kr2, kc, kl, fe, n_active, ff, f, p(3), x,y,z, nx, ny, nz
        logical :: recompute_opac, cross_coupling=.false.
        integer, save :: itimer1, itimer2, itimer3, itimer4, itimer5, itimer6
        !------------------------------------------------------------------------------------
        fe = fs + size(Inu, 4) - 1 
        nx = size(Inu,1)
        ny = size(Inu,2)
        nz = size(Inu,3)
  
        allocate(phi, mold=Inu(:,:,:,1:1))
        allocate(SumU, mold=self%n)
        allocate(SumChiUp, mold=self%n)
        allocate(SumChiDn, mold=self%n)
        allocate(blzfac, vij, vji, chi, Ieff, SumEta, hk_over_temp, mold=temp)
        allocate(Gij, mold=gamma_ij)
        allocate(Gji, mold=gamma_ji)
  
        n_active = maxval(self%sp%nac(fs:fe) + self%sp%nal(fs:fe))
        allocate(wnu(n_active), fac(n_active))
        allocate(lower(n_active), upper(n_active), trans(n_active), cross(n_active))
        allocate(up_done(self%nlevel), lo_done(self%nlevel))
  
        hk_over_temp = h_k / temp
  
        do ff=fs,fe
           call trace%begin('administration', itimer=itimer2)
           f = ff - fs + 1
           n_active = self%sp%nac(ff) + self%sp%nal(ff)
  
           hnu_inv = 1 / (hh*self%sp%nu(ff))
  
           if (self%sp%nac(ff) > 0) then
              twohnu3_c2 = twoh_c2 * self%sp%nu(ff)**3
              blzfac = exp(- hk_over_temp * self%sp%nu(ff))
           endif
           
           do kr=1, n_active
              if (kr <= self%sp%nac(ff)) then
                 kc = self%sp%ac(kr,ff)
                 lower(kr) = self%cont(kc)%i
                 upper(kr) = self%cont(kc)%j
                 trans(kr) = self%cont(kc)%ntrans
                 fac  (kr) = twohnu3_c2
                 wnu  (kr) = self%cont(kc)%wnu(ff) * hnu_inv ! wm needs to be multiplied later
              else
                 kl = self%sp%al(kr - self%sp%nac(ff),ff)
                 lower(kr) = self%line(kl)%i
                 upper(kr) = self%line(kl)%j
                 trans(kr) = self%line(kl)%ntrans
                 fac  (kr) = self%line(kl)%Aji / self%line(kl)%Bji
                 wnu  (kr) = self%line(kl)%wnu(ff) * hnu_inv ! wm needs to be multiplied later
              endif
           end do
  
           do kr=1, n_active
              do kr2=1, n_active
                 if (lower(kr) == upper(kr2)) then
                    cross_coupling = .true.
                    cross(kr2) = .true.
                 endif
              end do
           end do
  
           SumEta = 0
  
           call zero_sums(SumU, upper, up_done, n_active)
           call zero_sums(SumChiUp, lower, lo_done, n_active)
  
           if (cross_coupling) then
              call zero_sums(SumU, lower, lo_done, n_active)
              call zero_sums(SumChiDn, upper, up_done, n_active)
           endif
  
           call trace%end(itimer2)
           call trace%begin('vij', itimer=itimer3)
  
           do kr=1, n_active
              i = lower(kr)
              j = upper(kr)
              k = trans(kr)
  
              vij => Gij(:,:,:,k)
              vji => Gji(:,:,:,k)
  
              if (kr <= self%sp%nac(ff)) then
                 kc = self%sp%ac(kr,ff)
                 alpha = self%cont(kc)%alpha(ff) * wnu(kr)
                 ufac = fac(kr) / wnu(kr)
                 do z=1,nz
                 do y=1,ny
                 do x=1,nx
                    vij(x,y,z) = alpha
                    vji(x,y,z) = self%cont(kc)%ni_over_nj_lte(x,y,z) * blzfac(x,y,z) * alpha
                    chi(x,y,z) = vij(x,y,z) * self%n(x,y,z,i) - vji(x,y,z) * self%n(x,y,z,j)
                    U = vji(x,y,z) * ufac
  
                    SumU    (x,y,z,j) = SumU    (x,y,z,j) + U
                    SumChiUp(x,y,z,i) = SumChiUp(x,y,z,i) + chi(x,y,z)
     
                    SumEta(x,y,z) = SumEta(x,y,z) + U * self%n(x,y,z,j)
                 end do
                 end do
                 end do
              else
                 kl = self%sp%al(kr - self%sp%nac(ff),ff)
                 kl = idx_map(kl)
                 call lines(kl)%interpolate_phi(self%sp%nu(ff:ff), phi)
                 Bij = lines(kl)%Bij * wnu(kr)
                 Bji = lines(kl)%Bji * wnu(kr)
                 ufac = fac(kr) / wnu(kr)
                 do z=1,nz
                 do y=1,ny
                 do x=1,nx
                    vij(x,y,z) = Bij * phi(x,y,z,1)
                    vji(x,y,z) = Bji * phi(x,y,z,1)
                    chi(x,y,z) = vij(x,y,z) * self%n(x,y,z,i) - vji(x,y,z) * self%n(x,y,z,j)
                    U = vji(x,y,z) * ufac
  
                    SumU    (x,y,z,j) = SumU    (x,y,z,j) + U
                    SumChiUp(x,y,z,i) = SumChiUp(x,y,z,i) + chi(x,y,z)
     
                    SumEta(x,y,z) = SumEta(x,y,z) + U * self%n(x,y,z,j)
                 end do
                 end do
                 end do
              endif
                 
              if (cross(kr)) then
                 SumChiDn(:,:,:,j) = SumChiDn(:,:,:,j) + chi
              endif
  
              !if (ff==self%fcheck) then
              !   p(:) = self%pcheck
              !   f = self%fcheck
              !   print *, k, i, j, 'rate'
              !   print '(8(ES15.7))', wnu(kr), vij(p(1),p(2),p(3)), vji(p(1),p(2),p(3))
              !endif
           end do
  
           call trace%end(itimer3)
           call trace%begin('Gij', itimer=itimer4)
  
           Psi => all_Psi(:,:,:,f)
           Ieff = Inu(:,:,:,f) - Psi * SumEta
  
           up_done = .false.
           lo_done = .false.
           
           do kr=1, n_active
              i = lower(kr)
              j = upper(kr)
  
              if (.not. lo_done(i)) then
                 SumChiUp(:,:,:,i) = SumChiUp(:,:,:,i) * Psi
                 lo_done(i) = .true.
              endif
  
              if (cross(kr) .and. .not. up_done(j)) then
                 SumChiDn(:,:,:,j) = SumChiDn(:,:,:,j) * Psi
                 up_done(j) = .true.
              endif
           end do
  
           do kr=1, n_active
              i = lower(kr)
              j = upper(kr)
              k = trans(kr)
  
              !if (ff==self%fcheck) then
              !   p(:) = self%pcheck
              !   f = self%fcheck
              !   print *, k, i, j, 'cross coupling'
              !   print '(8(ES15.7))', SumChiUp(p(1),p(2),p(3),i), SumChiDn(p(1),p(2),p(3),j), SumU(p(1),p(2),p(3),j), &
              !   Ieff(p(1),p(2),p(3)), SumEta(p(1),p(2),p(3))
              !endif
              
              Gji(:,:,:,k) = Gji(:,:,:,k) * (fac(kr) + Ieff) - SumChiUp(:,:,:,i) * SumU(:,:,:,j)
  
              if (cross(kr)) then
                 Gij(:,:,:,k) = Gij(:,:,:,k) * Ieff + SumChiDn(:,:,:,j) * SumU(:,:,:,i)
              else
                 Gij(:,:,:,k) = Gij(:,:,:,k) * Ieff
              endif
  
              ! add to upward rate coefficient Gamma_ij  
              ! add to downward rate coefficient Gamma_ji 
              gamma_ij(:,:,:,k) = gamma_ij(:,:,:,k) + Gij(:,:,:,k)
              gamma_ji(:,:,:,k) = gamma_ji(:,:,:,k) + Gji(:,:,:,k)
           end do
  
           call trace%end(itimer4)
        end do
     END SUBROUTINE add_all_gammas2
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE solve_stat_equilibrium (self, temp, ne, totn, n, nstar, nh, gamma_ij, gamma_ji, cswitch)
        class(atom_type), target :: self
        class(transition_type), pointer :: trans
        real(4), intent(in) :: temp, ne, totn
        real(4), dimension(:), intent(in) :: nh
        real(8), dimension(:), intent(in) :: gamma_ij, gamma_ji
        real(4), dimension(:), intent(in) :: nstar
        real(4), dimension(:), intent(inout) :: n
        integer, intent(in) :: cswitch
        integer :: kr, pos, i
        real(8) :: tmax_update =1.0d7
        
        if (temp > tmax_update) then
           return !self%e = 1.0d0
        endif
  
        self%w=0.0d0
        self%e=0.0d0
        
        ! Add radiative rates to W
        do kr=1, self%nrad
            trans => self%get_trans(kr)
            self%w(trans%j,trans%i) = gamma_ij(kr)
            self%w(trans%i,trans%j) = gamma_ji(kr)
        enddo
  
        select case(trim(self%crout))
        case('gencol')
            call gencol_new(self%coll,self%id,self%nlevel,self%g,self%ev,self%ion,&
                        self%c,dble(temp),dble(ne),nstar,dble(nh))
        case default
            call stop('read_atom: collisional routine keyword: "' // self%crout // '" not supported.')
        end select
  
        self%w = self%w + self%c * cswitch
  
        ! fill diagonal elements of rate coefficient matrix
        do i=1, self%nlevel
            self%w(i,i) = -1.d0*sum(self%w(:,i))
        end do
  
        ! go from rate coefficients to actual rates   
        do i=1,self%nlevel
            self%w(i,:) = self%w(i,:) * n 
        enddo
  
        ! replace row with largest population with particle conservation
        ! equation.
        pos = maxloc(n,1)
        self%w(pos,:) = n
        self%e(pos)   = totn
        
        ! solve w n . x = e
        ! e contains the solution vector x on exit
        call lineq(self%w, self%e, self%nlevel)
  
        ! set limits on corrections
        where(self%e < 0.001d0) self%e = 0.001d0
  
        ! update the populations
        n = n * self%e
  
        ! set e to the relative change in the populations 
        self%e = self%e - 1.0d0
     END SUBROUTINE solve_stat_equilibrium
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     function get_trans(self, kr) result(trans)
        class(atom_type), intent(inout), target :: self
        integer, intent(in) :: kr
        class(transition_type), pointer :: trans
  
        if (kr > self%nline) then
           trans => self%cont(kr-self%nline)
        else
           trans => self%line(kr)
        endif
     end function
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE compute_lte_populations(self, temp, ne)
        class(atom_type), intent(inout) :: self
        integer :: i
        logical :: molec
        real(4), dimension(:,:,:), intent(in) :: temp, ne
        real(4), dimension(:,:,:), allocatable :: partition
        integer :: nx, ny, nz
        !---------------------------------------------------------------------------------
        nx = size(temp, 1)
        ny = size(temp, 2)
        nz = size(temp, 3)
        
        allocate(partition, mold=temp)
        if (.not. associated(self%nstar)) allocate(self%nstar(nx, ny, nz, self%nlevel))
        if (.not. associated(self%n    )) allocate(self%n    (nx, ny, nz, self%nlevel))
  
        select case (self%id(1:len_trim(self%id)))
        case ('ch','co','cn')
           molec=.true.
           !partition = totn * self%g(1) * exp(-e_k/temp*self%ev(1)) / qcalc(elemid, temp)
           print *, 'molecule partition functions need to be updated. Stopping...'
           stop
        case default
           molec=.false.
           partition = self%partition_function(temp, ne, self%totn)
        end select
  
        do i=1, self%nlevel
           self%nstar(:,:,:,i) = self%lte_pop(i, temp, ne, partition)
  
           if(any(self%nstar(:,:,:,i) <= 0.0)) then
              print *, ' ltepop: level ',i,' has zero lte population. Stopping...'
              stop
           end if
        end do
  
        ! scale total abundances if treated model is a molecule
        if(molec) then
           self%totn = partition * self%nstar(:,:,:,1)
        endif
  
        do i=1, self%nlevel
           self%n(:,:,:,i) = self%nstar(:,:,:,i) 
        end do
     END SUBROUTINE compute_lte_populations
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     elemental function lte_pop(self, i_lvl, temp, ne, partition) result(nstar)
        class(atom_type), intent(in) :: self
        real(4), intent(in) :: temp, ne, partition
        integer, intent(in) :: i_lvl
        real(8) :: evi, dxi, tns, nstar, saha
        logical :: molec
        integer :: i, j, l
        !---------------------------------------------------------------------------------
  
        if (i_lvl == 1) then
           nstar = partition
           return
        endif
  
        dxi=4.98e-4*5040./temp*sqrt(ne*bk*temp)
        !saha = temp*sqrt(temp) / (saha_const*ne)
        saha = log(ne*saha_const)-1.5*log(temp)
  
        j = self%ion(i_lvl)
        evi  = (self%ev(i_lvl) - self%ev(1)) - dxi * 0.5*(j-1)*j
  
        !tns = self%g(i_lvl)/self%g(1) * exp(-e_k/temp*evi)
        tns = log(self%g(i_lvl)) - log(self%g(1)) - e_k/temp*evi
  
        if(j > self%ion(1)) then
           l = j - self%ion(1)
           !tns = tns * saha**l
           tns = tns-float(l)*saha
        endif
  
        nstar = exp(tns) * partition
  
     end function lte_pop
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     elemental function partition_function(self, temp, ne, totn) result(sumn)
        class(atom_type), intent(in) :: self
        real(4), intent(in) :: temp, ne, totn
        real(8) :: sumn, evi, dxi, tns, saha
        integer :: i, j, l
        !--------------------------------------------------------------------------------
        dxi = 4.98e-4*5040./temp*sqrt(ne*bk*temp)
        !saha = temp*sqrt(temp) / (saha_const*ne)
        saha = log(ne*saha_const)-1.5*log(temp)
  
        sumn = 1.
  
        do i=2, self%nlevel
           j = self%ion(i)
           evi  = (self%ev(i) - self%ev(1)) - dxi * 0.5*(j-1)*j
  
           !tns = self%g(i)/self%g(1)*exp(-e_k/temp*evi)
           tns = log(self%g(i)) - log(self%g(1)) - e_k/temp*evi
  
           if(j > self%ion(1)) then
              l = j - self%ion(1)
              !tns = tns * saha**l
              tns = tns - float(l)*saha
           endif
  
           sumn = sumn + exp(tns)
        enddo 
  
        sumn = totn / sumn
     end function
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     pure elemental function compute_doppler_factor(self, temp) result(doppler_factor)
        class(atom_type), intent(in) :: self
        real(4),  intent(in)  :: temp
        real(4) :: doppler_factor
        !---------------------------------------------------------------------------------
        doppler_factor = sqrt( 2.0d0*bk / self%awgt * temp ) / cc
     end function compute_doppler_factor
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE compute_line_damping(self, ne, temp, nh) ! Compute the damping array
        class(atom_type), target, intent(inout) :: self
        real(4),dimension(:,:,:,:),intent(in)  :: nh
        real(4),dimension(:,:,:),  intent(in)  :: ne, temp
        real(4), dimension(:,:,:), pointer, contiguous :: tothi, unsold_temp, stark_pe
        integer :: i, kr, nterm, iterm
        class(line_type), pointer :: line
        !---------------------------------------------------------------------------------
        call os%allocate(tothi, shape(temp))
        call os%allocate(unsold_temp, shape(temp))
        call os%allocate(stark_pe, shape(temp))
  
        tothi = sum(nh(:,:,:,1:5), 4) ! total neutral hydrogen density
        unsold_temp =  8.0*BK / PI * (1.0/MH + 1.0/self%awgt) * temp
        stark_pe = ne * bk * temp**one_sixth
     
        do kr=1,self%nline
           line => self%line(kr)
           call os%allocate(line%adamp, shape(temp))
           line%adamp = 0
           call line%damping(line%adamp, temp, ne, tothi, unsold_temp, stark_pe)
           if (self%id == 'H') then
            call line%dampingH(line%adamp, temp, ne, nh(:,:,:,6))
           endif
        enddo
  
        call os%deallocate(tothi)
        call os%deallocate(unsold_temp)
        call os%deallocate(stark_pe)
     END SUBROUTINE compute_line_damping
     !------------------------------------------------------------------------------------
     
     !====================================================================================
     SUBROUTINE compute_cont_factors(self)
        class(atom_type), target, intent(inout) :: self
        integer :: i, j, kr
        class(cont_type), pointer :: cont
        !---------------------------------------------------------------------------------
        do kr=1,self%ncont
           cont => self%cont(kr)
           call os%allocate(cont%ni_over_nj_lte, shape(self%nstar(:,:,:,1)))
           cont%ni_over_nj_lte = self%nstar(:,:,:,cont%i) / self%nstar(:,:,:,cont%j)
        enddo
     END SUBROUTINE compute_cont_factors
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE zero_sums(sum_arr, indices, clean, n)
        real(4), dimension(:,:,:,:), intent(inout) :: sum_arr
        integer, dimension(:), intent(in) :: indices
        logical, dimension(:), intent(inout) :: clean
        integer, intent(in) :: n
        integer :: i
        !---------------------------------------------------------------------------------
        clean = .false.
        do i=1, n
           if (.not. clean(indices(i))) then
              sum_arr(:,:,:, indices(i)) = 0
              clean(indices(i)) = .true.
           endif
        end do
     END SUBROUTINE zero_sums
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE read (self, atom_file, qnorm, abnd, logf_cut, red_cut, blue_cut, in_lines, teff)
        !---------------------------------------------------------------------
        !
        !  reads an atomic model and calculates einstein a,b coefficients.
        !
        !---------------------------------------------------------------------
        use mclu
        use coll_mod, only: get_ncoll, read_collisions
        class(atom_type), target :: self
        character(len=*), intent(in) :: atom_file
        real(8), intent(in) :: qnorm
        real(8), intent(in), optional :: abnd, logf_cut, red_cut, blue_cut, teff
        integer, intent(in), optional, allocatable :: in_lines(:)
  
        ! Declaration of the local variables
        integer             :: i, j, io, latom, strl, kr, krs, iterm, nterm, ic
        integer             :: ncoll, max_n_temp, doubles
        character(len = 80) :: text, mulfmt
        character(len = 9), parameter :: rname='read_atom'
        character(len = :), allocatable :: myfmt ! temporary format string
        type(line_type), pointer :: lin => null()
        type(cont_type), pointer :: cnt => null()
        type(line_type),dimension(:),allocatable :: scratch_line
        real(8), allocatable, dimension(:) :: nu, alpha
        real(8) :: scr, gv, gq, cnt_ev
  
        open(newunit=latom, file=atom_file, form='formatted')
  
        call nextentry(latom)
        read(latom,fmt='(a)') self%id
        self%id=lowercase(adjustl(self%id),len_trim(self%id))
  
        call nextentry(latom)
        read(latom,*) self%abnd, self%awgt
        self%awgt=self%awgt*uu ! atomic weight [g]
        ! AJG: quickly change abundance without having to edit self
        if (present(abnd)) self%abnd = abnd
  
        call nextentry(latom)
        read(latom,*) self%nlevel,self%nline,self%ncont,self%nrfix
        self%nrad = self%nline + self%ncont
  
        self%qnorm = qnorm
  
        allocate(  self%line  ( self%nline ) ) 
        allocate(  self%cont  ( self%ncont ) )
  
        allocate(  self%ev    ( self%nlevel) ) ; self%ev = 0.0d0
        allocate(  self%g     ( self%nlevel) ) ; self%g = 0
        allocate(  self%label ( self%nlevel) ) ; self%label = ' '
        allocate(  self%ion   ( self%nlevel) ) ; self%ion = 0
  
        allocate(  self%ilin  ( self%nlevel,self%nlevel) ) ; self%ilin=0
        allocate(  self%icon  ( self%nlevel,self%nlevel) ) ; self%icon=0
  
        do i=1, self%nlevel
            call nextentry(latom)
            read(latom,*) self%ev(i),self%g(i), self%label(i), self%ion(i)
            self%ev(i)=self%ev(i)*hc_e  ! [eV]  
        end do
  
        ! sanity check on level energies
        do i = 2, self%nlevel
            if(self%ev(i) < self%ev(i - 1)) then
            myfmt = "('self: levels are: ', a, 2x, a, ' energy not monotonically &
                &increasing.')"
            write(*, myfmt) self%label(i), self%label(i - 1)
            call stop(' ')
            endif
        end do
  
        do i=1,self%nline
            lin => self%line(i)
  
            call nextentry(latom) 
  
            ! io WAS present for historical reasons, in MULTI io determines
            ! whether a line is wide or not. In the RH formulation all lines are wide
  
            read(latom,*) lin%j   , lin%i , lin%f , &
                lin%nnu, lin%qmax, lin%q0, lin%io, &
                lin%ga, gv, gq, lin%profile_type
  
            if (present(in_lines)) then
               if (allocated(in_lines)) then
                 if (findloc(in_lines, i, 1) <= 0) then
                    lin%active = .false.
                    cycle
                 endif
              endif
            endif
  
            if (present(logf_cut)) then
               scr = log10(self%g(lin%i) * lin%f)
               if (scr < logf_cut) then
                  lin%active = .false.
                  cycle
               endif
            endif
  
            ! io is now modified so that multicomponent lines can be read
            ! io = 0: Normal line
            ! io = 1: For now is the same as io = 0
            ! io > 1: Number of components in term
            ! Like in Multi, thereafter follows one extra input line per component:
  
            !   determ  wterm  fterm  gaterm  gwterm  gqterm
            ! determ    energy difference for component i of transition kr. unit cm-1.
            ! wterm     weight of component kterm in weighted average
            ! fterm     oscillator strength
            ! gaterm    damping constants
            ! gwterm
            ! gqterm
            ! adterm    adamp for component, calculated in damp
  
            nterm = lin%io
  
            if (nterm .ge. 2) then
  
                mulfmt = ' '
                write(mulfmt, fmt = '(a)') uppercase(rname, len_trim(rname)) &
                    // ' found ' // string2(nterm) // &
                    'components in transition KR = ' // string3(i)
                call messg(mulfmt)
  
                allocate(lin%determ (nterm) )
                allocate(lin%wterm  (nterm) )
                allocate(lin%fterm  (nterm) )
                allocate(lin%gaterm (nterm) )
                allocate(lin%gwterm (nterm) )
                allocate(lin%gqterm (nterm) )
                !allocate(lin%adterm, xs, xe, ys, ye, zs1, ze1, 1, nterm, rname)
  
                do iterm = 1, nterm
                
                read(latom,*) lin%determ(iterm)          , lin%wterm(iterm) , &
                        lin%fterm(iterm) , lin%gaterm(iterm), lin%gwterm(iterm), &
                        lin%gqterm(iterm)
  
                end do
  
            end if
  
            lin%lower_ev = self%ev(lin%i)
            lin%upper_ev = self%ev(lin%j)
            lin%hydrogen = trim(self%id) == 'h'
  
            lin%lambda0 = hc_e / ( lin%upper_ev - lin%lower_ev ) ! [cm]
            lin%nu0   = cc / lin%lambda0    ! [Hz]
  
            lin%Grat  = self%g(lin%i) / self%g(lin%j)
            lin%Aji   = f_to_Aji * lin%f * lin%Grat / lin%lambda0**2
            lin%Bji   = c2_2h / lin%nu0**3 * lin%Aji 
            lin%Bij   = lin%Bji / lin%Grat

            if (present(teff)) then
               if (teff > 0) then
                  lin%qmax = lin%estimate_qmax(teff)
                  lin%q0 = lin%qmax
               endif
            endif
  
            if (present(red_cut)) then
              scr = lin%lambda0 / (1.0-lin%qmax*qnorm*KM_TO_CM/cc) * 1e8
              if (scr > red_cut) then
                 lin%active = .false.
                 cycle
              endif
            endif
            
            if (present(blue_cut)) then
              scr = lin%lambda0 / (1.0+lin%qmax*qnorm*KM_TO_CM/cc) * 1e8
              if (scr < blue_cut) then
                 lin%active = .false.
                 cycle
              endif
            endif
  
            ! make sure we always have odd number of frequency points for a line
            if (mod(lin%nnu,2).eq.0) lin%nnu=lin%nnu+1
            allocate( lin%nu(  lin%nnu ) )
            lin%nq = lin%nnu
            
            ! find continuum level
            do ic=lin%j+1, self%nlevel
               if (self%ion(ic) .eq. self%ion(lin%j)+1) then
                 cnt_ev = self%ev(ic)
                 exit
               endif
            end do

            call lin%prepare_Stark(gq)
            call lin%prepare_VanderWaals(gv, self%ion(lin%i), cnt_ev)

            lin%profile_type = &
                adjustl(lowercase(lin%profile_type,len_trim(lin%profile_type)))
  
            lin%ired = 1
            lin%iblue = lin%nnu
  
            lin%ntrans=i
            self%ilin(lin%j,lin%i) = i
            self%ilin(lin%i,lin%j) = i
  
            if ( (trim(lin%profile_type) .ne. 'voigt'  ) .and. &
                (trim(lin%profile_type) .ne. 'prd'    ) .and. &
                (trim(lin%profile_type) .ne. 'doppler') .and. &
                (trim(lin%profile_type) .ne. 'hlinop' )         ) then 
  
                if ( self%id .eq. 'h' ) then
  
                 write (*,fmt = '(a,1x,a,a)') 'Profile not recognised:', &
                          uppercase(lin%profile_type, len_trim(lin%profile_type)), &
                          ' Setting it to HLINOP'
  
                 lin%profile_type = adjustl(lowercase('hlinop', len_trim('hlinop')))
                
                else
  
                 write (*, fmt = '(a,1x,a,a)') 'Profile not recognised:', &
                          uppercase(lin%profile_type, len_trim(lin%profile_type)), &
                          ' Setting it to VOIGT'
  
                 lin%profile_type = adjustl(lowercase('voigt', len_trim('voigt')))
                
                end if
            end if
        end do
  
        j = 0
        do i=1,self%nline
           if (self%line(i)%active) j = j + 1
        end do
        
        if (j /= self%nline) then
           allocate(scratch_line(j))
        
           j = 0
           do i=1,self%nline
              if (self%line(i)%active) then
                 j = j + 1
                 scratch_line(j) = self%line(i)
                 scratch_line(j)%ntrans=j
              endif
           end do
        
           deallocate(self%line)
           allocate(self%line(j))
           self%line = scratch_line
           deallocate(scratch_line)
           self%nline = j
           self%nrad = self%nline + self%ncont
        endif
  
        !  bound-free transitions in detail
        !  if numax.lt.0.0 frequency points in angstrom (starting with
        !  threshold and decreasing) and crossections in cm2 are read.
        do i=1,self%ncont
            cnt => self%cont(i)
            call nextentry(latom) 
  
            read(latom,*) cnt%j,cnt%i,cnt%alpha0,cnt%nnu,cnt%numax
  
            cnt%lower_ev = self%ev(cnt%i)
            cnt%upper_ev = self%ev(cnt%j)
        
            cnt%nu0   = ( cnt%upper_ev - cnt%lower_ev ) / h_e ! [Hz]
            cnt%numax = cc/(cnt%numax*AA_TO_CM) ! [from AA to Hz]
            cnt%ntrans=self%nline+i ! transition identifier for continua
  
            cnt%ired=1
            cnt%iblue=cnt%nnu
            self%icon( cnt%j, cnt%i ) = i
            self%icon( cnt%i, cnt%j ) = i
  
            if ( cnt%numax .lt. 0.0 ) then
  
                allocate( nu(cnt%nnu))
                allocate( alpha(cnt%nnu))
  
                ! cross sections given in detail
                cnt%bf_type='detail'
                do j=1,cnt%nnu
                 call nextentry(latom) 
                 read(latom,*)  nu(j), alpha(j)
                end do
                ! sanity check
                doubles = 0
                do j=2,cnt%nnu
                 if ( nu(j) .gt. nu(j-1) ) then
                       write(text,*) 'read_atom: continuum transition nr', i, ' wavelengths not decreasing'
                       call stop(text)
                 elseif ( nu(j) == nu(j-1) ) then
                    doubles = doubles + 1
                 endif
                end do
  
                ! making sure here that there are no double wavelenghts
                cnt%nnu = cnt%nnu - doubles
                allocate(cnt%nu(cnt%nnu), cnt%alpha(cnt%nnu))
                cnt%nu(1) = nu(1)
                cnt%alpha(1) = alpha(1)
  
                doubles = 1
                do j=2, cnt%iblue
                 if ( nu(j) < nu(j-1) ) then
                    doubles = doubles + 1
                    cnt%nu(doubles) = nu(j)
                    cnt%alpha(doubles) = alpha(j)
                 endif
                end do
  
                cnt%iblue = cnt%nnu
  
                deallocate( nu)
                deallocate( alpha)
  
                ! and convert to frequency [Hz]
                cnt%nu = cc / ( cnt%nu * AA_TO_CM )
                cnt%numax = cnt%nu(cnt%nnu)
  
            else
                allocate( cnt%nu(1:cnt%nnu)    )
                allocate( cnt%alpha(1:cnt%nnu) )
  
                ! hydrogenic transitions
                cnt%bf_type='hydrogenic'
  
            endif
  
        end do
  
        ! read choice of collisional SUBROUTINE
        ! GENCOL and HCOL are currently the only collisions subroutines.
        call nextentry(latom)
        read(latom,fmt='(a)') self%crout
  
        ! sanity check
        self%crout = adjustl(self%crout)
        strl=len_trim(self%crout)
        self%crout = lowercase(self%crout,strl)
  
        ! read collisional data
        call get_ncoll(latom, self%crout, ncoll, max_n_temp)
        allocate( self%coll(ncoll) )
        call read_collisions(latom, self%coll, max_n_temp, self%nlevel)
  
        close(latom)

        ! frequency grid in lines
        do i=1,self%nline
          call self%line(i)%freql(self%qnorm, symmetric=.true.)
        end do

        !  frequency grid in continua
        do i=1,self%ncont
            call self%cont(i)%freqc()
        end do
     END SUBROUTINE read
     !------------------------------------------------------------------------------------

    !===============================================================================
    ! initialises the spectrum structure
    !===============================================================================
    SUBROUTINE create_spectrum (self, sp)
      class(atom_type), target :: self
      class(spectrum_type), intent(out) :: sp
      type(line_type), pointer :: lin
      type(cont_type), pointer :: cnt
      integer :: i,nnu,js,je
      real(8),dimension(:),allocatable :: nu
      integer, save:: itimer=0
      !-----------------------------------------------------------------------------
      call trace%begin('atom_type%create_spectrum', itimer=itimer)
      ! compute total number of frequencies
      nnu=0
      do i=1,self%nline
         nnu=nnu+self%line(i)%nnu
      end do
      do i=1,self%ncont
         nnu=nnu+self%cont(i)%nnu
      end do
  
      ! save all frequencies
      allocate(nu(nnu))
      je=0
      do i=1,self%nline
         js=je+1
         je=js+self%line(i)%nnu-1
         nu(js:je)=self%line(i)%nu
      end do
      do i=1,self%ncont
         js=je+1
         je=js+self%cont(i)%nnu-1
         nu(js:je)=self%cont(i)%nu
      end do
  
      ! sort frequencies into ascending order 
      call refsor(nu)
   
     ! discard duplicates, nnu contains number of unique elements on exit
      call unista(nu,nnu)
   
     ! now set overall frequency grid
      sp%nnu=nnu
      allocate(sp%nu(nnu)) 
      sp%nu=nu(1:nnu)
      call trace%end(itimer)
    END SUBROUTINE create_spectrum

    !===============================================================================
    SUBROUTINE set_active_transitions (self, sp)
      class(atom_type), target :: self
      class(atomic_spectrum_type), intent(inout) :: sp
      type(line_type), pointer :: lin
      type(cont_type), pointer :: cnt
      real(8), dimension(:), allocatable :: dalpha_dnu,alpha_fromfile,nu_fromfile, alpha
      integer :: i, j, nnu, ired, iblue
      integer, save:: itimer=0
      !-----------------------------------------------------------------------------
      call trace%begin('atom_type%set_active_transitions', itimer=itimer)
     
      ! allocate list of active continua and lines
      allocate(sp%ac(self%ncont,sp%nnu)) ; sp%ac=0
      allocate(sp%nac(sp%nnu)) ; sp%nac=0
  
      allocate(sp%al(self%nline,sp%nnu)) ; sp%al=0
      allocate(sp%nal(sp%nnu)) ; sp%nal=0
      
      do i=1,self%ncont
        cnt => self%cont(i)

        ! find minimum and maximum index of this continuum in the sp%nu grid
        ired  = searchsorted(sp%nu, cnt%nu(1), outbounds_ok=.true.)
        iblue = searchsorted(sp%nu, cnt%nu(cnt%nnu), outbounds_ok=.true.)
        if (iblue == 0 .or. ired == sp%nnu) then
         cnt%ired  = 0
         cnt%iblue = 0
         cycle
        endif

        if (ired == 0) ired = 1
        ! total number of frequencies where this continuum is active
        nnu=iblue-ired+1

        ! min and max frequency where this continuum is active
        cnt%ired  = ired
        cnt%iblue = iblue
      
        ! add to number of active continua
        do j=ired,iblue
          sp%nac(j) = sp%nac(j) + 1 ! number of active continua at frequency index j
          sp%ac( sp%nac(j), j ) = i ! the index of this active continuum
        end do
        
        ! recompute the bf cross sections
        select case(cnt%bf_type)

        case('hydrogenic')
          ! use nu**-3 dependency
          cnt%nnu = nnu

          deallocate( cnt%nu      )
          allocate(   cnt%nu(ired:iblue) )
          cnt%nu=sp%nu(ired:iblue)

          deallocate( cnt%alpha      )
          allocate(   cnt%alpha(ired:iblue) )
          call cnt%freqc()

        case('detail')
          ! interpolate from given values
          if (allocated(dalpha_dnu)) deallocate(dalpha_dnu)
          allocate( dalpha_dnu(cnt%nnu) )
          call hermite_derivative(dalpha_dnu, cnt%nu, dble(cnt%alpha))
          
          if (allocated(alpha_fromfile)) deallocate(alpha_fromfile)
          allocate( alpha_fromfile(cnt%nnu) )
          alpha_fromfile = cnt%alpha

          if (allocated(nu_fromfile)) deallocate(nu_fromfile)
          allocate( nu_fromfile(cnt%nnu) )
          nu_fromfile = cnt%nu
          
          ! set new frequency grid        
          cnt%nnu   = nnu
          deallocate( cnt%nu      )
          allocate(   cnt%nu(ired:iblue) )
          cnt%nu=sp%nu(ired:iblue)
        
          ! allocate space for new cross sections
          deallocate( cnt%alpha             )
          allocate(   cnt%alpha(ired:iblue) )

          ! fill new cross sections 
          allocate(alpha(cnt%nnu))
          alpha = dble(cnt%alpha)
          call hermite_interpolation(alpha, cnt%nu, nu_fromfile, alpha_fromfile, dalpha_dnu)
          cnt%alpha = real(alpha)
          deallocate(alpha)

        case default
          print *, 'sp_from_atom : cont(i)%bf_type'// cnt%bf_type // 'not recognized'
          stop

        end select

        call sp%set_weights(cnt%wnu, cnt%nu, ired, iblue)
      end do

      do i=1,self%nline
          lin => self%line(i)

          ! find minimum and maximum index of this line in the sp%nu grid
          ired  = searchsorted(sp%nu, lin%nu(1), outbounds_ok=.true.)
          iblue = searchsorted(sp%nu, lin%nu(lin%nnu), outbounds_ok=.true.)
          if (iblue == 0 .or. ired == sp%nnu) then
           lin%ired  = 0
           lin%iblue = 0
           cycle
          endif

          if (ired == 0) ired = 1
          ! total number of frequencies where this line is active
          nnu=iblue-ired+1

          ! min and max frequency where this line is active
          lin%ired  = ired
          lin%iblue = iblue
        
          ! add to number of active lines
          do j=ired,iblue
            sp%nal(j) = sp%nal(j) + 1
            sp%al( sp%nal(j), j ) = i
          end do

          ! set new frequency grid for the line      
          lin%nnu   = nnu
          deallocate( lin%nu      )
          allocate(   lin%nu(ired:iblue) )
          lin%nu=sp%nu(ired:iblue) ! [Hz]
          allocate(   lin%q(ired:iblue) )
          lin%q = (lin%nu/lin%nu0-1.0d0) * CC / (self%qnorm * KM_TO_CM) ! [Doppler units]

          ! fill frequency integration weights for the line
          call sp%set_weights(lin%wnu, lin%nu, ired, iblue)
          call sp%set_weights(lin%wq , lin%q , ired, iblue)
      end do

      sp%maxac=maxval(sp%nac)
      sp%maxal=maxval(sp%nal)

      if (sp%verbose>0) then
          write(*,'(a,I0,a)') 'sp_from_atom: frequency grid has ', sp%nnu, ' points.'
      endif
      call trace%end(itimer)
    END SUBROUTINE set_active_transitions
  
  END MODULE atom_mod
  
  ! *******************************************************************
  