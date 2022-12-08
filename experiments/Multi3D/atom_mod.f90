! *******************************************************************
MODULE atom_mod
    USE io_mod
    USE m3io_mod
    USE prec_mod, only: inf8
    USE mconst
    USE coll_mod
    USE trans_mod
    USE mmath, only: refsor, unista
    USE op_math_mod, only: lineq, lineq_no_pivot, get_interval
    USE ion_mod, only: ion_type
    USE element_mod, only: element_type
    USE molecule_mod, only: molec_type
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
        integer, dimension(:,:), contiguous, pointer :: ac,al
        ! how many continua (nac) and lines (nal) are active at each frequency
        integer, dimension(:), contiguous, pointer :: nac,nal
        ! has this frequency an active prd line
        !logical,dimension(:),allocatable :: hasPRDline
    end type atomic_spectrum_type
    !---------------------------------------------------------------------
  
    type :: atom_type
        integer :: nrad   ! total number of radiative transitions
        integer :: nrfix  ! number of fixed transitions
        integer :: ncont  ! number of bound-free transitions
        integer :: nline  ! number of bound-bound transitions
        integer :: nlevel ! number of levels
        integer :: nprd=0 ! number of PRD transitions
        
        character(len=3) :: id    ! atom id string 
        character(len=20) :: crout ! which collisional SUBROUTINE to use 
   
        character(len=20), dimension(:), allocatable :: label ! level label
   
        integer, dimension(:), allocatable :: ion ! ionization stage of each level
   
        real(8) :: abnd ! abundance ( in the [hydrogen]=12 scale)
        real(8) :: awgt ! atomic weight [g]
   
        real(8), dimension(:), allocatable :: ev ! energy of level above ground state [eV]
        real(8), dimension(:), allocatable :: g ! multiplicity of level
  
        type(line_type), dimension(:), allocatable :: line
        type(cont_type), dimension(:), allocatable :: cont
        type(coll_type), dimension(:), pointer, contiguous :: coll
        type(atomic_spectrum_type), pointer :: sp
  
        ! rather constant doppler width factor, i.e. dnud = line(i)%nu0 * doppler_factor [Hz]
        real(4), dimension(:,:,:), pointer, contiguous :: doppler_factor
  
        ! LTE, NLTE and total population of active atom
        real(4), dimension(:,:,:,:), pointer, contiguous :: nstar => null()
        real(4), dimension(:,:,:,:), pointer, contiguous :: n     => null()
        real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij => null()
        real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ji => null()
        real(4), dimension(:,:,:), pointer, contiguous :: totn
        real(8), dimension(:,:), pointer, contiguous :: tgrid, cgrid, dcdt
        real(4), dimension(:), pointer, contiguous :: cdata
        logical, dimension(:), pointer, contiguous :: ecol, hcol
  
        real(8), dimension(:),   allocatable  :: e
        real(8), dimension(:,:), allocatable  :: C, W
        
        real(8) :: qnorm
        integer :: fcheck, pcheck(3), linecheck
        logical :: stiff 
        
     contains
        !procedure :: add_line_opacities
        procedure :: add_continuum_opacities
        procedure :: full_preconditioning
        procedure :: compute_doppler_factor
        procedure :: partition_function
        procedure :: lte_pop
        procedure :: compute_lte_populations
        procedure :: compute_cont_factors
        procedure :: solve_stat_equilibrium
        procedure :: create_spectrum
        procedure :: set_active_transitions
        procedure :: read
        procedure :: copy_to
     
    end type atom_type

    logical :: lse_no_pivot=.True.
  
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
        !------------------------------------------------------------------------------------
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
  
!     !====================================================================================
!     SUBROUTINE add_line_opacities(self, fs_patch, fe_patch, nu, temp, vproj, chi, eta)
!        class(atom_type), intent(in), target :: self
!        integer, intent(in) :: fs_patch, fe_patch
!        real(8), dimension(:), intent(in) :: nu
!        real(4), dimension(:,:,:), intent(in) :: temp, vproj
!        real(4), dimension(:,:,:,:), intent(inout) :: chi, eta
!        real(4), dimension(:,:,:,:), allocatable, target :: scratch_phi
!        real(4), dimension(:,:,:,:), pointer :: phi
!        real(4), dimension(:,:,:), allocatable :: vel_dop
!        real(8), dimension(:), allocatable :: hnu_4pi
!        class(line_type), pointer :: line
!        integer :: i, fcs, fce, fs, fe, nx, ny, nz, ffs, ffe, f, p(3)
!        logical :: skip
!  
!        ! all_chi and all_eta only range from fs:fe
!        ! nu however contains all frequencies from atom
!  
!        fcs = max(minval(self%line%ired),  fs_patch)
!        fce = min(maxval(self%line%iblue), fe_patch)
!  
!        if (fcs >= fe_patch) return
!        if (fce <= fs_patch) return
!  
!        nx = size(temp,1)
!        ny = size(temp,2)
!        nz = size(temp,3)
!  
!        allocate(hnu_4pi(fcs:fce))
!        allocate(scratch_phi(nx,ny,nz,fcs:fce))
!  
!        hnu_4pi = h_4pi * nu(fcs:fce)
!        vel_dop = (1 + vproj) / self%doppler_factor
!  
!        do i=1, self%nline
!           line => self%line(i)
!           call line%check_range(fcs,fce,fs,fe,skip)
!           if (skip) cycle
!        
!           ffs = fs-fs_patch+1
!           ffe = fe-fs_patch+1
!  
!           phi => scratch_phi(:,:,:,fs:fe)
!  
!           call line%opacities(nu(fs:fe), hnu_4pi(fs:fe), self%n, vel_dop, self%doppler_factor, phi, eta(:,:,:,ffs:ffe), chi(:,:,:,ffs:ffe))
!  
!           !if (self%fcheck >= fs .and. self%fcheck <= fe) then
!           !   p(:) = self%pcheck
!           !   f = self%fcheck
!           !   print *, i, line%i, line%j, 'chi line'
!           !   print '(8(ES15.7))', chi(p(1),p(2),p(3),f-fs_patch+1), &
!           !             scratch_phi(p(1),p(2),p(3), f) &
!           !          * (line%Bij * self%n(p(1),p(2),p(3),line%i) &
!           !          -  line%Bji * self%n(p(1),p(2),p(3),line%j)), &
!           !          scratch_phi(p(1),p(2),p(3), f), &
!           !          self%n(p(1),p(2),p(3),line%i), &
!           !          self%n(p(1),p(2),p(3),line%j), nu(f)
!           !endif
!        end do
!     END SUBROUTINE add_line_opacities
!     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE full_preconditioning(self, fs, temp, Ieff, Psi, gamma_ij, gamma_ji, idx_map, lines, lin_idx_map)
        class(atom_type), target, intent(in) :: self
        integer, intent(in) :: fs
        real(4), dimension(:,:,:), intent(in), contiguous :: temp
        real(4), dimension(:,:,:,:), intent(in), contiguous :: Psi
        real(4), dimension(:,:,:,:), intent(inout), contiguous :: gamma_ij, gamma_ji, Ieff
        class(line_type), dimension(:), intent(in) :: lines
        integer, dimension(:), intent(in) :: idx_map, lin_idx_map
        real(4), dimension(:,:,:), allocatable :: blzfac, hk_over_temp
        real(4), dimension(:,:,:), pointer, contiguous :: ni, nj
        real(4), dimension(:,:,:,:), allocatable :: SumChiDn, SumChiUp, phi, SumU
        real(4), dimension(:,:,:,:), allocatable :: vij, vji
        real(4), dimension(:), allocatable :: wnu
        integer, dimension(:), allocatable :: lower, upper, trans, ii, jj
        logical, dimension(:), allocatable :: cross
        real(4) :: chi, twohnu3_c2, alpha, hnu_inv, Bij, Bji, ufac, alpha2
        integer :: i, j, k, kr, kr2, kc, kl, fe, n_active, ff, f, p(3), x,y,z, nx, ny, nz, nupper, nlower, i2
        logical :: cross_coupling
        integer, save :: itimer1, itimer2, itimer3, itimer4, itimer5, itimer6
        !------------------------------------------------------------------------------------
        fe = fs + size(Psi, 4) - 1 
        nx = size(Psi,1)
        ny = size(Psi,2)
        nz = size(Psi,3)
  
        n_active = maxval(self%sp%nac(fs:fe) + self%sp%nal(fs:fe))

        allocate(wnu(n_active))
        allocate(lower(n_active), upper(n_active), trans(n_active), cross(n_active))
        allocate(ii(n_active), jj(n_active))
        allocate(Vij(nx,ny,nz,n_active))
        allocate(Vji(nx,ny,nz,n_active))
        
        nlower = 0
        nupper = 0
        do ff=fs,fe
          n_active = self%sp%nac(ff) + self%sp%nal(ff)
          do kr=1, n_active
            if (kr <= self%sp%nac(ff)) then
              kc = self%sp%ac(kr,ff)
              lower(kr) = self%cont(kc)%i
              upper(kr) = self%cont(kc)%j
            else
              kl = self%sp%al(kr - self%sp%nac(ff),ff)
              lower(kr) = self%line(kl)%i
              upper(kr) = self%line(kl)%j
            endif
          end do
          
          call unique_int(lower(1:n_active), lower(1:n_active), i)
          call unique_int(upper(1:n_active), upper(1:n_active), j)
          if (i > nlower) nlower = i
          if (j > nupper) nupper = j
        end do
  
        allocate(phi, mold=Psi(:,:,:,1:1))
        allocate(blzfac, hk_over_temp, mold=temp)
        allocate(SumChiUp(nx,ny,nz,nlower))
        allocate(SumChiDn(nx,ny,nz,nupper))
        allocate(SumU  (nx,ny,nz,nupper))
  
        hk_over_temp = h_k / temp

        do ff=fs,fe
           !call trace%begin('administration', itimer=itimer1)
           f = ff - fs + 1
           n_active = self%sp%nac(ff) + self%sp%nal(ff)
  
           hnu_inv = 1 / (hh*self%sp%nu(ff))
           twohnu3_c2 = twoh_c2 * self%sp%nu(ff)**3
  
           if (self%sp%nac(ff) > 0) then
              blzfac = exp(- hk_over_temp * real(self%sp%nu(ff)))
           endif
           
           do kr=1, n_active
              if (kr <= self%sp%nac(ff)) then
                 kc = self%sp%ac(kr,ff)
                 lower(kr) = self%cont(kc)%i
                 upper(kr) = self%cont(kc)%j
                 trans(kr) = idx_map(self%cont(kc)%ntrans)
                 wnu  (kr) = self%cont(kc)%wnu(ff) * hnu_inv ! wm needs to be multiplied later
              else
                 kl = self%sp%al(kr - self%sp%nac(ff),ff)
                 lower(kr) = self%line(kl)%i
                 upper(kr) = self%line(kl)%j
                 trans(kr) = idx_map(self%line(kl)%ntrans)
                 wnu  (kr) = self%line(kl)%wnu(ff) * hnu_inv ! wm needs to be multiplied later
              endif
           end do

           cross(1:n_active) = .false.
           cross_coupling = .false.
  
           do kr=1, n_active
              do kr2=1, n_active
                 if (lower(kr) == upper(kr2)) then
                    cross_coupling = .true.
                    cross(kr) = .true.
                 endif
              end do
           end do

           call unique_int(lower(1:n_active), lower(1:n_active), nlower)
           call unique_int(upper(1:n_active), upper(1:n_active), nupper)

           do kr=1, n_active
            if (kr <= self%sp%nac(ff)) then
              kc = self%sp%ac(kr,ff)
              ii(kr) = findloc(lower(1:nlower), self%cont(kc)%i, 1)
              jj(kr) = findloc(upper(1:nupper), self%cont(kc)%j, 1)
            else
              kl = self%sp%al(kr - self%sp%nac(ff),ff)
              ii(kr) = findloc(lower(1:nlower), self%line(kl)%i, 1)
              jj(kr) = findloc(upper(1:nupper), self%line(kl)%j, 1)
            endif
           end do

           do i=1,nlower
            SumChiUp(:,:,:,i) = 0
           end do

           do j=1,nupper
            SumU(:,:,:,j) = 0
            if (cross_coupling) SumChiDn(:,:,:,j) = 0
           end do
  
           !call trace%end(itimer1)
           !call trace%begin('vij', itimer=itimer2)
  
           do kr=1, n_active
              i = ii(kr)
              j = jj(kr)
              ni => self%n(:,:,:,lower(i))
              nj => self%n(:,:,:,upper(j))

              ufac = twohnu3_c2 / wnu(kr)
  
              if (kr <= self%sp%nac(ff)) then
                 kc = self%sp%ac(kr,ff)
                 ! alpha0 is part of ni_over_nj_lte for overflow reasons
                 alpha = self%cont(kc)%alpha(ff) * wnu(kr)
                 alpha2 = alpha * self%cont(kc)%alpha0
                 do z=1,nz
                 do y=1,ny
                 do x=1,nx
                    vij(x,y,z,kr) = alpha2
                    vji(x,y,z,kr) = alpha * self%cont(kc)%ni_over_nj_lte(x,y,z) * blzfac(x,y,z)
                    chi = vij(x,y,z,kr) * ni(x,y,z) - vji(x,y,z,kr) * nj(x,y,z)
  
                    SumU(x,y,z,j) = SumU(x,y,z,j) + vji(x,y,z,kr) * ufac
                    SumChiUp(x,y,z,i) = SumChiUp(x,y,z,i) + chi
                    if (cross(kr)) SumChiDn(x,y,z,j) = SumChiDn(x,y,z,j) + chi
                 end do
                 end do
                 end do
              else
                 kl = self%sp%al(kr - self%sp%nac(ff),ff)
                 kl = lin_idx_map(kl)
                 call lines(kl)%interpolate_phi(self%sp%nu(ff:ff), phi)
                 Bij = lines(kl)%Bij * wnu(kr)
                 Bji = lines(kl)%Bji * wnu(kr)
                 do z=1,nz
                 do y=1,ny
                 do x=1,nx
                    vij(x,y,z,kr) = phi(x,y,z,1) * Bij
                    vji(x,y,z,kr) = phi(x,y,z,1) * Bji
                    chi = vij(x,y,z,kr) * ni(x,y,z) - vji(x,y,z,kr) * nj(x,y,z)
  
                    SumU(x,y,z,j) = SumU(x,y,z,j) + vji(x,y,z,kr) * ufac
                    SumChiUp(x,y,z,i) = SumChiUp(x,y,z,i) + chi
                    if (cross(kr)) SumChiDn(x,y,z,j) = SumChiDn(x,y,z,j) + chi
                 end do
                 end do
                 end do
  
                 !if (ff==self%fcheck .and. all(self%pcheck>0)) then
                 !   p(:) = self%pcheck
                 !   print *, self%sp%al(kr - self%sp%nac(ff),ff), lower(i), upper(j), 'rate', p
                 !   print '(8(ES15.7))', wnu(kr), vij(p(1),p(2),p(3),kr), vji(p(1),p(2),p(3),kr)
                 !endif
              endif
           end do
  
           !call trace%end(itimer2)
           !call trace%begin('Ieff', itimer=itimer3)

           do j=1,nupper
            SumU(:,:,:,j) = SumU(:,:,:,j) * Psi(:,:,:,f)
           end do
  
           !Ieff = Inu(:,:,:,f) - sum(SumU(:,:,:,1:nupper), 4)
           !Ieff(:,:,:,f) = Ieff(:,:,:,f) - sum(SumU(:,:,:,1:nupper), 4)

           do kr=1, n_active
             i = ii(kr)
             j = jj(kr)
             if (cross(kr)) i2 = findloc(upper(1:nupper), lower(i), 1)
             k = trans(kr)

             do z=1,nz
             do y=1,ny
             do x=1,nx
              vji(x,y,z,kr) = vji(x,y,z,kr) * (twohnu3_c2 + Ieff(x,y,z,f)) - SumChiUp(x,y,z,i) * SumU(x,y,z,j)
              if (cross(kr)) then
               vij(x,y,z,kr) = vij(x,y,z,kr) * Ieff(x,y,z,f) + SumChiDn(x,y,z,j) * SumU(x,y,z,i2)
              else
               vij(x,y,z,kr) = vij(x,y,z,kr) * Ieff(x,y,z,f)
              endif
              gamma_ji(x,y,z,k) = gamma_ji(x,y,z,k) + vji(x,y,z,kr)
              gamma_ij(x,y,z,k) = gamma_ij(x,y,z,k) + vij(x,y,z,kr)
             end do
             end do
             end do

             !if (ff==self%fcheck .and. all(self%pcheck>0)) then
             !   p(:) = self%pcheck
             !   print *, lower(i), upper(j), 'Ieff', cross(kr), p
             !   print '(8(ES15.7))', Ieff(p(1),p(2),p(3)), vij(p(1),p(2),p(3),kr), vji(p(1),p(2),p(3),kr)
             !endif
           end do

           !call trace%end(itimer3)
        end do

        deallocate(Vij)
        deallocate(Vji)
     END SUBROUTINE full_preconditioning
     !------------------------------------------------------------------------------------
    
   !====================================================================================
   SUBROUTINE solve_stat_equilibrium (self, temp, ne, totn, n, nstar, nh, gamma_ij, gamma_ji, cswitch, check_NaN)
      class(atom_type), target :: self
      class(transition_type), pointer :: trans
      real(4), intent(in) :: temp, ne, totn
      real(4), dimension(:), intent(in) :: nh
      real(8), dimension(:), intent(in) :: gamma_ij, gamma_ji
      real(4), dimension(:), intent(in) :: nstar
      real(4), dimension(:), intent(inout) :: n
      integer, intent(in) :: cswitch
      logical, intent(in) :: check_NaN
      integer :: kr, pos, i, j
      real(8) :: tmax_update =1.0d7
      integer, save :: itimer=0, itimer2=0
      !------------------------------------------------------------------------------------
      if (temp > tmax_update) then
         return !self%e = 1.0d0
      endif

      self%w=0.0d0
      self%e=0.0d0
      
      ! Add radiative rates to W
      do kr=1, self%nrad
         if (kr <= self%nline) then
            i = self%line(kr)%i
            j = self%line(kr)%j
         else
            i = self%cont(kr-self%nline)%i
            j = self%cont(kr-self%nline)%j
         endif
         self%w(j,i) = gamma_ij(kr)
         self%w(i,j) = gamma_ji(kr)
      enddo

      call trace%begin('gencol', itimer=itimer)
      self%c = 0.d0

      !call classic_col(self%g,self%ev,self%ion,dble(temp),dble(ne),dble(nh), nstar, &
      !                 self%cdata,self%c,self%ecol,self%hcol,self%awgt)

      !print '(4(4(ES17.8), /))', self%c
      !self%c = 0

      call gencol(self%coll,self%id,self%nlevel,self%g,self%ev,self%ion,&
                  self%c,dble(temp),dble(ne),nstar,dble(nh),self%tgrid,self%cgrid,self%dcdt)

      !print '(4(4(ES17.8), /))', self%c
      !stop

      call trace%end(itimer=itimer)

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

      if (check_NaN) then
         do i=1,self%nlevel
            if (self%w(i,i) == 0.d0) then
               print *, 'zero in rate diagonal at i = ', i
               print *, 'Temp = ', temp
               stop
            endif
            do j=1,self%nlevel
               if (self%w(i,j) /= self%w(i,j)) then
                  print *, 'NaN in rate Matrix at i,j = ', i, j
                  print *, 'Temp = ', temp
                  stop
               elseif (self%w(i,j) > inf8) then
                  print *, 'Infinity in rate Matrix at i,j = ', i, j
                  print *, 'Temp = ', temp
                  stop
               endif
            end do
         end do
      endif

      call trace%begin('lineq_no_pivot', itimer=itimer2)

      ! solve w n . x = e
      ! e contains the solution vector x on exit
      if (lse_no_pivot) then
        call lineq_no_pivot(self%w, self%e, self%nlevel)
      else
        call lineq(self%w, self%e, self%nlevel)
      endif

      call trace%end(itimer2)

      if (check_NaN) then
         do i=1,self%nlevel
            if (self%e(i) /= self%e(i)) then
               print *, 'NaN in error vector at i = ', i
               print *, 'Temp = ', temp
               stop
            endif
         end do
      endif

      ! set limits on corrections
      where(self%e < 0.001d0) self%e = 0.001d0

      ! update the populations
      n = n * self%e

      ! set e to the relative change in the populations 
      self%e = self%e - 1.0d0
   END SUBROUTINE solve_stat_equilibrium
   !------------------------------------------------------------------------------------
  
   !====================================================================================
   SUBROUTINE compute_lte_populations(self, temp, ne, dxi, element, mol)
      class(atom_type), intent(inout) :: self
      class(element_type), pointer :: element
      class(ion_type), pointer :: ion
      class(molec_type), pointer :: mol
      integer :: i, j, x,y,z
      real(4), dimension(:,:,:), intent(in) :: temp, ne, dxi
      real(4), dimension(:,:,:), allocatable :: partition, saha
      integer :: nx, ny, nz
      logical :: zero_pop
      !---------------------------------------------------------------------------------
      nx = size(temp, 1)
      ny = size(temp, 2)
      nz = size(temp, 3)
      
      allocate(partition, mold=temp)
      allocate(saha, mold=temp)
      if (.not. associated(self%nstar)) allocate(self%nstar(nx, ny, nz, self%nlevel))
      if (.not. associated(self%n    )) allocate(self%n    (nx, ny, nz, self%nlevel))

      !dxi = 4.98e-4*5040./temp*sqrt(ne*bk*temp)

      if (associated(element)) then
         partition = self%totn
         saha = 0.0
         do i=1, self%nlevel
            j = self%ion(i)
            ion => element%ions(j)
            self%nstar(:,:,:,i) = self%lte_pop(i, temp, dxi, saha, partition)
            self%nstar(:,:,:,i) = self%nstar(:,:,:,i) * ion%anjon * ion%part
            self%n(:,:,:,i)     = self%nstar(:,:,:,i) 
         end do
      else
         !saha = temp*sqrt(temp) / (saha_const*ne)
         saha = 1.5*log(temp) - log(saha_const*ne) ! actually log(saha)

         if (associated(mol)) then
            partition = self%totn / mol%irwin_partf(dble(temp))
         else
            !select case (self%id(1:len_trim(self%id)))
            !case ('ch','co','cn')
            !   molec=.true.
            !   !partition = totn * self%g(1) * exp(-e_k/temp*self%ev(1)) / qcalc(elemid, temp)
            !   print *, 'molecule partition functions need to be updated. Stopping...'
            !   stop
            !case default
            !   molec=.false.
               partition = self%partition_function(temp, dxi, saha, self%totn)
            !end select
         endif

         do i=1, self%nlevel
            where(self%totn == 0)
               self%nstar(:,:,:,i) = 0
            elsewhere
               self%nstar(:,:,:,i) = self%lte_pop(i, temp, dxi, saha, partition)
            endwhere

            zero_pop = .false.
            do x=1, nx
            do y=1, ny
            do z=1, nz
               if (self%nstar(x,y,z,i) <= tiny(0.0) .and. .not. self%totn(x,y,z)==0) then
                  zero_pop = .true.
               endif
            end do
            end do
            end do
            if(zero_pop) then
               print *, ' WARNING: ltepop: level ',i,' has zero lte population.'
               !stop
            end if
         end do

         ! scale total abundances if treated model is a molecule
         !if(molec) then
         !   self%totn = partition * self%nstar(:,:,:,1)
         !endif

         do i=1, self%nlevel
            self%n(:,:,:,i) = self%nstar(:,:,:,i) 
         end do
      endif
   END SUBROUTINE compute_lte_populations
   !------------------------------------------------------------------------------------

   !====================================================================================
   elemental function lte_pop(self, i_lvl, temp, dxi, saha, partition) result(nstar)
      class(atom_type), intent(in) :: self
      real(4), intent(in) :: temp, dxi, saha, partition
      integer, intent(in) :: i_lvl
      real(8) :: evi, tns
      real(4) :: nstar
      logical :: molec
      integer :: i, j, l
      !---------------------------------------------------------------------------------
      j = self%ion(i_lvl)
      evi  = (self%ev(i_lvl) - self%ev(1)) - dxi * 0.5*(j-1)*j
      tns = log(self%g(i_lvl)) - e_k/temp*evi

      if(j > self%ion(1)) then
         l = j - self%ion(1)
         tns = tns + float(l)*saha
      endif

      if (tns+log(partition) < -75) then
         !nstar = 0d0
         nstar = tiny(1.0)
      else
         nstar = exp(tns) * partition
      endif
   end function lte_pop
   !------------------------------------------------------------------------------------
  
   !====================================================================================
   elemental function partition_function(self, temp, dxi, saha, totn) result(sumn)
      class(atom_type), intent(in) :: self
      real(4), intent(in) :: temp, dxi, totn, saha
      real(8) :: evi, tns
      real(4) :: sumn
      integer :: i, j, l
      !--------------------------------------------------------------------------------
      sumn = 0.
      do i=1, self%nlevel
         j = self%ion(i)
         evi  = (self%ev(i) - self%ev(1)) - dxi * 0.5*(j-1)*j

         tns = log(self%g(i)) - e_k/temp*evi

         if(j > self%ion(1)) then
            l = j - self%ion(1)
            tns = tns + float(l)*saha
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
   SUBROUTINE compute_cont_factors(self, temp, dxi, ne)
      class(atom_type), target, intent(inout) :: self
      real(4), dimension(:,:,:) :: temp, dxi, ne
      real(4) :: dev, tns, test, saha, lalpha
      integer :: i, j, k, kr, n(3), x, y, z, ii, jj
      class(cont_type), pointer :: cont
      !---------------------------------------------------------------------------------
      n = shape(self%totn)
      do kr=1,self%ncont
         cont => self%cont(kr)
         call os%allocate(cont%ni_over_nj_lte, n, label='cont%ni_over_nj_lte')
         i = cont%i
         j = cont%j
         ii = 0.5*(self%ion(i)*(self%ion(i)-1))
         jj = 0.5*(self%ion(j)*(self%ion(j)-1))
         lalpha = log(cont%alpha0)
         do z=1, n(3)
         do y=1, n(2)
         do x=1, n(1)
         saha = 1.5*log(temp(x,y,z)) - log(saha_const*ne(x,y,z)) ! actually log(saha)
         dev = self%ev(i) - self%ev(j) - dxi(x,y,z)*(ii - jj)
         tns = log(self%g(i)/self%g(j)) - e_k/temp(x,y,z)*dev
         tns = tns + (self%ion(i) - self%ion(j))*saha
         cont%ni_over_nj_lte(x,y,z) = exp(tns + lalpha)
         !cont%ni_over_nj_lte(x,y,z) = self%nstar(x,y,z,i) / self%nstar(x,y,z,j)
         end do
         end do
         end do
      enddo
   END SUBROUTINE compute_cont_factors
   !------------------------------------------------------------------------------------

   !====================================================================================
   SUBROUTINE read (self, atom_file, qnorm, abnd, logf_cut, red_cut, blue_cut, in_lines, &
                     read_coll)
      !---------------------------------------------------------------------
      !
      !  reads an atomic model.
      !
      !---------------------------------------------------------------------
      class(atom_type), target :: self
      character(len=*), intent(in) :: atom_file
      real(8), intent(in) :: qnorm
      real(8), intent(in), optional :: abnd, logf_cut, red_cut, blue_cut
      integer, intent(in), optional, allocatable :: in_lines(:)
      logical, intent(in) :: read_coll

      ! Declaration of the local variables
      integer :: i, j, ios, latom, strl, kr, itemp, iterm, nterm, ic, nc, pos
      integer :: ncoll, max_n_temp, doubles, ngrids, ii, jj, tri, ntemp, l, u
      character(len = 80) :: text, mulfmt
      character(len = 16) :: key
      character(len = 9), parameter :: rname='read_atom'
      character(len = :), allocatable :: myfmt ! temporary format string
      type(line_type), pointer :: lin => null()
      type(cont_type), pointer :: cnt => null()
      type(coll_type), pointer :: col => null()
      type(line_type),dimension(:),allocatable :: scratch_line
      real(8), allocatable, dimension(:) :: nu, alpha, tgrid
      real(8) :: scr, gv, gq, cnt_ev
      character(len=20) :: txt
      character(len=26),parameter :: upcase_alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      logical :: abort
      !---------------------------------------------------------------------
      call os%check_file(atom_file)
      open(newunit=latom, file=atom_file, form='formatted')
      abort = .false.

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
      call os%mem_register(storage_size(self%line), self%nline, label='atom%line')
      call os%mem_register(storage_size(self%cont), self%ncont, label='atom%cont')

      allocate(  self%ev    ( self%nlevel) ) ; self%ev = 0.0d0
      allocate(  self%g     ( self%nlevel) ) ; self%g = 0
      allocate(  self%label ( self%nlevel) ) ; self%label = ' '
      allocate(  self%ion   ( self%nlevel) ) ; self%ion = 0

      call os%allocate(self%cdata, self%nlevel*(self%nlevel-1)/2, zero=.true.)
      call os%allocate(self%ecol , self%nlevel*(self%nlevel-1)/2, zero=.true.)
      call os%allocate(self%hcol , self%nlevel*(self%nlevel-1)/2, zero=.true.)

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
         abort = .true.
         !stop
         endif
      end do

      do i=1,self%nline
         lin => self%line(i)

         call nextentry(latom) 

         ! io WAS present for historical reasons, in MULTI io determines
         ! whether a line is wide or not. In the RH formulation all lines are wide

         read(latom,*) jj, ii, lin%f , lin%nnu, lin%qmax, lin%q0, lin%io, &
               lin%ga, gv, gq, lin%profile_type

         lin%i = min(ii, jj)
         lin%j = max(ii, jj)

         tri = (lin%j-1)*(lin%j-2)/2 + lin%i
         self%cdata(tri) = lin%f

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
               write(stdout, *) uppercase(rname, len_trim(rname)), &
               ' found ' ,nterm, 'components in transition KR = ',i

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
         lin%Aji   = f_to_Aji * lin%f * lin%Grat / lin%lambda0**2 * hh
         lin%Bji   = c2_2h / lin%nu0**3 * lin%Aji
         lin%Bij   = lin%Bji / lin%Grat

         !if (present(teff)) then
         !   if (teff > 0) then
         !      lin%qmax = lin%estimate_qmax(teff)
         !      lin%q0 = lin%qmax
         !   endif
         !endif

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
         if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1
         !allocate( lin%nu(lin%nnu) )
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

         lin%ntrans = i

         if((trim(lin%profile_type) .ne. 'voigt'  ) .and. &
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
      !  if nu_max.lt.0.0 frequency points in angstrom (starting with
      !  threshold and decreasing) and crossections in cm2 are read.
      do i=1,self%ncont
         cnt => self%cont(i)
         call nextentry(latom) 

         read(latom,*) jj,ii,cnt%alpha0,cnt%nnu,cnt%nu_max

         cnt%i = min(ii, jj)
         cnt%j = max(ii, jj)

         if (self%ion(cnt%i) == self%ion(cnt%j)) then
            print *, 'bound-free transition', i, 'is between two levels of the same ion', cnt%i, cnt%j
            abort = .true.
            !stop
         endif

         tri = (cnt%j-1)*(cnt%j-2)/2 + cnt%i
         self%cdata(tri) = cnt%alpha0

         cnt%lower_ev = self%ev(cnt%i)
         cnt%upper_ev = self%ev(cnt%j)
      
         cnt%nu0   = ( cnt%upper_ev - cnt%lower_ev ) / h_e ! [Hz]
         cnt%nu_max = cc/(cnt%nu_max*AA_TO_CM) ! [from AA to Hz]
         cnt%ntrans=self%nline+i ! transition identifier for continua

         cnt%ired=1
         cnt%iblue=cnt%nnu

         if ( cnt%nu_max .lt. 0.0 ) then
               allocate( nu(cnt%nnu))
               allocate( alpha(cnt%nnu))

               ! cross sections given in detail
               cnt%bf_type='detail'
               do j=1,cnt%nnu
               call nextentry(latom) 
               read(latom,*)  nu(j), alpha(j)
               end do
               ! for overflow reasons
               alpha = alpha / cnt%alpha0
               ! sanity check
               doubles = 0
               do j=2,cnt%nnu
                  if ( nu(j) .gt. nu(j-1) ) then
                        write(stdout,*) 'read_atom: continuum transition nr', i, &
                                       ' wavelengths not decreasing'
                        abort = .true.
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
               cnt%nu_max = cnt%nu(cnt%nnu)

         else
               ! hydrogenic transitions
               cnt%bf_type='hydrogenic'
         endif
      end do

      ! read collisional data
      if (.not. read_coll) then
         close(latom)
         return
      endif

      ! read choice of collisional SUBROUTINE
      ! GENCOL is currently the only collision subroutine.
      call nextentry(latom)
      read(latom,fmt='(a)') self%crout

      ! sanity check
      self%crout = adjustl(self%crout)
      strl=len_trim(self%crout)
      self%crout = lowercase(self%crout,strl)
      if (self%crout /= 'gencol') then
         write(stdout, *) 'Collision routine ', trim(self%crout), ' not supported. Stopping...'
         abort=.true.
      endif


      write(stdout, *) ' - Reading collision data'
      ncoll = 0
      max_n_temp = 8
      ngrids = 0
   
      ! first run through to determine
      ! ncoll = number of keywords
      ! mtgrd = size of the largest temperature grid
   
      do while (.true.)
         call nextentry(latom)
         read(latom, '(a)', iostat=ios) key
         if (ios<0) then
            stop 'read_collisions: keyword "END" not found in collision file'
         endif
         key = uppercase(adjustl(key),len_trim(key))
   
         ! read next non-comment line if key doesnt start with a letter
         if (scan(key(1:1),upcase_alphabet)==0) cycle 
   
         ! we found a keyword, so interpret it
         select case(trim(key))
         case('END')
            exit
   
         case('TEMP')
            ! keep track of max size of temperature grid
            call nextentry(latom)
            read(latom,*) ntemp
            if (ntemp > max_n_temp) max_n_temp = ntemp
            ngrids = ngrids + 1
         case default
            if (.not.is_valid_process(key)) then
               write(stdout, *) 'read_collisions: unrecognized keyword '//key
               abort = .true.
            end if
            ncoll = ncoll + 1
         end select
      end do

      allocate( self%coll(ncoll) )
      call os%mem_register(storage_size(self%coll), ncoll, label='atom%coll')
      call os%allocate(self%tgrid, [max_n_temp, ngrids],   label='atom%tgrid')
      call os%allocate(self%cgrid, [max_n_temp, ncoll],    label='atom%cgrid')
      call os%allocate(self%dcdt,  [max_n_temp, ncoll],    label='atom%dcdt')
      allocate(tgrid(max_n_temp))
   
      rewind(latom)
      
      do while (.true.)
         call nextentry(latom)
         read(latom,*) txt
         if (lowercase(adjustl(txt), len_trim(txt)) == self%crout) exit
      enddo


      itemp = 0
      i=0
      do while (i<=ncoll) !le to catch the end keyword
         !  read the keyword 
         call nextentry(latom)
         read(latom,fmt='(a8)',iostat=ios) key
         if (ios<0) then
            write(stdout, *) 'read_collisions: keyword "END" not found in collision file'
            stop
         endif
   
         key = uppercase(adjustl(key),len_trim(key))

         if (trim(key) == 'END') then
            write(*,*) "read_collisions: read all collisions."
            exit
         else
            call nextentry(latom)
         endif
         
         select case(trim(key))
         case('TEMP')
            itemp = itemp + 1
            read(latom,*) ntemp, tgrid(1:ntemp)
            self%tgrid(1:ntemp, itemp) = tgrid(1:ntemp)
            if (ntemp < 2) then
               stop "read_collisions: ntemp < 2"
            endif
   
         ! MB  5.04.2018 - added CH0 keyword
         ! RH 25.03.2022 - added CHI keyword
         case('OHMEGA','CE','CH','CI','FIXRAT','CH0','CHI')
            nc = ntemp
         case('AR85-CEA','BURGESS') 
            nc = 1
         case('AR85-RR') 
            nc = 2
         case('AR85-CH','AR85-CH+','AR85-CHE','AR85-CHP','AR85-CHE+') 
            nc = 6
         case('SHULL82') 
            nc = 8
   
         case default
            if (.not.is_valid_process(key)) then
               write(stdout,*) 'read_collisions: unrecognized keyword ', key
               stop
            endif
         end select

         ! sanity check on levels
         if (is_valid_process(key)) then
            i=i+1
            col => self%coll(i)
            col%key = key
            col%nc = nc
            read(latom,*) l, u, self%cgrid(1:nc,i)
            
            if (l==u) then
               write(stdout, *) "ll==ul for process ", key, l, u
               abort = .true.
            endif

            select case(trim(key))
            case('OHMEGA','CE','CH','CI','FIXRAT','CH0','CHI')
               ! set derivatives for 3rd order hermite interpolation
               call hermite_derivative(self%dcdt(1:ntemp,i), tgrid(1:ntemp), self%cgrid(1:ntemp,i))
               col%itemp = itemp
            end select

            col%ll = min(l, u)
            col%ul = max(l, u)

            tri = (col%ul-1)*(col%ul-2)/2 + col%ll
            !self%cdata(tri) = -1.0 ! flag non-classical recipe

            select case(trim(key))
            case ('CE', 'CI')
               self%ecol(tri) = .true.
            case ('CH', 'CHI')
               self%hcol(tri) = .true.
            end select
   
            ! within level range
            if(col%ll < 0 .or. col%ul > self%nlevel) then
               write(stdout,*) 'key=', col%key
               write(stdout,*) 'll=',  col%ll 
               write(stdout,*) 'ul=',  col%ul 
               write(stdout,*) 'nk=',  self%nlevel
               write(stdout,*) 'read_atom:  level index outside range (1,nk)'
               abort = .true.
            endif
         end if
      end do
      close(latom)
      if (abort) stop
   END SUBROUTINE read

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
      real(8) :: dnu
      integer, save:: itimer=0
      !-----------------------------------------------------------------------------
      call trace%begin('atom_type%create_spectrum', itimer=itimer)
      ! compute total number of frequencies
      nnu = 0

      do i=1,self%nline
         lin => self%line(i)
         call lin%freql(self%qnorm, symmetric=.true.)
         nnu = nnu + lin%nnu
      end do
      do i=1,self%ncont
         cnt => self%cont(i)
         if (cnt%bf_type=='hydrogenic') then
            allocate( cnt%nu(1:cnt%nnu)    )
            allocate( cnt%alpha(1:cnt%nnu) )
            dnu = (cnt%nu_max - cnt%nu0) / (cnt%nnu-1)
            cnt%nu = (/ (cnt%nu0 + i*dnu, i=0, cnt%nnu-1) /)
            cnt%alpha = (cnt%nu / cnt%nu0 )**-3
         endif
         nnu = nnu + cnt%nnu
      end do

      ! save all frequencies
      allocate(nu(nnu))
      je=0
      do i=1,self%nline
         js=je+1
         je=js+self%line(i)%nnu-1
         nu(js:je)=self%line(i)%nu
         deallocate(self%line(i)%nu)
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
      call os%allocate(sp%nac, sp%nnu, zero=.true.)
      call os%allocate(sp%nal, sp%nnu, zero=.true.)

      do i=1,self%ncont
         cnt => self%cont(i)
         ! find minimum and maximum index of this continuum in the sp%nu grid
         call get_interval(sp%nu, cnt%nu0, cnt%nu_max, ired, iblue, limit=.false.)
         if (iblue == 0 .or. ired == sp%nnu) cycle
      
         ! add to number of active continua
         do j=ired,iblue
            sp%nac(j) = sp%nac(j) + 1 ! number of active continua at frequency index j
         end do
      end do

      do i=1,self%nline
         lin => self%line(i)
         ! find minimum and maximum index of this line in the sp%nu grid
         call get_interval(sp%nu, lin%nu_min, lin%nu_max, ired, iblue, limit=.false.)
         if (iblue == 0 .or. ired == sp%nnu) cycle
      
         ! add to number of active lines
         do j=ired,iblue
            sp%nal(j) = sp%nal(j) + 1
         end do
      end do

      call os%allocate(sp%ac, [maxval(sp%nac,1), sp%nnu], zero=.true., label='sp%ac')
      call os%allocate(sp%al, [maxval(sp%nal,1), sp%nnu], zero=.true., label='sp%al')
      sp%nac = 0
      sp%nal = 0
      
      do i=1,self%ncont
         cnt => self%cont(i)

         ! find minimum and maximum index of this continuum in the sp%nu grid
         call get_interval(sp%nu, cnt%nu0, cnt%nu_max, ired, iblue, limit=.false.)

         ! total number of frequencies where this continuum is active
         nnu=iblue-ired+1

         if (iblue == 0 .or. ired == sp%nnu .or. nnu == 1) then
            cnt%ired  = 0
            cnt%iblue = 0
            if (allocated(cnt%nu)) deallocate(cnt%nu)
            if (allocated(cnt%alpha)) deallocate(cnt%alpha)
            cnt%active = .false.
            cycle
         endif

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

            deallocate(cnt%nu)
            deallocate(cnt%alpha)
            allocate(cnt%nu(ired:iblue))
            allocate(cnt%alpha(ired:iblue))
            cnt%nu = sp%nu(ired:iblue)
            cnt%alpha = (cnt%nu / cnt%nu0 )**-3

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
            call get_interval(sp%nu, lin%nu_min, lin%nu_max, ired, iblue, limit=.false.)

            ! total number of frequencies where this line is active
            nnu=iblue-ired+1

            if (iblue == 0 .or. ired == sp%nnu .or. nnu==1) then
               lin%ired  = 0
               lin%iblue = 0
               if (allocated(lin%nu)) deallocate(lin%nu)
               if (allocated(lin%q)) deallocate(lin%q)
               lin%active = .false.
               cycle
            endif

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
            !deallocate( lin%nu      )
            allocate(lin%nu(nnu))
            allocate(lin%q(nnu))
            !allocate(   lin%nu(ired:iblue) )
            !allocate(   lin%q(ired:iblue) )
            lin%nu=sp%nu(ired:iblue) ! [Hz]
            !lin%q = (lin%nu/lin%nu0-1.0d0) * CC / (self%qnorm * KM_TO_CM) ! [Doppler units]
            lin%q = lin%nu/lin%nu0 - 1.0d0

            ! fill frequency integration weights for the line
            call sp%set_weights(lin%wnu, lin%nu, ired, iblue)
            call sp%set_weights(lin%wq , lin%q , ired, iblue)
      end do

      sp%maxac=maxval(sp%nac)
      sp%maxal=maxval(sp%nal)

      write(*,'(a,I0,a)') 'sp_from_atom: frequency grid has ', sp%nnu, ' points.'
      call trace%end(itimer)
   END SUBROUTINE set_active_transitions

   !===============================================================================
   SUBROUTINE unique_int(arr, out, nout)
      integer, intent(in), dimension(:) :: arr
      integer, intent(out), dimension(:) :: out
      integer, intent(out) :: nout
      integer :: i, j, n
      logical :: new
      !-----------------------------------------------------------------------------
      n = size(arr)
      nout = 0
      do i=1, n
         new = .true.
         do j=1, nout
            if (arr(i) == out(j)) new = .false.
         end do

         if (new) then
            nout = nout+1
            out(nout) = arr(i)
         endif
      end do
   END SUBROUTINE unique_int

   !===============================================================================
   SUBROUTINE copy_to (self, other)
      class(atom_type) :: self
      type(atom_type) :: other
      real(8), dimension(:,:), pointer :: tgrid, cgrid, dcdt
      type(coll_type), dimension(:), pointer :: col
      !-----------------------------------------------------------------------------
      tgrid => self%tgrid
      cgrid => self%cgrid
      dcdt  => self%dcdt
      col   => self%coll
      nullify(self%tgrid)
      nullify(self%cgrid)
      nullify(self%dcdt)
      nullify(self%coll)

      ! deep copy
      other = self

      self%tgrid => tgrid
      self%cgrid => cgrid
      self%dcdt  => dcdt 
      self%coll  => col
      other%tgrid => tgrid
      other%cgrid => cgrid
      other%dcdt  => dcdt 
      other%coll  => col
      
      call os%mem_register(storage_size(self%line), self%nline, label='atom%line')
      call os%mem_register(storage_size(self%cont), self%ncont, label='atom%cont')
   END SUBROUTINE

END MODULE atom_mod
  
  ! *******************************************************************
  