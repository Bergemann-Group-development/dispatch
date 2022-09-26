! *******************************************************************
MODULE atom_mod
    USE m3io_mod
    USE mconst
    USE coll_mod
    USE trans_mod
    USE mmath, only: refsor, unista
    USE op_math_mod, only: searchsorted, lineq, lineq_no_pivot
    USE ion_mod, only: ion_type
    USE element_mod, only: element_type
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
        integer, dimension(:,:),  allocatable :: ac,al
        ! how many continua (nac) and lines (nal) are active at each frequency
        integer, dimension(:),  allocatable :: nac,nal
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
   
        !integer,dimension(:,:), allocatable :: ilin,icon ! ilin(l1,l2) is the index of the line
        !                                                   ! between l1 and l2, if any, 
        !                                                   ! likewise for icon
   
        real(8) :: abnd ! abundance ( in the [hydrogen]=12 scale)
        real(8) :: awgt ! atomic weight [g]
   
        real(8), dimension(:), allocatable :: ev ! energy of level above ground state [eV]
        real(8), dimension(:), allocatable :: g ! multiplicity of level
  
        type(line_type), dimension(:), allocatable :: line
        type(cont_type), dimension(:), allocatable :: cont
        type(coll_type), dimension(:), allocatable :: coll
        type(atomic_spectrum_type), pointer :: sp
  
        ! rather constant doppler width factor, i.e. dnud = line(i)%nu0 * doppler_factor [Hz]
        real(4), dimension(:,:,:), pointer, contiguous :: doppler_factor
  
        ! LTE, NLTE and total population of active atom
        real(4), dimension(:,:,:,:), pointer, contiguous :: nstar => null()
        real(4), dimension(:,:,:,:), pointer, contiguous :: n     => null()
        real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij => null()
        real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ji => null()
        real(4), dimension(:,:,:), pointer, contiguous :: totn
  
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
     
    end type atom_type
  
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
        real(4), dimension(:,:,:,:), intent(in), contiguous :: Psi, Ieff
        real(4), dimension(:,:,:,:), intent(inout), contiguous :: gamma_ij, gamma_ji
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
                 alpha2 = alpha! * self%cont(kc)%alpha0
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
      integer, save :: itimer=0
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

      if (check_NaN) then
         do i=1,self%nlevel
         do j=1,self%nlevel
            if (self%w(i,j) /= self%w(i,j)) then
               print *, 'NaN in rate Matrix at i,j = ', i, j
               print *, 'Temp = ', temp
               stop
            endif
         end do
         end do
      endif

      ! solve w n . x = e
      ! e contains the solution vector x on exit
      !call lineq(self%w, self%e, self%nlevel)
      call lineq_no_pivot(self%w, self%e, self%nlevel)

      ! set limits on corrections
      where(self%e < 0.001d0) self%e = 0.001d0

      ! update the populations
      n = n * self%e

      ! set e to the relative change in the populations 
      self%e = self%e - 1.0d0
   END SUBROUTINE solve_stat_equilibrium
   !------------------------------------------------------------------------------------
  
   !====================================================================================
   SUBROUTINE compute_lte_populations(self, temp, ne, dxi, element)
      class(atom_type), intent(inout) :: self
      class(element_type), pointer :: element
      class(ion_type), pointer :: ion
      integer :: i, j
      logical :: molec
      real(4), dimension(:,:,:), intent(in) :: temp, ne, dxi
      real(4), dimension(:,:,:), allocatable :: partition, saha
      integer :: nx, ny, nz
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

         select case (self%id(1:len_trim(self%id)))
         case ('ch','co','cn')
            molec=.true.
            !partition = totn * self%g(1) * exp(-e_k/temp*self%ev(1)) / qcalc(elemid, temp)
            print *, 'molecule partition functions need to be updated. Stopping...'
            stop
         case default
            molec=.false.
            partition = self%partition_function(temp, dxi, saha, self%totn)
         end select

         do i=1, self%nlevel
            self%nstar(:,:,:,i) = self%lte_pop(i, temp, dxi, saha, partition)

            if(any(self%nstar(:,:,:,i) <= tiny(0.0))) then
               print *, ' WARNING: ltepop: level ',i,' has zero lte population.'
               !stop
            end if
         end do

         ! scale total abundances if treated model is a molecule
         if(molec) then
            self%totn = partition * self%nstar(:,:,:,1)
         endif

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
            !saha = 1.5*log(temp(x,y,z)) - log(saha_const*ne(x,y,z)) ! actually log(saha)
            !dev = self%ev(i) - self%ev(j) - dxi(x,y,z)*(ii - jj)
            !tns = log(self%g(i)/self%g(j)) - e_k/temp(x,y,z)*dev
            !tns = tns + (self%ion(i) - self%ion(j))*saha
            !cont%ni_over_nj_lte(x,y,z) = exp(tns + lalpha)
            cont%ni_over_nj_lte(x,y,z) = self%nstar(x,y,z,i) / self%nstar(x,y,z,j)
           end do
           end do
           end do
        enddo
     END SUBROUTINE compute_cont_factors
     !------------------------------------------------------------------------------------
  
     !====================================================================================
     SUBROUTINE read (self, atom_file, qnorm, abnd, logf_cut, red_cut, blue_cut, in_lines)
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
        real(8), intent(in), optional :: abnd, logf_cut, red_cut, blue_cut
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
  
        !allocate(  self%ilin  ( self%nlevel,self%nlevel) ) ; self%ilin=0
        !allocate(  self%icon  ( self%nlevel,self%nlevel) ) ; self%icon=0
  
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
            lin%Aji   = f_to_Aji * lin%f * lin%Grat / lin%lambda0**2! * hh
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
            !self%ilin(lin%j,lin%i) = i
            !self%ilin(lin%i,lin%j) = i
  
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
            !self%icon( cnt%j, cnt%i ) = i
            !self%icon( cnt%i, cnt%j ) = i
  
            if ( cnt%numax .lt. 0.0 ) then
  
                allocate( nu(cnt%nnu))
                allocate( alpha(cnt%nnu))
  
                ! cross sections given in detail
                cnt%bf_type='detail'
                do j=1,cnt%nnu
                 call nextentry(latom) 
                 read(latom,*)  nu(j), alpha(j)
                end do
                ! for overflow reasons
                !alpha = alpha / cnt%alpha0
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
  
  END MODULE atom_mod
  
  ! *******************************************************************
  