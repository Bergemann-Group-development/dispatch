! *******************************************************************
MODULE atom_mod
    USE io_mod
    USE m3io_mod
    USE handyfuncs_mod
    USE prec_mod, only: inf8
    USE mconst
    USE coll_mod
    USE trans_mod
    USE mmath, only: refsor, unista
    USE op_math_mod, only: lineq, get_interval, eqsyst, lineq_no_pivot
    USE ion_mod, only: ion_type, massive_type
    USE element_mod, only: element_type
    USE molecule_mod, only: molec_type
    USE os_mod
    USE trace_mod
    USE spectrum_mod, only: spectrum_type, set_weights
    USE hprofl_mod, only:prepare_hydrogen_damping
    ! *******************************************************************
    !
    !     Contains variables connected to the atom(s) 
    !     included in the calculations
    !
    ! *******************************************************************
    implicit none
    private

    !-------------------------------------------------------------------------------
    type, public, extends(spectrum_type) :: atomic_spectrum_type
      ! maximum number of active lines or continua at any given frequency
      integer :: maxal=0, maxac=0, nmax=0, nlower, nupper
      ! which continua (ac) and lines (al) are active at each frequency
      integer, dimension(:,:), contiguous, pointer :: ac,al
      ! how many continua (nac) and lines (nal) are active at each frequency
      integer, dimension(:), contiguous, pointer :: nac,nal
      ! has this frequency an active prd line
      !logical,dimension(:),allocatable :: hasPRDline
    end type atomic_spectrum_type
    !---------------------------------------------------------------------

    type, public :: local_data_t
      ! LTE, NLTE and total population of active atom
      real(4), dimension(:,:,:,:), pointer, contiguous :: nstar => null()
      real(4), dimension(:,:,:,:), pointer, contiguous :: n     => null()
      real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij => null()
      real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ji => null()
      real(4), dimension(:,:,:,:), pointer, contiguous :: ni_over_nj_lte => null()
      real(4), dimension(:,:,:,:), pointer, contiguous :: wphi => null()
      real(4), dimension(:,:,:),   pointer, contiguous :: totn => null() 
    end type local_data_t
  
    type, public, extends(massive_type) :: atom_type
      integer :: nrad   ! total number of radiative transitions
      integer :: nrfix  ! number of fixed transitions
      integer :: ncont  ! number of bound-free transitions
      integer :: nline  ! number of bound-bound transitions
      integer :: nlevel ! number of levels
      integer :: nprd=0 ! number of PRD transitions
      
      character(len=128) :: name
      character(len=3)   :: id    ! atom id string 
      character(len=20)  :: crout ! which collisional SUBROUTINE to use 
  
      character(len=20), dimension(:), allocatable :: label ! level label
  
      integer, dimension(:), allocatable :: ion ! ionization stage of each level
  
      real(8) :: abnd ! abundance ( in the [hydrogen]=12 scale)
      !real(8) :: awgt ! atomic weight [g]
      real(8) :: tmax_update=1.0d7
  
      real(8), dimension(:), allocatable :: ev ! energy of level above ground state [eV]
      real(8), dimension(:), allocatable :: g ! multiplicity of level

      type(line_type), dimension(:), pointer, contiguous :: line
      type(cont_type), dimension(:), pointer, contiguous :: cont
      type(coll_type), dimension(:), pointer, contiguous :: coll
      type(atomic_spectrum_type), pointer :: sp
      type(local_data_t), pointer :: data

      real(8), dimension(:,:), pointer, contiguous :: tgrid, cgrid, dcdt
      real(4), dimension(:), pointer, contiguous :: cdata
      logical, dimension(:), pointer, contiguous :: ecol, hcol

      real(4), dimension(:),   allocatable  :: oldn
      real(8), dimension(:),   allocatable  :: e
      real(8), dimension(:,:), allocatable  :: C, W
      
      real(8) :: qnorm
      integer :: fcheck, pcheck(3), linecheck
      logical :: stiff 

      real(8) :: drawin, mu_e
        
    contains
      procedure :: add_continuum_opacities
      procedure :: partition_function
      procedure :: lte_pop
      procedure :: compute_lte_populations
      procedure :: compute_cont_factors
      procedure :: solve_stat_equilibrium
      procedure :: create_spectrum
      procedure :: set_ired_iblue
      procedure :: set_active_transitions
      procedure :: get_active_nlevels
      procedure :: read
      procedure :: prepare_Drawin
      procedure :: classic_col
      procedure :: copy_to
      procedure :: link_transitions
    end type atom_type

    logical, public :: use_classic_col = .False., m1d_legacy_mode=.False.
    logical, public :: lse_no_pivot = .True., kramer_bf=.False., debye_shield=.True.
    real(8), public :: red_thresh_aa  = 20000
    real(8), public :: blue_thresh_aa = 500
    real(8), public :: cont_step_aa = 50
    real(8), public :: nu_tolerance = 1d6
    real(8), public :: logf_cut = -99.9, red_cut=huge(1.0), blue_cut=tiny(1.0)
    real(8), parameter :: ame=5.485799e-4, amh=1.0
    real(8), parameter :: mu_e=ame/(ame+amh)
    character(len=32), public :: lse_solver='lineq_no_pivot'
    integer, public :: nq_fac=1

    public :: unique_int
  
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
      if (.not. cnt%active) cycle
      call cnt%check_range(fcs,fce,fs,fe,skip)
      if (skip) cycle
  
      ffs = fs-fs_patch+1
      ffe = fe-fs_patch+1
      fs = fs-cnt%ired+1
      fe = fe-cnt%ired+1

      call cnt%opacities(cnt%alpha(fs:fe), self%data%ni_over_nj_lte(:,:,:,i), &
      self%data%n, self%data%nstar, eta(:,:,:,ffs:ffe), chi(:,:,:,ffs:ffe))

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
    integer, save :: itimer=0, itimer2=0
    logical :: zero_pop, abort
    logical, save :: first_time=.true.
    !------------------------------------------------------------------------------------
    if (temp > self%tmax_update) then
        self%e = 0.0d0
        return
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

    if (cswitch > 0) then
      call trace%begin('gencol', itimer=itimer)
      self%c = 0.d0

      if (use_classic_col) then
        call self%classic_col(temp,ne,nh,nstar)
      endif

      call gencol(self%coll,self%id,self%nlevel,self%g,self%ev,self%ion,&
                  self%c,dble(temp),dble(ne),nstar,dble(nh),self%tgrid,self%cgrid,self%dcdt)

      call trace%end(itimer=itimer)

      !open(newunit=kr, file='debug.bin', form='unformatted')
      !write(kr) nstar
      !write(kr) n
      !write(kr) totn
      !write(kr) self%w
      !write(kr) self%c

      do j=1, self%nlevel
      do i=1, self%nlevel
        self%w(i,j) = self%w(i,j) + self%c(i,j) * cswitch
      end do
      end do
    endif

    ! fill diagonal elements of rate coefficient matrix
    do i=1, self%nlevel
      self%w(i,i) = -sum(self%w(:,i))
    end do

    ! replace row with largest population with particle conservation
    ! equation.
    pos = maxloc(n,1)
    self%e(pos) = totn
    !zero_pop = any(n==0)

    !if (zero_pop) then
    !  self%oldn = n
    !  self%w(pos,:) = 1
    !else
      ! go from rate coefficients to actual rates   
      do i=1,self%nlevel
        self%w(i,:) = self%w(i,:) * n 
      enddo
      self%w(pos,:) = n
    !endif

    if (check_NaN .or. first_time) then
      first_time = .false.
      abort = .false.
      do i=1,self%nlevel
        if (self%w(i,i) == 0.d0) then
          print *, 'zero in rate diagonal at i = ', i
          abort = .true.
        endif
        do j=1,self%nlevel
          if (self%w(i,j) /= self%w(i,j)) then
            print *, 'NaN in rate Matrix at i,j = ', i, j
            abort = .true.
          elseif (self%w(i,j) > inf8) then
            print *, 'Infinity in rate Matrix at i,j = ', i, j
            abort = .true.
          endif
        end do
      end do
      if (abort) then
        print *, 'Temp = ', temp
        print *, 'Totn  = ', totn
        stop
      endif
    endif

    call trace%begin('lse_solver', itimer=itimer2)

    ! solve w n . x = e
    ! e contains the solution vector x on exit
    select case (lse_solver)
    case('eqsyst')
      call eqsyst(self%w, self%e, self%nlevel)
    case('lineq_no_pivot')
      call lineq_no_pivot(self%w, self%e, self%nlevel)
    case('lineq')
      call lineq(self%w, self%e, self%nlevel)
    end select

    call trace%end(itimer2)

    if (check_NaN) then
      do i=1,self%nlevel
        if (self%e(i) /= self%e(i)) then
          print *, 'NaN in error vector at i = ', i
          print *, 'Temp  = ', temp
          stop
        endif
      end do
    endif

    !if (zero_pop) then
    !  self%oldn = n
    !  n = self%e
    !  do i=1,self%nlevel
    !    if (self%oldn(i) > 0) then
    !      self%e(i) = n(i) / self%oldn(i)
    !    elseif (n(i) == 0) then
    !      self%e(i) = 1.d0
    !    else
    !      self%e(i) = huge(1.0)
    !    endif
    !  end do
    !  n = max(n, 0.)
    !else
      ! set limits on corrections
      self%e = max(self%e, 0.001d0)
      !self%e = max(self%e, tiny(1.0))
      ! update the populations
      n = n * self%e
    !endif

    if (check_NaN) then
      do i=1,self%nlevel
        if (n(i) < 0) then
          print *, 'Below zero in populations at i = ', i
          print *, 'Temp  = ', temp
          print *, n(i)
          print *, self%e(i)
          !print *, zero_pop, self%oldn(i)
          stop
        endif
      end do
    endif

    !write(kr) self%e
    !write(kr) n
    !close(kr)
    !stop

    ! set e to the relative change in the populations 
    self%e = self%e - 1.0d0
  END SUBROUTINE solve_stat_equilibrium
  !------------------------------------------------------------------------------------
  
  !====================================================================================
  SUBROUTINE compute_lte_populations(self, nstar, totn, temp, saha, dxi, element, mol)
    class(atom_type), intent(inout) :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion
    class(molec_type), pointer :: mol
    integer :: i, j, x,y,z
    real(4), intent(inout), dimension(:,:,:,:) :: nstar
    real(4), intent(inout), dimension(:,:,:) :: totn
    real(4), intent(in), dimension(:,:,:) :: temp, saha, dxi
    real(4), dimension(:,:,:), allocatable :: partition, f
    integer :: ijk(3)
    logical :: zero_pop
    integer, save :: itimer=0
    !---------------------------------------------------------------------------------
    call trace%begin('atom_mod%compute_lte_populations', itimer=itimer)
    allocate(partition, f, mold=temp)
    if (associated(element)) then
      do i=1, self%nlevel
        j = self%ion(i)
        if (j <= ubound(element%ions,1)) then
          ion => element%ions(j)
          ! ion%anjon already takes saha factor into account
          ! it also takes into account reductions in abundances due to molecules
          partition = ion%n * ion%part
          nstar(:,:,:,i) = self%lte_pop(i, temp, dxi, 1.0, partition)
        else
          print *, 'ion stage larger than max ion stage from UOP!'
        endif
      end do
    else
      partition = self%partition_function(temp, dxi, saha, totn)
      if (associated(mol)) then
        nstar(:,:,:,1) = partition ! store partition function from levels
        partition = totn / qjorg(temp)
        !partition = totn / mol%irwin_partf(dble(temp))
        partition = partition * self%g(1)*exp(-self%ev(1)*e_k/temp)
         ! correct for the use of a different partition function
        totn = partition / nstar(:,:,:,1) * totn
      endif

      do i=1, self%nlevel
        nstar(:,:,:,i) = self%lte_pop(i, temp, dxi, saha, partition)
      end do

      do i=1, self%nlevel
        where(totn == 0) nstar(:,:,:,i) = 0
        ijk = findloc(nstar(:,:,:,i), 0., mask=totn>0)
        if (any(ijk > 0)) then
          print *, ' WARNING: ltepop: level ',i,' has zero lte population at', ijk
          !stop
        end if
      end do
      if (associated(mol)) then
        associate(ndis => nstar(:,:,:,self%nlevel))
          totn = totn - ndis
          ndis = huge(1.e0)
          do j=1, mol%nions
            ion => mol%ions(j)%p
            ndis = min(ndis, ion%n*ion%part)
          end do
          totn = totn + ndis
        end associate
      endif
    endif
    call trace%end(itimer)
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
    evi = self%ev(i_lvl) - self%ev(1) - dxi * 0.5*(j-1)*j
    tns = log(self%g(i_lvl)/self%g(1)) - e_k/temp*evi

    if(j > self%ion(1)) then
      l = j - self%ion(1)
      tns = tns + float(l)*log(saha)
    endif

    tns = tns + log(partition)

    if (tns < -85) then
      nstar = tiny(1.0)
    else
      nstar = exp(tns)
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
      evi = self%ev(i) - self%ev(1) - dxi * 0.5*(j-1)*j

      tns = log(self%g(i)/self%g(1)) - e_k/temp*evi

      if(j > self%ion(1)) then
        l = j - self%ion(1)
        tns = tns + float(l)*log(saha)
      endif

      sumn = sumn + exp(tns)
    enddo 
    sumn = totn / sumn
  end function
  !------------------------------------------------------------------------------------
   
  !====================================================================================
  SUBROUTINE compute_cont_factors(self, temp, dxi, saha)
    class(atom_type), target, intent(inout) :: self
    real(4), dimension(:,:,:) :: temp, dxi, saha
    real(4) :: dev, tns, test, lalpha
    integer :: i, j, k, kr, n(4), x, y, z, ii, jj
    class(cont_type), pointer :: cont
    !---------------------------------------------------------------------------------
    n = [shape(temp), self%ncont]
    call os%allocate(self%data%ni_over_nj_lte, n, label='data%ni_over_nj_lte')

    do kr=1,self%ncont
      cont => self%cont(kr)
      if (.not. cont%active) cycle
      i = cont%i
      j = cont%j
      !ii = 0.5*(self%ion(i)*(self%ion(i)-1))
      !jj = 0.5*(self%ion(j)*(self%ion(j)-1))
      do z=1, n(3)
      do y=1, n(2)
      do x=1, n(1)
        !dev = self%ev(i) - self%ev(j) - dxi(x,y,z)*(ii - jj)
        !tns = log(self%g(i)/self%g(j)) - e_k/temp(x,y,z)*dev
        !tns = tns + (self%ion(i) - self%ion(j))*log(saha(x,y,z))
        !self%data%ni_over_nj_lte(x,y,z,kr) = exp(tns)
        self%data%ni_over_nj_lte(x,y,z,kr) = self%data%nstar(x,y,z,i) / self%data%nstar(x,y,z,j)
      end do
      end do
      end do
    enddo
  END SUBROUTINE compute_cont_factors
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE read (self, atom_file, qnorm, in_lines, read_coll)
    !---------------------------------------------------------------------
    !
    !  reads an atomic model.
    !
    !---------------------------------------------------------------------
    class(atom_type), target :: self
    character(len=*), intent(in) :: atom_file
    real(8), intent(in) :: qnorm
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
    real(8), allocatable, dimension(:) :: nu, alpha, tgrid
    real(8) :: scr, gv, gq, cnt_ev, nu_max, alpha0, thresh_aa, dnu
    real(8) :: gij, red_fcut, blu_fcut, vbar2e4
    integer :: ncont, nnu
    character(len=20) :: txt
    character(len=26),parameter :: upcase_alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    logical :: abort
    !---------------------------------------------------------------------
    call os%check_file(atom_file)
    open(newunit=latom, file=atom_file, form='formatted')
    abort = .false.

    call nextentry(latom)
    read(latom,fmt='(a)') self%id
    self%id=uppercase(adjustl(self%id))

    call nextentry(latom)
    read(latom,*) self%abnd, self%awgt
    call self%prepare_Drawin()

    call nextentry(latom)
    read(latom,*) self%nlevel, self%nline, self%ncont, self%nrfix

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

    allocate(self%line(self%nline)) 
    call os%mem_register(storage_size(self%line), self%nline, label='atom%line')

    do i=1,self%nline
      lin => self%line(i)
      call nextentry(latom) 

      ! io WAS present for historical reasons, in MULTI io determines
      ! whether a line is wide or not. In the RH formulation all lines are wide

      read(latom,*) jj, ii, lin%f , lin%nnu, lin%qmax, lin%q0, lin%io, &
            lin%ga, gv, gq, lin%profile_type

      lin%i = min(ii, jj)
      lin%j = max(ii, jj)
      lin%q0   = lin%q0  
      lin%qmax = lin%qmax

      ! make sure we always have an odd number of frequency points for a line
      if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1

      if (nq_fac /= 0) then
        lin%nnu = nq_fac * (lin%nnu-1) + 1
      endif

      tri = (lin%j-1)*(lin%j-2)/2 + lin%i
      self%cdata(tri) = lin%f

      if (present(in_lines)) then
        if (allocated(in_lines)) then
          if (findloc(in_lines, i, 1) <= 0) then
            call lin%deactivate()
            cycle
          endif
        endif
      endif

      if (logf_cut > -99) then
        scr = log10(self%g(lin%i)*lin%f)
        if (scr < logf_cut) then
          call lin%deactivate()
          cycle
        endif
      endif

      lin%upper_g  = self%g(lin%j)
      lin%lower_ev = self%ev(lin%i)
      lin%upper_ev = self%ev(lin%j)
      lin%hydrogen = trim(self%id) == 'h'

      lin%lambda0 = hc_e / ( lin%upper_ev - lin%lower_ev ) ! [cm]
      lin%nu0 = cc / lin%lambda0    ! [Hz]

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
        write(stdout, *) uppercase(rname), &
        ' found ' ,nterm, 'components in transition KR = ',i

        allocate(lin%determ (nterm) )
        allocate(lin%wterm  (nterm) )
        allocate(lin%fterm  (nterm) )
        allocate(lin%gaterm (nterm) )
        allocate(lin%gwterm (nterm) )
        allocate(lin%gqterm (nterm) )
        !allocate(lin%adterm, xs, xe, ys, ye, zs1, ze1, 1, nterm, rname)

        do iterm = 1, nterm
          read(latom,*) lin%determ(iterm), lin%wterm(iterm),  lin%fterm(iterm), &
                        lin%gaterm(iterm), lin%gwterm(iterm), lin%gqterm(iterm)
        end do

        lin%determ = cc/lin%nu0 * lin%determ
      end if

      !allocate( lin%nu(lin%nnu) )
      lin%nq = lin%nnu

      call os%allocate(lin%nu_input, lin%nnu, 'line%nu_input')
      call os%allocate(lin%wq, lin%nnu, 'line%wq')
      call lin%get_nu(qnorm, lin%nu_input)
      lin%nu_min = lin%nu_input(1)
      lin%nu_max = lin%nu_input(lin%nnu)
      call set_weights(lin%wq, lin%nu_input / lin%nu0 - 1)

      gij = self%g(lin%i) / self%g(lin%j)
      lin%Aji = f_to_Aji * lin%f * gij / lin%lambda0**2 * hh
      lin%Bji = c2_2h / lin%nu0**3 * lin%Aji
      lin%Bij = lin%Bji / gij

      if (red_cut < huge(1.0) .or. blue_cut > tiny(1.0)) then
        red_fcut = cc*1e8 / red_cut
        blu_fcut = cc*1e8 / blue_cut
        if (lin%nu_max < red_fcut .or. lin%nu_min > blu_fcut) then
          call lin%deactivate()
          cycle
        endif
      endif
      
      ! find continuum level
      do ic=lin%j+1, self%nlevel
        if (self%ion(ic) .eq. self%ion(lin%j)+1) then
          cnt_ev = self%ev(ic)
          exit
        endif
      end do

      lin%ntrans = i
      lin%profile_type = trim(lowercase(lin%profile_type))

      select case (lin%profile_type)
      case ('voigt', 'hlinop')
        continue
      case default
        write (*, fmt = '(a,1x,a,a)') 'Profile not recognised:', &
                lin%profile_type, ' Setting it to VOIGT'
        lin%profile_type = 'voigt'
      end select

      select case (lin%profile_type)
      case ('hlinop')
        call prepare_hydrogen_damping(lin)
      case ('voigt')
        call lin%prepare_Stark(gq)
        vbar2e4 = self%compute_vbar2(10000.)
        call lin%prepare_VanderWaals(gv, self%ion(lin%i), cnt_ev, vbar2e4)
      end select
    end do

    ncont = count(self%ion < maxval(self%ion))
    if (ncont < self%ncont) then
        print *, 'Unrealisticly many bound-free transitions. Something is wrong!'
    endif
    allocate(self%cont(ncont))
    call os%mem_register(storage_size(self%cont), ncont, label='atom%cont')

    do i=1, ncont
      cnt => self%cont(i)
      cnt%bf_type='Kramer'
      ! find continuum level
      do ic=i, self%nlevel
        if (self%ion(ic) .eq. self%ion(i)+1) exit
      end do
      cnt%i = i
      cnt%j = ic
      cnt%ntrans = self%nline + i ! transition identifier for continua
    end do

    !  bound-free transitions in detail
    !  if nu_max.lt.0.0 frequency points in angstrom (starting with
    !  threshold and decreasing) and crossections in cm2 are read.
    do i=1,self%ncont
      call nextentry(latom) 
      read(latom,*) jj, ii, alpha0, nnu, nu_max

      if (self%ion(ii) == self%ion(jj)) then
        print *, 'Bound-free transition', i, ' is between two levels of the same ion', ii, jj
        abort = .true.
      endif

      j  = max(ii,jj)
      ii = min(ii,jj)

      cnt => self%cont(ii)

      if (cnt%bf_type /= 'Kramer') then
        print *, 'Multiple bound-free transitions from lower state', ii, 'Not implemented!'
        abort = .true.
      endif
      if (j /= cnt%j) then
        print *, 'Bound-free transition', i, ' has exited state as upper level. Worth checking!'
        cnt%j = j
      endif

      cnt%alpha0 = alpha0
      cnt%nnu    = nnu

      if ( nu_max .lt. 0.0 ) then
        ! cross sections given in detail
        cnt%bf_type='detail'

        call os%allocate(cnt%nu_input, cnt%nnu, 'cont%nu_input')
        call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')

        do j=1,cnt%nnu
          call nextentry(latom) 
          read(latom,*)  cnt%nu_input(j), cnt%alpha_input(j)
        end do
        
        do j=2,cnt%nnu
          if ( cnt%nu_input(j) > cnt%nu_input(j-1) ) then
            write(stdout,*) 'read_atom: continuum transition nr', i, &
                          ' wavelengths not decreasing'
            abort = .true.
          endif
        end do

        ! and convert to frequency [Hz]
        cnt%nu_input = cc / ( cnt%nu_input * AA_TO_CM )
      else
        ! hydrogenic transitions
        cnt%bf_type='hydrogenic'
        cnt%nu_max = cc / ( nu_max * AA_TO_CM )
      endif
    end do

    self%ncont = ncont
    self%nrad = self%nline + self%ncont

    j = 0
    do i=1, self%ncont
      cnt => self%cont(i)

      cnt%lower_ev = self%ev(cnt%i)
      cnt%upper_ev = self%ev(cnt%j)
      cnt%nu0 = ( cnt%upper_ev - cnt%lower_ev ) / h_e ! [Hz]

      if (cnt%bf_type == 'hydrogenic') then
        call os%allocate(cnt%nu_input, cnt%nnu, 'cont%nu_input')
        call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')
        dnu = (cnt%nu_max - cnt%nu0) / (cnt%nnu-1)
        cnt%nu_input = [(cnt%nu0 + i*dnu, i=0, cnt%nnu-1)]
        cnt%alpha_input = (cnt%nu0 / cnt%nu_input)**3.d0 * cnt%alpha0

      elseif (cnt%bf_type == 'Kramer') then
        if (kramer_bf) then
          thresh_aa = cc / ( cnt%nu0 * AA_TO_CM )
          thresh_aa = min(thresh_aa, red_thresh_aa)
          nnu = ceiling((thresh_aa - blue_thresh_aa) / cont_step_aa) + 1
          cnt%nnu = nnu
          call os%allocate(cnt%nu_input, cnt%nnu, 'cont%nu_input')
          call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')
          cnt%nu_input(nnu:2:-1) = (/ (blue_thresh_aa + j*cont_step_aa, j=0, nnu-2) /)
          cnt%nu_input(1) = thresh_aa
          cnt%nu_input = cc / ( cnt%nu_input * AA_TO_CM )
          cnt%alpha_input = Kramer(cnt%nu_input, self%ev(cnt%j) - self%ev(cnt%i))
          cnt%alpha0 = cnt%alpha_input(1)
          j = j+1
        else
          call cnt%deactivate()
          cycle
        endif
      endif

      tri = (cnt%j-1)*(cnt%j-2)/2 + cnt%i
      self%cdata(tri) = cnt%alpha0

      cnt%nu_min = cnt%nu_input(1)
      cnt%nu_max = cnt%nu_input(cnt%nnu)
    end do

    if (j>0) then
      print *, j, " bound-free transitions auto-completed using Kramer's formula!"
    endif

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
    !self%crout = adjustl(self%crout)
    !strl=len_trim(self%crout)
    self%crout = lowercase(adjustl(self%crout))

    select case (self%crout)
    case ('end')
      close(latom)
      return
    case ('gencol')
      ! continue
    case default
      write(stdout, *) 'Collision routine ', trim(self%crout), ' not supported. Stopping...'
      stop
    end select

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
      key = uppercase(adjustl(key))

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
      if (lowercase(adjustl(txt)) == self%crout) exit
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

      key = uppercase(adjustl(key))

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

  !====================================================================================
  SUBROUTINE prepare_Drawin(self)
    class(atom_type), intent(inout) :: self
    real(8) :: drmu, mu_h, ama
    !---------------------------------------------------------------------
    ama = self%awgt
    mu_h = amh/(ama+amh)
    drmu = PI*MH*ama*mu_h
    self%drawin = 1.40758E-15 * sqrt(2.0*BK/drmu) * mu_e/mu_h
  END SUBROUTINE prepare_Drawin

  !====================================================================================
  SUBROUTINE classic_col(self,temp,ne,nh,nstar)
    class(atom_type), intent(inout) :: self
    real(4), intent(in) :: temp, ne
    real(4), dimension(:), intent(in) :: nh
    real(4), dimension(:), intent(in) :: nstar
    real(8) :: ni, nj, nij, cup, cdn, EK_over_temp, divu
    real(8) :: dxi, dev, ex, u0, f, tmp, sqt, fac, nu0
    real(8) :: seaton, allen, regemorter
    integer :: i,j, ii, jj, tri, n
    logical :: cont
    logical, save :: first_time=.true.
    !---------------------------------------------------------------------
    if (debye_shield) then
      dxi = 4.98e-4*5040./temp*sqrt(ne*bk*temp)
    else
      dxi = 0
    endif
    sqt = sqrt(temp)
    allen = 5.465e-11 * Rydk * ne / sqt
    regemorter = 14.5 * allen * Rydk / e_k
    seaton = 1.55e13 * ne / sqt
    EK_over_temp = e_k / temp

    n = 0
    do i=1, self%nlevel
      ni = max(nstar(i), tiny(1.0))
      do j=i+1, self%nlevel
        tri = (j-1)*(j-2)/2 + i

        if (.not. self%ecol(tri) .or. .not. self%hcol(tri)) then
          fac = self%cdata(tri)! / self%g(j)
          cont = self%ion(j) > self%ion(i)
          nj = max(nstar(j), tiny(1.0))
          nij = ni / nj

          dev = self%ev(j) - self%ev(i)
          !if (dev == 0) cycle
          dev = max(dev, 1d-8) ! to avoid division by zero

          if (cont) then
            jj = self%ion(j) * (self%ion(j)-1)
            ii = self%ion(i) * (self%ion(i)-1)
            dev = dev - dxi*0.5*(jj-ii)
          endif

          u0 = dev * EK_over_temp
          ex = exp(-u0)

          ! collisions with electrons
          if (.not. self%ecol(tri)) then
            if (cont) then
              ! ionization by electron collisions
              if (fac == 0) then
                nu0 = dev * EE/HH
                fac = Kramer(nu0, dev)! / self%g(j)
              endif
              tmp = min(0.1*self%ion(i), 0.3)
              cup = seaton * fac * tmp / u0 * ex
            else
              ! excitation of neutrals by electron collisions 
              if (fac > 0) then
                ! Van Regemorter
                f = max(fac, 1e-5)
                tmp = max(0.2*ex, 0.276*expin(u0))
                cup = regemorter * f * tmp / dev
              else
                ! Allen
                cup = allen * ex / self%g(i)
              endif
            endif
            !cup = cup * self%g(j)
            cdn = cup * nij
            self%c(j,i) = self%c(j,i) + cup ! upward rate
            self%c(i,j) = self%c(i,j) + cdn ! downward rate
            n = n + 1
          endif

          ! collisions with hydrogen (Drawin)
          if (.not. self%hcol(tri)) then
            divu = (RydV/dev)**2
            tmp = 2.0/u0
            tmp = (1.0 + tmp) / (1.0 + (mu_e*tmp)**2)
            cup = nh(1) * self%drawin * divu * tmp * ex * sqt
            
            if (.not. cont) then
              f = max(fac, 1e-5)
              cup = cup * f
            endif
            !cup = cup * self%g(j)
            cdn = cup * nij
            self%c(j,i) = self%c(j,i) + cup ! upward rate
            self%c(i,j) = self%c(i,j) + cdn ! downward rate
            n = n + 1
          endif
        endif
        if (cont) exit
      end do
    end do
    if (first_time .and. n>0) then
      first_time = .false.
      write(stdout, *) n, ' collisional transitions filled with classical recipes.'
    endif
  END SUBROUTINE classic_col

  !===============================================================================
  ! initialises the spectrum structure
  !===============================================================================
  SUBROUTINE create_spectrum (self, sp, bb, bf, sort)
    class(atom_type), target :: self
    class(spectrum_type) :: sp
    logical, intent(in), optional :: bb, bf, sort
    type(line_type), pointer :: lin
    type(cont_type), pointer :: cnt
    integer :: i,nnu,js,je
    real(8), dimension(:), allocatable :: nu
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_type%create_spectrum', itimer=itimer)
    ! compute total number of frequencies
    nnu = sp%nnu
    if (option(bb)) nnu = nnu + sum(self%line%nnu, self%line%active)
    if (option(bf)) nnu = nnu + sum(self%cont%nnu, self%cont%active)

    ! save all frequencies
    allocate(nu(nnu))
    if (allocated(sp%nu)) nu(1:sp%nnu) = sp%nu

    je=sp%nnu
    if (option(bb)) then
      do i=1,self%nline
        lin => self%line(i)
        if (.not. lin%active) cycle
        js=je+1
        je=js+lin%nnu-1
        nu(js:je)=lin%nu_input
        lin%ired = js
        lin%iblue = je
      end do
    endif
    if (option(bf)) then
      do i=1,self%ncont
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        js=je+1
        je=js+cnt%nnu-1
        nu(js:je)=cnt%nu_input
        cnt%ired = js
        cnt%iblue = je
      end do
    endif

    ! sort frequencies into ascending order 
    if (option(sort, .true.)) then
      sp%sorted = .true.
      call refsor(nu)

      ! discard duplicates, nnu contains number of unique elements on exit
      !call unista(nu,nnu)
      call remove_close_entries(nu, nu_tolerance)
    endif

    ! now set overall frequency grid
    if (allocated(sp%nu)) deallocate(sp%nu)
    allocate(sp%nu, source=nu)
    sp%nnu = size(nu)

    call trace%end(itimer)
  END SUBROUTINE create_spectrum

  !===============================================================================
  SUBROUTINE set_ired_iblue (self, sp, bb, bf, set_active)
    class(atom_type), target :: self
    class(atomic_spectrum_type), target :: sp
    logical, intent(in), optional :: bb, bf, set_active
    type(line_type), pointer :: lin
    type(cont_type), pointer :: cnt
    integer :: i, j, nnu, ired, iblue
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_type%set_ired_iblue', itimer=itimer)
    if (option(bf)) then
      do i=1, size(self%cont)
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        ! find minimum and maximum index of this continuum in the sp%nu grid
        call sp%get_interval(cnt, ired, iblue)
        if (iblue == 0 .or. ired >= sp%nnu) then
          call cnt%deactivate()
        else
          if (iblue==ired) print *, 'WARNING, single frequency for cont', i
          ! total number of frequencies where this continuum is active
          nnu=iblue-ired+1
          cnt%nnu = nnu
          ! min and max frequency where this continuum is active
          cnt%ired  = ired
          cnt%iblue = iblue

          ! set new frequency grid
          call os%allocate(cnt%alpha, cnt%nnu)
          cnt%nu => sp%nu(ired:iblue)
          
          ! recompute the bf cross sections
          select case(cnt%bf_type)
          case('hydrogenic', 'Kramer')
            ! use nu**-3 dependency
            cnt%alpha = (cnt%nu0 / cnt%nu)**3.d0 * cnt%alpha0
          case('detail')
            ! fill new cross sections 
            call cnt%hermite_interp(cnt%alpha, cnt%nu)
          end select

          call os%allocate(cnt%wnu, cnt%nnu, 'cont%wnu')
          call set_weights(cnt%wnu, cnt%nu)
        endif
      end do
    endif

    if (option(bb)) then
      do i=1, size(self%line)
        lin => self%line(i)
        if (.not. lin%active) cycle
        ! find minimum and maximum index of this continuum in the sp%nu grid
        call sp%get_interval(lin, ired, iblue)
        if (iblue == 0 .or. ired == sp%nnu) then
          call lin%deactivate()
        else
          if (iblue==ired) print *, 'WARNING, single frequency for line', i
          ! total number of frequencies where this line is active
          nnu=iblue-ired+1
          lin%nnu = nnu
          ! min and max frequency where this line is active
          lin%ired  = ired
          lin%iblue = iblue

          lin%nu => sp%nu(ired:iblue)
          call os%allocate(lin%q, lin%nnu)
          lin%q = lin%nu/lin%nu0 - 1.0d0

          ! fill frequency integration weights for the line
          call os%allocate(lin%wnu, lin%nnu, 'line%wnu')
          call set_weights(lin%wnu, lin%nu)
        endif
      end do
    endif

    if (option(set_active, .false.)) then
      call self%set_active_transitions(sp, bb, bf)
    endif
    call trace%end(itimer)
  END SUBROUTINE set_ired_iblue

  !===============================================================================
  SUBROUTINE set_active_transitions (self, sp, bb, bf)
    class(atom_type), target :: self
    class(atomic_spectrum_type), target :: sp
    logical, intent(in), optional :: bb, bf
    type(line_type), pointer :: lin
    type(cont_type), pointer :: cnt
    integer :: i, j, nnu, ired, iblue
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_type%set_active_transitions', itimer=itimer)
    
    ! allocate list of active continua and lines
    call os%allocate(sp%nac, sp%nnu, zero=.true.)
    call os%allocate(sp%nal, sp%nnu, zero=.true.)

    if (option(bf)) then
      do i=1, size(self%cont)
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        sp%nac(cnt%ired:cnt%iblue) = sp%nac(cnt%ired:cnt%iblue) + 1 ! number of active continua at frequency index j
      end do
    endif

    if (option(bb)) then
      do i=1, size(self%line)
        lin => self%line(i)
        if (.not. lin%active) cycle
        sp%nal(lin%ired:lin%iblue) = sp%nal(lin%ired:lin%iblue) + 1
      end do
    endif

    sp%maxac = maxval(sp%nac)
    sp%maxal = maxval(sp%nal)
    sp%nmax  = maxval(sp%nac+sp%nal)

    call os%allocate(sp%ac, [sp%maxac, sp%nnu], zero=.true., label='sp%ac')
    call os%allocate(sp%al, [sp%maxal, sp%nnu], zero=.true., label='sp%al')
    sp%nac = 0
    sp%nal = 0

    if (option(bf)) then    
      do i=1, size(self%cont)
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        ! add to number of active continua
        do j=cnt%ired,cnt%iblue
          sp%nac(j) = sp%nac(j) + 1 ! number of active continua at frequency index j
          sp%ac( sp%nac(j), j ) = i ! the index of this active continuum
        end do
      end do
    endif

    if (option(bb)) then
      do i=1, size(self%line)
        lin => self%line(i)
        if (.not. lin%active) cycle
        ! find minimum and maximum index of this continuum in the sp%nu grid
        call sp%get_interval(lin, ired, iblue)
        ! add to number of active lines
        do j=lin%ired,lin%iblue
          sp%nal(j) = sp%nal(j) + 1
          sp%al( sp%nal(j), j ) = i
        end do
      end do
    endif

    call self%get_active_nlevels(sp)

    write(*,'(a,I0,a)') 'sp_from_atom: frequency grid has ', sp%nnu, ' points.'
    write(*,'(a,I0,a)') 'sp_from_atom: ', sum(sp%nac+sp%nal), ' active transition-frequency points'
    call trace%end(itimer)
  END SUBROUTINE set_active_transitions

  !===============================================================================
  SUBROUTINE get_active_nlevels(self, sp)
    class(atom_type), target :: self
    class(atomic_spectrum_type), intent(inout) :: sp
    integer, dimension(:), allocatable :: lower, upper
    integer :: f, i, j, n, kr, kc, kl
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_type%get_active_nlevels', itimer=itimer)
    sp%nlower = 0
    sp%nupper = 0

    allocate(lower(sp%nmax))
    allocate(upper(sp%nmax))

    do f=1, sp%nnu
      if (f > 1) then
        if (all(sp%al(:,f) == sp%al(:,f-1))) then
          if (all(sp%ac(:,f) == sp%ac(:,f-1))) then
            cycle
          endif
        endif
      endif

      n = sp%nac(f) + sp%nal(f)
      
      do kr=1, n
        if (kr <= sp%nac(f)) then
          kc = sp%ac(kr,f)
          lower(kr) = self%cont(kc)%i
          upper(kr) = self%cont(kc)%j
        else
          kl = sp%al(kr - sp%nac(f),f)
          lower(kr) = self%line(kl)%i
          upper(kr) = self%line(kl)%j
        endif
      end do
      
      call unique_int(lower(1:n), i)
      call unique_int(upper(1:n), j)
      sp%nlower = max(i, sp%nlower)
      sp%nupper = max(j, sp%nupper)
    end do

    deallocate(lower)
    deallocate(upper)
    call trace%end(itimer)
  END SUBROUTINE get_active_nlevels

  !===============================================================================
  ! returns arr without duplicate elements (assumes that all entries are positive!)
  !===============================================================================
  SUBROUTINE unique_int(arr, nout, idx)
    integer, intent(inout), dimension(:) :: arr
    integer, intent(out) :: nout
    integer, intent(out), optional, dimension(:) :: idx
    integer, dimension(:), allocatable, save :: scratch
    !$omp threadprivate (scratch)
    integer :: i, j, n, last
    logical :: new, return_index
    !-----------------------------------------------------------------------------
    return_index = present(idx)
    n = size(arr)

    if (.not. allocated(scratch)) then
      allocate(scratch(n))
    elseif (n /= size(scratch)) then
      deallocate(scratch)
      allocate(scratch(n))
    endif

    nout = 0
    do i=1, n
      new = .true.
      do j=1, nout
        new = arr(i)/=scratch(j)
        if (.not. new) exit
      end do

      if (new) then
        nout = nout+1
        scratch(nout) = arr(i)
        j = nout
      endif

      if (return_index) then
        idx(i) = j
      endif
    end do

    arr(1:nout) = scratch(1:nout)
  END SUBROUTINE unique_int

  !===============================================================================
  SUBROUTINE copy_to (self, other)
    class(atom_type) :: self
    type(atom_type) :: other
    !-----------------------------------------------------------------------------
    ! deep copy
    other = self

    allocate(other%line, source=self%line)
    call os%mem_register(storage_size(self%line), self%nline, label='atom%line')
    
    allocate(other%cont, source=self%cont)
    call os%mem_register(storage_size(self%cont), self%ncont, label='atom%cont')
  END SUBROUTINE

  !===============================================================================
  SUBROUTINE link_transitions (self, cont, line)
    class(atom_type) :: self
    type(cont_type), dimension(:), pointer :: cont
    type(line_type), dimension(:), pointer :: line
    !-----------------------------------------------------------------------------
    self%cont => cont
    self%line => line
    !self%ncont = size(cont)
    !self%nline = size(line)
    !self%nrad  = self%ncont + self%nline
  END SUBROUTINE

  !===============================================================================
  !     calculates partition function for ch
  !     data from u.g. jorgensen et al. a&a 1996, 315, 204
  !
  !: qjorg  98-08-13  new routine: (rouppe van der voort)
  !:        calculates partition functions for ch
  !:
  !:        00-10-25  modifications: (martin asplund)
  !:        constant extrapolation for t<1000 and t>8000
  !:        avoided exceeding array limits for t>8000
  !===============================================================================
  PURE ELEMENTAL FUNCTION qjorg(temp)  
    !     data from jorgensen, table 1.
    !     temperature range t={1000, 8000}
    !     note that function extapolates for temp>8000 and temp<1000
    !     for extrapolation for temp>8500 or temp<1000 the results
    !     differ significant from a polynomial extrapolation (polint)
    real(4), intent(in) :: temp
    real(4), parameter :: pfdata(15) = &
                [203.5,  324.7,  471.2,  646.2,  850.8, 1085.1, 1348.1, &
                1638.7, 1955.3, 2296.1, 2659.3, 3043.2, 3446.2, 3866.8, 4303.7]
    real(4), parameter :: tdata(15) = &
                [1000., 1500., 2000., 2500., 3000., 3500., 4000., &
                 4500., 5000., 5500., 6000., 6500., 7000., 7500., 8000.]
    real(4) :: qjorg, x0, x1, y0, y1, b, a
    integer :: i
    !-----------------------------------------------------------------------------
    !     simple linear interpolation between two nearest points
    do i=1, 14
      if (temp .eq. tdata(i)) then
        qjorg = pfdata(i)
        return
      elseif (temp .gt. tdata(i) .and. temp .lt. tdata(i+1)) then 
        x0 = tdata(i)
        x1 = tdata(i+1)
        y0 = pfdata(i)
        y1 = pfdata(i+1)
      elseif (temp .le. tdata(1)) then
        qjorg = pfdata(1)
        return
      elseif (temp .ge. tdata(15)) then
        qjorg = pfdata(15)
        return
      endif
    end do
    b = y0
    a = (y1 - y0)/(x1 - x0)
    qjorg = a*(temp - x0) + b
  END FUNCTION qjorg

END MODULE atom_mod
  
! *******************************************************************
  