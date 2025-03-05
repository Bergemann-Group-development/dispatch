! *******************************************************************
MODULE atom_mod
    USE io_mod
    USE handyfuncs_mod
    USE prec_mod, only: inf8
    USE const_mod
    USE coll_mod
    USE trans_mod
    USE mmath_mod, only: refsor, unista
    USE op_math_mod, only: lineq, get_interval, eqsyst, lineq_no_pivot
    USE ion_mod, only: ion_type, massive_type, debye_shield, massive_type
    USE element_mod, only: element_type, isotope_type
    USE molecule_mod, only: molec_type
    USE os_mod
    USE trace_mod
    USE spectrum_mod, only: spectrum_type, set_weights, nu_tolerance
    USE hprofl_mod, only:prepare_hydrogen_damping
    USE, intrinsic :: ieee_arithmetic
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

    type, public :: atom_data_t
      ! LTE, NLTE and total population of active atom
      real(8), dimension(:,:,:,:), pointer, contiguous :: nstar => null()
      real(8), dimension(:,:,:,:), pointer, contiguous :: n     => null()
      real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij => null()
      real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ji => null()
      real(4), dimension(:,:,:,:), pointer, contiguous :: nibj => null()
      real(4), dimension(:,:,:,:), pointer, contiguous :: wphi => null()
      real(8), dimension(:,:,:),   pointer, contiguous :: totn => null() 
      real(8), dimension(:,:,:,:,:), pointer, contiguous :: oldn => null()
      real(8), dimension(:),   allocatable  :: e
      real(8), dimension(:,:), allocatable  :: C, W
      logical :: converged=.false.
    end type atom_data_t
  
    type, public, extends(massive_type) :: atom_type
      integer :: nrad   ! total number of radiative transitions
      integer :: nrfix  ! number of fixed transitions
      integer :: ncont  ! number of bound-free transitions
      integer :: nline  ! number of bound-bound transitions
      integer :: nlevel ! number of levels
      integer :: nprd=0 ! number of PRD transitions
      integer :: verbose=0
      
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

      type(line_type), dimension(:), pointer, contiguous :: line => null()
      type(cont_type), dimension(:), pointer, contiguous :: cont => null()
      type(coll_type), dimension(:), pointer, contiguous :: coll => null()
      type(atomic_spectrum_type), pointer :: sp => null()
      type(atom_data_t), pointer :: data => null()
      type(atom_data_t), dimension(:), pointer :: rotated_data => null()

      real(8), dimension(:,:), pointer, contiguous :: tgrid => null()
      real(8), dimension(:,:), pointer, contiguous :: cgrid => null()
      real(8), dimension(:,:), pointer, contiguous :: dcdt  => null()
      real(4), dimension(:), pointer, contiguous :: cdata => null()
      logical, dimension(:), pointer, contiguous :: ecol => null()
      logical, dimension(:), pointer, contiguous :: hcol => null()
      
      real(8) :: qnorm
      real(8) :: drawin, mu_e
        
    contains
      procedure :: check_ionisation_thresholds
      procedure :: lte_pop
      procedure :: compute_lte_populations
      procedure :: compute_cont_factors
      procedure :: solve_stat_equilibrium
      procedure :: create_spectrum
      procedure :: set_ired_iblue
      procedure :: set_active_transitions
      procedure :: get_active_nlevels
      procedure :: read
      procedure :: read_pandas
      procedure :: prepare_damping
      procedure :: prepare_Drawin
      procedure :: classic_col
      procedure :: copy_to
      procedure :: link_transitions
    end type atom_type

    logical, public :: use_classic_col = .False., m1d_legacy_mode=.False.
    logical, public :: lse_no_pivot = .True., kramer_bf=.False.
    logical, public :: disable_bf=.False., disable_bb=.False.
    logical, public :: level_partf=.true.
    real(8), public :: red_thresh_aa  = 20000
    real(8), public :: blue_thresh_aa = 500
    real(8), public :: cont_step_aa = 50
    real(8), public :: logf_cut = -99.9, red_cut=huge(1.0), blue_cut=tiny(1.0)
    real(8), parameter :: ame=5.485799e-4, amh=1.0
    real(8), parameter :: mu_e=ame/(ame+amh)
    real(8), public :: fmin =1e-6
    character(len=32), public :: lse_solver='lineq_no_pivot'
    integer, public :: nq_min=0, nq_max=huge(1)
    real(4), public :: nq_fac=1.0

    public :: unique_int
  
  contains
  
  !====================================================================================
  SUBROUTINE solve_stat_equilibrium (self, temp, ne, totn, n, nstar, nh1, gamma_ij, gamma_ji, cswitch, check_NaN, k)
    class(atom_type), target :: self
    class(transition_type), pointer :: trans
    real(8), intent(in) :: temp, ne, totn, nh1
    real(8), dimension(:), intent(in) :: gamma_ij, gamma_ji
    real(8), dimension(:), intent(in) :: nstar
    real(8), dimension(:), intent(inout) :: n
    integer, intent(in) :: cswitch
    logical, intent(in) :: check_NaN
    integer :: kr, pos, i, j, k
    integer, save :: itimer=0, itimer2=0
    logical :: zero_pop, abort
    logical, save :: first_time=.true.
    !------------------------------------------------------------------------------------
    if (temp > self%tmax_update) then
        self%data%e = 0.0d0
        return
    endif

    self%data%w=0.0d0
    self%data%e=0.0d0
    
    ! Add radiative rates to W
    do kr=1, self%nrad
      if (kr <= self%nline) then
        trans => self%line(kr)
      else
        trans => self%cont(kr-self%nline)
      endif
      if (.not. trans%active) cycle
      i = trans%i
      j = trans%j
      self%data%w(j,i) = self%data%w(j,i) + gamma_ij(kr)
      self%data%w(i,j) = self%data%w(i,j) + gamma_ji(kr)
    enddo

    if (cswitch > 0) then
      call trace%begin('gencol', itimer=itimer)
      self%data%c = 0.d0

      if (use_classic_col) then
        call self%classic_col(temp,ne,nh1,nstar)
      endif

      if (associated(self%coll)) then
        call gencol(self%coll,self%id,self%nlevel,self%g,self%ev,self%ion,&
                    self%data%c,temp,ne,nstar,nh1,self%tgrid,self%cgrid,self%dcdt)
      endif

      call trace%end(itimer=itimer)

      !open(newunit=kr, file='debug.bin', form='unformatted')
      !write(kr) nstar
      !write(kr) n
      !write(kr) totn
      !write(kr) self%data%w
      !write(kr) self%data%c

      do j=1, self%nlevel
      do i=1, self%nlevel
        self%data%w(i,j) = self%data%w(i,j) + self%data%c(i,j) * cswitch
      end do
      end do
    endif

    ! fill diagonal elements of rate coefficient matrix
    do i=1, self%nlevel
      self%data%w(i,i) = -sum(self%data%w(:,i))
    end do

    ! replace row with largest population with particle conservation
    ! equation.
    pos = maxloc(n,1)
    self%data%e(pos) = totn
    !zero_pop = any(n==0)

    !if (zero_pop) then
    !  self%data%oldn = n
    !  self%data%w(pos,:) = 1
    !else
      ! go from rate coefficients to actual rates   
      do i=1,self%nlevel
        self%data%w(i,:) = self%data%w(i,:) * n
      enddo
      self%data%w(pos,:) = n
    !endif

    if (check_NaN .or. first_time) then
      first_time = .false.
      abort = .false.
      do i=1,self%nlevel
        if (self%data%w(i,i) == 0.d0) then
          print *, 'zero in rate diagonal at i = ', i
          abort = .true.
        endif
        do j=1,self%nlevel
          if (self%data%w(i,j) /= self%data%w(i,j)) then
            print *, 'NaN in rate Matrix at i,j = ', i, j
            abort = .true.
          elseif (self%data%w(i,j) > inf8) then
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
      call eqsyst(self%data%w, self%data%e, self%nlevel)
    case('lineq_no_pivot')
      call lineq_no_pivot(self%data%w, self%data%e, self%nlevel)
    case('lineq')
      call lineq(self%data%w, self%data%e, self%nlevel)
    end select

    call trace%end(itimer2)

    if (check_NaN) then
      do i=1,self%nlevel
        if (self%data%e(i) /= self%data%e(i)) then
          print *, 'NaN in error vector at i = ', i
          print *, 'Temp  = ', temp
          stop
        endif
      end do
    endif

    !if (zero_pop) then
    !  self%data%oldn = n
    !  n = self%data%e
    !  do i=1,self%nlevel
    !    if (self%data%oldn(i) > 0) then
    !      self%data%e(i) = n(i) / self%data%oldn(i)
    !    elseif (n(i) == 0) then
    !      self%data%e(i) = 1.d0
    !    else
    !      self%data%e(i) = huge(1.0)
    !    endif
    !  end do
    !  n = max(n, 0.)
    !else
      ! set limits on corrections
      self%data%e = max(self%data%e, 0.001d0)
      !self%data%e = max(self%data%e, tiny(1.0))
      ! update the populations
      n = n * self%data%e
    !endif

    if (check_NaN) then
      do i=1,self%nlevel
        if (n(i) < 0) then
          print *, 'Below zero in populations at i = ', i
          print *, 'Temp  = ', temp
          print *, n(i)
          print *, self%data%e(i)
          !print *, zero_pop, self%data%oldn(i)
          stop
        endif
      end do
    endif

    !write(kr) self%data%e
    !write(kr) n
    !close(kr)
    !stop

    ! set e to the relative change in the populations 
    self%data%e = self%data%e - 1.0d0
  END SUBROUTINE solve_stat_equilibrium
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE check_ionisation_thresholds(self, element)
    class(atom_type), intent(inout) :: self
    class(element_type), pointer :: element
    integer :: nions, j1, j2, i, j
    real(4), dimension(:), allocatable :: xion
    real(4) :: prev
    logical :: first_time=.true.
    !---------------------------------------------------------------------------------
    if (.not. first_time) return
    j1 = minval(self%ion)
    j2 = maxval(self%ion)
    allocate(xion(j1:j2))
    i = 1
    do j=j1, j2
      do i=i, self%nlevel
        if (self%ion(i) == j) then
          xion(j) = self%ev(i)
          exit
        endif
      end do
    end do
    xion(j1:j2-1) = xion(j1+1:j2) - xion(j1:j2-1)
    do j=j1, min(j2-1, element%nion-1)
      if (abs(xion(j)/element%ions(j)%xion - 1) > 0.1) then
        print *, 'WARNING: ionisation threshold for this ion differs from EOS value!'
        print '(I3, 2(F12.4))', j, xion(j), element%ions(j)%xion
      end if
    end do
    first_time = .false.
  END SUBROUTINE check_ionisation_thresholds
  
  !====================================================================================
  SUBROUTINE compute_lte_populations(self, nstar, totn, temp, lsaha, dxi, element, mol)
    class(atom_type), intent(inout) :: self
    class(element_type), pointer :: element
    class(ion_type), pointer :: ion
    class(molec_type), pointer :: mol
    integer :: i, j, x,y,z
    real(8), intent(inout), dimension(:,:,:,:) :: nstar
    real(8), intent(inout), dimension(:,:,:) :: totn
    real(4), intent(in), dimension(:,:,:) :: temp, lsaha, dxi
    real(4), dimension(:,:,:), allocatable :: partition, lvl_part, e_kT
    real(4), dimension(:,:,:), allocatable :: ev_off, ex, norm
    real(4) :: xion
    integer :: ijk(3), idis
    logical :: molec
    integer, save :: itimer=0
    !---------------------------------------------------------------------------------
    call trace%begin('atom_mod%compute_lte_populations', itimer=itimer)
    allocate(partition, lvl_part, norm, mold=temp)
    allocate(e_kT, ev_off, ex, mold=temp)
    e_kT = e_k / temp
    molec = associated(mol)
    if (molec) idis = findloc(self%ion, 2, 1)
  
    nstar = 0
    do i=1, self%nlevel
      if (molec .and. i==idis) cycle
      nstar(:,:,:,i) = self%lte_pop(i, e_kT, dxi, lsaha)
    end do
    lvl_part = sum(nstar, dim=4)
  
    if (.not. level_partf) then
      if (associated(element)) then
        partition = 0
        ev_off = 0
        do j=1, element%nion
          ion => element%ions(j)
          ex = float(j-1)*lsaha - ev_off*e_kT
          where (ex < -80.)
            ex = 0
          elsewhere
            ex = exp(ex)
          endwhere
          partition = partition + ex / ion%part
          ev_off = ev_off + ion%xion - j*dxi
        end do
      elseif (molec) then
        partition = 1 / mol%part
      endif
    elseif (self%id == 'CH') then
      partition = qjorg(temp)
    else
      partition = lvl_part
    endif
    
    norm = totn / partition
    
    do i=1, self%nlevel
      nstar(:,:,:,i) = norm * nstar(:,:,:,i)
    end do
    
    do i=1, self%nlevel
      where(totn == 0) nstar(:,:,:,i) = 0 !ieee_value(nstar(:,:,:,i), ieee_quiet_nan)
      ijk = findloc(nstar(:,:,:,i), tiny(1.d0), mask=totn>0)
      if (any(ijk > 0)) then
        print '(A, I5, A, 3(I4))', ' WARNING: ltepop: level ',i,' has zero lte population at', ijk
        !stop
      end if
    end do
    if (molec) then
      ! replace dissociated state with the ground state population of rarest component
      ! a factor of g0 is missing, which is the statistical weight of the ground state
      ! for carbon this happens to be 1, but for other molecules it may not be
      associate(ndis => nstar(:,:,:,idis))
        ndis = huge(1.e0)
        do j=1, mol%nions
          ion => mol%ions(j)%p
          ndis = min(ndis, ion%n*ion%part)
        end do
      end associate
    endif
  
    ! this is important for number conservation in the SE solver
    totn = sum(nstar, dim=4)
    call trace%end(itimer)
  END SUBROUTINE compute_lte_populations
  !------------------------------------------------------------------------------------

  !====================================================================================
  elemental function lte_pop(self, i, e_kT, dxi, lsaha) result(nstar)
    class(atom_type), intent(in) :: self
    real(4), intent(in) :: e_kT, dxi, lsaha
    integer, intent(in) :: i
    real(8) :: evi, tns
    real(8) :: nstar
    integer :: j, l
    !---------------------------------------------------------------------------------
    j = self%ion(i)
    evi = self%ev(i) - dxi * 0.5*(j-1)*j
    tns = log(self%g(i)) - e_kT*evi

    if(j > self%ion(1)) then
      l = j - self%ion(1)
      tns = tns + float(l)*lsaha
    endif

    nstar = exp(tns)
  end function lte_pop
  !------------------------------------------------------------------------------------
   
  !====================================================================================
  SUBROUTINE compute_cont_factors(self, adata)
    class(atom_type), target, intent(inout) :: self
    class(atom_data_t) :: adata
    real(4) :: bj
    integer :: i, j, k, kr, n(4), x, y, z
    class(cont_type), pointer :: cont
    !---------------------------------------------------------------------------------
    n = [shape(adata%totn), self%ncont]
    call os%allocate(adata%nibj, n, label='data%nibj', zero=.true.)

    do kr=1,self%ncont
      cont => self%cont(kr)
      if (.not. cont%active) cycle
      i = cont%i
      j = cont%j
      do z=1, n(3)
      do y=1, n(2)
      do x=1, n(1)
        if (adata%nstar(x,y,z,i) /= 0) then
          bj = adata%n(x,y,z,j) / adata%nstar(x,y,z,j)
          adata%nibj(x,y,z,kr) = bj * adata%nstar(x,y,z,i)
        endif
      end do
      end do
      end do
    enddo
  END SUBROUTINE compute_cont_factors
  !------------------------------------------------------------------------------------

  !====================================================================================
  SUBROUTINE read_pandas (self, atom_file, qnorm, in_lines)
    class(atom_type), target :: self
    type(line_type), pointer :: lin => null()
    type(cont_type), pointer :: cnt => null()
    type(coll_type), pointer :: col => null()
    character(len=*), intent(in) :: atom_file
    real(8), intent(in) :: qnorm
    integer, intent(in), optional, allocatable :: in_lines(:)
    real(8), pointer, contiguous, dimension(:) :: bf_nu, bf_alpha
    integer, pointer, contiguous, dimension(:) :: icont
    integer :: i, j, ic, tri, idx, latom, nnu_bf, ncol, ntri, ii, jj
    integer :: iostat, nnu, ncont, ntemp, l, u
    real(8) :: scr, gv, cnt_ev, nu_max, alpha0, thresh_aa, dnu
    real(8) :: gij, red_fcut, blu_fcut, vbar2e4, loggf, lam_air
    character(512) :: dummy_string, text
    character(8) :: key, str1, str2, str3
    logical :: abort
    integer, save :: itimer=0
    !---------------------------------------------------------------------------------
    call trace%begin('atom_mod%read_pandas', itimer=itimer)
    call os%check_file(atom_file)
    abort = .false.
    open(newunit=latom, file=atom_file, form='formatted')
    call nextentry(latom)
    read(latom,fmt='(a)') self%id
    self%id=uppercase(adjustl(self%id))

    call nextentry(latom)
    read(latom,*) self%abnd, self%awgt
    call self%prepare_Drawin()

    call nextentry(latom)
    read(latom,*) self%nlevel, self%nline, nnu_bf, ncol

    allocate(self%ev    (self%nlevel))
    allocate(self%g     (self%nlevel))
    allocate(self%label (self%nlevel))
    allocate(self%ion   (self%nlevel))

    ntri = self%nlevel*(self%nlevel-1)/2
    call os%allocate(self%cdata, ntri, zero=.true.)
    call os%allocate(self%ecol , ntri, zero=.true.)
    call os%allocate(self%hcol , ntri, zero=.true.)

    call nextentry(latom)
    read(latom,*) ! skip level header

    do i=1, self%nlevel
      call nextentry(latom)
      read(latom,*) idx, self%ev(i), self%g(i), self%ion(i), self%label(i)
    end do

    if (self%nline > 0) then
      call nextentry(latom)
      read(latom,*) ! skip bbtrans header
    endif

    allocate(self%line(self%nline)) 
    call os%mem_register(storage_size(self%line), self%nline, label='atom%line')

    do i=1, self%nline
      lin => self%line(i)
      call nextentry(latom)
      if (size(in_lines) > 0) then
        if (findloc(in_lines, i, 1) <= 0) then
          call lin%deactivate()
          read(latom,*) ! skip line
          cycle
        endif
      endif

      read(latom, '(A)') text
      read(text, *, iostat=iostat) idx, lin%i, lin%j, lam_air, lin%lower_ev, loggf, &
            lin%ga, lin%gv, lin%gq, lin%nnu, lin%qmax, lin%q0, lin%profile_type, lin%iso

      if (iostat > 0) then
        write(stdout,*) 'Error reading following line '
        write(stdout,*) text
        abort = .true.
      endif

      if (disable_bb) call lin%deactivate()
            
      lin%i = lin%i + 1
      lin%j = lin%j + 1

      lin%lam0 = air2vac(lam_air) * AA_TO_CM ! [cm]
      lin%nu0  = CC / lin%lam0               ! [Hz]
      lin%upper_ev = lin%lower_ev + H_E*lin%nu0
      lin%hydrogen = uppercase(trim(self%id)) == 'H'

      lin%f = 10**loggf / self%g(lin%i)
      gij = self%g(lin%i) / self%g(lin%j)
      call lin%init_einstein_coefficients(lin%f, gij)

      ! make sure we always have an odd number of frequency points for a line
      if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1
      lin%nnu = nq_fac * (lin%nnu-1) + 1
      if (lin%io <= 1 .and. m1d_legacy_mode) then
        lin%nnu = 2 * (lin%nnu-1) + 1
      endif

      if (nq_min < 0) then
        lin%nnu = max(lin%nnu, int(abs(nq_min) + 10*(loggf+5)))
        if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1
      else
        lin%nnu = max(lin%nnu, nq_min)
      endif
      lin%nnu = min(lin%nnu, nq_max)

      if (lin%iso == 0) then
        tri = (lin%j-1)*(lin%j-2)/2 + lin%i
        self%cdata(tri) = lin%f
      endif

      if (logf_cut > -99) then
        scr = log10(self%g(lin%i)*lin%f)
        if (scr < logf_cut) then
          call lin%deactivate()
          cycle
        endif
      endif

      call os%allocate(lin%nu_input, lin%nnu, 'line%nu_input')
      call lin%get_nu(qnorm, lin%nu_input, legacy=m1d_legacy_mode)
      lin%nu_min = lin%nu_input(1)
      lin%nu_max = lin%nu_input(lin%nnu)

      if (red_cut < huge(1.0) .or. blue_cut > tiny(1.0)) then
        red_fcut = cc*1e8 / red_cut
        blu_fcut = cc*1e8 / blue_cut
        if (lin%nu_max < red_fcut .or. lin%nu_min > blu_fcut) then
          call lin%deactivate()
          cycle
        endif
      endif

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
    end do

    ncont = self%nlevel
    allocate(self%cont(ncont))
    call os%mem_register(storage_size(self%cont), ncont, label='atom%cont')

    do i=1, ncont
      cnt => self%cont(i)
      cnt%bf_type='Kramer'
      ! find continuum level
      ic = findloc(self%ion, self%ion(i)+1, 1)
      cnt%i = i
      cnt%j = ic
      cnt%ntrans = self%nline + i ! transition identifier for continua
      if (ic == 0) call cnt%deactivate()
    end do

    call os%allocate(icont, nnu_bf)
    call os%allocate(bf_nu, nnu_bf)
    call os%allocate(bf_alpha, nnu_bf)

    if (nnu_bf > 0) then
      call nextentry(latom)
      read(latom,*) ! skip bftrans header

      do i=1, nnu_bf
        call nextentry(latom) 
        read(latom,*) icont(i), bf_nu(i), bf_alpha(i)
        icont(i) = icont(i)+1
      end do

      ! and convert to frequency [Hz]
      bf_nu = cc / ( bf_nu * AA_TO_CM )

      !  bound-free transitions in detail
      do i=1, ncont
        cnt => self%cont(i)
        ii = findloc(icont, i, 1)
        jj = findloc(icont, i, 1, back=.true.)

        if (jj-ii > 0) then
          cnt%nnu = jj-ii+1
          cnt%nu_input => bf_nu(ii:jj)
          cnt%alpha_input => bf_alpha(ii:jj)
          cnt%alpha0 = cnt%alpha_input(1)
          ! cross sections given in detail
          cnt%bf_type='detail'
        endif
      end do
    endif

    self%ncont = ncont
    self%nrad = self%nline + self%ncont

    j = 0
    do i=1, self%ncont
      cnt => self%cont(i)
      if (disable_bf) call cnt%deactivate()
      if (.not. cnt%active) cycle
      cnt%lower_ev = self%ev(cnt%i)
      cnt%upper_ev = self%ev(cnt%j)
      cnt%nu0 = ( cnt%upper_ev - cnt%lower_ev ) / h_e ! [Hz]

      if (cnt%bf_type == 'Kramer') then
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
      self%cdata(tri) = self%cdata(tri) + cnt%alpha0

      cnt%nu_min = cnt%nu_input(1)
      cnt%nu_max = cnt%nu_input(cnt%nnu)
    end do

    if (j>0) then
      print *, j, " bound-free transitions auto-completed using Kramer's formula!"
    endif

    if (ncol > 0) then
      read(latom, '(A)') dummy_string
      call count_words(dummy_string, ntemp)
      ntemp = ntemp - 3

      allocate( self%coll(ncol) )
      call os%mem_register(storage_size(self%coll), ncol, label='atom%coll')
      call os%allocate(self%tgrid, [ntemp, 1],  label='atom%tgrid')
      call os%allocate(self%cgrid, [ntemp, ncol],    label='atom%cgrid')
      call os%allocate(self%dcdt,  [ntemp, ncol],    label='atom%dcdt')

      read(dummy_string, *), str1, str2, str3, self%tgrid(:, 1)

      do i=1, ncol
        call nextentry(latom) 
        read(latom,*) idx, key, l, u, self%cgrid(1:ntemp,i)
        l = l+1
        u = u+1
        
        if (l==u) then
          write(stdout, *) "ll==ul for process ", key, l, u
          abort = .true.
        endif

        ! within level range
        if(l < 0 .or. u > self%nlevel) then
          write(stdout,*) 'key=', key
          write(stdout,*) 'll=',  l 
          write(stdout,*) 'ul=',  u 
          write(stdout,*) 'nk=',  self%nlevel
          write(stdout,*) 'read_atom:  level index outside range (1,nk)'
          abort = .true.
        endif

        select case(trim(key))
        case('OHM','OHMEGA','CE','CH','CI','FIXRAT','CH0','CH_CE','CHI')
          col => self%coll(i)
          col%key = key
          col%ll = l
          col%ul = u

          select case(trim(key))
          case('OHM','OHMEGA','CE','CI','CH0','CH_CE')
            self%cgrid(1:ntemp,i) = log(self%cgrid(1:ntemp,i))
          end select
          
          ! set derivatives for 3rd order hermite interpolation
          call hermite_derivative(self%dcdt(1:ntemp,i), self%tgrid(1:ntemp,1), self%cgrid(1:ntemp,i))
          col%nc = ntemp
          col%itemp = 1
          tri = (col%ul-1)*(col%ul-2)/2 + col%ll
          !self%cdata(tri) = -1.0 ! flag non-classical recipe

          select case(trim(key))
          case ('CE', 'CI')
            if (self%ecol(tri)) then
              print '(A, I4, I4, A4)', 'Collisional rate already set for transition ', col%ll, col%ul, trim(key)
            endif
            self%ecol(tri) = .true.
          case ('CH', 'CHI')
            if (self%hcol(tri)) then
              print '(A, I4, I4, A4)', 'Collisional rate already set for transition ', col%ll, col%ul, trim(key)
            endif
            self%hcol(tri) = .true.
          end select
        case('FOSC')
          self%coll(i)%nc = 0
          tri = (u-1)*(u-2)/2 + l
          self%cdata(tri) = self%cgrid(1,i)
        end select
      end do
    endif
    close(latom)
    if (abort) stop
    call trace%end(itimer)
  END SUBROUTINE read_pandas

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
    logical, dimension(:,:), allocatable :: ij_check
    real(8), allocatable, dimension(:) :: nu, alpha, tgrid
    real(8) :: scr, gv, cnt_ev, nu_max, alpha0, thresh_aa, dnu
    real(8) :: gij, red_fcut, blu_fcut, vbar2e4, loggf
    integer :: ncont, nnu
    character(len=20) :: txt
    character(len=26),parameter :: upcase_alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    logical :: abort
    integer, save :: itimer=0
    !---------------------------------------------------------------------
    call trace%begin('atom_mod%read', itimer=itimer)
    call os%check_file(atom_file)
    open(newunit=latom, file=atom_file, form='formatted')
    abort = .false.

    call nextentry(latom)
    read(latom,fmt='(a)') self%id
    self%id=uppercase(adjustl(self%id))

    call nextentry(latom)
    read(latom,*,iostat=ios) self%abnd, self%awgt
    call print_failed_line(latom, ios)
    call self%prepare_Drawin()

    call nextentry(latom)
    read(latom,*,iostat=ios) self%nlevel, self%nline, self%ncont, self%nrfix
    call print_failed_line(latom, ios)

    allocate(  self%ev    ( self%nlevel) ) ; self%ev = 0.0d0
    allocate(  self%g     ( self%nlevel) ) ; self%g = 0
    allocate(  self%label ( self%nlevel) ) ; self%label = ' '
    allocate(  self%ion   ( self%nlevel) ) ; self%ion = 0

    call os%allocate(self%cdata, self%nlevel*(self%nlevel-1)/2, zero=.true.)
    call os%allocate(self%ecol , self%nlevel*(self%nlevel-1)/2, zero=.true.)
    call os%allocate(self%hcol , self%nlevel*(self%nlevel-1)/2, zero=.true.)

    do i=1, self%nlevel
      call nextentry(latom)
      read(latom,*,iostat=ios) self%ev(i),self%g(i), self%label(i), self%ion(i)
      self%ev(i)=self%ev(i)*hc_e  ! [eV]  
      call print_failed_line(latom, ios)
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

    allocate(ij_check(self%nlevel, self%nlevel), source=.false.)
    allocate(self%line(self%nline)) 
    call os%mem_register(storage_size(self%line), self%nline, label='atom%line')

    do i=1,self%nline
      lin => self%line(i)
      call nextentry(latom) 

      ! io WAS present for historical reasons, in MULTI io determines
      ! whether a line is wide or not. In the RH formulation all lines are wide

      read(latom,*,iostat=ios) jj, ii, lin%f , lin%nnu, lin%qmax, lin%q0, lin%io, &
            lin%ga, lin%gv, lin%gq, lin%profile_type
      call print_failed_line(latom, ios)

      lin%i = min(ii, jj)
      lin%j = max(ii, jj)
      lin%q0   = lin%q0  
      lin%qmax = lin%qmax

      if (ij_check(lin%i, lin%j)) then
        print *, 'Transition assigned twice!', lin%i, lin%j
      else
        ij_check(lin%i, lin%j) = .true.
      endif

      if (lin%j > self%nlevel) then
        print *, 'Transition to non-existing level!', i, lin%i, lin%j
        abort = .true.
      endif

      ! make sure we always have an odd number of frequency points for a line
      if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1

      lin%nnu = nq_fac * (lin%nnu-1) + 1
      if (lin%io <= 1 .and. m1d_legacy_mode) then
        lin%nnu = 2 * (lin%nnu-1) + 1
      endif

      if (nq_min < 0) then
        loggf = log10(lin%f*self%g(lin%i))
        lin%nnu = max(lin%nnu, int(abs(nq_min) + 10*(loggf+5)))
        if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1
      else
        lin%nnu = max(lin%nnu, nq_min)
      endif
      lin%nnu = min(lin%nnu, nq_max)

      tri = (lin%j-1)*(lin%j-2)/2 + lin%i
      self%cdata(tri) = lin%f

      lin%upper_g  = self%g(lin%j)
      lin%lower_ev = self%ev(lin%i)
      lin%upper_ev = self%ev(lin%j)
      lin%hydrogen = uppercase(trim(self%id)) == 'H'

      lin%lam0 = hc_e / ( lin%upper_ev - lin%lower_ev ) ! [cm]
      lin%nu0 = cc / lin%lam0    ! [Hz]

      ! io is now modified so that multicomponent lines can be read
      ! io = 0: Normal line
      ! io = 1: For now is the same as io = 0
      ! io > 1: Number of components in term
      ! Like in Multi, thereafter follows one extra input line per component:

      !   dq  wt  fterm  da  gwterm  gqterm
      ! dq    energy difference for component i of transition kr. unit cm-1.
      ! wt     weight of component kterm in weighted average
      ! fterm     oscillator strength
      ! da    damping constants
      ! gwterm
      ! gqterm
      ! adterm    adamp for component, calculated in damp

      nterm = lin%io
      if (nterm .ge. 2) then
        write(stdout, *) uppercase(rname), &
        ' found ' ,nterm, 'components in transition KR = ',i

        allocate(lin%dq (nterm) )
        allocate(lin%wt  (nterm) )
        allocate(lin%fterm  (nterm) )
        allocate(lin%da (nterm) )
        allocate(lin%gwterm (nterm) )
        allocate(lin%gqterm (nterm) )
        !allocate(lin%adterm, xs, xe, ys, ye, zs1, ze1, 1, nterm, rname)

        do iterm = 1, nterm
          read(latom,iostat=ios) lin%dq(iterm), lin%wt(iterm),  lin%fterm(iterm), &
                        lin%da(iterm), lin%gwterm(iterm), lin%gqterm(iterm)
          call print_failed_line(latom, ios)
        end do

        lin%dq = cc/lin%nu0 * lin%dq
        lin%da = (lin%da - lin%ga) * over_4pi / lin%nu0
      end if

      ! in gfortran the size of an unallocated allocatable array can be larger than 0
      if (allocated(in_lines) .and. size(in_lines) > 0) then
        if (findloc(in_lines, i, 1) <= 0) then
          call lin%deactivate()
          cycle
        endif
      endif

      if (logf_cut > -99) then
        scr = log10(self%g(lin%i)*lin%f)
        if (scr < logf_cut) then
          call lin%deactivate()
          cycle
        endif
      endif

      call os%allocate(lin%nu_input, lin%nnu, 'line%nu_input')
      call lin%get_nu(qnorm, lin%nu_input, legacy=nterm<=1 .and. m1d_legacy_mode)
      lin%nu_min = lin%nu_input(1)
      lin%nu_max = lin%nu_input(lin%nnu)

      gij = self%g(lin%i) / self%g(lin%j)
      call lin%init_einstein_coefficients(lin%f, gij)

      if (red_cut < huge(1.0) .or. blue_cut > tiny(1.0)) then
        red_fcut = cc*1e8 / red_cut
        blu_fcut = cc*1e8 / blue_cut
        if (lin%nu_max < red_fcut .or. lin%nu_min > blu_fcut) then
          call lin%deactivate()
          cycle
        endif
      endif

      lin%ntrans = i
      lin%profile_type = trim(lowercase(lin%profile_type))

      select case (lin%profile_type)
      case ('voigt', 'hlinop')
        continue
      case default
        if (index(lin%profile_type, '.')>0) then
          continue
        else
          write (*, fmt = '(a,1x,a,a)') 'Profile not recognised:', &
                  lin%profile_type, ' Setting it to VOIGT'
        endif
        lin%profile_type = 'voigt'
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
      read(latom,*,iostat=ios) jj, ii, alpha0, nnu, nu_max
      call print_failed_line(latom, ios)

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
          read(latom,*,iostat=ios)  cnt%nu_input(j), cnt%alpha_input(j)
          call print_failed_line(latom, ios)
        end do
        
        do j=2,cnt%nnu
          if ( cnt%nu_input(j) > cnt%nu_input(j-1) ) then
            write(stdout,*) 'read_atom: continuum transition nr', i, &
                          ' wavelengths not decreasing!'
            abort = .true.
          elseif ( cnt%nu_input(j) == cnt%nu_input(j-1) ) then
            write(stdout,*) 'read_atom: continuum transition nr', i, &
                          ' wavelengths are not unique!'
            cnt%nu_input(j:nnu) = cnt%nu_input(j:nnu) - 1.0d-6
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
      if (disable_bf) then
        call cnt%deactivate()
        cycle
      endif
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
      key = trim(uppercase(adjustl(key)))

      ! read next non-comment line if key doesnt start with a letter
      if (scan(key(1:1),upcase_alphabet)==0) cycle 
      

      ! we found a keyword, so interpret it
      select case(key)
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

      key = trim(uppercase(adjustl(key)))

      if (key == 'END') then
        !write(*,*) "read_collisions: read all collisions."
        exit
      else
        call nextentry(latom)
      endif
      
      select case(key)
      case('TEMP')
        itemp = itemp + 1
        read(latom,*,iostat=ios) ntemp, tgrid(1:ntemp)
        call print_failed_line(latom, ios)
        self%tgrid(1:ntemp, itemp) = tgrid(1:ntemp)
        if (ntemp < 2) then
            stop "read_collisions: ntemp < 2"
        endif

      ! MB  5.04.2018 - added CH0 keyword
      ! RH 25.03.2022 - added CHI keyword
      case('OHM','OHMEGA','CE','CH','CI','FIXRAT','CH0','CH_CE','CHI')
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
        read(latom,*,iostat=ios) l, u, self%cgrid(1:nc,i)
        call print_failed_line(latom, ios)
        
        if (l==u) then
          write(stdout, *) "ll==ul for process ", key, l, u
          abort = .true.
        endif

        if (any(self%cgrid(1:nc,i) <= 0)) then
          write(stdout, *) "Non-positive cross sections for process ", key, l, u
          abort = .true.
        else
          select case(trim(key))
          case('OHM','OHMEGA','CE','CI','CH0','CH_CE')
            self%cgrid(1:ntemp,i) = log(self%cgrid(1:ntemp,i))
          end select
        endif

        select case(trim(key))
        case('OHM','OHMEGA','CE','CH','CI','FIXRAT','CH0','CH_CE','CHI')
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
          if (self%ecol(tri)) then
            print '(A, I4, I4, A4)', 'Collisional rate already set for transition ', col%ll, col%ul, trim(key)
          endif
          self%ecol(tri) = .true.
        case ('CH', 'CHI')
          if (self%hcol(tri)) then
            print '(A, I4, I4, A4)', 'Collisional rate already set for transition ', col%ll, col%ul, trim(key)
          endif
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
    call trace%end(itimer)
  END SUBROUTINE read

  !====================================================================================
  SUBROUTINE prepare_damping(self, species)
    class(atom_type), intent(inout) :: self
    class(massive_type), intent(in), pointer :: species
    class(ion_type), pointer :: ions(:)
    class(isotope_type), pointer :: isotopes(:)
    class(line_type), pointer :: line
    logical :: is_element
    integer :: i, j, ic, iso, tri
    real(8) :: vbar2e4, cnt_ev, gij
    !---------------------------------------------------------------------
    select type(species)
    class is (element_type)
      is_element = .true.
      ions => species%ions
      isotopes => species%isotope
    class default
      is_element = .false.
    end select

    do i=1, self%nline
      line => self%line(i)
      if (.not. line%active) cycle
      j = self%ion(line%j)

      if (is_element .and. line%iso > 0) then
        iso = findloc(isotopes%iso, line%iso, 1)
        if (iso > 0) then
          line%f = line%f * isotopes(iso)%frac
          gij = self%g(line%i) / self%g(line%j)
          call line%init_einstein_coefficients(line%f, gij)
          tri = (line%j-1)*(line%j-2)/2 + line%i
          self%cdata(tri) = self%cdata(tri) + line%f
        else
          print *, 'isotope', line%iso, ' not found. Deactivating line ', i
          call line%deactivate()
          cycle
        endif
      endif

      if (line%hydrogen) then
        call prepare_hydrogen_damping(line)
      else
        call line%prepare_Stark(line%gq)
        vbar2e4 = self%compute_vbar2(10000.)
        if (associated(species)) then
          if (is_element) then
            cnt_ev = ions(j)%xion
          else
            cnt_ev = species%xion
          endif
        else
          ! find continuum level
          ic = findloc(self%ion, j+1, 1)
          cnt_ev = merge(self%ev(ic), 999.99, ic > 0)
        endif
        call line%prepare_VanderWaals(line%gv, self%ion(line%i), cnt_ev, vbar2e4)
      endif
    end do
  END SUBROUTINE prepare_damping

  !====================================================================================
  SUBROUTINE prepare_Drawin(self)
    class(atom_type), intent(inout) :: self
    real(8) :: drmu, mu_h, ama
    real(8), parameter :: const=16*PI*Rbohr**2 ! 1.40758E-15
    !---------------------------------------------------------------------
    ama = self%awgt
    mu_h = amh/(ama+amh)
    drmu = PI*MH*ama*mu_h
    self%drawin = const * sqrt(2.0*BK/drmu) * mu_e/mu_h
  END SUBROUTINE prepare_Drawin

  !====================================================================================
  SUBROUTINE classic_col(self,temp,ne,nh1,nstar)
    class(atom_type), intent(inout) :: self
    real(8), intent(in) :: temp, ne, nh1
    real(8), dimension(:), intent(in) :: nstar
    real(8) :: ni, nj, nij, cup, cdn, e_kT, divu
    real(8) :: dxi, dev, ex, u0, f, tmp, sqt, fac, nu0
    real(8) :: seaton, allen, regemorter
    integer :: i,j, ii, jj, tri, n
    logical :: cont
    logical, save :: first_time=.true.
    real(8), parameter :: sq32 = sqrt(3.d0) / (2.d0)
    real(8), parameter :: sqr3_2pi = sq32 / PI
    real(8), parameter :: sq8bmp = sqrt(8*BK/EM/PI)
    real(8), parameter :: c0 = PI*Rbohr**2 * sq8bmp
    real(8), parameter :: cseaton = sq8bmp/sq32/(1.d0/137.d0)*Rydk !1.55e13
    !---------------------------------------------------------------------
    if (debye_shield) then
      dxi = 4.98e-4*5040./temp*sqrt(ne*bk*temp)
    else
      dxi = 0
    endif
    sqt = sqrt(temp)
    ! See Hubeny % Mihalas 2014, p. 276 for these recipes
    allen = c0 * Rydk * ne / sqt
    regemorter = 14.5 * allen * Rydk / e_k
    seaton = 1.55e13 * ne / sqt
    e_kT = e_k / temp

    n = 0
    do i=1, self%nlevel
      ni = nstar(i)
      do j=i+1, self%nlevel
        tri = (j-1)*(j-2)/2 + i

         ! inverted bound-free transitions are unclear at the moment
        if (self%ion(j) < self%ion(i)) cycle

        cont = self%ion(j) > self%ion(i)
        if (.not. self%ecol(tri) .or. .not. self%hcol(tri)) then
          fac = self%cdata(tri)
          nj = nstar(j)
          nij = ni / nj

          dev = self%ev(j) - self%ev(i)
          dev = max(dev, 1d-8) ! to avoid division by zero

          if (cont) then
            jj = self%ion(j) * (self%ion(j)-1)
            ii = self%ion(i) * (self%ion(i)-1)
            dev = dev - dxi*0.5*(jj-ii)
          else
            f = max(fac, fmin*self%g(j))
          endif

          u0 = dev * e_kT
          ex = exp(-u0)

          ! collisions with electrons
          if (.not. self%ecol(tri)) then
            if (cont) then
              ! ionization by electron collisions
              if (fac == 0) then
                nu0 = dev * EE/HH
                fac = Kramer(nu0, dev)
              endif
              tmp = min(0.1*self%ion(i), 0.3)
              cup = seaton * fac * tmp / u0 * ex
              cup = cup * cbfesc
            else
              ! excitation of neutrals by electron collisions 
              !if (fac > 0) then
                ! Van Regemorter
                tmp = max(0.2*ex, sqr3_2pi*expin(u0))
                cup = regemorter * f * tmp / dev
              !else
                ! Allen
                ! This recipe is flawed as it does not work with superlevels
                !cup = allen * ex / self%g(i)
              !endif
              cup = cup * cbbesc
            endif
            cdn = cup * nij
            self%data%c(j,i) = self%data%c(j,i) + cup ! upward rate
            self%data%c(i,j) = self%data%c(i,j) + cdn ! downward rate
            n = n + 1
          endif

          ! collisions with hydrogen (Drawin)
          if (.not. self%hcol(tri)) then
            divu = (RydV/dev)**2
            tmp = 2.0/u0
            tmp = (1.0 + tmp) / (1.0 + (mu_e*tmp)**2)
            cup = nh1 * self%drawin * divu * tmp * ex * sqt
            
            if (cont) then
              cup = cup * cbfhsc
            else
              cup = cup * f
              cup = cup * cbbhsc
            endif
            cdn = cup * nij
            self%data%c(j,i) = self%data%c(j,i) + cup ! upward rate
            self%data%c(i,j) = self%data%c(i,j) + cdn ! downward rate
            n = n + 1
          endif
        endif
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
    if (option(bf) .and. associated(self%cont)) then
      do i=1, size(self%cont)
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        ! find minimum and maximum index of this continuum in the sp%nu grid
        call sp%get_interval(cnt, ired, iblue)
        if (iblue == 0 .or. ired >= sp%nnu) then
          call cnt%deactivate()
        elseif (iblue==ired) then
          print *, 'WARNING, single frequency for cont', i
          call cnt%deactivate()
        else
          ! total number of frequencies where this continuum is active
          nnu=iblue-ired+1
          cnt%nnu = nnu
          ! min and max frequency where this continuum is active
          cnt%ired  = ired
          cnt%iblue = iblue

          ! set new frequency grid
          allocate(cnt%alpha(ired:iblue))
          call os%mem_register(storage_size(cnt%alpha), cnt%nnu)
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

          allocate(cnt%wnu(ired:iblue))
          call os%mem_register(storage_size(cnt%wnu), cnt%nnu, 'cont%wnu')
          call set_weights(cnt%wnu, cnt%nu)
          cnt%wnu = cnt%wnu / (h_4pi*cnt%nu)
        endif
      end do
    endif

    if (option(bb) .and. associated(self%line)) then
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
          !call os%allocate(lin%q, lin%nnu)
          allocate(lin%q(ired:iblue))
          call os%mem_register(storage_size(lin%q), nnu)
          lin%q(ired:iblue) = lin%nu/lin%nu0 - 1.0d0

          ! fill frequency integration weights for the line
          allocate(lin%wnu(ired:iblue))
          call os%mem_register(storage_size(lin%wnu), lin%nnu, 'line%wnu')
          call set_weights(lin%wnu, lin%nu)
          lin%wnu = lin%wnu / (h_4pi*lin%nu)

          allocate(lin%wq(ired:iblue))
          call os%mem_register(storage_size(lin%wq), nnu)
          call set_weights(lin%wq, dble(lin%q))
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

    !call self%get_active_nlevels(sp)

    write(*,'(a,I0,a)') 'sp_from_atom: frequency grid has ', sp%nnu, ' points.'
    if (self%verbose > 2) then
      write(*,'(a,I0,a)') 'sp_from_atom: ', sum(sp%nac+sp%nal), ' active transition-frequency points'
    endif
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
  