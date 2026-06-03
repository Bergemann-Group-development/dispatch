! *******************************************************************
MODULE atom_mod
    USE io_mod
    USE handyfuncs_mod
    USE m3d_debug_mod
    USE prec_mod, only: inf8
    USE const_mod
    USE trans_mod
    USE op_math_mod, only: lineq, get_interval, eqsyst, lineq_no_pivot
    USE level_mod, only: state_type
    USE ion_mod, only: ion_type, massive_type, debye_shield, massive_type, particle_type
    USE element_mod, only: element_type, isotope_type
    USE molecule_mod, only: molec_type
    USE composition_mod, only: composition_type
    USE os_mod
    USE trace_mod
    USE spectrum_mod, only: spectrum_type, set_weights
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

    !---------------------------------------------------------------------
    type, public, extends(spectrum_type) :: atomic_spectrum_type
      ! maximum number of active lines or continua at any given frequency
      integer :: maxal=0, maxac=0, nmax=0, nlower, nupper
      ! which continua (ac) and lines (al) are active at each frequency
      integer, dimension(:,:), contiguous, pointer :: ac,al
      ! how many continua (nac) and lines (nal) are active at each frequency
      integer, dimension(:), contiguous, pointer :: nac,nal
      ! how many different atoms are active at each frequency
      integer, dimension(:), contiguous, pointer :: na
      ! has this frequency an active prd line
      !logical,dimension(:),allocatable :: hasPRDline
    end type atomic_spectrum_type
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    type, public :: atom_data_t
      ! LTE, NLTE and total population of active atom
      real(8), dimension(:,:,:),   allocatable :: totn
      real(8), dimension(:,:,:,:), allocatable :: nstar
      real(8), dimension(:,:,:,:), allocatable :: b
      real(8), dimension(:,:,:,:,:), allocatable :: oldb
      real(8), dimension(:,:,:,:), allocatable :: gamma_ij
      real(8), dimension(:,:,:,:), allocatable :: gamma_ji
      real(4), dimension(:,:,:,:), allocatable :: wphi
      real(4), dimension(:,:,:),   allocatable :: inv_df
      real(8), dimension(:),   allocatable  :: e
      real(8), dimension(:,:), allocatable  :: C, W
      real(8), dimension(:), allocatable :: lvlerr
      real(8) :: maxerr
      logical :: converged=.false.
      integer :: elvl, epos(3)
    end type atom_data_t
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    type, public :: atom_trans_t
      type(line_type), dimension(:), pointer, contiguous :: line => null()
      type(cont_type), dimension(:), pointer, contiguous :: cont => null()
      real(4), dimension(:,:,:,:), allocatable :: adamp
      real(4), dimension(:,:,:,:), allocatable :: gamma_ij
      real(4), dimension(:,:,:,:), allocatable :: gamma_ji
    contains
      procedure :: allocate => allocate_trans
      procedure :: set_ired_iblue
      procedure :: set_active_transitions
      procedure :: append_transitions
    end type atom_trans_t
    !---------------------------------------------------------------------
  
    !---------------------------------------------------------------------
    type, public, extends(massive_type) :: atom_type
      integer :: nrad   ! total number of radiative transitions
      integer :: nrfix  ! number of fixed transitions
      integer :: ncont  ! number of bound-free transitions
      integer :: nline  ! number of bound-bound transitions
      integer :: nlevel ! number of levels
      integer :: nprd=0 ! number of PRD transitions
      integer :: verbose=0
      integer :: idep=0 ! index of atom in departure coefficients array
      logical :: precomputed = .false. ! whether departure coeffs are precomputed
      
      character(len=128) :: name
      character(len=3)   :: id    ! atom id string 
      character(len=20)  :: crout ! which collisional SUBROUTINE to use 
  
      character(len=20), dimension(:), allocatable :: label ! level label
  
      integer, dimension(:), allocatable :: ion ! ionization stage of each level
  
      real(8) :: abnd ! abundance ( in the [hydrogen]=12 scale)
      !real(8) :: awgt ! atomic weight [g]
  
      real(8), dimension(:), allocatable :: ev ! energy of level above ground state [eV]
      real(8), dimension(:), allocatable :: g ! multiplicity of level
      integer, dimension(:), allocatable :: irad, jrad ! indices of levels in radiative transition

      type(atom_trans_t), pointer :: trans => null()
      type(atom_data_t), pointer :: data => null()
      type(atom_data_t), dimension(:), pointer :: rotated_data => null()
      type(coll_type), dimension(:), pointer, contiguous :: coll => null()
      class(massive_type), pointer :: species => null()

      real(8), dimension(:,:), pointer, contiguous :: tgrid => null()
      real(8), dimension(:,:), pointer, contiguous :: cgrid => null()
      real(8), dimension(:,:), pointer, contiguous :: dcdt  => null()
      real(4), dimension(:), pointer, contiguous :: cdata => null()
      logical, dimension(:), pointer, contiguous :: ecol => null()
      logical, dimension(:), pointer, contiguous :: hcol => null()
      
      real(8) :: qnorm
      real(8) :: drawin, mu_e
      logical :: check=.false.

      real(8), pointer :: ng_a => null()
      real(8), pointer :: ng_b => null()
      real(8), dimension(:,:), pointer, contiguous :: ng_d => null()
        
    contains
      procedure :: init_atom_data
      procedure :: check_ionisation_thresholds
      procedure :: lte_pop
      procedure :: compute_lte_populations
      procedure :: normalize_profiles
      procedure :: update_populations
      procedure :: setup_transition_matrix
      procedure :: read
      procedure :: read_pandas
      procedure :: prepare_damping
      procedure :: prepare_Drawin
      procedure :: classic_col
      procedure :: gencol
      procedure :: set_colliders
      procedure :: ng_acceleration
      procedure :: ng_acceleration_accumulate
      procedure :: ng_acceleration_coeff
      procedure :: apply_ng_acceleration
    end type atom_type
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    type coll_type
      character(len=16) :: key ! transition type 
      integer :: nc ! size of coefficient grid
      integer :: ll,ul ! lower/upper level for transition
      integer :: itemp
      integer :: dir=0
      logical :: log = .false.
      integer :: np=0, ion=0, mol=0 ! collider identification
    end type coll_type
    !---------------------------------------------------------------------

    character(len=8), dimension(7), parameter :: theprocesses = &
    (/'OHMEGA  ','CI      ','CH      ','CE      ', &
      'OHM     ','CH0     ','CHI     '/)

    real(8), public :: cbbesc=1., cbfesc=1., cbbhsc=1., cbfhsc=1., colscale=1.
    logical, public :: use_classic_col = .False., m1d_legacy_mode=.False.
    logical, public :: lse_no_pivot = .True., kramer_bf=.False.
    logical, public :: disable_bf=.False., disable_bb=.False.
    logical, public :: level_partf=.false.
    logical, public :: ch0_hmip=.false.
    real(8), public :: red_thresh_aa  = 20000
    real(8), public :: blue_thresh_aa = 500
    real(8), public :: cont_step_aa = 50
    real(8), public :: logf_cut = -99.9, red_cut=huge(1.0), blue_cut=tiny(1.0)
    real(8), parameter :: ame=5.485799e-4, amh=1.0
    real(8), parameter :: mu_e=ame/(ame+amh)
    real(8), public :: fmin=1e-6
    real(8) :: tmax_update=1.0d7
    character(len=32), public :: lse_solver='lineq_no_pivot'
    integer, public :: nq_min=0, nq_max=huge(1)
    real(4), public :: nq_fac=1.0
    integer, public :: cswitch=1
    integer, public :: ng_step=-1, ng_start=8
    logical, public :: ng_global=.true.
    logical, public :: save_rate_matrix=.false.

    public :: unique_int
  
  contains

  !====================================================================================
  SUBROUTINE allocate_trans(self, nb, nf, bb, bf)
    class(atom_trans_t), target :: self
    integer, intent(in), optional :: nb, nf
    type(line_type), intent(in), dimension(:), optional :: bb
    type(cont_type), intent(in), dimension(:), optional :: bf
    !---------------------------------------------------------------------------------
    if (present(nb)) then
      allocate(self%line(nb))
    endif
    if (present(nf)) then
      allocate(self%cont(nf))
    endif
    if (present(bb)) then
      allocate(self%line(size(bb)))
      self%line = bb
    endif
    if (present(bf)) then
      allocate(self%cont(size(bf)))
      self%cont = bf
    endif
    call os%mem_register(storage_size(self%line), size(self%line), label='atom%line')
    call os%mem_register(storage_size(self%cont), size(self%cont), label='atom%cont')
  END SUBROUTINE allocate_trans

  !===============================================================================
  SUBROUTINE init_atom_data (self, comp, adata, iterate, found)
    class(atom_type), target :: self
    class(composition_type), target :: comp
    class(atom_data_t), target :: adata
    class(element_type), pointer :: element
    class(molec_type), pointer :: molec
    logical, intent(in) :: iterate
    logical, intent(out), optional :: found
    real(4) :: nmax, nlim
    integer, save:: itimer=0
    integer :: i, j, m(4), nl
    !-----------------------------------------------------------------------------
    call trace%begin('atom_type%init_atom_data', itimer=itimer)
    nl = self%nlevel
    m = [shape(comp%temp), nl]
    call os%allocate(adata%nstar, m, label='nstar')
    call os%allocate(adata%totn, m(1:3), label='totn')
    call os%allocate(adata%inv_df, m(1:3), label='inv_df')
    
    call comp%get_totn(self%id, self%abnd, adata%totn, element, molec)
    if (present(found)) found = associated(element) .or. associated(molec)
    
    if (associated(molec)) tmax_update = comp%tmolim
    if (associated(element)) call self%check_ionisation_thresholds(element)
    call self%compute_lte_populations(adata%nstar, adata%totn, real(comp%temp), real(comp%lsaha), real(comp%dxi), element, molec)
    adata%inv_df = self%inverse_doppler_factor(real(comp%temp), comp%vmic)
    
    if (iterate) then
      call os%allocate(adata%C, [nl, nl])
      call os%allocate(adata%W, [nl, nl])
      call os%allocate(adata%e, nl)
      call os%allocate(adata%lvlerr, nl)
      call os%allocate(adata%wphi, [m(1:3), self%nline], label='wphi', zero=.true.)
      call os%allocate(adata%b, m, label='b')
      do i=1, self%nlevel
        adata%b(:,:,:,i) = 1.d0
      end do
      if (ng_step > 0) then
        call os%allocate(adata%oldb, [m, 3], label='oldb')
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE init_atom_data

  !===============================================================================
  SUBROUTINE normalize_profiles(self, iter)
    class(atom_type), target :: self
    class(atom_data_t), pointer :: adata
    class(transition_type), pointer :: line, cont
    integer, intent(in) :: iter
    real(8), dimension(:,:,:,:), pointer, contiguous :: gamma_ij, gamma_ji, nstar
    real(4), dimension(:,:,:,:), pointer, contiguous :: wphi
    real(8) :: gij
    integer :: i
    ! ----------------------------------------------------------------------------
    adata => self%data
    wphi => adata%wphi
    nstar => adata%nstar
    gamma_ij => adata%gamma_ij
    gamma_ji => adata%gamma_ji

    if (iter == 1) then
      ! Normalize wphi
      do i=1, self%nline
        if (.not. self%trans%line(i)%active) cycle
        wphi(:,:,:,i) = 1 / wphi(:,:,:,i)
      end do
    endif

    do i=1, self%nline
      line => self%trans%line(i)
      if (.not. line%active) cycle
      gij = self%g(line%i) / self%g(line%j)
      gamma_ij(:,:,:,i) = gamma_ij(:,:,:,i) * wphi(:,:,:,i)
      gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * wphi(:,:,:,i)
      gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * gij
    end do
    
    do i=i, self%nrad
      cont => self%trans%cont(i-self%nline)
      if (.not. cont%active) cycle
      gamma_ji(:,:,:,i) = gamma_ji(:,:,:,i) * nstar(:,:,:,cont%i) / nstar(:,:,:,cont%j)
    end do

    if (allocated(adata%oldb)) then
      do i=1, self%nlevel
        adata%oldb(:,:,:,i,1) = adata%oldb(:,:,:,i,2)
        adata%oldb(:,:,:,i,2) = adata%oldb(:,:,:,i,3)
        adata%oldb(:,:,:,i,3) = adata%b(:,:,:,i)
      end do
    endif
  END SUBROUTINE normalize_profiles

  !===============================================================================
  SUBROUTINE update_populations (self, comp, local, pcheck, cij, cji, w, c)
    class(atom_type), target :: self
    class(composition_type) :: comp
    class(atom_data_t), pointer :: adata
    real(8), intent(out) :: local
    integer, intent(in) :: pcheck(3)
    real(4), intent(out), optional, dimension(:,:,:) :: cij, cji
    real(8), intent(out), optional, dimension(:,:) :: w, c
    integer :: i, m(3), x, y, z, nk
    real(8) :: norm
    real(8), parameter :: edamp = 1e-3
    logical :: relax, ok
    integer, save :: itimer=0, itimer2=0
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_type%update_populations', itimer=itimer)
    adata => self%data
    m = shape(adata%totn)
    adata%maxerr = 0
    adata%lvlerr = 0
    nk = self%nlevel
    
    do z=1, m(3)
      do y=1, m(2)
        do x=1, m(1)
          if (adata%totn(x,y,z) > 0 .and. comp%temp(x,y,z) < tmax_update) then
            call self%setup_transition_matrix(comp, x,y,z)

            ! e contains the correction n/n_old 
            associate(e=>adata%e(:), b=>adata%b(x,y,z,:), nstar=>adata%nstar(x,y,z,:), &
                      maxerr=>adata%maxerr, lvlerr=>adata%lvlerr)

              if (self%check .and. save_rate_matrix) then
                if (all([x, y, z] == pcheck)) then
                  forall(i=1:nk) w(i,:) = adata%w(i,:) / (b*nstar)
                  c = adata%c
                end if
              endif

              call trace%begin('lse_solver', itimer=itimer2)
              ! solve w n . x = e
              ! e contains the solution vector x on exit
              select case (lse_solver)
              case('eqsyst')
                call eqsyst(adata%w, e, nk)
              case('lineq_no_pivot')
                call lineq_no_pivot(adata%w, e, nk)
              case('lineq')
                call lineq(adata%w, e, nk)
              end select
              call trace%end(itimer2)
              
              if (check_NaN) then
                call check_array(e, 'lineq solution', check_zero=.true., check_negative=.true., abort=.false., ok=ok)
                if (.not.ok) then
                  print *, 'ERROR: in lineq solution at x,y,z = ', x, y, z
                  stop
                endif
              endif
              
              ! compute the exponential update factor
              where (e < edamp)
                e = exp(max(e, 0.0d0)/edamp - 1d0 + log(edamp))
              endwhere

              ! apply the update
              b = b * e

              ! compute the normalisation scalar
              norm = adata%totn(x,y,z) / sum(b * nstar)

              ! enforce strict macroscopic particle conservation
              b = b * norm
              
              if (check_NaN) then
                call check_array(e, 'population correction', check_zero=.true., abort=.false., ok=ok)
                if (.not.ok) then
                  print *, 'ERROR: in population correction vector at x,y,z = ', x, y, z
                  stop
                endif
              endif

              if (check_NaN) then
                call check_array(b, 'departure coefficients', check_zero=.true., abort=.false., ok=ok)
                if (.not.ok) then
                  print *, 'ERROR: in departure coefficients at x,y,z = ', x, y, z
                  stop
                endif
              endif

              ! compute the true relative change in the populations
              e = (e * norm) - 1.0d0

              ! keep track of maximum population change
              lvlerr = merge(e, lvlerr, abs(e) > abs(lvlerr))
              i = maxloc(abs(lvlerr), dim=1)

              if (abs(lvlerr(i)) > abs(maxerr)) then
                maxerr = lvlerr(i)
                adata%epos = [x, y, z]
                adata%elvl = i
              endif

              if (self%check) then
                if (all([x, y, z] == pcheck)) then
                  local = e(lvlcheck)
                end if
                if (present(cij)) then
                  cij(x,y,z) = adata%c(self%jrad(linecheck), self%irad(linecheck))
                  cji(x,y,z) = adata%c(self%irad(linecheck), self%jrad(linecheck))
                endif
              endif
            end associate
          endif
        end do
      end do
    end do

    call trace%end(itimer)
  END SUBROUTINE update_populations
  
  !====================================================================================
  SUBROUTINE setup_transition_matrix (self, comp, x,y,z)
    class(atom_type), target :: self
    class(composition_type) :: comp
    class(transition_type), pointer :: trans
    real(8), dimension(:), pointer :: b, nstar, rij, rji
    real(8) :: totn
    integer :: kr, pos, i, j, x,y,z
    integer, save :: itimer=0, itimer2=0
    logical :: abort, ok
    logical, save :: first_time=.true.
    !------------------------------------------------------------------------------------
    totn = self%data%totn(x,y,z)
    b => self%data%b(x,y,z,:)
    nstar => self%data%nstar(x,y,z,:)
    rij => self%data%gamma_ij(x,y,z,:)
    rji => self%data%gamma_ji(x,y,z,:)

    self%data%w=0.0d0
    self%data%e=0.0d0

    ! Add radiative rates to W
    do kr=1, self%nrad
      if (kr <= self%nline) then
        trans => self%trans%line(kr)
      else
        trans => self%trans%cont(kr-self%nline)
      endif
      if (.not. trans%active) cycle
      i = trans%i
      j = trans%j
      self%data%w(j,i) = self%data%w(j,i) + rij(kr)
      self%data%w(i,j) = self%data%w(i,j) + rji(kr)
    enddo

    ! one can switch off the collisional rates by setting cswitch to 0
    if (cswitch > 0) then
      self%data%c = 0.d0

      if (use_classic_col) then
        call self%classic_col(comp,x,y,z)
      endif

      if (associated(self%coll)) then
        call self%gencol(comp,x,y,z)
      endif

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

    ! replace row with largest population with particle conservation equation.
    pos = maxloc(b * nstar, dim=1)
    self%data%e(pos) = totn

    ! go from rate coefficients to actual rates   
    do i=1, self%nlevel
      self%data%w(i,:) = self%data%w(i,:) * b * nstar
    enddo
    self%data%w(pos,:) = b * nstar

    if (check_NaN .or. first_time) then
      first_time = .false.
      call check_array([(self%data%w(i,i), i=1, self%nlevel)], 'rate diagonal', check_zero=.true., abort=.false., ok=ok)
      call check_array(self%data%w, 'rate matrix', check_zero=.false., abort=.false., ok=ok)
      if (.not. ok) then
        print *, 'ERROR: in rate matrix at x,y,z = ', x, y, z
        stop
      endif
    endif
  END SUBROUTINE setup_transition_matrix
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
          ex = float(j-1)*lsaha - ev_off*e_kT + ion%part
          ex = exp(ex)
          partition = partition + ex
          ev_off = ev_off + ion%xion - j*dxi
        end do
      elseif (molec) then
        partition = exp(-mol%part)
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

    if (molec) then
      ! replace dissociated state with the ground state population of rarest component
      ! a factor of g0 is missing, which is the statistical weight of the ground state
      ! for carbon this happens to be 1, but for other molecules it may not be
      associate(ndis => nstar(:,:,:,idis))
        ndis = huge(1.e0)
        partition = exp(-mol%part)
        do j=1, mol%nions
          ion => mol%ions(j)%p
          ndis = min(ndis, ion%n * partition)
        end do
      end associate
    endif
    
    do i=1, self%nlevel
      where(totn == 0) nstar(:,:,:,i) = 0 !ieee_value(nstar(:,:,:,i), ieee_quiet_nan)
      if (check_NaN) then
        call check_array(nstar(:,:,:,i), 'LTE populations', check_zero=.true., abort=.true.)
      endif
    end do
  
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
    real(8) :: gij, red_fcut, blu_fcut, vbar2e4, loggf
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

    ncont = self%nlevel
    allocate(self%trans)
    call self%trans%allocate(self%nline, ncont)

    do i=1, self%nline
      lin => self%trans%line(i)
      call nextentry(latom)
      if (size(in_lines) > 0) then
        if (findloc(in_lines, i, 1) <= 0) then
          call lin%deactivate()
          read(latom,*) ! skip line
          cycle
        endif
      endif

      read(latom, '(A)') text
      read(text, *, iostat=iostat) idx, lin%i, lin%j, lin%lam_air, lin%lower_ev, loggf, &
            lin%ga_in, lin%gv_in, lin%gq_in, lin%nnu, lin%qmax, lin%q0, lin%profile_type, lin%iso

      if (iostat > 0) then
        write(stdout,*) iostat, ' Error reading following line '
        write(stdout,*) text
        abort = .true.
      endif

      if (disable_bb) call lin%deactivate()
            
      lin%i = lin%i + 1
      lin%j = lin%j + 1

      lin%lam0 = air2vac(lin%lam_air) * AA_TO_CM ! [cm]
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
      elseif (lin%nnu > 0) then
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

      call lin%get_nu(qnorm, lin%nu_input, legacy=m1d_legacy_mode)
      call os%mem_register(storage_size(lin%nu_input), lin%nnu, 'line%nu_input')
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

    do i=1, ncont
      cnt => self%trans%cont(i)
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

      !  bound-free transitions in detail
      do i=1, ncont
        cnt => self%trans%cont(i)
        ii = findloc(icont, i, 1)
        jj = findloc(icont, i, 1, back=.true.)

        if (jj-ii > 0) then
          cnt%nnu = jj-ii+1
          cnt%nu_input => bf_nu(ii:jj)
          cnt%alpha_input => bf_alpha(ii:jj)
          cnt%alpha0 = cnt%alpha_input(1)

          do j=2,cnt%nnu
            if ( cnt%nu_input(j) > cnt%nu_input(j-1) ) then
              write(stdout,*) 'read_atom: continuum transition nr', i, &
                            ' wavelengths not decreasing!'
              abort = .true.
            elseif ( cnt%nu_input(j) == cnt%nu_input(j-1) ) then
              write(stdout,*) 'read_atom: continuum transition nr', i, &
                            ' wavelengths are not unique!'
              cnt%nu_input(j:cnt%nnu) = cnt%nu_input(j:cnt%nnu) - 1.0d-6
            endif
          end do
          
          ! cross sections given in detail
          cnt%bf_type='detail'
        endif
      end do

      ! and convert to frequency [Hz]
      bf_nu = cc / ( bf_nu * AA_TO_CM )
    endif

    self%ncont = ncont
    self%nrad = self%nline + self%ncont

    j = 0
    do i=1, self%ncont
      cnt => self%trans%cont(i)
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

    allocate(self%irad(self%nrad), self%jrad(self%nrad))
    self%irad(1:self%nline) = self%trans%line%i
    self%jrad(1:self%nline) = self%trans%line%j
    self%irad(self%nline+1:self%nrad) = self%trans%cont%i
    self%jrad(self%nline+1:self%nrad) = self%trans%cont%j

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
        read(latom, *, iostat=iostat) idx, key, l, u, self%cgrid(1:ntemp,i)
        call print_failed_line(latom, iostat)
        l = l+1
        u = u+1
        if (u > l) then
          tri = (u-1)*(u-2)/2 + l
        else
          tri = (l-1)*(l-2)/2 + u
        endif
        
        if (l==u) then
          write(stdout, '(A, A, I4, I4, A)') "ll==ul for process ", trim(key), l, u, ". Ignoring this one."
          cycle
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

        col => self%coll(i)
        col%key = key

        select case(trim(key))
        case('FOSC')
          col%nc = 0
          self%cdata(tri) = self%cgrid(1,i)
        case default
          col%ll = l
          col%ul = u

          self%cgrid(1:ntemp,i) = log(self%cgrid(1:ntemp,i))
          col%log = .true.
          
          ! set derivatives for 3rd order hermite interpolation
          call hermite_derivative(self%dcdt(1:ntemp,i), self%tgrid(1:ntemp,1), self%cgrid(1:ntemp,i))
          col%nc = ntemp
          col%itemp = 1
          !self%cdata(tri) = -1.0 ! flag non-classical recipe

          select case(trim(key))
          case ('CE', 'CI')
            if (self%ecol(tri)) then
              print '(A, I4, I4, A4)', 'WARNING: Collisional rate already set for transition ', col%ll, col%ul, trim(key)
            endif
            self%ecol(tri) = .true.
          case ('CH', 'CHI')
            if (self%hcol(tri)) then
              print '(A, I4, I4, A4)', 'WARNING: Collisional rate already set for transition ', col%ll, col%ul, trim(key)
            endif
            self%hcol(tri) = .true.
          end select
        end select
      end do
    endif
    close(latom)
    if (abort) stop 'abort in atom_mod%read_pandas'
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
    character(512) :: text
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

    ncont = count(self%ion < maxval(self%ion))
    if (ncont < self%ncont) then
        print *, 'WARNING: Unrealisticly many bound-free transitions. Something is wrong!'
    endif

    allocate(ij_check(self%nlevel, self%nlevel), source=.false.)

    allocate(self%trans)
    call self%trans%allocate(self%nline, ncont)

    do i=1,self%nline
      lin => self%trans%line(i)
      call nextentry(latom) 

      ! io WAS present for historical reasons, in MULTI io determines
      ! whether a line is wide or not. In the RH formulation all lines are wide

      read(latom,*,iostat=ios) jj, ii, lin%f , lin%nnu, lin%qmax, lin%q0, lin%io, &
            lin%ga_in, lin%gv_in, lin%gq_in, lin%profile_type
      call print_failed_line(latom, ios)

      lin%i = min(ii, jj)
      lin%j = max(ii, jj)
      lin%q0   = lin%q0  
      lin%qmax = lin%qmax

      if (ij_check(lin%i, lin%j)) then
        print '(A, I4, I4)', 'WARNING: Transition assigned twice!', lin%i, lin%j
      else
        ij_check(lin%i, lin%j) = .true.
      endif

      if (lin%j > self%nlevel) then
        print '(A, I8, I4, I4)', 'ERROR: Transition to non-existing level!', i, lin%i, lin%j
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
      elseif (lin%nnu > 0) then
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
          read(latom, '(A)') text
          read(text, *, iostat=ios) lin%dq(iterm), lin%wt(iterm),  lin%fterm(iterm), &
                        lin%da(iterm), lin%gwterm(iterm), lin%gqterm(iterm)
          call print_failed_line(latom, ios)
        end do

        lin%dq = cc/lin%nu0 * lin%dq
        lin%da = (lin%da - lin%ga_in) * over_4pi / lin%nu0
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

      call lin%get_nu(qnorm, lin%nu_input, legacy=nterm<=1 .and. m1d_legacy_mode)
      call os%mem_register(storage_size(lin%nu_input), lin%nnu, 'line%nu_input')
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

    do i=1, ncont
      cnt => self%trans%cont(i)
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
        print '(A, I4, A, I4, I4)', 'ERROR: Bound-free transition', i, ' is between two levels of the same ion', ii, jj
        abort = .true.
      endif

      j  = max(ii,jj)
      ii = min(ii,jj)

      cnt => self%trans%cont(ii)

      if (cnt%bf_type /= 'Kramer') then
        print '(A, I4, A)', 'WARNING: Multiple bound-free transitions from lower state', ii, '. Ignoring the latter!'
        do j=1, nnu
          call nextentry(latom) 
          read(latom, '(A)') ! skipping
        end do
        cycle
      endif
      if (j /= cnt%j) then
        print '(A, I4, A)', 'WARNING: Bound-free transition', i, ' has exited state as upper level. Worth checking!'
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
      cnt => self%trans%cont(i)
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

    allocate(self%irad(self%nrad), self%jrad(self%nrad))
    self%irad(1:self%nline) = self%trans%line%i
    self%jrad(1:self%nline) = self%trans%line%j
    self%irad(self%nline+1:self%nrad) = self%trans%cont%i
    self%jrad(self%nline+1:self%nrad) = self%trans%cont%j

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
      case('OHM','OHMEGA','CE','CH','CI','CH0','CHI')
        nc = ntemp

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
          write(stdout, '(A, A, I4, I4, A)') "ll==ul for process ", trim(key), l, u, ". Ignoring this one."
          cycle
        endif

        if (any(self%cgrid(1:nc,i) <= 0)) then
          write(stdout, *) "Non-positive cross sections for process ", key, l, u
          abort = .true.
        else
          select case(trim(key))
          case('OHM','OHMEGA','CE','CI','CH0')
            self%cgrid(1:ntemp,i) = log(self%cgrid(1:ntemp,i))
            col%log = .true.
          end select
        endif

        select case(trim(key))
        case('OHM','OHMEGA','CE','CH','CI','CH0','CHI')
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
            print '(A, I4, I4, A4)', 'WARNING: Collisional rate already set for transition ', col%ll, col%ul, trim(key)
          endif
          self%ecol(tri) = .true.
        case ('CH', 'CHI')
          if (self%hcol(tri)) then
            print '(A, I4, I4, A4)', 'WARNING: Collisional rate already set for transition ', col%ll, col%ul, trim(key)
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
    if (abort) stop 'abort in atom_mod%read'
    call trace%end(itimer)
  END SUBROUTINE read

  !====================================================================================
  SUBROUTINE prepare_damping(self, species)
    class(atom_type), intent(inout) :: self
    class(massive_type), intent(in), pointer :: species
    class(ion_type), pointer :: ions(:)
    class(isotope_type), pointer :: isotopes(:)
    class(line_type), pointer :: line
    logical :: is_molecule
    integer :: i, j, ic, iso, tri
    real(8) :: vbar2e4, cnt_ev, gij
    !---------------------------------------------------------------------
    nullify(ions)
    select type(species)
    class is (element_type)
      ions => species%ions
      isotopes => species%isotope
    class is (molec_type)
      isotopes => species%isotope
    end select
    is_molecule = .not. associated(ions)

    do i=1, self%nline
      line => self%trans%line(i)
      if (.not. line%active) cycle
      j = self%ion(line%j)

      if (line%iso > 0) then
        iso = findloc(isotopes%iso, line%iso, 1)
        if (iso > 0) then
          line%f = line%f * isotopes(iso)%frac
          gij = self%g(line%i) / self%g(line%j)
          call line%init_einstein_coefficients(line%f, gij)
          tri = (line%j-1)*(line%j-2)/2 + line%i
          self%cdata(tri) = self%cdata(tri) + line%f
        else
          print *, 'WARNING: isotope', line%iso, ' not found. Deactivating line ', i
          call line%deactivate()
          cycle
        endif
      endif

      if (line%hydrogen) then
        call prepare_hydrogen_damping(line)
      else
        call line%prepare_Radiative(line%ga_in)
        vbar2e4 = self%compute_vbar2(10000.)
        ! the following commented section is generally more correct, but not backwards compatible.
        ! mul23 uses the ionisation threshold from the atom, not from the EOS.
        ! I don't know why many of our atoms, e.g. Mn, have wrong ionisation thresholds...
        !if (associated(species)) then
        !  if (.not. is_molecule) then
        !    cnt_ev = ions(j)%xion
        !  else
        !    cnt_ev = species%xion ! will not be used
        !  endif
        !else
        ! find continuum level
        ic = findloc(self%ion, j+1, 1)
        cnt_ev = merge(self%ev(ic), 999.99d0, ic > 0)
        !endif
        call line%prepare_Stark(line%gq_in, self%ion(line%i), cnt_ev)
        call line%prepare_VanderWaals(line%gv_in, self%ion(line%i), cnt_ev, vbar2e4, is_mol=is_molecule)
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
  SUBROUTINE classic_col(self, comp, x, y, z)
    class(atom_type), intent(inout) :: self
    class(composition_type) :: comp
    real(8) :: temp, ne, nh1
    real(8) :: ni, nj, nij, cup, cdn, e_kT, divu
    real(8) :: dxi, dev, ex, u0, f, tmp, sqt, fac, nu0
    real(8) :: seaton, allen, regemorter
    integer :: x, y, z, i, j, ii, jj, tri, n
    logical :: cont
    logical, save :: first_time=.true.
    real(8), parameter :: sq32 = sqrt(3.d0) / (2.d0)
    real(8), parameter :: sqr3_2pi = sq32 / PI
    real(8), parameter :: sq8bmp = sqrt(8*BK/EM/PI)
    real(8), parameter :: c0 = PI*Rbohr**2 * sq8bmp
    real(8), parameter :: cseaton = sq8bmp/sq32/(1.d0/137.d0)*Rydk !1.55e13
    integer, save :: itimer=0
    !---------------------------------------------------------------------
    call trace%begin('classic_col', itimer=itimer)
    dxi = comp%dxi(x,y,z)
    temp = comp%temp(x,y,z)
    ne = comp%ne(x,y,z)
    e_kT = comp%e_kT(x,y,z)
    nh1 = comp%elements(1)%ions(1)%levels(1)%n(x,y,z)
    sqt = sqrt(temp)
    ! See Hubeny % Mihalas 2014, p. 276 for these recipes
    allen = c0 * Rydk * ne / sqt
    regemorter = 14.5 * allen * Rydk / e_k
    seaton = 1.55e13 * ne / sqt

    n = 0
    do i=1, self%nlevel
      ni = self%data%nstar(x,y,z,i)
      do j=i+1, self%nlevel
        tri = (j-1)*(j-2)/2 + i

         ! inverted bound-free transitions are unclear at the moment
        if (self%ion(j) < self%ion(i)) cycle

        cont = self%ion(j) > self%ion(i)
        if (.not. self%ecol(tri) .or. .not. self%hcol(tri)) then
          fac = self%cdata(tri)
          nj = self%data%nstar(x,y,z,j)
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
            if (cup /= 0) then
              cup = cup * colscale
              cdn = cup * nij
              self%data%c(j,i) = self%data%c(j,i) + cup ! upward rate
              self%data%c(i,j) = self%data%c(i,j) + cdn ! downward rate
              n = n + 1
            endif
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
            if (cup /= 0) then
              cup = cup * colscale
              cdn = cup * nij
              self%data%c(j,i) = self%data%c(j,i) + cup ! upward rate
              self%data%c(i,j) = self%data%c(i,j) + cdn ! downward rate
              n = n + 1
            endif
          endif
        endif
      end do
    end do
    if (first_time .and. n>0) then
      first_time = .false.
      write(stdout, *) n, ' collisional transitions filled with classical recipes.'
    endif
    call trace%end(itimer)
  END SUBROUTINE classic_col

  !===============================================================================
  ! initialises the spectrum structure
  !===============================================================================
  SUBROUTINE append_transitions (self, sp, bb, bf, sort)
    class(atom_trans_t), target :: self
    class(spectrum_type) :: sp
    logical, intent(in), optional :: bb, bf, sort
    type(line_type), pointer :: lin
    type(cont_type), pointer :: cnt
    integer :: i,nnu,js,je
    real(8), dimension(:), allocatable :: nu
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_trans_t%append_transitions', itimer=itimer)
    ! compute total number of frequencies
    nnu = sp%nnu
    if (option(bb)) nnu = nnu + sum(self%line%nnu, self%line%active)
    if (option(bf)) nnu = nnu + sum(self%cont%nnu, self%cont%active)

    ! save all frequencies
    allocate(nu(nnu))
    if (allocated(sp%nu)) nu(1:sp%nnu) = sp%nu

    je=sp%nnu
    if (option(bb)) then
      do i=1, size(self%line)
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
      do i=1, size(self%cont)
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        js=je+1
        je=js+cnt%nnu-1
        nu(js:je)=cnt%nu_input
        cnt%ired = js
        cnt%iblue = je
      end do
    endif

    ! now set overall frequency grid
    if (allocated(sp%nu)) deallocate(sp%nu)
    allocate(sp%nu, source=nu)
    sp%nnu = size(nu)

    call trace%end(itimer)
  END SUBROUTINE append_transitions

  !===============================================================================
  SUBROUTINE set_ired_iblue (self, sp, bb, bf, set_active)
    class(atom_trans_t), target :: self
    class(atomic_spectrum_type), target :: sp
    logical, intent(in), optional :: bb, bf, set_active
    type(line_type), pointer :: lin
    type(cont_type), pointer :: cnt
    integer :: i, j, nnu, ired, iblue
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_trans_t%set_ired_iblue', itimer=itimer)
    if (option(bf) .and. associated(self%cont)) then
      do i=1, size(self%cont)
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        ! find minimum and maximum index of this continuum in the sp%nu grid
        call sp%get_interval(cnt, ired, iblue)
        if (iblue == 0 .or. ired >= sp%nnu) then
          call cnt%deactivate()
        elseif (iblue==ired) then
          print '(A, I6)', 'WARNING: single frequency for cont', i
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
          if (iblue==ired) print '(A, I8)', 'WARNING: single frequency for line', i
          ! total number of frequencies where this line is active
          nnu=iblue-ired+1
          lin%nnu = nnu
          ! min and max frequency where this line is active
          lin%ired  = ired
          lin%iblue = iblue

          lin%nu => sp%nu(ired:iblue)
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
    class(atom_trans_t), target :: self
    class(atomic_spectrum_type), target :: sp
    logical, intent(in), optional :: bb, bf
    type(line_type), pointer :: lin
    type(cont_type), pointer :: cnt
    integer :: i, j, nnu, ired, iblue
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_trans_t%set_active_transitions', itimer=itimer)

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

    call trace%end(itimer)
  END SUBROUTINE set_active_transitions

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

  !====================================================================================
  SUBROUTINE set_colliders(self, comp)
    class(atom_type) :: self
    class(composition_type) :: comp
    class(particle_type), pointer :: collider
    integer :: i, s, l
    integer, save :: itimer=0
    !---------------------------------------------------------------------
    call trace%begin('atom_type%set_colliders', itimer=itimer)
    do i=1, size(self%coll)
      associate(col => self%coll(i))
          select case (trim(col%key))
          case ('CH','CHI','CH0','OHM','OHMEGA','CE','CI','FOSC')
            ! no collider needed
          case default
            l = len(col%key)
            s = index(lowercase(col%key), 'up')
            if (s > 0) then
                s = s + 2
                col%dir = +1
            else
                s = index(lowercase(col%key), 'dn')
                if (s > 0) then
                  s = s + 2
                  col%dir = -1
                end if
            endif
            call comp%particle_by_id(trim(col%key(s:l)), collider, verbose=.true., &
                                     np=col%np, ion=col%ion, mol=col%mol)
          end select
      end associate
    end do
    call trace%end(itimer)
  END SUBROUTINE set_colliders

  !====================================================================================
  SUBROUTINE gencol(self, comp, x, y, z)
    class(atom_type) :: self
    class(composition_type) :: comp
    character(8) :: key
    real(8) :: temp, ne
    real(8) :: ct, cup, cdn, dekt, de_k, de
    real(8) :: stg, nji, fac, nh1, n_collider
    real(8) :: lsaha, dxi, dev, tns, e_kT
    integer :: x,y,z,ul,ll,i,j,k,iz, nhlevels, ii, jj, nc, itemp
    real(8), parameter :: c0=((Eryd/sqrt(em)) * pi*Rbohr**2) * sqrt(8.0/(pi*bk))
    integer, save :: itimer=0
    !---------------------------------------------------------------------
    call trace%begin('gencol', itimer=itimer)
    ! get local copies for ease of notation
    temp = comp%temp(x,y,z)
    ne = comp%ne(x,y,z)
    dxi = comp%dxi(x,y,z)
    e_kT = comp%e_kT(x,y,z)
    lsaha = comp%lsaha(x,y,z)
    stg=sqrt(temp)
    nh1 = comp%elements(1)%ions(1)%levels(1)%n(x,y,z)

    ! loop over all the transitions
    do i=1, size(self%coll, 1)
      associate(col => self%coll(i))
        if (col%nc == 0) cycle
        cup=0.0d0
        cdn=0.0d0
        ll=col%ll
        ul=col%ul

        dev = self%ev(ul) - self%ev(ll)
        if (self%ion(ul) > self%ion(ll)) then
          ii = self%ion(ll) * (self%ion(ll)-1)
          jj = self%ion(ul) * (self%ion(ul)-1)
          dev = dev - dxi*0.5*(jj-ii)
        endif

        if (self%data%nstar(x,y,z,ll) > 0) then
          nji = self%data%nstar(x,y,z,ul) / self%data%nstar(x,y,z,ll)
        else
          tns = - dev*e_kT + (self%ion(ul) - self%ion(ll))*lsaha
          nji = self%g(ul)/self%g(ll) * exp(tns)
        endif

        nc = col%nc
        if (nc > 0) then
          itemp = col%itemp
          call hermite_interpolation(ct,temp,self%tgrid(1:nc,itemp),self%cgrid(1:nc,i),self%dcdt(1:nc,i))

          if (col%log) then
              ct = exp(ct)
          endif
        endif

        if (col%np /= 0) then
          n_collider = comp%elements(col%np)%ions(col%ion)%n(x,y,z)
        elseif (col%mol /= 0) then
          n_collider = comp%molecules(col%mol)%n(x,y,z)
        endif

        key = trim(adjustl(col%key))
        select case (key)
        case('OHM','OHMEGA')
          !  excitation of electron collisons with positive ions
          cdn = c0 / self%g(ul) * ct * ne / stg
          cup = cdn * nji

        case('CE')
          ! excitation of neutrals by electron collisions 
          cdn = self%g(ll)/self%g(ul) * ct * ne * stg
          cdn = cdn * cbbesc
          cup = cdn * nji

        case('CH')
          ! b-b collisions with neutral hydrogen
          ! these are downward rates (cdn); cup are upward rates
          cdn = ct * nh1
          cdn = cdn * cbbhsc
          cup = cdn * nji

        case('CHI')
          ! b-f collisions with neutral hydrogen
          ! these are upwards rates (cup); cdn are downward rates
          cup = ct * nh1
          cup = cup * cbfhsc
          cdn = cup / nji

        case('CI')
          ! ionization by electron collisions
          cup = ct * ne * stg * exp(-(self%ev(ul) - self%ev(ll))*e_kT)
          cup = cup * cbfesc
          cdn = cup / nji

        case('CH0')
          ! MB 5.04.2018 - charge transfer with H-, A+H <--> A+ + H-
          ! these are upwards rates (CT ionisation or ion-pair formation)
          ! also assumed for CI; but the difference is that CI downwards rates 
          ! are computed from detailed balance (they are bi-directional);
          ! CH0 are one-directional and the levels must be in the correct order
          !  in the atom (, i.e. ll  --> lu for ioninisation)
          ! in the QM papers, the inverse exothermic rates (neutralisation) 
          ! are derived from the detailed balance
          ! N_kj(T,E_k) = N_jk(t,E_k) * exp(-E_jk/ (k_B * T))
          cup = ct * nh1
          if (ch0_hmip) then
            cup = cup * exp(-0.754 * e_kT)
          endif
          cup = cup * cbfhsc
          cdn = cup / nji

        case default
          if (col%dir == -1) then ! downwards transition
            cdn = ct * n_collider
            cup = cdn * nji

          elseif (col%dir == 1) then ! upwards transition
            cup = ct * n_collider
            cup = cdn / nji
          end if
        end select

        !********************************************************************

        cup = cup * colscale
        cdn = cdn * colscale
        
        self%data%c(ul,ll) = self%data%c(ul,ll) + cup ! upward rate Cll->ul 
        self%data%c(ll,ul) = self%data%c(ll,ul) + cdn ! downward rate Cul->ll
      end associate
    end do ! transition loop 
    call trace%end(itimer)
  END SUBROUTINE gencol

  !====================================================================================
  logical function is_valid_process(key)
    character(len=8),intent(in) :: key
    integer :: i
    !---------------------------------------------------------------------
    is_valid_process=.false.
    do i=1,size(theprocesses)
        if ( key(1:len_trim(key)) == trim(theprocesses(i)))  then
          is_valid_process=.true.
          exit
        end if
    end do
  end function is_valid_process

  !====================================================================================
  ELEMENTAL FUNCTION expin (x) result(ex1)
    real(8), intent(in) :: x
    real(8) :: ex, ex1
    real(8), parameter, dimension(6) :: &
    a = [-44178.5471728217,57721.7247139444,9938.31388962037,    &
          1842.11088668000,101.093806161906,5.03416184097568]
    real(8), parameter, dimension(5) :: &
    b = [76537.3323337614,32597.1881290275,6106.10794245759,     &
          635.419418378382,37.2298352833327]
    real(8), parameter, dimension(7) :: &
    c = [ 4.65627107975096e-7,.999979577051595,9.04161556946329, &
        24.3784088791317,23.0192559391333,6.90522522784444,     &
          .430967839469389]
    real(8), parameter, dimension(6) :: &
    d = [10.0411643829054,32.4264210695138, 41.2807841891424,    &
        20.4494785013794,3.31909213593302, .103400130404874]
    real(8), parameter, dimension(7) :: &
    e = [-.999999999998447,-26.6271060431811,-241.055827097015,  &
        -895.927957772937,-1298.85688746484,-545.374158883133,  &
          -5.66575206533869]
    real(8), parameter, dimension(6) :: &
    f = [28.6271060422192,292.310039388533,1332.78537748257,     &
        2777.61949509163,2404.01713225909, 631.657483280800]
    !---------------------------------------------------------------------
    ex = exp(-x)
    if (x > 4) then
      ex1 = (ex+ex*(e(1)+(e(2)+(e(3)+(e(4)+(e(5)+(e(6)+e(7)/x)/x)/x)/x)/x)/x)/  &
            (x+f(1)+(f(2)+(f(3)+(f(4)+(f(5)+f(6)/x)/x)/x)/x)/x))/x
    elseif (x > 1) then
      ex1 = ex*(c(7)+(c(6)+(c(5)+(c(4)+(c(3)+(c(2)+c(1)*x)*x)*x)*x)*x)*x)/      &
            (d(6)+(d(5)+(d(4)+(d(3)+(d(2)+(d(1)+x)*x)*x)*x)*x)*x)
    elseif (x > 0) then
      ex1 = (a(1)+(a(2)+(a(3)+(a(4)+(a(5)+a(6)*x)*x)*x)*x)*x)/                  &
            (b(1)+(b(2)+(b(3)+(b(4)+(b(5)+x)*x)*x)*x)*x) - log(x)
    else
      ex1 = 0.d0
    endif
  END FUNCTION expin

  !====================================================================================
  ELEMENTAL FUNCTION Kramer (nu, dev) result(hyd)
    real(8), intent(in) :: nu, dev
    real(8) :: nl, hyd
    real(8), parameter :: kramer_const = 64*PI**4 * EC**10 * EM / (3*sqrt(3.d0) * CC * HH**6) ! 2.815e29
    !---------------------------------------------------------------------
    if (dev > 0) then
      nl = sqrt(RydV/dev)
      hyd = kramer_const * gaunt(nu,ceiling(nl)) / (nl**5 * nu**3)
    else
      hyd = 0.
    endif
  END FUNCTION Kramer

  !====================================================================================
  ELEMENTAL FUNCTION gaunt (nu, n)
    real(8), intent(in) :: nu
    integer, intent(in) :: n
    real(8) :: gaunt
    real(8), parameter, dimension(6) :: &
    a1 = [ 1.078e00, 1.0926e0, 1.0983e0, 1.0954e0, 1.0912e0, 1.0925e0], &
    b1 = [-8.754e14,-2.019e14,-9.45e13, -5.188e13,-3.2e13,  -2.331e13], &
    c1 = [-1.791e29, 1.836e28, 9.177e27, 3.552e27, 1.576e27, 9.325e26], &
    a2 = [ 0.798e00, 0.768e00, 0.793e00, 0.831e00, 0.758e00, 0.79e00],  &
    b2 = [ 5.358e15, 6.242e15, 5.48e15,  4.094e15, 6.633e15, 5.808e15], &
    c2 = [-3.484e31,-3.208e31,-2.318e31,-1.43e31, -3.32e31, -2.844e31]
    !---------------------------------------------------------------------
    if (n > 6) then
      gaunt = 1
    else
      if (nu <= 1e16) then
          gaunt = a1(n) + (b1(n) + c1(n)/nu) / nu
      else
          gaunt = a2(n) + (b2(n) + c2(n)/nu) / nu
      endif
    endif
  END FUNCTION gaunt

  !===============================================================================
  SUBROUTINE ng_acceleration (self, ng_iter, emax, epos)
    class(atom_type):: self
    integer, intent(in) :: ng_iter
    real(8), intent(out) :: emax
    integer, intent(out) :: epos(4)
    integer :: i, x, y, z, k, nvec, m(3)
    real(8), dimension(:), allocatable :: n1, n2, n3, n4, d0, d1, d2, e, w
    real(8) :: d01, d02, d11, d12, d21, d22, a, b, den, enew
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_type%ng_acceleration', itimer=itimer)
    emax = 0
    if (ng_global) then
      if (ng_iter==0) then
        call self%ng_acceleration_accumulate()
      elseif (ng_iter==1) then
        ! the first one to get here will convert d into a and b
        ! and will reset d to zero; the others will reuse a and b
        !$omp critical (ng_global)
        call self%ng_acceleration_coeff()
        !$omp end critical (ng_global)
        call self%apply_ng_acceleration(emax, epos)
      endif
    elseif (ng_iter == 0) then
      allocate(self%ng_d(0:2,1:2), source=0.d0)
      allocate(self%ng_a, self%ng_b)
      call self%ng_acceleration_accumulate()
      call self%ng_acceleration_coeff()
      call self%apply_ng_acceleration(emax, epos)
      deallocate(self%ng_d)
      deallocate(self%ng_a, self%ng_b)
    endif
    call trace%end(itimer)
  END SUBROUTINE ng_acceleration

  !===============================================================================
  SUBROUTINE ng_acceleration_accumulate (self)
    class(atom_type):: self
    integer :: i, x, y, z, k, nvec, m(3)
    real(8), dimension(:), allocatable :: n1, n2, n3, n4, d0, d1, d2, e
    real(8) :: w
    ! ----------------------------------------------------------------------------
    m = shape(self%data%totn)
    nvec = self%nlevel * product(m(1:3))
    allocate(e(nvec))
    allocate(n1, n2, n3, n4, d0, d1, d2, mold=e)
    i = 0
    do z=1, m(3)
    do y=1, m(2)
    do x=1, m(1)
    do k=1, self%nlevel
      i = i+1
      ! Extract natural logarithms; explicit weighting is no longer required
      n1(i) = log(self%data%oldb(x,y,z,k,1))
      n2(i) = log(self%data%oldb(x,y,z,k,2))
      n3(i) = log(self%data%oldb(x,y,z,k,3))
      n4(i) = log(self%data%b(x,y,z,k))
    end do
    end do
    end do
    end do

    d0 = n4 - n3
    d1 = n3 - n2
    d2 = n2 - n1

    d1 = d0 - d1
    d2 = d0 - d2

    associate(d => self%ng_d)
      d(0,1) = d(0,1) + dot_product(d0, d1)
      d(0,2) = d(0,2) + dot_product(d0, d2)
      d(1,1) = d(1,1) + dot_product(d1, d1)
      d(1,2) = d(1,2) + dot_product(d1, d2)
      d(2,1) = d(2,1) + dot_product(d2, d1)
      d(2,2) = d(2,2) + dot_product(d2, d2)
    end associate
  END SUBROUTINE ng_acceleration_accumulate

  !===============================================================================
  SUBROUTINE ng_acceleration_coeff (self)
    class(atom_type):: self
    real(8) :: den
    ! ----------------------------------------------------------------------------
    associate(d => self%ng_d)
      if (any(d /= 0.d0)) then
        den = d(1,1)*d(2,2) - d(1,2)*d(2,1)
        if (abs(den) > 0) then
          self%ng_a = (d(0,1)*d(2,2) - d(0,2)*d(2,1)) / den
          self%ng_b = (d(0,2)*d(1,1) - d(0,1)*d(1,2)) / den
        else
          self%ng_a = 0
          self%ng_b = 0
        endif
        d = 0 ! reset for next ng update
      endif
    end associate
  END SUBROUTINE ng_acceleration_coeff

  !===============================================================================
  SUBROUTINE apply_ng_acceleration (self, emax, epos)
    class(atom_type):: self
    integer :: i, x, y, z, k, epos(4), m(3)
    real(8) :: n1, n2, n3, n4, e
    real(8) :: emax, enew
    ! ----------------------------------------------------------------------------
    m = shape(self%data%totn)
    i = 0
    associate(a => self%ng_a, b => self%ng_b)
      do z=1, m(3)
      do y=1, m(2)
      do x=1, m(1)
        do k=1, self%nlevel
          i = i+1
          n1 = log(self%data%oldb(x,y,z,k,1))
          n2 = log(self%data%oldb(x,y,z,k,2))
          n3 = log(self%data%oldb(x,y,z,k,3))
          n4 = log(self%data%b(x,y,z,k))
          
          ! Extrapolate in logarithmic space
          e = (1.d0 - a - b)*n4 + a*n3 + b*n2
          
          ! Exponentiate to recover strictly positive departure coefficients
          self%data%b(x,y,z,k) = exp(e)
          
          ! Track relative error
          enew = exp(e - n4) - 1.d0

          if (abs(enew) > abs(emax)) then
            emax = enew
            epos = [k, x, y, z]
          endif
        end do
        ! Enforce strict macroscopic particle conservation after exponentiation
        self%data%b(x,y,z,:) = self%data%b(x,y,z,:) * &
            (self%data%totn(x,y,z) / sum(self%data%b(x,y,z,:) * self%data%nstar(x,y,z,:)))
      end do
      end do
      end do
    end associate
  END SUBROUTINE apply_ng_acceleration

END MODULE atom_mod
  
! *******************************************************************
  