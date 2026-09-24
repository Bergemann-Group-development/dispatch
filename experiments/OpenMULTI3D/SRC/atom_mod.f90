! *******************************************************************
MODULE atom_mod
    USE io_mod
    USE handyfuncs_mod
    USE m3d_debug_mod
    USE prec_mod, only: inf8
    USE const_mod
    USE trans_mod
    USE op_math_mod, only: lineq, get_interval, eqsyst, lineq_no_pivot
    USE interpol_mod, only: hermite_interpolation, hermite_derivative
    USE level_mod, only: state_type
    USE ion_mod, only: ion_type, massive_type, debye_shield, massive_type, particle_type
    USE element_mod, only: element_type, isotope_type
    USE molecule_mod, only: molec_type
    USE composition_mod, only: composition_type
    USE line_types_mod, only: linelist_type
    USE linelist_mod, only: append_line, finalize_linelist, &
                            reserve_linelist, ll_bij_const => bij_const
    USE os_mod
    USE trace_mod
    USE console_mod
    USE spectrum_mod, only: spectrum_type, set_weights
    USE hprofl_mod, only:prepare_hydrogen_damping
    USE ng_mod, only: ng_state_t, ng_promote, ng_should_fire, ng_coefficients, &
                      ng_accumulate, ng_extrapolate, ng_count_clamped, ng_condition
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
      ! SoA copy of the (non-hydrogen) bound-bound transitions for the SoA line kernel
      ! path; built by atom_type%build_ll (globally for the linelist_out.txt export,
      ! per rt-patch for the opacity/rate kernels)
      type(linelist_type) :: ll
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

      ! Ng acceleration state, shared by every patch copy of this atom: the
      ! target is allocated once on the master atom and the pointer survives the
      ! copy, so all patches accumulate into one object.
      type(ng_state_t), pointer :: ng => null()


    contains
      procedure :: init_atom_data
      procedure :: check_ionisation_thresholds
      procedure :: ev_keys
      procedure :: lte_pop
      procedure :: compute_lte_populations
      procedure :: normalize_profiles
      procedure :: update_populations
      procedure :: setup_transition_matrix
      procedure :: read_classic
      procedure :: read_pandas
      procedure :: prepare_damping
      procedure :: build_ll
      procedure :: prepare_Drawin
      procedure :: classic_col
      procedure :: gencol
      procedure :: set_colliders
      procedure :: ng_acceleration
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
    logical, public :: use_classic_col = .False., legacy_freql=.False.
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
    !---------------------------------------------------------------------------
    ! Guard rails on the statistical-equilibrium update (see update_populations).
    ! max_update is the largest factor by which a level population may change in
    ! one iteration, in either direction; min_pop_frac keeps a level's population
    ! from ratcheting down into a range where its column of the rate matrix is
    ! numerically empty; conv_pop_floor is the population fraction above which a
    ! level carries full weight in the convergence decision, so that convlim and
    ! convmax are not held hostage by levels with no effect on any opacity.
    !---------------------------------------------------------------------------
    real(8), public :: max_update=1d2, min_pop_frac=1d-30, conv_pop_floor=1d-8
    !---------------------------------------------------------------------------
    ! Under-relaxation of the SE update.  The update is multiplicative (b = b*e),
    ! so the natural damping is geometric: e -> e**sor_w, which is exactly
    ! ln b -> ln b + sor_w*correction and keeps b positive for any sor_w > 0.
    ! If the undamped error map has eigenvalue lambda, the damped one has
    ! mu = 1 + sor_w*(lambda - 1), so a mode that alternates sign and grows at
    ! |lambda| = 1.85 (see fe607_convergence.md) needs sor_w < 2/2.85 = 0.70.
    ! sor_w = 1 is the undamped scheme, i.e. bit-identical to before.
    !---------------------------------------------------------------------------
    real(8), public :: sor_w=1d0
    integer, public :: nq_min=0, nq_max=huge(1)
    real(4), public :: nq_fac=1.0
    ! Bound-bound wavelength grid read from a file (vacuum AA, one value per record,
    ! '#' or '!' comments) instead of being assembled from the per-line q0/qmax
    ! quadratures; those then only fix each line's lam_min/lam_max, i.e. the interval
    ! of the shared grid its rates are integrated over. tumult/line_quadrature.py
    ! writes the grid together with a matching atom copy.
    character(len=256), public :: line_grid_file=''
    real(8), dimension(:), allocatable, private, save :: line_grid
    integer, public :: cswitch=1
    ! Ng acceleration schedule:
    !   < 0  off (default)
    !   = 0  automatic -- accumulate every iteration and let the gate in ng_mod
    !        decide, using the contraction rate, its stationarity, the fraction
    !        of corrections at the update limiter, and the conditioning of the
    !        least-squares fit.  ng_start is then unused.
    !   > 0  fixed schedule, applying on the iteration after each
    !        mod(iter-ng_start, ng_step)==0.
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
      ! iterate history for Ng, needed by the automatic trigger (ng_step==0) as
      ! much as by the fixed schedule
      if (ng_step >= 0) then
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
      ! Normalize wphi; with the new-linelist rate kernel, lines below the LTE
      ! strength cutoff accumulate no wphi (and no gamma) -- leave those at 1 so
      ! the zero radiative rates pass through instead of producing 0*Inf
      do i=1, self%nline
        if (.not. self%trans%line(i)%active) cycle
        where (wphi(:,:,:,i) > 0.)
          wphi(:,:,:,i) = 1 / wphi(:,:,:,i)
        elsewhere
          wphi(:,:,:,i) = 1.
        end where
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
    real(8) :: norm, elo, eup
    real(8), dimension(:), allocatable :: bold, ework
    integer(8) :: ng_nclamp, ng_ntot
    logical :: relax, ok
    integer, save :: itimer=0, itimer2=0
    ! -- diag_pops bookkeeping ---------------------------------------------------
    real(8), dimension(:,:), allocatable :: w_save
    real(8), dimension(:), allocatable :: e_save, x_sol
    integer :: d_nfloor, d_ncell, d_nfcell, d_wpos(3), d_wlvl, d_nneg, d_ipos, j
    integer :: d_npos(3), d_nlvl
    real(8) :: d_resid, d_wresid, d_emax, d_range, d_wrange, d_bmin, d_wb, d_wfrac
    real(8) :: d_nmin, d_nfrac
    ! worst (most negative) off-diagonal: value, its level pair, and where it sits,
    ! so a persistent M-matrix violation can be tied to a region of the atmosphere
    integer :: d_gpos(3), d_glvl(2)
    real(8) :: d_gneg, d_gr
    character(len=256) :: msg
    ! ----------------------------------------------------------------------------
    call trace%begin('atom_type%update_populations', itimer=itimer)
    adata => self%data
    m = shape(adata%totn)
    adata%maxerr = 0
    adata%lvlerr = 0
    nk = self%nlevel
    eup = max(max_update, 1.0d0)
    elo = 1d0/eup
    allocate (bold(nk), ework(nk))
    ! Corrections cut by the limiter, counted unconditionally: the Ng gate needs
    ! it to tell an asymptotic iteration from a clamped one, and it is a count()
    ! over a vector already in cache.
    ng_nclamp = 0
    ng_ntot = 0
    if (diag_pops) then
      allocate (w_save(nk,nk), e_save(nk), x_sol(nk))
      d_nfloor = 0; d_ncell = 0; d_nfcell = 0; d_nneg = 0
      d_wresid = 0; d_emax = 0; d_wrange = 0; d_bmin = huge(1d0)
      d_wpos = 0; d_wlvl = 0; d_wb = 0; d_wfrac = 0
      d_nmin = huge(1d0); d_npos = 0; d_nlvl = 0; d_nfrac = 0
      d_gneg = 0; d_gpos = 0; d_glvl = 0
    endif

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

              if (diag_pops) then
                w_save = adata%w
                e_save = e
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

              if (diag_pops) then
                x_sol = e
                ! backward error of the LU solve, scaled the usual way:
                ! |W x - rhs|_inf / (|W|_inf |x|_inf)
                d_resid = maxval(abs(matmul(w_save, x_sol) - e_save))
                d_range = maxval(abs(w_save)) * maxval(abs(x_sol))
                if (d_range > 0) d_resid = d_resid / d_range
                ! dynamic range of the column scaling n_j that the no-pivot LU sees
                d_range = maxval(b*nstar) / max(minval(b*nstar), tiny(1d0))
                d_ncell = d_ncell + 1
                i = count(x_sol < elo .or. x_sol > eup)
                d_nfloor = d_nfloor + i
                if (i > 0) d_nfcell = d_nfcell + 1
                d_bmin = min(d_bmin, minval(b))
                d_wresid = max(d_wresid, d_resid)
                ! an M-matrix has no negative off-diagonal entries; negative
                ! (preconditioned) rates break the positivity of the solution
                d_ipos = maxloc(b*nstar, dim=1)
                do j=1, nk
                  do i=1, nk
                    if (i /= j .and. i /= d_ipos) then
                      if (w_save(i,j) < 0) then
                        d_nneg = d_nneg + 1
                        ! rank by the entry relative to its own column diagonal,
                        ! so the worst violation is the one that actually matters
                        ! for the solve, not merely the largest raw number
                        d_gr = w_save(i,j) / max(abs(w_save(j,j)), 1d-99)
                        if (d_gr < d_gneg) then
                          d_gneg = d_gr
                          d_gpos = [x, y, z]
                          d_glvl = [i, j]
                        endif
                      endif
                    endif
                  end do
                end do
                ! track the cell with the largest upward correction
                i = maxloc(x_sol, dim=1)
                if (x_sol(i) > d_emax) then
                  d_emax = x_sol(i)
                  d_wrange = d_range
                  d_wpos = [x, y, z]
                  d_wlvl = i
                  d_wb = b(i)
                  d_wfrac = b(i)*nstar(i) / adata%totn(x,y,z)
                endif
              endif

              if (check_NaN) then
                ! only NaN/Inf are fatal here; zero or negative corrections are
                ! absorbed by the exponential damping below
                call check_array(e, 'lineq solution', abort=.false., ok=ok)
                if (.not.ok) then
                  write(msg,'(a,3(1x,i0))') 'invalid lineq solution at x,y,z =', x, y, z
                  call console_error (msg)
                  stop
                endif
              endif
              
              ! Limit how far one iteration may move a population.  Corrections in
              ! [elo, eup] pass through untouched; outside that range both branches
              ! join on smoothly (value and slope), and saturate at 0.37*elo below and
              ! 2*eup above.  The lower branch is also what maps a non-positive
              ! solution onto a small positive factor, i.e. what keeps populations
              ! positive.  Without the upper branch a single cell can multiply a level
              ! by 1e9 in one iteration, and the next formal solution sees that opacity.
              ng_nclamp = ng_nclamp + count(e < elo .or. e > eup)
              ng_ntot   = ng_ntot + nk
              where (e < elo)
                e = elo * exp(max(e, 0.0d0)/elo - 1d0)
              elsewhere (e > eup)
                e = eup * (2d0 - eup/e)
              endwhere

              ! Under-relaxation, applied after the limiter so it acts on a
              ! strictly positive factor.  e**sor_w scales the correction in
              ! log-population space; the conservation renormalisation below then
              ! restores the particle budget that the damped step no longer
              ! satisfies exactly.  sor_w = 1 leaves e untouched.
              if (sor_w /= 1d0) then
                e = e**sor_w
              endif

              ! apply the update
              bold = b
              b = b * e

              ! compute the normalisation scalar
              norm = adata%totn(x,y,z) / sum(b * nstar)

              if (diag_pops) then
                ! how much of the particle budget the update handed to a single
                ! level, and by how much conservation then rescales every level
                if (norm < d_nmin) then
                  d_nmin = norm
                  d_npos = [x, y, z]
                  d_nlvl = maxloc(b*nstar, dim=1)
                  d_nfrac = maxval(b*nstar) / adata%totn(x,y,z)
                endif
              endif

              ! enforce strict macroscopic particle conservation
              b = b * norm

              ! Keep every level numerically alive.  A population that ratchets down
              ! iteration after iteration empties its column of the rate matrix, and
              ! can never come back, since the update is multiplicative (b = b*e).
              ! min_pop_frac of the atoms in the cell is far below anything that
              ! contributes to an opacity or a rate.
              where (nstar > 0) b = max(b, min_pop_frac*adata%totn(x,y,z)/nstar)

              if (check_NaN) then
                call check_array(e, 'population correction', check_zero=.true., abort=.false., ok=ok)
                if (.not.ok) then
                  write(msg,'(a,3(1x,i0))') 'invalid population correction vector at x,y,z =', x, y, z
                  call console_error (msg)
                  stop
                endif
              endif

              if (check_NaN) then
                call check_array(b, 'departure coefficients', check_zero=.true., abort=.false., ok=ok)
                if (.not.ok) then
                  write(msg,'(a,3(1x,i0))') 'invalid departure coefficients at x,y,z =', x, y, z
                  call console_error (msg)
                  stop
                endif
              endif

              ! the true relative change, taken from the populations that were
              ! actually stored, so that the limiter and the floor above are included
              where (bold > 0)
                e = b/bold - 1.0d0
              elsewhere
                e = 0.0d0
              endwhere

              ! keep the full per-level record for the convergence report
              lvlerr = merge(e, lvlerr, abs(e) > abs(lvlerr))

              ! The convergence decision, on the other hand, weights each level by the
              ! share of the atoms it actually holds.  The top levels of the last ion
              ! are radiatively controlled by a near-zero radiation field, so they swing
              ! by orders of magnitude while holding 1e-10 of the iron or less; on an
              ! unweighted maximum they own convlim and convmax outright and the run
              ! never terminates.  Levels above conv_pop_floor keep their full weight,
              ! so this cannot mask a level that contributes to an opacity or a rate.
              ework = abs(e) * min(1.0d0, b*nstar/(conv_pop_floor*adata%totn(x,y,z)))
              i = maxloc(ework, dim=1)
              if (ework(i) > abs(maxerr)) then
                maxerr = sign(ework(i), e(i))
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

    if (diag_pops) then
      if (d_ncell > 0) then
        !$omp critical (diag_pops_print)
        print '(A,A3,A,I6,A,I7,A,I4,A,I9,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,3I4,A,I4)', &
          'DIAG ', self%id, ' clamp=', d_nfloor, '/', d_ncell*nk, ' cells=', d_nfcell, &
          ' negoff=', d_nneg, ' bmin=', d_bmin, ' resid=', d_wresid, &
          ' | emax=', d_emax, ' b=', d_wb, ' n/ntot=', d_wfrac, ' at', d_wpos, ' lvl=', d_wlvl
        print '(A,A3,A,ES9.2,A,ES9.2,A,3I4,A,I4)', &
          'DIAGN ', self%id, ' norm_min=', d_nmin, ' top_frac=', d_nfrac, &
          ' at', d_npos, ' lvl=', d_nlvl
        ! the worst M-matrix violation and where it lives, in global z
        if (d_nneg > 0) then
          print '(A,A3,A,ES9.2,A,I4,A,I4,A,3I4)', &
            'DIAGG ', self%id, ' worst_negoff/diag=', d_gneg, &
            ' i=', d_glvl(1), ' j=', d_glvl(2), ' at', &
            d_gpos + [0, 0, diag_zs-1]
        endif
        !$omp end critical (diag_pops_print)
      endif
      deallocate (w_save, e_save, x_sol)
    endif
    deallocate (bold, ework)

    if (associated(self%ng)) call ng_count_clamped(self%ng, ng_nclamp, ng_ntot)

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
    logical :: abort, ok, db_here
    logical, save :: first_time=.true.
    character(len=256) :: msg
    !------------------------------------------------------------------------------------
    ! resolved once per cell: the transition loop below is the hottest in the
    ! solver and must not carry an array comparison per transition
    db_here = diag_db
    if (db_here) db_here = all([x,y,z] == pcheck)
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
      ! Detailed balance check: deep down the radiation field is Planckian, so the
      ! radiative rates alone have to satisfy n_i* R_ij = n_j* R_ji for every
      ! transition.  Any transition that fails this is out of LTE for a numerical
      ! reason, not a physical one.
      if (db_here) then
        if (rji(kr) /= 0) then
          print '(A,F8.1,A,I6,A,L2,A,I4,A,I4,A,ES11.4,A,ES11.4,A,ES11.4)', &
            'DIAGDB T=', comp%temp(x,y,z), ' kr=', kr, ' cont=', kr>self%nline, &
            ' i=', i, ' j=', j, &
            ' Rij=', rij(kr), ' Rji=', rji(kr), &
            ' ratio=', nstar(i)*rij(kr)/(nstar(j)*rji(kr))
        endif
      endif
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
        write(msg,'(a,3(1x,i0))') 'invalid rate matrix at x,y,z =', x, y, z
        call console_error (msg)
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
    character(len=256) :: msg
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
        write(msg,'(a,a,a,i0,a,f8.4,a,f8.4,a)') 'ionisation threshold of ', trim(self%id), &
          ' ion ', j, ' differs from the EOS value (', xion(j), ' vs ', element%ions(j)%xion, ' eV)'
        call console_warn (ch_m3d, msg)
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
  ! Level energies in the canonical departure-key convention: excitation energy above the
  ! ground state of the level's own ionisation stage, tagged with 100*stage so that levels
  ! of different stages cannot collide in find_closest.  This mirrors the linelist
  ! lower_ev convention (energy above the ground state of that ion), whereas atom%ev is
  ! cumulative from the ground state of the neutral.
  !------------------------------------------------------------------------------------
  function ev_keys(self) result(keys)
    class(atom_type), intent(in) :: self
    real(8), dimension(self%nlevel) :: keys
    real(8) :: ground_ev
    integer :: l
    !---------------------------------------------------------------------------------
    ground_ev = self%ev(1)
    do l=1, self%nlevel
      if (l > 1) ground_ev = merge(self%ev(l), ground_ev, self%ion(l) > self%ion(l-1))
      keys(l) = 100*self%ion(l) + self%ev(l) - ground_ev
    end do
  end function ev_keys
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
    real(8), pointer, contiguous, dimension(:) :: bf_lam, bf_alpha
    integer, pointer, contiguous, dimension(:) :: icont
    integer :: i, j, ic, tri, idx, latom, nnu_bf, ncol, ntri, ii, jj
    integer :: iostat, nnu, ncont, ntemp, l, u
    real(8) :: scr, gv, cnt_ev, nu_max, alpha0, thresh_aa, dlam
    real(8) :: gij, vbar2e4, loggf
    character(512) :: dummy_string, text
    character(8) :: key, str1, str2, str3
    logical :: abort
    integer, save :: itimer=0
    character(len=256) :: msg
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
      if (allocated(in_lines)) then
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
        write(msg,'(a,i0,a)') 'read_atom: cannot parse line (iostat=', iostat, '):'
        call console_error (trim(msg)//' '//trim(text))
        abort = .true.
      endif

      if (disable_bb) call lin%deactivate()
            
      lin%i = lin%i + 1
      lin%j = lin%j + 1

      lin%lam0 = air2vac(lin%lam_air) ! [AA]
      lin%nu0  = HZ_TO_AA / lin%lam0  ! [Hz]
      lin%upper_ev = lin%lower_ev + H_E*lin%nu0
      lin%hydrogen = uppercase(trim(self%id)) == 'H'

      lin%f = 10**loggf / self%g(lin%i)
      gij = self%g(lin%i) / self%g(lin%j)
      call lin%init_einstein_coefficients(lin%f, gij)

      ! make sure we always have an odd number of frequency points for a line
      if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1
      lin%nnu = nq_fac * (lin%nnu-1) + 1
      if (lin%io <= 1 .and. legacy_freql) then
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

      call lin%get_lam(qnorm, lin%lam_input, legacy=legacy_freql)
      call os%mem_register(storage_size(lin%lam_input), lin%nnu, 'line%lam_input')
      lin%lam_min = lin%lam_input(1)
      lin%lam_max = lin%lam_input(lin%nnu)

      if (red_cut < huge(1.0) .or. blue_cut > tiny(1.0)) then
        if (lin%lam_max < blue_cut .or. lin%lam_min > red_cut) then
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
    call os%allocate(bf_lam, nnu_bf)
    call os%allocate(bf_alpha, nnu_bf)

    if (nnu_bf > 0) then
      call nextentry(latom)
      read(latom,*) ! skip bftrans header

      do i=1, nnu_bf
        call nextentry(latom) 
        read(latom,*) icont(i), bf_lam(i), bf_alpha(i)
        icont(i) = icont(i)+1
      end do

      !  bound-free transitions in detail
      do i=1, ncont
        cnt => self%trans%cont(i)
        ii = findloc(icont, i, 1)
        jj = findloc(icont, i, 1, back=.true.)

        if (jj-ii > 0) then
          cnt%nnu = jj-ii+1
          call os%allocate(cnt%lam_input, cnt%nnu, 'cont%lam_input')
          call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')
          
          ! file wavelengths are threshold-decreasing; reverse to increasing
          cnt%lam_input = bf_lam(jj:ii:-1)
          cnt%alpha_input = bf_alpha(jj:ii:-1)
          cnt%alpha0 = cnt%alpha_input(cnt%nnu)

          do j=2,cnt%nnu
            if ( cnt%lam_input(j) < cnt%lam_input(j-1) ) then
              write(msg,'(a,i0,a)') 'read_atom: continuum transition ', i, &
                            ' wavelengths not decreasing'
              call console_error (msg)
              abort = .true.
            elseif ( cnt%lam_input(j) == cnt%lam_input(j-1) ) then
              write(msg,'(a,i0,a)') 'read_atom: continuum transition ', i, &
                            ' wavelengths are not unique; shifted by 1e-6'
              call console_warn (ch_m3d, msg)
              cnt%lam_input(j:nnu) = cnt%lam_input(j:nnu) + 1.0d-6
            endif
          end do

          ! cross sections given in detail
          cnt%bf_type='detail'
        endif
      end do
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
          thresh_aa = HZ_TO_AA / cnt%nu0
          thresh_aa = min(thresh_aa, red_thresh_aa)
          nnu = ceiling((thresh_aa - blue_thresh_aa) / cont_step_aa) + 1
          cnt%nnu = nnu
          call os%allocate(cnt%lam_input, cnt%nnu, 'cont%lam_input')
          call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')
          cnt%lam_input(1:nnu-1) = (/ (blue_thresh_aa + j*cont_step_aa, j=0, nnu-2) /)
          cnt%lam_input(nnu) = thresh_aa
          cnt%alpha_input = Kramer(cnt%lam_input, self%ev(cnt%j) - self%ev(cnt%i))
          cnt%alpha0 = cnt%alpha_input(nnu)
          j = j+1
        else
          call cnt%deactivate()
          cycle
        endif
      endif

      tri = (cnt%j-1)*(cnt%j-2)/2 + cnt%i
      self%cdata(tri) = self%cdata(tri) + cnt%alpha0

      cnt%lam_min = cnt%lam_input(1)
      cnt%lam_max = cnt%lam_input(cnt%nnu)
    end do

    if (j>0) then
      if (console_on(ch_m3d, 1)) write(console_unit,'(i0,a)') j, " bound-free transitions auto-completed using Kramers' formula"
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
          write(msg,'(a,a,2(1x,i0),a)') "ll==ul for collision process ", trim(key), l, u, "; ignored"
          call console_warn (ch_m3d, msg)
          cycle
        endif

        ! within level range
        if(l < 0 .or. u > self%nlevel) then
          write(msg,'(3a,3(a,i0))') 'read_atom: level index outside range (1,nk) for ', &
            trim(key), ':', ' ll=', l, ' ul=', u, ' nk=', self%nlevel
          call console_error (msg)
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
            if (self%ecol(tri) .and. console_on(ch_m3d, 2)) then
              write(msg,'(a,2(1x,i0),1x,a)') 'collisional rate already set for transition', col%ll, col%ul, trim(key)
              call console_warn (ch_m3d, msg)
            endif
            self%ecol(tri) = .true.
          case ('CH', 'CHI')
            if (self%hcol(tri) .and. console_on(ch_m3d, 2)) then
              write(msg,'(a,2(1x,i0),1x,a)') 'collisional rate already set for transition', col%ll, col%ul, trim(key)
              call console_warn (ch_m3d, msg)
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
  SUBROUTINE read_classic (self, atom_file, qnorm, in_lines, read_coll)
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
    real(8) :: scr, gv, cnt_ev, lam_min, alpha0, thresh_aa, dlam
    real(8) :: gij, red_fcut, blu_fcut, vbar2e4, loggf
    integer :: ncont, nnu
    character(len=20) :: txt
    character(len=26),parameter :: upcase_alphabet="ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    character(512) :: text
    character(len=256) :: msg
    logical :: abort
    integer, save :: itimer=0
    !---------------------------------------------------------------------
    call trace%begin('atom_mod%read_classic', itimer=itimer)
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
        call console_warn (ch_m3d, 'more bound-free transitions than ionisable levels in '//trim(self%id))
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
        write(msg,'(a,2(1x,i0))') 'bound-bound transition assigned twice:', lin%i, lin%j
        call console_warn (ch_m3d, msg)
      else
        ij_check(lin%i, lin%j) = .true.
      endif

      if (lin%j > self%nlevel) then
        write(msg,'(a,3(1x,i0))') 'transition to a non-existing level:', i, lin%i, lin%j
        call console_error (msg)
        abort = .true.
      endif

      ! make sure we always have an odd number of frequency points for a line
      if (mod(lin%nnu,2) == 0) lin%nnu=lin%nnu+1

      lin%nnu = nq_fac * (lin%nnu-1) + 1
      if (lin%io <= 1 .and. legacy_freql) then
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

      lin%nu0 = ( lin%upper_ev - lin%lower_ev ) / h_e ! [Hz]
      lin%lam0 = HZ_TO_AA / lin%nu0 ! [AA]

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
        if (console_on(ch_m3d, 2)) write(console_unit,'(2a,i0,a,i0)') uppercase(rname), &
          ': found ', nterm, ' components in transition KR = ', i

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

        lin%dq = lin%lam0 * AA_TO_CM * lin%dq
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

      call lin%get_lam(qnorm, lin%lam_input, legacy=nterm<=1 .and. legacy_freql)
      call os%mem_register(storage_size(lin%lam_input), lin%nnu, 'line%lam_input')
      lin%lam_min = lin%lam_input(1)
      lin%lam_max = lin%lam_input(lin%nnu)

      gij = self%g(lin%i) / self%g(lin%j)
      call lin%init_einstein_coefficients(lin%f, gij)

      if (red_cut < huge(1.0) .or. blue_cut > tiny(1.0)) then
        if (lin%lam_max < blue_cut .or. lin%lam_min > red_cut) then
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
    !  if lam_min .lt.0.0 frequency points in angstrom (starting with
    !  threshold and decreasing) and crossections in cm2 are read.
    do i=1,self%ncont
      call nextentry(latom) 
      read(latom,*,iostat=ios) jj, ii, alpha0, nnu, lam_min
      call print_failed_line(latom, ios)

      if (self%ion(ii) == self%ion(jj)) then
        write(msg,'(a,i0,a,2(1x,i0))') 'bound-free transition ', i, ' connects two levels of the same ion:', ii, jj
        call console_error (msg)
        abort = .true.
      endif

      j  = max(ii,jj)
      ii = min(ii,jj)

      cnt => self%trans%cont(ii)

      if (cnt%bf_type /= 'Kramer') then
        write(msg,'(a,i0,a)') 'multiple bound-free transitions from lower state ', ii, '; ignoring the later one'
        call console_warn (ch_m3d, msg)
        do j=1, nnu
          call nextentry(latom) 
          read(latom, '(A)') ! skipping
        end do
        cycle
      endif
      if (j /= cnt%j) then
        write(msg,'(a,i0,a)') 'bound-free transition ', i, ' has an excited state as upper level'
        call console_warn (ch_m3d, msg)
        cnt%j = j
      endif

      cnt%alpha0 = alpha0
      cnt%nnu    = nnu

      if ( lam_min .lt. 0.0 ) then
        ! cross sections given in detail
        cnt%bf_type='detail'

        call os%allocate(cnt%lam_input, cnt%nnu, 'cont%lam_input')
        call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')

        do j=cnt%nnu,1,-1
          call nextentry(latom) 
          read(latom,*,iostat=ios)  cnt%lam_input(j), cnt%alpha_input(j)
          call print_failed_line(latom, ios)
        end do
        
        do j=2,cnt%nnu
          if ( cnt%lam_input(j) < cnt%lam_input(j-1) ) then
            write(msg,'(a,i0,a)') 'read_atom: continuum transition ', i, &
                          ' wavelengths not decreasing'
            call console_error (msg)
            abort = .true.
          elseif ( cnt%lam_input(j) == cnt%lam_input(j-1) ) then
            write(msg,'(a,i0,a)') 'read_atom: continuum transition ', i, &
                          ' wavelengths are not unique; shifted by 1e-6'
            call console_warn (ch_m3d, msg)
            cnt%lam_input(j:nnu) = cnt%lam_input(j:nnu) + 1.0d-6
          endif
        end do
      else
        ! hydrogenic transitions
        cnt%bf_type='hydrogenic'
        cnt%lam_min = lam_min
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
      cnt%lam0 = HZ_TO_AA / cnt%nu0 ! [AA]

      if (cnt%bf_type == 'hydrogenic') then
        call os%allocate(cnt%lam_input, cnt%nnu, 'cont%lam_input')
        call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')
        dlam = (cnt%lam0 - cnt%lam_min) / (cnt%nnu-1)
        cnt%lam_input = [(cnt%lam0 - i*dlam, i=cnt%nnu-1, 0, -1)]
        cnt%alpha_input = (cnt%lam_input / cnt%lam0)**3.d0 * cnt%alpha0

      elseif (cnt%bf_type == 'Kramer') then
        if (kramer_bf) then
          thresh_aa = min(cnt%lam0, red_thresh_aa)
          nnu = ceiling((thresh_aa - blue_thresh_aa) / cont_step_aa) + 1
          cnt%nnu = nnu
          call os%allocate(cnt%lam_input, cnt%nnu, 'cont%lam_input')
          call os%allocate(cnt%alpha_input, cnt%nnu, 'cont%alpha_input')
          cnt%lam_input(2:nnu) = (/ (blue_thresh_aa + j*cont_step_aa, j=0, nnu-2) /)
          cnt%lam_input(1) = thresh_aa
          cnt%alpha_input = Kramer(cnt%lam_input, self%ev(cnt%j) - self%ev(cnt%i))
          cnt%alpha0 = cnt%alpha_input(nnu)
          j = j+1
        else
          call cnt%deactivate()
          cycle
        endif
      endif

      tri = (cnt%j-1)*(cnt%j-2)/2 + cnt%i
      self%cdata(tri) = cnt%alpha0

      cnt%lam_min = cnt%lam_input(1)
      cnt%lam_max = cnt%lam_input(cnt%nnu)
    end do

    if (j>0) then
      if (console_on(ch_m3d, 1)) write(console_unit,'(i0,a)') j, " bound-free transitions auto-completed using Kramers' formula"
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
      call console_error ('collision routine '//trim(self%crout)//' not supported')
      stop
    end select

    if (console_on(ch_m3d, 1)) write(console_unit,'(a)') 'Reading collision data'
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
          call console_error ('read_collisions: unrecognized keyword '//trim(key))
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
        call console_error ('read_collisions: keyword "END" not found in collision file')
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
          call console_error ('read_collisions: unrecognized keyword '//trim(key))
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
          write(msg,'(a,a,2(1x,i0),a)') "ll==ul for collision process ", trim(key), l, u, "; ignored"
          call console_warn (ch_m3d, msg)
          cycle
        endif

        if (any(self%cgrid(1:nc,i) <= 0)) then
          write(msg,'(2a,2(1x,i0))') 'non-positive cross sections for collision process ', trim(key), l, u
          call console_error (msg)
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
          if (self%ecol(tri) .and. console_on(ch_m3d, 2)) then
            write(msg,'(a,2(1x,i0),1x,a)') 'collisional rate already set for transition', col%ll, col%ul, trim(key)
              call console_warn (ch_m3d, msg)
          endif
          self%ecol(tri) = .true.
        case ('CH', 'CHI')
          if (self%hcol(tri) .and. console_on(ch_m3d, 2)) then
            write(msg,'(a,2(1x,i0),1x,a)') 'collisional rate already set for transition', col%ll, col%ul, trim(key)
              call console_warn (ch_m3d, msg)
          endif
          self%hcol(tri) = .true.
        end select

        ! within level range
        if(col%ll < 0 .or. col%ul > self%nlevel) then
          write(msg,'(3a,3(a,i0))') 'read_collisions: level index outside range (1,nk) for ', &
            trim(col%key), ':', ' ll=', col%ll, ' ul=', col%ul, ' nk=', self%nlevel
          call console_error (msg)
          abort = .true.
        endif
      end if
    end do
    close(latom)
    if (abort) stop 'abort in atom_mod%read'
    call trace%end(itimer)
  END SUBROUTINE read_classic

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
    character(len=256) :: msg
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
          write(msg,'(a,i0,a,i0)') 'isotope ', line%iso, ' not found; deactivating line ', i
          call console_warn (ch_m3d, msg)
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
  ! Convert the (non-hydrogen) bound-bound transitions of trans%line(:) into the SoA
  ! linelist trans%ll for the SoA line kernels (linelist_mod).  Broadening is
  ! copied as prepared by prepare_damping (same over_4pi/nu0 scaling), populations use
  ! the composition ion densities via lnbase/ipart (= ion stage; lower_ev re-referenced
  ! to the stage ground level), and multi-component lines (io>1) are expanded into one
  ! entry per component sharing the same lower/upper levels and transition column.
  !
  ! f_offset present -> per-rt-patch build: the per-line frequency windows
  ! [iblue,ired] (already clamped to the patch range) are stored 1-based relative to
  ! f_offset, and itrans is the subset column (index in trans%line, matching
  ! trans%gamma_ij).  Absent -> global build for the linelist_out.txt export, sorted
  ! by wavelength.
  !------------------------------------------------------------------------------------
  SUBROUTINE build_ll (self, trans, f_offset)
    class(atom_type), intent(in) :: self
    class(atom_trans_t), target :: trans
    integer, intent(in), optional :: f_offset
    class(line_type), pointer :: line
    class(isotope_type), pointer :: isotopes(:)
    real(8) :: gf, gpref0, llam0, lower_ev, dej, frac, loggf, wt
    integer :: i, ih, s, ig, iso, off, k0, n
    !----------------------------------------------------------------------------------
    nullify(isotopes)
    if (associated(self%species)) then
      select type(species => self%species)
      class is (element_type)
        isotopes => species%isotope
      class is (molec_type)
        isotopes => species%isotope
      end select
    endif

    off = 0
    if (present(f_offset)) off = f_offset - 1

    ! pre-size the arrays so the per-line appends do not reallocate
    n = 0
    do i=1, size(trans%line)
      line => trans%line(i)
      if (.not. line%active .or. line%hydrogen) cycle
      n = n + max(1, line%io)
    end do
    if (n == 0) return
    call reserve_linelist(trans%ll, n)

    do i=1, size(trans%line)
      line => trans%line(i)
      if (.not. line%active) cycle
      if (line%hydrogen) cycle     ! hydrogen keeps the old (Stark-profile) formalism

      s  = self%ion(line%i)
      ig = findloc(self%ion, s, 1) ! stage ground level: energy reference matching ion%n
      lower_ev = line%lower_ev - self%ev(ig)

      ! line%f already contains the isotope fraction (prepare_damping); take it back
      ! out of the exported loggf so the record matches a plain input linelist
      frac = 1
      if (line%iso > 0 .and. associated(isotopes)) then
        iso = findloc(isotopes%iso, line%iso, 1)
        if (iso > 0) frac = isotopes(iso)%frac
      endif

      gf     = self%g(line%i) * line%f
      gpref0 = log(gf * ll_bij_const / line%nu0)
      loggf  = log10(gf / frac)
      llam0  = log(line%lam0)
      dej    = H_E * line%nu0
      k0     = trans%ll%nline

      if (line%io > 1) then
        do ih=1, line%io
          wt = line%wt(ih)
          if (wt <= 0) cycle
          call append_line(trans%ll, llam0 - line%dq(ih), line%nu0, lower_ev, &
            gpref0 + log(wt), line%ga + line%da(ih), line%gw1, line%gw2, line%gq, &
            line%unsold, line%i, line%j, i, s, line%Aij*wt, wt, dej, &
            line%iso, loggf + log10(wt), line%gv_in, line%ga_in, line%gq_in)
        end do
      else
        call append_line(trans%ll, llam0, line%nu0, lower_ev, gpref0, &
          dble(line%ga), dble(line%gw1), dble(line%gw2), dble(line%gq), line%unsold, &
          line%i, line%j, i, s, line%Aij, 1.d0, dej, line%iso, &
          loggf, line%gv_in, line%ga_in, line%gq_in)
      endif

      if (present(f_offset)) then
        ! all entries just appended (components) share the line's frequency window;
        ! the window grid is symmetric about the line centre, so its midpoint seeds f0
        trans%ll%iws(k0+1:trans%ll%nline)  = line%iblue - off
        trans%ll%iwe(k0+1:trans%ll%nline)  = line%ired - off
        trans%ll%icen(k0+1:trans%ll%nline) = (line%iblue + line%ired) / 2 - off
      endif
    end do

    if (.not. present(f_offset)) then
      call finalize_linelist(trans%ll)
      ! the old bb formalism evaluates every active line, so flag them all as used
      ! for the linelist_out.txt export (the SoA kernels, which would set the flag
      ! per line, are currently disabled in rt_patch_mod)
      trans%ll%used(1:trans%ll%nline) = .true.
    endif
  END SUBROUTINE build_ll

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
    real(8) :: dxi, dev, ex, u0, f, tmp, sqt, fac, lam0
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
                lam0 = HZ_TO_AA*h_e / dev
                fac = Kramer(lam0, dev)
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
      if (console_on(ch_m3d, 1)) write(console_unit,'(i0,a)') n, ' collisional transitions filled with classical recipes'
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
    real(8), dimension(:), allocatable :: lam
    logical :: use_grid
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('atom_trans_t%append_transitions', itimer=itimer)
    ! compute total number of frequencies
    use_grid = option(bb) .and. line_grid_file /= ''
    if (use_grid .and. .not. allocated(line_grid)) call read_line_grid()
    nnu = sp%nnu
    if (use_grid) then
      nnu = nnu + size(line_grid)
    elseif (option(bb)) then
      nnu = nnu + sum(self%line%nnu, self%line%active)
    endif
    if (option(bf)) nnu = nnu + sum(self%cont%nnu, self%cont%active)

    ! save all frequencies
    allocate(lam(nnu))
    if (allocated(sp%lam)) lam(1:sp%nnu) = sp%lam

    je=sp%nnu
    if (use_grid) then
      ! the caller sorts the grid, and set_ired_iblue then locates every line on it
      js=je+1
      je=js+size(line_grid)-1
      lam(js:je)=line_grid
    elseif (option(bb)) then
      do i=1, size(self%line)
        lin => self%line(i)
        if (.not. lin%active) cycle
        js=je+1
        je=js+lin%nnu-1
        lam(js:je)=lin%lam_input
        lin%iblue = js
        lin%ired = je
      end do
    endif
    if (option(bf)) then
      do i=1, size(self%cont)
        cnt => self%cont(i)
        if (.not. cnt%active) cycle
        js=je+1
        je=js+cnt%nnu-1
        lam(js:je)=cnt%lam_input
        cnt%iblue = js
        cnt%ired = je
      end do
    endif

    ! now set overall wavelength grid
    call sp%reset(lam=lam)
    call trace%end(itimer)
  END SUBROUTINE append_transitions

  !===============================================================================
  ! Read line_grid_file into line_grid: one vacuum wavelength [AA] per record,
  ! blank lines and lines starting with '#' or '!' skipped.
  !===============================================================================
  SUBROUTINE read_line_grid
    character(len=256) :: text
    integer :: lunit, iostat, n
    !-----------------------------------------------------------------------------
    call os%check_file(line_grid_file)
    open(newunit=lunit, file=trim(line_grid_file), form='FORMATTED', status='old')
    n = 0
    do
      read(lunit, '(A)', iostat=iostat) text
      if (iostat /= 0) exit
      if (is_grid_record(text)) n = n + 1
    end do
    if (n == 0) then
      call console_error ('no wavelengths in line_grid_file '//trim(line_grid_file))
      stop 'abort in atom_mod%read_line_grid'
    endif
    allocate(line_grid(n))
    rewind(lunit)
    n = 0
    do
      read(lunit, '(A)', iostat=iostat) text
      if (iostat /= 0) exit
      if (.not. is_grid_record(text)) cycle
      n = n + 1
      read(text, *) line_grid(n)
    end do
    close(lunit)
    if (console_on(ch_m3d, 1)) write(console_unit,'(A,I0,A,A)') 'Line grid: ', n, ' wavelengths from ', trim(line_grid_file)
  contains
    logical function is_grid_record(text)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: t
      t = adjustl(text)
      is_grid_record = len_trim(t) > 0 .and. t(1:1) /= '#' .and. t(1:1) /= '!'
    end function is_grid_record
  END SUBROUTINE read_line_grid

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
        if (ired == 0 .or. iblue >= sp%nnu) then
          call cnt%deactivate()
        elseif (iblue == ired) then
          if (console_on(ch_m3d, 2)) write(console_unit,'(a,i0,a)') 'continuum ', i, ' covers a single frequency; deactivated'
          call cnt%deactivate()
        else
          ! total number of frequencies where this continuum is active
          nnu = ired - iblue + 1
          cnt%nnu = nnu
          ! min and max frequency where this continuum is active
          cnt%ired  = ired
          cnt%iblue = iblue

          ! set new frequency grid
          allocate(cnt%alpha(iblue:ired))
          call os%mem_register(storage_size(cnt%alpha), cnt%nnu)
          cnt%lam => sp%lam(iblue:ired)
          cnt%nu => sp%nu(iblue:ired)
          
          ! recompute the bf cross sections
          select case(cnt%bf_type)
          case('hydrogenic', 'Kramer')
            ! use nu**-3 dependency
            cnt%alpha = (cnt%lam / cnt%lam0)**3.d0 * cnt%alpha0
          case('detail')
            ! fill new cross sections (interpolated in frequency, see hermite_interp)
            call cnt%hermite_interp(cnt%alpha, cnt%nu)
          end select

          allocate(cnt%wnu(iblue:ired))
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
        if (ired == 0 .or. iblue == sp%nnu) then
          call lin%deactivate()
        else
          if (iblue == ired) then
            if (console_on(ch_m3d, 2)) write(console_unit,'(a,i0,a)') 'line ', i, ' covers a single frequency'
          endif
          ! total number of frequencies where this line is active
          nnu = ired - iblue + 1
          lin%nnu = nnu
          ! min and max frequency where this line is active
          lin%ired  = ired
          lin%iblue = iblue

          lin%lam => sp%lam(iblue:ired)
          lin%nu => sp%nu(iblue:ired)
          allocate(lin%q(iblue:ired))
          call os%mem_register(storage_size(lin%q), nnu)
          ! frequency-oriented offset (~ +dnu/nu0) so the Doppler sign matches nu space
          lin%q(iblue:ired) = 1.0d0 - lin%lam/lin%lam0

          ! fill frequency integration weights for the line
          allocate(lin%wnu(iblue:ired))
          call os%mem_register(storage_size(lin%wnu), lin%nnu, 'line%wnu')
          call set_weights(lin%wnu, lin%nu)
          lin%wnu = lin%wnu / (h_4pi*lin%nu)

          allocate(lin%wq(iblue:ired))
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
        sp%nac(cnt%iblue:cnt%ired) = sp%nac(cnt%iblue:cnt%ired) + 1 ! number of active continua at frequency index j
      end do
    endif

    if (option(bb)) then
      do i=1, size(self%line)
        lin => self%line(i)
        if (.not. lin%active) cycle
        sp%nal(lin%iblue:lin%ired) = sp%nal(lin%iblue:lin%ired) + 1
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
  ELEMENTAL FUNCTION Kramer (lam, dev) result(hyd)
    real(8), intent(in) :: lam, dev
    real(8) :: nu, nl, hyd
    real(8), parameter :: kramer_const = 64*PI**4 * EC**10 * EM / (3*sqrt(3.d0) * CC * HH**6) ! 2.815e29
    !---------------------------------------------------------------------
    nu = HZ_TO_AA / lam
    if (dev > 0) then
      nl = sqrt(RydV/dev)
      hyd = kramer_const * gaunt(lam, ceiling(nl)) / (nl**5 * nu**3)
    else
      hyd = 0.
    endif
  END FUNCTION Kramer

  !====================================================================================
  ELEMENTAL FUNCTION gaunt (lam, n)
    real(8), intent(in) :: lam
    integer, intent(in) :: n
    real(8) :: inv_nu, gaunt
    real(8), parameter, dimension(6) :: &
    a1 = [ 1.078e00, 1.0926e0, 1.0983e0, 1.0954e0, 1.0912e0, 1.0925e0], &
    b1 = [-8.754e14,-2.019e14,-9.45e13, -5.188e13,-3.2e13,  -2.331e13], &
    c1 = [-1.791e29, 1.836e28, 9.177e27, 3.552e27, 1.576e27, 9.325e26], &
    a2 = [ 0.798e00, 0.768e00, 0.793e00, 0.831e00, 0.758e00, 0.79e00],  &
    b2 = [ 5.358e15, 6.242e15, 5.48e15,  4.094e15, 6.633e15, 5.808e15], &
    c2 = [-3.484e31,-3.208e31,-2.318e31,-1.43e31, -3.32e31, -2.844e31]
    !---------------------------------------------------------------------
    inv_nu = lam * AA_TO_HZ
    if (n > 6) then
      gaunt = 1
    else
      if (inv_nu >= 1e-16) then
          gaunt = a1(n) + (b1(n) + c1(n) * inv_nu) * inv_nu
      else
          gaunt = a2(n) + (b2(n) + c2(n) * inv_nu) * inv_nu
      endif
    endif
  END FUNCTION gaunt

  !===============================================================================
  ! Ng acceleration driver.  Accumulation is unconditional and cheap (a handful
  ! of inner products against an LU solve per cell), so the schedule is not what
  ! decides whether to accelerate -- ng_should_fire is, from the quantities the
  ! accumulation already produces.
  !
  ! ng_step > 0 keeps the fixed schedule: accumulate when mod(iter-ng_start,
  ! ng_step)==0, apply on the iteration after.  ng_step == 0 selects automatic
  ! mode, where every iteration accumulates and the gate decides.  ng_step < 0
  ! disables Ng entirely and this is never reached.
  !
  ! In both modes application trails accumulation by one iteration, because a
  ! patch cannot know when the other patches have finished contributing to the
  ! Gram matrix.
  !===============================================================================
  SUBROUTINE ng_acceleration (self, iter, err, convlim, emax, epos)
    class(atom_type):: self
    integer, intent(in) :: iter
    real(8), intent(in) :: err                ! current convergence measure
    real(8), intent(in) :: convlim            ! its target, for the headroom test
    real(8), intent(out) :: emax
    integer, intent(out) :: epos(4)
    logical :: fire, applied, apply_here
    character(len=15) :: ng_reason
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin ('atom_type%ng_acceleration', itimer=itimer)
    emax = 0
    epos = 0
    if (.not. associated(self%ng)) then
      call trace%end(itimer)
      return
    endif

    associate (ng => self%ng, adata => self%data)
      ! Promote last iteration's partial sums, decide, and read the decision back
      ! -- all under one lock, so a patch cannot act on a half-updated state.
      ! Patches advance their own iteration counters, so both the roll and the
      ! decision are monotone in iter: a straggler still on an earlier iteration
      ! neither redoes them nor undoes them, and simply does not accelerate.
      apply_here = .false.
      !$omp critical (ng_shared)
      call ng_promote(ng, iter)
      if (iter > ng%iter_decided) then
        ng%iter_decided = iter
        if (ng_step > 0) then
          fire = mod(iter - ng_start, ng_step) == 1
          if (fire) fire = ng%s_done > 0d0
          ng_reason = merge('FIRE     ', 'off-cycle', fire)
        else
          fire = ng_should_fire(ng, iter, err, convlim, ng_reason)
        endif
        if (fire) then
          call ng_coefficients(ng)
          ng%iter_fired = iter
        endif
        if (diag_pops) then
          ! why the gate did what it did: the contraction rate and its
          ! predecessor (stationarity), the limiter fraction (linearity), and the
          ! scaled Gram determinant (conditioning)
          write(stdout,'(A,I4,A,A15,A,F8.4,A,F8.4,A,ES10.3,A,ES10.3,A,ES10.3)') &
            'DIAGNG iter=', iter, ' verdict=', ng_reason, ' r=', ng%r, &
            ' r_prev=', ng%r_prev, ' clamp=', &
            real(ng%nclamp_done,8)/max(real(ng%ntot_done,8), 1d0), &
            ' cond=', ng_condition(ng), ' err=', err
          flush(stdout)
        endif
      endif
      apply_here = ng%coeff_valid .and. ng%iter_fired == iter
      !$omp end critical (ng_shared)

      if (apply_here) then
        call ng_extrapolate(ng, adata%b, adata%oldb, adata%nstar, adata%totn, &
                            self%nlevel, emax, epos, applied)
      endif

      ! Feed this iteration's differences in for the next decision.  Skipped on
      ! an iteration that extrapolated: b no longer belongs to the same
      ! geometric sequence as oldb, so the differences would be meaningless.
      if (ng%iter_fired /= iter .and. allocated(adata%oldb)) then
        call ng_accumulate(ng, adata%b, adata%oldb, self%nlevel)
      endif
    end associate
    call trace%end(itimer)
  END SUBROUTINE ng_acceleration

END MODULE atom_mod
  
! *******************************************************************
  