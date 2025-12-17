MODULE composition_patch_mod
  USE io_mod
  USE os_mod
  USE omp_mod
  USE task_mod, only: task_t
  USE link_mod, only: nbor_t
  USE gpatch_mod
  USE handyfuncs_mod
  USE optab_mod, only: optab
  USE osmet_mod, only: osmet
  USE ion_mod, only: particle_type, ion_type, debye_shield
  USE molecule_mod, only: molec_type
  USE element_mod, only: element_type
  USE m3d_patch_mod,   only: m3d_patch_t, composition, comoving_frame, ilambd, dlam, &
  use_rho, use_ne, overlapping_bb, overlapping_bf, binned_opac, &
  maxiter, rotate_atmos, long_vec, decouple_continuum, n_nu, use_vel, &
  linelist_in_output, make_eos, overlapping_bb_in_output, atmos_format, atmos_file
  USE m3d_debug_mod
  USE composition_mod, only: composition_type, exclude_elements, only_elements, switch, &
  irwin_dis, ref_lam, debug_opac, abundances, relative_abundances, mhd_eos, helium_scattering, &
  departure_pointer_t, read_departures, departure_indices, sample_back_opac, &
  sample_back_opac_in_output, sample_linelist, opac_table_vmic
  USE atom_mod, only: m1d_legacy_mode, use_classic_col, level_partf
  USE nlte_patch_mod,  only: nlte_patch_t
  USE rt_patch_mod, only:rt_patch_t
  USE bftrans_mod, only: H_occupation_probability
  USE bbtrans_mod, only: tolerance, bbtrans_type, disable_ga, disable_stark, disable_vdw, ga_fallback
  USE hprofl_mod, only: hydrogen_BPO, stark_convolution
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE ray_mod, only: set_carlson_quadratures
  USE list_mod, only: list_t
  USE link_mod, only: nbor_t
  USE bits_mod
  USE atom_mod, only: atom_type, atomic_spectrum_type
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  type, public, extends(m3d_patch_t) :: composition_patch_t
    class(gpatch_t), pointer :: mhd ! pointer to MHD patch
    class(nlte_patch_t), pointer :: nlte_patch
  contains
    procedure :: pre_init => composition_pre_init
    procedure :: execute
    procedure :: create_sub_patches
    procedure :: print_post_update_info
    procedure :: save_abundances_to_file
    procedure :: init_tabulated_comp
    procedure :: init_eos_comp
    procedure :: save_eos_variables
    procedure :: save_number_densities
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type
  
  real(8), public :: FeH=-99.9, dtemp=250., ldne=0.25, ldrho=0.25, ldtemp=0.1
  logical, public :: save_atmos=.true., save_eos=.false.
  real(4), public, dimension(:), allocatable :: external_abund
  real(4), public :: alpha_enhancement=0.0, r_process_enhancement=0.0
  character(len=256), save :: save_species=''
  character(len=256), save :: absdat_file='./input_multi3d/TS_absdat.dat'
  character(len=256), save :: abund_file='./input_multi3d/abund/abund_magg'
  character(len=256), save :: isotope_file='./input_multi3d/isotopes'
  character(len=256), save :: molec_file='./input_multi3d/molec_data_BC.txt'
  character(len=256), save :: opac_table=''
  character(len=256), save :: absmet_file=''
  character(len=512), save :: linelist_folder=''
  character(len=256), save, dimension(:), allocatable :: line_lists
  
  public :: read_composition_params

contains

!===============================================================================
SUBROUTINE read_composition_params ()
  integer :: i, iostat
  real(4) :: tmolim=-99
  integer, save:: itimer=0, verbose=0
  namelist /composition_params/ verbose, absdat_file, abund_file, absmet_file, &
  isotope_file, molec_file, H_occupation_probability, opac_table, m1d_legacy_mode, &
  sample_back_opac, sample_back_opac_in_output, switch, mhd_eos, debye_shield, &
  dtemp, ldne, ldrho, ldtemp, tmolim, irwin_dis, opac_table_vmic, save_species, &
  ref_lam, debug_opac, abundances, relative_abundances, alpha_enhancement, r_process_enhancement, &
  disable_ga, disable_stark, disable_vdw, helium_scattering, sample_linelist, hydrogen_BPO, &
  stark_convolution, ga_fallback
  namelist /linelist_params/ dlam, exclude_elements, line_lists, linelist_folder, &
  linelist_in_output, tolerance, only_elements
  ! ----------------------------------------------------------------------------
  call trace%begin('composition%init', itimer=itimer)

  allocate(line_lists(128))
  line_lists = ''

  rewind (io%input)
  read (io%input, composition_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, composition_params)
  end if
  call io%namelist_warning ('composition_params', iostat)

  if (m1d_legacy_mode) then
    decouple_continuum = .true.
    overlapping_bb = .false.
    overlapping_bf = .false.
    use_rho = .false.
    use_vel(1:3) = .false.
    sample_back_opac = .false.
    h_occupation_probability = .false.
    use_classic_col = .false.
    irwin_dis = .false.
    level_partf = .true.
    helium_scattering = .false.
    !comoving_frame = .true.
    !preconditioning = .false.
  endif

  write (io%namelist, composition_params)

  rewind (io%input)
  read (io%input, linelist_params, iostat=iostat)
  if (iostat > 0) then
    rewind (io%input)
    read (io%input, linelist_params)
  end if
  call io%namelist_warning ('linelist_params', iostat)

  if (linelist_folder /= '') then
    deallocate(line_lists)
    call os%get_filenames(linelist_folder, line_lists)
    line_lists = [(trim(linelist_folder)//'/'//trim(line_lists(i)), i=1, size(line_lists))]
  else
    call truncate_array(line_lists, count(line_lists/=''))
  endif
  write (io%namelist, linelist_params)

  if (verbose > 0) write (*,*) ' - Initializing the NEW opacity package'
  composition%verbose = verbose

  if (absmet_file /= '') then
      call osmet%init(absmet_file)
  endif

  call composition%init(absdat_file, abund_file, FeH, &
  isotope_file=isotope_file, molec_file=molec_file, tmolim=tmolim, &
  external_abund=external_abund, alpha_enhancement=alpha_enhancement, &
  r_process_enhancement=r_process_enhancement)
  call trace%end(itimer)
END SUBROUTINE read_composition_params

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  class(nlte_patch_t), pointer :: nlte_patch
  class(atomic_spectrum_type), save, pointer :: output_sp
  integer :: i, mu, nw, n_mu
  integer, save:: itimer=0
  logical, save :: first_time=.true., external_sp
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%composition_pre_init', itimer=itimer)
  self%kind = 'composition'
  call self%m3d_patch_t%pre_init()
  nw = merge(1+size(long_vec,2), 1, rotate_atmos)
  call self%set_nv_and_nw(8,nw)
  self%static_mem=.true.

  if (first_time) then
    allocate(output_sp)
    binned_opac = opac_table /= ''
      
    if (binned_opac) then
      call optab%init(opac_table, output_sp)
    else
      call output_sp%init()
      output_sp%name = 'output'
    endif

    external_sp = output_sp%nnu > 0

    if (n_nu < 0) then
      n_nu = omp%nthreads
    endif
  endif

  allocate(nlte_patch)
  call nlte_patch%pre_init_atom()
  call nlte_patch%pre_init_spectrum(output_sp)

  if (associated(nlte_patch%atoms) .or. read_departures) then
    call nlte_patch%pre_init_and_link(self)
    self%nlte_patch => nlte_patch
  endif

  if (first_time) then
    first_time = .false.
    do i=1, size(line_lists)
      call composition%init_linelist(line_lists(i), dlam, output_sp%nu)
    end do
  endif

  ! make sure there are no more f-domain splits than threads
  if (output_sp%nnu > 0) n_nu = min(n_nu, output_sp%nnu)
  if (decouple_continuum) then
    if (nlte_patch%bb_sp%nnu > 0) n_nu = min(n_nu,  nlte_patch%bb_sp%nnu)
    if (nlte_patch%bf_sp%nnu > 0) n_nu = min(n_nu,  nlte_patch%bf_sp%nnu)
  else
    if (nlte_patch%sp%nnu > 0) n_nu = min(n_nu,  nlte_patch%sp%nnu)
  endif

  self%n_nu = n_nu
  n_mu = size(long_vec,2)
  if (associated(nlte_patch%atoms)) then
    if (maxiter > 0) then
      if (decouple_continuum) then
        call self%create_sub_patches([1, n_nu],      [1,2], nlte_patch%bb_sp, nlte_patch%bb_atoms)
        call self%create_sub_patches([1, n_nu]+n_nu, [1,2], nlte_patch%bf_sp, nlte_patch%bf_atoms)
      else
        call self%create_sub_patches([1, n_nu], [1,2], nlte_patch%sp, nlte_patch%atoms)
      endif
    endif
  endif

  if (output_sp%nnu > 0) then
    if (n_mu > 0 .or. (make_eos .and. self%key(1)==1)) then
      if (ilambd > 0) then
        call self%create_sub_patches([-n_nu, -1], [1,2], output_sp, nlte_patch%output_atoms)
      endif
      call self%create_sub_patches([-n_nu, -1], [-1,-1], output_sp, nlte_patch%output_atoms)
    endif
  endif

  if (.not. associated(self%nlte_patch)) then
    deallocate(nlte_patch)
  endif

  if (save_species /= '') then
    call self%save_number_densities(composition, pre_init=.true.)
  endif

  call self%save_abundances_to_file()
  call self%post_pre_init()
  call trace%end(itimer)
END SUBROUTINE composition_pre_init

!===============================================================================
SUBROUTINE create_sub_patches (self, nu, mu, sp, atoms)
  class(composition_patch_t) :: self
  integer, intent(in) :: nu(2), mu(2)
  class(atomic_spectrum_type), intent(in), pointer :: sp
  type(atom_type), intent(in), pointer, optional :: atoms(:)
  class(rt_patch_t), pointer :: rt_patch
  integer :: n, m
  logical :: new_nlte
  !-----------------------------------------------------------------------------
  if (.not. associated(sp)) return
  if (sp%nnu == 0) return
  do m=mu(1), mu(2)
    do n=nu(1), nu(2)
      allocate(rt_patch)
      call rt_patch%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp, atoms=atoms)
    end do
  end do
END SUBROUTINE create_sub_patches

!===============================================================================
SUBROUTINE execute (self)
  class(composition_patch_t) :: self
  class(composition_type), pointer :: comp
  real(4), dimension(:,:,:), pointer :: d, e, pe, chi
  real(4), dimension(:,:,:,:), pointer, contiguous :: authority
  real(4), dimension(:,:,:), pointer, contiguous :: ne, rho, temp, vx, vy, vz
  real(4), dimension(:,:,:), pointer, contiguous :: ux, uy, uz, vmic
  real(8), dimension(:,:,:), allocatable :: ne8, rho8, temp8
  real(8) :: vec(3)
  integer :: i, n(3), li(3), ui(3)
  integer, save:: itimer=0
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%execute', itimer=itimer)
  n = self%n

  if (self%iter == 0) then
    li = self%li
    ui = self%ui

    allocate(self%composition)
    if (rotate_atmos) allocate(self%rotated_composition(self%n_out))

    allocate(temp8(n(1),n(2),n(3)))
    allocate(rho8 (n(1),n(2),n(3)))
    allocate(ne8  (n(1),n(2),n(3)))

    if (atmos_format == 'save') then
      call self%load5('atmos', run_dir=atmos_file, mem=self%mem(:,:,:,:,self%it,:))
    else
      do i=1, self%nw
        authority => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,self%it,i)
        e  => authority(:,:,:,self%patch%idx%e)
        d  => authority(:,:,:,self%patch%idx%d)
        ux => authority(:,:,:,self%patch%idx%px)
        uy => authority(:,:,:,self%patch%idx%py)
        uz => authority(:,:,:,self%patch%idx%pz)
        pe => authority(:,:,:,6)
        vmic => authority(:,:,:,7)
      
        authority => self%mem(:,:,:,:,self%it,i)
        ne   => authority(:,:,:,self%i_ne)
        rho  => authority(:,:,:,self%i_rho)
        temp => authority(:,:,:,self%i_temp)
        vx   => authority(:,:,:,self%i_vx)
        vy   => authority(:,:,:,self%i_vy)
        vz   => authority(:,:,:,self%i_vz)
        chi  => authority(:,:,:,7)
        authority(:,:,:,8) = vmic
    
        temp = e
        rho  = d
        ne   = pe / (cgs%k_b * temp)
        vx   = ux / cgs%c
        vy   = uy / cgs%c
        vz   = uz / cgs%c
      end do
    endif

    do i=1, self%nw
      authority => self%mem(:,:,:,:,self%it,i)
      ne   => authority(:,:,:,self%i_ne)
      rho  => authority(:,:,:,self%i_rho)
      temp => authority(:,:,:,self%i_temp)
      chi  => authority(:,:,:,7)
      vmic => authority(:,:,:,8)
      
      temp8 = dble(temp)
      rho8  = dble(rho)
      ne8   = dble(ne)

      ! copy composition with all the information from the abund
      ! and absdat file to make it threadprivate
      if (i==1) then
        self%composition = composition
        comp => self%composition
      else
        self%rotated_composition(i-1) = composition
        comp => self%rotated_composition(i-1)
      endif

      ! this is needed, otherwise all pointers will still point
      ! to ions and molecules in original composition
      call comp%reassign_pointer()
      ! loading temp, ne and rho in to composition objekt
      ! ion and molecule equilibrium are solved immediately
      call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_rho, use_ne=use_ne, vmic=vmic)

      if (.not. use_ne) ne = real(comp%ne)
      if (.not. use_rho) rho = real(comp%rho)
      if (save_atmos) then
        chi = comp%chi500
      endif
    end do
    
    if (save_atmos) then
      call self%save5(self%mem(:,:,:,:,self%it,:), 'atmos')
    endif

    if (save_eos) then
      call self%save_eos_variables()
    endif

    if (save_species /= '') then
      call self%save_number_densities(self%composition, pre_init=.false.)
    endif

    if (sample_back_opac .or. sample_back_opac_in_output) then
      call self%init_tabulated_comp(self%mem(:,:,:,self%i_temp,self%it,:), &
                                    self%mem(:,:,:,self%i_ne,  self%it,:), &
                                    dtemp=dtemp, ldne=ldne)
    elseif (make_eos) then
      call self%init_eos_comp(self%composition%temp, self%composition%rho, ldtemp=ldtemp, ldrho=ldrho)
    endif

    self%dtime = max(self%dtime, io%end_time - self%dtime)
    self%out_ray = .true. ! so that it is still executed after convergence
  else
    self%dtime = self%patch%dtime
    !print *, self%m3d_patch_info()
    !print *, 'updating composition more than once is not implemented yet. Stopping...'
    !stop
  endif
  call trace%end(itimer)
END SUBROUTINE execute

!===============================================================================
SUBROUTINE init_tabulated_comp (self, temp, ne, dtemp, ldne)
    class(composition_patch_t), target :: self
    real(4), dimension(:,:,:,:), intent(in) :: temp, ne
    real(8), intent(in) :: dtemp, ldne
    integer :: i, j, ntemp, nne, nnu
    real(8), dimension(:,:,:), allocatable :: tab_temp, tab_ne, tab_rho
    real(4), dimension(:,:,:), allocatable :: tab_vmic
    real(8) :: mint, maxt, mine, maxe
    !---------------------------------------------------------------------------
    allocate(self%tab_comp)
    self%tab_comp = composition
    call self%tab_comp%reassign_pointer()
    mint = minval(temp)
    maxt = maxval(temp)
    mine = log(minval(ne))
    maxe = log(maxval(ne))
    ntemp = ceiling(maxt/dtemp) - floor(mint/dtemp) + 1
    nne   = ceiling(maxe/ldne)  - floor(mine/ldne)  + 1
    mint = floor(mint/dtemp)*dtemp
    mine = floor(mine/ldne)*ldne
    allocate(tab_temp(ntemp, nne, 1))
    allocate(tab_rho (ntemp, nne, 1))
    allocate(tab_ne  (ntemp, nne, 1))
    allocate(tab_vmic(ntemp, nne, 1))
    do j=1, nne
    do i=1, ntemp
      tab_temp(i,j,1) = mint + (i-1)*dtemp
      tab_ne  (i,j,1) = mine + (j-1)*ldne
    end do
    end do
    tab_ne = exp(tab_ne)
    tab_vmic = opac_table_vmic
    call self%tab_comp%load_state_variables(tab_temp, tab_ne, tab_rho, vmic=tab_vmic, &
                                            use_ne=.true., use_rho=.false.)
END SUBROUTINE init_tabulated_comp

!===============================================================================
SUBROUTINE init_eos_comp (self, temp, rho, ldtemp, ldrho)
  class(composition_patch_t), target :: self
  real(8), dimension(:,:,:), intent(in) :: temp, rho
  real(8), intent(in) :: ldtemp, ldrho
  real(8) :: minTemp, maxTemp, minRho, maxRho
  integer :: i, j, ntemp, nrho, nnu
  real(8), dimension(:,:,:), allocatable :: tab_temp, tab_ne, tab_rho
  real(4), dimension(:,:,:), allocatable :: eos, tab_vmic
  real(8) :: minLnRho, maxLnRho, minLnT, maxLnT
  !---------------------------------------------------------------------------
  allocate(self%tab_comp)
  self%tab_comp = composition
  call self%tab_comp%reassign_pointer()

  minTemp = minval(temp)
  maxTemp = maxval(temp)
  minRho = minval(rho)
  maxRho = maxval(rho)

  minLnRho = log(minRho)
  maxLnRho = log(maxRho)
  minLnT = log(minTemp)
  maxLnT = log(maxTemp)

  ntemp = ceiling((maxLnT - minLnT) / ldtemp) 
  nrho  = ceiling((maxLnRho - minLnRho) / ldrho) 

  allocate(tab_temp(ntemp, nrho, 1))
  allocate(tab_rho (ntemp, nrho, 1))
  allocate(tab_ne  (ntemp, nrho, 1))
  allocate(tab_vmic (ntemp, nrho, 1))

  do j=1, nrho
  do i=1, ntemp
    tab_temp(i,j,1) = minLnT + (i-1)*ldtemp
    tab_rho (i,j,1) = minLnRho + (j-1)*ldrho
  end do
  end do

  tab_rho = exp(tab_rho)
  tab_temp = exp(tab_temp)
  tab_vmic = opac_table_vmic
  call self%tab_comp%load_state_variables(tab_temp, tab_ne, tab_rho, use_rho=.true., use_ne=.false., vmic=tab_vmic)

  allocate(eos(6, ntemp, nrho))
  eos(1, :, :) = safe_to_real4(self%tab_comp%temp(:, :, 1))
  eos(2, :, :) = safe_to_real4(self%tab_comp%rho(:, :, 1))
  eos(3, :, :) = safe_to_real4(self%tab_comp%ne(:, :, 1))
  eos(4, :, :) = safe_to_real4(self%tab_comp%pg(:, :, 1))
  eos(5, :, :) = safe_to_real4(self%tab_comp%E(:, :, 1))
  eos(6, :, :) = self%tab_comp%chi500(:, :, 1)
  call self%save3(eos, 'eos')

contains
  elemental function safe_to_real4(d_val) result(s_val)
    real(8), intent(in) :: d_val
    real(4) :: s_val
    s_val = real(max(min(d_val, dble(huge(1.0_4))), -dble(huge(1.0_4))), 4)
  end function safe_to_real4
END SUBROUTINE init_eos_comp

!===============================================================================
SUBROUTINE save_eos_variables (self)
  class(composition_patch_t), target :: self
  class(composition_type), pointer :: comp
  real(4), dimension(:,:,:,:), allocatable :: eos
  integer :: nx, ny, nz
  integer, save :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('composition_patch%save_eos_variables', itimer=itimer)
  comp => self%composition
  nx = size(comp%temp,1)
  ny = size(comp%temp,2)
  nz = size(comp%temp,3)
  allocate(eos(nx,ny,nz,2))
  eos(:,:,:,1) = comp%pg
  eos(:,:,:,2) = comp%E
  call self%save4(eos, 'eos')
  call trace%end(itimer)
END SUBROUTINE save_eos_variables

!===============================================================================
SUBROUTINE print_post_update_info (self)
  class(composition_patch_t):: self
  integer :: p(3)
  !-----------------------------------------------------------------------------
  if (self%verbose > 1) then
    p = self%pcheck
    memcheck(1) = self%mem(p(1),p(2),p(3),self%i_temp,1,1)
    memcheck(2) = self%mem(p(1),p(2),p(3),self%i_rho,1,1)
    memcheck(3) = self%mem(p(1),p(2),p(3),self%i_ne,1,1)
    call self%print_to_stdout(arr=memcheck(1:3))
  endif
END SUBROUTINE print_post_update_info

!===============================================================================
SUBROUTINE save_abundances_to_file (self)
  class(composition_patch_t), target :: self
  class(element_type), pointer :: element
  real(8) :: abnd
  integer :: lunit, i
  !-----------------------------------------------------------------------------
  if (self%ipos(3)==0) then
    open(newunit=lunit, file=trim(io%run_dir)//'abund.txt', form='FORMATTED')
    do i=1, composition%nel
      element => composition%elements(i)
      abnd = log10(element%abund) + 12.0
      abnd = max(abnd, -99.9)
      write(lunit, '(A3, F8.4)') element%id, abnd
    end do
    close(lunit)
  endif
END SUBROUTINE save_abundances_to_file

!===============================================================================
SUBROUTINE save_number_densities (self, comp, pre_init)
  class(composition_patch_t), target :: self
  class(composition_type) :: comp
  class(ion_type), pointer :: ion
  class(molec_type), pointer :: mol
  class(particle_type), pointer :: particle
  integer :: i, nspecies
  logical, intent(in) :: pre_init
  character(len=256) :: local_text
  character(len=6), dimension(:), allocatable :: species
  character(len=6) :: key
  !-----------------------------------------------------------------------------
  local_text = save_species
  call count_words(local_text, nspecies, sep=',')
  allocate(species(nspecies))
  call replace_char(local_text, ' ', '*')
  call replace_char(local_text, ',', ' ')
  read(local_text, *) species

  do i=1, nspecies
    key = species(i)
    call replace_char(key, '*', ' ')
    key = adjustl(key)
    call comp%particle_by_id(trim(key), particle, verbose=.true.)
    if (associated(particle)) then
      if (pre_init) then
        particle%in_opac = .true.
      else
        call self%save3(real(particle%n), trim(key))
      endif
    endif
  end do
END SUBROUTINE save_number_densities

!===============================================================================
SUBROUTINE composition_patch_nbor_relations(self, another, needed, needs_me, download)
  class(composition_patch_t) :: self
  class(task_t), pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  !-----------------------------------------------------------------------------
  call trace%begin ('composition_patch_t%nbor_relations')

  call self%m3d_patch_t%nbor_relations(another, needed, needs_me, download)

  if (needed .or. needs_me) then
    if (another%id /= self%patch%id) then
      select case (another%kind)
      case ('nlte_patch')
        ! is needed after convergence and needs me
      case ('rt_patch')
        needed = .false.
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
  endif
  call trace%end()
END SUBROUTINE composition_patch_nbor_relations

END MODULE composition_patch_mod