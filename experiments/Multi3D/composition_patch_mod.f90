MODULE composition_patch_mod
  USE io_mod
  USE os_mod
  USE task_mod, only: task_t
  USE link_mod, only: nbor_t
  USE gpatch_mod
  USE handyfuncs_mod
  USE optab_mod, only: optab
  USE element_mod, only: element_type
  USE m3d_patch_mod,   only: m3d_patch_t, composition, comoving_frame, ilambd, dlam, &
  use_density, use_ne, overlapping_bb, overlapping_bf, line_blanketing, binned_opac, &
  save_Qrad, maxiter, rotate_atmos, long_vec, sample_back_opac, decouple_continuum, &
  linelist_in_output, sample_back_opac_in_output, make_eos
  USE m3d_debug_mod
  USE composition_mod, only: composition_type, exclude_elements
  USE atom_patch_mod,  only: atom_patch_t
  USE nlte_patch_mod, only:nlte_patch_t
  USE bftrans_mod, only: H_occupation_probability
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE ray_mod, only: set_carlson_quadratures
  USE list_mod, only: list_t
  USE link_mod, only: nbor_t
  USE bits_mod
  USE atom_mod, only: atom_type, atomic_spectrum_type
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  type, public, extends(m3d_patch_t) :: composition_patch_t
    class(gpatch_t), pointer :: mhd ! pointer to MHD patch
    class(atom_patch_t), pointer :: atom_patch
  contains
    procedure :: pre_init => composition_pre_init
    procedure :: execute
    procedure :: create_sub_patches
    procedure :: print_post_update_info
    procedure :: save_abundances_to_file
    procedure :: init_tabulated_comp
    procedure :: init_eos_comp
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type
  
  real(8), public :: FeH=0, dtemp=250., ldne=0.25, ldrho=0.25, ldtemp=0.1
  logical, public :: save_atmos=.true.

contains

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  class(atom_patch_t), pointer :: atom_patch
  class(atomic_spectrum_type), save, pointer :: output_sp
  integer :: i, iostat, n_nu, nlists, mu, nw, n_mu
  integer, save:: itimer=0, itimer2=0
  integer, save :: verbose=0
  character(len=256), save :: absdat_file='./input_multi3d/TS_absdat.dat'
  character(len=256), save :: abund_file='./input_multi3d/abund_magg'
  character(len=256), save :: isotope_file=''!'./input_multi3d/isotopes_abund.txt'
  character(len=256), save :: molec_file='./input_multi3d/molec_data.txt'
  character(len=256), save :: opac_table=''
  character(len=256), save :: absmet_file=''
  character(len=512), save :: linelist_file=''
  character(len=256), save :: line_lists(10)=''
  logical, save :: first_time=.true., absmet_big_end=.false., external_sp
  namelist /composition_params/ verbose, absdat_file, abund_file, absmet_file, &
  absmet_big_end, isotope_file, molec_file, H_occupation_probability, opac_table, &
  sample_back_opac, sample_back_opac_in_output, &
  dtemp, ldne, ldrho, ldtemp
  namelist /linelist_params/ dlam, exclude_elements, linelist_file, line_lists, &
  linelist_in_output
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%composition_pre_init', itimer=itimer)
  self%kind = 'composition'
  call self%m3d_patch_t%pre_init()
  nw = merge(1+size(long_vec,2), 1, rotate_atmos)
  call self%set_nv_and_nw(7,nw)
  self%static_mem=.true.

  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, composition_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, composition_params)
    end if
    call io%namelist_warning ('composition_params', iostat)
    write (io%namelist, composition_params)

    rewind (io%input)
    read (io%input, linelist_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, linelist_params)
    end if
    call io%namelist_warning ('linelist_params', iostat)
    write (io%namelist, linelist_params)

    if (verbose > 0) write (*,*) ' - Initializing the NEW opacity package'
    composition%verbose = verbose

    line_blanketing = absmet_file /= ''
    call trace%begin('compositiont%init', itimer=itimer2)
    call composition%init(absdat_file, abund_file, absmet_file, FeH, &
    isotope_file=isotope_file, molec_file=molec_file, absmet_big_end=absmet_big_end)

    if (linelist_file /= '') then
      call composition%init_linelist(linelist_file)
    else
      i = 1
      do while (line_lists(i) /= '' .and. i <= size(line_lists))
        call composition%init_linelist(line_lists(i))
        i = i + 1
      end do
    endif
    
    allocate(output_sp)
    binned_opac = opac_table /= ''
      
    if (binned_opac) then
      call optab%init(opac_table, output_sp)
    else
      call output_sp%init()
      output_sp%name = 'output'
    endif

    external_sp = output_sp%nnu > 0
    call trace%end(itimer2)
  endif

  allocate(atom_patch)
  call atom_patch%pre_init_atom()
  call atom_patch%pre_init_spectrum(output_sp)

  n_nu = self%n_nu
  n_mu = size(long_vec,2)
  if (associated(atom_patch%atom)) then
    call atom_patch%pre_init_and_link(self)
    self%atom_patch => atom_patch

    if (maxiter > 0) then
      do mu=1,2
        if (decouple_continuum) then
          call self%create_sub_patches([1, n_nu],      mu, atom_patch%bb_sp, atom_patch%bb_atom)
          call self%create_sub_patches([1, n_nu]+n_nu, mu, atom_patch%bf_sp, atom_patch%bf_atom)
        else
          call self%create_sub_patches([1, n_nu], mu, atom_patch%sp, atom_patch%atom)
        endif
      end do
    endif
    
    if (n_mu>0) then
      if (ilambd > 0) then
        do mu=1,2
          call self%create_sub_patches([-n_nu, -1], mu, output_sp, atom_patch%output_atom)
        end do
      endif
      call self%create_sub_patches([-n_nu, -1], -1, output_sp, atom_patch%output_atom)
    endif
  else
    deallocate(atom_patch)
    if (external_sp .and. n_mu>0) then
      if (ilambd > 0) then
        do mu=1,2
          call self%create_sub_patches([-n_nu, -1], mu, output_sp)
        end do
      endif
      call self%create_sub_patches([-n_nu, -1], -1, output_sp)
      !stop 'No atom nor output spectrum specified! Stopping...'
    endif
  endif

  
  call self%save_abundances_to_file()
  call self%post_pre_init()
  call trace%end(itimer)
END SUBROUTINE composition_pre_init

!===============================================================================
SUBROUTINE create_sub_patches (self, nu, mu, sp, atom)
  class(composition_patch_t) :: self
  integer, intent(in) :: nu(2), mu
  class(atomic_spectrum_type), intent(in), pointer :: sp
  type(atom_type), intent(in), pointer, optional :: atom
  class(nlte_patch_t), pointer :: nlte
  integer :: n
  logical :: new_nlte
  !-----------------------------------------------------------------------------
  if (.not. associated(sp)) return
  if (sp%nnu == 0) return
  do n=nu(1), nu(2)
    allocate(nlte)
    if (present(atom)) nlte%atom => atom
    call nlte%pre_init_and_link(self, i_mu=mu, i_nu=n, sp=sp)
  end do
END SUBROUTINE create_sub_patches

!===============================================================================
SUBROUTINE execute (self)
  class(composition_patch_t) :: self
  class(composition_type), pointer :: comp
  real(4), dimension(:,:,:), pointer :: d, e, pe
  real(4), dimension(:,:,:,:), pointer, contiguous :: authority, chi
  real(4), dimension(:,:,:), pointer, contiguous :: ne, rho, temp, vx, vy, vz, px, py, pz
  real(8), dimension(:,:,:), allocatable :: ne8, rho8, temp8
  real(8) :: vec(3), nu(1)
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

    do i=1, self%nw
      authority => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,self%it,i)
      e  => authority(:,:,:,self%patch%idx%e)
      d  => authority(:,:,:,self%patch%idx%d)
      px => authority(:,:,:,self%patch%idx%px)
      py => authority(:,:,:,self%patch%idx%py)
      pz => authority(:,:,:,self%patch%idx%pz)
      pe => authority(:,:,:,6)
    
      authority => self%mem(:,:,:,:,self%it,i)
      ne   => authority(:,:,:,self%i_ne)
      rho  => authority(:,:,:,self%i_rho)
      temp => authority(:,:,:,self%i_temp)
      vx   => authority(:,:,:,self%i_vx)
      vy   => authority(:,:,:,self%i_vy)
      vz   => authority(:,:,:,self%i_vz)
      chi(1:1,1:n(1),1:n(2),1:n(3)) => authority(:,:,:,7)
  
      temp = e
      rho  = d
      ne   = pe / (cgs%k_b * temp)
      vx   = px / (cgs%c * rho)
      vy   = py / (cgs%c * rho)
      vz   = pz / (cgs%c * rho)

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
      call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_density, use_ne=use_ne)

      if (.not. use_ne) ne = real(comp%ne)
      if (.not. use_density) rho = real(comp%rho)
      if (save_atmos) then
        nu = 1d8 * cgs%c / 5000.d0
        call comp%absorption(nu, chi, sorted=.true.)
      endif
    end do
    
    if (save_atmos) then
      call self%save5(self%mem(:,:,:,:,self%it,:), 'atmos')
    endif

    if (sample_back_opac) then
      call self%init_tabulated_comp(self%composition%temp, self%composition%ne, dtemp=dtemp, ldne=ldne)
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
    real(8), dimension(:,:,:), intent(in) :: temp, ne
    real(8), intent(in) :: dtemp, ldne
    integer :: i, j, ntemp, nne, nnu
    real(8), dimension(:,:,:), allocatable :: tab_temp, tab_ne, tab_rho
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
    do j=1, nne
    do i=1, ntemp
      tab_temp(i,j,1) = mint + (i-1)*dtemp
      tab_ne  (i,j,1) = mine + (j-1)*ldne
    end do
    end do
    tab_ne = exp(tab_ne)
    call self%tab_comp%load_state_variables(tab_temp, tab_ne, tab_rho, use_rho=.false., use_ne=.true.)
END SUBROUTINE init_tabulated_comp

!===============================================================================
SUBROUTINE init_eos_comp (self, temp, rho, ldtemp, ldrho)
  class(composition_patch_t), target :: self
  real(8), dimension(:,:,:), intent(in) :: temp, rho
  real(8), intent(in) :: ldtemp, ldrho
  real(8) :: minTemp, maxTemp, minRho, maxRho
  integer :: i, j, ntemp, nrho, nnu
  real(8), dimension(:,:,:), allocatable :: tab_temp, tab_ne, tab_rho
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

  do j=1, nrho
  do i=1, ntemp
    tab_temp(i,j,1) = minLnT + (i-1)*ldtemp
    tab_rho  (i,j,1) = minLnRho + (j-1)*ldrho
  end do
  end do

  tab_rho = exp(tab_rho)
  tab_temp = exp(tab_temp)

  call self%tab_comp%load_state_variables(tab_temp, tab_ne, tab_rho, use_rho=.true., use_ne=.false.)
END SUBROUTINE init_eos_comp

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
      case ('atom_patch')
        ! is needed after convergence and needs me
      case ('nlte_patch')
        needed = .false.
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
  endif
  call trace%end()
END SUBROUTINE composition_patch_nbor_relations

END MODULE composition_patch_mod