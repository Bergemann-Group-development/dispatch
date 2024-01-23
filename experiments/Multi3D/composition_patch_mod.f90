MODULE composition_patch_mod
  USE io_mod
  USE os_mod
  USE task_mod, only: task_t
  USE link_mod, only: nbor_t
  USE gpatch_mod
  USE handyfuncs_mod
  USE optab_mod, only: optab
  USE element_mod, only: element_type
  USE m3d_patch_mod,   only: m3d_patch_t, comoving_frame, ilambd, decouple_continuum, dlam, &
  use_density, use_ne, overlapping_bb, overlapping_bf, line_blanketing, binned_opac, &
  save_Qrad, maxiter, rotate_atmos, long_vec
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
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type
  
  real(8), public :: FeH=0
  logical, public :: save_atmos=.true.
  type(composition_type), public, save :: composition

contains

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  class(atom_patch_t), pointer :: atom_patch
  class(atomic_spectrum_type), save, pointer :: output_sp
  integer :: i, iostat, n_nu, nlists, mu, nw
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
  absmet_big_end, isotope_file, molec_file, H_occupation_probability, opac_table
  namelist /linelist_params/ dlam, exclude_elements, linelist_file, line_lists
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

  ! copy composition with all the information from the abund
  ! and absdat file into dispackage to make it threadprivate
  allocate(self%composition)
  self%composition = composition
  ! this is needed, otherwise all pointers will still point
  ! to ions and molecules in multi3d%composition
  call self%composition%reassign_pointer()

  allocate(atom_patch)
  atom_patch%composition => self%composition
  call atom_patch%pre_init_atom()
  call atom_patch%pre_init_spectrum(output_sp)

  n_nu = self%n_nu
  if (associated(atom_patch%atom)) then
    call atom_patch%pre_init_and_link(self)
    self%atom_patch => atom_patch

    if (ilambd + maxiter > 0) then
      do mu=1,2
        if (decouple_continuum) then
          call self%create_sub_patches([1, n_nu],      mu, atom_patch%bb_sp, atom_patch%bb_atom)
          call self%create_sub_patches([1, n_nu]+n_nu, mu, atom_patch%bf_sp, atom_patch%bf_atom)
        else
          call self%create_sub_patches([1, n_nu], mu, atom_patch%sp, atom_patch%atom)
        endif
      end do
    endif
    
    call self%create_sub_patches([-n_nu, -1], -1, output_sp, atom_patch%output_atom)
  else
    if (.not. external_sp) then
      stop 'No atom nor output spectrum specified! Stopping...'
    endif
    deallocate(atom_patch)
    call self%create_sub_patches([-n_nu, -1], -1, output_sp)
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

    if (rotate_atmos) then
      allocate(self%rotated_composition(self%n_out))
    endif

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

      if (i==1) then
        allocate(temp8, source=dble(temp))
        allocate(rho8, source=dble(rho))
        allocate(ne8, source=dble(ne))
  
        ! loading temp, ne and rho in to composition objekt
        ! ion and molecule equilibrium are solved immediately
        call self%composition%load_state_variables(temp8, ne8, rho8, use_rho=use_density, use_ne=use_ne)
    
        if (save_atmos) then
          nu = 1d8 * cgs%c / 5000.d0
          call self%composition%absorption(nu, chi, sorted=.true.)
          call self%save4(authority, 'atmos')
        endif
      else
        temp8 = authority(:,:,:,self%i_temp)
        rho8  = authority(:,:,:,self%i_rho)
        ne8   = authority(:,:,:,self%i_ne)
        associate(comp => self%rotated_composition(i-1))
          comp = self%composition
          call comp%reassign_pointer()
          call comp%load_state_variables(temp8, ne8, rho8, use_rho=use_density, use_ne=use_ne)
        end associate
      endif
    end do


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
SUBROUTINE print_post_update_info (self)
  class(composition_patch_t):: self
  integer :: p(3)
  !-----------------------------------------------------------------------------
  if (self%verbose > 1) then
    p = self%pcheck
    call self%print_to_stdout(arr=dble(self%mem(p(1),p(2),p(3),:self%i_temp,1,1)))
  endif
END SUBROUTINE print_post_update_info

!===============================================================================
SUBROUTINE save_abundances_to_file (self)
  class(composition_patch_t), target :: self
  class(element_type), pointer :: element
  integer :: lunit, i
  !-----------------------------------------------------------------------------
  if (self%ipos(3)==0) then
    open(newunit=lunit, file=trim(io%run_dir)//'abund.txt', form='FORMATTED')
    do i=1, self%composition%nel
      element => self%composition%elements(i)
      write(lunit, '(A3, F8.4)') element%id, log10(element%abund) + 12.0
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