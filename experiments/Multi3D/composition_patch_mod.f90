MODULE composition_patch_mod
  USE io_mod
  USE os_mod
  USE mpi_mod
  USE omp_mod
  USE gpatch_mod
  USE handyfuncs_mod
  USE element_mod, only: element_type
  USE m3d_patch_mod,   only: m3d_patch_t, comoving_frame, ilambd, decouple_continuum, &
  use_density, overlapping_bb, overlapping_bf, line_blanketing
  USE composition_mod, only: composition_type, dlam, exclude_elements
  USE back_cont_patch_mod,  only: back_cont_patch_t
  USE angle_rad_patch_mod,   only: rad_patch_t
  USE atom_patch_mod,  only: atom_patch_t
  USE atom_rate_patch_mod,  only: atom_rate_t
  USE net_opac_patch_mod,  only: net_opac_t
  USE back_line_patch_mod, only:back_line_patch_t
  USE out_rad_patch_mod, only:out_rad_patch_t
  USE proj_patch_mod, only:proj_patch_t
  USE mean_rad_patch_mod, only:mean_rad_patch_t
  USE task_mod
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
    procedure :: update => composition_update
    procedure :: create_sub_patches
    procedure :: print_post_update_info
    procedure :: save_abundances_to_file
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type
  
  real(8), public :: FeH=0

contains

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  type(composition_type), save :: composition
  class(atom_patch_t), pointer :: atom_patch
  class(atomic_spectrum_type), save, pointer :: output_sp
  integer :: i, iostat, n_nu, n_mu, nlists
  integer, save:: itimer=0, itimer2=0
  integer, save :: verbose=0
  character(len=128), save :: absdat_file='./input_multi3d/absdat'
  character(len=128), save :: abund_file='./input_multi3d/abund'
  character(len=128), save :: isotope_file=''!'./input_multi3d/isotopes_abund.txt'
  character(len=128), save :: molec_file='./input_multi3d/molec_data.txt'
  character(len=128), save :: absmet_file=''
  character(len=512), save :: linelist_file=''
  character(len=128), save :: line_lists(10)=''
  logical, save :: first_time=.true., absmet_big_end=.false., external_sp
  namelist /composition_params/ verbose, absdat_file, abund_file, absmet_file, &
  absmet_big_end, isotope_file, molec_file
  namelist /linelist_params/ dlam, exclude_elements, linelist_file, line_lists
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%composition_pre_init', itimer=itimer)
  self%kind = 'composition'
  call self%m3d_patch_t%pre_init()
  call self%set_nv_and_nw(9,1)
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

    if (verbose) write (*,*) ' - Initializing the NEW opacity package'
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

    call trace%end(itimer2)
    
    allocate(output_sp)
    call output_sp%init()
    output_sp%name = 'output'
    external_sp = output_sp%nnu > 0
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
  call atom_patch%pre_init_and_link(self)

  n_nu = self%n_nu
  n_mu = self%n_mu

  if (atom_patch%atom_allocated) then
    self%atom_patch => atom_patch
    if (decouple_continuum) then
      call self%create_sub_patches([1, n_nu],      [0, n_mu], self%composition, atom_patch%bb_sp)
      call self%create_sub_patches([1, n_nu]+n_nu, [0, n_mu], self%composition, atom_patch%bf_sp)
    else
      call self%create_sub_patches([1, n_nu], [0, n_mu], self%composition, atom_patch%sp)
    endif
  else
    deallocate(atom_patch)
    if (.not. external_sp) then
      stop 'No atom nor output spectrum specified! Stopping...'
    endif
  endif

  call self%create_sub_patches([-n_nu, 0], [-self%n_out, 0], self%composition, output_sp)
  
  call self%save_abundances_to_file()
  call self%post_pre_init()
  call trace%end(itimer)
END SUBROUTINE composition_pre_init

!===============================================================================
SUBROUTINE create_sub_patches (self, nu, mu, composition, sp)
  class(composition_patch_t) :: self
  integer, intent(in) :: nu(2), mu(2)
  class(composition_type), intent(in), target :: composition
  class(atomic_spectrum_type), intent(in), pointer :: sp
  class(back_cont_patch_t), pointer :: cont
  class(back_line_patch_t), pointer :: line
  class(out_rad_patch_t), pointer :: out
  class(net_opac_t), pointer :: net
  class(proj_patch_t), pointer :: proj
  class(rad_patch_t), pointer :: rad
  class(mean_rad_patch_t), pointer :: mean
  integer :: n, m
  logical :: cf, new_cont, new_line, new_rad, new_out, new_net, new_proj, new_mean
  !-----------------------------------------------------------------------------
  if (.not. associated(sp)) return
  if (sp%nnu == 0) return
  do n=nu(1), nu(2)
    cf = comoving_frame .and. n>0
    nullify(mean)
  do m=mu(1), mu(2)
    new_mean = n>0  .and. m==0
    new_cont = n/=0 .and. m==0
    new_line = n/=0 .and. (m<0 .or. merge(m==0, m/=0, cf)) .and. composition%ll_initialised
    new_net  = n/=0 .and. (m<0 .or. merge(m==0, m/=0, cf))
    new_proj = m/=0 .and. cf
    new_rad  = n/=0 .and. m/=0
    new_out  = n==0 .and. m<0 .and. self%patch%boundaries%is_set(bits%zu)

    if (new_mean) then
      allocate(mean)
      call mean%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp)
    endif
    if (new_cont) then
      allocate(cont)
      cont%composition => composition
      if (associated(mean)) then
        cont%mean_rad => mean
        mean%back_cont => cont
      endif
      call cont%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp)
    endif
    if (new_line) then
      allocate(line)
      line%composition => composition
      call line%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp)
    endif
    if (new_net) then
      allocate(net)
      call net%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp)
      if (associated(self%atom_patch)) then
        call self%atom_patch%add_to_key_list(net)
      endif
      if (associated(mean)) then
        mean%opac => net
      endif
    endif
    if (new_proj) then
      allocate(proj)
      call proj%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp)
    endif
    if (new_rad) then
      allocate(rad)
      if (associated(mean)) rad%mean_rad => mean
      call rad%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp)
    endif
    if (new_out) then
      allocate(out)
      call out%pre_init_and_link(self, i_mu=m, i_nu=n, sp=sp)
    endif
  end do
end do
END SUBROUTINE create_sub_patches

!===============================================================================
SUBROUTINE composition_update (self)
  class(composition_patch_t) :: self
  real(4), dimension(:,:,:), pointer :: d, e, pe
  real(4), dimension(:,:,:), pointer, contiguous :: ne, rho, temp
  real(4), dimension(:,:,:,:), pointer :: nh
  integer, save:: itimer=0
  integer :: li(3), ui(3)
  ! ----------------------------------------------------------------------------
  if (self%is_finished()) return
  call trace%begin('composition_patch_t%composition_update', itimer=itimer)

  if (self%iter == 0) then
    li = self%patch%li
    ui = self%patch%ui
  
    e  => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%e,self%it,1)
    d  => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),self%patch%idx%d,self%it,1)
    pe => self%patch%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),6               ,self%it,1)
  
    ne   => self%mem(:,:,:,self%i_ne,  self%it,1)
    rho  => self%mem(:,:,:,self%i_rho, self%it,1)
    temp => self%mem(:,:,:,self%i_temp,self%it,1)
    nh   => self%mem(:,:,:,self%i_nh1:self%i_nh6,self%it,1)
  
    temp = e
    rho  = d * scaling%d
    ne   = pe * scaling%p / (cgs%k_b * temp)
  
    ! loading temp, ne and rho in to composition objekt
    ! ion and molecule equilibrium are solved immediately
    call self%composition%load_state_variables(temp, ne, rho, use_rho=use_density)
    call self%composition%get_nhyd(nh)

    self%dtime = max(self%dtime, io%end_time - self%dtime)
  else
    self%dtime = self%patch%dtime
    !print *, self%m3d_patch_info()
    !print *, 'updating composition more than once is not implemented yet. Stopping...'
    !stop
  endif
  call trace%end(itimer)
END SUBROUTINE composition_update

!===============================================================================
SUBROUTINE print_post_update_info (self)
  class(composition_patch_t):: self
  integer :: p(3)
  !-----------------------------------------------------------------------------
  p = self%pcheck
  call self%print_to_stdout(arr=self%mem(p(1),p(2),p(3),:self%i_temp,1,1))
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
      case ('back_cont', 'back_line')
        needed = .false.
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
  endif
  call trace%end()
END SUBROUTINE composition_patch_nbor_relations

END MODULE composition_patch_mod