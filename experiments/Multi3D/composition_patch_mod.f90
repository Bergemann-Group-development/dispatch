module composition_patch_mod
  USE io_mod
  USE os_mod
  USE mpi_mod
  USE omp_mod
  USE gpatch_mod
  USE element_mod, only: element_type
  USE m3d_patch_mod,   only: m3d_patch_t, hash_string
  USE composition_mod, only: composition_type
  USE back_opac_patch_mod,  only: back_opac_patch_t
  USE rad_patch_mod,   only: rad_patch_t
  USE atom_patch_mod,  only: atom_patch_t
  USE rate_patch_mod,  only: rate_patch_t
  USE net_opac_patch_mod,  only: net_opac_t
  USE back_line_patch_mod, only:back_line_patch_t
  USE out_rad_patch_mod, only:out_rad_patch_t
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
    type(composition_type) :: composition
    real(8) :: FeH
  contains
    procedure :: pre_init => composition_pre_init
    procedure :: update => composition_update
    procedure :: print_post_update_info
    procedure :: save_abundances_to_file
    procedure :: nbor_relations => composition_patch_nbor_relations
    procedure :: init_nbors_hash
  end type

contains

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  type(composition_type), save :: composition
  class(atom_patch_t), pointer :: atom_patch
  class(back_opac_patch_t), pointer :: back_opac
  class(rad_patch_t), pointer :: rad
  class(out_rad_patch_t), pointer :: out_rad
  class(net_opac_t), pointer :: net_opac
  class(back_line_patch_t), pointer :: back_line
  class(atomic_spectrum_type), save, pointer :: output_sp
  integer :: iostat, i_out
  integer, save:: itimer=0, itimer2=0
  integer, save :: n_nu=1, verbose=0
  character(len=128), save :: absdat_file='./input_multi3d/absdat'
  character(len=128), save :: abund_file='./input_multi3d/abund'
  character(len=128), save :: isotope_file='./input_multi3d/isotopes_abund.txt'
  character(len=128), save :: molec_file='./input_multi3d/molec_data.txt'
  character(len=128), save :: absmet_file=''
  character(len=128), save :: linelist_file=''
  integer :: i, dir, p, m, n, m2, n2
  logical, save :: first_time=.true., absmet_big_end=.false.
  namelist /composition_params/ verbose, absdat_file, abund_file, absmet_file, &
  absmet_big_end, linelist_file, isotope_file, molec_file
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%composition_pre_init', itimer=itimer)
  self%kind = 'composition'
  call self%m3d_patch_t%pre_init()
  call self%set_nv_and_nw(9,1)
  self%static_mem=.true.

  if (first_time) then
    rewind (io%input)
    read (io%input, composition_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, composition_params)
    end if
    call io%namelist_warning ('composition_params', iostat)
    write (io%namelist, composition_params)

    if (verbose) write (*,*) ' - Initializing the NEW opacity package'
    composition%verbose = verbose

    call trace%begin('compositiont%init', itimer=itimer2)
    call composition%init(absdat_file, abund_file, absmet_file, self%FeH, &
    isotope_file=isotope_file, molec_file=molec_file, absmet_big_end=absmet_big_end)
    call composition%init_linelist(linelist_file)
    call trace%end(itimer2)
    
    allocate(output_sp)
    call output_sp%init()
    first_time = .false.
  endif

  ! copy composition with all the information from the abund
  ! and absdat file into dispackage to make it threadprivate
  self%composition = composition
  ! this is needed, otherwise all pointers will still point
  ! to ions and molecules in multi3d%composition
  call self%composition%reassign_pointer()

  allocate(atom_patch)
  atom_patch%composition => self%composition
  call atom_patch%pre_init_and_link(self, sp=output_sp)
  if (atom_patch%atom_allocated) then
    if (output_sp%nnu == 0) then
      output_sp => atom_patch%sp
    endif

    do n=1, self%n_nu
      allocate(back_opac)
      back_opac%composition => self%composition
      call back_opac%composition%reassign_pointer()
      call back_opac%pre_init_and_link(self, i_mu=0, i_nu=n, sp=atom_patch%sp)
    end do
  
    do n=1, self%n_nu
      do m=1, self%n_omega
        allocate(net_opac)
        call net_opac%pre_init_and_link(self, i_mu=m, i_nu=n, sp=atom_patch%sp)
      end do
    end do
  
    do n=1, self%n_nu
      do m=0, self%n_omega
        allocate(rad)
        call rad%pre_init_and_link(self, i_mu=m, i_nu=n, sp=atom_patch%sp)
      end do
    end do
  else
    deallocate(atom_patch)
    if (output_sp%nnu == 0) then
      stop 'No atom nor output spectrum specified! Stopping...'
    endif
  endif

  do m=0, self%n_out
    m2 = m+self%n_omega
    do n=1, self%n_nu
      n2 = n+self%n_nu
      if (m==0) then
        allocate(back_opac)
        back_opac%composition => self%composition
        call back_opac%composition%reassign_pointer()
        call back_opac%pre_init_and_link(self, i_mu=0, i_nu=n2, sp=output_sp)
      else
        if (self%composition%ll_initialised) then
          allocate(back_line)
          back_line%composition => self%composition
          call back_line%pre_init_and_link(self, i_mu=m2, i_nu=n2, sp=output_sp)
        endif

        allocate(net_opac)
        call net_opac%pre_init_and_link(self, i_mu=m2, i_nu=n2, sp=output_sp)

        allocate(rad)
        call rad%pre_init_and_link(self, i_mu=m2, i_nu=n2, sp=output_sp)
      endif
    end do

    if (m > 0 .and. self%boundaries%upper) then
      allocate(out_rad)
      call out_rad%pre_init_and_link(self, i_mu=m2, sp=output_sp)
    endif
  end do
  call self%save_abundances_to_file(trim(io%run_dir) // 'abund.txt')
  call self%post_pre_init()
  call trace%end(itimer)
END SUBROUTINE composition_pre_init

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
    call self%composition%load_state_variables(temp, ne, rho)
    call self%composition%get_nhyd(nh)

    self%dtime = io%end_time - self%dtime
  else
    self%dtime = self%patch%dtime
    !print *, self%m3d_patch_info()
    !print *, 'updating composition more than once is not implemented yet. Stopping...'
    !stop
  endif

  call self%post_update()
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
SUBROUTINE save_abundances_to_file (self, file)
  class(composition_patch_t), target :: self
  class(element_type), pointer :: element
  character(len=*) :: file
  integer :: lunit, i
  !-----------------------------------------------------------------------------
  open(newunit=lunit, file=file, form='FORMATTED')
  do i=1, self%composition%nel
    element => self%composition%elements(i)
    write(lunit, '(A3, F8.4)') element%id, log10(element%abund) + 12.0
  end do
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
      case ('back_opac', 'back_line')
        needed = .false.
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
  endif
  call trace%end()
END SUBROUTINE composition_patch_nbor_relations

!===============================================================================
!> Fast method to generate nbor lists in cases with very long task lists
!===============================================================================
SUBROUTINE init_nbors_hash (self, nbors)
  class(composition_patch_t):: self
  class(nbor_t), pointer:: nbors
  !.............................................................................
  class(task_t), pointer:: task
  class(nbor_t), pointer:: nbor
  integer, save:: itimer=0
  integer :: i, i_nu, i_mu, key(4)
  !-----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%init_nbors_hash', itimer=itimer)
  call self%m3d_patch_t%init_nbors_hash(nbors)

  key(1:3) = self%key(1:3)
  key(4) = hash_string('atom_patch')
  call self%check_nbor(key, nbors)

  do i_nu=1, 2*self%n_nu
    key(2) = i_nu
    key(4) = hash_string('back_opac')
    call self%check_nbor(key, nbors)
  end do

  do i_nu=self%n_nu+1, 2*self%n_nu
    key(2) = i_nu
    do i_mu=1, self%n_out
      key(3) = i_mu+self%n_omega
      key(4) = hash_string('back_line')
      call self%check_nbor(key, nbors)
    end do
  end do
  call trace%end (itimer)
END SUBROUTINE init_nbors_hash

end module composition_patch_mod