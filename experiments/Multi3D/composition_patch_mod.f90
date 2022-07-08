module composition_patch_mod
  USE io_mod
  USE mpi_mod
  USE omp_mod
  USE gpatch_mod
  USE m3d_patch_mod,   only: m3d_patch_t
  USE composition_mod, only: composition_type
  USE opac_patch_mod,  only: opacity_t
  USE rad_patch_mod,   only: radiation_t
  USE atom_patch_mod,  only: atom_patch_t
  USE rate_patch_mod,  only: rate_patch_t
  USE net_opac_patch_mod,  only: net_opac_t
  USE back_line_patch_mod, only:back_line_patch_t
  USE task_mod
  USE os_mod
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE ray_mod, only: set_carlson_quadratures
  USE list_mod, only: list_t
  USE bits_mod
  USE atom_mod, only: atom_type, atomic_spectrum_type
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  type, public, extends(m3d_patch_t) :: composition_patch_t
    class(gpatch_t), pointer :: mhd ! pointer to MHD patch
    type(composition_type) :: composition
  contains
    procedure :: pre_init => composition_pre_init
    procedure :: update => composition_update
    procedure :: post_update => composition_patch_post_update
    procedure :: finish => composition_finish
    procedure :: deallocate_mem => composition_patch_deallocate_mem
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type

contains

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  type(composition_type), save :: composition
  class(atom_patch_t), pointer :: atom_patch
  class(opacity_t), pointer :: opac
  class(radiation_t), pointer :: rad
  class(net_opac_t), pointer :: net_opac
  class(back_line_patch_t), pointer :: back_line
  class(atomic_spectrum_type), save, pointer :: output_sp
  integer :: iostat, i_out
  integer, save:: itimer=0, itimer2=0
  integer, save :: n_nu=1, verbose=0
  character(len=80), save :: absdat_file='./input_multi3d/absdat'
  character(len=80), save :: abund_file='./input_multi3d/abund'
  character(len=80), save :: absmet_file=''
  character(len=80), save :: linelist_file=''
  integer :: i, dir, p, m, n, m2, n2
  logical, save :: first_time=.true., absmet_big_end=.false.
  namelist /composition_params/ verbose, absdat_file, abund_file, absmet_file, absmet_big_end, linelist_file
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%composition_pre_init', itimer=itimer)
  self%kind = 'composition'
  call self%m3d_patch_t%pre_init()
  call self%append_to_task_list()
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
    call composition%init(absdat_file, abund_file, absmet_file, absmet_big_end=absmet_big_end)
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
  atom_patch%kind = 'atom_patch'
  atom_patch%composition => self%composition
  call atom_patch%pre_init_and_link(self, sp=output_sp)
  if (output_sp%nnu == 0) then
    output_sp => atom_patch%sp
  endif

  do n=1, self%n_nu
    do m=0, self%n_omega
      if (m==0) then
        allocate(opac)
        opac%kind = 'back_opac'
        opac%composition => self%composition
        call opac%composition%reassign_pointer()
        call opac%pre_init_and_link(self, i_mu=m, i_nu=n, sp=atom_patch%sp)
      else
        allocate(net_opac)
        net_opac%kind = 'net_opac'
        call net_opac%pre_init_and_link(self, i_mu=m, i_nu=n, sp=atom_patch%sp)
      endif

      allocate(rad)
      rad%kind = merge('angle_rad', 'mean_rad ', m>0)
      call rad%pre_init_and_link(self, i_mu=m, i_nu=n, sp=atom_patch%sp)
    end do
  end do

  do m=0, self%n_out
    m2 = m+self%n_omega
    do n=1, self%n_nu
      n2 = n+self%n_nu
      if (m==0) then
        allocate(opac)
        opac%kind = 'back_opac'
        opac%composition => self%composition
        call opac%composition%reassign_pointer()
        call opac%pre_init_and_link(self, i_mu=0, i_nu=n2, sp=output_sp)
      else
        if (self%composition%ll_initialised) then
          allocate(back_line)
          back_line%kind = 'back_line'
          back_line%composition => self%composition
          call back_line%pre_init_and_link(self, i_mu=m2, i_nu=n2, sp=output_sp)
        endif

        allocate(net_opac)
        net_opac%kind = 'net_opac'
        call net_opac%pre_init_and_link(self, i_mu=m2, i_nu=n2, sp=output_sp)

        allocate(rad)
        rad%kind = 'angle_rad'
        call rad%pre_init_and_link(self, i_mu=m2, i_nu=n2, sp=output_sp)
      endif
    end do

    if (m > 0 .and. self%boundaries%upper) then
      allocate(rad)
      rad%kind = 'out_rad'
      call rad%pre_init_and_link(self, i_mu=m2, sp=output_sp)
    endif
  end do
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
  call trace%begin('composition_patch_t%composition_update', itimer=itimer)
  call self%pre_update()

  if (self%finished) then
    return
  elseif (self%iter > 0) then
    print *, 'updating composition more than once is not implemented yet. Stopping...'
    stop
  endif

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

  call self%post_update()
  call trace%end(itimer)
END SUBROUTINE composition_update

!===============================================================================
SUBROUTINE composition_patch_post_update (self)
  class(composition_patch_t) :: self
  integer :: p(3), z
  logical :: skip
  !-----------------------------------------------------------------------------
  call trace%begin ('composition_patch_t%composition_patch_post_update')
  if (self%verbose > 1) then
    p = self%pcheck
    if (all(p > 0)) then
      write(stdout,'(A, 8(ES15.7))') self%m3d_patch_info(), &
      self%mem(p(1),p(2),p(3),:self%i_temp,1,:), self%composition%grph
    elseif (self%verbose > 2) then
      write(stdout,'(A,A)') self%m3d_patch_info()
    endif
  end if

  call self%m3d_patch_t%post_update()
  call trace%end()
END SUBROUTINE composition_patch_post_update

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
    if (another%id /= self%master_id) then
      select case (another%kind)
      case ('back_opac', 'back_line', 'atom_patch')
        needed = .false.
      case default
        needed = .false.; needs_me = .false.
      end select
    endif
  endif
  call trace%end()
END SUBROUTINE composition_patch_nbor_relations

!===============================================================================
SUBROUTINE composition_finish (self, rotated)
  class(composition_patch_t):: self
  logical, optional :: rotated
  !-----------------------------------------------------------------------------
  call self%m3d_patch_t%finish(rotated=rotated)
  call self%set(bits%frozen)
END SUBROUTINE composition_finish

!===============================================================================
SUBROUTINE composition_patch_deallocate_mem(self, deallocator)
  class(composition_patch_t) :: self
  class(m3d_patch_t), intent(in), target :: deallocator
  ! ----------------------------------------------------------------------------
  ! no memory to be deallocated
END SUBROUTINE composition_patch_deallocate_mem


end module composition_patch_mod