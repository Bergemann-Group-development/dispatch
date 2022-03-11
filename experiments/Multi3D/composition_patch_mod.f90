module composition_patch_mod
  USE io_mod
  USE gpatch_mod
  USE m3d_patch_mod,   only: m3d_patch_t
  USE composition_mod, only: composition_type
  USE opac_patch_mod,  only: opacity_t
  USE rad_patch_mod,   only: radiation_t
  USE atom_patch_mod,  only: atom_patch_t
  USE rate_patch_mod,  only: rate_patch_t
  USE task_mod
  USE os_mod
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE ray_mod, only: set_carlson_quadratures
  USE list_mod, only: list_t
  USE bits_mod
  implicit none
  private
  
  !-----------------------------------------------------------------------------
  type, public, extends(m3d_patch_t) :: composition_patch_t
    class(gpatch_t), pointer :: mhd ! pointer to MHD patch
    type(composition_type) :: composition
    class(opacity_t), pointer, dimension(:) :: back_opac
    class(opacity_t), pointer, dimension(:,:) :: opac
    class(radiation_t), pointer, dimension(:) :: out_rad
    class(radiation_t), pointer, dimension(:,:) :: rad
    class(atom_patch_t), pointer :: atom_patch
  contains
    procedure :: pre_init => composition_pre_init
    procedure :: update => composition_update
    procedure :: finish => composition_finish
    procedure :: add_to_task_list
    procedure :: nbor_relations => composition_patch_nbor_relations
  end type

contains

!===============================================================================
SUBROUTINE composition_pre_init (self)
  class(composition_patch_t), target:: self
  type(composition_type), save :: composition
  integer :: iostat, i_out
  integer, save:: itimer=0
  integer, save :: n_nu=1, verbose=0
  character(len=80), save :: absdat_file='./input_multi3d/absdat'
  character(len=80), save :: abund_file='./input_multi3d/abund'
  integer :: i, dir, p, m, n
  logical, save :: first_time=.true.
  namelist /composition_params/ verbose, absdat_file, abund_file
  ! ----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%composition_pre_init', itimer=itimer)
  self%kind = 'composition'
  call self%m3d_patch_t%pre_init()
  call self%set_nv_and_nw(9,1)

  if (first_time) then
    rewind (io%input)
    read (io%input, composition_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, composition_params)
    end if
    call io%namelist_warning ('composition_params', iostat)
    write (io%namelist, composition_params)

    write (*,*) ' - Initializing the NEW opacity package'
    call composition%init(absdat_file, abund_file)

    first_time = .false.
  endif

  self%verbose = verbose

  ! copy composition with all the information from the abund
  ! and absdat file into dispackage to make it threadprivate
  self%composition = composition
  ! this is needed, otherwise all pointers will still point
  ! to ions and molecules in multi3d%composition
  call self%composition%reassign_pointer()

  allocate(self%atom_patch)
  self%atom_patch%kind = 'atom_patch'
  call self%atom_patch%pre_init_and_link(self)
  self%sp => self%atom_patch%sp
  self%atom_patch%composition => self%composition

  allocate(self%back_opac(self%n_nu))
  allocate(self%opac(self%n_mu,  self%n_nu))
  allocate(self%rad(0:self%n_mu, self%n_nu))

  do n=1, self%n_nu
    self%back_opac(n)%kind = 'back_opac'
    self%back_opac(n)%composition => self%composition
    self%back_opac(n)%mean_rad => self%rad(0,n)
    call self%back_opac(n)%composition%reassign_pointer()
    call self%back_opac(n)%pre_init_and_link(self, i_nu=n)

    self%rad(0,n)%kind = 'mean_rad'
    call self%rad(0,n)%pre_init_and_link(self, i_nu=n)
  end do

  do m=1, self%n_mu
    do n=1, self%n_nu
      self%opac(m,n)%kind = 'net_opac'
      self%opac(m,n)%rad  => self%rad(m,n)
      call self%opac(m,n)%pre_init_and_link(self, i_mu=m, i_nu=n)

      self%rad(m,n)%kind = 'angle_rad'
      self%rad(m,n)%opac => self%opac(m,n)
      self%rad(m,n)%mean_rad => self%rad(0, n)
      call self%rad(m,n)%pre_init_and_link(self, i_mu=m, i_nu=n)
    end do
  end do

  if (self%boundaries%upper) then
    allocate(self%out_rad(self%n_out))
    do i=1, self%n_out
      self%out_rad(i)%kind = 'out_rad'
      self%out_rad(i)%sheep => self%rad(i+self%n_omega, :)
      do n=1, self%n_nu
        self%out_rad(i)%sheep(n)%out_rad => self%out_rad(i)
      end do
      call self%out_rad(i)%pre_init_and_link(self, i_mu=i+self%n_omega)
    end do
  endif
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
SUBROUTINE add_to_task_list(self, task_list)
  class(composition_patch_t) :: self
  class(list_t), pointer:: task_list
  integer :: m,n
  !-----------------------------------------------------------------------------
  call trace%begin('composition_patch_t%add_to_task_list')
  call self%append_to_task_list(task_list)
  call self%atom_patch%add_to_task_list(task_list)
  
  do n=1, self%n_nu
    call self%back_opac(n)%append_to_task_list(task_list)
    call self%rad(0,n)%append_to_task_list(task_list)
  end do

  do m=1, self%n_mu
    do n=1, self%n_nu
      call self%opac(m,n)%append_to_task_list(task_list)
      call self%rad(m,n)%append_to_task_list(task_list)
    end do
  end do

  if (self%boundaries%upper) then
    do m=1, self%n_out
      call self%out_rad(m)%append_to_task_list(task_list)
    end do
  endif
  call trace%end()
END SUBROUTINE add_to_task_list

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
      if (another%kind=='back_opac') then
        needed = .false.
      else
        needed = .false.; needs_me = .false.
      endif
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

end module composition_patch_mod