!===============================================================================
!> Extension of short charadcteristis RT module with Multi3D specifics
!===============================================================================
MODULE m3d_rt_mod
  USE io_mod
  use omp_lib
  USE trace_mod
  USE rt_mod
  USE rt_integral_mod
  USE multi3d_mod
  USE radau_mod
  use mpi_mod
  USE boundaries_mod
  USE math_mod
  USE units_mod, only: cgs
  USE bits_mod
  USE gpatch_mod
  USE scs_mod, only: dis_package_t
  use scaling_mod, only: scaling
  USE task_mod
  USE list_mod
  USE link_mod
  USE halo_mod
  USE kinds_mod
  implicit none
  private
  type, public, extends(rt_t):: m3d_rt_t
    !---------------------------------------------------------------------------
    ! Before adding things here, check of the variable already exists in rt_t,
    ! and consider if it should be shared, or Multi3D specific (new name)
    !---------------------------------------------------------------------------
    real:: task_variable
    !---------------------------------------------------------------------------
    ! Add arrays here to store population numbers are whatever else
    !---------------------------------------------------------------------------
    !real, pointer, contiguous:: Jint(:,:,:,:) => null()
    !real, pointer, contiguous:: Npop(:,:,:,:) => null()
    integer:: fbin                                     ! Frequency bin
    integer,dimension(:,:),allocatable :: nu_regions   ! Frequency index limits
    integer :: nu_s,nu_e, n_updates
    type(dis_package_t)  :: dispackage    ! Dispatch <-> Multi communicator
    integer :: i_n,i_n_lte,i_ne           ! Index range of pops and ne in mem arr
    integer :: i_chi,i_snu,i_inu,n_nu     ! Index range of chi, snu and inu 
    real(8),dimension(:,:,:),pointer :: temp8 => null()  ! Input atmosphere (double prec.cast)
    real(8),dimension(:,:,:),pointer :: rho8  => null()
    real(8),dimension(:,:,:),pointer :: ne8   => null()
    real(8),dimension(:,:,:),pointer :: vx8   => null()
    real(8),dimension(:,:,:),pointer :: vy8   => null()
    real(8),dimension(:,:,:),pointer :: vz8   => null()

    class(rt_t), pointer :: opac  => null()  ! pointer to head of opac sheep list, similarly to omeg

    logical :: is_linked = .false.
    character(len=64) :: meta_file, ray_name
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: bndries
    procedure:: update
    procedure:: after_update
    procedure:: output
    procedure:: get_nu_regions
    procedure:: reduce_omegaX
    procedure:: link_dispackage
    procedure:: nbor_relations
    procedure:: init_level
    procedure:: init_task_list
    procedure:: task_info
    procedure:: load_initial_conditions
    procedure:: reduce_angle_quadrature
    procedure:: diagnose
    procedure:: convert_mhd_variables
    procedure:: save_meta4
    procedure:: save_meta8
    procedure:: setup_metafile
    procedure:: pre_init_sheep_list
  end type

  type, public, extends(halo_t):: halo_m3d_t
  contains
    procedure:: info
  end type


  integer, save:: verbose=0                           ! verbosity
  real,    save:: rt_freq=0                           ! how often RT should be updated WRT MHD -- 0=unlimited
  logical, save:: dnload_mhd=.true.
  logical, save:: save_output=.true.
  integer, save:: maxiter=1
  character(len=64), save:: out_folder = './out/'
CONTAINS

!===============================================================================
!> Pre-initialize: read parameters, if any.  Probably we can rely on the
!> existing parameters inside Multi3D and DISPATCH.
!===============================================================================
recursive SUBROUTINE pre_init (self)
  class(m3d_rt_t), target:: self
  !-----------------------------------------------------------------------------
  class(m3d_rt_t), pointer:: omega
  character(128), save:: folder='input_multi3d'
  integer:: n_omega, iostat
  integer, save:: n_mu=1, n_phi=1, n_bin=2, nt=2, n_warmup=3, ng(3)=[4,4,1]
  real,    save:: rt_grace=0.01
  real(8), save:: rt_timestep=0.01_8
  logical, save:: first_time=.true.
  logical, save:: on=.true.
  !-----------------------------------------------------------------------------
  namelist /multi3d_params/ folder, on, verbose, n_mu, n_phi, n_bin, save_output, &
    out_folder, maxiter, rt_grace, rt_timestep
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%pre_init')

  if (first_time) then
    first_time = .false.
    call io%lock%set ('multi3d_t')
    rewind (io%input)
    read (io%input, multi3d_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, multi3d_params)
    end if
    call io%namelist_warning ('multi3d_params', iostat)
 
    multi3d%dir = folder
    call multi3d%pre_init()                                                  ! Pre-init m3d
    call io%lock%unset ('multi3d_t')
  end if

  ! Multi3D pre init is ran only once per OMP thread
  ! This means that another dispackages of the same thread can not be filled 
  ! within the same routine. Fill them here with thread global variables
  self%dispackage%ref_atom     => multi3d%ref_atom
  self%dispackage%ref_spectrum => multi3d%ref_spectrum

  ! Set the scaling pointer in the dispackage
  self%dispackage%rho_scale => scaling%d
  self%dispackage%p_scale   => scaling%p
  self%dispackage%vel_scale => scaling%u

  self%dispackage%ng           = self%patch%ng
  self%dispackage%geometry%nx  = self%patch%n(1)
  self%dispackage%geometry%ny  = self%patch%n(2)
  self%dispackage%geometry%nz  = self%patch%n(3)

  if (self%pre_initialized_rt) call io%abort ('double rt pre_init')
  self%pre_initialized_rt = .true.
  if (self%pre_initialized) print *, 'pre_init issued somewhere'
  call io%lock%unset ('rt_t%pre_init')

  ! ----------------------------------------------------------------------------
  ! Mute a warning about bad nw, nv, and ng - these are set later.
  ! nw and ng - just below,
  ! nv - in the init procedure, where we will know, how many RT directions there
  ! will be. 
  ! ----------------------------------------------------------------------------
  !self%nw = self%patch%nw
  !self%nv = self%patch%nv
  self%ng = self%patch%ng
  call self%gpatch_t%pre_init()

  ! ----------------------------------------------------------------------------
  ! Fill namelist variables
  ! ----------------------------------------------------------------------------
  self%on = on
  self%nt = nt
  self%ng = self%patch%ng
  self%nw = 1                                                                  ! no time derivatives
  self%n  = self%patch%n
  self%n_warmup = n_warmup
  self%n_bin = n_bin
  self%grace = rt_grace
  self%rt_timestep = rt_timestep
  self%n_mu  = n_mu
  self%n_phi = n_phi
  self%in_task_list = .false.

  if (self%on) then 
    call self%task_t%init_unique (same=.true.) ! obtain an unique ID
    self%type    = 'm3d_rt_t'
    !======================================================================================
    ! Omega0: > Collects the transition rates from omegaX tasks
    !         > Stores population numbers + rates 
    !
    ! OmegaX: > Computes radiative transfer
    !         > Stores chi, snu and inu in its frequency range
    !
    ! The storage that needs to be communicated between the tasks is chi,snu and inu !!!
    ! The angle tasks do not need to know the population numbers in the ghost zones
    ! if the opacity is passed!
    !=====================================================================================
    if (self%i_omega==0) then ! i_omega is set to zero by default in rt_mod
      self%kind    = 'm3d_rt_sc0'
      self%main => self         ! convenience
      if (self%n_mu > 0) then
        !-------------------------------------------------------------------------
        ! Number of ray directions
        !-------------------------------------------------------------------------
        !n_omega = merge (1+self%n_phi*(self%n_mu-1), 0, self%n_mu > 0)
        n_omega = self%n_phi*self%n_mu ! only for comparison with 'set_a2'
        self%n_omega = 2*n_omega*n_bin
        call self%pre_init_sheep_list(self%omeg, 'omeg_sheep') ! this calls self%pre_init recursively
        !call self%pre_init_sheep_list(self%opac, 'opac_sheep') ! this calls self%pre_init recursively
      end if
      
      self%n_nu = self%dispackage%ref_spectrum%nnu
      self%dispackage%nu_s  = 1
      self%dispackage%nu_e  = self%n_nu

      self%nv      = 2   ! population numbers (+ heating)
      self%i_n     = 1
      self%i_n_lte = 2
      self%nw      = self%dispackage%ref_atom%nlevel
    elseif(self%kind == 'omeg_sheep') then
      self%nv = 3
      self%nw = self%n_nu  ! chi, snu and inu (for every nu in this patch)
      self%i_chi = 1
      self%i_snu = 2
      self%i_inu = 3
    end if
    self%mv = self%nv
    self%nt = 1
    self%n_updates = 0
    if (save_output) call self%setup_metafile()
  end if
  call trace%end ()
END SUBROUTINE pre_init

SUBROUTINE pre_init_sheep_list(self, first_sheep, kind)
  class(m3d_rt_t), target :: self
  class(rt_t), pointer :: first_sheep
  character(len=*), intent(in) :: kind
  class(m3d_rt_t), pointer :: sheep
  real, dimension(:), allocatable:: mu, w_mu, mux, muy, muz
  integer :: i_omega, dir, i_mu, i_phi, n, fbin

  call self%get_nu_regions(self%n_bin)            ! Get frequency domains for all tasks
  
  !-------------------------------------------------------------------------
  ! Assign angles and integration weights
  !-------------------------------------------------------------------------
  allocate (mu  (self%n_mu))
  allocate (w_mu(self%n_mu))
  !call radau%init (self%n_mu, mu, w_mu)
  mu = [0.57735026] ! only for comparison with 'set_a2'
  mux = [-0.57735026,  0.57735026,  0.57735026, -0.57735026]
  muy = [ 0.57735026,  0.57735026, -0.57735026, -0.57735026]
  w_mu = 1
  !-----------------------------------------------------------------------
  ! Assign angles and initialize clones
  !-----------------------------------------------------------------------
  i_omega = 0
  do dir=-1,1,2
    do i_mu=1,self%n_mu
      n = self%n_phi
      !if (i_mu == self%n_mu) n = 1
      do i_phi=1,n
        do fbin=1, self%n_bin
          allocate(sheep)
          sheep%main => self
          sheep%master_id = self%id
          if (associated(first_sheep)) then
            sheep%prev => self%tail
            sheep%prev%next => sheep
            self%tail => sheep
          else
            first_sheep => sheep
            self%tail => sheep
          end if

          i_omega = i_omega+1
          sheep%i_omega = i_omega
          sheep%kind = kind

          !sheep%mu = dir*mu(i_mu)
          sheep%mu = -dir*mu(i_mu)

          !sheep%rt_phi = math%pi*(modulo(real(2*(i_phi-1))/self%n_phi + &
          !  (dir+1)/2.0, 2.0) - 1.0)
          !if (i_mu == self%n_mu) sheep%rt_phi = 0.0
          sheep%rt_phi = mux(i_phi) ! only for comparison with 'set_a2'

          ! ------------------------------------------------------------------------
          ! Determine closest axis direction for omega/=0
          ! ------------------------------------------------------------------------
          !sheep%rt_hat(1) = sin(acos(sheep%mu))*cos(sheep%rt_phi)   ! sin(theta)*cos(phi), dx
          !sheep%rt_hat(2) = sin(acos(sheep%mu))*sin(sheep%rt_phi)   ! sin(theta)*sin(phi), dy
          !sheep%rt_hat(3) = sheep%mu                               ! cos(theta), dz

          sheep%rt_hat(1) = mux(i_phi)
          sheep%rt_hat(2) = muy(i_phi)
          sheep%rt_hat(3) = sheep%mu
          sheep%axis = maxloc(abs(sheep%rt_hat), 1) ! direction of max propagation

          sheep%dispackage%mux   = dble(sheep%rt_hat(1))
          sheep%dispackage%muy   = dble(sheep%rt_hat(2))
          sheep%dispackage%muz   = dble(sheep%rt_hat(3))

          sheep%w       = w_mu(i_mu)*math%pi4/n/2.
          sheep%patch   => self%patch

          sheep%fbin    = fbin
          sheep%nu_s    = self%nu_regions(fbin,1)  ! Starting frequency index in spectrum freq. grid
          sheep%nu_e    = self%nu_regions(fbin,2)  ! Ending frequency index in spectrum freq. grid

          sheep%dispackage%fbin = sheep%fbin
          sheep%dispackage%nu_s = sheep%nu_s
          sheep%dispackage%nu_e = sheep%nu_e
          sheep%dispackage%nfint = self%n_bin
          sheep%n_nu = sheep%nu_e-sheep%nu_s + 1

          call sheep%pre_init
        end do
      end do
    end do
  end do
END SUBROUTINE pre_init_sheep_list

SUBROUTINE setup_metafile(self)
  class(m3d_rt_t),intent(inout):: self
  integer :: meta_unit

  if (self%i_omega == 0) then
    write(self%ray_name,'(A1, I3.3, A9)') 'p', self%patch%id, '_shepherd'
  else
    write(self%ray_name,'(A1, I3.3, A6, I3.3)') 'p', self%patch%id, '_sheep', self%i_omega
  endif

  self%meta_file = trim(out_folder) // trim(self%ray_name) // '_meta.txt'
  open(newunit=meta_unit, file=self%meta_file, form='formatted')

  write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'position', self%patch%position

  if (self%i_omega == 0) then
    write(meta_unit, '(A8, 1X, I4)') 'nmu',   self%n_mu
    write(meta_unit, '(A8, 1X, I4)') 'nphi',  self%n_phi
    write(meta_unit, '(A8, 1X, I4)') 'nfbin', self%n_bin
  else
    write(meta_unit, '(A8, 1X, SP, 3(F12.8))') 'rt_hat', self%rt_hat
    write(meta_unit, '(A8, 1X, I4)') 'fbin', self%fbin
  endif

  close(meta_unit)
END SUBROUTINE


! PE: obtain frequency splitted regions
subroutine get_nu_regions(self,n_splits) 
    class(m3d_rt_t),intent(inout):: self
    integer,intent(in) :: n_splits
    integer :: i,nu_bin_size

    ! Simple approach, can be adjusted in the future if needed (based on active lines etc.)
    if (allocated(self%nu_regions)) deallocate(self%nu_regions)
    allocate(self%nu_regions(n_splits,2))

    nu_bin_size =  floor(self%dispackage%ref_spectrum%nnu / real(n_splits))
    self%nu_regions(1,1)        = 1                         ! Starpoint
    self%nu_regions(1,2)        = nu_bin_size               ! Starpoint
    do i = 2, n_splits
        self%nu_regions(i,1) = self%nu_regions(i-1,2) + 1
        self%nu_regions(i,2) = self%nu_regions(i,1) + nu_bin_size - 1
    enddo
    self%nu_regions(n_splits,2) = self%dispackage%ref_spectrum%nnu  ! Endpoint
end subroutine get_nu_regions

!===============================================================================
!> Initialize the rt_t data types.  To setup for n_warmup initial solutions, we
!> set the time rt_t times = -n_warmup * rt_t%timestep, making sure to call
!> rt_t%rt_eos for the initial MHD state. This subroutine is recursive, as
!> omega0 initializes all the omegaX. 
!===============================================================================
SUBROUTINE init (self)
  class (m3d_rt_t):: self
  ! ............................................................................
  class (gpatch_t), pointer :: p
  class (rt_t), pointer:: omega
  integer :: i, shape_patch(6)
  logical, save :: first_time = .true.
  ! ----------------------------------------------------------------------------
  if (self%initialized) then
    print *, self%id, self%i_omega
    call io%abort ('double rt init')
  end if 
  call trace%begin ('rt_t%init')
  self%initialized = .true.
  p => self%patch
  if (self%on) then
    do i=0, self%n_omega
      if (i==0) omega => self%main
      if (i==1) omega => self%omeg
      if (i>=2) omega => omega%next
      ! ----------------------------------------------------------------------------
      ! relevant synchronization items for cartesian setup
      ! ----------------------------------------------------------------------------
      omega%ipos       = p%ipos
      omega%rank       = p%rank
      omega%box        = p%box
      omega%size       = p%size
      omega%position   = p%position
      omega%ds         = p%ds
      omega%origin     = p%origin
      omega%llc_cart   = p%llc_cart
      omega%centre_nat = p%centre_nat
      omega%llc_nat    = p%llc_nat
      omega%status     = p%status
      omega%periodic   = p%periodic
      omega%ng         = p%ng
      omega%n          = p%n

      if (verbose > 0) &
        write (io_unit%mpi,'(3i4,3f8.3,i4,a)') p%id, self%id, self%i_omega, p%position, p%rank,'pID, sID, omega, pos, rank'

      if (verbose == -1) &
        write (io_unit%log,*) 'rt_t%init: id, mhd%id i_omega =', self%id, p%id, self%i_omega

      
      call omega%gpatch_t%init
      omega%mem = tiny(1.)
      omega%it = 1
      omega%new = omega%it

      if (verbose > 0) write (io_unit%log,*) 'call self%patch_t%setup'

      !-----------------------------------------------------------------------------
      ! Synchronize with MHD patch
      ! ----------------------------------------------------------------------------
      omega%level      = 1 ! p%level
      omega%restart    = p%restart
      omega%warmup_end_time = p%time
      omega%t(omega%it) = omega%time
      ! -- add an `eps` to avoid round-off errors
      omega%out_next = (int(omega%time/io%out_time+1.0e-4)+1)*io%out_time
      omega%dt_fixed = omega%rt_timestep
      omega%dtime    = omega%rt_timestep
      omega%dt       = omega%rt_timestep
      omega%time     = p%time - (omega%n_warmup-omega%grace)*omega%dtime
      omega%t        = omega%time
      omega%qr_time  = omega%time - 2.0*omega%grace*omega%dtime

      ! ----------------------------------------------------------------------------
      ! Set omega0 time to be marginally after omegaX, so radiative rates are
      ! computed first and populations are updated after
      ! ----------------------------------------------------------------------------
      if (i == 0) then
        omega%time  = omega%time + 2.0*omega%grace*omega%dtime
      end if
      omega%iout    = omega%restart+1
      omega%verbose = verbose

      ! ----------------------------------------------------------------------------
      ! Do not check density or validate with it. This is rather important to be set
      ! ----------------------------------------------------------------------------
      call omega%set (bits%no_density)
      call omega%set (bits%internal)
      
      ! ----------------------------------------------------------------------------
      ! Allocate and set unsigned and pervolume. Opacity, source function,
      ! intensity and populations are all positive definite
      ! Although quantities are generally per volume, but density is not stored
      ! locally, thus this flag must not be set.
      ! ----------------------------------------------------------------------------
      allocate (omega%unsigned(self%nv), omega%pervolume(self%nv))
      omega%pervolume = .false.
      omega%unsigned  = .true.

      ! ----------------------------------------------------------------------------
      ! Init boundary conditions
      ! ----------------------------------------------------------------------------  
      !call omega%init_bdries

      select type (o=>omega)
      class is (m3d_rt_t)
        shape_patch = shape(o%patch%mem)
        if (.not.allocated(o%dispackage%j)) allocate(o%dispackage%j(shape_patch(1),shape_patch(2),shape_patch(3),&
                                                          o%dispackage%nu_s:o%dispackage%nu_e))
        if (o%i_omega.eq.0) then
          ! Init the real8 cast for atmos
          if (.not.associated(o%temp8)) allocate(o%temp8(shape_patch(1),shape_patch(2),shape_patch(3)))
          if (.not.associated(o%rho8))  allocate(o%rho8(shape_patch(1),shape_patch(2),shape_patch(3)))
          if (.not.associated(o%ne8))   allocate(o%ne8(shape_patch(1),shape_patch(2),shape_patch(3)))
          if (.not.associated(o%vx8))   allocate(o%vx8(shape_patch(1),shape_patch(2),shape_patch(3)))
          if (.not.associated(o%vy8))   allocate(o%vy8(shape_patch(1),shape_patch(2),shape_patch(3)))
          if (.not.associated(o%vz8))   allocate(o%vz8(shape_patch(1),shape_patch(2),shape_patch(3)))

          o%dispackage%nstar => o%mem(:,:,:,o%i_n_lte,o%it,:)
          o%dispackage%n     => o%mem(:,:,:,o%i_n    ,o%it,:)
          o%dispackage%jold  => o%dispackage%j
    
          o%boundaries%topbot => o%patch%boundaries%topbot
          o%boundaries%temp   => o%patch%boundaries%temp
          o%boundaries%rho    => o%patch%boundaries%rho
          o%boundaries%pe     => o%patch%boundaries%pe

        else
          select type (om=>o%main)
          class is (m3d_rt_t)
            o%temp8 => om%temp8
            o%rho8  => om%rho8
            o%ne8   => om%ne8
            o%vx8   => om%vx8
            o%vy8   => om%vy8
            o%vz8   => om%vz8

            o%dispackage%n     => om%dispackage%n
            o%dispackage%nstar => om%dispackage%nstar
            o%dispackage%jold(1:,1:,1:,o%dispackage%nu_s:) => &
                          om%dispackage%j(:,:,:,o%dispackage%nu_s:o%dispackage%nu_e)    ! The omegax task refers to its master J as Jold    
                          
            o%boundaries%topbot => om%boundaries%topbot
            o%boundaries%temp   => om%boundaries%temp
            o%boundaries%rho    => om%boundaries%rho
            o%boundaries%pe     => om%boundaries%pe

            !allocate(o%dispackage%ups)
            !o%dispackage%ups = om%dispackage%ups
            !o%dispackage%ups => om%dispackage%ups
          end select
        endif

        o%dispackage%temp => o%temp8
        o%dispackage%rho  => o%rho8
        o%dispackage%ne   => o%ne8
        o%dispackage%vx   => o%vx8
        o%dispackage%vy   => o%vy8
        o%dispackage%vz   => o%vz8
        
        o%dispackage%geometry%nx = o%n(1)
        o%dispackage%geometry%ny = o%n(2)
        o%dispackage%geometry%nz = o%n(3)
        o%dispackage%ng          = o%ng

        o%dispackage%wm    = dble(o%w)
        o%dispackage%ds    = dble(o%ds * scaling%l)

        call multi3d%init(o%dispackage,skip_ray_in=(o%i_omega.eq.0))
        o%dispackage%task_id = o%task_t%id                                            ! Save the task id as the id of the dispackage for cross referencing
      end select

      if (.not. omega%on) then
        !---------------------------------------------------------------------------
        ! Deallocate the self%mem for this task, but keep the task, for simplicity
        !---------------------------------------------------------------------------
        call omega%dealloc()
        omega%n_omega = 0
      end if
    end do
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
SUBROUTINE convert_mhd_variables(self)
  class(m3d_rt_t),target :: self

  !print *, 'temp',self%patch%mem(5,5,:,self%patch%idx%e,  1, 1)
  !print *, 'shape',shape(self%patch%mem(:,:,:,self%patch%idx%e,  1, 1)),shape(self%temp8)
  !print *, 'assoc',associated(self%temp8)

  self%temp8 = dble(self%patch%mem(:,:,:,self%patch%idx%e,  1, 1))
  self%rho8  = dble(self%patch%mem(:,:,:,self%patch%idx%d,  1, 1)) * scaling%d
  where (self%temp8 > 0) 
    self%ne8 = dble(self%patch%mem(:,:,:,6,                 1, 1)) * scaling%p / (cgs%k_b * self%temp8)
    self%vx8 = dble(self%patch%mem(:,:,:,self%patch%idx%px, 1, 1)) / self%rho8 *scaling%d*scaling%u
    self%vy8 = dble(self%patch%mem(:,:,:,self%patch%idx%py, 1, 1)) / self%rho8 *scaling%d*scaling%u
    self%vz8 = dble(self%patch%mem(:,:,:,self%patch%idx%pz, 1, 1)) / self%rho8 *scaling%d*scaling%u
  end where
END subroutine convert_mhd_variables
!===============================================================================

!===============================================================================
SUBROUTINE load_initial_conditions(self)
    class(m3d_rt_t),target :: self
    class(rt_t),pointer    :: omega_ptr

    call self%convert_mhd_variables()

    self%mesh%lower_boundary = self%patch%mesh%lower_boundary
    self%mesh%upper_boundary = self%patch%mesh%upper_boundary
    self%boundaries%topbot   = self%patch%boundaries%topbot
    self%boundaries%temp     = self%patch%boundaries%temp
    self%boundaries%rho      = self%patch%boundaries%rho
    self%boundaries%pe       = self%patch%boundaries%pe

    self%dispackage%atom_rad%wphi => self%dispackage%atom_rad%wphi_contribution          

    !call self%init_bdries

    call multi3d%pre_update (self%dispackage, skip_ray_in = .true.)
    call self%dispackage%update_populations_lte ()                                     ! compute lte populations
    self%dispackage%n = self%dispackage%nstar                                          ! Set the initial population to LTE
    self%dispackage%j = 0.0d0                                                          ! Set the initial Jnu to 0. This is only relevant for scattering
                                                                                       ! The omega0 task does not need a jold. 
                                                                                       ! Jold from the sub-task point to the right spot in j
    omega_ptr => self%omeg                                                             ! Loop over all sub tasks and initialize
    do while (associated(omega_ptr))
        if (omega_ptr%i_omega .eq. 0) then
            omega_ptr => omega_ptr%next
            cycle
        endif

        select type (om=>omega_ptr)
        class is (m3d_rt_t)                                                             ! If this is a m3d_rt_t omegaX, initialize the source function
            call om%link_dispackage   ()                                                ! Link the dispackage source function to the current slice
            call multi3d%pre_update   (om%dispackage, skip_ray_in = .false.)            ! Sync dispackage and set up MULTI arrays 
            om%dispackage%atom_rad%wphi => self%dispackage%atom_rad%wphi_contribution   ! Set the pointer to the quadrature array of the master
            call om%dispackage%update_source_function ()                                ! Compute the source function for every height coordinate
            om%dispackage%inu_mem = 0 !!!om%dispackage%snu_mem                        ! Set the initial intensity equal to the source function

        end select
        
        omega_ptr => omega_ptr%next
    enddo
    call self%reduce_angle_quadrature ()                                                ! collect the angle quadrature
END SUBROUTINE

subroutine reduce_angle_quadrature(self)
  class(m3d_rt_t),intent(inout) :: self
  class(rt_t),pointer :: omega

  ! Collect the wphi arrays from all omegaX tasks
  self%dispackage%atom_rad%wphi = 0.0d0
  omega => self%main%omeg
  do while (associated(omega))
    if ((omega%master_id.ne.self%id).or.(omega%i_omega.eq.0)) then
      omega => omega%next
      cycle                                         
    endif
    
    select type (o=>omega)                            
    class is (m3d_rt_t)
      self%dispackage%atom_rad%wphi = self%dispackage%atom_rad%wphi + o%dispackage%wm * o%dispackage%atom_rad%wphi_contribution
    end select

    omega => omega%next
  enddo
  self%dispackage%atom_rad%wphi = 4.0d0 * cgs%pi / self%dispackage%atom_rad%wphi  ! wphi contains the sum of wphi from each angle
end subroutine reduce_angle_quadrature

!===============================================================================
!> Update the population numbers and other dynamic task data.  Consider spliting
!> this into three steps; put data to Multi3D in pre_update, update, and get
!> data back from Mulit3D in post_update.
!===============================================================================
SUBROUTINE update (self)
  use timer_mod
  class(m3d_rt_t) :: self
  logical, save :: first_time = .true.
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  integer :: lvl, nz
  character(len=3) :: solver
  character(len=48) :: text
  !-----------------------------------------------------------------------------
  call trace%begin ('m3d_rt_t%update')

  self%it = 1
  lvl = 3
  nz = size(self%dispackage%n, 3)
  self%n_updates = self%n_updates+1
  self%dispackage%topbot = self%patch%boundaries%topbot

  if (self%n_updates > maxiter) then
    print *, 'Iteration: ', self%n_updates
    stop
  endif

  call timer%tic ()

  if (self%i_omega /= 0) then                                                          ! rad. transfer solver
    solver = 'rad'
    self%dispackage%position = self%patch%position
    self%dnload_time         = self%time - 2.0 * self%grace * self%dtime               ! Force download at the beginning of this omega-nu task
    call self%dnload          ()                                                       ! it will assume that the upwind tasks are already done and 
                                                                                       ! that the mem object is downloaded! (So Chi,Inu and Snu are available in the ghost zones)
    call self%link_dispackage ()                                                       ! Link the dispackage source function to the current slice
    call multi3d%pre_update   (self%dispackage, skip_ray_in = .false.)                 ! Sync dispackage and set up MULTI arrays 
    call self%dispackage%update_radiation()                
    if (save_output) call self%output()
      
  else
    solver = 'pop'
    call self%dispackage%sync         ()
    call self%reduce_angle_quadrature ()                                                    ! collect the angle quadrature
    call self%reduce_omegaX           ()                                                    ! Reduce the rates etc. obtained from angle tasks

    !print *, 'old population sum: ', sum(self%dispackage%n(5,5,:,:))
    !print *, 'old population numbers (level ',lvl,'): ', self%dispackage%n(5,5,2:nz-1,lvl)
    call self%dispackage%update_populations()
    !call self%after_update ()                                                              ! Shift the populations
    !print *, 'new population sum: ', sum(self%dispackage%n(5,5,:,:))
    print *, 'new population numbers (level ',lvl,'): ', self%dispackage%n(5,5,2:nz-1,lvl)
    print *, ""
  endif

  write(text, *), 'time to solve ' // solver // ' task ', self%task_t%id
  call timer%toc (text, 1)

  !print *, '                  ...task', self%task_t%id, 'updated.'
  !print *, '====================================================='

  !call multi3d%post_update ! (... argument list ...)
  !self%dnload_time = self%time - 2.0 * self%grace * self%dtime

  call trace%end ()
END SUBROUTINE update

subroutine diagnose(self)
  class(m3d_rt_t),intent(in),target  :: self
  integer,dimension(:),      pointer :: is,ie

  is => self%dispackage%irp(:,1)
  ie => self%dispackage%irp(:,2)

  !$omp critical

  print *, ''
  print *, 'Task with id', self%task_t%id, 'has pre-updated.'
  print *, 'The dispackage in question is assigned to the task',self%dispackage%task_id
  print *, 'dim. parameters are (is,ie)', self%dispackage%irp
  print *, 'atmos parameters are (t,ne,r)', self%dispackage%temp(is(1),is(2),is(3)),&
                                            self%dispackage%ne(is(1),is(2),is(3)),&
                                            self%dispackage%rho(is(1),is(2),is(3))
  print *, ''

  !$omp end critical
end subroutine diagnose

subroutine after_update(self) ! Can this safely be called post_update?
  class(m3d_rt_t),intent(inout) :: self
  integer :: it_l,in_l,in_lte_l

  if (self%i_omega.eq.0) then
    ! Shift the arrays one time slice ahead
    it_l     = self%it
    in_l     = self%i_n
    in_lte_l = self%i_n_lte
  
    self%mem(:,:,:,in_l,    it_l + 1,:) = self%dispackage%n
    self%mem(:,:,:,in_lte_l,it_l + 1,:) = self%dispackage%nstar
    self%dispackage%n     => self%mem(:,:,:,in_l,    it_l + 1,:)
    self%dispackage%nstar => self%mem(:,:,:,in_lte_l,it_l + 1,:)
  endif

end subroutine after_update

subroutine reduce_omegaX(self)
  class(m3d_rt_t),target :: self
  class(rt_t),pointer    :: omeg
  integer,dimension(:,:),pointer :: ir
  real(8), dimension(:,:,:,:), pointer :: wphi, Jnu
  real(4), dimension(:,:,:,:), pointer :: Inu
  real(8) :: wmu_4pi
  integer :: nu_s, nu_e

  ! Index range of the patch (no ghost zones)
  ir => self%dispackage%irp

  self%dispackage%j = 0.0d0                           ! Zero the old radiation field
  self%dispackage%atom_rad%gamma_ij = 0.0d0
  self%dispackage%atom_rad%gamma_ji = 0.0d0

  omeg => self%main%omeg
  do while (associated(omeg))                         ! Loop through all the angle tasks and save the rates etc.
    if ((omeg%master_id.ne.self%id).or.(omeg%i_omega.eq.0)) then
      omeg => omeg%next
      cycle                                           ! This is not a task for this master
    endif
    
    select type (o=>omeg)                             ! Add the neighbor quants to the quants of this task
    class is (m3d_rt_t)

      self%dispackage%atom_rad%gamma_ij   = self%dispackage%atom_rad%gamma_ij   + o%dispackage%atom_rad%gamma_ij
      self%dispackage%atom_rad%gamma_ji   = self%dispackage%atom_rad%gamma_ji   + o%dispackage%atom_rad%gamma_ji
      
      wmu_4pi= o%dispackage%wm / (4.0d0*cgs%pi)
      nu_s = o%dispackage%nu_s
      nu_e = o%dispackage%nu_e
      Jnu => self%dispackage%j   (ir(1,1):ir(1,2),ir(2,1):ir(2,2),ir(3,1):ir(3,2),nu_s:nu_e)
      Inu => o%dispackage%inu_mem(ir(1,1):ir(1,2),ir(2,1):ir(2,2),ir(3,1):ir(3,2),nu_s:nu_e)

      Jnu = Jnu + wmu_4pi * dble(Inu)
    end select  
    omeg => omeg%next
  enddo

  wphi => self%dispackage%atom_rad%wphi(ir(1,1):ir(1,2),ir(2,1):ir(2,2),ir(3,1):ir(3,2),:)

  self%dispackage%atom_rad%gamma_ij(:,:,:,1:size(wphi,4)) =  self%dispackage%atom_rad%gamma_ij(:,:,:,1:size(wphi,4)) * wphi
  self%dispackage%atom_rad%gamma_ji(:,:,:,1:size(wphi,4)) =  self%dispackage%atom_rad%gamma_ji(:,:,:,1:size(wphi,4)) * wphi
end subroutine reduce_omegaX

subroutine link_dispackage(self)
  class(m3d_rt_t), intent(inout) :: self
  class(rt_t),pointer            :: omeg

  self%dispackage%chi_mem(1:,1:,1:,self%nu_s:) => self%mem(:,:,:,self%i_chi,self%it,1:self%n_nu)
  self%dispackage%snu_mem(1:,1:,1:,self%nu_s:) => self%mem(:,:,:,self%i_snu,self%it,1:self%n_nu)
  self%dispackage%inu_mem(1:,1:,1:,self%nu_s:) => self%mem(:,:,:,self%i_inu,self%it,1:self%n_nu)

  call self%dispackage%sync ()
  
  !if ((self%i_omega .eq. 1 ).and. (self%id == 13)) then
  !  omeg => self%main%omeg
  !  do while (associated(omeg))                         
  !   
  !    if ((omeg%i_omega.ne.3)) then
  !      omeg => omeg%next
  !      cycle    
  !    
  !    else
  !      select type (o=>omeg)
  !      class is (m3d_rt_t)
  !        self%dispackage%ups = o%dispackage%ups
  !      end select 
  !      exit   
  !    
  !    endif
  !  enddo
  !endif


end subroutine link_dispackage

!===============================================================================
!> Do Multi3D specific output
!===============================================================================
SUBROUTINE output (self)
  use spectrum_var
  use parallel_mod
  use mconst
  class(m3d_rt_t):: self
  character(len=64) :: filename
  integer :: lunit
  !-----------------------------------------------------------------------------
  !call trace%begin ('m3d_rt_t%output')

  ! if at top-boundary
  if (self%patch%boundaries%topbot(1) .and. self%mu > 0) then
    write(*,*) 'Writing output for angle ', self%ray_name

    filename = trim(out_folder) // trim(self%ray_name) // '_lam.bin'
    open(newunit=lunit, file=filename, access='stream')
    write(lunit) 1e8 * cc / sp%nu(fs:fe)
    call self%save_meta8(filename, sp%nu)
    close(lunit)

    filename = trim(out_folder) // trim(self%ray_name) // '_inu.bin'
    open(newunit=lunit, file=filename, access='stream')
    write(lunit) self%mem(:,:,self%ui(3),self%i_inu,self%it,:)
    call self%save_meta4(filename, self%mem(:,:,self%ui(3),self%i_inu,self%it,:))
    close(lunit)
  endif

  !call multi3d%output
  !call trace%end ()
END SUBROUTINE output

!===============================================================================
!> The nbor relations among the omega
!
! The omega0 task has a different role than in the default setup!
! It will be called AFTER the omegaX tasks have solved the radiative transfer,
! because it will collect the rates and sets + solves the stat. eq. equations.
!===============================================================================
subroutine nbor_relations(self, another, needed, needs_me, download)
  class(m3d_rt_t)       :: self            ! current M3d rt task
  class(task_t),pointer :: another         ! another task in the list
  logical :: needed, needs_me, download    ! Is the other task: needed, 
                                           !    does it need me, should I download?
  class(m3d_rt_t),pointer :: rt_another
  class(task_t),pointer   :: other
  

  call trace%begin ('rt_t%nbor_relations')

  ! If this RT task has no RT on, it should have no nbor relations
  if (.not.(self%on)) then
    needed = .false.; needs_me = .false.; download = .false.
    call trace%end()
    return
  end if

  ! First check the main relations, based on location and level of refinement
  call self%task_t%nbor_relations(another, needed, needs_me, download)

  ! Now, if another is RT, add it appropriately, if it is MHD, it has to overlap
  ! and have different kind of relations
  select type (another)
  class is (m3d_rt_t)

    rt_another => another
    if (.not.(rt_another%on)) then
      ! if nbor RT is not on, we don't need any relations to it.
      needed = .false.; needs_me = .false.; download = .false.
      call trace%end()
      return
    end if

    if ((rt_another%level<self%level) .and. (rt_another%i_omega==self%i_omega)) then
      if (all(rt_another%size>(abs(self%position-rt_another%position)))) then
        print *, 'found RT parent', self%i_omega, self%level, rt_another%level
        self%parent => rt_another
        self%link%parent => rt_another%link
        self%parentid = rt_another%id
        return
      end if
    end if
    
    if ((rt_another%i_omega == 0).and.(self%i_omega==0)) then                   ! Both are omega0
      ! omega0 need all the relations with another omega0, thus defaut setting
      !if (dnload_mhd) download = .false.
      call trace%end()
      return
    else if ((rt_another%i_omega == 0).and.self%i_omega/=0) then                ! Self: X, other: 0
      if (self%main%id == rt_another%id) then
        ! non-omega0 needs the master omega0 for its level populations
        needed   = .true.
        needs_me = .true.
        download = .false. 
      else 
        ! non-omega0 doesn't need to know about non-overlapping nbor omega0 
        needed = .false.; needs_me = .false.; download = .false.
      end if
    else if (self%i_omega==0) then                                              ! Self: 0, other: X
      if (self%id==rt_another%main%id) then
        ! omega0 needs other omegax tasks
        needed   = .true.
        needs_me = .true.
        download = .false. 
      else 
        ! non-omega0 doesn't need to know about non-overlapping nbor omega0 
        needed = .false.; needs_me = .false.; download = .false.
      end if
  
    else if (self%needs(rt_another)) then              
      ! The omega task needs the upstream task,
      ! but the upstream task does not need the downstream omega task.
      if ((self%rank==mpi%rank).and.(rt_another%rank/=mpi%rank)) needed=.true.
      needs_me=.false.
    else if (rt_another%needs(self)) then
      ! Downstream RT task will need me, but it is not needed and I should not
      ! download data from it.
      if ((self%rank==mpi%rank).and.(rt_another%rank/=mpi%rank)) needs_me=.true.
      needed = .false.
      download = .false.
    else
      needed = .false.; needs_me = .false.; download = .false.
    end if
    
  class is (gpatch_t)
    other => another
    if (all(self%position==other%position)) then
      ! Check if this is the "master" MHD patch. It might not be associated yet, i.e.
      ! after moving to another MPI rank or changing level of refinement. in that 
      ! case, associate it.
      needed   = .true.
      needs_me = .false.
      download = .false.
      if (self%i_omega .eq. 0) needs_me = .true.   ! For returning the heating rate
      if (.not.(associated(self%patch))) self%patch => gpatch%cast2gpatch(other)
    else
      needed = .false.; needs_me = .false.; download = .false.
    end if 
  class default
    needed = .false.; needs_me = .false.; download = .false.
  end select

  call trace%end()
end subroutine nbor_relations

!===============================================================================
!> Setup initial states in the patches
!===============================================================================
SUBROUTINE bndries (self)
  class(m3d_rt_t):: self
  real(var_k), dimension(:,:,:,:,:), pointer :: mem
  !----------------------------------------------------------------------------
  call trace_begin('m3d_rt_t%bndries')
  mem => self%mem(:,:,:,:,self%it,:)
  if (size(self%boundaries%topbot) /= 2) then
    print *, 'boundary conditions not initialised for task ', self%id 
  elseif (any(self%boundaries%topbot)) then
    call self%dispackage%sync () ! needed to set global "fs" correctly
    call self%boundaries%boundary(mem)
  endif
  call trace%end()
END SUBROUTINE bndries


!===============================================================================
SUBROUTINE init_level (self)
  class(m3d_rt_t):: self
  !-----------------------------------------------------------------------------
  if (any(self%n > 1)) then
    ! The level does not depend on the patch size in non-regular grids
    self%level = 1
  end if
END SUBROUTINE init_level

!===============================================================================
!> Add the (1+n_omega) RT tasks to the task list
!===============================================================================
SUBROUTINE init_task_list(self, task_list)
  class(m3d_rt_t) :: self
  class(task_t), pointer:: rt_task
  class(list_t), pointer:: task_list
  class(rt_t),   pointer:: omega
  class(link_t), pointer:: link
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%add_tasks_to')
  !-----------------------------------------------------------------------------
  if (self%on) then
    call append_link(self, task_list)
    omega => self%omeg
    do while(associated(omega))
      call append_link(omega, task_list)
      omega => omega%next
    end do
  end if
  call trace%end()
END SUBROUTINE init_task_list

!===============================================================================
!> Append RT tasks to the task list
!===============================================================================
SUBROUTINE append_link(patch, task_list)
  class(rt_t),   target :: patch
  class(link_t), pointer:: link
  class(list_t), pointer:: task_list

  patch%task_list => task_list
  if (.not.associated(patch%link)) then
    allocate (halo_m3d_t:: link)
    link%task => patch
    patch%link => link
  endif
  if (.not.patch%in_task_list) then
    call task_list%append (patch%link)
  endif
  patch%in_task_list = .true.
END SUBROUTINE append_link

!===============================================================================
!> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
!===============================================================================
SUBROUTINE task_info (self)
  class(m3d_rt_t):: self
  !-----------------------------------------------------------------------------
  write(io_unit%mpi,'(1x,a,12x,i6,2i4,7x,2l1,2x,3f9.3,i4,2f7.3,i4)') &
    'rt_t%task_info: id, rank, level, BV, pos, i_omega, mu, phi, i_omega =', &
    self%id, self%rank, self%level, &
    self%is_set(bits%boundary), self%is_set(bits%virtual), self%position, &
    self%i_omega, self%mu, self%rt_phi, self%i_omega
END SUBROUTINE task_info


!===============================================================================
!> Initialize the link, and possibly the nbor nodes
!===============================================================================
SUBROUTINE info (self, label)
  class(halo_m3d_t):: self
  character(*), optional:: label
  class(nbor_t), pointer:: nbor
  real(8):: distance(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('halo_t%info')
  if (verbose < 2) then
    call self%task%task_info()
    nbor => self%nbor
    do while (associated(nbor))
      distance = nbor%task%distance(self%task)
      select type(nt => nbor%task)
      class is(m3d_rt_t)        
        write(io_unit%mpi,'(3x,a16,a,i6,":",i2.2,i4,2x,3l1,2x,2l1,3x,3f9.3, a10, i4)') &
          nt%kind, ' nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
          nt%id, nt%level, nt%rank, nbor%needed, nbor%needs_me, &
          nbor%download, nt%is_set(bits%boundary), &
          nt%is_set(bits%virtual), distance/ &
          (0.5d0*(self%task%size+nt%size)), ' i_omega: ', nt%i_omega
      class default
        write(io_unit%mpi,'(3x,a16,a,i6,":",i2.2,i4,2x,3l1,2x,2l1,3x,3f9.3)') &
          nt%kind, ' nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
          nt%id, nt%level, nt%rank, nbor%needed, nbor%needs_me, &
          nbor%download, nt%is_set(bits%boundary), &
          nt%is_set(bits%virtual), distance/ &
          (0.5d0*(self%task%size+nt%size))
      end select
      nbor => nbor%next
    end do
  else
    call self%halo_t%info(label)
  end if
  call trace%end ()
END SUBROUTINE info

SUBROUTINE save_meta8(self, filename, arr)
  class(m3d_rt_t):: self
  character(len=*), intent(in) :: filename
  real(8), intent(in) :: arr(..)
  integer :: meta_unit

  open(newunit=meta_unit, file=self%meta_file, form='formatted', access='append')
  write(meta_unit, '(A64, I4, 8(I5))') filename, storage_size(arr), shape(arr)
  close(meta_unit)
END SUBROUTINE save_meta8

SUBROUTINE save_meta4(self, filename, arr)
  class(m3d_rt_t):: self
  character(len=*), intent(in) :: filename
  real(4), intent(in) :: arr(..)
  integer :: meta_unit

  open(newunit=meta_unit, file=self%meta_file, form='formatted', access='append')
  write(meta_unit, '(A64, I4, 8(I5))') filename, storage_size(arr), shape(arr)
  close(meta_unit)
END SUBROUTINE save_meta4

END MODULE m3d_rt_mod
