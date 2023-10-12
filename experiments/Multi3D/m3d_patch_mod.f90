MODULE m3d_patch_mod
    !USE rt_integral_mod
    !USE radau_mod
    !USE boundaries_mod
    !USE eos_mod
    USE os_mod
    USE list_mod
    USE halo_mod
    USE gpatch_mod
    USE link_mod
    USE task_mod
    USE guards_mod
    USE snapshots_mod
    !USE index_mod
    USE mesh_mod
    USE io_mod
    USE data_io_mod
    !USE aux_mod
    USE timer_mod
    !USE scaling_mod
    USE trace_mod
    USE mpi_mod
    USE omp_mod
    USE omp_lock_mod
    USE shared_mod
    USE omp_timer_mod
    USE math_mod
    USE bits_mod
    USE mconst
    USE handyfuncs_mod
    USE prec_mod, only: inf4
    !USE units_mod
    !USE scratch_mod
    USE m3d_halo_mod, only: m3d_halo_t, print_finished, print_download, lock_mem, unlock_mem, print_lock
    USE spectrum_mod, only: nu_tolerance
    USE atom_mod, only: atom_type, atomic_spectrum_type, lse_no_pivot, lse_solver, m1d_legacy_mode, nq_fac
    USE composition_mod, only: composition_type
    USE ray_mod
    USE m3d_debug_mod
    USE bbtrans_mod, only: disable_stark, disable_vdw, disable_ga
    implicit none
    private

    type, public :: key_link_t
      integer :: key(4)
      class(key_link_t), pointer :: next => null()
    end type key_link_t
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(gpatch_t) :: m3d_patch_t
      class(gpatch_t),    pointer :: patch => null()         ! pointer to mhd patch
      class(m3d_patch_t), pointer :: m3d_parent => null()    ! pointer to parent patch
      type(composition_type), pointer :: composition
      class(atomic_spectrum_type), pointer :: sp => null()
      type(atom_type), pointer :: atom => null()
      class(key_link_t), pointer :: key_list => null()
      class(m3d_patch_t), pointer :: last_allocator, last_deallocator
      real(4), dimension(:,:,:,:,:,:), pointer :: authority_mem => null()
      real(4), dimension(:,:,:,:,:), pointer :: out_mem => null() 
      integer, dimension(:,:), pointer :: periodic_idx
      integer, dimension(:), pointer :: periodic_i, periodic_j
      logical :: on=.false.       ! allows local turn-off
      logical :: in_task_list = .false.
      logical :: initialized  = .false., nullified = .false.
      real(8) :: t_offset=0, t_push=0
      real(8) :: allocation_time, deallocation_time
      integer :: iter=0, maxiter, master_maxiter, dir=0
      integer :: n_mu, n_nu, n_omega, n_out, zs, ze
      integer :: fcheck, pcheck(3), linecheck, lvlcheck
      logical :: out_ray=.false., static_mem=.false.
      logical :: compressed=.false., save_output=.false., saved=.false.
      logical :: print_update_info, checked_bounds=.false.
      logical :: free_self=.false.
      integer :: i_nu=0, i_mu=0
      integer :: i_ne=1, i_rho=2, i_temp=3, i_nh1=4, i_nh6=9
      integer :: i_chi=1, i_eta=2, i_snu=2
      integer :: i_Gij=1, i_Gji=2
      integer :: i_Inu=1
      character(3) :: mem_label='mem'
      logical :: prepend=.true.
    contains
      procedure :: pre_init
      procedure :: post_pre_init
      procedure :: init
      procedure :: pre_update
      procedure :: update
      procedure :: execute
      procedure :: post_update
      procedure :: finish
      procedure :: is_finished
      procedure :: nbor_relations => m3d_nbor_relations
      procedure :: pre_init_and_link
      procedure :: add_to_key_list
      procedure :: set_nv_and_nw
      procedure :: set_initial_time
      procedure :: set_maxiter
      procedure :: init_bdries
      procedure :: m3d_name
      procedure :: init_level
      procedure :: task_info
      procedure :: append_to_task_list
      procedure :: extrapolate_guards
      procedure :: is_ahead_of
      procedure :: periodic_grid
      procedure :: init_periodic_idx
      procedure :: rotate
      procedure :: distance_to
      procedure :: free_unused_memory
      procedure :: nbors_finished
      procedure :: skip_iterations
      procedure :: claim_nbor_place
      procedure :: push
      procedure :: allocate_mem
      procedure :: deallocate_mem
      procedure :: compress_mem
      procedure :: decompress_mem
      procedure :: check_mem
      procedure :: check_mem_NaN
      procedure :: check_access
      procedure :: m3d_patch_info
      procedure :: m3d_patch_info2
      procedure :: print_to_stdout
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: init_nbors_hash
      procedure :: check_nbor
      procedure :: hash_key
      procedure :: save
      procedure :: load_mem
      procedure :: set_place
      procedure :: queue_before
      procedure :: sort_nbors
    end type m3d_patch_t
  
    !-----------------------------------------------------------------------------
    integer, save, public:: verbose=0
    integer, public :: maxiter, ilambd, nz_total
    integer, public :: cswitch=1
    integer, public, save :: save_id=-1
    real(8), public, save, dimension(:), allocatable:: wmu, long_wmu
    real(8), public, save, dimension(:,:), allocatable:: vec, long_vec
    real(8), public:: convlim=1d-3
    logical, public, save :: converged=.false., use_vel=.true., auto_preconditioning=.false.
    logical, public, save :: preconditioning=.true., binned_opac=.false., save_Qrad=.false.
    logical, public, save :: check_NaN=.false., interp_log=.false.
    logical, public, save :: save_back_cont=.false., save_atom_line=.false.
    logical, public, save :: linear_long=.true., decouple_continuum=.false.
    logical, public, save :: detailed_prints=.false., comoving_frame=.false.
    logical, public, save :: use_density=.true., use_ne=.true., line_blanketing=.false.
    logical, public, save :: overlapping_bb=.true., overlapping_bf=.true.
    logical, public, save :: absmet_in_output=.false.
    logical, public, save :: print_start=.false., print_push=.false., print_send=.false.
    real(8), public, save :: dlam=100.
    character(len=128), save :: save_patch_kind=''
    character(len=64), public, save:: custom_mu=''
    real(8), parameter :: t_push=0.1d0
    real(4), public :: vel_max=0.0
    real(4), public :: ray_offset(2)=[0.,0.]
    integer, public :: pivot_z=0
    type(mpi_counter_t), public :: atom_update_counter, global_maxerr, global_pos
      
    public :: hash_string, pre_update, update, print_execution_history

  contains
  
  !===============================================================================
  SUBROUTINE pre_init (self)
    class(m3d_patch_t), target:: self
    !-----------------------------------------------------------------------------
    logical, save :: first_time=.true.
    logical, save :: on=.true.
    integer :: iostat, m(6), i
    character(len=256) :: command
    integer, save :: n_nu=1, itimer=0, ntheta=3, nphi=4
    character(len=32), save :: quad_scheme='radau2'
    character(len=32), save :: long_scheme='lobatto'
    real(8), dimension(:), allocatable:: cust_mu
    namelist /m3d_params/ on, verbose, fcheck, linecheck, pcheck, quad_scheme, &
    long_scheme, n_nu, maxiter, ilambd, lvlcheck, comoving_frame, save_id, &
    ntheta, nphi, check_NaN, custom_mu, lse_no_pivot, save_patch_kind, cswitch, &
    save_back_cont, detailed_prints, print_finished, print_start, print_push, &
    print_send, linear_long, interp_log, ray_offset, print_download, save_atom_line, &
    lse_solver, decouple_continuum, m1d_legacy_mode, nu_tolerance, disable_stark, &
    disable_vdw, disable_ga, pivot_z, use_vel, preconditioning, save_Qrad, &
    absmet_in_output, print_lock, auto_preconditioning
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_init', itimer=itimer)
    self%type = trim(self%kind)
    !-----------------------------------------------------------------------------
    ! Read parameters
    !-----------------------------------------------------------------------------
    call io%lock%set ('m3d_patch_t%pre_init')
    if (first_time) then
      first_time = .false.
      rewind (io%input)
      read (io%input, m3d_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, m3d_params)
      end if
      call io%namelist_warning ('m3d_params', iostat)
      write (io%namelist, m3d_params)
      if (mpi%master) then
        call os%mkdir(trim(io%run_dir)//'save')
        write (command, *) 'rm ', trim(io%run_dir), '*.txt ', trim(io%run_dir), '*.bin'
        call execute_command_line (trim(command), exitstat=iostat)
        write (command, *) 'rm ', trim(io%run_dir), 'save/*.txt ', trim(io%run_dir), 'save/*.bin'
        call execute_command_line (trim(command), exitstat=iostat)
      endif
      one_dimensional = all(self%patch%n(1:2) == 1)
      !-------------------------------------------------------------------------
      ! Assign angles and integration weights
      !-------------------------------------------------------------------------
      select case (quad_scheme)
      case ('disk_center', 'set_a2', 'set_a4', 'set_a6', 'set_a8', 'set_b4', 'set_b6', 'set_b8')
        call set_carlson_quadratures (quad_scheme, vec, wmu, up_and_down=.true.)
      case ('radau1', 'radau2', 'radau3', 'radau4', 'radau5')
        call set_radau_quadratures (quad_scheme, vec, wmu, up_and_down=.true.)
      case default
        print *, 'Quadrature keyword "'// trim(quad_scheme) // '" not recognized for SC solver. Stopping...'
        stop
      end select

      select case (long_scheme)
      case ('lobatto')
        call set_lobatto_quadratures (ntheta, nphi, long_vec, long_wmu)
      case ('disk_center', 'set_a2', 'set_a4', 'set_a6', 'set_a8', 'set_b4', 'set_b6', 'set_b8')
        call set_carlson_quadratures (long_scheme, long_vec, long_wmu, up_and_down=.false.)
      case ('radau1', 'radau2', 'radau3', 'radau4', 'radau5')
        call set_radau_quadratures (long_scheme, long_vec, long_wmu, up_and_down=.false.)
      case ('custom')
        call count_words(custom_mu, ntheta)
        allocate(cust_mu(ntheta))
        read(custom_mu, *) cust_mu 
        call set_custom_quadratures (cust_mu, nphi, long_vec, long_wmu)
      case ('grid+', 'gridx', 'grid*')
        call set_grid_quadratures(long_scheme, self%patch%ds, long_vec, long_wmu, up_and_down=.false.)
      case default
        print *, 'Quadrature keyword "'// trim(long_scheme) // '" not recognized for LC solver. Stopping...'
        stop
      end select

      if (m1d_legacy_mode) then
        decouple_continuum = .true.
        overlapping_bb = .false.
        overlapping_bf = .false.
        use_density = .false.
        use_vel = .false.
        comoving_frame = .true.
        preconditioning = .false.
        nq_fac = 2
      endif

      if (pivot_z == 0) pivot_z = nint(self%patch%box(3) / self%patch%ds(3))
    end if
  
    self%verbose = verbose

    self%n_mu = size(vec, 2)
    self%n_out = size(long_vec, 2)
    self%n_nu = n_nu
    self%out_ray = self%i_nu < 0 .or. self%i_mu < 0
  
    if (self%pre_initialized) print *, 'pre_init issued somewhere'
    call io%lock%unset ('m3d_patch_t%pre_init')
  
    call self%patch_t%pre_init()

    self%id = 0 ! this is assigned later
    self%on = on
    self%in_task_list = .false.
  
    ! ----------------------------------------------------------------------------
    ! relevant synchronization items for cartesian setup
    ! ----------------------------------------------------------------------------
    self%n          = self%patch%n
    self%ng         = 0
    self%ipos       = self%patch%ipos
    self%rank       = self%patch%rank
    self%box        = self%patch%box
    self%size       = self%patch%size
    self%position   = self%patch%position
    self%ds         = self%patch%ds
    self%origin     = self%patch%origin
    self%llc_cart   = self%patch%llc_cart
    self%centre_nat = self%patch%centre_nat
    self%llc_nat    = self%patch%llc_nat
    self%status     = self%patch%status ! be aware that this copies all bits
    self%periodic   = self%patch%periodic
    self%mesh       => self%patch%mesh
    self%task_list  => self%patch%task_list

    self%ze = self%zs + self%n(3)

    select case (self%kind)
    case ('atom_patch', 'out_rad')
      self%save_output = .true.
    case default
      self%save_output = string_in_text(save_patch_kind, self%kind)
    end select

    self%save_output = self%save_output .or. self%id == save_id

    ! ----------------------------------------------------------------------------
    ! Do not check density or validate with it. This is rather important to be set
    ! ----------------------------------------------------------------------------
    call self%set (bits%no_density)
    call self%set (bits%no_filled)
    
    call self%set_place(self, -2)
    call self%init_level()

    self%nt = 1
    self%time  = self%patch%time
    self%dtime = self%patch%dtime
    self%dt_fixed = self%patch%dtime
    self%out_next = self%patch%out_next
    call self%task_t%init()
    self%lock%kind = 'm3d'
    self%lock%id = self%id
    self%key = self%hash_key() ! overwrite usual hash key

    if (.not. associated(self%mem_lock)) then
      allocate (self%mem_lock(self%nt))
      do i=1,self%nt
        call self%mem_lock(i)%init ('meml', self%id)
      end do
    endif

    call self%set_maxiter()
    call trace%end (itimer)
  END SUBROUTINE pre_init

  !===============================================================================
  SUBROUTINE post_pre_init (self)
    class(m3d_patch_t) :: self
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin('m3d_patch_t%post_pre_init', itimer=itimer)
    !-----------------------------------------------------------------------------
    ! Preliminary region information, relative to the task position.  This may
    ! change in patch_t%init, but the settings below allow overlap() calls to work
    !-----------------------------------------------------------------------------
    self%authority%l = - self%size*0.5_8
    self%authority%u = + self%size*0.5_8
    self%interest%l  = self%authority%l - self%ng*self%size/self%n
    self%interest%u  = self%authority%u + self%ng*self%size/self%n

    self%gn = self%n + 2*self%ng
    if (self%ng(3)==0) then
      call self%set(bits%no_download)
    endif

    call self%init_bdries()
    call self%set_initial_time()
    call self%append_to_task_list()
    call trace%end(itimer)
  END SUBROUTINE post_pre_init

  !===============================================================================
  recursive SUBROUTINE pre_init_and_link (self, m3d_parent, i_mu, i_nu, sp, atom)
    class(m3d_patch_t), target :: self, m3d_parent
    class(atomic_spectrum_type), pointer, optional :: sp
    class(atom_type), target, optional :: atom
    integer, optional :: i_mu, i_nu
    !-----------------------------------------------------------------------------
    if (present(i_mu)) self%i_mu = i_mu
    if (present(i_nu)) self%i_nu = i_nu
    if (present(sp)) self%sp => sp
    if (present(atom)) self%atom => atom
    self%composition => m3d_parent%composition
    self%patch => m3d_parent%patch
    self%m3d_parent => m3d_parent
    self%zs = m3d_parent%zs
    call self%pre_init()
    call m3d_parent%add_to_key_list(self)
  END SUBROUTINE pre_init_and_link

  !===============================================================================
  SUBROUTINE set_nv_and_nw (self, nv, nw)
    class(m3d_patch_t) :: self
    integer :: nv, nw
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_nv_and_nw')
    self%nv = nv
    self%nw = nw
    self%mv = self%nv
    call trace%end()
  END SUBROUTINE set_nv_and_nw
  
  !===============================================================================
  SUBROUTINE init (self)
    class(m3d_patch_t) :: self
    integer :: li(3), ui(3), status
    integer, save:: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%init', itimer=itimer)
    ! internal is set in task%pre_init, which is called by patch%init
    status = self%status
    call self%patch_t%init()
    call data_io%init (self)
    self%initialized=.true.
    self%status = status

    if (check_NaN) then
      self%mem = tiny(1.)
    endif
    ! ----------------------------------------------------------------------------
    ! Set unsigned and pervolume.
    ! Although quantities are generally per volume, but density is not stored
    ! locally, thus this flag must not be set.
    ! opacities can become negative in rare cases (masers and lasers)!
    ! ----------------------------------------------------------------------------
    self%pervolume = .false.
    self%unsigned  = .false.

    li = self%li
    ui = self%ui

    self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
    if (self%save_output) self%out_mem => self%authority_mem(:,:,:,:,self%it,:)

    self%pcheck = pcheck
    self%fcheck = fcheck
    self%linecheck = linecheck
    self%lvlcheck = lvlcheck
    
    self%pcheck(3) = self%pcheck(3) - self%zs
    self%pcheck = self%pcheck + self%ng

    where (pcheck < 0)
      self%pcheck = self%pcheck + nz_total + 1
    endwhere
    call trace%end(itimer)
  END SUBROUTINE init

  !===============================================================================
  SUBROUTINE set_initial_time (self)
    class(m3d_patch_t) :: self
    real(8) :: tnu
    integer :: thread
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_initial_time')
    !-----------------------------------------------------------------------------
    ! Synchronize with MHD patch
    ! ----------------------------------------------------------------------------
    !self%nt = 1
    !allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))

    !self%time     = self%patch%time
    !self%dt_fixed = self%patch%dtime
    !self%dtime    = self%patch%dtime
    !self%dt       = self%patch%dtime
    !self%out_next = self%patch%out_next
    self%t_offset = 0.5d0

    !if (self%i_mu == 0) then
    !  thread = mod(self%ipos(3),omp%nthreads)
    !  call self%set_place(self, thread2place(thread))
    !endif


    !select case (self%kind)
    !case ('composition', 'atom_patch')
    !thread = mod(self%ipos(3),omp%nthreads)
    !call self%set_place(self, thread2place(thread, omp%ncores, omp%nthreads))
    !end select


    select case (self%kind)
    !case ('composition', 'atom_patch')
    !  self%t_offset = 0.4d0
    !case ('back_cont')
      !if (save_back_cont .and. .not. self%out_ray) then
      !  self%t_offset = 0.4d0
      !endif
    case ('out_rad')
      self%t_offset = 0.9d0

      !if (self%boundaries%is_set(bits%zl + bits%zu)) then
      !  call self%set_place(self, self%place, push=.true.)
      !endif
      !if (self%ipos(3) < omp%nthreads) then
      !  call self%set_place(self, self%place, push=.true.)
      !endif
    end select

    tnu = abs(self%i_nu) / real(self%n_nu)
    self%t_offset = self%t_offset + 0.1*tnu

    if (self%boundaries%is_clear(bits%zl + bits%zu)) then
      self%t_offset = self%t_offset + 0.001
    endif

    if (self%out_ray) then
      self%time = self%time + ilambd + maxiter
    endif

    self%time = self%time + self%t_offset * self%dtime

    self%t = self%time            ! in case already set
    self%iit = 1                   ! initial slots
    self%it = 1                    ! slot for current values
    self%new = 1                   ! slot for new values
    self%iit(self%nt) = self%new   ! slot for new values
    self%time = self%t(self%it)    ! initial time

    select case (self%kind)
    case ('atom_patch', 'atom_cont')
      ! this is set back to 1 in update
      self%dtime = self%dtime + self%dtime*ilambd
    case ('atom_rate')
      self%time = self%time + self%dtime*ilambd
    end select
    call trace%end()
  END SUBROUTINE set_initial_time

  !===============================================================================
  SUBROUTINE set_maxiter (self)
    class (m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_maxiter')
    self%master_maxiter = maxiter + ilambd

    select case (self%kind)
    case ('composition')
      self%maxiter = 1
    case ('atom_line')
      if (save_atom_line) then
        self%maxiter = 0
      else
        self%maxiter = self%m3d_parent%maxiter - 1
      endif
    case ('atom_patch')
      self%maxiter = maxiter
    case ('atom_cont')
      self%maxiter = self%m3d_parent%maxiter-1
    case ('atom_rate')
      self%maxiter = self%m3d_parent%maxiter + ilambd - 1
      self%iter = ilambd
    case ('back_cont')
      if (save_back_cont .and. .not. line_blanketing) then
        self%maxiter = min(ilambd, maxiter - 1)
      else
        self%maxiter = ilambd + maxiter - 1
      endif
    case ('mean_rad')
      !self%maxiter = ilambd - 1
      !if (comoving_frame) then
      !  self%maxiter = self%maxiter + maxiter
      !endif
      self%maxiter = ilambd + maxiter - 1
    case ('angle_rad', 'net_opac', 'proj_patch')
      self%maxiter = ilambd + maxiter - 1
    end select

    if (self%out_ray) self%maxiter = 0
    call trace%end()
  END SUBROUTINE set_maxiter
  
  !===============================================================================
  !> Initialize boundary bits
  !===============================================================================
  SUBROUTINE init_bdries (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    if (.not. associated (self%boundaries)) then
      allocate (self%boundaries)
      call self%boundaries%pre_init (self)
    endif
    if (self%out_ray .and. maxiter==0) then
      call self%boundaries%clear(bits%zu)
    endif
  END SUBROUTINE init_bdries

  !===============================================================================
  SUBROUTINE pre_update (self)
    class(m3d_patch_t), target :: self
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_update', itimer=itimer)
    call lock_mem (self, 'update')
    if (print_start) then
      write(stdout, *) 'Update on', self%id, 'has begun', self%finished, self%time
    endif
    if (.not. self%is_finished()) then
      if (converged .and. .not. self%out_ray) then
        call self%finish()
        call self%free_unused_memory()
        self%print_update_info = .false.
        call print_execution_history(self, 'SKP')
      else
        call print_execution_history(self, 'PRE')
        call self%check_mem(locked=.true.)
      endif
    endif
    if (self%verbose > 1) then
      if (.not. self%checked_bounds) then
        self%print_update_info = self%check_in_bounds()
        self%checked_bounds = .true.
      endif
    endif
    self%place = 0
    call trace%end(itimer)
  END SUBROUTINE pre_update

  !===============================================================================
  SUBROUTINE update (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    if (.not. self%is_finished()) then
      call self%check_access(self)
      call self%execute()
      call self%post_update()
    endif
    call unlock_mem (self, 'update')
  END SUBROUTINE update

  !===============================================================================
  SUBROUTINE execute (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    ! fallback
  END SUBROUTINE execute
  
  !===============================================================================
  SUBROUTINE post_update (self)
    class(m3d_patch_t) :: self
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%post_update', itimer=itimer)
    if (self%kind == 'atom_patch' .and. converged) then
      call print_execution_history(self, 'CON')
    else
      call print_execution_history(self, 'END')
    endif
    if (self%verbose > 1) then
      if (self%print_update_info) then
        call self%print_post_update_info()
      elseif (self%verbose > 2) then
        call self%print_to_stdout()
      endif
    endif
    if (check_NaN) then
      call self%check_mem_NaN()
    endif
    if (self%save_output .and. self%finished) then
      call self%save(self%out_mem, self%kind)
    endif
    call trace%end(itimer)
  END SUBROUTINE post_update

  !===============================================================================
  SUBROUTINE rotate (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    self%iter = self%iter + 1
    if (self%iter > self%maxiter .or. self%nbors_finished()) then
      call self%finish()
    endif
    call self%task_t%rotate
    call self%free_unused_memory()
    !call self%claim_nbor_place()
  END SUBROUTINE rotate

  !===============================================================================
  FUNCTION m3d_patch_info (self) result(text)
    class(m3d_patch_t) :: self
    character(len=merge(48, 32, self%verbose>3)) :: text
    integer :: i
    !-----------------------------------------------------------------------------
    i = omp%place
    if (self%verbose > 3) then
      write(text, '(I3, 1x, I3, 1x, F8.4, A12, I8, 3(I4))') i, self%place, self%time, &
                      trim(self%kind), self%id, self%ipos(3)+1, self%i_nu, self%i_mu
    else
      write(text, '(A12, I8, 3(I4))')  trim(self%kind), self%id, &
      self%ipos(3)+1, self%i_nu, self%i_mu
    endif
  END FUNCTION m3d_patch_info

  !===============================================================================
  FUNCTION m3d_patch_info2 (self) result(text)
    class(m3d_patch_t) :: self
    character(len=38) :: text
    !-----------------------------------------------------------------------------
    write(text, '(A12, I8)') trim(self%kind), self%id
  END FUNCTION m3d_patch_info2
  
  !===============================================================================
  SUBROUTINE finish (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%finish')
    self%time = io%end_time
    self%t = io%end_time
    call self%set_place(self, -2)
    call self%set(bits%no_download)
    call self%set(bits%ignore)
    if (self%save_output) call self%save(self%out_mem, self%kind)
    call trace%end()
  END SUBROUTINE finish

  !===============================================================================
  LOGICAL FUNCTION is_finished (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    is_finished = self%time >= io%end_time
  END FUNCTION is_finished
  
  !===============================================================================
  SUBROUTINE free_unused_memory (self)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbor, nbnbor
    class(task_t), pointer:: nbtask, nbnbtask
    logical :: mem_not_needed
    !-----------------------------------------------------------------------------
    if (self%free_self) then
      call self%deallocate_mem(self, print=.true.)
    endif
    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed .and. nbor%task%is_clear(bits%external)) then
        nbtask => nbor%task
        select type (nbtask)
        class is (m3d_patch_t)
          if (.not. nbtask%static_mem) then
            nbnbor => nbtask%link%nbor
            mem_not_needed=.true.
            do while (associated(nbnbor) .and. mem_not_needed)
              nbnbtask => nbnbor%task
              if (nbnbtask%is_clear(bits%virtual)) then
                if (nbnbor%needs_me) then
                  if (nbnbtask%key(1)==nbtask%key(1) .or. nbnbtask%is_clear(bits%no_download)) then
                    if (nbnbtask%time < nbtask%time) then
                      mem_not_needed=.false.
                    endif
                  endif
                endif
              endif
              nbnbor => nbnbor%next
            end do
            if (mem_not_needed .and. nbtask%initialized) then
              call lock_mem(nbtask, 'free_unused_memory')
              call nbtask%deallocate_mem(self, print=.true.)
              call unlock_mem(nbtask, 'free_unused_memory')
            endif
          endif
        end select
      endif
      nbor => nbor%next
    end do
  END SUBROUTINE free_unused_memory
  
  !===============================================================================
  FUNCTION nbors_finished (self)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbor, nbnbor
    class(task_t), pointer:: nbtask, nbnbtask
    logical :: nbors_finished
    !-----------------------------------------------------------------------------
    nbors_finished=.true.
    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed) then
        nbtask => nbor%task
        select type (nbtask)
        class is (m3d_patch_t)
          if (.not. nbtask%finished) then
            nbors_finished=.false.
          endif
        end select
      endif
      nbor => nbor%next
    end do
  END FUNCTION nbors_finished

  !===============================================================================
  SUBROUTINE m3d_nbor_relations(self, another, needed, needs_me, download)
    class(m3d_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    class(task_t), pointer :: other
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_nbor_relations')

    if (.not.associated(another)) then
      write (stderr,*) mpi%rank, 'another not associated'
      call io%abort ('task_t another not associated')
    end if

    needed=.true. ; needs_me=.true. ; download=any(self%ng>0)

    if (another%is_set (bits%frozen)) then
      needs_me = .false.
      needed = .false.
    endif

    if (self%rank /= mpi%rank) then
      !---------------------------------------------------------------------------
      ! virtual tasks don't need guard zones
      !---------------------------------------------------------------------------
      download = .false.
      if (another%rank /= mpi%rank) then
        !---------------------------------------------------------------------------
        ! task is virtual, nbor is virtual
        !---------------------------------------------------------------------------
        needed = .false.
        needs_me = .false.
      endif
    endif

    if (needed .or. needs_me) then
      select type(another)
      class is(m3d_patch_t)
        if (another%kind == self%kind) then
          if (another%patch%id == self%patch%id) then
            needed=.false. ; needs_me=.false. ; download=.false.
          elseif((another%i_mu /= self%i_mu) .or. (another%i_nu /= self%i_nu)) then
            needed=.false. ; needs_me=.false. ; download=.false.
          endif
        else
          download=.false.
          if (another%patch%id /= self%patch%id) then
            needed=.false. ; needs_me=.false.
          else
            if(self%i_mu /= 0 .and. another%i_mu /= 0) then
              if (another%i_mu /= self%i_mu) then
                needed=.false. ; needs_me=.false. ; download=.false.
              endif
            endif

            if(self%i_nu /= 0 .and. another%i_nu /= 0) then
              if (another%i_nu /= self%i_nu) then
                needed=.false. ; needs_me=.false. ; download=.false.
              endif
            endif
          endif
        endif
      class default
        if (another%id == self%patch%id) then
          needs_me=.false. ; download=.false.
        else
          needed=.false. ; needs_me=.false. ; download=.false.
        endif
      end select
    endif
  
    call trace%end()
  END SUBROUTINE m3d_nbor_relations
  
  !===============================================================================
  SUBROUTINE append_to_task_list(self)
    class (m3d_patch_t), target :: self
    class (m3d_halo_t), pointer :: link
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%add_to_task_list', itimer=itimer)
    if (associated(self%task_list)) then
      if (self%maxiter >= 0) then
        allocate (link)
        call os%mem_register(storage_size(link), 1, label='link')
        call link%init()
        link%task => self
        self%link => link
        !call link%append_to(self%task_list)
        call self%task_list%append_link(self%link)
        call self%task_list%hash_table%add (self%key, self%link)
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE append_to_task_list
  
  !===============================================================================
  !> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
  !===============================================================================
  SUBROUTINE task_info (self, unit)
    class(m3d_patch_t):: self
    integer, optional:: unit
    integer:: unit_L
    !-----------------------------------------------------------------------------
    unit_L = shared%option (unit, io_unit%mpi)
    write(unit_L,'(1x,a,2x,a10,i6,4i4,7x,2l1,2x,3f9.3)') &
      'm3d_patch_t%task_info: kind, id, rank, level, nu, mu, BV, pos', &
      self%kind, self%id, self%rank, self%level, self%i_nu, self%i_mu, &
      self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
  END SUBROUTINE task_info
  
  !===============================================================================
  SUBROUTINE init_level (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    self%level = self%patch%level
    !self%level = 1
  END SUBROUTINE init_level
  
  !===============================================================================
  SUBROUTINE periodic_grid (self, only)
    class(m3d_patch_t):: self
    integer, optional:: only
    integer :: iv, it, iw, i,j,k,t, ks, ke
    integer, save :: itimer=0
    !----------------------------------------------------------------------------
    if (.not. self%is_periodic()) return
    if (all(self%ng(:2) == 0) .or. .not. associated(self%mem)) return
    call trace%begin ('m3d_patch_t%periodic_grid', itimer=itimer)

    if (.not. associated(self%periodic_idx)) then
      call self%init_periodic_idx()
    endif

    if (present(only)) then
      ks = only
      ke = only
    else
      ks = 1
      ke = self%gn(3)
    endif

    do iw=1, self%nw
    do it=1, self%nt
    do iv=1, self%nv
    do k=ks, ke
      do t=1, size(self%periodic_idx, 2)
        i = self%periodic_idx(1,t)
        j = self%periodic_idx(2,t)
        self%mem(i,j,k,iv,it,iw) = self%mem(self%periodic_i(i),self%periodic_j(j),k,iv,it,iw)
      end do
    end do
    end do
    end do
    end do
    call trace%end(itimer)
  END SUBROUTINE periodic_grid 

  !===============================================================================
  SUBROUTINE init_periodic_idx (self)
    class(m3d_patch_t):: self
    integer, dimension(3) :: nn
    integer, dimension(:), allocatable :: ii, jj
    integer :: i, j, k, t, nidx, li(3)
    !----------------------------------------------------------------------------
    call os%allocate(self%periodic_i, self%gn(1))
    call os%allocate(self%periodic_j, self%gn(2))

    li = self%ng + 1

    do i=1, self%gn(1)
      self%periodic_i(i) = modulo(i-li(1),self%n(1)) + li(1)
    end do
    do j=1, self%gn(2)
      self%periodic_j(j) = modulo(j-li(2),self%n(2)) + li(2)
    end do

    nidx = product(self%gn(1:2)) - product(self%n(1:2))

    nn = self%ng + self%n
    call os%allocate(self%periodic_idx, [2, nidx])
    t = 1
    do j=1, self%gn(2)
    do i=1, self%gn(1)
      if ((i <= self%ng(1) .or. i > nn(1)) .or. (j <= self%ng(2) .or. j > nn(2))) then
        self%periodic_idx(1, t) = i
        self%periodic_idx(2, t) = j
        t = t + 1
      endif
    end do
    end do
  END SUBROUTINE init_periodic_idx
  
  !===============================================================================
  character(len=17) FUNCTION m3d_name (self)
    class(m3d_patch_t):: self
    !----------------------------------------------------------------------------
    write(m3d_name, '(A12,I5)') trim(self%kind), self%id
  END FUNCTION
  
  !===============================================================================
  SUBROUTINE extrapolate_guards (self, target_mem)
    class(m3d_patch_t):: self
    real(4), pointer, contiguous, optional:: target_mem(:,:,:,:)
    !.............................................................................
    ! do nothing
  END SUBROUTINE
  
  !===============================================================================
  !> Check if self (which is an nbor task) is ahead of target (which is the one to
  !> possibly move to the ready queue), using self%dtime*target%grace as the grace
  !> period, since we want to limit the extrapolation in the nbor task to at most
  !> target%grace*self%dtime.
  !===============================================================================
  LOGICAL FUNCTION is_ahead_of (self, target)
    class (m3d_patch_t) :: self
    class (task_t), pointer :: target
    !.............................................................................
    real:: grace
    real(8):: nbtime, nbdtime, tgtime
    logical:: needed, needs_me, download
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%is_ahead_of')
    grace = self%grace + shared%extra_grace
                                            call self%lock%set ('is_ahead_of')
    !---------------------------------------------------------------------------
    ! dt(self%dt) is always consistent with t(:), while dtime changes in
    ! courant_condition(), before other attributes change in rotate()
    !---------------------------------------------------------------------------
    nbtime = self%time
    tgtime = target%time
    is_ahead_of = .false.
    if (self%is_set(bits%frozen)) then
      is_ahead_of = .true.
    elseif (nbtime == tgtime) then
      call self%nbor_relations(target, needed, needs_me, download)
      if (needed .eqv. needs_me) then
        select type(target)
        class is (m3d_patch_t)
          if (associated(target%m3d_parent)) then
            is_ahead_of = .not. target%m3d_parent%id == self%id
          else
            is_ahead_of = .true.
          endif
        end select
      endif
    elseif (nbtime > tgtime) then
      is_ahead_of = .true.
    end if
                                            call self%lock%unset ('is_ahead_of')
    call trace%end()
  END FUNCTION is_ahead_of

  !===============================================================================
  !> Disambiguation: distance from self to task -- applies (only) to Cartesian
  !> meshes
  !===============================================================================
  FUNCTION distance_to (self, task) RESULT (d)
    class(m3d_patch_t):: self
    class(task_t):: task
    real(8):: d(3), p(3)
    !-----------------------------------------------------------------------------
    select type (task)
    class is (m3d_patch_t)
      if (self%i_mu /= 0 .and. task%i_mu /= 0) then
        if (self%i_mu /= task%i_mu) then
          d = huge(1.0)
          return
        endif
      endif
      if (self%i_nu /= 0 .and. task%i_nu /= 0) then
        if (self%i_nu /= task%i_nu) then
          d = huge(1.0)
          return
        endif
      endif
    end select
    d = self%task_t%distance_to(task)
  END FUNCTION distance_to

  !===============================================================================
  SUBROUTINE allocate_mem (self, allocator, m_in, label)
    class(m3d_patch_t):: self
    class(m3d_patch_t), intent(in), target :: allocator
    integer, intent(in), optional :: m_in(6)
    character(len=*), optional :: label
    integer :: m(6), li(3), ui(3)
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    if (.not. associated(self%mem)) then
      call trace%begin ('m3d_patch_t%allocate_mem', itimer=itimer)
      if (present(m_in)) then
        m = m_in
      else
        m = [self%gn, self%nv, self%nt, self%nw]
        self%compressed = .false.
      endif
      if (present(label)) then
        self%mem_label = label
      else
        self%mem_label = 'mem'
      endif
      call os%allocate(self%mem, m, label=trim(self%kind)//'%'//self%mem_label)
      if (check_NaN) then
        self%mem = tiny(1.0)
      endif
      self%nullified = .false.
      li = self%li
      ui = self%ui
      if (self%compressed) then
        li(3) = 1
        ui(3) = 1 
      endif
      self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
      if (self%save_output) then
        self%out_mem => self%authority_mem(:,:,:,:,self%it,:)
      endif
      self%last_allocator => allocator
      self%allocation_time = wallclock()
      call trace%end (itimer)
    endif
  END SUBROUTINE allocate_mem

  !===============================================================================
  SUBROUTINE deallocate_mem (self, deallocator, print)
    class(m3d_patch_t):: self
    class(m3d_patch_t), intent(in), target :: deallocator
    logical, optional :: print
    integer :: m(6), li(3), ui(3)
    logical :: free
    !-----------------------------------------------------------------------------
    if (associated(self%mem)) then
      if (self%is_clear(bits%busy) .or. deallocator%id==self%id) then
        call os%deallocate(self%mem, label=trim(self%kind)//'%'//self%mem_label)
        self%nullified = .false.
        self%last_deallocator => deallocator
        self%deallocation_time = wallclock()
        if (present(print)) then
          if (print .and. (self%print_update_info .or. self%verbose > 2) .and. detailed_prints) then
            write(stdout, *) '  deallocated', self%m3d_patch_info2()
          endif
        endif
      endif
    endif
  END SUBROUTINE deallocate_mem

  !===============================================================================
  SUBROUTINE check_mem (self, locked)
    class(m3d_patch_t), target :: self
    logical, intent(in), optional :: locked
    integer :: orig_place
    logical :: push
    !-----------------------------------------------------------------------------
    if (.not. associated(self%mem)) then
      if (.not. self%initialized) then
        if (self%is_set(bits%virtual)) then
          self%kind = trim(self%kind)//'V'
          self%static_mem = .true.
        endif
        !push = self%place >= 100
        orig_place = self%place
        call self%init()
        self%place = orig_place
        !if (push) self%place = self%place + 100
        !if (mod(orig_place,100) /= mod(self%place,100)) then
        !  print *, self%m3d_patch_info2(), ' changed place from ', orig_place
        !endif
        self%last_allocator => self
        self%allocation_time = wallclock()
      else
        call self%allocate_mem(self)
      endif
      !if (self%place /= omp%place .and. self%i_mu > 0) then
      !  self%place = omp%place ! memory affinity
      !endif
    elseif (self%compressed) then
      call self%decompress_mem(self, locked=option(locked, .true.))
    endif
  END SUBROUTINE check_mem

  !===============================================================================
  SUBROUTINE compress_mem (self, compressor, static, locked)
    class(m3d_patch_t):: self
    class(m3d_patch_t), intent(in), target :: compressor
    logical, intent(in) :: static, locked
    integer :: i, k, m(6)
    real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: dn
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    if (.not. self%compressed) then
      call trace%begin ('ray_patch_t%compress_mem', itimer=itimer)
      if (.not. locked) call lock_mem(self, 'compress_mem')
      self%static_mem = static
      self%compressed = .true.
      m = shape(self%mem)
      m(3) = 1
      k = merge(self%ui(3), self%li(3),self%dir>0)
      call os%allocate(dn, m)
      do i=1, self%nv
        dn(:,:,1:1,i,:,:) = self%mem(:,:,k:k,i,:,:)
      end do
      call self%deallocate_mem(compressor)
      call self%allocate_mem(compressor, m, 'zip')
      self%allocation_time = self%time
      do i=1, self%nv
        self%mem(:,:,1:1,i,:,:) = dn(:,:,1:1,i,:,:)
      end do
      call os%deallocate(dn)
      if (.not. locked) call unlock_mem (self, 'compress_mem')
      if (self%print_update_info .and. detailed_prints) then
        write(stdout,*) self%m3d_patch_info2(), ' has been compressed'
      endif
      call trace%end(itimer)
    endif
  END SUBROUTINE compress_mem

  !===============================================================================
  SUBROUTINE decompress_mem (self, decompressor, locked)
    class(m3d_patch_t) :: self
    class(m3d_patch_t), intent(in), target :: decompressor
    logical, intent(in), optional :: locked
    integer :: i, k, m(6)
    real(4), dimension(:,:,:,:,:,:), pointer, contiguous :: up
    integer, save:: itimer=0
    !-----------------------------------------------------------------------------
    if (self%compressed) then
      call trace%begin ('ray_patch_t%decompress_mem', itimer=itimer)
      if (.not. locked) call lock_mem(self, 'decompress_mem')
      self%static_mem = .false.
      self%compressed = .false.
      m = shape(self%mem)
      call os%allocate(up, m)
      do i=1, self%nv
        up(:,:,1:1,i,:,:) = self%mem(:,:,1:1,i,:,:)
      end do
      call self%deallocate_mem(decompressor)
      call self%allocate_mem(decompressor)
      k = merge(self%ui(3), self%li(3),self%dir>0)
      do i=1, self%nv
        self%mem(:,:,k:k,i,:,:) = up(:,:,1:1,i,:,:)
      end do
      call os%deallocate(up)
      if (.not. locked) call unlock_mem (self, 'decompress_mem')
      if (self%print_update_info .and. detailed_prints) then
        write(stdout,*) self%m3d_patch_info2(), ' has been decompressed'
      endif
      call trace%end(itimer)
    endif
  END SUBROUTINE decompress_mem

  !===============================================================================
  SUBROUTINE check_mem_NaN (self)
    class(m3d_patch_t):: self
    integer :: i,j,k,v,t,w
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%check_mem_NaN', itimer=itimer)
    if (associated(self%mem)) then
      do w=1, self%nw
      do t=1, self%nt
      do v=1, self%nv
      do k=1, self%gn(3)
      do j=1, self%gn(2)
      do i=1, self%gn(1)
        if (self%mem(i,j,k,v,t,w) /= self%mem(i,j,k,v,t,w)) then
          write(stdout, '(A,A,3(I4),I5,I3,I3)') self%m3d_patch_info(), '  NaN', i,j,k,v,t,w
          stop
        elseif (self%mem(i,j,k,v,t,w) > inf4) then
          write(stdout, '(A,A,3(I4),I5,I3,I3)') self%m3d_patch_info(), '  Infinity', i,j,k,v,t,w
          stop
        !elseif (self%unsigned(w) .and. self%mem(i,j,k,v,t,w) < 0) then
        !  write(stdout, '(A,A,3(I4),I5,I3,I3)') self%m3d_patch_info(), '  negative', i,j,k,v,t,w
        !  stop
        endif
      end do
      end do
      end do
      end do
      end do
      end do
    endif
    call trace%end(itimer)
  END SUBROUTINE check_mem_NaN

  !===============================================================================
  SUBROUTINE claim_nbor_place (self)
    class(m3d_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    integer, save :: itimer=0
    integer :: key(4), topscore
    logical :: pushed, need
    !-----------------------------------------------------------------------------
    pushed = .false.
    key = self%key
    nbor => self%link%nbor
    do while (associated (nbor))
      nbtask => nbor%task
      if (nbor%needs_me .and. (.not. nbtask%finished)) then
        if (count(key(1:3)==nbtask%key(1:3)) == 3) then
          if(is_ready(self, nbor%link)) then
            select type (nbtask)
            class is (m3d_patch_t)
              !nbtask%prepend = .true.
              call self%push(nbtask)
              pushed = .true.
            end select
          endif
        endif
      endif
      nbor => nbor%next
    end do

  contains
    !-----------------------------------------------------------------------------
    function is_ready(self, link) result(ok)
      class(task_t) :: self
      class(link_t) :: link
      class(nbor_t), pointer :: nbor
      class(task_t), pointer :: task, nbtask
      logical :: ok
      integer :: ignore
      !-----------------------------------------------------------------------------
      ok = .true.
      ignore = bits%frozen + bits%remove + bits%external
      task => link%task
      nbor => link%nbor
      do while (associated (nbor))                        ! keep going until end
        nbtask => nbor%task
        if (nbor%needed .and. nbtask%is_clear(ignore)) then
          if (abs(task%level-nbtask%level) <= 1) then
            if (.not. nbor%task%is_ahead_of (task)) then
              ok = .false.
            endif
          endif
        endif
        nbor => nbor%next    
      end do
    end function is_ready

    !-----------------------------------------------------------------------------
    subroutine key_test (task, choice, topscore, key)
      class(task_t), target :: task
      class(task_t), pointer :: choice
      integer :: key(:), topscore
      integer :: i, score
      !-----------------------------------------------------------------------------
      score = 0
      if(key(1) == task%key(1)) then
        score = count(key == task%key(1:3))
        if (score == 3) then
          score = 10
        elseif (key(2)==0) then
          score = count(task%key(1:3) /= 0)
        endif
      elseif (any(key < 0)) then
        score = count(key == task%key(1:3))
        score = merge(score, 0, all(key(2:3) == task%key(2:3)))
      endif
      if (score > topscore) then
        topscore = score
        choice => task
      endif
      !print '(3x, I3, I5, 1x, A)', score, task%id, trim(task%kind)
    end subroutine key_test
  END SUBROUTINE claim_nbor_place

  !===============================================================================
  SUBROUTINE set_place (self, other, place, push)
    class(m3d_patch_t) :: self
    class(task_t):: other
    integer :: place
    logical, optional :: push
    !-----------------------------------------------------------------------------
    !if (option(push) .and. place < 100) then
    !  other%place = place+100
    !else
      other%place = place
    !endif
    if (self%id /= other%id .and. print_push) then
      print '(A, A, A, A, I4, A, I4)', trim(self%m3d_patch_info2()), ' set place of ', &
      trim(other%kind), ' ', other%id, ' to ', other%place
    endif
  END SUBROUTINE set_place

  !===============================================================================
  SUBROUTINE push (self, other)
    class(m3d_patch_t) :: self, other
    !-----------------------------------------------------------------------------
    if (other%is_set(bits%ready)) return
    other%t_push = t_push
    other%time = other%time - other%t_push
    if (self%id /= other%id .and. print_push) then
      print '(A, A, A, 1x, I6)', trim(self%m3d_patch_info2()), ' pushed ', &
      trim(other%kind), other%id
    endif
  END SUBROUTINE push
  
  !===============================================================================
  SUBROUTINE skip_iterations (self, n)
    class(m3d_patch_t):: self
    integer, intent(in) :: n
    !-----------------------------------------------------------------------------
    self%time = self%time + n * self%dtime
    self%iter = self%iter + n

    if (self%iter > self%maxiter) then
      call self%finish()
    endif
  END SUBROUTINE skip_iterations

  !===============================================================================
  SUBROUTINE check_access (self, other)
    class(m3d_patch_t):: self
    class(task_t):: other
    !-----------------------------------------------------------------------------
    select type(other)
    class is (m3d_patch_t)
      if (.not. associated(other%mem)) then
        call io%lock%set('check_access')
        print *, 'memory of ', other%m3d_patch_info()
        print *, 'cannot be accessed from ', self%m3d_patch_info()
        print *, 'static mem =', other%static_mem
        if (associated(other%last_deallocator)) then
          print '(a, a, F10.6)', 'it was last deallocated by', other%last_deallocator%m3d_patch_info2(), other%deallocation_time
        endif
        if (associated(other%last_allocator)) then
          print '(a, a, F10.6)', 'it was last   allocated by', other%last_allocator%m3d_patch_info2(), other%allocation_time
        endif
        call io%lock%unset('check_access')
        call shared%abort('memory array not accessible')
      endif
    end select
  END SUBROUTINE check_access

  !===============================================================================
  SUBROUTINE print_to_stdout (self, arr, text)
    class(m3d_patch_t) :: self
    real(4), dimension(:), optional :: arr
    character(len=*), optional :: text
    !-----------------------------------------------------------------------------
    !$omp critical (cr_print_to_stdout)
    if (present(arr)) then
      if (size(arr)==0) then
        write(stdout,'(A, A)') self%m3d_patch_info(), 'invalid index!'
      else
        write(stdout,'(A, 6(ES15.7))') self%m3d_patch_info(), arr
      endif
    elseif (present(text)) then
      write(stdout,'(A, A, A)') self%m3d_patch_info(), ' ', trim(text)
    else
      write(stdout,'(A)') self%m3d_patch_info()
    endif
    !$omp end critical (cr_print_to_stdout)
  END SUBROUTINE print_to_stdout

  !===============================================================================
  SUBROUTINE print_post_update_info (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    call self%print_to_stdout()
  END SUBROUTINE print_post_update_info

  !===============================================================================
  FUNCTION check_in_bounds (self) result (cib)
    class(m3d_patch_t) :: self
    logical :: cib
    !-----------------------------------------------------------------------------
    cib = all(self%pcheck > 0) .and. all(self%pcheck <= self%gn)
  END FUNCTION check_in_bounds

  !===============================================================================
  SUBROUTINE print_execution_history (self, status)
    class(task_t) :: self
    character(len=3) :: status
    !-----------------------------------------------------------------------------
    !$omp critical (cr_print_execution_history)
    write(io_unit%dispatcher, &
    '(F17.8, 2(I3, 1x), 1x, A3, 1x, A12, 1x, I8, 1x, 3(I3,1x), F17.8)') &
    wallclock(), mpi%rank, omp%thread, status, trim(self%kind), self%id, self%key(1:3), self%time
    !$omp end critical (cr_print_execution_history)
  END SUBROUTINE print_execution_history

  !===============================================================================
  !> Return a unique hash key for the patch
  !===============================================================================
  FUNCTION hash_key (self, position, size, level, box, origin) RESULT (key)
    class(m3d_patch_t):: self
    real(8), dimension(3), optional:: position, origin, size, box
    integer, optional:: level
    integer:: key(4)
    !-----------------------------------------------------------------------------
    key(1) = self%ipos(3) + 1 
    key(2) = self%i_nu
    key(3) = self%i_mu
    key(4) = hash_string(self%kind)
  END FUNCTION hash_key

  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
    class(key_link_t), pointer :: klink
    integer, save:: itimer=0
    integer :: i, i_nu, key(4)
    !-----------------------------------------------------------------------------
    key(1) = self%key(1)
    key(2:3) = 0
    key(4) = hash_string(self%patch%kind)
    call self%check_nbor(key, nbors)

    if (associated(self%m3d_parent)) then
      key(4) = hash_string(self%m3d_parent%kind)
      call self%check_nbor(key, nbors)
    endif

    klink => self%key_list
    do while (associated(klink))
      call self%check_nbor(klink%key, nbors)
      klink => klink%next
    end do
  END SUBROUTINE init_nbors_hash

  !===============================================================================
  SUBROUTINE check_nbor (self, key, nbors)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbors
    class(link_t), pointer:: nblink
    integer:: key(4)
    logical:: added
    !-----------------------------------------------------------------------------
    call self%task_list%hash_table%get (key, nblink)
    if (associated(nblink)) then
      select type (link => self%link)
      class is (halo_t)
        call link%conditionally_add_nbor (nblink, nbors, added)
        if (added) then
          !$omp atomic update
          self%n_nbors = self%n_nbors + 1
        end if
      end select
    end if
  END SUBROUTINE check_nbor

  !===============================================================================
  SUBROUTINE save (self, mem, label)
    class(m3d_patch_t):: self
    real(4), dimension(:,:,:,:,:) :: mem 
    character(len=*) :: label
    character(len=128) :: meta_file, folder
    character(len=32) :: filename
    integer :: lunit, iw, iv
    integer, save :: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin('m3d_patch_t%save', itimer=itimer)
    if (size(mem) == 0) then
      call shared%abort('mem not associated in save')
    endif
    write(folder, '(A, A)') trim(io%run_dir), 'save/'
    if (self%out_ray) then
      write(filename, '(A, A, 3(I3.3), A)') trim(label), '_out_', abs(self%key(1:3)), '.bin'
    elseif (.not. converged .or. self%kind=='atom_patch') then
      write(filename, '(A, A, 3(I3.3), A)') trim(label), '_', abs(self%key(1:3)), '.bin'
    else
      call trace%end(itimer)
      return
    endif
    meta_file = trim(label) // '_meta.txt'
    call io%lock%set('save')
    if (os%exists(trim(folder)//meta_file)) then
      open(newunit=lunit, file=trim(folder)//meta_file, form='formatted', position='append')
    else
      call os%mkdir(folder)
      open(newunit=lunit, file=trim(folder)//meta_file, form='formatted')
      write(lunit, '(A, A8, A22, A5, 8(A10))') '*', ' filename', '', 'bytes', &
      'nx', 'ny', 'nz', 'nv', 'nw'
    endif
    if (.not. self%saved) write(lunit, '(A32, I4, 8(I10))') filename, storage_size(mem), shape(mem)
    close(lunit)
    open(newunit=lunit, file=trim(folder)//filename, access='stream', status='replace')
    write(lunit) mem
    close(lunit)
    call io%lock%unset('save')
    self%saved = .true.
    call trace%end(itimer)
  END SUBROUTINE save

  !===============================================================================
  SUBROUTINE load_mem (self, folder, mem, label)
    class(m3d_patch_t), intent(inout) :: self
    character(len=128), intent(in) :: folder
    character(len=*), intent(in) :: label
    real(4), intent(inout), dimension(:,:,:,:,:) :: mem
    character(len=128) :: filename, sfolder, metafile, file
    integer :: lunit, file_size, mem_size, iostat, prec, oldm(5), newm(5)
    ! ----------------------------------------------------------------------------
    call self%check_access(self)
    write(sfolder, '(A, A)') trim(folder), 'save/'
    write(filename, '(A, A, 3(I3.3), A)') trim(label), '_', self%key(1:3), '.bin'
    write(metafile, '(A, A, 3(I3.3), A)') trim(label), '_meta.txt'
    call os%check_file(trim(sfolder)//filename)
    call os%check_file(trim(sfolder)//metafile)
    open(newunit=lunit, file=trim(sfolder)//metafile, form='formatted', iostat=iostat)
    read(lunit,*) ! skip header
    do while (iostat==0)
      read(lunit, *, iostat=iostat) file, prec, oldm
      if (file == filename) exit
    end do
    close(lunit)
    if (iostat /= 0) then
      print *, 'file not found in metafile! Stopping...'
      stop
    endif
    newm = shape(mem)
    open(newunit=lunit, file=trim(sfolder)//filename, access='stream')
    if (all(oldm == newm)) then
      read(lunit) mem
    else
      write(stdout, *) 'Cannot load patch save ', trim(sfolder)//trim(filename), &
      'Memory sizes do not match!'
      print *, oldm
      print *, newm
      stop
    endif
    close(lunit)
  END SUBROUTINE load_mem

  !===============================================================================
  SUBROUTINE add_to_key_list (self, other)
    class(m3d_patch_t), intent(inout) :: self, other
    class(key_link_t), pointer :: next, link, prev
    ! ----------------------------------------------------------------------------
    nullify(prev)
    next => self%key_list
    do while(associated(next))
      prev => next
      next => next%next
    end do
    allocate(link)
    link%key = other%key
    if (associated(prev)) then
      prev%next => link
    else
      self%key_list => link
    endif
  END SUBROUTINE add_to_key_list

  !===============================================================================
  FUNCTION queue_before (self, other) result(queue)
    class(m3d_patch_t):: self
    class(task_t):: other
    logical :: queue
    !-----------------------------------------------------------------------------
    if (self%prepend) then
      queue = other%time >= self%time
      !self%prepend = .false.
    else
      queue = self%task_t%queue_before(other)
    endif
    if (self%t_push > 0.) then
      self%time = self%time + self%t_push
      self%t_push = 0.
    endif
  END FUNCTION queue_before

  !===============================================================================
  SUBROUTINE sort_nbors (self)
    class(m3d_patch_t) :: self
    class(nbor_t), pointer:: nbor
    class(task_t), pointer:: nbtask
    integer:: n_nbors
    !-----------------------------------------------------------------------------
    nbor => self%link%nbor
    do while(associated(nbor))
      nbtask => nbor%link%task
      nbtask%istep = nbtask%level
      nbtask%level = -count(nbtask%key == self%key)
      nbor => nbor%next
    end do  

    nbor => self%link%nbor
    call sort_nbors_by_level(self%link, nbor, nbor)

    do while(associated(nbor))
      nbtask => nbor%link%task
      nbtask%level = nbtask%istep
      nbtask%istep = 0
      nbor => nbor%next
    end do  
  END SUBROUTINE sort_nbors

END MODULE m3d_patch_mod