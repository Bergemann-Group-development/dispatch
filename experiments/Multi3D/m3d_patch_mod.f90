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
    USE prec_mod, only: inf4
    !USE units_mod
    !USE scratch_mod
    USE m3d_halo_mod, only: m3d_halo_t, print_finished
    USE atom_mod, only: atomic_spectrum_type, lse_no_pivot
    USE ray_mod, only: set_carlson_quadratures, set_lobatto_quadratures, set_custom_quadratures
    implicit none
    private
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(gpatch_t) :: m3d_patch_t
      class(gpatch_t),      pointer :: patch => null()         ! pointer to mhd patch
      class(m3d_patch_t),   pointer :: m3d_parent => null()    ! pointer to parent patch
      class(atomic_spectrum_type), pointer :: sp => null()
      class(m3d_patch_t),   pointer :: last_allocator, last_deallocator
      logical                  :: on=.false.       ! allows local turn-off
      logical                  :: in_task_list = .false.
      logical                  :: initialized  = .false.
      real(4), dimension(:,:,:,:,:,:), pointer :: authority_mem => null()
      real(4), dimension(:,:,:,:,:), pointer :: out_mem => null() 
      integer, dimension(:,:), pointer :: periodic_idx
      integer, dimension(:), pointer :: periodic_i, periodic_j
      real(8) :: t_offset=0, t_push=0
      real(8) :: allocation_time, deallocation_time
      integer :: iter=0, maxiter, master_maxiter, dir=0
      integer :: n_mu, n_nu, n_omega, n_out, zs, ze
      integer :: fcheck, pcheck(3), linecheck, lvlcheck
      logical :: out_ray=.false., static_mem=.false.
      logical :: compressed=.false., save_output=.false.
      logical :: print_update_info
      integer :: i_nu=0, i_mu=0
      integer :: i_ne=1, i_rho=2, i_temp=3, i_nh1=4, i_nh6=9
      integer :: i_chi=1, i_eta=2, i_snu=2
      integer :: i_Gij=1, i_Gji=2
      integer :: i_Inu=1
    contains
      procedure :: pre_init => m3d_pre_init
      procedure :: post_pre_init => m3d_post_pre_init
      procedure :: init => m3d_init 
      procedure :: pre_update => m3d_pre_update
      procedure :: post_update => m3d_post_update
      procedure :: finish
      procedure :: is_finished
      procedure :: nbor_relations => m3d_nbor_relations
      procedure :: pre_init_and_link
      procedure :: set_nv_and_nw
      procedure :: set_initial_time
      procedure :: set_maxiter
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
      procedure :: push_nbors_for_update
      procedure :: push
      procedure :: allocate_mem
      procedure :: deallocate_mem
      procedure :: check_mem
      procedure :: check_mem_NaN
      procedure :: check_access
      procedure :: m3d_patch_info
      procedure :: m3d_patch_info2
      procedure :: print_to_stdout
      procedure :: print_post_update_info
      procedure :: check_in_bounds
      procedure :: print_execution_history
      procedure :: init_nbors_hash
      procedure :: check_nbor
      procedure :: hash_key
      procedure :: output
      procedure :: save
    end type m3d_patch_t
  
    !-----------------------------------------------------------------------------
    integer, save:: verbose=0
    integer, public :: maxiter, ilambd
    integer, public :: cswitch=1
    integer, public, save :: recl=0, save_id=-1
    real(8), public, save, dimension(:), allocatable:: wmu, long_wmu
    real(8), public, save, dimension(:,:), allocatable:: vec, long_vec
    real(8), public:: convlim=1d-3
    logical, public, save :: converged=.false.
    logical, public, save :: check_NaN=.false., save_back_cont=.true.
    logical, public, save :: detailed_prints=.false., comoving_frame=.false.
    logical, public, save :: print_start=.false., print_push=.false., print_send=.false.
    character(len=128), save :: save_patch_kind=''
    character(len=64), public, save:: custom_mu
    real(8), parameter :: t_push=0.1d0
    real(4), public :: vel_max=0.0
    type(mpi_counter_t), public :: atom_update_counter, global_maxerr, global_pos
      
    public :: hash_string, count_words

  contains
  
  !===============================================================================
  SUBROUTINE m3d_pre_init (self)
    class(m3d_patch_t), target:: self
    !-----------------------------------------------------------------------------
    logical, save :: first_time=.true.
    logical, save :: on=.true.
    integer :: iostat, m(6)
    character(len=256) :: command
    integer, save :: fcheck=1, linecheck=3, pcheck(3)=[1,1,1], lvlcheck=1
    integer, save :: n_nu=1, itimer=0, ntheta=3, nphi=4, total_n(3)
    character(len=32), save :: quad_scheme='set_a2'
    character(len=32), save :: long_scheme='lobatto'
    real(8), dimension(:), allocatable:: cust_mu
    namelist /m3d_params/ on, verbose, fcheck, linecheck, pcheck, quad_scheme, &
    long_scheme, n_nu, maxiter, ilambd, lvlcheck, comoving_frame, save_id, &
    ntheta, nphi, check_NaN, custom_mu, lse_no_pivot, save_patch_kind, cswitch, &
    save_back_cont, detailed_prints, print_finished, print_start, print_push, print_send
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_init', itimer=itimer)
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
      !-------------------------------------------------------------------------
      ! Assign angles and integration weights
      !-------------------------------------------------------------------------
      call set_carlson_quadratures (quad_scheme, vec, wmu, up_and_down=.true.)

      if (long_scheme == 'lobatto') then
        call set_lobatto_quadratures (ntheta, nphi, long_vec, long_wmu)
      elseif (long_scheme == 'custom') then
        call count_words(custom_mu, ntheta)
        allocate(cust_mu(ntheta))
        read(custom_mu, *) cust_mu 
        call set_custom_quadratures (cust_mu, nphi, long_vec, long_wmu)
      else
        call set_carlson_quadratures (long_scheme, long_vec, long_wmu, up_and_down=.false.)
      endif

      total_n = nint(self%patch%box / self%patch%ds)
  
      where (pcheck < 0)
        pcheck = total_n + pcheck + 1
      end where
    end if
  
    self%verbose = verbose

    self%n_omega = size(vec, 2)
    self%n_out = size(long_vec, 2)
    self%n_mu = self%n_omega + self%n_out
    self%n_nu = n_nu
    self%out_ray = self%i_nu > self%n_nu .or. self%i_mu > self%n_omega
  
    if (self%pre_initialized) print *, 'pre_init issued somewhere'
    call io%lock%unset ('m3d_patch_t%pre_init')
  
    call self%patch_t%pre_init()
    if (.not.associated(self%lock)) then
      allocate (self%lock)
      call self%lock%init('m3d')
    end if
  
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
    self%place      = self%patch%place
    self%mesh       => self%patch%mesh
    self%task_list  => self%patch%task_list

    self%zs = nint(self%origin(3) / self%ds(3)) + 1 
    self%ze = self%zs + self%n(3) - 1

    self%pcheck = pcheck
    self%fcheck = fcheck
    self%linecheck = linecheck
    self%lvlcheck = lvlcheck

    self%save_output = self%kind == 'atom_patch' .or. self%id == save_id &
    .or. string_in_text(save_patch_kind, self%kind)

    ! ----------------------------------------------------------------------------
    ! Do not check density or validate with it. This is rather important to be set
    ! ----------------------------------------------------------------------------
    call self%set (bits%no_density)
    call self%set (bits%no_filled)

    call self%set_maxiter()
    call trace%end (itimer)
  END SUBROUTINE m3d_pre_init

  !===============================================================================
  SUBROUTINE m3d_post_pre_init (self)
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

    self%pcheck(3) = self%pcheck(3) - self%zs + 1
    self%pcheck = self%pcheck + self%ng
    self%gn = self%n + 2*self%ng
    self%print_update_info = self%check_in_bounds()
    if (self%ng(3)==0) then
      call self%set(bits%no_download)
    endif
    call self%set_initial_time()
    call self%append_to_task_list()
    call trace%end(itimer)
  END SUBROUTINE m3d_post_pre_init

  !===============================================================================
  recursive SUBROUTINE pre_init_and_link (self, m3d_parent, i_mu, i_nu, sp)
    class(m3d_patch_t), target :: self, m3d_parent
    class(atomic_spectrum_type), pointer, optional :: sp
    integer, optional :: i_mu, i_nu
    !-----------------------------------------------------------------------------
    if (present(i_mu)) self%i_mu = i_mu
    if (present(i_nu)) self%i_nu = i_nu
    if (present(sp)) self%sp => sp
    self%patch => m3d_parent%patch
    self%m3d_parent => m3d_parent
    call self%pre_init()
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
  SUBROUTINE m3d_init (self)
    class(m3d_patch_t) :: self
    integer :: li(3), ui(3), status
    integer, save:: itimer=0
    ! ----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_init', itimer=itimer)
    ! internal is set in task%pre_init, which is called by patch%init
    status = self%status
    call self%patch_t%init()
    call data_io%init (self)
    self%status = status
    call self%init_level()

    self%mem = tiny(1.)
    ! ----------------------------------------------------------------------------
    ! Set unsigned and pervolume. Opacity, source function,
    ! intensity and populations are all positive definite
    ! Although quantities are generally per volume, but density is not stored
    ! locally, thus this flag must not be set.
    ! ----------------------------------------------------------------------------
    self%pervolume = .false.
    self%unsigned  = .true.

    li = self%li
    ui = self%ui

    self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
    if (self%save_output) self%out_mem => self%authority_mem(:,:,:,:,self%it,:)
  
    self%initialized=.true.
    call trace%end(itimer)
  END SUBROUTINE m3d_init

  !===============================================================================
  SUBROUTINE set_initial_time (self)
    class(m3d_patch_t) :: self
    real(8) :: tnu
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_initial_time')
    !-----------------------------------------------------------------------------
    ! Synchronize with MHD patch
    ! ----------------------------------------------------------------------------
    self%nt = 1
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))

    self%time     = self%patch%time
    self%dt_fixed = self%patch%dtime
    self%dtime    = self%patch%dtime
    self%dt       = self%patch%dtime
    self%out_next = self%patch%out_next

    select case (self%kind)
    case ('composition')
      self%t_offset = 0.1d0
    case ('atom_patch')
      self%t_offset = 0.2d0
    case ('back_line')
      self%t_offset = 0.3d0
    case ('out_rad')
      self%t_offset = 0.9d0
    case default
      self%t_offset = 0.5d0
      if (.not. self%out_ray) then
        tnu = t_push * (self%i_nu-1) / self%n_nu
        self%t_offset = self%t_offset + tnu
        if (self%kind == 'back_cont' .and. save_back_cont) then
          self%t_offset = 0.3d0
        endif
      endif
    end select

    if (self%out_ray) then
      self%time = self%time + ilambd + maxiter
      if (self%patch%boundaries%is_set(bits%zl)) then
        self%time = self%time - t_push / 2
      endif
    endif

    self%time = self%time + self%t_offset * self%dtime

    self%t  = self%time            ! in case already set
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
      self%maxiter = min(0, self%m3d_parent%maxiter-1)
    case ('atom_patch')
      self%maxiter = maxiter
    case ('line_opac', 'atom_cont')
      self%maxiter = self%m3d_parent%maxiter - 1
    case ('atom_rate')
      self%maxiter = self%m3d_parent%maxiter + ilambd - 1
      self%iter = ilambd
    case ('back_line')
      self%maxiter = min(0, maxiter - 1)
    case ('back_cont')
      if (save_back_cont) then
        self%maxiter = min(ilambd, maxiter - 1)
      else
        self%maxiter = ilambd + maxiter - 1
      endif
    case ('mean_rad')
      self%maxiter = ilambd - 1
      if (comoving_frame) then
        self%maxiter = self%maxiter + maxiter
      endif
    case ('angle_rad', 'net_opac', 'proj_patch')
      self%maxiter = ilambd + maxiter - 1
    end select

    if (self%out_ray) self%maxiter = 0
    call trace%end()
  END SUBROUTINE set_maxiter

  !===============================================================================
  SUBROUTINE m3d_pre_update (self)
    class(m3d_patch_t), target :: self
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_update', itimer=itimer)
    if (print_start) then
      write(stdout, *) 'Update on', self%id, 'has begun', self%finished, self%time
    endif
    if (.not. self%is_finished()) then
      if (converged .and. .not. self%out_ray) then
        call self%finish()
        call self%free_unused_memory()
      else
        call self%check_mem()
      endif
      if (.not. self%is_finished()) then
        call self%print_execution_history('PRE')
      endif
    endif
    call trace%end(itimer)
  END SUBROUTINE m3d_pre_update
  
  !===============================================================================
  SUBROUTINE m3d_post_update (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_post_update')
    if (self%kind == 'atom_patch' .and. converged) then
      call self%print_execution_history('CON')
    else
      call self%print_execution_history('END')
    endif
    call io%lock%set
    if (self%verbose > 1) then
      if (self%print_update_info) then
        call self%print_post_update_info()
      elseif (self%verbose > 2) then
        call self%print_to_stdout()
      endif
    endif
    if (self%iter == self%maxiter .or. self%nbors_finished()) then
      call self%finish()
    endif
    call self%free_unused_memory()
    call self%push_nbors_for_update()
    call io%lock%unset
    if (check_NaN) then
      call self%check_mem_NaN()
    endif
    self%iter = self%iter + 1
    call trace%end()
  END SUBROUTINE m3d_post_update

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
    call self%set(bits%no_download)
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
                if (nbnbor%task%id == self%id) then
                  if (.not. self%is_finished()) then
                    if (self%time + self%t_push + self%dtime < nbtask%time) then
                      mem_not_needed=.false.
                    endif
                  endif
                elseif (nbnbor%needs_me) then
                  if (nbnbtask%time < nbtask%time) then
                    mem_not_needed=.false.
                  endif
                endif
              endif
              nbnbor => nbnbor%next
            end do
            if (mem_not_needed .and. nbtask%initialized) then
              call nbtask%mem_lock(1)%set ('free_unused_memory')
              call nbtask%deallocate_mem(self, print=.true.)
              call nbtask%mem_lock(1)%unset ('free_unused_memory')
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

    needed=.true. ; needs_me=.true. ; download=.true.

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
            if(self%i_mu > 0 .and. another%i_mu > 0) then
              if (another%i_mu /= self%i_mu) then
                needed=.false. ; needs_me=.false. ; download=.false.
              endif
            endif

            if(self%i_nu > 0 .and. another%i_nu > 0) then
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
        call link%init()
        link%task => self
        self%link => link
        !call link%append_to(self%task_list)
        call self%task_list%append_link(self%link)
        self%key = self%hash_key()
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
      'm3d_patch_t%task_info: kind, id, rank, level, mu, nu, BV, pos', &
      self%kind, self%id, self%rank, self%level, self%i_mu, self%i_nu, &
      self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
  END SUBROUTINE task_info
  
  !===============================================================================
  SUBROUTINE init_level (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    self%level = 1
  END SUBROUTINE init_level

  !===============================================================================
  SUBROUTINE rotate (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    call self%task_t%rotate
  END SUBROUTINE rotate
  
  !===============================================================================
  SUBROUTINE periodic_grid (self, only)
    class(m3d_patch_t):: self
    integer, optional:: only
    integer :: iv, it, iw, i,j,k,t
    integer, save :: itimer=0
    logical :: only_z
    !----------------------------------------------------------------------------
    if (.not. self%is_periodic()) return
    if (all(self%ng(:2) == 0) .or. .not. associated(self%mem)) return
    call trace%begin ('m3d_patch_t%periodic_grid', itimer=itimer)
    call self%check_access(self)

    only_z = present(only)
    if (only_z) k = only

    if (.not. associated(self%periodic_idx)) then
      call self%init_periodic_idx()
    endif

    do iw=1, self%nw
    do it=1, self%nt
    do iv=1, self%nv
      if (only_z) then
        do t=1, size(self%periodic_idx, 2)
          i = self%periodic_idx(1,t)
          j = self%periodic_idx(2,t)
          self%mem(i,j,k,iv,it,iw) = self%mem(self%periodic_i(i),self%periodic_j(j),k,iv,it,iw)
        end do
      else
        do k=1, self%gn(3)
          do t=1, size(self%periodic_idx, 2)
            i = self%periodic_idx(1,t)
            j = self%periodic_idx(2,t)
            self%mem(i,j,k,iv,it,iw) = self%mem(self%periodic_i(i),self%periodic_j(j),k,iv,it,iw)
          end do
        end do
      endif
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
    real(4), pointer, optional:: target_mem(:,:,:,:)
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
    class (m3d_patch_t):: self
    class (task_t):: target
    !.............................................................................
    real:: grace
    real(8):: nbtime, nbdtime, tgtime
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%is_ahead_of')
    grace = self%grace + shared%extra_grace
    if (omp_lock%tasks) then
                                              call self%lock%set ('is_ahead_of')
      !---------------------------------------------------------------------------
      ! dt(self%dt) is always consistent with t(:), while dtime changes in
      ! courant_condition(), before other attributes change in rotate()
      !---------------------------------------------------------------------------
      nbtime = self%time
      tgtime = target%time
      if (self%is_set(bits%frozen)) then
        is_ahead_of = .true.
      else
        is_ahead_of = nbtime > tgtime
      end if
                                              call self%lock%unset ('is_ahead_of')
    else
      !$omp atomic read
      nbtime = self%time
      !$omp atomic read
      nbdtime = self%dtime
      !$omp atomic read
      tgtime = target%time
      ! -- if a task is frozen (in time) assume it is forever ahead of other tasks
      if (self%is_set(bits%frozen)) then
        is_ahead_of = .true.
      ! -- a negative nbor time signals an invalid (e.g. initial RT) state
    !  else if (nbtime < 0.0) then
    !    is_ahead_of = .false.
      ! -- do not use a grace interval for the first few steps
      else if (self%istep < 3) then
        is_ahead_of = nbtime > tgtime
      ! -- use a grace interval that is fraction of the nbor time step
      else
        is_ahead_of = nbtime + nbdtime*grace > tgtime 
      end if
    end if
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
      if (self%i_mu > 0 .and. task%i_mu > 0) then
        if (self%i_mu /= task%i_mu) then
          d = huge(1.0)
          return
        endif
      endif
      if (self%i_nu > 0 .and. task%i_nu > 0) then
        if (self%i_nu /= task%i_nu) then
          d = huge(1.0)
          return
        endif
      endif
    end select
    d = self%task_t%distance_to(task)
  END FUNCTION distance_to

  !===============================================================================
  SUBROUTINE allocate_mem (self, allocator, m_in)
    class(m3d_patch_t):: self
    class(m3d_patch_t), intent(in), target :: allocator
    integer, intent(in), optional :: m_in(6)
    integer :: m(6), li(3), ui(3)
    integer, save :: itimer
    !-----------------------------------------------------------------------------
    if (.not. associated(self%mem)) then
      call trace%begin ('m3d_halo_t%allocate_mem', itimer=itimer)
      if (present(m_in)) then
        m = m_in
      else
        m(:3) = self%gn
        m( 4) = self%nv
        m( 5) = self%nt
        m( 6) = self%nw
      endif
      call os%allocate(self%mem, m, label=trim(self%kind)//'%mem')
      if (.not. self%compressed) then
        li = self%li
        ui = self%ui
        self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
      endif
      self%last_allocator => allocator
      self%allocation_time = allocator%time
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
        call os%deallocate(self%mem, label=trim(self%kind)//'%mem')
        self%last_deallocator => deallocator
        self%deallocation_time = deallocator%time
        if (present(print)) then
          if (print .and. (self%print_update_info .or. self%verbose > 2) .and. detailed_prints) then
            write(stdout, *) '  deallocated', self%m3d_patch_info2()
          endif
        endif
      endif
    endif
  END SUBROUTINE deallocate_mem

  !===============================================================================
  SUBROUTINE check_mem (self)
    class(m3d_patch_t):: self
    integer :: orig_place
    !-----------------------------------------------------------------------------
    if (.not. associated(self%mem)) then
      if (.not. self%initialized) then
        if (self%is_set(bits%virtual)) then
          self%kind = trim(self%kind)//'V'
          self%static_mem = .true.
        endif
        orig_place = self%place
        call self%init()
        if (orig_place == -1) self%place = -1
      else
        call self%allocate_mem(self)
      endif
      if (self%place /= omp%place .and. self%i_mu > 0) then
        self%place = omp%place ! memory affinity
      else
      endif
    endif
  END SUBROUTINE check_mem

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
  SUBROUTINE push_nbors_for_update (self)
    class(m3d_patch_t) :: self
    class(nbor_t), pointer:: nbor, nbnbor
    class(task_t), pointer:: nbtask, nbnbtask
    integer :: n_ready, n_needed, n_done
    integer, save :: itimer=0
    logical :: ok
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%push_nbors_for_update', itimer=itimer)
    select case (self%kind)
    case ('composition', 'atom_patch', 'atom_rate')
      ok = .false.
    case ('back_cont', 'atom_cont')
      ok = self%out_ray
    case default
      ok = .true.
    end select
    if (ok) then
      nbor => self%link%nbor
      do while (associated (nbor))
        nbtask => nbor%task
        n_needed = 0
        n_ready = 0
        n_done = 0
        select type (nbtask)
        class is (m3d_patch_t)
        if (nbor%needs_me .and. nbtask%iter==self%iter .and. nbtask%t_push==0) then
            nbnbor => nbtask%link%nbor
            do while (associated(nbnbor))
              if (nbnbor%needed .and. .not. nbnbor%needs_me) then
                nbnbtask => nbnbor%task
                n_needed = n_needed + 1
                if (nbnbtask%id == self%id) then
                  n_done = n_done + 1
                elseif (nbnbtask%time > nbtask%time) then
                  n_done = n_done + 1
                elseif (nbnbtask%is_set(bits%ready)) then
                  n_ready = n_ready + 1
                endif
              endif
              nbnbor => nbnbor%next
            end do
            if (n_done == n_needed) then
              call self%push(nbtask)
            elseif ((n_done + n_ready) == n_needed) then
              nbnbor => nbtask%link%nbor
              do while (associated(nbnbor))
                if (nbnbor%needed) then
                  nbnbtask => nbnbor%task
                  if (nbnbtask%is_set(bits%ready) .and. nbnbtask%id /= self%id) then
                    select type (nbnbtask)
                    class is (m3d_patch_t)
                      if (nbnbtask%iter == self%iter .and. nbnbtask%is_clear(bits%busy) .and. nbnbtask%t_push == 0) then
                        call self%push(nbnbtask)
                        call self%task_list%qlock%set ('remove')
                        call self%task_list%remove_from_queue(nbnbor%link)
                        call self%task_list%qlock%unset ('remove')
                        call self%task_list%queue_by_time(nbnbor%link, nbnbtask%istep)
                      endif
                    end select
                  endif
                endif
                nbnbor => nbnbor%next
              end do
            endif
          end if
        end select
        nbor => nbor%next
      end do
    endif

    if (self%t_push > 0) then
      self%time = self%time + self%t_push
      self%t_push = 0
    endif
    call trace%end(itimer)
  END SUBROUTINE push_nbors_for_update

  !===============================================================================
  SUBROUTINE push (self, other)
    class(m3d_patch_t):: self, other
    !-----------------------------------------------------------------------------
    other%t_push = other%t_push + t_push
    other%time = other%time - t_push
    if (other%i_mu > 0) then
      other%place = self%place
    endif
    if ((self%print_update_info .or. self%verbose > 2) .and. print_push) then
      write(stdout, *) '  pushed', other%m3d_patch_info2()
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
      !call other%mem_lock(1)%set('check_access')
      if (.not. associated(other%mem)) then
        print *, 'memory of ', other%m3d_patch_info()
        print *, 'cannot be accessed from ', self%m3d_patch_info()
        print *, 'static mem =', other%static_mem
        print *, 'it was last deallocated by', other%last_deallocator%id, other%deallocation_time
        print *, 'it was last   allocated by', other%last_allocator%id, other%allocation_time
      endif
      !call other%mem_lock(1)%unset('check_access')
    end select
  END SUBROUTINE check_access

  !===============================================================================
  SUBROUTINE print_to_stdout (self, arr, text)
    class(m3d_patch_t) :: self
    real(4), dimension(:), optional :: arr
    character(len=*), optional :: text
    !-----------------------------------------------------------------------------
    if (present(arr)) then
      write(stdout,'(A, 6(ES15.7))') self%m3d_patch_info(), arr
    elseif (present(text)) then
      write(stdout,'(A, A, A)') self%m3d_patch_info(), ' ', trim(text)
    else
      write(stdout,'(A)') self%m3d_patch_info()
    endif
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
  SUBROUTINE count_words (string, nwords)
    character(len=*), intent(in) :: string
    integer, intent(out) :: nwords
    integer :: i, n
    logical :: new_word=.true.
    !-----------------------------------------------------------------------------
    n = len(string)
    nwords = 0
    do i=1, n
      if (string(i:i) == ' ') then
        new_word = .true.
      elseif(new_word) then
        nwords = nwords + 1
        new_word = .false.
      endif
    end do
  END SUBROUTINE count_words

  !===============================================================================
  SUBROUTINE print_execution_history (self, status)
    class(m3d_patch_t) :: self
    character(len=3) :: status
    !-----------------------------------------------------------------------------
    write(io_unit%dispatcher, '(F17.8, I3, 1x, A3, 1x, A)') wallclock(), mpi%rank, status, self%m3d_patch_info()
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
  FUNCTION hash_string(text) result(hashed) 
    character(*), intent(in) :: text 
    integer :: hashed ! 32-bit integers assumed 
    integer, parameter :: magic_numb = Z'5D7A9F43'  ! hexadecimal
    integer :: i, j 
    !-----------------------------------------------------------------------------
    hashed = 0
    do i = 1, len_trim(text) 
      j = mod(i-1, 4) * 8 
      hashed = IEOR( hashed, ISHFT( ICHAR( text(i:i) ), j ) ) 
    end do
    hashed = abs( IEOR( hashed, magic_numb ) )
  END FUNCTION hash_string 

  !===============================================================================
  !> Fast method to generate nbor lists in cases with very long task lists
  !===============================================================================
  SUBROUTINE init_nbors_hash (self, nbors)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbors
    !.............................................................................
    class(task_t), pointer:: task
    class(nbor_t), pointer:: nbor
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
  SUBROUTINE output (self)
    class(m3d_patch_t):: self
    integer:: snap, rec, unit, id, level, ipos(3), n(3), m(5), nv, nw, rank
    integer:: i_nu, i_mu
    real(8), dimension(3) :: ds, position, size, box
    character(len=128) filename, knd
    class(snapshot_t), pointer:: snapshot
    namelist /m3d_out_nml/ snap, box, ds
    namelist /patch_nml/ id, position, size, level, rank, ipos, i_nu, i_mu, knd, &
    n, nv, nw, m
    !-----------------------------------------------------------------------------
    call io%lock%set ('m3d_patch_t%output')
    snap = self%iout
    snapshot => snapshots%get (snap)
    !-----------------------------------------------------------------------------
    ! Open a direct access file if not already open
    !-----------------------------------------------------------------------------
    if (snapshot%data_unit == 0) then
      snapshot%record = 0
      write (filename,'(a,i5.5,"/snapshot_",i5.5,".dat")') &
        trim(io%run_dir), snap, mpi%rank
      open (newunit=snapshot%data_unit, file=trim(filename), form='unformatted', &
          recl=recl,  access='direct', status='unknown')
    end if
    unit = snapshot%data_unit
    snapshot%record = snapshot%record + 1
    self%record = snapshot%record
    !-----------------------------------------------------------------------------
    ! Write out
    !-----------------------------------------------------------------------------
    write (unit, rec=snapshot%record) self%authority_mem
    flush (unit)
    !-----------------------------------------------------------------------
    ! Output namelist file with metadata
    !-----------------------------------------------------------------------
    unit = snapshot%meta_unit
    if (unit == 0) then
      box = self%box
      ds  = self%ds
      write(filename,'(a,i5.5,"/snapshot.nml")') trim(io%outputname), snap
      open (io_unit%nml1, file=trim(filename), form='formatted', &
            delim='apostrophe', status='unknown')
      write (io_unit%nml1, m3d_out_nml)
      close (io_unit%nml1)
      !-------------------------------------------------------------------------
      ! Open a new "data/run/rank_rrrrr_patches.nml"
      !-------------------------------------------------------------------------
      write (filename,'(a,i5.5,"/rank_",i5.5,"_patches.nml")') &
        trim(io%outputname), snap, mpi%rank
      open (newunit=unit, file=trim(filename), form='formatted', &
            delim='apostrophe', status='unknown')
      snapshot%meta_unit = unit
    endif
    id       = self%id
    position = self%position
    size     = self%size
    level    = self%level
    rank     = self%rank 
    ipos     = self%ipos 
    i_nu     = self%i_nu
    i_mu     = self%i_mu
    knd      = self%kind 
    n        = self%n    
    nv       = self%nv   
    nw       = self%nw   
    m        = [n, nv, nw]
    write (unit, patch_nml)
    flush (unit)
    call io%lock%unset ('m3d_patch_t%output')
  END SUBROUTINE output

  !===============================================================================
  SUBROUTINE save (self, mem, label)
    class(m3d_patch_t):: self
    real(4), dimension(:,:,:,:,:) :: mem 
    character(len=*) :: label
    character(len=128) :: meta_file, folder
    character(len=32) :: filename
    integer :: lunit, iw, iv
    ! ----------------------------------------------------------------------------
    write(folder, '(A, A)') trim(io%run_dir), 'save/'
    write(filename, '(A, A, 3(I3.3), A)') trim(label), '_', self%key(1:3), '.bin'
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
    write(lunit, '(A32, I4, 8(I10))') filename, storage_size(mem), shape(mem)
    close(lunit)
    call io%lock%unset()
    open(newunit=lunit, file=trim(folder)//filename, access='stream')
    write(lunit) mem
    close(lunit)
  END SUBROUTINE save

  !===============================================================================
  FUNCTION string_in_text (text, str) result(in_text)
    character(len=*) :: str
    character(len=*) :: text
    character(len=128) :: dummy
    integer :: i, j 
    logical :: in_text
    ! ----------------------------------------------------------------------------
    in_text = .false.
    dummy = ''
    j = 0
    do i=1, len(text)
      if (text(i:i) == ' ') then
        if (dummy(1:j) == str) exit
        j = 0
      else
        j = j + 1
        dummy(j:j) = text(i:i)
      endif
    end do
    if (dummy(1:j) == str) in_text = .true.
  END FUNCTION

END MODULE m3d_patch_mod