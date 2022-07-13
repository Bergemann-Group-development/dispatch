module m3d_patch_mod

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
    !USE index_mod
    USE mesh_mod
    USE io_mod
    !USE aux_mod
    USE timer_mod
    !USE scaling_mod
    USE trace_mod
    USE mpi_mod
    USE omp_mod
    USE omp_lock_mod
    USE shared_mod
    !USE omp_timer_mod
    USE math_mod
    USE bits_mod
    USE mconst
    !USE units_mod
    !USE scratch_mod
    USE atom_mod, only: atomic_spectrum_type
    USE ray_mod, only: set_carlson_quadratures, set_lobatto_quadratures
    implicit none
    private
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(gpatch_t) :: m3d_patch_t
      class(gpatch_t),      pointer :: patch => null()         ! pointer to mhd patch
      class(m3d_patch_t),   pointer :: m3d_parent => null()    ! pointer to parent patch
      class(atomic_spectrum_type), pointer :: sp => null()
      class(m3d_patch_t),   pointer :: last_allocator, last_deallocator
      integer                  :: master_id
      logical                  :: on=.false.       ! allows local turn-off
      logical                  :: in_task_list = .false.
      logical                  :: initialized  = .false.
      real(4), dimension(:,:,:,:,:,:), pointer :: authority_mem => null()
      real(8) :: t_offset=0, t_push=0
      real(8) :: allocation_time, deallocation_time
      integer :: iter=0, maxiter, master_maxiter
      integer :: n_mu, n_nu, n_omega, n_out, zs, ze
      integer :: fcheck, pcheck(3), linecheck, lvlcheck
      logical :: out_ray=.false., static_mem=.false.
      logical :: compressed=.false.
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
      procedure :: finish => m3d_finish
      procedure :: nbor_relations => m3d_nbor_relations
      procedure :: pre_init_and_link
      procedure :: bndries
      procedure :: set_nv_and_nw
      procedure :: set_initial_time
      procedure :: set_maxiter
      procedure :: m3d_name
      procedure :: init_level
      procedure :: task_info
      procedure :: append_to_task_list
      procedure :: extrapolate_guards
      procedure :: is_ahead_of
      procedure :: periodic_grid => m3d_periodic_grid
      procedure :: index_only => m3d_index_only
      procedure :: distance_to
      procedure :: free_unused_memory
      procedure :: skip_iterations
      procedure :: push_nbors_for_update
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
    end type m3d_patch_t
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(halo_t):: halo_m3d_t
    contains
      procedure:: info
      procedure:: download
    end type
  
    !-----------------------------------------------------------------------------
    integer, save:: verbose=0
    integer, public :: maxiter, ilambd
    real(8), public, save, dimension(:), allocatable:: wmu, long_wmu
    real(8), public, save, dimension(:,:), allocatable:: vec, long_vec
    logical, public, save :: converged=.false., dynamic_memory=.true.
    logical, public, save :: check_NaN=.false., save_back_opac=.true.
    logical, public, save :: detailed_prints=.false., print_finished=.false.
    logical, public, save :: print_start=.false.
    character(4), public, save :: conditioning='full'
    character(len=64), public, save:: out_folder = './out/'
      
  contains
  
  !===============================================================================
  SUBROUTINE m3d_pre_init (self)
    class(m3d_patch_t), target:: self
    !-----------------------------------------------------------------------------
    logical, save :: first_time=.true.
    logical, save :: on=.true.
    integer :: iostat
    character(len=256) :: command
    integer, save :: fcheck=1, linecheck=3, pcheck(3)=[1,1,1], lvlcheck=1
    integer, save :: n_nu=1, itimer=0, ntheta=3, nphi=4, total_n(3)
    character(len=32), save :: quad_scheme='set_a2'
    character(len=32), save :: long_scheme='disk_center'
    namelist /m3d_params/ on, verbose, fcheck, linecheck, pcheck, quad_scheme, &
    long_scheme, n_nu, maxiter, ilambd, dynamic_memory, lvlcheck, &
    ntheta, nphi, check_NaN, conditioning, out_folder, save_back_opac, &
    detailed_prints, print_finished, print_start
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
      call os%mkdir(trim(out_folder))
      !write (command, *) 'mkdir -p ', trim(out_folder)
      !call execute_command_line (trim(command), exitstat=iostat)
      write (command, *) 'rm ', trim(out_folder), '*'
      call execute_command_line (trim(command), exitstat=iostat)

      !-------------------------------------------------------------------------
      ! Assign angles and integration weights
      !-------------------------------------------------------------------------
      call set_carlson_quadratures (quad_scheme, vec, wmu, up_and_down=.true.)

      if (long_scheme == 'lobatto') then
        call set_lobatto_quadratures (ntheta, nphi, long_vec, long_wmu)
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
  
    call self%gpatch_t%pre_init()
    if (.not.associated(self%lock)) then
      allocate (self%lock)
      call self%lock%init('m3d')
    end if
    call self%set_initial_time()
    call self%set_maxiter()
  
    self%on = on
    self%in_task_list = .false.
  
    ! ----------------------------------------------------------------------------
    ! relevant synchronization items for cartesian setup
    ! ----------------------------------------------------------------------------
    self%master_id  = self%patch%id
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
    self%status     = self%patch%status
    self%periodic   = self%patch%periodic
    self%mesh       => self%patch%mesh
    self%task_list  => self%patch%task_list

    select case (self%kind)
    case ('net_opac', 'angle_rad')
     self%ng(3) = 1 ! needed for nbor_relations
    end select

    !-----------------------------------------------------------------------------
    ! Preliminary region information, relative to the task position.  This may
    ! change in patch_t%init, but the settings below allow overlap() calls to work
    !-----------------------------------------------------------------------------
    self%authority%l = - self%size*0.5_8
    self%authority%u = + self%size*0.5_8
    self%interest%l  = self%authority%l - self%ng*self%size/self%n
    self%interest%u  = self%authority%u + self%ng*self%size/self%n

    self%zs = nint(self%origin(3) / self%ds(3)) + 1 
    self%ze = self%zs + self%n(3) - 1

    self%pcheck = pcheck
    self%fcheck = fcheck
    self%linecheck = linecheck
    self%lvlcheck = lvlcheck

    ! ----------------------------------------------------------------------------
    ! Do not check density or validate with it. This is rather important to be set
    ! ----------------------------------------------------------------------------
    call self%set (bits%no_density)
    call self%set (bits%internal)
    call self%set (bits%external)
    call self%set (bits%static)

    allocate (self%boundaries)
    call self%boundaries%pre_init (self)

    call trace%end (itimer)
  END SUBROUTINE m3d_pre_init

  !===============================================================================
  SUBROUTINE m3d_post_pre_init (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    self%pcheck(3) = self%pcheck(3) - self%zs + 1
    self%pcheck = self%pcheck + self%ng
    self%gn = self%n + 2*self%ng
    self%print_update_info = self%check_in_bounds()
    if (all(self%ng==0)) then
      self%dnload_time = io%end_time
    endif
    self%place = mod((self%i_nu-1)*self%n_mu + self%i_mu, omp%nthreads)
  END SUBROUTINE m3d_post_pre_init

  !===============================================================================
  SUBROUTINE pre_init_and_link (self, m3d_parent, i_mu, i_nu, sp)
    class(m3d_patch_t), target :: self, m3d_parent
    class(atomic_spectrum_type), pointer, optional :: sp
    integer, optional :: i_mu, i_nu
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_init_and_link', itimer=itimer)
    if (present(i_mu)) self%i_mu = i_mu
    if (present(i_nu)) self%i_nu = i_nu
    if (present(sp)) self%sp => sp

    self%patch => m3d_parent%patch
    self%m3d_parent => m3d_parent
    call self%pre_init()
    call trace%end(itimer)
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
    
    ! ----------------------------------------------------------------------------
    ! Allocate and set unsigned and pervolume. Opacity, source function,
    ! intensity and populations are all positive definite
    ! Although quantities are generally per volume, but density is not stored
    ! locally, thus this flag must not be set.
    ! ----------------------------------------------------------------------------
    allocate (self%unsigned(self%nv), self%pervolume(self%nv))
    self%pervolume = .false.
    self%unsigned  = .true.
    call trace%end()
  END SUBROUTINE set_nv_and_nw
  
  !===============================================================================
  SUBROUTINE m3d_init (self)
    class(m3d_patch_t) :: self
    integer :: li(3), ui(3)
    integer, save:: itimer=0
    ! ----------------------------------------------------------------------------
    if (self%is_set(bits%external)) then
      self%mem => null()
      return
    endif
    call trace%begin ('m3d_patch_t%m3d_init', itimer=itimer)
    call self%gpatch_t%init()
    call self%init_level()

    self%mem = tiny(1.)

    li = self%li
    ui = self%ui

    self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
  
    call self%set (bits%internal)
    self%initialized=.true.
    call trace%end(itimer)
  END SUBROUTINE m3d_init

  !===============================================================================
  SUBROUTINE set_initial_time (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_initial_time')
    !-----------------------------------------------------------------------------
    ! Synchronize with MHD patch
    ! ----------------------------------------------------------------------------
    self%nt = 1
    allocate (self%t(self%nt), self%dt(self%nt), self%iit(self%nt))

    self%dt_fixed = self%patch%dtime
    self%dtime    = self%patch%dtime
    self%dt       = self%patch%dtime
    self%out_next = self%patch%out_next

    select case (self%kind)
    case ('composition')
      self%t_offset = 0d0
    case ('atom_patch')
      self%t_offset = 0d0
    case ('lines_patch')
      self%t_offset = 0.3
    case ('back_opac')
      if (save_back_opac) then
        self%t_offset = 0.3d0
      else
        self%t_offset = 0.5d0
      endif
    !case ('cont_opac', 'line_opac')
    !  self%t_offset = 0.6
    !case ('net_opac')
    !  self%t_offset = 0.4
    !case ('angle_rad')
    !  self%t_offset = 0.3
    !case ('mean_rad', 'out_rad')
    !  self%t_offset = 0.8
    !case ('rate_patch')
    !  self%t_offset = 0.4
    case default
      self%t_offset = 0.5d0
    end select

    if (self%out_ray) then
      self%t_offset = self%t_offset + ilambd + maxiter
    endif
  
    self%time = self%patch%time + self%t_offset * self%dtime
    self%t  = self%time            ! in case already set
    self%iit = 1                   ! initial slots
    self%it = 1                    ! slot for current values
    self%new = 1                   ! slot for new values
    self%iit(self%nt) = self%new   ! slot for new values
    self%time = self%t(self%it)    ! initial time
    call trace%end()
  END SUBROUTINE set_initial_time

  !===============================================================================
  SUBROUTINE set_maxiter (self)
    class (m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%set_maxiter')
    self%master_maxiter = maxiter + ilambd

    select case (self%kind)
    case ('composition', 'lines_patch')
      self%maxiter = 0
    case ('atom_patch')
      self%maxiter = maxiter
    case ('rate_patch', 'line_opac', 'cont_opac')
      self%maxiter = self%m3d_parent%maxiter - 1
    case ('back_line')
      self%maxiter = min(0, maxiter - 1)
    case ('back_opac')
      if (save_back_opac) then
        self%maxiter = min(ilambd, maxiter-1)
      else
        self%maxiter = ilambd + maxiter - 1
      endif
    case ('mean_rad')
      if (conditioning=='self') then
        self%maxiter = ilambd + maxiter - 1
      else
        self%maxiter = ilambd - 1
      endif
    case ('angle_rad', 'net_opac')
      self%maxiter = ilambd + maxiter - 1
    end select

    if (self%out_ray) self%maxiter = ilambd
    call trace%end()
  END SUBROUTINE set_maxiter

  !===============================================================================
  SUBROUTINE m3d_pre_update (self)
    class(m3d_patch_t), target :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_update')
    if (print_start) then
      print *, 'Update on', self%id, 'has begun', self%finished, self%time
    endif
    call self%check_mem()
    if (converged .and. .not. self%out_ray) then
      call self%finish()
    endif
    call trace%end()
  END SUBROUTINE m3d_pre_update
  
  !===============================================================================
  SUBROUTINE m3d_post_update (self)
    class(m3d_patch_t) :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_post_update')
    call io%lock%set
    if (self%verbose > 1) then
      if (self%print_update_info) then
        call self%print_post_update_info()
      elseif (self%verbose > 4) then
        call self%print_to_stdout()
      endif
    endif
    if (self%iter == self%maxiter) then
      call self%finish()
    endif
    if (dynamic_memory) then
      call self%free_unused_memory()
      call self%push_nbors_for_update()
    endif
    if (check_NaN) then
      call self%check_mem_NaN()
    endif
    self%iter = self%iter + 1
    call io%lock%unset
    call trace%end()
  END SUBROUTINE m3d_post_update

  !===============================================================================
  FUNCTION m3d_patch_info (self) result(text)
    class(m3d_patch_t) :: self
    character(len=42) :: text
    integer :: i
    !-----------------------------------------------------------------------------
    i = omp_get_thread_num()
    write(text, '(I3, 1x, F5.2, A12, I8, 3(I4))') i, self%time, trim(self%kind), &
    self%id, self%ipos(3)+1, self%i_nu, self%i_mu
  END FUNCTION m3d_patch_info

  !===============================================================================
  FUNCTION m3d_patch_info2 (self) result(text)
    class(m3d_patch_t) :: self
    character(len=38) :: text
    !-----------------------------------------------------------------------------
    write(text, '(A12, I8)') trim(self%kind), self%id
  END FUNCTION m3d_patch_info2
  
  !===============================================================================
  RECURSIVE SUBROUTINE m3d_finish (self, rotated)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbor, nbnbor
    class(task_t), pointer:: nbtask, nbnbtask
    logical, optional :: rotated
    logical :: nbors_finished
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_finish')
    self%finished = .true.
    if (print_finished) print *, self%m3d_patch_info(), 'has finished!'
    self%time = io%end_time + self%t_offset * self%dtime
    call self%set(bits%no_download)
    if (present(rotated)) self%rotated = rotated

    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed .and. nbor%task%is_clear(bits%external)) then
        nbtask => nbor%task
        select type (nbtask)
        class is (m3d_patch_t)
          nbors_finished=.true.
          nbnbor => nbtask%link%nbor
          do while (associated(nbnbor))
            if (nbnbor%needs_me) then
              nbnbtask => nbnbor%task
              select type (nbnbtask)
              class is (m3d_patch_t)
                if (.not. nbnbtask%finished) then
                  nbors_finished=.false.
                endif
              end select
            endif
            nbnbor => nbnbor%next
          end do
          if (nbors_finished) then
            if (nbtask%finished) then
              if (.not. nbtask%static_mem) then
                call nbtask%mem_lock(1)%set ('m3d_finish')
                call nbtask%deallocate_mem(self, print=.true.)
                call nbtask%mem_lock(1)%unset ('m3d_finish')
              endif
            else
              call nbtask%finish()
            endif
          endif
        end select
      endif
      nbor => nbor%next
    end do
    call trace%end()
  END SUBROUTINE m3d_finish

  !===============================================================================
  SUBROUTINE m3d_nbor_relations(self, another, needed, needs_me, download)
    class(m3d_patch_t) :: self
    class(task_t), pointer :: another         ! another task in the list
    logical :: needed, needs_me, download    ! Is the other task: needed, 
                                             !    does it need me, should I download?
    class(task_t), pointer :: other
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_nbor_relations')
    !call self%task_t%nbor_relations(another, needed, needs_me, download)

    if (.not.associated(another)) then
      write (stderr,*) mpi%rank, 'another not associated'
      call io%abort ('task_t another not associated')
    end if

    needed=.true. ; needs_me=.true. ; download=.true.

    if (another%is_set (bits%frozen)) then
      needs_me = .false.
      needed = .false.
    endif

    if (self%rank == mpi%rank) then
      call self%clear(bits%external)
      if (another%rank /= mpi%rank) then
        call self%clear(bits%internal)
      end if
    else
      call self%clear(bits%internal)
      if (another%rank==mpi%rank) then
        call self%set(bits%virtual)
        call self%clear(bits%external)
      !---------------------------------------------------------------------------
      ! task is virtual, nbor is virtual
      !---------------------------------------------------------------------------
      else
        needed = .false.
        needs_me = .false.
      end if
      !---------------------------------------------------------------------------
      ! virtual tasks don't need guard zones
      !---------------------------------------------------------------------------
      download = .false.
    end if

    if (needed .or. needs_me) then
      select type(another)
      class is(m3d_patch_t)
        if (another%kind == self%kind) then
          if (another%master_id == self%master_id) then
            needed=.false. ; needs_me=.false. ; download=.false.
          elseif((another%i_mu /= self%i_mu) .or. (another%i_nu /= self%i_nu)) then
            needed=.false. ; needs_me=.false. ; download=.false.
          endif
        else
          download=.false.
          if (another%master_id /= self%master_id) then
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
        if (another%id == self%master_id) then
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
    class (link_t), pointer :: link
    class (task_t), pointer :: task
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%add_to_task_list')
    if (associated(self%task_list)) then
      if (self%maxiter >= 0) then
        call self%task_t%init_unique (same=.true.) ! obtain an unique ID
        allocate (halo_m3d_t:: link)
        call link%init()
        link%task => self
        self%link => link
        call self%task_list%append (link)
      endif
    endif
    call trace%end()
  END SUBROUTINE append_to_task_list
  
  !===============================================================================
  SUBROUTINE info (self, label)
    class(halo_m3d_t):: self
    character(*), optional:: label
    class(nbor_t), pointer:: nbor
    real(8):: distance(3)
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%info')
    if (verbose < 5) then
      call self%task%task_info()
      nbor => self%nbor
      do while (associated(nbor))
        distance = nbor%task%distance(self%task)
        select type(nt => nbor%task)
        class is(m3d_patch_t)        
          write(io_unit%mpi,'(3x,a16,a,i6,":",i2.2,i4,2x,3l1,2x,2l1,3x,3f9.3, a10, i4)') &
            nt%kind, ' nbor_info: id, rank, level, needed, needs_me, download, BV, pos =', &
            nt%id, nt%level, nt%rank, nbor%needed, nbor%needs_me, &
            nbor%download, nt%is_set(bits%boundary), &
            nt%is_set(bits%virtual), distance/ &
            (0.5d0*(self%task%size+nt%size))
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
  
  !===============================================================================
  !> For diagnostic output of tasks and nbor lists to data/run/rank_RRRRR.log
  !===============================================================================
  SUBROUTINE task_info (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    write(io_unit%mpi,'(1x,a,2x,a10,i6,4i4,7x,2l1,2x,3f9.3)') &
      'm3d_patch_t%task_info: kind, id, rank, level, mu, nu, BV, pos', &
      self%kind, self%id, self%rank, self%level, self%i_mu, self%i_nu, &
      self%is_set(bits%boundary), self%is_set(bits%virtual), self%position
  END SUBROUTINE task_info
  
  !===============================================================================
  SUBROUTINE init_level (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%init_level')
    self%level = 1
    call trace%end()
  END SUBROUTINE init_level
  
  !===============================================================================
  SUBROUTINE bndries (self)
    class(m3d_patch_t):: self
    real(4), dimension(:,:,:,:,:), pointer :: mem
    !----------------------------------------------------------------------------
    call trace_begin('m3d_patch_t%bndries')
    mem => self%mem(:,:,:,:,self%it,:)
    if (self%boundaries%lower .or. self%boundaries%upper) then
      call self%boundaries%boundary(mem)
    endif
    call trace%end()
  END SUBROUTINE bndries
  
  !===============================================================================
  SUBROUTINE m3d_periodic_grid (self, only)
    class(m3d_patch_t):: self
    integer, optional:: only
    integer, dimension(3) :: li, ui, lb, ub
    real(4), dimension(:,:,:,:,:), pointer :: mem
    integer :: z, dir, off, is, ie, v, w
    !----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_periodic_grid')
    if (any(self%ng(:2) > 0)) then
      li = self%mesh%li; ui = self%mesh%ui
      lb = self%mesh%lb; ub = self%mesh%ub
      mem => self%mem(:,:,:,:,self%it,:)
  
      if (all(self%n(:2) == 1)) then
        do w=1, self%nw
          do v=1, self%nv
            do z=1, self%gn(3)
              mem(:,:,z,v,w) = mem(li(1),li(2),z,v,w)
            end do
          end do 
        end do
      else
        do dir=-1,1,2
          is = dir*max(dir*li(1), dir*ui(1))
          ie = dir*min(dir*li(1), dir*ui(1))
  
          do off=dir, dir*self%ng(1), dir 
            mem(is+off,li(2):ui(2),:,:,:) = mem(ie-dir+off,li(2):ui(2),:,:,:)
          end do
        end do
  
        do dir=-1,1,2
          is = dir*max(dir*li(2), dir*ui(2))
          ie = dir*min(dir*li(2), dir*ui(2))
  
          do off=dir, dir*self%ng(2), dir
            mem(:,is+off,:,:,:) = mem(:,ie-dir+off,:,:,:)
          end do
        end do
      endif
    endif
    call trace%end()
  END SUBROUTINE m3d_periodic_grid
  
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
  !> Comprehensive entry point
  !===============================================================================
  SUBROUTINE download (self, nbors, only, all_tgt, all_src)
    class(halo_m3d_t)          :: self
    class(nbor_t), pointer :: nbors
    integer,      optional :: only
    logical,      optional :: all_tgt, all_src
    !.............................................................................
    class(nbor_t), pointer :: tnbors
    class(nbor_t), pointer :: nbor
    class(task_t), pointer :: target, source
    integer, save          :: itimer=0
    !-----------------------------------------------------------------------------
    ! Decide whether to use a coarse clone or not
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_halo_t%download', itimer=itimer)
                                            call self%lock%set ('m3d_halo_t%download')
    call self%copy_nbor_list (self%nbor, tnbors)
                                          call self%lock%unset ('m3d_halo_t%download')
    target => self%task
    !.............................................................................
    ! important for dynamic initialisation of patches
    select type (target)
    class is (m3d_patch_t)
      if (target%finished) then
        call trace%end (itimer)
        return
      endif
      call target%check_mem()
    end select
    !-----------------------------------------------------------------------------
    ! loop over nbors
    !-----------------------------------------------------------------------------
    nbor => nbors
    do while (associated(nbor))
      if (nbor%download) then
        if (nbor%task%is_clear (bits%remove)) then
          source => nbor%task
          select type (target)
          class is (m3d_patch_t)
            call target%check_access(source)
          end select
          call guards%source_to_target (source, target, n_op=self%n_op, &
            only=only, all_src=all_src, all_tgt=all_tgt)
        end if
      end if
      nbor => nbor%next
    end do
    !-----------------------------------------------------------------------------
    ! Boundary conditions
    !-----------------------------------------------------------------------------
    call self%task%bndries()
    call self%task%periodic_grid()
    call self%remove_nbor_list (tnbors)
    call trace%end (itimer)
  END SUBROUTINE download

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
  !> Return only the index in self corresponding to a point in native space.
  !> "pos" is a relative position in target coordinate space, and since we are
  !> not informed about the target grid spacing it must be adjusted for staggering
  !> before the call.
  !===============================================================================
  FUNCTION m3d_index_only (self, pos, iv, roundup, rounddn, nearest, nowrap, local) &
    RESULT (ii)
  class(m3d_patch_t)      :: self
  real(8), intent(in) :: pos(3)
  integer, optional   :: iv
  logical, optional   :: roundup, rounddn, nearest, nowrap, local
  integer             :: ii(3)
  !.............................................................................
  integer             :: i
  real(8)             :: p(3), h(3), md
  class(mesh_t), pointer :: mesh
  !-----------------------------------------------------------------------------
  if (present(iv)) &
    h = self%index_stagger(iv)
  do i=1,3
    if (size(self%mem, i) == 1) then
      ii(i) = 1
      cycle
    end if
    mesh => self%mesh(i)
    if (.not.present(local)) then
      p(i) = pos(i) - self%position(i)                          ! relative position
    else
      p(i) = pos(i)
    end if
    if (self%periodic(i)) then
      if (.not.present(nowrap)) then
        md = modulo(p(i)+0.5_8*mesh%b, mesh%b)
        if (md == 0) md = min(p(i)+0.5_8*mesh%b, mesh%b)
        p(i) = md-0.5_8*mesh%b   ! periodic wrap
      end if
    end if
    p(i) = p(i)/mesh%d + self%offset(i)                         ! float index
    if (present(iv)) then
      p(i) = p(i) + mesh%h(iv)                                  ! staggering
    end if
    if (present(roundup)) then
      ii(i) = ceiling(p(i))                                     ! round up
    else if (present(rounddn)) then
      ii(i) = floor(p(i))                                       ! round dn
    else if (present(nearest)) then
      ii(i) = nint(p (i))                                       ! round nearest
    else
      ii(i) = floor(p(i) + 0.0001)                              ! floor w margin
    end if
  end do
  END FUNCTION m3d_index_only
  
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
            do while (associated(nbnbor))
              if (nbnbor%needs_me .and. nbnbor%task%id /= self%id) then
                nbnbtask => nbnbor%task
                select type (nbnbtask)
                class is (m3d_patch_t)
                  if (nbnbtask%time < nbtask%time) then
                    mem_not_needed=.false.
                  endif
                end select
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
  SUBROUTINE allocate_mem (self, allocator, m_in)
    class(m3d_patch_t):: self
    class(m3d_patch_t), intent(in), target :: allocator
    integer, intent(in), optional :: m_in(6)
    integer :: m(6), li(3), ui(3)
    !-----------------------------------------------------------------------------
    if (.not. associated(self%mem)) then
      if (present(m_in)) then
        m = m_in
      else
        m(:3) = self%gn
        m( 4) = self%nv
        m( 5) = self%nt
        m( 6) = self%nw
      endif
      call os%allocate(self%mem, m, label=trim(self%kind)//'%mem')
      self%thread = omp%thread ! thread affinity
      self%place  = omp%place  ! memory affinity
      if (.not. self%compressed) then
        li = self%li
        ui = self%ui
        self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
      endif
      self%last_allocator => allocator
      self%allocation_time = allocator%time
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
          if (print .and. (self%print_update_info .or. self%verbose > 4) .and. detailed_prints) then
            write(stdout, *) '  deallocated', self%m3d_patch_info2()
          endif
        endif
      endif
    endif
  END SUBROUTINE deallocate_mem

  !===============================================================================
  SUBROUTINE check_mem (self)
    class(m3d_patch_t):: self
    !-----------------------------------------------------------------------------
    if (.not. self%initialized) then
      call self%init()
    elseif (.not. associated(self%mem)) then
      call self%allocate_mem(self)
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
    real(8) :: t_push = 0.05d0
    integer :: n_ready, n_needed, n_done
    integer, save :: itimer=0
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%push_nbors_for_update', itimer=itimer)
    select case (self%kind)
    case ('composition', 'atom_patch', 'rate_patch', 'lines_patch')
      ! do not push others
    case default
      nbor => self%link%nbor
      do while (associated (nbor))
        nbtask => nbor%task
        n_needed = 0
        n_ready = 0
        n_done = 0
        select type (nbtask)
        class is (m3d_patch_t)
        if (nbor%needs_me .and. nbtask%iter==self%iter) then
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
              nbtask%t_push = nbtask%t_push + t_push
              nbtask%time = nbtask%time - t_push
              if ((self%print_update_info .or. self%verbose > 4) .and. detailed_prints) then
                write(stdout, *) '  pushed', nbtask%m3d_patch_info2()
              endif
            elseif ((n_done + n_ready) == n_needed) then
              nbnbor => nbtask%link%nbor
              do while (associated(nbnbor))
                if (nbnbor%needed) then
                  nbnbtask => nbnbor%task
                  if (nbnbtask%is_set(bits%ready) .and. nbnbtask%id /= self%id) then
                    select type (nbnbtask)
                    class is (m3d_patch_t)
                      if (nbnbtask%iter == self%iter .and. nbnbtask%is_clear(bits%busy) .and. nbnbtask%t_push == 0) then
                        nbnbtask%t_push = nbnbtask%t_push + t_push
                        nbnbtask%time = nbnbtask%time - t_push
                        call self%task_list%qlock%set ('remove')
                        call self%task_list%remove_from_queue(nbnbor%link)
                        call self%task_list%qlock%unset ('remove')
                        call self%task_list%queue_by_time(nbnbor%link)
                        if ((self%print_update_info .or. self%verbose > 4) .and. detailed_prints) then
                          write(stdout, *) '  pushed', nbnbtask%m3d_patch_info2()
                        endif
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
    end select

    if (self%t_push > 0) then
      self%time = self%time + self%t_push
      self%t_push = 0
    endif
    call trace%end(itimer)
  END SUBROUTINE push_nbors_for_update
  
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
      call other%mem_lock(1)%set('check_access')
      if (.not. associated(other%mem)) then
        print *, 'memory of ', other%m3d_patch_info()
        print *, 'cannot be accessed from ', self%m3d_patch_info()
        print *, 'it was last deallocated by', other%last_deallocator%id, other%deallocation_time
        print *, 'it was last   allocated by', other%last_allocator%id, other%allocation_time
        print *, 'compressed =', other%compressed
      endif
      call other%mem_lock(1)%unset('check_access')
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

  end module m3d_patch_mod
  