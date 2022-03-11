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
    !USE mesh_mod
    USE io_mod
    !USE aux_mod
    !USE timer_mod
    !USE scaling_mod
    USE trace_mod
    !USE mpi_mod
    !USE omp_mod
    USE omp_lock_mod
    USE shared_mod
    !USE omp_timer_mod
    !USE math_mod
    USE bits_mod
    !USE units_mod
    !USE scratch_mod
    USE spectrum_mod, only: spectrum_type
    USE ray_mod, only:set_carlson_quadratures
    implicit none
    private
  
    !-----------------------------------------------------------------------------
  
    type, public, extends(gpatch_t) :: m3d_patch_t
      class(gpatch_t),      pointer :: patch => null()         ! pointer to mhd patch
      class(m3d_patch_t),   pointer :: m3d_parent => null()    ! pointer to parent patch
      class(spectrum_type), pointer :: sp => null()
      integer                  :: master_id
      logical                  :: on=.false.       ! allows local turn-off
      logical                  :: in_task_list = .false.
      logical                  :: initialized  = .false.
      real(4), dimension(:,:,:,:,:,:), pointer :: authority_mem => null()
      real(8) :: t_offset=0
      integer :: iter=0, maxiter, master_maxiter
      integer :: n_mu, n_nu, n_omega, n_out
      integer :: fcheck, pcheck(3), linecheck
      logical :: out_ray=.false.
      integer :: i_nu=0, i_mu=0
      integer :: i_ne=1, i_rho=2, i_temp=3, i_nh1=4, i_nh6=9
      integer :: i_chi=1, i_eta=2, i_snu=2
      integer :: i_Gij=1, i_Gji=2
      integer :: i_Inu=1
    contains
      procedure :: pre_init => m3d_pre_init
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
      
  contains
  
  !===============================================================================
  SUBROUTINE m3d_pre_init (self)
    class(m3d_patch_t), target:: self
    !-----------------------------------------------------------------------------
    logical, save :: first_time=.true.
    logical, save :: on=.true.
    integer :: iostat
    integer, save :: fcheck=1, linecheck=3, pcheck(3)=[1,1,1]
    integer, save :: n_nu=1
    character(len=32), save :: quad_scheme='set_a2'
    character(len=32), save :: long_scheme='disk_center'
    namelist /m3d_params/ on, verbose, fcheck, linecheck, pcheck, quad_scheme, &
    long_scheme, n_nu, maxiter, ilambd
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_init')
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

      !-------------------------------------------------------------------------
      ! Assign angles and integration weights
      !-------------------------------------------------------------------------
      call set_carlson_quadratures (quad_scheme, vec, wmu, up_and_down=.true.)
      call set_carlson_quadratures (long_scheme, long_vec, long_wmu, up_and_down=.false.)
    end if
  
    self%pcheck = pcheck
    self%fcheck = fcheck
    self%linecheck = linecheck
    self%verbose = verbose

    self%n_omega = size(vec, 2)
    self%n_out = size(long_vec, 2)
    self%n_mu = self%n_omega + self%n_out
    self%n_nu = n_nu
    self%out_ray = self%i_mu > self%n_omega  
  
    where (self%pcheck < 0)
      self%pcheck = self%patch%n + self%pcheck + 1
    end where
  
    if (self%pre_initialized) print *, 'pre_init issued somewhere'
    call io%lock%unset ('m3d_patch_t%pre_init')
  
    call self%gpatch_t%pre_init()
    if (.not.associated(self%lock)) then
      allocate (self%lock)
      call self%lock%init('m3d')
    end if
    call self%set_initial_time()
  
    self%on = on
    self%in_task_list = .false.
    if (self%on) then 
      call self%task_t%init_unique (same=.true.) ! obtain an unique ID
    end if
  
    ! ----------------------------------------------------------------------------
    ! relevant synchronization items for cartesian setup
    ! ----------------------------------------------------------------------------
    self%n          = self%patch%n
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

    self%ng         = 0

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

    ! ----------------------------------------------------------------------------
    ! Do not check density or validate with it. This is rather important to be set
    ! ----------------------------------------------------------------------------
    call self%set (bits%no_density)

    allocate (self%boundaries)
    call self%boundaries%pre_init (self)
    call trace%end ()
  END SUBROUTINE m3d_pre_init

  !===============================================================================
  SUBROUTINE pre_init_and_link (self, m3d_parent, i_mu, i_nu)
    class(m3d_patch_t), target :: self, m3d_parent
    integer, optional :: i_mu, i_nu
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_init_and_link')
    if (present(i_mu)) self%i_mu = i_mu
    if (present(i_nu)) self%i_nu = i_nu

    self%patch  => m3d_parent%patch
    self%m3d_parent => m3d_parent
    if (self%kind /= 'atom_patch') self%sp => m3d_parent%sp
    call self%pre_init()
    call trace%end()
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
    ! ----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_init')
    call self%gpatch_t%init()
    call self%init_level()
    call self%set_maxiter()

    self%mem = tiny(1.)

    li = self%li
    ui = self%ui

    self%authority_mem => self%mem(li(1):ui(1),li(2):ui(2),li(3):ui(3),:,:,:)
  
    call self%set (bits%internal)
    self%initialized=.true.
    call trace%end()
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
      self%t_offset = 0.1
    case ('atom_patch')
      self%t_offset = 0.2
    case ('back_opac', 'cont_opac', 'line_opac')
      self%t_offset = 0.3
    case ('net_opac')
      self%t_offset = 0.4
    case ('angle_rad')
      self%t_offset = 0.6
    case ('mean_rad', 'out_rad')
      self%t_offset = 0.8
    case ('rate_patch')
      self%t_offset = 0.9
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
    case ('composition')
      self%maxiter = 0
    case ('atom_patch')
      self%maxiter = maxiter
    case ('rate_patch', 'line_opac')
      self%maxiter = self%m3d_parent%maxiter - 1
    case ('cont_opac')
      self%maxiter = self%m3d_parent%maxiter
    case ('back_opac')
      self%maxiter = ilambd + 1
    case ('mean_rad')
      self%maxiter = ilambd
    case ('angle_rad', 'net_opac')
      self%maxiter = ilambd + maxiter - 1
    end select

    if (self%out_ray) self%maxiter = 0
    call trace%end()
  END SUBROUTINE set_maxiter

  !===============================================================================
  SUBROUTINE m3d_pre_update (self)
    class(m3d_patch_t), target :: self
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%pre_update')
    if (self%finished) then
      print '(A, I4, 3(A))', 'task ', self%id, ' of kind ', trim(self%kind), ' is updating despite finished!'
      print *, 'task%time = ', self%time
    elseif(.not. self%initialized) then
      call self%init()
    endif
    call trace%end()
  END SUBROUTINE m3d_pre_update
  
  !===============================================================================
  SUBROUTINE m3d_post_update (self)
    class(m3d_patch_t) :: self
    character(len=38) :: text
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_post_update')
    if (verbose > 1) then
      write(text, '(F5.1, A12, 4(I4))') self%time, trim(self%kind), self%id, self%ipos(3), &
                                        self%i_nu, self%i_mu
      select case (self%kind)
      case('rate_patch')
        write(*,'(A, 6(ES15.7))') text, &
        self%authority_mem(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%linecheck,1,:)
      case('out_rad')
        write(*,'(A, 6(ES15.7))') text, &
        self%authority_mem(self%pcheck(1),self%pcheck(2),1,self%fcheck,1,:)
      case('composition')
        write(*,'(A, 6(ES15.7))') text, &
        self%authority_mem(self%pcheck(1),self%pcheck(2),self%pcheck(3),:self%i_temp,1,:)
      case default
        write(*,'(A, 6(ES15.7))') text, &
        self%authority_mem(self%pcheck(1),self%pcheck(2),self%pcheck(3),self%fcheck,1,:)
      end select
    endif
  
    if (self%iter == self%maxiter) then
      call self%finish()
    endif
  
    self%iter = self%iter + 1
    call trace%end()
  END SUBROUTINE m3d_post_update
  
  !===============================================================================
  SUBROUTINE m3d_finish (self, rotated)
    class(m3d_patch_t):: self
    class(nbor_t), pointer:: nbor, nbnbor
    class(task_t), pointer:: nbtask, nbnbtask
    logical, optional :: rotated
    logical :: nbors_finished
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%m3d_finish')
    self%finished = .true.
    self%time = io%end_time + self%t_offset * self%dtime
    call self%set(bits%no_download)
    if (present(rotated)) self%rotated = rotated

    nbor => self%link%nbor
    do while (associated(nbor))
      if (nbor%needed) then
        nbtask => nbor%task
        select type (nbtask)
        class is (m3d_patch_t)
        if (nbtask%finished) then
          nbnbor => nbtask%link%nbor
          nbors_finished=.true.
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
            call os%deallocate(nbtask%mem)
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
    call self%task_t%nbor_relations(another, needed, needs_me, download)

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
  SUBROUTINE append_to_task_list(self, task_list)
    class (m3d_patch_t), target :: self
    class (list_t), pointer :: task_list
    class (link_t), pointer :: link
    class (task_t), pointer :: task
    !-----------------------------------------------------------------------------
    call trace%begin ('m3d_patch_t%add_to_task_list')
    if (self%maxiter >= 0) then
      allocate (halo_m3d_t:: link)
      call link%init()
      link%task => self
      self%link => link
      call task_list%append (link)
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
    if (verbose < 3) then
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
    write(io_unit%mpi,'(1x,a,2x,a10,i6,2i4,7x,2l1,2x,3f9.3)') &
      'm3d_patch_t%task_info: kind, id, rank, level, BV, pos', &
      self%kind, self%id, self%rank, self%level, &
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
    call trace%begin ('halo_t%download', itimer=itimer)
                                            call self%lock%set ('halo_t%download')
    call self%copy_nbor_list (self%nbor, tnbors)
                                          call self%lock%unset ('halo_t%download')
    target => self%task
    !.............................................................................
    ! important for dynamic initialisation of patches
    select type (target)
    class is (m3d_patch_t)
      if (.not. target%initialized) call target%init()
    end select
    !.............................................................................
    call self%set_filled (target)
    !-----------------------------------------------------------------------------
    ! loop over nbors
    !-----------------------------------------------------------------------------
    nbor => nbors
    do while (associated(nbor))
      if (nbor%download) then
        if (nbor%task%is_clear (bits%remove)) then
          source => nbor%task
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
    call self%get_filled (target)
    call self%remove_nbor_list (tnbors)
    call trace%end (itimer)
  END SUBROUTINE download
      
  end module m3d_patch_mod
  