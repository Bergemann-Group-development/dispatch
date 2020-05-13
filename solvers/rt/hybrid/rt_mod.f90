!===============================================================================
!> A skeleton for RT

!>omega(0) should also be in the task list, independent on any tasks, just updates
!eos whenever the time comes.

! noot done: pre_update, update
! pre_update -- just take the Q values
! update - if omega==0, just update the rk and src, solve otherwise

! fix boundary conditions

!write a paper about RT and other extras in dispatch

!===============================================================================
MODULE rt_mod
  USE boundaries_mod
  USE io_mod
  USE trace_mod
  USE mpi_mod
  USE mesh_mod
  USE omp_mod
  USE omp_timer_mod
  USE timer_mod
  USE index_mod
  USE gpatch_mod
  !USE radau_mod
  USE link_mod
  USE eos_mod
  USE bits_mod
  USE math_mod
  USE task_mod
  USE rt_components_mod
  USE rt_solver_mod
  USE rt_boundaries_mod
  USE ray_geometry_v3_mod
  USE scaling_mod
  USE kinds_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! RT patch data type, extended from gpatch, with extra RT memory and parameters
  !-----------------------------------------------------------------------------
  type, public, extends(gpatch_t):: rt_t
    class(gpatch_t), pointer   :: patch                                         ! pointer to MHD patch
    class(omega_t), dimension(:), pointer :: omega => null()                    ! different ray directions
    type(rt_boundaries_t)      :: rt_boundaries
    integer                    :: n_lambda = 0                                  ! number of frequency bins
    integer                    :: i_omega=0                                     ! for cross-compatibility with SC solver
    integer                    :: n_omega=0                                     ! number of rt space-angles
    real(8)                    :: eos_time
    real, pointer              :: w_lambda(:)
    real   :: rt_grace=0.0
    real(8):: rt_llc(3)    = 0d0
    real(8):: rt_extent(3) = 0d0
    real(8):: rt_timestep=0.01_8
    real   :: epsilon=1.0
    integer:: n_mu=3, n_phi=4
    logical:: on = .false.
    logical:: vertical_propagation=.false.
    logical:: constant_q = .true.
    procedure(optical_depth_staggered), pointer :: optical_depth => null()
    character(len=16)                 :: solver ='integral'             ! choice of solver
    procedure(solve_integral), pointer:: solve => null()
    real(kind=KindRTVar), pointer              :: tt(:,:,:)    => null()
    real(kind=KindRTVar), pointer              :: rk(:,:,:,:)  => null()
    real(kind=KindRTVar), pointer              :: src(:,:,:,:) => null()
    real(kind=KindRTVar), pointer              :: jj(:,:,:)    => null()
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_omega
    procedure:: init_nbors
    procedure:: pre_update
    procedure:: post_update
    procedure:: update
    procedure:: rt_eos
    procedure:: setup_boundaries
    procedure:: warm_up
    !
    procedure :: get_ds
    procedure :: get_q
    procedure :: optical_depth_centered
    procedure :: optical_depth_staggered
    procedure :: optical_depth_staggered_weighted
    procedure :: ray_vol_map3
    procedure :: ray_vol_map4
    procedure :: solve_integral
    procedure :: solve_feautrier
    procedure :: update_q
    procedure :: nbor_relations
  end type
  integer, save:: verbose=0
  type(rt_t), public:: rt
CONTAINS

!===============================================================================
!> Read parameters
!===============================================================================
SUBROUTINE pre_init (self, link)
  class (rt_t):: self
  class (link_t), pointer :: link
  !-----------------------------------------------------------------------------
  integer, save:: n_mu=3, n_phi=4, n_lambda=0, nt=3
  real,    save:: rt_grace=0.0
  real,    save:: epsilon = 1.0
  real(8), save:: rt_llc(3)    = 0d0
  real(8), save:: rt_extent(3) = 0d0
  real(8), save:: rt_timestep=0.01_8
  integer:: dir, i_mu, i_phi, i_omega, n_omega, iostat, n
  character(len=16), save :: solver = 'integral'
  character(len=16), save :: dtau = 'weighted'
  logical, save:: vertical_propagation=.false.
  logical, save:: constant_q = .true.
  logical, save:: first_time=.true.
  logical, save:: on =.false.
  !-----------------------------------------------------------------------------
  namelist /hybrid_rt_params/ on, nt, rt_llc, rt_extent, rt_timestep, n_mu, n_lambda, n_phi, verbose, &
                       rt_grace, vertical_propagation, constant_q, dtau, solver, epsilon
                       
  !-----------------------------------------------------------------------------
  ! Read parameters
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%pre_init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, hybrid_rt_params)
    if (io%master) write (io%output, hybrid_rt_params)
  end if
  call io%lock%unset
  if (on) then
    if (.not. io%master) verbose=0
    self%on = on
    self%nt = nt
    self%n_lambda = n_lambda
    self%rt_grace = rt_grace
    self%rt_llc = rt_llc
    self%rt_extent = rt_extent
    self%rt_timestep = rt_timestep
    self%vertical_propagation = vertical_propagation
    self%constant_q = constant_q
    self%n_mu  = n_mu
    self%n_phi = n_phi
    self%epsilon = epsilon
    ! If EOS n_lambda is not overridden by input, use the default.
    if (self%n_lambda == 0) self%n_lambda=eos%n_lambda
    !
    self%solver = solver                                                  ! solver selection
    select case(self%solver)
    case('integral')
      self%solve => solve_integral
    case('feautrier')
      self%solve => solve_feautrier
    case default
      print *, 'WARNING: rt solver not recognised. Using integral solver.'
      self%solve => solve_integral
    end select

    select case(trim(dtau))
    case('centered')
      self%optical_depth => optical_depth_centered
    case('staggered')
      self%optical_depth => optical_depth_staggered
    case('weighted')
      self%optical_depth => optical_depth_staggered_weighted
    end select
    !---------------------------------------------------------------------------
    ! Assign the MHD patch and request for two more mem indices
    !---------------------------------------------------------------------------
    self%patch => gpatch%cast2gpatch(link%task)
    if (self%patch%idx%qr <1) then
      call self%patch%idx%request_index(self%patch%idx%qr, self%patch%nv)
      call self%patch%idx%request_index(self%patch%idx%tt, self%patch%nv)
    end if
  end if
  call trace%end()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize the object
!===============================================================================
SUBROUTINE init (self)
  USE link_mod
  class (rt_t):: self
  ! ............................................................................
  class (gpatch_t), pointer :: p
  class (rt_t),     pointer :: omega
  class (link_t),   pointer :: link
  class (mesh_t),   pointer :: mesh
  class (*), pointer :: dummy
  integer:: dir, i_mu, i_phi, iostat, n, nv, nd, ng(3)
  real, dimension(:), allocatable:: mu, w_mu
  logical, save:: rg_time=.true.
  ! ----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('rt_t%init')
  p => self%patch
  !-----------------------------------------------------------------------------
  ! interconnect itself
  !-----------------------------------------------------------------------------
  nullify(self%link)
  allocate(link)
  dummy => p%connect%rt
  select type(dummy)
  class is (task_t)
    link%task => dummy
  end select
  self%link => link
  !-----------------------------------------------------------------------------
  ! Initialize the patch
  !-----------------------------------------------------------------------------
  call self%init_unique                                                         ! get an unique ID
  self%track   = .false.
  self%kind    = 'rt_t'                                                         ! to be used in data_io_mod
  ng=1                                                                          ! rt patches have only one guard zone
  if (ray_geometry%head%ray_list%directionality(1) == 0) then
    nd = 2
  else
    nd = 1
  end if
  nv = ray_geometry%head%n_ray_dirs*nd*self%n_lambda                            ! memory slot stores all heating rates per direction per freq bin
  call self%patch_t%setup(position=p%position, size=p%size, n=p%n, ng=ng, &
                          nt=self%nt, nv = nv, nw = 1)
  !-----------------------------------------------------------------------------
  ! Synchronize with MHD patch
  ! ---------------------------------------------------------------------------- 
  self%box        = p%box
  self%size       = p%size
  self%origin     = p%origin
  self%llc_cart   = p%llc_cart
  self%llc_nat    = p%llc_nat
  self%restart    = p%restart
  self%time       = p%time
  self%it         = max(self%nt-1,1)
  self%t(self%it) = self%time
  self%new        = self%nt
  self%nw         = 1
  self%out_next   = (int(self%time/io%out_time+1.0e-4)+1)*io%out_time           ! add an `eps` to avoid round-off errors
  self%dt_fixed   = self%rt_timestep                                            ! Fixed timestep for RT
  self%iout       = self%restart+1
  self%grace      = self%rt_grace
  self%periodic   = p%periodic
  !-----------------------------------------------------------------------------
  ! Check if this MHD patch needs RT
  !-----------------------------------------------------------------------------
  call self%setup_boundaries()
  if(self%on) then
    !---------------------------------------------------------------------------
    ! Allocate heating for MHD
    !---------------------------------------------------------------------------
    if (.not.associated(self%patch%heating_per_unit_volume)) then
      allocate(self%patch%heating_per_unit_volume(p%gn(1), p%gn(2), p%gn(3)))
      self%patch%heating_per_unit_volume = 0.0
    end if
    !---------------------------------------------------------------------------
    ! Allocate wights and omegas. 
    ! Integration weights (w_bin is baked into the source function in this case)
    !---------------------------------------------------------------------------
    allocate(self%w_lambda(self%n_lambda))
    if (self%n_lambda == 1) then
      self%w_lambda = 1.0
    else
      self%w_lambda = eos%w_lambda
    end if
    self%n_omega = ray_geometry%head%n_ray_dirs
    allocate (self%omega(self%n_omega))
    call self%init_omega()
    self%dt_fixed = self%rt_timestep
    !
    allocate (self%unsigned(self%nv), self%pervolume(self%nv))
    self%unsigned  = .false.                                                    ! heating rates can be both positive and negative
    self%pervolume = .false.                                                    ! although heating is per unit volume, there is no density memory slot in RT patches, thus this has to be .false.
    call self%set (bits%no_density)
    ! --------------------------------------------------------------------------
    ! Allocate memory
    ! --------------------------------------------------------------------------
    if (associated(self%mem)) deallocate (self%mem)
    nullify (self%mem)
    allocate (self%mem (self%gn(1),&
                        self%gn(2),&
                        self%gn(3),&
                        self%nv,self%nt,1))
    call io%bits_mem(storage_size(self%mem)  ,product(shape(self%mem))   ,'q')
    self%mem_allocated = .true.                                                 ! mark
    self%mem = 0.0 
    !---------------------------------------------------------------------------
    ! Allocate EOS data
    !---------------------------------------------------------------------------
    allocate (self%rk (self%gn(1),self%gn(2),self%gn(3),self%n_lambda))
    allocate (self%tt (self%gn(1),self%gn(2),self%gn(3)))
    allocate (self%src(self%gn(1),self%gn(2),self%gn(3),self%n_lambda))
    allocate (self%jj (self%gn(1),self%gn(2),self%gn(3)))
    call io%bits_mem(storage_size(self%tt  ),product(shape(self%tt )),'tt' )
    call io%bits_mem(storage_size(self%rk  ),product(shape(self%rk )),'rk' )
    call io%bits_mem(storage_size(self%src ),product(shape(self%src)),'src')
    call io%bits_mem(storage_size(self%jj  ),product(shape(self%jj )),'jj' )
  else
    call self%dealloc()
  end if
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Init different ray directions
!===============================================================================
SUBROUTINE init_omega(self)
  class(rt_t) :: self
  ! ............................................................................
  class (ray_dirs_t), pointer :: ray_dir
  class (ray_bndl_t), pointer :: bundle
  class (omega_t),    pointer :: omega
  integer :: a_id, ifl, irl
  ! ----------------------------------------------------------------------------
  a_id = 1                                                                      ! omega personal ID in the omega array
  ifl  = 1
  irl  = ray_geometry%head%n_ray_dirs*eos%n_lambda + 1                          ! where reverse ray data starts in mem()
  ray_dir => ray_geometry%head%ray_list%head                                    ! pick a ray direction
  do while(associated(ray_dir))                                                 ! loop over all ray directions
    bundle          => ray_dir%bundle(1)                                        ! there is a single bundle for this solver
    omega           => self%omega(a_id)
    omega%link      => self%link                                                ! might be useful
    omega%direction =  0                                                        ! bi-directional
    omega%id        =  a_id                                                     ! omega personal ID in the omega array
    omega%bundle    => bundle
    omega%ifl       =  ifl
    omega%ifu       =  ifl + self%n_lambda - 1
    omega%irl       =  irl
    omega%iru       =  irl + self%n_lambda - 1
    call self%get_ds(bundle,omega%ds)                                           ! get distance between ray points
    !
    a_id    = a_id+1
    ifl     = ifl + self%n_lambda
    irl     = irl + self%n_lambda
    ray_dir => ray_dir%next                                                     ! next direction
  end do
  !
  nullify(ray_dir, bundle, omega)
END SUBROUTINE init_omega


!===============================================================================
!> Add upstream and downstream tasks to an nbor list.  The upstream nbors should
!> have the settings needed=T, needs_me=F, download=T, while the downstream 
!> tasks should have needed=F, needs_me=T, download=F.
!>
!> Each omega subtask of an MHD task has, for each MHD nbor, either
!> 1) no connection, or
!> 2) an upstream connection, or
!> 3) a downstream connection
!>
!> The nbor's sub-task in the omega direction should not be added to then MHD
!> task's nbor list, and should, in the local sub-tasks nbor list, either
!> 1) not be added
!> 2) be added with needed=T, needs_me=F, download=T, or
!> 3) be added with needed=F, needs_me=t, download=F
!>
!===============================================================================
SUBROUTINE init_nbors (self)
  class (rt_t):: self
  !.............................................................................
  type(rt_t), pointer:: omega
  class(*), pointer:: nb_rt_task
  class(link_t), pointer:: nbor, new
  class(task_t), pointer:: nbtask
  class(gpatch_t), pointer:: patch, nb_patch
  integer:: n_nbor
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('rt_t%init_nbors')
  patch => self%patch
  call self%clear(bits%init_nbors)
  n_nbor = 0
  nbor => patch%link%nbor
  do while (associated(nbor))
    !-----------------------------------------------------------------------
    ! Add a new nbor, a new nbor%link, and pointers from link to task and
    ! from task to link.  The omega-task on nbor already has its own link.
    !-----------------------------------------------------------------------
    nbtask => nbor%task
    select type (nbtask)
    class is (gpatch_t)
      nb_patch   => nbtask
      nb_rt_task => nbtask%connect%rt
      if (associated(nb_rt_task)) then
        select type(nb_rt_task)
        class is (rt_t)
          if (nb_rt_task%on) then
            !-------------------------------------------------------------------
            ! Create an "upstrean type" nbor node pointing to the nbor sub-task
            !-------------------------------------------------------------------
            allocate (new)
            new%link => nb_rt_task%link
            new%task  => nb_rt_task
            ! -- add to the sub-tasks nbor list, in MPI rank order
            !call self%link%add_nbor_by_rank (self%link%nbors_by_level, new)
            call self%link%add_nbor_by_rank (self%link%nbor, new)
            !-------------------------------------------------------------------
            ! Create a "downstrean type" nbor node pointing back to omega
            !-------------------------------------------------------------------
            allocate (new)
            new%link => self%link
            new%task  => self%link%task
            !! -- add it to the nbor's nbor chain, in rank order
            !call nb_rt_task%link%add_nbor_by_rank (nb_rt_task%link%nbors_by_level, new)
            call nb_rt_task%link%add_nbor_by_rank (nb_rt_task%link%nbor, new)
            n_nbor = n_nbor+1
          end if
        end select
      end if
    end select
    nbor => nbor%next
  end do
    if (verbose > 2) &
      write(io%output,*) self%id, 'rt_t%init_nbors: n_nbor=', n_nbor
  call trace%end()
END SUBROUTINE init_nbors


!===============================================================================
!> Check if patch belongs to an RT box or is at its boundary.
!===============================================================================
SUBROUTINE setup_boundaries(self)
  class (rt_t) :: self
  ! ----------------------------------------------------------------------------
  real(8) :: llc(3), urc(3), rt_urc(3)
  integer :: lbit(3), ubit(3), i
  ! ............................................................................
  call self%init_bdries
  llc = self%position-0.5_8*self%size
  urc = self%position+0.5_8*self%size
  rt_urc = self%rt_llc+self%rt_extent
  !
  lbit = [bits%xl, bits%yl, bits%zl]                                            ! lower boundary bits -- convenience
  ubit = [bits%xu, bits%yu, bits%zu]                                            ! upper boundary bits -- convenience
  !
  do i = 1,3
    if (self%rt_extent(i) > 0.) then
      if ((urc(i) < self%rt_llc(i))) then
        call self%patch%set(bits%no_rt)
        self%on = .false.
      else if (llc(i) <= self%rt_llc(i)) then
        call self%boundaries%set (lbit(i))
      end if
      !
      if ((llc(i) > rt_urc(i))) then
        call self%patch%set(bits%no_rt)
        self%on = .false.
      else if (urc(i) >= rt_urc(i)) then
        call self%boundaries%set (ubit(i))
      end if
    end if
  end do
  !
  call self%rt_boundaries%init(self%mesh)
END SUBROUTINE setup_boundaries

!===============================================================================
!>
!===============================================================================
SUBROUTINE pre_update (self)
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(rt_t), pointer:: omega
  real, dimension(:,:,:), pointer:: q
  integer, save:: itimer=0
  real(8):: time
  integer:: l(3), u(3)
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  !
  call trace%begin ('rt_t%pre_update', itimer=itimer)
  l = self%patch%mesh%lo
  u = self%patch%mesh%uo
  q => self%patch%heating_per_unit_volume(l(1):u(1),l(2):u(2),l(3):u(3))
  if (self%constant_q) then
    call self%get_q(time, q)
    !print *, 'minmax heating', minval(q), maxval(q)
    !print *, 'minmax loc heating', minloc(q), maxloc(q)
    !print *, 'after Q retrieval for heating'
    !print *, self%patch%heating_per_unit_volume(4,l(2):u(2),4)
  else
    call io%abort("time-interpolated q in rt_t%pre_update is not implemented yet")
  end if
  call trace%end(itimer)
END SUBROUTINE pre_update

!=======================================================================
!> Get heating rates for given time in terms of the already averaged <Q>.
!> Scattering reduces the heat transfer in proportion to epsilon.
!=======================================================================
SUBROUTINE get_q (self, time, q)
  class (rt_t) :: self
  real(8):: time
  real(kind=KindScalarVar), dimension(:,:,:) :: q
  !.....................................................................
  integer :: iv
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rt_t%get_q', itimer=itimer)
  associate (l => self%li, u => self%ui)
  !q = q + (self%jj (l(1):u(1),l(2):u(2),l(3):u(3)) &
  !         - self%src(l(1):u(1),l(2):u(2),l(3):u(3),1)) &
  !        * self%rk(l(1):u(1),l(2):u(2),l(3):u(3),1) * math%pi4
  do iv = 1, self%n_lambda
    q = q + (self%jj - self%src(:,:,:,iv)) * self%rk(:,:,:,iv) * math%pi4
  end do
  end associate
  call trace%end (itimer)
END SUBROUTINE get_q

!===============================================================================
!>
!===============================================================================
SUBROUTINE post_update (self)
  class (rt_t) :: self
  !-----------------------------------------------------------------------------
  class(rt_t), pointer:: omega
  real, dimension(:,:,:), pointer:: q
  integer, save:: itimer=0
  real(8):: time
  integer:: l(3), u(3), lp(3), up(3)
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  !
  call trace%begin ('rt_t%post_update', itimer=itimer)
  l=self%li
  u=self%ui
  lp=self%patch%mesh%lo
  up=self%patch%mesh%uo
  if (self%constant_q) then
    q => self%patch%mem(lp(1):up(1),lp(2):up(2),lp(3):up(3),self%patch%idx%qr,self%patch%new,1)
    q = 0.0
    call self%get_q(time, q)
    !print *, 'after Q retrieval for printing'
    !print *, q(2,:,2)
    !print *, 'minmax qr', minval(q), maxval(q)
    !print *, 'minmax qr loc', minloc(q), maxloc(q)
    self%patch%mem(lp(1):up(1),&
                   lp(2):up(2),&
                   lp(3):up(3),&
                   self%patch%idx%tt,self%patch%new,1) = self%tt
  else
    call io%abort("rt_t:: post_update requested non constant q")
  end if
  call trace%end(itimer)
END SUBROUTINE post_update

!===============================================================================
!> Run through all angular solvers and attempt to update them.  Solvers that
!> lack upstream info (based on internal criteria) will do nothing, while
!> solvers that have enough information will be updated.  The angular solvers
!> are so lightweight (taking of the order 100-200 ns) that it would be too
!> much overhead to consider them separate tasks.
!>
!> It would be easy to spawn independent OMP tasks for the ones that are
!> ready, and then wait on them, but for now we avoid complications.
!>
!> The update step computes the current RT patch time, based on the current
!> time in the attached RT solvers for each angle
!===============================================================================
SUBROUTINE update (self)
  class (rt_t):: self
  !.............................................................................
  class(omega_t),    pointer :: omega
  class(ray_bndl_t), pointer :: bundle
  class(packet_t),   pointer :: packet
  integer:: i_omega, i_packet, i_bin, ir, nrays, npt, iv, c(3), l, u
  real(kind=KindRTVar), dimension(:,:  ), pointer :: q_fw0, q_re0
  real(kind=KindRTVar), dimension(:,:,:), pointer :: q_fw , q_re
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  !
  call trace%begin ('rt_t%update', itimer=itimer)
  !
  if (self%id == io%id_debug) &
    write (io%output,*) 'DBG rt_t%udpate: time =', self%time, &
    self%is_set(bits%boundary), self%is_set(bits%virtual)
  !
  call self%rt_eos(self%time)                                                   ! Update EOS for RT
  !
  associate(mem => self%mem(:,:,:,:,self%it,1))
  !
  !print *, 'mem before BC'
  !print *, self%mem(1,:,2,1,self%it,1)
  !
  ! ----------------------------------------------------------------------------
  ! reset <Q>
  ! ----------------------------------------------------------------------------
  self%jj = 0.0
  !
  !print *, 'mem after BC'
  !print *, self%mem(1,:,2,1,self%it,1)
  ! ----------------------------------------------------------------------------
  ! Go through omega and update RT
  ! ---------------------------------------------------------------------------- 
  do i_omega = 1, self%n_omega
    omega => self%omega(i_omega)
    bundle => omega%bundle
    ! --------------------------------------------------------------------------
    ! Update boundary conditions
    ! --------------------------------------------------------------------------
    ! forward direction
    l = omega%ifl
    u = omega%ifu
    call self%rt_boundaries%condition (mem=mem(:,:,:,l:u), src=self%src, &
                        rk=self%rk, tt=self%tt, mu=bundle%mu,   w=bundle%weight)
    ! reverse direction
    l = omega%irl
    u = omega%iru
    call self%rt_boundaries%condition (mem=mem(:,:,:,l:u), src=self%src, &
                        rk=self%rk, tt=self%tt, mu=bundle%mu,   w=bundle%weight)
    ! --------------------------------------------------------------------------
    ! go through packets in this bundle
    ! --------------------------------------------------------------------------
    do i_packet = 1, bundle%n_packets
      packet => bundle%packet(i_packet)                                         ! pick up a ray packet
      nrays  =  packet%n_rays                                                   ! number of rays in this packet
      npt    =  packet%n_points+2                                               ! +2 boundary points
      ! ------------------------------------------------------------------------
      ! allocate ray scratches and set them to 0
      ! ------------------------------------------------------------------------
      allocate(q_fw0 (nrays,eos%n_lambda))
      allocate(q_re0 (nrays,eos%n_lambda))
      allocate(q_fw (nrays,eos%n_lambda,npt))
      allocate(q_re (nrays,eos%n_lambda,npt))
      q_fw0=0.0
      q_re0=0.0
      q_fw =0.0
      q_re =0.0
      !

      !print *, 'update: reading mapped qfw0'
      ! ------------------------------------------------------------------------
      ! map boundary heating rates in forward direction
      ! ------------------------------------------------------------------------
      do ir = 1, nrays
        !q_fw0(ir,:) = omega%pdata(ir)%in(:,self%it)
        c= packet%rev_coords(:,ir)
        !print *, ir, omega%pdata(ir)%in(:,:),q_fw0(ir,:) , 'mem',self%mem(c(1),c(2),c(3),:,:,1)
        q_fw0(ir,:) = mem(c(1),c(2),c(3),omega%ifl:omega%ifu)
        !print *, ir, omega%pdata(ir)%in(:,:),q_fw0(ir,:) , 'mem',self%mem(c(1),c(2),c(3),:,:,1)
      end do
      ! ------------------------------------------------------------------------
      ! map boundary heating rates in reverse direction
      ! ------------------------------------------------------------------------
      do ir = 1, nrays
        !q_re0(tr,:) = omega%pdata(ir)%in(:,self%it)
        c= packet%fwd_coords(:,ir)
        q_re0(ir,:) = mem(c(1),c(2),c(3),omega%irl:omega%iru)
      end do
      !
      !print *, 'qfw0'
      !print *, q_fw0(33:65,:)
      !print *, 'before update, RT mem bnd'
      !print *, self%mem(1,:,2,omega%irl:omega%iru,self%it,1)
      !
      ! ------------------------------------------------------------------------
      ! solve the radiative transfer equation for this particular link
      ! ------------------------------------------------------------------------
      call self%solve(omega, packet,q_in_fw=q_fw0,  q_in_re=q_re0, &
                           q_out_fw=q_fw,  q_out_re=q_re)          
      !print *, 'omega q out fw', shape(q_fw)
      !print *, q_fw(33:65,:,2)
      !
      ! ------------------------------------------------------------------------
      ! update the heating rates, obtained from this packet
      ! ------------------------------------------------------------------------
      call self%update_q(omega, packet, q_fw=q_fw, q_re=q_re, it = self%new)
      !print *, 'after update, RT mem'
      !print *, self%mem(2,:,2,1,self%it,1)
      deallocate (q_fw, q_re, q_fw0, q_re0)
    end do
    ! --------------------------------------------------------------------------
    ! Compute the updated QQ = <Q>
    ! --------------------------------------------------------------------------
    do iv=1, self%n_lambda
      self%jj = self%jj + mem(:,:,:,omega%ifl+iv-1) &
                        * bundle%weight * self%w_lambda(iv) / math%pi4
      self%jj = self%jj + mem(:,:,:,omega%irl+iv-1) &
                        * bundle%weight * self%w_lambda(iv) / math%pi4
    end do
  end do
  do iv=1, self%n_lambda
    self%jj = self%jj + self%src(:,:,:,iv)
  end do
  end associate
  !
  ! ----------------------------------------------------------------------------
  ! Update timestep
  ! ----------------------------------------------------------------------------
  self%dtime = self%dt_fixed
  call self%rotate()
  !
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Compute the RT EOS quantities for the patch, using self%eos_time to guard
!> againts duplicate calculations when multi-threading.  Check both before and
!> after entering a critical region, so most threads never hit the critical region.
!===============================================================================
SUBROUTINE rt_eos (self, time)
  class (rt_t):: self
  real(8):: time
  !.............................................................................
  class (gpatch_t), pointer :: patch
  integer, save:: itimer=0
  real, dimension(:,:,:), pointer:: d, e, s, e1, k
  integer:: m(4), ib, i3, jt(2), l(3), u(3), i, pl(3), pu(3)
  real:: g1
  real:: z, pt(2)
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('rt_t%rt_eos', itimer=itimer)
  patch => self%patch
  ! ----------------------------------------------------------------------------
  ! interpolate to the local time
  ! ----------------------------------------------------------------------------
  call patch%time_interval (time, jt, pt, all=.true.)
  !  
  m = shape(self%rk)
  l = self%mesh%lb
  u = self%mesh%ub
  pl = patch%mesh%lo
  pu = patch%mesh%uo
  ! ----------------------------------------------------------------------------
  ! allocate scratches
  ! ----------------------------------------------------------------------------
  allocate (d(m(1),m(2),m(3)))                                                  ! density scratch
  allocate (e(m(1),m(2),m(3)))                                                  ! internal energy scratch
  d = patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(1),1)*pt(1) &
    + patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%d,jt(2),1)*pt(2)
  ! ----------------------------------------------------------------------------
  ! Take internal energy and density from MHD patch
  ! ----------------------------------------------------------------------------
  if (patch%kind(1:9)=='stagger2e') then
    e = patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%e,jt(1),1)*pt(1) &
      + patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%e,jt(2),1)*pt(2)
    e = e/d
  else if (patch%kind(1:9)=='stagger2_') then
    associate(s1 => patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%s,jt(1),1), &
              s2 => patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%s,jt(2),1))
    allocate (s(m(1),m(2),m(3)))
    s = s1*pt(1) + s2*pt(2)
    g1 = patch%gamma-1.0
    e = d**g1*exp(s/d*g1)/g1
    deallocate(s)
    end associate
  else if (patch%kind(1:7)=='ramses_') then
    associate(e1 => patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%e           ,jt(1),1), &
              e2 => patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%e           ,jt(2),1), &
              p1 => patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%px:patch%idx%pz,jt(1),1),&
              p2 => patch%mem(pl(1):pu(1),pl(2):pu(2),pl(3):pu(3),patch%idx%px:patch%idx%pz,jt(2),1))
    allocate(k(m(1),m(2),m(3)))
    k = 0.0
    e = e1*pt(1)+e2*pt(2)
    do i=1,size(p1,4)
      k = k+(p1(:,:,:,i)*pt(1)+p2(:,:,:,i)*pt(2))**2
    end do
    k = 0.5*k/d
    e = (e-k)/d
    deallocate(k)
    end associate
  else
    print *,trim(patch%kind)
    error stop "rt_mod (1): MHD solver type not known"
  end if
  ! ----------------------------------------------------------------------------
  ! Lookup the needed RT EOS quantities
  ! ----------------------------------------------------------------------------
  call eos%lookup_table (m, d=d, ee=e, rk=self%rk, src=self%src, tt=self%tt)
  !
  if (verbose > 3) then
    write(io%output,*) patch%id,'rt_t%rt_eos:  d',minval(d),maxval(d)
    write(io%output,*) patch%id,'rt_t%rt_eos:  e',minval(e),maxval(e)
    write(io%output,*) patch%id,'rt_t%rt_eos: rk',minval(self%rk),maxval(self%rk)
    write(io%output,*) patch%id,'rt_t%rt_eos: tt',minval(self%tt),maxval(self%tt),scaling%temp
  end if
  !
  deallocate (d, e)
  self%eos_time = time
  !> -------------------------------------------------------------------
  !> Combine scattering mean <Q> with thermal source part
  !> -------------------------------------------------------------------
  do i = 1, self%n_lambda
    self%src(:,:,:,i) = (1.-self%epsilon)*self%jj + self%epsilon*self%src(:,:,:,i)
  end do 
  call trace%end (itimer)
END SUBROUTINE rt_eos

!===============================================================================
!>
!===============================================================================
SUBROUTINE update_q (self, omega, packet, q_fw, q_re, it)
  class(rt_t) :: self
  class(omega_t), pointer  :: omega
  class(packet_t), pointer :: packet
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_fw, q_re
  integer:: it
  ! ....................................................................
  real, dimension(:,:,:,:), pointer :: qqf, qqr
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), wl, ray, pt, ib, ifl, ifu, irl, iru, nwb
  integer, save:: itimer=0
  real(8):: wc
  !---------------------------------------------------------------------
  call trace_begin ('rt_t%update_q', 4, itimer=itimer)
  ifl = omega%ifl
  ifu = omega%ifu
  irl = omega%irl
  iru = omega%iru
  nwb = self%n_lambda
  !
  qqf  => self%mem(:,:,:,ifl:ifu ,it,1)                                          ! point to the forward direction memory slot
  qqr => self%mem(:,:,:,irl:iru,it,1)                                           ! point to the reverse direction memory slot
  cell => packet%coords                                                         ! point to coordinates
  do ray = 1, packet%n_rays                                                     ! loop over rays
    do pt = 1, packet%n_points                                                  ! loop over points
      p = cell(:,pt,ray)                                                        ! get the cell position
      do wl=1,nwb
        qqf(p(1),p(2),p(3),wl) = q_fw(ray,wl,pt+1)
        qqr(p(1),p(2),p(3),wl) = q_re(ray,wl,pt+1)
      end do
    end do
  end do
  
  call trace_end (itimer)
END SUBROUTINE update_q

!===============================================================================
!> Prepare and solve the RTE using the integral formulation
!===============================================================================
SUBROUTINE solve_integral (self, omega, packet, q_in_fw, q_in_re, q_out_fw, q_out_re)
  class(rt_t) :: self
  class(omega_t) :: omega
  class(packet_t):: packet
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_fw
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_re
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_fw
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_re
  !.....................................................................
  real(kind=KindRTVar), dimension(:,:,:), pointer :: s, rk, dtau, exp0
  integer                         :: nrays, npt, nwb, ir
  integer, save                   :: itimer=0
  ! ----------------------------------------------------------------------------
  call trace_begin ('rt_t%solve_integral', itimer=itimer)
  nwb   =  self%n_lambda                                                        ! number of wl bins
  nrays =  packet%n_rays                                                        ! number of rays in this packet
  
  ! ---------------------------------------------------------------------------
  ! RT solvers need at least 3 points to work with, as they use
  ! derivatives. Rays, shorter, than 3 points, have to be extended
  ! to patch boundaries to have accurate solutions.
  ! ----------------------------------------------------------------------------
  npt = packet%n_points + 2
  ! ----------------------------------------------------------------------------
  ! allocate all the scratches
  ! ----------------------------------------------------------------------------
  allocate(s     (nrays, nwb, npt))                                             ! source function
  allocate(rk    (nrays, nwb, npt))                                             ! opacity
  allocate(dtau  (nrays, nwb, npt))                                             ! optical depth between the points
  allocate(exp0  (nrays, nwb, npt))                                             ! first exponential factor
  dtau = 0.0
  ! ----------------------------------------------------------------------------
  ! map the quantities onto rays
  ! ----------------------------------------------------------------------------
  call self%ray_vol_map4(packet, self%src, s)                                   ! map source function to rays
  call self%ray_vol_map4(packet, self%rk, rk)                                   ! map opacities to rays
  call self%optical_depth(omega%ds, rk, dtau)                                   ! calculate dtau
  !print *, 'solve integral:'
  !print *, 'nwb, nrays, npt, ds, minmax dtau: ', nwb, nrays, npt, omega%ds, minval(dtau), maxval(dtau)
#if defined __INTEL_COMPILER
  exp0(:,:,2:) = 2**(-dtau(:,:,2:)*1.44269504)
#else
  exp0(:,:,2:) = exp(-dtau(:,:,2:))                                             ! expopential factors
#endif
  ! ----------------------------------------------------------------------------
  ! Solve the RTE using the integral formulation
  ! ----------------------------------------------------------------------------
  call rt_solver%integral_double(dtau, exp0, S, q_in_fw, q_in_re, q_out_fw, q_out_re)
  !
  deallocate (s, rk, dtau, exp0)
  call trace_end (itimer)
END SUBROUTINE solve_integral

!=======================================================================
!> Same as solve_integral, but for the Feautrier formulation
!=======================================================================
SUBROUTINE solve_feautrier (self, omega, packet, q_in_fw, q_in_re, q_out_fw, q_out_re)
  class(rt_t) :: self
  class(omega_t) :: omega
  class(packet_t), pointer :: packet
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_fw
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_re
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_fw
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_re
  !.....................................................................
  call io%abort("rt_t::feautrier solver not updated yet")
END SUBROUTINE solve_feautrier

!===============================================================================
!> "warm up" the simulation by only running the RT tasks without updating their
!> timestep. The only difference from "update" is that warm_up doesn't update EOS,
!> doesn't update the RT task time and doesn't rotate the memory slots.
!===============================================================================
SUBROUTINE warm_up(self)
  class (rt_t):: self
  !.............................................................................
  class(omega_t),    pointer :: omega
  class(ray_bndl_t), pointer :: bundle
  class(packet_t),   pointer :: packet
  integer:: i_omega, i_packet, i_bin, ir, nrays, npt, iv, c(3), l, u
  real(kind=KindRTVar), dimension(:,:  ), pointer :: q_fw0, q_re0
  real(kind=KindRTVar), dimension(:,:,:), pointer :: q_fw , q_re
  real(8):: time
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('rt_t%warm_up', itimer=itimer)
  !
  associate(mem => self%mem(:,:,:,:,self%it,1))
  ! ----------------------------------------------------------------------------
  ! reset <Q>
  ! ----------------------------------------------------------------------------
  self%jj = 0.0
  ! ----------------------------------------------------------------------------
  ! Go through omega and update RT
  ! ---------------------------------------------------------------------------- 
  do i_omega = 1, self%n_omega
    omega => self%omega(i_omega)
    bundle => omega%bundle
    ! --------------------------------------------------------------------------
    ! Update boundary conditions
    ! --------------------------------------------------------------------------
    ! forward direction
    l = omega%ifl
    u = omega%ifu
    call self%rt_boundaries%condition (mem=mem(:,:,:,l:u), src=self%src, &
                        rk=self%rk, tt=self%tt, mu=bundle%mu,   w=bundle%weight)
    ! reverse direction
    l = omega%irl
    u = omega%iru
    call self%rt_boundaries%condition (mem=mem(:,:,:,l:u), src=self%src, &
                        rk=self%rk, tt=self%tt, mu=bundle%mu,   w=bundle%weight)
    ! --------------------------------------------------------------------------
    ! go through packets in this bundle
    ! --------------------------------------------------------------------------
    do i_packet = 1, bundle%n_packets
      packet => bundle%packet(i_packet)                                         ! pick up a ray packet
      nrays  =  packet%n_rays                                                   ! number of rays in this packet
      npt    =  packet%n_points+2                                               ! +2 boundary points
      ! ------------------------------------------------------------------------
      ! allocate ray scratches and set them to 0
      ! ------------------------------------------------------------------------
      allocate(q_fw0 (nrays,eos%n_lambda))
      allocate(q_re0 (nrays,eos%n_lambda))
      allocate(q_fw (nrays,eos%n_lambda,npt))
      allocate(q_re (nrays,eos%n_lambda,npt))
      q_fw0=0.0
      q_re0=0.0
      q_fw =0.0
      q_re =0.0
      ! ------------------------------------------------------------------------
      ! map boundary heating rates in forward direction
      ! ------------------------------------------------------------------------
      do ir = 1, nrays
        c= packet%rev_coords(:,ir)
        q_fw0(ir,:) = mem(c(1),c(2),c(3),omega%ifl:omega%ifu)
      end do
      ! ------------------------------------------------------------------------
      ! map boundary heating rates in reverse direction
      ! ------------------------------------------------------------------------
      do ir = 1, nrays
        c= packet%fwd_coords(:,ir)
        q_re0(ir,:) = mem(c(1),c(2),c(3),omega%irl:omega%iru)
      end do
      ! ------------------------------------------------------------------------
      ! solve the radiative transfer equation for this particular link
      ! ------------------------------------------------------------------------
      call self%solve(omega, packet,q_in_fw=q_fw0,  q_in_re=q_re0, &
                           q_out_fw=q_fw,  q_out_re=q_re)          
      ! ------------------------------------------------------------------------
      ! update the heating rates, obtained from this packet
      ! ------------------------------------------------------------------------
      call self%update_q(omega, packet, q_fw=q_fw, q_re=q_re, it = self%it)
      deallocate (q_fw, q_re, q_fw0, q_re0)
    end do
    ! --------------------------------------------------------------------------
    ! Compute the updated QQ = <Q>
    ! --------------------------------------------------------------------------
    do iv=1, self%n_lambda
      self%jj = self%jj + mem(:,:,:,omega%ifl+iv-1) &
                        * bundle%weight * self%w_lambda(iv) / math%pi4
      self%jj = self%jj + mem(:,:,:,omega%irl+iv-1) &
                        * bundle%weight * self%w_lambda(iv) / math%pi4
    end do
  end do
  do iv=1, self%n_lambda
    self%jj = self%jj + self%src(:,:,:,iv)
  end do
  end associate
  !
  call trace%end (itimer)
END SUBROUTINE warm_up

!===============================================================================
!> calculate the distance between two ray points
!===============================================================================
SUBROUTINE get_ds(self, bundle, ds)
  class(rt_t) :: self
  class(ray_bndl_t), pointer :: bundle
  real(kind=8) :: cs(3), ds
  !.....................................................................
  real(kind=8) :: dst(3), c(3)
  !---------------------------------------------------------------------
  dst = real(bundle%dstep,kind=8)
  c   = self%mesh%d
  ds  = ((c(1)*dst(1))**2+(c(2)*dst(2))**2+(c(3)*dst(3))**2)**0.5
END SUBROUTINE get_ds

!=======================================================================
!> Compute the optical depth along rays; cell-centred version
!=======================================================================
SUBROUTINE optical_depth_centered (self, ds, rk, dtau)
  class (rt_t) :: self
  real(kind=8) :: ds
  real(kind=KindRTVar), dimension(:,:,:), pointer :: rk
  real(kind=KindRTVar), dimension(:,:,:), pointer :: dtau
  !.....................................................................
  integer                         :: pt
  !---------------------------------------------------------------------
  do pt = 1, size(rk,3)
    dtau(:,:,pt) = ds*rk(:,:,pt)
  end do
END SUBROUTINE optical_depth_centered

!=======================================================================
!> Compute the optical depth along rays; face-centred version
!=======================================================================
SUBROUTINE optical_depth_staggered (self, ds, rk, dtau)
  class (rt_t) :: self
  real(kind=8) :: ds
  real(kind=KindRTVar), dimension(:,:,:), pointer :: rk
  real(kind=KindRTVar), dimension(:,:,:), pointer :: dtau
  !.....................................................................
  integer                         :: pt
  !---------------------------------------------------------------------
  dtau(:,:,1) = 0.0
  do pt = 2, size(rk,3)
    dtau(:,:,pt) = ds * 0.5 * (rk(:,:,pt-1) + rk(:,:,pt))
  end do
END SUBROUTINE optical_depth_staggered

!=======================================================================
!> Compute the optical depth along rays. Dtau is staggered, therefore it
!> has only a half value on the lower edge of the patch.
!> Based on fig.4 in Howell et al. (2001) or Ramsey & Dullemond (2014), eq.9
!=======================================================================
SUBROUTINE optical_depth_staggered_weighted (self, ds, rk, dtau)
  class (rt_t) :: self
  real(kind=8) :: ds
  real(kind=KindRTVar), dimension(:,:,:), pointer :: rk
  real(kind=KindRTVar), dimension(:,:,:), pointer :: dtau,dt
  !.....................................................................
  integer                         :: pt, n
  !---------------------------------------------------------------------
  dtau(:,:,1) = 0.0
  do pt=2,size(rk,3)
    dtau(:,:,pt) = min((rk(:,:,pt-1)+rk(:,:,pt))*0.5, &
                   max(2*rk(:,:,pt-1)*rk(:,:,pt)/(rk(:,:,pt-1)+rk(:,:,pt)), &
                   1.33333/ds)) *ds
  end do
END SUBROUTINE optical_depth_staggered_weighted

!=======================================================================
!> Map cell data onto rays
!=======================================================================
SUBROUTINE ray_vol_map4 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t):: packet
  real(kind=KindRTVar), dimension(:,:,:)   , pointer :: out
  real(kind=KindRTVar), dimension(:,:,:,:) , pointer :: in
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), p2(3), ray, pt, nwb, wl, lp
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rt_t%ray_vol_map4', itimer=itimer)
  cell => packet%coords
  nwb = size(out,2)
  lp  = size(out,3)
  !> -------------------------------------------------------------------
  !> we need to consider the first and last points separately,
  !> as they do not belong to these rays implicitly.
  !> -------------------------------------------------------------------
  do ray = 1, packet%n_rays                                             ! loop over all rays
    do pt = 1, packet%n_points                                          ! loop over all points
      p = cell(:,pt,ray)                                                ! get the right cell coordinate
      do wl =1, nwb                                                     ! loop over wl bins
        out(ray,wl,pt+1) = in(p(1),p(2),p(3),wl)                        ! first point is missing
      end do
    end do
    p  = packet%rev_coords(:,ray)
    p2 = packet%fwd_coords(:,ray)
    do wl=1, nwb                                                        ! fill in the first and last points
      out(ray,wl,1)  = in(p(1),p(2),p(3), wl)                           ! Get the first point
      out(ray,wl,lp) = in(p2(1),p2(2),p2(3), wl)                        ! Get the last point
    end do
  end do
  if (io%verbose > 1) then
    if (any(in<0)) then
      print *, in
      error stop "ray vol map4 in<0"
    end if
    if (any(out<0)) then
      do pt = 1, size(out,3)
      print *,"pt",pt, out(:,:,pt)
      end do
      error stop "ray vol map4 out<0"
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE ray_vol_map4

!=======================================================================
!> Map cell data onto rays. In this case, it is wavelenght independent.
!=======================================================================
SUBROUTINE ray_vol_map3 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t):: packet
  real(kind=KindRTVar), dimension(:,:)   , pointer, intent(out) :: out
  real(kind=KindRTVar), dimension(:,:,:) , pointer, intent(in)  :: in
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), ray, pt, lp
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rt_t%ray_vol_map3', itimer=itimer)
  cell => packet%coords
  lp  = size(out,2)
  do ray = 1, packet%n_rays                                             ! loop over rays
    do pt = 1, packet%n_points                                          ! loop over points
      p = cell(:,pt,ray)                                                ! get the right cell
      out(ray,pt+1) = in(p(1),p(2),p(3))
    end do
    p  = packet%rev_coords(:,ray)
    out(ray,1) = in(p(1),p(2),p(3))                                     ! Get the first point
    p = packet%fwd_coords(:,ray)
    out(ray,lp) = in(p(1),p(2),p(3))                                    ! Get the last point
  end do
  call trace%end (itimer)
END SUBROUTINE ray_vol_map3

SUBROUTINE nbor_relations (self, another, needed, needs_me, download)
  class(rt_t):: self
  class(task_t), pointer:: another
  logical:: needed, needs_me, download
  !.............................................................................
  class(rt_t), pointer:: rt_another
  class(task_t), pointer:: other
  !-----------------------------------------------------------------------------
  call trace%begin ('rt_t%nbor_relations')
  !
  ! ----------------------------------------------------------------------------
  ! First check the main relations, based on location and level of refinement
  ! ----------------------------------------------------------------------------
  call self%task_t%nbor_relations(another, needed, needs_me, download)
  ! ----------------------------------------------------------------------------
  ! Now, if another is RT, add it appropriately, if it is MHD, it has to overlap
  ! and have different kind of relations
  ! ----------------------------------------------------------------------------
  select type (another)
  class is (rt_t)
  rt_another => another
  ! ----------------------------------------------------------------------------
  ! for no chains solver nothing has to be done here - all logicals are true, 
  ! unless vertical_propagation = .true. or some other reason.
  ! ----------------------------------------------------------------------------
  if(self%vertical_propagation) then
    if (any(self%position(1:2) /= rt_another%position(1:2))) &
      needed = .false.
  end if
  if (self%i_omega /= rt_another%i_omega) download = .false.
  !
  class default
  other => another
  ! ----------------------------------------------------------------------------
  ! Check if this is the "master" MHD patch. It might not be associated yet, i.e.
  ! after moving to another MPI rank or changing level of refinement. in that 
  ! case, associate it.
  ! ----------------------------------------------------------------------------
  if (self%overlaps(other)) then
    download = .false.  ! don't download 
    if (.not.(associated(self%patch))) self%patch => gpatch%cast2gpatch(other)
  end if
  end select

  !
  call trace%end()
END SUBROUTINE nbor_relations

END MODULE rt_mod
