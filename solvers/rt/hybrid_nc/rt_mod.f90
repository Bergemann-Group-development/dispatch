MODULE rt_mod
  USE bits_mod
  USE eos_mod
  USE io_mod
  USE link_mod
  USE kinds_mod
  USE math_mod
  USE mesh_mod
  USE ray_geometry_v3_mod
  USE rt_components_mod
  USE rt_boundaries_mod
  USE rt_solver_mod
  USE patch_mod
  USE task_mod
  USE trace_mod
  USE download_mod
  USE data_hub_mod
  USE scaling_mod
  USE units_mod
  USE omp_timer_mod
  USE planet_mod
  implicit none
  private
!=======================================================================
!> rt_t is the main task.
!> every task has a set of chain links, belonging to this patch. Every
!> item in the set is then connected to appropriate upstream/downstream
!> nbors. Also, any patch might end up being the chain head for some
!> bundle. In that case, it nudges it's own link.
!=======================================================================
  type, extends(patch_t), public :: rt_t
    class(chain_link_t), dimension(:), pointer :: my_links => null()     ! list of link, which belongs to me
    logical, dimension(:), pointer :: updated => null()                 ! list of flags, indicates if particular direction has been updated
    type(link_t), pointer             :: mhd_link                       ! link to my mhd_patch
    type(rt_outer_bndry_t)            :: rt_boundaries                  ! RT outer boundaries
    type(data_hub_t) :: data_hub
    real(kind=8)                      :: upd_time                       ! time when the last internal update was done
    integer(kind=8)                   :: n_packet                       ! number of bundles
    real(kind=KindRTVar), dimension(:), pointer       :: w_lambda                       ! wavelength weights
    real(kind=KindRTVar), dimension(:), pointer       :: w_bundle                       ! space angle weights
    real(kind=KindRTVar)                              :: cr_heating
    real(kind=KindScalarVar), dimension(:,:,:),   pointer :: d     => null()
    real(kind=KindRTVar), dimension(:,:,:),   pointer :: e     => null()
    real(kind=KindRTVar), dimension(:,:,:,:), pointer :: rk    => null()
    real(kind=KindRTVar), dimension(:,:,:,:), pointer :: src   => null()
    real(kind=KindRTVar), dimension(:,:,:),   pointer :: jj    => null()
    procedure(optical_depth_staggered), pointer :: optical_depth => null()
    character(len=16)                 :: solver ='integral'             ! choice of solver
    procedure(solve_integral), pointer:: solve => null()
    logical :: do_output      = .false.
    logical :: data_ready     = .false.
    logical :: do_chains      = .false.
    logical :: dry = .false.
    logical :: do_rt_sph_bc  = .true.
    real    :: epsilon=1.0
  contains
    procedure :: absolute_point_position
    procedure :: boundary_conditions
    procedure :: cell2pos
    procedure :: dnload
    procedure :: execute_no_chain
    procedure :: finish
    procedure :: get_ds
    procedure :: get_q
    procedure :: has_OB
    procedure :: has_IB
    procedure :: init => init_rt
    procedure :: init_chain
    procedure :: init_no_chain
    procedure :: map_bndry
    procedure :: nudge
    procedure :: optical_depth_centered
    procedure :: optical_depth_staggered
    procedure :: optical_depth_staggered_weighted
    procedure :: output => output_legacy_rt
    procedure :: pos2cell
    procedure :: pos2cell_interp
    procedure :: ray_vol_map3
    procedure :: ray_vol_map4
    procedure :: link_connects_to
    procedure :: solve_integral
    procedure :: solve_feautrier
    procedure :: update => update_rt
    procedure :: update_q
  end type
  logical:: debug=.false., detailed_timer=.false.
CONTAINS

!=======================================================================
!> The rt_patch initialization subroutine. The patch is responsible for
!> all the directions and their bundles. Al this is stored in chain links.
!> FLOW:
!> 1) read the namelist
!> 2) initialize the patch, set bits
!> 3) set
!=======================================================================
SUBROUTINE init_rt (self)
  class (rt_t) :: self
  !.....................................................................
  logical, save :: first_time= .true.
  logical, save :: do_output  = .false.
  logical, save :: use_data_hub = .true.
  logical, save :: periodic(3) = .false.
  logical, save :: do_rt_sph_bc = .false.
  real(8), save :: dtime=0.01_8, grace=0.00_8, init_time=0d0
  real(kind=KindRTVar), save :: cr_heating=0.0
  character(len=16) :: solver = 'integral'
  character(len=16) :: dtau = 'weighted'
  integer :: ignore, nd, i
  integer, save:: itimer = 0
  namelist /rt_params/ init_time, dtime, grace, solver, do_output, &
    cr_heating, use_data_hub, detailed_timer, debug, periodic, dtau, &
    do_rt_sph_bc
  ! --------------------------------------------------------------------

  call trace_begin ('rt_t%init', itimer=itimer)
  call io%lock%set
  if (first_time) then
    grace    = self%mhd_link%task%grace
    periodic = self%mhd_link%task%periodic
    rewind (io%input); read (io%input, rt_params, iostat=ignore)
    if (io%master) write (*, rt_params)
    call rt_solver%init (detailed_timer)
    first_time=.false.
  end if
  call io%lock%unset
  self%do_output = do_output
  if (do_output) error stop "rt_mod::init:: do_output=true, but IO for RT is not rewritten for the master branch"
  if (ray_geometry%head%ray_list%directionality(1) == 1) then
    nd = 2
  else
    nd = 1
  end if
  self%nv = ray_geometry%head%n_ray_dirs*nd*eos%n_lambda                ! memory slot stores all heating rates per direction per freq bin
  self%nw = 1                                                           ! no time deriv
  self%ng = 1
  self%periodic = periodic
  call self%patch_t%init                                                ! initialize patch
  self%use_data_hub = use_data_hub
  self%mem = 0.0
  call self%set (bits%static+bits%no_download)                          ! special bits
  self%dtime = dtime                                                    ! 0.1 second default
  self%time = init_time                                                 ! 0.0 default
  self%dt_fixed = dtime                                                 ! timestep is fixed for radiative transfer
  self%grace = grace                                                    ! grace on by default, but very small, so nboring rt_patches would nod do anything stupid
  self%kind = 'rt_solver'                                               ! set the patch kind for output post-processing purposes
  self%do_chains = ray_geometry%head%do_chains                          ! are there chains in the system, or patches are effectively isolated
  if (self%do_chains) then
    self%n_packet = ray_geometry%head%total_n_packets * nd              ! number of packets
  else
    self%n_packet = ray_geometry%head%total_n_packets
  end if
  self%cr_heating = cr_heating                                          ! store cr heating, if it is present
  do i=1,3
    self%mesh(i)%h(:)=0
  end do
  !> -------------------------------------------------------------------
  !> set up rt-specific parameters : weights for frequency bins and
  !> space-angles. Also, select the appropriate solver.
  !> -------------------------------------------------------------------
  allocate (self%w_lambda(self%nv))                                     ! frequency bin weights
  allocate (self%w_bundle(self%nv))                                     ! space angle weights
  self%w_bundle = math%pi4/(ray_geometry%head%n_ray_dirs * nd)               ! space angle weights FIXME: could/should be different weights for different directions

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

  allocate(self%my_links(self%n_packet))                                ! allocate all the links this patch contains
  if (self%do_chains) then
    allocate(self%updated(self%n_packet))                               ! flags, indicating update status of all the bundles
    self%updated(:) = .false.                                           ! none are updated initially
    call self%init_chain(nd)
  else
    call self%init_no_chain
  end if
  allocate (self%jj(self%gn(1), self%gn(2), self%gn(3)))                ! mean intensity
  self%jj = 0.0
  call trace%end(itimer)
END SUBROUTINE init_rt

!=======================================================================
!> Set up all the links, belonging to this patch. Every link represents
!> a separate bundle. One rt patch manages all the directions.
!=======================================================================
SUBROUTINE init_chain(self, nd)
  class(rt_t) :: self
  integer :: nd
  ! ....................................................................
  class (ray_dirs_t),   pointer :: ray_dir
  class (ray_bndl_t),   pointer :: ray_bnd
  class (packet_t),     pointer :: packet
  class (chain_link_t), pointer :: this_link
  integer :: a_id, ivl, dd, ibnd, ipck
  ! --------------------------------------------------------------------
  a_id = 1                                                              ! link personal ID in the link array
  ivl = 1
  do dd = 1,nd                                                          ! separate links for forward/reverse directions
    ray_dir => ray_geometry%head%ray_list%head                            ! pick a ray direction
    do while(associated(ray_dir))                                         ! loop over all ray directions
      do ibnd =1, ray_dir%n_bundles                                       ! loop over ray bundles
        ray_bnd => ray_dir%bundle(ibnd)
        do ipck = 1, ray_bnd%n_packets
          packet => ray_bnd%packet(ipck)
          this_link => self%my_links(a_id)
          this_link%link => self%link                                   ! links should know, what rt patch hosts them
          this_link%direction = dd                                      ! forward/reverse direction links
          this_link%id = a_id                                           ! link personal ID in the link array
          this_link%dir_id = ray_bnd%dir_id
          this_link%packet => packet                                     ! link has to know, which packet it corresponds to
          this_link%ivl = ivl
          this_link%ivu = ivl + eos%n_lambda-1
          allocate(this_link%pdata(packet%n_rays))                      ! allocate Q_in data for the link. Array consists of all the rays
          call self%get_ds(packet,this_link%ds)
          a_id = a_id + 1                                               ! next ID
        end do
      end do
      self%w_lambda(ivl:ivl+eos%n_lambda-1) = eos%w_lambda
      ray_dir => ray_dir%next                                           ! next direction
      ivl = ivl + eos%n_lambda
    end do
  end do
  nullify(ray_dir, ray_bnd, packet)
END SUBROUTINE init_chain

!=======================================================================
!> Set up all the links, belonging to this patch. Every link represents
!> a separate bundle. One rt patch manages all the directions.
!=======================================================================
SUBROUTINE init_no_chain(self)
  class(rt_t) :: self
  ! ....................................................................
  class (ray_dirs_t),   pointer :: ray_dir
  class (ray_bndl_t),   pointer :: ray_bnd
  class (packet_t),     pointer :: packet
  class (chain_link_t), pointer :: this_link
  integer :: a_id, ivl, ivl2, ray, pt, ipk
  ! --------------------------------------------------------------------
  a_id = 1                                                              ! link personal ID in the link array
  ivl = 1
  ivl2 = ray_geometry%head%n_ray_dirs*eos%n_lambda + 1
  ray_dir => ray_geometry%head%ray_list%head                            ! pick a ray direction
  do while(associated(ray_dir))                                         ! loop over all ray directions
    ray_bnd => ray_dir%bundle(1)
    do ipk = 1, ray_bnd%n_packets
      packet => ray_bnd%packet(ipk)
      this_link => self%my_links(a_id)
      this_link%link => self%link                                       ! links should know, what rt patch hosts them
      this_link%direction = 3                                           ! bi-directional
      this_link%id = a_id                                               ! link personal ID in the link array
      this_link%dir_id = ray_bnd%dir_id
      this_link%packet => packet                                        ! link has to know, which packet it corresponds to
      this_link%ivl = ivl
      this_link%ivu = ivl + eos%n_lambda-1
      this_link%ivl2 = ivl2
      this_link%ivu2 = ivl2 + eos%n_lambda-1
      this_link%master = .true.
      allocate(this_link%pdata(packet%n_rays*2))                        ! allocate Q_in data for the link. Array consists of all the rays
      call self%get_ds(packet,this_link%ds)
      call self%map_bndry(this_link)
      a_id = a_id + 1                                                   ! next ID
    end do
    self%w_lambda(ivl:ivl+eos%n_lambda-1)   = eos%w_lambda
    self%w_lambda(ivl2:ivl2+eos%n_lambda-1) = eos%w_lambda
    ray_dir => ray_dir%next                                             ! next direction
    ivl = ivl + eos%n_lambda
    ivl2 = ivl2 + eos%n_lambda
  end do
  nullify(ray_dir, ray_bnd, packet)
END SUBROUTINE init_no_chain

!=======================================================================
!> Type-specific version of `dnload`. This is a *duplicate* of what's in
!> `gpatch_t` because `rt_t` extends `patch_t`.
!=======================================================================
SUBROUTINE dnload (self, only)
  class (rt_t) :: self
  integer, optional:: only
  !---------------------------------------------------------------------
  if (self%use_data_hub) then
    call self%data_hub%update (self%link)
  else
    call download%download_link (self%link)
  end if
END SUBROUTINE dnload

!=======================================================================
!> This subroutine maps the incoming heating rates to the outer boundary
!> conditions
!=======================================================================
SUBROUTINE map_bndry(self, link)
  class(rt_t) :: self
  class(chain_link_t) :: link
  !.....................................................................
  integer :: r, i, c(3), lv, uv, lt, ut
  integer, dimension(:,:), pointer:: up, dn
  !---------------------------------------------------------------------
  select case(link%direction)
  case(1)                                                               ! forward rays
    lv = link%ivl; uv = link%ivu                                        ! memory slot indices
    up => link%packet%rev_coords                                        ! upstream coordinates
    do r = 1, link%packet%n_rays
      c = up(:,r)
      link%pdata(r)%in => self%mem(c(1),c(2),c(3),lv:uv,:,1)            ! map to guard cells
    end do
  case(2)                                                               ! reverse rays
    lv = link%ivl; uv = link%ivu                                        ! memory slot indices
    up => link%packet%fwd_coords                                        ! dnstream coordinates
    do r = 1, link%packet%n_rays
      c = up(:,r)
      link%pdata(r)%in => self%mem(c(1),c(2),c(3),lv:uv,:,1)            ! map to guard cells
    end do
  case(3)                                                               ! both directions
    lv = link%ivl ; uv = link%ivu                                       ! memory slot forward
    lt = link%ivl2; ut = link%ivu2                                      ! memory slot reverse
    up => link%packet%rev_coords                                        ! upstream coordinates
    dn => link%packet%fwd_coords                                        ! dnstream coordinates
    i = link%packet%n_rays
    do r = 1, link%packet%n_rays
      i = i + 1
      c = up(:,r)
      link%pdata(r)%in => self%mem(c(1),c(2),c(3),lv:uv,:,1)            ! map forward to mem
      c = dn(:,r)
      link%pdata(i)%in => self%mem(c(1),c(2),c(3),lt:ut,:,1)            ! map reverse to mem
    end do
   case default
     error stop "rt_t::map_bndry:: direction not understood"
   end select
END SUBROUTINE map_bndry

!=======================================================================
!> Solve the RT problem, by first obtaining the thermodynamics at the
!> designated RT time, and using them to compute the opacity and source
!> function values. Then nudge the chain heads.
!=======================================================================
SUBROUTINE update_rt (self)
  class (rt_t) :: self
  class(task_t),  pointer:: task
  class(patch_t), pointer:: mhd
  integer, save:: itimer=0
  integer :: ii, iy, l(3), u(3), i, iv
  real(kind=KindRTVar) :: g1
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: k
  !---------------------------------------------------------------------
  if (self%data_ready) return
  call trace%begin ('rt_t%update', itimer=itimer)
 ! if (io%verbose > 1) print *, 'RADIATIVE TRANSFER TIME', self%time
 ! if (self%time < 0d0) then
 !   self%time  = 0d0                                                    ! set current time, if it was negative (to do the first RT run before MHD updates)
 ! end if
  task => self%mhd_link%task                                            ! pick the rmhd task
  select type (task)
  class is (patch_t)
    mhd => task
  end select
  l = mhd%mesh%lo
  u = mhd%mesh%uo
  self%time = mhd%time
  self%dtime= mhd%dtime
  self%upd_time = self%time
  !> -------------------------------------------------------------------
  !> allocate all the scratch arrays
  !> -------------------------------------------------------------------
  allocate (self%d (self%gn(1), self%gn(2), self%gn(3)))                ! density scratch
  allocate (self%e (self%gn(1), self%gn(2), self%gn(3)))                ! entropy scratch
  allocate (self%rk  (self%gn(1), self%gn(2), self%gn(3), eos%n_lambda))! opacity scratch
  allocate (self%src (self%gn(1), self%gn(2), self%gn(3), eos%n_lambda))! source function scratch

  !> -------------------------------------------------------------------
  !> Take entropy and density from MHD patch
  !> -------------------------------------------------------------------
  self%d = mhd%mem(l(1):u(1),l(2):u(2),l(3):u(3),mhd%idx%d,mhd%it,1)
  if (mhd%kind(1:9)=='stagger2e') then
    self%e = mhd%mem(l(1):u(1),l(2):u(2),l(3):u(3),mhd%idx%e,mhd%it,1)
    self%e = self%e/self%d
  else if (mhd%kind(1:9)=='stagger2_') then
    associate(d => mhd%mem(l(1):u(1),l(2):u(2),l(3):u(3),mhd%idx%d,mhd%it,1), &
              s => mhd%mem(l(1):u(1),l(2):u(2),l(3):u(3),mhd%idx%s,mhd%it,1))
    g1 = mhd%gamma-1.0
    self%e = d**g1*exp(s/d*g1)/g1
    end associate
  else if (mhd%kind(1:7)=='ramses_') then
    associate(d => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%d           ,mhd%it,1), &
              e => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%e           ,mhd%it,1), &
              p => self%mem(l(1):u(1),l(2):u(2),l(3):u(3),self%idx%px:self%idx%pz,mhd%it,1))
    allocate(k(size(d,1),size(d,2),size(d,3)))
    k = 0.0
    do i=1,size(p,4)
      k = k+p(:,:,:,i)**2
    end do
    k = 0.5*k/d
    self%e = (e-k)/d
    deallocate(k)
    end associate
  else
    print *,trim(mhd%kind)
    error stop "rt_mod (1): MHD solver type not known"
  end if
  !> -------------------------------------------------------------------
  !> Get the source function and rk
  !> -------------------------------------------------------------------
  call eos%lookup (e=self%e, d=self%d, rk=self%rk, src=self%src)
  !> -------------------------------------------------------------------
  !> Combine scattering mean <Q> with thermal source part
  !> -------------------------------------------------------------------
  do iv = 1, eos%n_lambda
    self%src(:,:,:,iv) = (1.-self%epsilon)*self%jj + self%epsilon*self%src(:,:,:,iv)
  end do

  if (self%do_rt_sph_bc) then
    call rt_boundaries%spherical_jagged(self%src,self%rk,self%mesh)
  end if
   !   print *, 'minmax rk, src', minval(self%rk), maxval(self%rk), minval(self%src), maxval(self%src)
   ! print *,'along x'
   ! print *, self%src(:,self%n(2)/2,self%n(3)/2,:)
   ! print *, 'while self%e'
   ! print *, self%e(:,self%n(2)/2,self%n(3)/2)
   ! print *, 'but s'
   ! print *, mhd%mem(:,self%n(2)/2,self%n(3)/2,mhd%i_s,mhd%it,1)
  self%data_ready =.true.                                               ! flag, that internal data is ready
  if(io%verbose >1) then
    print *, self%id, 'rk', minval(self%rk), maxval(self%rk), 'cgs', &
      minval(self%rk(:,:,:,1)/self%d/scaling%m*scaling%l**2), &
      maxval(self%rk(:,:,:,1)/self%d/scaling%m*scaling%l**2)
    print *, self%id, 'src', minval(self%src), maxval(self%src)
  end if
  if (self%do_chains) then
    call self%rt_boundaries%condition (self%mem(:,:,:,:,self%it,1), self%src, self%rk)
    call self%nudge                                                     ! nudge the chain heads
  else
    call self%dnload
    call self%rt_boundaries%condition (self%mem(:,:,:,:,self%it,1), self%src, self%rk)
    call self%execute_no_chain
  end if
  !> -------------------------------------------------------------------
  !> Compute the updated QQ = <Q>
  !> -------------------------------------------------------------------
  self%jj = 0.0
  do iv = 1, self%nv
    self%jj = self%jj + &
            + self%mem(:,:,:,iv,self%it,1) &
            * self%w_bundle(iv) * self%w_lambda(iv)/math%pi4
  end do
  self%jj = self%jj + self%src(:,:,:,1)
  !
  call trace%end (itimer)
END SUBROUTINE update_rt

!=======================================================================
!> Nudge the chain heads. If all patches in the chain become ready, then
!> solve for the net heating rate, storing it in a slot that rotates
!> with time for each patch in the chain.
!=======================================================================
SUBROUTINE nudge (self)
  class(rt_t) :: self
  !.....................................................................
  class(chain_link_t), pointer :: link
  integer :: ip
  !---------------------------------------------------------------------
  do ip = 1, self%n_packet
    link => self%my_links(ip)
    !$omp atomic update
    link%head%n_ready = link%head%n_ready+1                             ! nudge the master link
  if ( link%head%nudged(link%chain_nm)) then                            ! mark myself as done in the ready array
    error stop "Whoops: nudging the nudged"
  end if
    !$omp atomic write
    link%head%nudged(link%chain_nm) = .true.                            ! mark myself as done in the ready array
    if (all(link%head%nudged)) call execute_chain(link)                  ! execute the chain, when all links are ready
  end do
END SUBROUTINE nudge

!=======================================================================
!> Loop through all links in the chain and solve the actual rt
!=======================================================================
SUBROUTINE execute_chain (self)
  class(chain_link_t), pointer :: self
  !.....................................................................
  class(rt_t),          pointer :: rt_patch
  class(task_t),        pointer :: task
  class(chain_link_t),  pointer :: link, prev
  real(kind=KindRTVar), dimension(:,:  ), pointer :: q0
  real(kind=KindRTVar), dimension(:,:,:), pointer :: q
  integer :: ray, it
  integer, save:: itimer=0
  logical :: first
  !---------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('rmhd_t%execute_chain', itimer=itimer)
  link => self%head                                                     ! take the master link
  nullify(prev)                                                         ! just in case - there has to be no "prev" for the master link
  task => link%link%task                                                ! get the patch of current chain link
  select type (task)
  class is (rt_t)
  rt_patch => task
  end select
  allocate(q0 (link%packet%n_rays,eos%n_lambda))
  first = .true.
  do while (associated(link))
    allocate(q (link%packet%n_rays,eos%n_lambda,link%packet%n_points+2))
    if (associated(prev)) then
      task => prev%link%task                                            ! get the patch of upstream chain link
      select type (task)
      class is (rt_t)
      rt_patch => task
      end select
      if  (rt_patch%upd_time < rt_patch%time) then
        it = rt_patch%iit(rt_patch%nt-2)
      else
        it = rt_patch%it
      end if
      do ray=1, link%packet%n_rays
        q0(ray,:) = link%pdata(ray)%in(:,it)
      end do
      task => link%link%task                                            ! get the patch of current chain link
      select type (task)
      class is (rt_t)
      rt_patch => task
      end select
    else
      task => link%link%task                                            ! get the patch of current chain link
      select type (task)
      class is (rt_t)
      rt_patch => task
      end select
      if  (rt_patch%upd_time < rt_patch%time) then
        it = rt_patch%iit(rt_patch%nt-2)
      else
        it = rt_patch%it
      end if
      do ray=1, link%packet%n_rays
        q0(ray,:) = link%pdata(ray)%in(:,it)
      end do
    end if
    !
    select case(self%direction)
      case(1)                                                           ! forward rays
        call rt_patch%solve(link,q_in_fw=q0,q_out_fw=q)                 ! do the radiative transfer for this particular link
        call rt_patch%update_q(link, q_fw=q)                            ! update the memory slot
      case(2)                                                           ! reverse rays
        call rt_patch%solve(link,q_in_re=q0,q_out_re=q)                 ! do the radiative transfer for this particular link
        call rt_patch%update_q(link, q_re=q)                            ! update the memory slot
      case default
      error stop "rt_t::execute_chain:: direction not understood"
    end select
    rt_patch%updated(link%id) = .true.                                  ! mark, that this link has been updated
    if(all(rt_patch%updated)) call rt_patch%finish                      ! if all links are done, this rt_patch is complete
    deallocate(q)
    prev => link
    link => link%next
  end do
  deallocate(q0)                                                        ! deallocate the scratch
  self%head%nudged = .false.                                            ! reset the nudge logicals
  self%head%n_ready = 0                                                 ! no links are ready in this chain
  nullify(rt_patch, task, prev)
  if (detailed_timer) call trace%end(itimer)
END SUBROUTINE execute_chain

!=======================================================================
!> Loop through all links in the chain and solve the actual rt
!=======================================================================
SUBROUTINE execute_no_chain (self)
  class(rt_t) :: self
  ! ....................................................................
  class(chain_link_t),  pointer :: link
  class(packet_t),      pointer :: packet
  real(kind=KindRTVar), dimension(:,:  ), pointer :: q_fw0, q_re0
  real(kind=KindRTVar), dimension(:,:,:), pointer :: q_fw , q_re
  integer :: ip, ipck, nrays, npt, ir,tr
  integer, save:: itimer=0
  ! --------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('rt_t%execute_no_chain', itimer=itimer)
  do ip = 1, self%n_packet
    link => self%my_links(ip)
    packet => link%packet
    nrays  =  packet%n_rays                                             ! number of rays in this packet
    npt    =  packet%n_points+2
    allocate(q_fw0 (nrays,eos%n_lambda))
    allocate(q_re0 (nrays,eos%n_lambda))
    allocate(q_fw (nrays,eos%n_lambda,npt))
    allocate(q_re (nrays,eos%n_lambda,npt))
    do ir = 1, nrays
      q_fw0(ir,:) = link%pdata(ir)%in(:,self%it)
    end do
    tr=1
    do ir = nrays+1, nrays*2
      q_re0(tr,:) = link%pdata(ir)%in(:,self%it)
      tr=tr+1
    end do
    call self%solve(link,q_in_fw=q_fw0, q_in_re=q_re0, &
                         q_out_fw=q_fw, q_out_re=q_re)          ! do the radiative transfer for this particular link
   ! if (self%id == 147) then
   !   print *, 'fw', q_fw0(12,:), 're', q_re0(12,:)
   !   print *, 'ray 12 fw',q_fw(12,:,:)
   !   print *, 'ray 12 re', q_re(12,:,:)
   ! end if
    call self%update_q(link, q_fw=q_fw, q_re=q_re)
    deallocate (q_fw, q_re, q_fw0, q_re0)
  end do
  if (detailed_timer) call trace%end(itimer)
END SUBROUTINE execute_no_chain

!=======================================================================
!> Check if the rays in self match up in the downstream with the rays in other
!> INPUT: self and nbor are rt_patches. mine is my link i am testing
!>   against nbor's link other. logical connects - true if chain can be
!>   made, logical is_upstream - if the link is in upstream from me, but
!>   is on a different refinement level, it could still potentially
!>   give me the intensities i need.
!> RESULT: two logicals: if can be connected and if one is upstream of
!> another, but on different refinement level.
!> SUBROUTINES CALLED FROM WITHIN:
!> rt_mod_t::rt_t%cell2pos
!> rt_mod_t::rt_t%pos2cell
!> FLOW:
!> 1) both logicals are set to false by default.
!> 2) both links are checked whether they correspond to same directions and
!>   direction IDs (spatial direction).
!> 3) rt tasks, corresponding to these links, are checked for different
!>   refinement level. If that is the case, appropriate flag is set.
!> 4) otherwise the connections are tested. Only one ray in a packet is
!>   enough to see, whether connection is good - ray_geometry module makes
!>   sure all packets are fully interconnected.
!> 5) the last point is taken from other link and the upstream position
!>   is taken from me. if the points coincide via cell2pos and pos2cell,
!>   there is a connection.
!=======================================================================
SUBROUTINE link_connects_to (self, nbor, mine, other, connects, is_upstream)
  class (rt_t) :: self
  class (rt_t), pointer :: nbor
  class (chain_link_t) :: mine
  class (chain_link_t) :: other
  logical :: connects, is_upstream
  !.....................................................................
  integer :: ray, upc(3), ir, oc(3), pt, i
  !---------------------------------------------------------------------
  connects = .false.
  is_upstream = .false.

  if (mine%dir_ID /= other%dir_ID)       return                         ! different ray spatial direction
  if (mine%direction /= other%direction) return                         ! different direction

  if (self%level /= nbor%level) then
    if (abs(self%level-nbor%level) < 2) then
    !  select case(mine%direction)
    !  case(1)
    !    if (mine%bundle%lwall(1) /= 0) then
    !      a = 2; b = 3; c=1
    !    else if (mine%bundle%lwall(2) /= 0) then
    !      a = 3; b = 1; c=2
    !    else if (mine%bundle%lwall(3) /= 0) then
    !      a = 1; b = 2; c=3
    !    end if
    !  case(2)
    !    if (mine%bundle%uwall(1) /= 0) then
    !      a = 2; b = 3
    !    else if (mine%bundle%uwall(2) /= 0) then
    !      a = 3; b = 1
    !    else if (mine%bundle%uwall(3) /= 0) then
    !      a = 1; b = 2
    !    end if
    !  case default
    !    error stop "rt_t::link_connects_to:: direction not understood(1)"
    !  end select

     ! mrn = 1
     ! do ip = 1, mine%bundle%n_packets
     !   do ray = 1, mine%bundle%packet(ip)%n_rays
     !     first_pt = mine%get_upstream_coord(ip,ray)                    ! take the first point
     !     ppos = self%myposition(first_pt)                              ! determine the absolute point position
     !     cp = nbor%cell2pos(ppos)                                      ! locate where the point could be in the trial patch
     !     do i=1,3                                                      ! is the ray originates in this patch?
     !       within(i) = .false.
     !       if ((cell(i)>=nbor%li(i)).and.(cell(i)<=nbor%ui(i))) within(i) = .true.
     !     end do
     !     if (all(within)) then
     !       ! determine for (a) axis


     error stop "rt_mod: error 10: different levels might not be ready yet"
     end if
     !   mrn = 1                                                         ! my global ray number in packet list
     !   do ip = 1, mine%bundle%n_packets                                ! loop over my packets
     !     do ray = 1, mine%bundle%packet(ip)%n_rays                     ! loop over my rays in the packet
     !       first_pt = mine%bundle%packet(ip)%rev_coords(:,ray)         ! take the first point
     !       do pl = 0,1
     !       call nbor%pos2cell_interp(self%myposition(first_pt),pl,cell,iw)
     !
     !       do i=1,3
     !         within(i) = .false.
     !         if ((cell(i)>=nbor%li(i)).and.(cell(i)<=nbor%ui(i))) within(i) = .true.
     !       end do
     !       if (all(within)) then
     !         call nbor%cell2shell(cell, s_coord, face)
     !         select case(face)                                  ! half of these go to other mine%direction case
     !         case(-1)
     !           wall => nbor%I_yz(:,:,:,2)
     !         case(-2)
     !           wall => nbor%I_xz(:,:,:,2)
     !         case(-3)
     !           wall => nbor%I_xy(:,:,:,2)
     !         end select
     !         print *, s_coord, size(wall,3), cell, pl, iw, first_pt, 'face',face, 'pos', self%myposition(first_pt)
     !         floop: do ics=1, size(wall,3)
     !           if (wall(s_coord(1), s_coord(2),ics)%dir_ID == mine%bundle%dir_ID) then
     !             cs => wall(s_coord(1), s_coord(2),ics)
     !             exit floop
     !           end if
     !         end do floop
     !         call mine%add_upstream(mrn, cs, iw)
     !         is_upstream = .true.
     !       end if
     !       end do
     !     end do
     !   end do
    !  case(2)
    !    mrn = 1                                                         ! my global ray number in packet list
    !    do ip = 1, mine%bundle%n_packets                                ! loop over my packets
    !      do ray = 1, mine%bundle%packet(ip)%n_rays                     ! loop over my rays in the packet
    !        first_pt = mine%bundle%packet(ip)%fwd_coords(:,ray)         ! take the first point
    !        do pl = 0,1
    !        call nbor%pos2cell_interp(self%myposition(first_pt),pl,cell,iw)
    !        within = .false.
    !        do i=1,3
    !          if ((cell(i)>=nbor%li(i)).and.(cell(i)<=nbor%ui(i))) within(i) = .true.
    !        end do
    !        if (all(within)) then
    !          call nbor%cell2shell(cell, s_coord, face)
    !          select case(face)                                  ! half of these go to other mine%direction case
    !          case(1)
    !            wall => nbor%I_yz(:,:,:,2)
    !          case(2)
    !            wall => nbor%I_xz(:,:,:,2)
    !          case(3)
    !            wall => nbor%I_xy(:,:,:,2)
    !          end select
    !          print *, s_coord, size(wall,3), cell, pl
    !          find_loop: do ics=1, size(wall,3)
    !            if (wall(s_coord(1), s_coord(2),ics)%dir_ID == mine%bundle%dir_ID) then
    !              cs => wall(s_coord(1), s_coord(2),ics)
    !              exit find_loop
    !            end if
    !          end do find_loop
    !          call mine%add_upstream(mrn, cs, iw)
    !          is_upstream = .true.
    !        end if
    !        end do
    !      end do
    !    end do
    !  end select
    !end if
  else
    ray = 1                                                             ! testing only one ray in the bundle
    upc = mine%get_upstream_coord_same(ray)
    select case(mine%direction)
    case(1)                                                             ! look at last points in other rays
      pt = other%packet%n_points
    case(2)                                                             ! look at first points in other rays
      pt = 1
    case default
      error stop "rt_t::link_connects_to:: direction not understood(2)"
    end select
    do ir = 1, other%packet%n_rays
      oc = other%packet%coords(:,pt,ir)
      if (all(oc==upc)) then
        connects = .true.
        exit
      end if
    end do
  end if
END SUBROUTINE link_connects_to

!=======================================================================
!> Find the absolute position of index 'ray' in ray packet 'rays'.
!> If the optional bn parameter is present, the answer depends on
!> the direction, otherwise 'pt' chooses the point along the ray.
!=======================================================================
FUNCTION absolute_point_position (self, link, ray, pt, bn) result (cell_pos)
  class (rt_t) :: self
  class (chain_link_t), pointer :: link
  class (packet_t), pointer :: packet
  integer, optional :: bn
  integer :: ray, pt, cp(3)
  real(kind=8) :: cell_pos(3)
!-----------------------------------------------------------------------
  packet => link%packet
  if (present(bn)) then
    select case(link%direction)
    case(1)
      cp = packet%rev_coords(:,ray)
    case(2)
      cp = packet%fwd_coords(:,ray)
    end select
  else
    cp = packet%coords(:,pt,ray)
  end if
  cell_pos = self%myposition(cp)
END FUNCTION absolute_point_position

LOGICAL FUNCTION has_OB (self, link)
  class (rt_t) :: self
  class (chain_link_t) :: link
  integer :: i(3), t
  real    :: d(3)
  real(8) :: p(3)
  !.....................................................................
  has_OB = .false.
  d = link%packet%dstep
  select case (link%direction)
  case(1)
    i = link%packet%rev_coords(:,1)
    p = self%myposition(i)
    do t = 1,3
      if (d(t)>0) then
        if (p(t)<rt_boundaries%llc(t)) has_OB = .true.
      else
        if (p(t)>rt_boundaries%urc(t)) has_OB = .true.
      end if
    end do
  case(2)
    i = link%packet%fwd_coords(:,1)
    p = self%myposition(i)
    do t = 1,3
      if (d(t)>0) then
        if (p(t)>rt_boundaries%urc(t)) has_OB = .true.
      else
        if (p(t)<rt_boundaries%llc(t)) has_OB = .true.
      end if
    end do
  case default
    error stop "rt_t::has_OB:: direction not understood"
  end select
END FUNCTION has_OB

LOGICAL FUNCTION has_IB (self, link)
  class (rt_t) :: self
  class (chain_link_t) :: link
  integer :: i(3), d(3),t
  real(8) :: p(3)
  !.....................................................................
  has_IB = .false.
END FUNCTION has_IB

SUBROUTINE finish(self)
  class(rt_t) :: self
  integer :: ii, iv, wb
  real, dimension(:,:,:), pointer :: qq
  class(task_t),  pointer:: task
  class(patch_t), pointer:: mhd
  !---------------------------------------------------------------------
  call trace%begin ('rt_t%finish')
  if (.not.self%data_ready) then
    error stop "rt_t%finish called unexpectedly"
  end if
  self%data_ready = .false.
  if (self%do_chains)  self%updated(:)=.false.
  task => self%mhd_link%task                                            ! pick the rmhd task
  select type (task)
  class is (patch_t)
    mhd => task
  end select
  if (.not.self%dry.and.self%do_output) then
    !call self%output
    error stop "rt_mod::finish:: output called, but IO for RT needs to be rewritten"
  end if
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)            ! predict
  deallocate (self%d, self%e, self%rk, self%src)                        ! deallocate scratch
  call trace%end()
END SUBROUTINE finish

SUBROUTINE solve_integral (self, link, q_in_fw, q_in_re, q_out_fw, q_out_re)
  class(rt_t) :: self
  class(chain_link_t),    pointer :: link
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_fw
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_re
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_fw
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_re
  !.....................................................................
  class(packet_t),        pointer :: packet
  real(kind=KindRTVar), dimension(:,:,:), pointer :: s, rk, dtau, exp0
  integer                         :: nrays, npt, nwb, ir
  integer, save                   :: itimer=0
  ! ------------- testing purposes
    class(task_t),  pointer:: task
  class(patch_t), pointer:: mhd
  logical :: tst(3)
  real:: fnd(3)
  ! --------------------------------

  ! --------------------------------------------------------------------
  if (detailed_timer) call trace_begin ('rt_t%solve_integral', itimer=itimer)
  packet => link%packet
  nwb   =  eos%n_lambda                                                 ! number of wl bins
  nrays =  packet%n_rays                                                ! number of rays in this packet
  !> -------------------------------------------------------------------
  !> RT solvers need at least 3 points to work with, as they use
  !> derivatives. Rays, shorter, than 3 points, have to be extended
  !> to patch boundaries to have accurate solutions.
  !> -------------------------------------------------------------------
  npt = packet%n_points + 2
  !> -----------------------------------------------------------------
  !> allocate all the scratches
  !> -----------------------------------------------------------------
  allocate(s     (nrays, nwb, npt))                                     ! source function
  allocate(rk    (nrays, nwb, npt))                                     ! opacity
  allocate(dtau  (nrays, nwb, npt))                                     ! optical depth between the points
  allocate(exp0  (nrays, nwb, npt))                                     ! first exponential factor
  dtau = 0.0
  !> -------------------------------------------------------------------
  !> map the quantities onto rays
  !> -------------------------------------------------------------------
  call self%ray_vol_map4(packet, self%src, s)                           ! map source function to rays
  call self%ray_vol_map4(packet, self%rk, rk)                           ! map opacities to rays
  call self%optical_depth(link%ds, rk, dtau)
#if defined __INTEL_COMPILER
  exp0(:,:,2:) = 2**(-dtau(:,:,2:)*1.44269504)
#else
  exp0(:,:,2:) = exp(-dtau(:,:,2:))                                                     ! expopential factors
#endif
  if (self%do_chains) then
    if (link%direction ==1) then
      call rt_solver%integral_forward (dtau, exp0, s, q_in_fw, q_out_fw)! solve forward
    else
      call rt_solver%integral_reverse (dtau, exp0, s, q_in_re, q_out_re)! solve reverse
    end if
  else
    if (debug) then
      write (31) shape(dtau)
      write (31) dtau, exp0, s, q_in_fw
      print *,'integral_double   dtau:',self%id,minval(dtau),maxval(dtau),shape(dtau)
      print *,'integral_double      s:',minval(s),maxval(s)
      !call rt_solver%integral_double_ref(dtau, exp0, S, q_in_fw, q_in_re, q_out_fw, q_out_re)
      print *,'integral_double   q_fw:',minval(q_out_fw),maxval(q_out_fw)
      print *,'integral_double   q_rw:',minval(q_out_re),maxval(q_out_re)
      write (31) q_out_fw,q_out_re
      call rt_solver%integral_double(dtau, exp0, S, q_in_fw, q_in_re, q_out_fw, q_out_re)
      print *,'integral_double   q_fw:',minval(q_out_fw),maxval(q_out_fw)
      print *,'integral_double   q_rw:',minval(q_out_re),maxval(q_out_re)
      write (31) q_out_fw,q_out_re
    else
      call rt_solver%integral_double(dtau, exp0, S, q_in_fw, q_in_re, q_out_fw, q_out_re)
    end if
  end if
  deallocate (s, rk, dtau, exp0)
  if (detailed_timer) call trace_end (itimer)
END SUBROUTINE solve_integral

!=======================================================================
!> Same as solve_integral, but for the Feautrier formulation
!=======================================================================
SUBROUTINE solve_feautrier (self, link, q_in_fw, q_in_re, q_out_fw, q_out_re)
  class(rt_t) :: self
  class(chain_link_t),    pointer :: link
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_fw
  real(kind=KindRTVar), dimension(:,:),   pointer, optional :: q_in_re
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_fw
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_out_re
  !.....................................................................
  class(packet_t),        pointer :: packet
  integer                         :: nrays, npt, nwb, ir
  real(kind=KindRTVar), dimension(:,:,:), pointer :: s, rk, dtau
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  if (detailed_timer) call trace_begin ('rt_t%solve_feautrier', itimer=itimer)
  error stop "rt_t::solve_featurier:: feautrier solver not consistent"
  !> TODO: if feautrier solver is used, there should be no forward/reverse
  !> directionality in memory slots, or forward/reverse heating rates have
  !> to be split in rt_solver_mod
  packet => link%packet
  nwb   =  eos%n_lambda                                                 ! number of wl bins
  nrays =  packet%n_rays                                                ! number of rays in this packet
  !> -------------------------------------------------------------------
  !> RT solvers need at least 3 points to work with, as they use
  !> derivatives. Rays, shorter, than 3 points, have to be extended
  !> to patch boundaries to have accurate solutions.
  !> -------------------------------------------------------------------
  npt = packet%n_points + 2
  !> -----------------------------------------------------------------
  !> allocate all the scratches
  !> -----------------------------------------------------------------
  allocate(s     (nrays, nwb, npt))                                     ! source function
  allocate(rk    (nrays, nwb, npt))                                     ! opacity
  allocate(dtau  (nrays, nwb, npt))                                     ! optical depth between the points
  !> -------------------------------------------------------------------
  !> map the quantities onto rays
  !> -------------------------------------------------------------------
  call self%ray_vol_map4(packet, self%src, s)                           ! map source function to rays
  call self%ray_vol_map4(packet, self%rk, rk)                           ! map opacities to rays
  call self%optical_depth(link%ds, rk, dtau)
  call rt_solver%feautrier(dtau, S, q_out_fw, q_in_fw, q_in_re, .true.)
  if (detailed_timer) call trace_end (itimer)
END SUBROUTINE solve_feautrier

!=======================================================================
!> Map cell data onto rays
!=======================================================================
SUBROUTINE ray_vol_map4 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real(kind=KindRTVar), dimension(:,:,:)   , pointer :: out
  real(kind=KindRTVar), dimension(:,:,:,:) , pointer :: in
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), p2(3), ray, pt, nwb, wl, lp
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('rt_t%ray_vol_map4', itimer=itimer)
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
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE ray_vol_map4

!=======================================================================
!> Map cell data onto rays. In this case, it is wavelenght independent.
!=======================================================================
SUBROUTINE ray_vol_map3 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real(kind=KindRTVar), dimension(:,:)   , pointer, intent(out) :: out
  real(kind=KindRTVar), dimension(:,:,:) , pointer, intent(in)  :: in
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), ray, pt, lp
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('rt_t%ray_vol_map3', itimer=itimer)
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
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE ray_vol_map3

SUBROUTINE update_q (self, link, q_fw, q_re)
  class(rt_t) :: self
  class(chain_link_t), pointer :: link
  real(kind=KindRTVar), dimension(:,:,:), pointer, optional :: q_fw, q_re
  ! ....................................................................
  class(packet_t), pointer :: packet
  real, dimension(:,:,:), pointer :: q
  real, dimension(:,:,:,:), pointer :: qq, qq2
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), wl, ray, pt, ib, ivl, ivu, it, ivl2, ivu2, nwb
  integer, save:: itimer=0
  real(8):: wc
  !---------------------------------------------------------------------
  if (detailed_timer) call trace_begin ('rt_t%update_q', 4, itimer=itimer)
  ivl = link%ivl
  ivu = link%ivu
  nwb = size(q_fw,2)
  packet => link%packet
  if (self%do_chains) then
    if  (self%upd_time < self%time) then
      it = self%iit(self%nt-2)
    else
      it = self%it
    end if
    if (present(q_fw)) q => q_fw
    if (present(q_re)) q => q_re
    qq => self%mem(:,:,:,ivl:ivu,it,1)                                  ! point to the memory slot

    cell => packet%coords                                                 ! point to coordinates
    do ray = 1, packet%n_rays                                             ! loop over rays
      do pt = 1, packet%n_points                                          ! loop over points
        p = cell(:,pt,ray)                                                ! get the cell position
        do wl=1,nwb
          qq(p(1),p(2),p(3),wl) = q(ray,wl,pt+1)
        end do
      end do
    end do
  else
    it = self%it
    ivl2 = link%ivl2
    ivu2 = link%ivu2
    qq  => self%mem(:,:,:,ivl :ivu ,it,1)                               ! point to the forward direction memory slot
    qq2 => self%mem(:,:,:,ivl2:ivu2,it,1)                               ! point to the reverse direction memory slot
    cell => packet%coords                                               ! point to coordinates
    !wc = wallclock()
    do ray = 1, packet%n_rays                                           ! loop over rays
      do pt = 1, packet%n_points                                        ! loop over points
        p = cell(:,pt,ray)                                              ! get the cell position
        do wl=1,nwb
          qq (p(1),p(2),p(3),wl) = q_fw(ray,wl,pt+1)
          qq2(p(1),p(2),p(3),wl) = q_re(ray,wl,pt+1)
        end do
      end do
    end do
    !print *, (wallclock()-wc)/(packet%n_rays*packet%n_points)           ! takes 0.25 mus/pt on laptop
  end if
  if (detailed_timer) call trace_end (itimer)
END SUBROUTINE update_q

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
!> Get heating rates for given time in terms of the already averaged <Q>.
!> Scattering reduces the heat transfer in proportion to epsilon.
!=======================================================================
SUBROUTINE get_q (self, time, de)
  class (rt_t) :: self
  real(8):: time
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: de
  !.....................................................................
  integer :: iv
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rt_t%get_q', itimer=itimer)
  associate (l => self%li, u => self%ui)
  de = de + (self%jj (l(1):u(1),l(2):u(2),l(3):u(3)) &
           - self%src(l(1):u(1),l(2):u(2),l(3):u(3),1)) &
          * self%rk(l(1):u(1),l(2):u(2),l(3):u(3),1) * math%pi4
  end associate
  call trace%end (itimer)
END SUBROUTINE get_q

!=======================================================================
!> Return the absolute position of integer index triple 'cell'
!=======================================================================
FUNCTION cell2pos (self, cell) result(p)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  p = self%position - 0.5_8*self%size + (cell-self%li+1)*self%ds
END FUNCTION cell2pos

!=======================================================================
!> Convert absolute position to integer index triple 'cell'
!=======================================================================
FUNCTION pos2cell (self, p) result (cell)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  !cell = nint((p-self%position+0.5_8*self%size)/self%ds)+self%li-1

  cell = nint((p-self%position)/self%ds+self%offset)

END FUNCTION pos2cell

!=======================================================================
!> Convert absolute position to integer index triple 'cell'
!> currently simple spatial interpolation is used and should be sufficient.
!> if it is not accurate enough, it should be changed to something else.
!=======================================================================
SUBROUTINE pos2cell_interp (self, p, fc, cell, iw)
  class(rt_t)  :: self
  real(kind=8) :: p(3), pp(3)
  real :: iw
  integer:: fc, cell(3)
!-----------------------------------------------------------------------
  !cell = nint((p-self%position+0.5_8*self%size)/self%ds)+self%li-1
  pp = (p-self%position)/self%ds+self%offset
  if (fc==0) then
    cell = floor(pp)
    iw   = product(1. - (pp-cell))
  elseif (fc==1) then
    cell = ceiling(pp)
    iw   = product(pp-floor(pp))
  end if
END SUBROUTINE pos2cell_interp

!=======================================================================
!> Apply boundary conditions
!=======================================================================
SUBROUTINE boundary_conditions (self)
  class(rt_t):: self
  !---------------------------------------------------------------------
END SUBROUTINE boundary_conditions

SUBROUTINE get_ds(self, packet, ds)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real(kind=8) :: cs(3), ds
  !.....................................................................
  real(kind=8) :: dst(3), c(3)
  !---------------------------------------------------------------------
  dst = real(packet%dstep,kind=8)
  c   = self%mesh%d
  ds  = ((c(1)*dst(1))**2+(c(2)*dst(2))**2+(c(3)*dst(3))**2)**0.5

END SUBROUTINE get_ds

!=======================================================================
!> Write results to disk
!=======================================================================
SUBROUTINE output_legacy_rt (self, experiment_name)
  USE iso_fortran_env, only : real32
  class(rt_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3), idir, il, ind, nd
  real:: pt, qt
  real(kind=KindRTVar), dimension(:,:,:), pointer:: buf, b2
  real(8):: time
  integer, save:: iodir=-1, itimer=0
  !-----------------------------------------------------------------------------
  if (.not. io%do_output) return
  call trace_begin('rt_t%output_legacy', itimer=itimer)
  do while (self%time >= self%out_next)
  !$omp critical (out_cr)
  if (io%verbose>2) print *, self%id,'RT legacy output, TIME:', self%time
  self%max_files = self%max_files-1
  if (self%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  
  error stop "ERROR: rt_mod::output_legacy_rt:: IO for RT has to be redone"
  
  
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is true,
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (io%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (io%time_derivs>0) then
    no = min(2,self%nw)
  else
    no = 1
  end if
  if (io%out_time==0.0) then
    time = self%time
  else if (io%time_derivs>0) then
    time = self%t(self%nt-2)
  else
    time = self%out_next
  end if
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  !if (self%no_mans_land) then
  !  self%ioformat = 1
  !else
  !  self%ioformat = 2
  !end if
  write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout+1
  write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), self%iout, self%id
  if (io%debug(0)) then
    if (io%debug(2) .or. self%track) &
      print'(a,3x,"time:",2f13.6)', ' patch_t%output_legacy: '//trim(filename), time, self%time
  end if
  open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
  !write (io%data_unit,'(i2.2)') self%ioformat
  write (io%data_unit,'(a)') trim(self%kind)
  write (io%data_unit,'(a)') trim(self%eos)
  write (io%data_unit,'(a)') trim(self%opacity)
  write (io%data_unit,'(1p,2e18.10,6(2x,3e18.10)1x,1x,2e10.3)') time, self%dtime, &
        self%position, self%ds, self%velocity, self%llc_nat, self%llc_cart, &
        self%mesh%centre_nat, self%quality, eos%gamma
  write (io%data_unit,'(6(3i10.1,2x),3i2.1,1x,1i2.1,i3)') int(self%box/self%ds+0.5), &
        u-l+1, self%li, self%ui, self%n, self%gn, self%ng, 3, self%level
  write (io%data_unit,'(3(1i8.1,1x),i2.1)') self%id, self%iout, self%istep, self%mesh_type
  close (io%data_unit)
  !-----------------------------------------------------------------------------
  ! Write a binary file with the data cube.
  !-----------------------------------------------------------------------------
  n = u-l+1
  open (io%data_unit, file=trim(filename)//'.dat', form='unformatted', access='direct', &
    status='unknown', recl=real32*product(n))
  do jo=1,no
    iv = 1
      allocate (buf(n(1),n(2),n(3)))
      buf = 0.0
      !-------------------------------------------------------------------------
      ! If we write out time derivatives we want a pair of consistent values
      ! and time derivatives, which we have at iit(nt-2).
      !-------------------------------------------------------------------------
      b2 => buf(self%li(1):self%ui(1),&
                self%li(2):self%ui(2),&
                self%li(3):self%ui(3))
      call self%get_q(time,b2)
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,2e11.3)') &
          trim(filename),1,self%it,jo,minval(buf),maxval(buf)
      write (io%data_unit,rec=iv+(jo-1)*2) buf
      if (io%verbose>1 .and. iv==self%idx%d) print *,'rho minmax',minval(buf), maxval(buf)
    iv = 2
      buf = self%src(:,:,:,1)
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,2e11.3)') &
          trim(filename),1,self%it,jo,minval(buf),maxval(buf)
      write (io%data_unit,rec=iv+(jo-1)*2) buf
      if (io%verbose>1 .and. iv==self%idx%d) print *,'rho minmax',minval(buf), maxval(buf)
    iv = 3
      buf = self%rk(:,:,:,1)
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,2e11.3)') &
          trim(filename),1,self%it,jo,minval(buf),maxval(buf)
      write (io%data_unit,rec=iv+(jo-1)*2) buf
      if (io%verbose>1 .and. iv==self%idx%d) print *,'rho minmax',minval(buf), maxval(buf)
      deallocate (buf)
  end do
  close (io%data_unit)
   self%iout = self%iout+1
   self%out_next = (nint(self%out_next/io%out_time)+1)*io%out_time
  !$omp end critical (out_cr)
  end do
  call trace_end (itimer)
END SUBROUTINE output_legacy_rt

END MODULE rt_mod
