MODULE rt_mod
  USE bits_mod
  USE eos_mod
  USE io_mod
  USE link_mod
  USE kinds_mod
  USE ray_geometry_v3_mod
  USE rt_components_mod
  USE rt_solver_mod
  USE patch_mod
  USE task_mod
  USE timestep_mod
  USE trace_mod
  USE scaling_mod
  implicit none
  private 
!> =====================================================================
!> rt_t is the main task. 
!> every task has a set of chain links, belonging to this patch. Every
!> item in the set is then connected to appropriate upstream/downstream
!> nbors. Also, any patch might end up being the chain head for some
!> bundle. In that case, it nudges it's own link.
!> =====================================================================
  type, extends(patch_t), public :: rt_t
    class(chain_link_t), dimension(:), pointer :: my_links => null()     ! list of link, which belongs to me
    logical, dimension(:), pointer :: updated => null()                 ! list of flags, indicates if particular direction has been updated
    integer :: wall_idx(6)                                              ! index, where data for particular wall starts
    type(link_t), pointer             :: mhd_link                       ! link to my mhd_patch
    type(rt_outer_bndry_t)            :: rt_boundaries                  ! RT outer boundaries
    real(kind=8)                      :: upd_time                       ! time when the last internal update was done
    integer                           :: n_lambda                       ! number of wavelengths
    integer                           :: n_bundle                       ! number of bundles
    real, dimension(:), pointer       :: w_lambda                       ! wavelength weights
    real, dimension(:), pointer       :: w_bundle                       ! space angle weights
    real                              :: cr_heating
    real, dimension(:,:,:),   pointer :: d     => null()
    real, dimension(:,:,:),   pointer :: e_th  => null()
    real, dimension(:,:,:,:), pointer :: rk    => null()
    real, dimension(:,:,:,:), pointer :: src   => null()
    character(len=16)                 :: solver ='integral'             ! choice of solver
    procedure(solve_integral), pointer:: solve => null()
    logical :: do_output      = .false.
    logical :: data_ready     = .false.
    logical :: do_chains      = .true.
  contains
    procedure :: absolute_point_position
    procedure :: boundary_conditions
    procedure :: cell2pos
    procedure :: cell2shell
    procedure :: density_weight_factor
    procedure :: dnload
    procedure :: execute_no_chain
    procedure :: finish
    procedure :: get_bndry
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
  real, parameter:: pi4 = 16.*atan(1.)
  
CONTAINS

!> =====================================================================
!> The rt_patch initialization subroutine. The patch is responsible for
!> all the directions and their bundles. Al this is stored in chain links.
!> FLOW:
!> 1) read the namelist
!> 2) initialize the patch, set bits
!> 3) set 
!> =====================================================================
SUBROUTINE init_rt (self)
  class (rt_t) :: self
  !.....................................................................
  logical, save :: first_time= .true.
  logical, save :: do_output  = .false.
  real(8), save :: dtime=0.01_8, grace=0.01_8, init_time=0d0
  real,    save :: cr_heating=0.0
  character(len=16) :: solver = 'integral'
  integer :: ignore, nd
  namelist /rt_params/ init_time, dtime, grace, solver, do_output, cr_heating
  ! --------------------------------------------------------------------
  
  call trace_begin ('rt_t%init')
  !$omp critical (init_cr)
  if (first_time) then
    rewind (io%input); read (io%input, rt_params, iostat=ignore)
    if (io%master) write (*, rt_params)
    call rt_solver%init
    first_time=.false.
  end if
  !$omp end critical (init_cr)
  self%periodic = .false.
  self%do_output = do_output
  if (ray_geometry%head%ray_list%directionality(1) == 1) then
    nd = 2
  else
    nd = 1
  end if
  self%nv = ray_geometry%head%n_ray_dirs*nd*eos%n_lambda                ! memory slot stores all heating rates per direction per freq bin
  self%nw = 1                                                           ! no time deriv  
  call self%patch_t%init                                                ! initialize patch
  call self%set (bits%static+bits%no_download)                          ! special bits
  self%dtime = dtime                                                    ! 0.1 second default
  self%time = init_time                                                 ! 0.0 default
  self%dt_fixed = dtime                                                 ! timestep is fixed for radiative transfer
  self%grace = 0.0                                                      ! grace on by default, but very small, so nboring rt_patches would nod do anything stupid
  self%kind = 'rt_solver'                                                ! set the patch kind for output post-processing purposes
  self%guard_zones = .true.
  call self%rt_boundaries%init (self%gn, ray_geometry%head%n_ray_dirs*nd, self%mesh)
  self%do_chains = ray_geometry%head%do_chains                          ! are there chains in the system, or patches are effectively isolated
  self%n_lambda = eos%n_lambda                                          ! number of wavelengths
  if (self%do_chains) then
    self%n_bundle = ray_geometry%head%total_n_bundles * nd                ! number of n_bundles
  else
    self%n_bundle = ray_geometry%head%total_n_bundles
  end if
  self%cr_heating = cr_heating                                          ! store cr heating, if it is present
  
  !> -------------------------------------------------------------------
  !> set up rt-specific parameters : weights for frequency bins and
  !> space-angles. Also, select the appropriate solver.
  !> -------------------------------------------------------------------
  allocate (self%w_lambda(self%nv))                                     ! frequency bin weights
  allocate (self%w_bundle(self%nv))                                     ! space angle weights
  self%w_bundle = pi4/2./(ray_geometry%head%n_ray_dirs * nd)            ! space angle weights FIXME: could/should be different weights for different directions
  
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
  
  allocate(self%updated(self%n_bundle))                                 ! flags, indicating update status of all the bundles
  self%updated(:) = .false.                                             ! none are updated initially
  allocate(self%my_links(self%n_bundle))                                ! allocate all the links this patch contains
  if (self%do_chains) then
    call self%init_chain(nd)
  else
    call self%init_no_chain
  end if
  call trace_end
END SUBROUTINE init_rt

!> =====================================================================
!> Set up all the links, belonging to this patch. Every link represents
!> a separate bundle. One rt patch manages all the directions.
!> =====================================================================
SUBROUTINE init_chain(self, nd)
  class(rt_t) :: self
  integer :: nd
  ! ....................................................................
  class (ray_dirs_t),   pointer :: ray_dir
  class (ray_bndl_t),   pointer :: ray_bnd
  class (packet_t),     pointer :: packet
  class (chain_link_t), pointer :: this_link
  integer :: a_id, ivl, dd, ibnd
  ! --------------------------------------------------------------------
  a_id = 1                                                              ! link personal ID in the link array
  ivl = 1
  do dd = 1,nd                                                          ! separate links for forward/reverse directions
    ray_dir => ray_geometry%head%ray_list%head                            ! pick a ray direction
    do while(associated(ray_dir))                                         ! loop over all ray directions
      do ibnd =1, ray_dir%n_bundles                                       ! loop over ray bundles
        ray_bnd => ray_dir%bundle(ibnd)                                                                              
        this_link => self%my_links(a_id)
        this_link%link => self%link                                     ! links should know, what rt patch hosts them
        this_link%direction = dd                                        ! forward/reverse direction links
        this_link%id = a_id                                             ! link personal ID in the link array
        this_link%dir_id = ray_bnd%dir_id
        this_link%bundle => ray_bnd                                     ! link has to know, which bundle it corresponds to
        this_link%directionality = ray_dir%directionality(1)
        this_link%ivl = ivl
        this_link%ivu = ivl + eos%n_lambda-1
        self%w_lambda(ivl:ivl+eos%n_lambda-1) = eos%w_lambda
        allocate(this_link%packet_data(ray_bnd%n_rays))                 ! allocate intensity in data for the link. Array consists of all the rays
        call this_link%get_ds(self%ds,this_link%ds)
        a_id = a_id + 1                                                 ! next ID
      end do
      ray_dir => ray_dir%next                                           ! next direction
      ivl = ivl + eos%n_lambda
    end do
  end do
  nullify(ray_dir, ray_bnd)
END SUBROUTINE init_chain

!> =====================================================================
!> Set up all the links, belonging to this patch. Every link represents
!> a separate bundle. One rt patch manages all the directions.
!> =====================================================================
SUBROUTINE init_no_chain(self)
  class(rt_t) :: self
  ! ....................................................................
  class (ray_dirs_t),   pointer :: ray_dir
  class (ray_bndl_t),   pointer :: ray_bnd
  class (packet_t),     pointer :: packet
  class (chain_link_t), pointer :: this_link
  integer :: a_id, ivl, ivl2
  ! --------------------------------------------------------------------
  a_id = 1                                                              ! link personal ID in the link array
  ivl = 1
  ivl2 = ray_geometry%head%n_ray_dirs*eos%n_lambda + 1
  ray_dir => ray_geometry%head%ray_list%head                            ! pick a ray direction
  do while(associated(ray_dir))                                         ! loop over all ray directions
    ray_bnd => ray_dir%bundle(1)                                                                              
    this_link => self%my_links(a_id)
    this_link%link => self%link                                         ! links should know, what rt patch hosts them
    this_link%direction = 3                                             ! bi-directional
    this_link%id = a_id                                                 ! link personal ID in the link array
    this_link%dir_id = ray_bnd%dir_id
    this_link%bundle => ray_bnd                                         ! link has to know, which bundle it corresponds to
    this_link%directionality = ray_dir%directionality(1)
    this_link%ivl = ivl
    this_link%ivu = ivl + eos%n_lambda-1
    this_link%ivl2 = ivl2
    this_link%ivu2 = ivl2 + eos%n_lambda-1
    self%w_lambda(ivl:ivl+eos%n_lambda-1)   = eos%w_lambda
    self%w_lambda(ivl2:ivl2+eos%n_lambda-1) = eos%w_lambda
    call this_link%get_ds(self%ds,this_link%ds)
    a_id = a_id + 1                                                     ! next ID
    ray_dir => ray_dir%next                                             ! next direction
    ivl = ivl + eos%n_lambda
    ivl2 = ivl2 + eos%n_lambda
  end do
  nullify(ray_dir, ray_bnd)
END SUBROUTINE init_no_chain

!> =====================================================================
!> =====================================================================
SUBROUTINE dnload (self)
  class (rt_t) :: self
END SUBROUTINE dnload

!> =====================================================================
!> This subroutine maps the incoming heating rates to the outer boundary
!> conditions
!> =====================================================================
SUBROUTINE map_bndry(self, link)
  class(rt_t) :: self
  class(chain_link_t) :: link
  !.....................................................................
  class(packet_t), pointer :: packet
  real, dimension(:,:,:), pointer:: wall
  integer :: ray, rid, c(3), pt, a,b, i, ipck, ng(3), it
  logical :: selected
  !---------------------------------------------------------------------
  rid=1
  do ipck =1, link%bundle%n_packets
    packet => link%bundle%packet(ipck)                                  ! loop over ray packets in the bundle
    do ray = 1, packet%n_rays                                 
      allocate(link%packet_data(rid)%in(eos%n_lambda,1))
      rid = rid + 1
    end do
  end do  
END SUBROUTINE map_bndry

!> =====================================================================
!> This subroutine maps the incoming heating rates to the outer boundary
!> conditions
!> =====================================================================
SUBROUTINE get_bndry(self, link, q_fw, q_re)
  class(rt_t) :: self
  class(chain_link_t) :: link
  real, dimension(:,:), pointer :: q_fw
  real, dimension(:,:), pointer, optional :: q_re
  !.....................................................................
  class(packet_t), pointer :: packet
  real, dimension(:,:,:), pointer:: wall
  integer :: ray, rid, c(3), pt, dir, ipck, ivl, ivu, ivl2, ivu2
  logical :: selected
  !---------------------------------------------------------------------
  if (self%do_chains) then
    ivl=link%ivl
    ivu=link%ivu
    select case(link%direction)
    case(1)
      rid=1
      do ipck =1, link%bundle%n_packets
        packet => link%bundle%packet(ipck)                              ! loop over ray packets in the bundle
        do ray = 1, packet%n_rays                                 
          c = packet%rev_coords(:,ray)
          q_fw(rid,:) = self%mem(c(1),c(2),c(3),ivl:ivu,self%it,1)
          rid = rid + 1
        end do
      end do
    case(2)
      rid=1
      do ipck =1, link%bundle%n_packets
        packet => link%bundle%packet(ipck)                              ! loop over ray packets in the bundle
        do ray = 1, packet%n_rays                                 
          c = packet%fwd_coords(:,ray)
          q_fw(rid,:) = self%mem(c(1),c(2),c(3),ivl:ivu,self%it,1)
          rid = rid + 1
        end do
      end do
    end select
  else
    ivl =link%ivl
    ivu =link%ivu
    ivl2=link%ivl2
    ivu2=link%ivu2
    rid=1
    do ipck =1, link%bundle%n_packets
      packet => link%bundle%packet(ipck)                                ! loop over ray packets in the bundle
      do ray = 1, packet%n_rays                                 
        c = packet%rev_coords(:,ray)
        q_fw(rid,:) = self%mem(c(1),c(2),c(3),ivl :ivu ,self%it,1)
        c = packet%fwd_coords(:,ray)
        q_re(rid,:) = self%mem(c(1),c(2),c(3),ivl2:ivu2,self%it,1)
        rid = rid + 1
      end do
    end do
  end if
END SUBROUTINE get_bndry

SUBROUTINE cell2shell(self, c_coord, s_coord, wall)
  class(rt_t) :: self
  integer :: c_coord(3), s_coord(2), wall
  !.....................................................................
  integer :: c(3)
  !---------------------------------------------------------------------
  c = c_coord !- self%li + 1
  if (c_coord(1) == self%li(1)) then
    s_coord(1) = c(2); s_coord(2) = c(3)
    wall = -1
  elseif (c_coord(1) == self%ui(1)) then
    s_coord(1) = c(2); s_coord(2) = c(3)
    wall = 1
  elseif (c_coord(2) == self%li(2)) then
    s_coord(1) = c(3); s_coord(2) = c(1)
    wall = -2
  elseif (c_coord(2) == self%ui(2)) then
    s_coord(1) = c(3); s_coord(2) = c(1)
    wall = 2
  elseif (c_coord(3) == self%li(3)) then
    s_coord(1) = c(1); s_coord(2) = c(2)
    wall = -3
  elseif (c_coord(3) == self%ui(3)) then
    s_coord(1) = c(1); s_coord(2) = c(2)
    wall = 3
  end if
END SUBROUTINE cell2shell

!> =====================================================================
!> Solve the RT problem, by first obtaining the thermodynamics at the
!> designated RT time, and using them to compute the opacity and source
!> function values. Then nudge the chain heads.
!> =====================================================================
SUBROUTINE update_rt (self)
  class (rt_t) :: self
  class(task_t),  pointer:: task
  class(patch_t), pointer:: mhd
  integer, save:: itimer=0
  integer :: ii, iy
  real    :: g1
  !---------------------------------------------------------------------
  if (self%data_ready) return
  
  call trace%begin ('rt_t%update', itimer=itimer)
 ! if (io%verbose > 1) print *, 'RADIATIVE TRANSFER TIME', self%time
  
  if (self%time < 0d0) then
    self%time  = 0d0 !+ self%dtime                                      ! set current time, if it was negative (to do the first RT run before MHD updates)
  end if
  
  task => self%mhd_link%task                                            ! pick the rmhd task
  select type (task)
  class is (patch_t)
    mhd => task
  end select
  
  self%upd_time = self%time
  !> -------------------------------------------------------------------
  !> allocate all the scratch arrays
  !> -------------------------------------------------------------------
  allocate (self%d (self%gn(1), self%gn(2), self%gn(3)))                ! density scratch
  allocate (self%e_th(self%gn(1), self%gn(2), self%gn(3)))              ! thermal energy scratch
  allocate (self%rk  (self%gn(1), self%gn(2), self%gn(3), self%n_lambda))! opacity scratch
  allocate (self%src (self%gn(1), self%gn(2), self%gn(3), self%n_lambda))! source function scratch
  
  !> -------------------------------------------------------------------
  !> interpolate in time and space to get densoties and internal energies
  !> -------------------------------------------------------------------
  call mhd%log_interpolate (self%time, mhd%i_d, self%d)               
  call mhd%log_interpolate (self%time, mhd%i_e, self%e_th)
  self%d = exp(self%d)
  self%e_th    = exp(self%e_th)               
  if (mhd%kind(1:9)=='stagger2e') then
    self%e_th = self%e_th / self%d
  else if (mhd%kind(1:9)=='stagger2_') then
    g1 = mhd%gamma-1.0
    self%e_th = self%d**g1*exp(self%e_th*g1/self%d)/g1
    self%e_th = self%e_th/self%d
  else
    print *,trim(mhd%kind)
    error stop "rt_mod (1): MHD solver type not known"
  end if
  print *, mhd%kind

  call eos%lookup (shape(self%e_th), x=self%e_th, y=self%d, &           ! call the table lookup, ...
                   rk=self%rk, src=self%src)  
  print *, self%src(:,8,8,:)
  print *, 'rk',self%rk(:,8,8,:)
  self%data_ready =.true.                                               ! flag, that internal data is ready
  
  if(io%verbose >1) then
    print *, self%id, 'rk', minval(self%rk), maxval(self%rk)
    print *, self%id, 'src', minval(self%src), maxval(self%src)
  end if
  self%mem(:,:,:,:,self%it,1) = 0.0
  call self%rt_boundaries%condition(self%mem(:,:,:,:,self%it,1))
  call trace%end (itimer)
  if (self%do_chains) then
    call self%nudge                                                     ! nudge the chain heads
  else
    call self%execute_no_chain
  end if
END SUBROUTINE update_rt 

!> =====================================================================
!> Nudge the chain heads. If all patches in the chain become ready, then
!> solve for the net heating rate, storing it in a slot that rotates
!> with time for each patch in the chain.
!> =====================================================================
SUBROUTINE nudge (self)
  class(rt_t) :: self
  !.....................................................................
  class(chain_link_t), pointer :: link
  class(task_t), pointer :: tmp
  integer :: b_id
  !---------------------------------------------------------------------
  do b_id = 1, self%n_bundle
    link => self%my_links(b_id)
    
  !  if (io%verbose > 4) then
  !    print *, 'link no.',r_id
  !    print *, 'before nudge, n_ready',link%head%n_ready, link%head%l_chain
  !    print *,'nudged by master', link%master,'rt time',self%time
  !    print *, 'other links updated',self%updated
  !    print *, 'my data is ready',self%data_ready
  !    print *, 'my ID is:',self%id,'link', link%id
  !  end if
    
  !  tmp => link%head%link%task
  !  select type(tmp)
  !  class is (rt_t)
  !  end select
    
  !  if (io%verbose > 4) print *, 'nudging task ID:', tmp%id,'link id', link%head%id, &
  !  'nudged allocated:',associated(link%head%nudged),'size'
  !  if (io%verbose > 4) print *, size(link%head%nudged),'its content:',link%head%nudged, 'head assoc', associated(link%head)
  !  if (io%verbose > 4) print *, 'l_chain',link%head%l_chain
    !$omp atomic update
    link%head%n_ready = link%head%n_ready+1                             ! nudge the master link
  !  if (link%head%nudged(link%chain_nm)) print *,'this already was nudged..'
  !  if (link%head%n_ready > link%head%l_chain) print *,'double nudging?'
  if ( link%head%nudged(link%chain_nm)) then                            ! mark myself as done in the ready array
    error stop "Whoops: nudging the nudged"
  end if
    !$omp atomic write
    link%head%nudged(link%chain_nm) = .true.                            ! mark myself as done in the ready array
        
    if (all(link%head%nudged)) call execute_chain(link)                  ! execute the chain, when all links are ready
  end do
END SUBROUTINE nudge

!> =====================================================================
!> Loop through all links in the chain and solve the actual rt
!> =====================================================================
SUBROUTINE execute_chain (self)
  class(chain_link_t), pointer :: self
  !.....................................................................
  class(rt_t),          pointer :: rt_patch
  class(task_t),        pointer :: task
  class(chain_link_t),  pointer :: link, prev
  real, dimension(:,:), pointer :: q_in
  integer :: ray, it
  integer, save:: itimer=0
  logical :: first
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%execute_chain', itimer=itimer)
  link => self%head                                                     ! take the master link
  nullify(prev)                                                         ! just in case - there has to be no "prev" for the master link
  task => link%link%task                                            ! get the patch of current chain link
  select type (task)
  class is (rt_t)
  rt_patch => task
  end select
  allocate(q_in (size(link%packet_data),rt_patch%n_lambda))                 ! allocate the intensity scratch
  q_in = 0
  first = .true.
  do while (associated(link))
    if (associated(prev)) then
      task => prev%link%task                                              ! get the patch of upstream chain link
      select type (task)
      class is (rt_t)
      rt_patch => task
      end select
      if  (rt_patch%upd_time < rt_patch%time) then
        it = rt_patch%iit(rt_patch%nt-2)
      else
        it = rt_patch%it
      end if
      do ray=1, link%bundle%n_rays
        q_in(ray,:) = link%packet_data(ray)%in(:,it)
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
      call rt_patch%get_bndry(link, q_in)
    end if
    !
    call rt_patch%solve(link, q_in)                                     ! do the radiative transfer for this particular link                                                              
    rt_patch%updated(link%id) = .true.                                  ! mark, that this link has been updated
    if(all(rt_patch%updated)) call rt_patch%finish                      ! if all links are done, this rt_patch is complete
    prev => link
    link => link%next
  end do
  deallocate(q_in)                                                      ! deallocate the scratch
  self%head%nudged = .false.                                                 ! reset the nudge logicals
  self%head%n_ready = 0                                                      ! no links are ready in this chain
  nullify(rt_patch, task, prev)  
  call trace%end(itimer)
END SUBROUTINE execute_chain

!> =====================================================================
!> Loop through all links in the chain and solve the actual rt
!> =====================================================================
SUBROUTINE execute_no_chain (self)
  class(rt_t) :: self
  ! ....................................................................
  class(chain_link_t),  pointer :: link
  real, dimension(:,:), pointer :: q_fw, q_re
  integer :: il
  ! --------------------------------------------------------------------
  call self%dnload
  do il = 1, self%n_bundle
    link => self%my_links(il)
    allocate(q_fw (link%bundle%n_rays,self%n_lambda))
    allocate(q_re (link%bundle%n_rays,self%n_lambda))
    call self%get_bndry(link, q_fw, q_re)                               ! get the boundary data if applicable
    call self%solve(link, q_fw, q_re)                                   ! do the radiative transfer for this particular link           
    deallocate (q_fw, q_re)
  end do
  call self%finish
END SUBROUTINE execute_no_chain

!> =====================================================================
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
!> =====================================================================
SUBROUTINE link_connects_to (self, nbor, mine, other, connects, is_upstream)
  class (rt_t) :: self
  class (rt_t), pointer :: nbor
  class (chain_link_t) :: mine
  class (chain_link_t) :: other
  logical :: connects, is_upstream
  !.....................................................................
  integer :: a, b, c, mrn, ip, ray, first_pt(3), cp(3), cl(3), cu(3)
  real(kind=8) :: ppos
  
  integer :: pp,  rev(3), rid, pid
  integer :: cell(3), i, s_coord(2),face, ics, pl
  real :: iw
  logical :: tout(3), within(3)
  !---------------------------------------------------------------------
  connects = .false.
  is_upstream = .false.
  tout=.false.
  
  if (mine%bundle%dir_ID /= other%bundle%dir_ID)     return             ! different ray spatial direction
  if (mine%direction /= other%direction)             return             ! different direction
    
  if (self%level /= nbor%level) then
    if (abs(self%level-nbor%level) < 2) then
      select case(mine%direction)
      case(1)
        if (mine%bundle%lwall(1) /= 0) then
          a = 2; b = 3; c=1
        else if (mine%bundle%lwall(2) /= 0) then
          a = 3; b = 1; c=2
        else if (mine%bundle%lwall(3) /= 0) then
          a = 1; b = 2; c=3
        end if

      case(2)
        if (mine%bundle%uwall(1) /= 0) then
          a = 2; b = 3
        else if (mine%bundle%uwall(2) /= 0) then
          a = 3; b = 1
        else if (mine%bundle%uwall(3) /= 0) then
          a = 1; b = 2
        end if
      end select

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
          

     error stop "rt_mod: error 10" 
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
    select case(mine%direction)                                           
    case(1)
      if (mine%bundle%b_id_rev  /= other%bundle%ID) return              ! wrong bundle to make a chain
        first_pt =  mine%bundle%packet(1)%coords(:,1,ray)                 ! taking the coordinates of the first point in the ray
        if (self%no_mans_land) then
          rid = mine%bundle%packet(1)%r_id_rev(ray)                         ! ray id in a packet list of a ray in the previous bundle 
          pid = mine%bundle%packet(1)%p_id_rev(ray)                         ! packet id in the previous bundle
          rev = other%bundle%packet(pid)%fwd_coords(:,rid)                ! zeroth point in the ray (location of incoming intensity)               
        else
          rid = mine%bundle%packet(1)%r_id_rev(ray)                         ! ray id in a packet list of a ray in the previous bundle 
          pid = mine%bundle%packet(1)%p_id_rev(ray)                         ! packet id in the previous bundle
          pp = other%bundle%packet(pid)%n_points
          rev = other%bundle%packet(pid)%coords(:,pp,rid)                  ! zeroth point in the ray (location of incoming intensity)               
        end if
      rid = mine%bundle%packet(1)%r_id_rev(ray)                         ! ray id in a packet list of a ray in the previous bundle 
      pid = mine%bundle%packet(1)%p_id_rev(ray)                         ! packet id in the previous bundle
      rev = other%bundle%packet(pid)%fwd_coords(:,rid)                ! zeroth point in the ray (location of incoming intensity)                                  
    case(2)
      if (mine%bundle%b_id_fw /= other%bundle%ID) return                ! wrong bundle to make a chain
      pp = mine%bundle%packet(1)%n_points
      first_pt =  mine%bundle%packet(1)%coords(:,pp,ray)                ! taking the coordinates of the last point in the ray
      if (self%no_mans_land) then
        rid = mine%bundle%packet(1)%r_id_fw(ray)                          ! ray id in a packet list of a ray in the previous bundle 
        pid = mine%bundle%packet(1)%p_id_fw(ray)                          ! packet id in the previous bundle
        rev = other%bundle%packet(pid)%rev_coords(:,rid)                ! zeroth point in the ray (location of incoming intensity)            
      else
        rid = mine%bundle%packet(1)%r_id_fw(ray)                          ! ray id in a packet list of a ray in the previous bundle 
        pid = mine%bundle%packet(1)%p_id_fw(ray)                          ! packet id in the previous bundle
        rev = other%bundle%packet(pid)%coords(:,1,rid)                    ! zeroth point in the ray (location of incoming intensity)            
      end if
                                                   
    end select
    first_pt = nbor%pos2cell(self%myposition(first_pt))

    do pp = 1,3
      if (first_pt(pp) == rev(pp)) tout(pp) = .true.
    end do
    
    if (all(tout)) connects = .true.
  end if
  
END SUBROUTINE link_connects_to

!> =====================================================================
!> Find the absolute position of index 'ray' in ray packet 'rays'.
!> If the optional bn parameter is present, the answer depends on
!> the direction, otherwise 'pt' chooses the point along the ray.
!> =====================================================================
FUNCTION absolute_point_position (self, link, packet, ray, pt, bn) result (cell_pos)
  class (rt_t) :: self
  class (chain_link_t), pointer :: link
  class (packet_t), pointer :: packet
  integer, optional :: bn
  integer :: ray, pt, cp(3)
  real(kind=8) :: cell_pos(3)
!-----------------------------------------------------------------------
  call trace_begin ('rt_t%absolute_point_position', 4)
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
  call trace_end
END FUNCTION absolute_point_position

LOGICAL FUNCTION has_OB (self, link)
  class (rt_t) :: self
  class (chain_link_t) :: link
  integer :: i(3), d(3),t
  real(8) :: p(3)
  !.....................................................................
  has_OB = .false.
  d = link%bundle%dstep
  if (link%direction == 1) then
    i = link%bundle%packet(1)%rev_coords(:,1)
  else 
    i = link%bundle%packet(1)%fwd_coords(:,1)
  end if
    p = self%myposition(i)
  
  do t = 1,3
    if (link%direction ==1) then
      if (d(t)>0) then
        if (p(t)<rt_boundaries%llc(t)) has_OB = .true.
      else
        if (p(t)>rt_boundaries%urc(t)) has_OB = .true.
      end if
    else
      if (d(t)>0) then
        if (p(t)>rt_boundaries%urc(t)) has_OB = .true.
      else
        if (p(t)<rt_boundaries%llc(t)) has_OB = .true.
      end if
    end if
  end do

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
  integer :: ii
  !---------------------------------------------------------------------
  self%updated(:)=.false.
  self%data_ready = .false.
  if (self%do_output) call self%output
  deallocate (self%d, self%e_th, self%rk, self%src)                     ! deallocate scratch
END SUBROUTINE finish

!> =====================================================================
!> This procedure should loop over all bundle directions, spawning two
!> tasks for each direction, with solvers running in the forward and
!> reverse directions.
!>
!> The procedure implies an implicit wait for the upstream net intensity
!> which arrives at different times for the two directions.  The first
!> call increments the counter from 0 to 1, and when the 2nd call makes
!> it 2, the procedure knows that temporary work arrays should be
!> deallocated.
!>
!> In this prototype case we use only one (vertical = y) direction, but
!> in the general case we loop over the included ray directions, doing
!> a waittask after the loop, so all tasks are active at the same time.
!> Here we spawn one task per wavelength, just to demonstrate tasks.
!>
!> The tasks that are missing upstream values need to to a taskwait,
!> and try again later.  However, in the current special case we always
!> have both boundart values, so there are no waits.
!> =====================================================================
SUBROUTINE solve_integral (self, link, q_in_fw, q_in_re)
  class (rt_t) :: self
  class(chain_link_t), pointer  :: link
  real, dimension(:,:), pointer, optional :: q_in_fw
  real, dimension(:,:), pointer, optional :: q_in_re
  !.....................................................................
  class(packet_t),        pointer :: packet
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q_fw, q_re
  real, dimension(:,:),   pointer :: dwf, q_fw0, q_re0
  integer                         :: nrays, npt, nwb, ipck, iray, ir
  integer, save                   :: itimer=0
  class(task_t), pointer :: tmp
  ! --------------------------------------------------------------------
  !call trace_begin ('rt_t%solve_integral', itimer=itimer)

  nwb   =  self%n_lambda                                                ! number of wl bins
  iray = 1
  do ipck = 1, link%bundle%n_packets                                    ! loop over packets
    packet => link%bundle%packet(ipck)                                  ! pick a packet
    nrays =  packet%n_rays                                              ! number of rays in this packet
    !> -----------------------------------------------------------------
    !> RT solvers need at least 3 points to work with, as they use
    !> derivatives. Rays, shorter, than 3 points, have to be extended
    !> to patch boundaries to have accurate solutions.
    !> TODO: inner boundaries modify npt
    !> -----------------------------------------------------------------
    npt = packet%n_points + 2

 
    !> -----------------------------------------------------------------
    !> allocate all the scratches
    !> -----------------------------------------------------------------
    allocate(s     (nrays, nwb, npt))                                   ! source function
    allocate(rk    (nrays, nwb, npt))                                   ! opacity
    allocate(dtau  (nrays, nwb, npt))                                   ! optical depth between the points
    allocate(exp0  (nrays, nwb, npt))                                   ! first exponential factor
    allocate(dwf   (nrays,      npt))                                   ! density weight factor
    if (present(q_in_fw)) then
      allocate(q_fw0(nrays, nwb))                                       ! incoming forward rays
      allocate(q_fw (nrays, nwb, npt))                                  ! forward heating rates
    end if
    if (present(q_in_re)) then
      allocate(q_re0(nrays, nwb))                                       ! incoming reverse rays
      allocate(q_re (nrays, nwb, npt))                                  ! reverse heating rates
    end if
       
    !> -----------------------------------------------------------------
    !> map the quantities onto rays
    !> -----------------------------------------------------------------
    call self%ray_vol_map4(packet, self%src, s)                         ! map source function to rays
    call self%ray_vol_map4(packet, self%rk, rk)                         ! map opacities to rays
    call self%optical_depth_staggered(link%ds, rk, dtau)                
    call self%ray_vol_map3(packet, self%d, dwf)                         ! density weight factor
    dwf = dwf * link%ds
    if (present(q_in_fw)) then
      do ir=1, nrays
        q_fw0(ir,:) = q_in_fw(iray+ir-1,:)/dwf(ir,1)                    ! take the appropriate heating rates
      end do
    end if
    if (present(q_in_re)) then
      do ir=1, nrays
        q_re0(ir,:) = q_in_re(iray+ir-1,:)/dwf(ir,npt)                  ! take the appropriate heating rates
      end do
    end if
#if defined __INTEL_COMPILER
      exp0 = 2**(-dtau*1.44269504) 
#else
      exp0 = exp(-dtau)                                                 ! expopential factors
#endif        
    if (self%do_chains) then
      if (link%direction ==1) then
        call rt_solver%integral_forward (dtau, exp0, s, q_fw0, q_fw)    ! solve forward
        call self%update_q(packet, dwf, link%id, q_fw=q_fw)             ! update the heating rates
        deallocate (q_fw0, q_fw)
      else
        call rt_solver%integral_reverse (dtau, exp0, s, q_re0, q_re)    ! solve reverse
        call self%update_q(packet, dwf, link%id, q_re=q_re)             ! update the heating rates
        deallocate (q_re0, q_re)
      end if
    else
      call rt_solver%integral_double(dtau, exp0, S, q_fw0, q_re0, q_fw, q_re)
      call self%update_q(packet, dwf, link%id,q_fw=q_fw, q_re=q_re)
      deallocate (q_fw0, q_re0, q_fw, q_re)
    end if
    deallocate (s, rk, dtau, exp0, dwf)
    iray = iray + nrays - 1
  end do

  !call trace_end (itimer)
END SUBROUTINE solve_integral

!> =====================================================================
!> Same as solve_integral, but for the Feautrier formulation
!> =====================================================================
SUBROUTINE solve_feautrier (self, link, i_in, i_in2)
  class (rt_t) :: self
  class(chain_link_t), pointer  :: link
  real, dimension(:,:), pointer, optional :: i_in
  real, dimension(:,:), pointer, optional :: i_in2
  !.....................................................................
  class(packet_t),        pointer :: packet
  integer                         :: nrays, npt, nwb, ipck
  real, dimension(:,:), pointer :: dwf
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve_feautrier', itimer=itimer)
  
  if (.not.present(i_in2)) error stop "solve Feautrier: needs 2 directions to act"
  error stop "Feautrier formulation not set up properly yet"
  do ipck=1, link%bundle%n_packets
    packet => link%bundle%packet(ipck)
    nrays =  packet%n_rays
    npt   =  packet%n_points
    allocate(s   (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(rk  (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(dtau(nrays, nwb, npt+2))
    allocate(exp0(nrays, nwb, npt+2))
    allocate(q   (nrays, nwb, npt+2))
    allocate(dwf (nrays,      npt+2))
    call self%ray_vol_map4(packet, self%src, s)                         ! map source function to rays
    call self%ray_vol_map4(packet, self%rk, rk)                         ! map opacities to rays
    call self%optical_depth_staggered(link%ds, rk, dtau)                ! optical depth is staggered for Feautrier formulation
    call self%ray_vol_map3(packet, self%d, dwf)                         ! density weight factor
    dwf = dwf * link%ds  
    if (io%verbose>1) &
        print *,'dtau',self%id, minval(dtau), maxval(dtau)
    
    call rt_solver%feautrier (dtau, s, q, I_in, I_in2,.true.)         ! solve bi-directional
    call self%update_q(packet, dwf, link%id, q)
    
    deallocate (s, rk, dtau, exp0, q, dwf)
  end do

  call trace_end (itimer)
END SUBROUTINE solve_feautrier

!> =====================================================================
!> Map cell data onto rays
!> =====================================================================
SUBROUTINE ray_vol_map4 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:,:) , pointer, intent(in)  :: in
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), p2(3), ray, pt, nwb, wl, lp
  !---------------------------------------------------------------------
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
      !DIR$ SIMD
      do wl =1, nwb                                                     ! loop over wl bins
        out(ray,wl,pt+1) = in(p(1),p(2),p(3),wl)                        ! first point is missing
      end do
      !DIR$ SIMD
    end do
    p  = packet%rev_coords(:,ray)
    p2 = packet%fwd_coords(:,ray)
    do wl=1, nwb                                                        ! fill in the first and last points
      out(ray,wl,1)  = in(p(1),p(2),p(3), wl)                           ! Get the first point
      out(ray,wl,lp) = in(p2(1),p2(2),p2(3), wl)                        ! Get the last point
    end do
  end do

    
  if (io%verbose > 1) then
    if (any(in==0)) then
      print *, in
      error stop "ray vol map4 in=0"
    end if
    if (any(out==0)) then
      do pt = 1, size(out,3)
      print *,"pt",pt, out(:,:,pt)
      end do
      error stop "ray vol map4 out=0"
    end if
  end if
  
END SUBROUTINE ray_vol_map4



!> =====================================================================
!> Map cell data onto rays. In this case, it is wavelenght independent.
!> =====================================================================
SUBROUTINE ray_vol_map3 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:) , pointer, intent(in)  :: in
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), ray, pt, lp
  !---------------------------------------------------------------------
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
  
END SUBROUTINE ray_vol_map3


SUBROUTINE update_q (self, packet, dwf, ib, q_fw, q_re)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:,:), pointer, optional :: q_fw, q_re
  real, dimension(:,:,:), pointer :: q
  real, dimension(:,:), pointer :: dwf
  real, dimension(:,:,:,:), pointer :: qq, qq2
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), wl, ray, pt, ib, ivl, ivu, it, ivl2, ivu2
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  !call trace_begin ('rt_t%update_q', 4, itimer=itimer)
  ivl = self%my_links(ib)%ivl
  ivu = self%my_links(ib)%ivu
  if  (self%upd_time < self%time) then
    it = self%iit(self%nt-2)
  else
    it = self%it
  end if
  if (self%do_chains) then
    if (present(q_fw)) q => q_fw
    if (present(q_re)) q => q_re
    qq => self%mem(:,:,:,ivl:ivu,it,1)                                    ! point to the memory slot
    cell => packet%coords                                                 ! point to coordinates
    do ray = 1, packet%n_rays                                             ! loop over rays
      do pt = 1, packet%n_points                                          ! loop over points
        p = cell(:,pt,ray)                                                ! get the cell position
        do wl =1, size(q,2)
          qq(p(1),p(2),p(3),wl) = &
          q(ray,wl,pt+1) * dwf(ray,pt+1)
        end do
      end do
    end do
  else
    ivl2 = self%my_links(ib)%ivl2
    ivu2 = self%my_links(ib)%ivu2
    qq  => self%mem(:,:,:,ivl :ivu ,it,1)                               ! point to the forward direction memory slot
    qq2 => self%mem(:,:,:,ivl2:ivu2,it,1)                               ! point to the reverse direction memory slot
    cell => packet%coords                                               ! point to coordinates
    do ray = 1, packet%n_rays                                           ! loop over rays
      do pt = 1, packet%n_points                                        ! loop over points
        p = cell(:,pt,ray)                                              ! get the cell position
        do wl =1, size(q_fw,2)
          qq (p(1),p(2),p(3),wl) = &
          q_fw(ray,wl,pt+1) * dwf(ray,pt+1)
          qq2(p(1),p(2),p(3),wl) = &
          q_re(ray,wl,pt+1) * dwf(ray,pt+1)
        end do
      end do
    end do
  end if
  !call trace_end (itimer)
END SUBROUTINE update_q

!> =====================================================================
!> Compute the optical depth along rays
!> =====================================================================
SUBROUTINE optical_depth_centered (self, ds, rk, dtau)
  class (rt_t) :: self
  real(kind=8) :: ds
  real, dimension(:,:,:), pointer :: rk
  real, dimension(:,:,:), pointer :: dtau
  !.....................................................................
  integer                         :: pt
  !---------------------------------------------------------------------
  do pt = 1, size(rk,3)
    dtau(:,:,pt) = ds*rk(:,:,pt)
  end do

END SUBROUTINE optical_depth_centered


!> =====================================================================
!> Compute the optical depth along rays. Dtau is staggered, therefore it
!> has only a half value on the lower edge of the patch.
!> =====================================================================
SUBROUTINE optical_depth_staggered (self, ds, rk, dtau)
  class (rt_t) :: self
  real(kind=8) :: ds
  real, dimension(:,:,:), pointer :: rk
  real, dimension(:,:,:), pointer :: dtau,dt
  !.....................................................................
  integer                         :: pt
  !---------------------------------------------------------------------
  allocate (dt(size(dtau,1),size(dtau,2),size(dtau,3)))
  do pt=1,size(rk,3)
      dtau(:,:,pt) = rk(:,:,pt)* ds
  end do
  dt=dtau
  do pt=2,size(rk,3)-1
    dtau(:,:,pt) = (dtau(:,:,pt) + dt(:,:,pt-1)) * .5
  end do
  deallocate(dt)
END SUBROUTINE optical_depth_staggered

!> =====================================================================
!> Compute the optical depth along rays. Dtau is staggered, therefore it
!> has only a half value on the lower edge of the patch.
!> =====================================================================
SUBROUTINE density_weight_factor (self, ds, dwf)
  class (rt_t) :: self
  real(kind=8) :: ds
  real, dimension(:,:,:), pointer :: dwf
  !.....................................................................
  integer                         :: iy, iz
  !---------------------------------------------------------------------
  do iz=1,size(dwf,3)
    do iy=1,size(dwf,2)
      dwf(:,iy,iz) = 1./(self%d(:,iy,iz) * ds)
    end do
  end do
END SUBROUTINE density_weight_factor

!> =====================================================================
!> Get heating rates for given time in terms of internal energy. Heating
!> rates are interpolated in time. 
!> =====================================================================
SUBROUTINE get_q (self, time, de)
  class (rt_t) :: self
  real(8):: time
  real(kind=KindScalarVar), dimension(:,:,:), pointer :: de
  !.....................................................................
  integer       :: jt(2), idir, nd, ind, il, iv
  real          :: pt(2)
  integer, save :: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rt_t%get_q', itimer=itimer)
  call self%time_interval (time, jt, pt)
  do iv = 1, self%nv
    de= de + (self%mem(:,:,:,iv,jt(1),1)*pt(1) + &
              self%mem(:,:,:,iv,jt(2),1)*pt(2)) * self%w_lambda(iv) * self%w_bundle(iv)
  end do
  call trace%end (itimer)
END SUBROUTINE get_q

!> =====================================================================
!> Return the absolute position of integer index triple 'cell'
!> =====================================================================
FUNCTION cell2pos (self, cell) result(p)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  p = self%position - 0.5_8*self%size + (cell-self%li+1)*self%ds
END FUNCTION cell2pos

!> =====================================================================
!> Convert absolute position to integer index triple 'cell'
!> =====================================================================
FUNCTION pos2cell (self, p) result (cell)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  !cell = nint((p-self%position+0.5_8*self%size)/self%ds)+self%li-1
  
  cell = nint((p-self%position)/self%ds+self%offset)  
  
END FUNCTION pos2cell

!> =====================================================================
!> Convert absolute position to integer index triple 'cell'
!> currently simple spatial interpolation is used and should be sufficient.
!> if it is not accurate enough, it should be changed to something else.
!> =====================================================================
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

!> =====================================================================
!> Apply boundary conditions
!> =====================================================================
SUBROUTINE boundary_conditions (self)
  class(rt_t):: self
  !---------------------------------------------------------------------
END SUBROUTINE boundary_conditions

!> =====================================================================
!> Write results to disk
!> =====================================================================
SUBROUTINE output_legacy_rt (self, experiment_name)
  USE iso_fortran_env, only : real32
  class(rt_t):: self
  character(len=64) filename, logname, experiment_name
  optional:: experiment_name
  !.............................................................................
  integer:: iv, li(3)=0, ui(3)=0, l(3), u(3), no, jo, it, i
  real:: ds(3)
  logical:: first_time=.true.
  integer:: it1, it2, n(3), idir, il, ind, nd
  real:: pt, qt
  real, dimension(:,:,:), pointer:: buf
  real(8):: time
  !-----------------------------------------------------------------------------
  call trace_begin('patch_mod::output_legacy')
  do while (self%time >= self%out_next)
  !$omp critical (out_cr)
  print *, self%id,'RT legacy output, TIME:', self%upd_time
  self%max_files = self%max_files-1
  if (self%max_files <= 0) then
    write (io_unit%log,*) 'TOO MANY FILES'
    stop
  end if
  !-----------------------------------------------------------------------------
  ! If guard_zones is true, include the guard zones.  If time_derivs is true,
  ! write out the time derivatives also.
  !-----------------------------------------------------------------------------
  if (self%guard_zones) then
    l = self%mesh%lb
    u = self%mesh%ub
  else
    l = self%mesh%li
    u = self%mesh%ui
  end if
  if (self%time_derivs) then
    no = min(2,self%nw)
  else
    no = 1
  end if
  if (self%out_time==0.0) then
    time = self%time
  else if (self%time_derivs) then
    time = self%t(self%nt-2)
  else
    time = self%out_next
  end if
  !-----------------------------------------------------------------------------
  ! Write a small text file with info
  !-----------------------------------------------------------------------------
  if (self%no_mans_land) then
    self%ioformat = 1
  else
    self%ioformat = 2
  end if
  write (filename,'(a,i5.5,"/")') trim(io%outputname), self%iout
  call system ('mkdir -p '//trim(filename))
  write (filename,'(a,i5.5,"/",i5.5)') trim(io%outputname), self%iout, self%id
  if (io%debug(0)) then
    if (io%debug(1) .or. self%track) &
      print'(a,3x,"time:",2f13.6)', ' patch_t%output_legacy: '//trim(filename), time, self%time
  end if
  open (io%data_unit, file=trim(filename)//'.txt', form='formatted', status='unknown')
  write (io%data_unit,'(i2.2)') self%ioformat
  write (io%data_unit,'(a)') trim(self%kind)
  write (io%data_unit,'(a)') trim(self%eos)
  write (io%data_unit,'(a)') trim(self%opacity)
  write (io%data_unit,'(1p,2e18.10,6(2x,3e18.10)1x,1x,2e10.3)') time, self%dtime, &
        self%position, self%ds, self%velocity, self%llc_nat, self%llc_cart, &
        self%mesh%centre_nat, self%quality, self%gamma
  write (io%data_unit,'(6(3i10.1,2x),3i2.1,1x,1i2.1,i3)') int(self%box/self%ds+0.5), &
        u-l+1, self%li, self%ui, self%n, self%gn, self%ng, 2, self%level
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
      call self%get_q(time,buf)
      !print *, buf(:,12,12)
      !buf = buf*self%dtime*scaling%f/scaling%l**3/scaling%t
          !buf = self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it2,jo)*self%dtime*scaling%f/scaling%l**3/scaling%t!*qt &
             ! + self%mem(l(1):u(1),l(2):u(2),l(3):u(3),iv,it2,jo)*pt
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),1,it,jo,minval(buf),maxval(buf),pt
      write (io%data_unit,rec=iv+(jo-1)*2) buf
      if (io%verbose>1 .and. iv==self%i_d) print *,'rho minmax',minval(buf), maxval(buf)
    iv = 2
      buf = buf + self%src(l(1):u(1),l(2):u(2),l(3):u(3),1)
      if (io%verbose>1) &
        write (io_unit%log,'(1x,a,2i5,i3,1p,3e11.3)') &
          trim(filename),1,it,jo,minval(buf),maxval(buf),pt
      write (io%data_unit,rec=iv+(jo-1)*2) buf
      if (io%verbose>1 .and. iv==self%i_d) print *,'rho minmax',minval(buf), maxval(buf)
      deallocate (buf)
  end do
  close (io%data_unit)
   self%iout = self%iout+1
   self%out_next = (nint(self%out_next/self%out_time)+1)*self%out_time
  !$omp end critical (out_cr)
  end do
  call trace_end
END SUBROUTINE output_legacy_rt

END MODULE rt_mod
