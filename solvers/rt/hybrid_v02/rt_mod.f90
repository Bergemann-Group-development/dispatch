MODULE rt_mod
  USE bits_mod
  USE eos_mod
  USE io_mod
  USE link_mod
  USE ray_geometry_v3_mod
  USE rt_boundary_mod
  USE rt_solver_mod
  USE patch_mod
  USE task_mod
  USE timestep_mod
  USE trace_mod
  implicit none
  private

!> =====================================================================
!> intensity data type is used to store intensities in the outer shell of
!> the patch as well as quickly retrieve incoming intensities from nbors
!> =====================================================================
  type :: I_data_t
    real, dimension(:), pointer :: in  => null()                        ! (wl_bin) incoming intensity.
    real, dimension(:), pointer :: out => null()                        ! leaving intensity
    integer :: dir_ID
  end type
!> =====================================================================
!> rt_chain_t is used to make a linked list of chain rings. Each rt_task
!> owns a set of rings, related to particular direction. Each ring in 
!> the set is then connected to a downstream nbor. 
!> =====================================================================
  type, public :: chain_ring_t
    class(link_t),                 pointer :: link => null()            ! link to rt_t task
    class(chain_ring_t),           pointer :: head => null()            ! pointer to the master chain ring
    class(chain_ring_t),           pointer :: tail => null()            ! pointer to the last ring
    class(chain_ring_t),           pointer :: next => null()
    class(chain_ring_t),           pointer :: prev => null()
    class(chain_ring_t),           pointer :: temp => null()            ! temporal chain connection to do any required modifications
    class(chain_ring_t), dimension(:), pointer :: upstream => null()    ! possible upstream of different refinement level
    class(i_data_t), dimension(:), pointer :: packet_data => null()     ! quick pointing to the right locations in the data sets
    class(ray_bndl_t),             pointer :: bundle => null()          ! ray bundle of this chain ring
    logical, dimension(:), pointer :: nudged                            ! which rings have already nudged me
    integer :: direction                                                ! 1 = forward / 2 = reverse
    integer :: a_id                                                     ! location id in the ring array of a particular patch
    logical :: master = .false.                                         ! Am I the ring to rule them all?
    logical :: BC(3) = .false.                                          ! outer/inner/virtual or level
    integer :: l_chain = 1                                              ! length of this chain, points to one number, set on the head of the chain
    integer :: n_ready                                                  ! number of ready patches in this chain. When n_ready == l_chain, the chain can be executed.
    integer :: chain_nm                                                 ! ring number in the chain
  contains
    procedure :: check_masters
    procedure :: connect_to_upstream_diff_level
    procedure :: finish_chain
    procedure :: store_I_out
  end type  

!> =====================================================================
!> rt_t is the main task. 
!> every task has a set of chain rings, belonging to this patch. Every
!> item in the set is then connected to appropriate upstream/downstream
!> nbors. Also, any patch might end up being the chain head for some
!> bundle. In that case, it nudges it's own ring.
!> =====================================================================
  type, extends(patch_t), public :: rt_t
    class(chain_ring_t), dimension(:), pointer :: my_rings => null()     ! list of rings, which belongs to me
    logical, dimension(:), pointer :: updated => null()                 ! list of flags, indicates if particular direction has been updated
    type(I_data_t), dimension(:,:,:,:), pointer:: I_yz => null()        ! intensity shell data for zy plane (-/+x wall), it is primaty, that is, the largest and search starts with this one.
    type(I_data_t), dimension(:,:,:,:), pointer:: I_xz => null()        ! (x dim | z dim | n_out_dirs | lower/upper). This is the secondary one, meaning z has -2 length. Tertiary has both dims -2.
    type(I_data_t), dimension(:,:,:,:), pointer:: I_xy => null()        ! n_out_dirs is the number of directions leaving any given surface (3 dirs = 1 n_out_dir, 13 dirs = 9 n_out_dirs). This is tertiary surface.
    type(link_t), pointer             :: mhd_link                       ! link to my mhd_patch
    integer                           :: dim(3)                         ! patch dimensions
    integer                           :: n_lambda                       ! number of wavelengths
    integer                           :: n_bundle                       ! number of bundles
    real, dimension(:), pointer       :: w_lambda                       ! wavelength weights
    real, dimension(:), pointer       :: w_bundle                       ! space angle weights
    real                              :: cr_heating
    real, dimension(:,:,:),   pointer :: lnd  => null()
    real, dimension(:,:,:),   pointer :: lne  => null()
    real, dimension(:,:,:,:), pointer :: rk   => null()
    real, dimension(:,:,:,:), pointer :: src  => null()
    real, dimension(:,:,:),   pointer :: dedt => null()
    character(len=16)                 :: solver ='integral'             ! choice of solver
    procedure(solve_integral), pointer:: solve => null()
    logical :: do_output  = .false.
    logical :: data_ready = .false.
    logical :: do_chains  = .true.
  contains
    procedure :: absolute_point_position
    procedure :: boundary_conditions
    procedure :: cell2pos
    procedure :: density_weight_factor
    procedure :: dnload
    procedure :: get_boundary
    procedure :: execute_chain
    procedure :: finish
    procedure :: has_OB
    procedure :: init => init_rt
    procedure :: interpolate
    procedure :: map_intensity_shell
    procedure :: nudge
    procedure :: optical_depth_centered
    procedure :: optical_depth_staggered
    procedure :: pos2cell
    procedure :: ray_vol_map3
    procedure :: ray_vol_map4
    procedure :: ring_connects_to
    procedure :: solve_integral
    procedure :: solve_feautrier
    procedure :: update => update_rt
    procedure :: update_q  
  end type
  real, parameter:: pi4 = 16.*atan(1.)
  
CONTAINS

!> =====================================================================
!> The rt_patch initialization subroutine. The patch is responsible for
!> all the directions and their bundles. Al this is stored in chain rings.
!> FLOW:
!> 1) read the namelist
!> 2) initialize the patch, set bits
!> 3) set 
!> =====================================================================
SUBROUTINE init_rt (self)
  class (rt_t) :: self
  class (ray_bndl_t), pointer :: ray_bnd
  class (ray_dirs_t), pointer :: ray_dir
  class (packet_t), pointer      :: packet
  class (chain_ring_t), pointer :: ring
  logical, save :: first_time= .true.
  logical, save :: store_rt  = .false.
  real(8), save :: dtime=0.01_8, grace=0.01_8, init_time=-999d0
  real,    save :: cr_heating=0.0
  character(len=16) :: solver = 'integral'
  integer :: a_id, dd, n_out
  integer:: ignore, ibnd
  namelist /rt_params/ init_time, dtime, grace, solver, store_rt, cr_heating
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
  
  self%do_output = store_rt
  self%nv = 1                                                           ! time slots for dedt
  self%nt = 1                                                           ! only one time slot FIXME 
  self%nw = 1                                                           ! no time deriv
  
  call self%patch_t%init                                                ! initialize patch
  call self%set (bits%static+bits%no_download)                          ! special bits
  
  self%i_e = 1                                                          ! dedt slot
  self%dedt => self%mem(:,:,:,self%i_e,self%it,1)                       ! map dedt to mem
  
  self%n_lambda = eos%n_lambda                                          ! number of wavelengths
  self%n_bundle = ray_geometry%head%total_n_bundles * 2                 ! double n_bundles since there are 2 directions
  self%cr_heating = cr_heating                                          ! store cr heating, if it is present
  self%do_chains = ray_geometry%head%do_chains                          ! are there chains in the system, or patches are effectively isolated
  
  
  allocate (self%w_lambda(self%n_lambda))                               ! frequency bin weights
  allocate (self%w_bundle(self%n_bundle))                               ! space angle weights
  self%w_lambda = eos%w_lambda
  self%w_bundle = pi4/2./self%n_bundle                                  ! space angle weights FIXME: could be different weights for different directions
  
  self%dtime = dtime                                                    ! 0.1 second default
  self%time = init_time                                                 ! -999.0 default
  self%dt_fixed = dtime
  self%grace = grace                                                    ! grace on by default, but very small
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
  self%updated = .false.                                                ! none are updated initially
  allocate(self%my_rings(self%n_bundle))                                ! allocate all the rings, this patch contains
  
  ! temporary workaround to set the number of leaving rays per wall. for 3
  ! directions it is 1, for 13 directions (45 degree space-angle separation)
  ! it is 9.
  if (ray_geometry%head%n_ray_dirs == 3) then
    n_out = 1
  else if (ray_geometry%head%n_ray_dirs == 13) then
    n_out = 9
  else
    error stop "init_rt (1) making Intensity shell problem: ray_geometry should give n_out"
  end if

  
  ! --- Allocate the intensity outer shell ---
  allocate(self%I_yz(self%dim(2), self%dim(3), n_out, 2))
  allocate(self%I_xz(self%dim(3)-2, self%dim(1), n_out, 2))
  allocate(self%I_xy(self%dim(1)-2, self%dim(2)-2, n_out, 2))
  ! FIXME: options for collapsed dimensions, small patches, etc.
  
  ! -- set up all the rings
  a_id = 1                                                              ! ring personal ID in the ring array
  ray_dir => ray_geometry%head%ray_list%head                            ! pick a ray direction
  do while(associated(ray_dir))                                         ! loop over all ray directions
    do ibnd =1, ray_dir%n_bundles                                       ! loop over ray bundles
      ray_bnd => ray_dir%bundle(ibnd)                        
      do dd=1,2                                                         ! separate rings for forward/reverse directions
        self%my_rings(a_id)%link => self%link                           ! rings should know, what rt patch hosts them
        self%my_rings(a_id)%direction = dd                              ! forward/reverse direction rings
        self%my_rings(a_id)%a_id = a_id                                 ! ring personal ID in the ring array
        self%my_rings(a_id)%bundle => ray_bnd                           ! ring has to know, which bundle it corresponds to
        allocate(self%my_rings(a_id)%packet_data(ray_bnd%n_rays))       ! allocate intensity in/out data for the ring. Array consists of all the rays
        call self%map_intensity_shell(self%my_rings(a_id))              ! map the rings to the intensity outer shell (their "out" position in packet_data)
        a_id = a_id + 1                                                 ! next ID
      end do
    end do
    ray_dir => ray_dir%next                                             ! next direction
  end do
  nullify(ray_dir, ray_bnd)
  call trace_end  
END SUBROUTINE init_rt   

!> =====================================================================
!> =====================================================================
SUBROUTINE dnload (self)
  class (rt_t) :: self
END SUBROUTINE dnload

!> =====================================================================
!> This subroutine goes though all the rays in the bundle, belonging to
!> particular chain ring and maps all the rays to particular cell in the
!> outer shell in the patch. This will completely remove the intensity
!> lookup time, unless the patches are moving.
!> INPUT: rt_patch as self and the ring to be mapped.
!> OUTPUT: none
!> RESULT: the ring knows where to put output intensities on the outer 
!>   shell of the rt patch.
!> FLOW: 
!> 1) the bundle lwall and uwall are checked to find, which faces the
!>   bundle cross. The -/+x face is the primary one, consisting of y*z
!>   cells, the -/+y face is a secondary one, consisting of x*(z-2) cells
!>   and the -/+z face is the tertiary one, consisting of (x-2)*(y-2) cells.
!>   the last indice in the shells is the lower(1)/upper(2) wall. Indices
!>   a and b are set to correspond to the right axes to obtain the cell
!>   coordinates later.
!> 2) a loop over all packets in the bundle is done.
!> 3) if direction is forward, one needs the coordinates of the last point
!>   in the ray, if it is reverse, one needs the coordinates for the first
!>   point.
!> 4) cell coordinate is taken for particular ray
!> 5) a loop is initialized to find unocuppied data entry point (different
!>   directions use same cells and nobody knows, which direction will do
!>   this subroutine first). As soon as available entry point is found,
!>   the quick lookup array point, related to particular ray, points to
!>   that location. Work for this ray is done, moving to next ray or packet.
!> =====================================================================
SUBROUTINE map_intensity_shell(self, ring)
  class(rt_t) :: self
  class(chain_ring_t) :: ring
  !.....................................................................
  class(packet_t), pointer :: packet
  class(I_data_t), dimension(:,:,:), pointer:: wall
  class(I_data_t), pointer :: this
  integer :: ray, rid, c(3), pt, dir, a,b, i, ipck
  logical :: selected
  !---------------------------------------------------------------------
  dir = ring%direction                                                  ! forward of reverse
    
  ! --- select the right wall and coordinate axis ---
  if ((ring%bundle%lwall(1) /= 0).and.(dir ==2)) then
    wall => self%I_yz(:,:,:, 1)
    a = 2;  b = 3
  else if ((ring%bundle%uwall(1) /= 0).and.(dir ==1)) then
    wall => self%I_yz(:,:,:, 2)
    a = 2;  b = 3
  else if ((ring%bundle%lwall(2) /= 0).and.(dir ==2)) then
    wall => self%I_xz(:,:,:, 1)
    a = 3;  b = 1
  else if ((ring%bundle%uwall(2) /= 0).and.(dir ==1)) then
    wall => self%I_xz(:,:,:, 2)
    a = 3;  b = 1
  else if ((ring%bundle%lwall(3) /= 0).and.(dir ==2)) then
    wall => self%I_xy(:,:,:, 1)
    a = 1;  b = 2
  else if ((ring%bundle%uwall(3) /= 0).and.(dir ==1)) then
    wall => self%I_xy(:,:,:, 2)
    a = 1;  b = 2
  else
    error stop "map_intensity_shell:: something went wrong..."
  end if
  
  rid=1
  do ipck =1, ring%bundle%n_packets
    packet => ring%bundle%packet(ipck)                                  ! loop over ray packets in the bundle                                        
    if (dir ==1) then
       pt = packet%n_points
    else
      pt = 1
    end if
    do ray = 1, packet%n_rays                                           ! loop over rays in the packet
      selected = .false.
      c = packet%coords(:,pt,ray)
      i=1
      do while (.not.selected)
        this => wall(c(a),c(b),i)
        if (associated(this%out)) then                                       ! allocate wall data.
          i = i+1
        else
          wall(c(a),c(b),i)%dir_ID = ring%bundle%dir_ID                 ! know, which direction this intensity is from
          allocate(wall(c(a),c(b),i)%out(self%n_lambda))                ! allocate the outgoing intensity storage
          ring%packet_data(rid)%out => wall(c(a),c(b),i)%out            ! point the fast lookup to the shell cell position
          selected = .true.
        end if
      end do
      rid=rid+1
    end do
  end do
END SUBROUTINE map_intensity_shell


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
  !---------------------------------------------------------------------
  if (self%data_ready) return
  
  call trace%begin ('rmhd_t%update', itimer=itimer)
  if (io%verbose > 1) print *, 'RADIATIVE TRANSFER TIME', self%time
  if (self%time < 0d0) then
    self%time  = 0d0                                                    ! set current time, if it was negative (to do the first RT run before MHD updates)
  end if
  self%dedt = 0.0                                                       ! reset net cooling
  task => self%mhd_link%task                                            ! pick the rmhd task
  select type (task)
  class is (patch_t)
    mhd => task
  end select
  ! --- allocate all the scratch arrays ---
  allocate (self%lnd (self%dim(1), self%dim(2), self%dim(3)))           ! log density scratch
  allocate (self%lne (self%dim(1), self%dim(2), self%dim(3)))           ! log internal energy scratch
  allocate (self%rk  (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  allocate (self%src (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  
  call mhd%log_interpolate (self%time, mhd%i_d, self%lnd)               ! interpolate in time and space to get densoties
  call mhd%log_interpolate (self%time, mhd%i_e, self%lne)               ! interpolate in time and space to get internal energies
  if (io%verbose > 1) then
    print *, self%id, 'lnd', minval(self%lnd), maxval(self%lnd)
    print *, self%id, 'lne', minval(self%lne), maxval(self%lne)
  end if
  self%lne = self%lne - self%lnd                                        ! internal energy per unit mass
  
  call eos%lookup (shape(self%lne), lne=self%lne, lnd=self%lnd, &       ! call the table lookup, ...
                   rk=self%rk, src=self%src)                            ! ... results go into rt_t    
  self%data_ready =.true.                                               ! flag, that internal data is ready
  call trace%end (itimer)
  call self%nudge                                                       ! nudge the chain heads
  call self%courant_condition
  call timestep%update (self%id, self%iit, self%dt, self%mem, self%mesh)
END SUBROUTINE update_rt 

!> =====================================================================
!> Nudge the chain heads. If all patches in the chain become ready, then
!> solve for the net heating rate, storing it in a slot that rotates
!> with time for each patch in the chain.
!> =====================================================================
SUBROUTINE nudge (self)
  class(rt_t) :: self
  class(chain_ring_t), pointer :: ring
  class(task_t), pointer :: tmp
  integer :: r_id
  !---------------------------------------------------------------------
  do r_id = 1, self%n_bundle
    ring => self%my_rings(r_id)
    !$omp critical (chain_cr)
    if (io%verbose > 1) then
      print *, 'ring no.',r_id
      print *, 'before nudge, n_ready',ring%head%n_ready, ring%head%l_chain
      print *,'nudged by master', ring%master,'rt time',self%time
      print *, 'other rings updated',self%updated
      print *, 'my data is ready',self%data_ready
      print *, 'my ID is:',self%id,'ring', ring%a_id
    end if
    
    tmp => ring%head%link%task
    select type(tmp)
    class is (rt_t)
    if (io%verbose > 1) print *, 'nudging task ID:', tmp%id,'ring id', ring%head%a_id, &
    'nudged allocated:',associated(ring%head%nudged),'size'
    if (io%verbose > 1) print *, size(ring%head%nudged),'its content:',ring%head%nudged, 'head assoc', associated(ring%head)
    if (io%verbose > 1) print *, 'l_chain',ring%head%l_chain
    end select
    
    ring%head%n_ready = ring%head%n_ready+1
    ring%head%nudged(ring%chain_nm) = .true.
    if (ring%head%n_ready > ring%head%l_chain) print *, "WARNING: l_chain exceeded in nudges.....!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    
    !if (ring%head%n_ready > ring%head%l_chain) call self%execute_chain(ring%head)
    if (all(ring%head%nudged)) call self%execute_chain(ring%head)

    !$omp end critical (chain_cr) 
  end do
END SUBROUTINE nudge

!> =====================================================================
!> Loop through all rings in the chain and solve the actual rt
!> =====================================================================
SUBROUTINE execute_chain (self, head)
  class(rt_t) :: self
  class(rt_t), pointer :: rt_patch
  class(task_t), pointer :: task
  class(chain_ring_t), pointer :: head, ring, prev
  real, dimension(:,:), pointer :: i_in
  integer :: ray
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  !call trace%begin ('rmhd_t%execute_chain', itimer=itimer)
  head%n_ready = 0
  ring => head
  nullify(prev)
  do while (associated(ring))
    task => ring%link%task                                              ! get the patch of this chain ring
    select type (task)
    class is (rt_t)
    rt_patch => task
    end select

    allocate(i_in(size(ring%packet_data),rt_patch%n_lambda))
    if(associated(prev)) then
      do ray=1, ring%bundle%n_rays
        I_in(ray,:) = ring%packet_data(ray)%in
      end do
    else
      call rt_patch%get_boundary(ring,i_in)
    end if    
    call rt_patch%solve(ring, i_in) !optionally, add second batch of intensities to do 2 dirs at once

    deallocate(i_in)
    rt_patch%updated(ring%a_id) = .true.
    if(all(rt_patch%updated)) call rt_patch%finish
    prev => ring
    ring => ring%next
  end do
 ! print *, 'CHAIN END'
  head%nudged = .false.
  nullify(rt_patch, task, prev)  
  !call trace%end(itimer)
END SUBROUTINE execute_chain

!> =====================================================================
!> Boundary data is on inner/outer boundary, different refinement level,
!> or on a virtual patch. Lookup through the use of cell2pos & pos2cell
!> subroutines. Slightly more costly, but necessary. Virtual patches only
!> store wall intensities per direction, whereas for different refinement
!> levels bundle chains are effectively broken. Outer and inner boundaries
!> have to modify the ray lenghs so that RT would not be solved beyond them.
!> this setting could be done only once, unless patches or boundaries are
!> moving or AMR is on.
!> =====================================================================
SUBROUTINE get_boundary (self, ring, I_in)
  class(rt_t) :: self
  class(chain_ring_t), pointer :: ring, upstream
  real, dimension(:,:), pointer :: I_in
  !---------------------------------------------------------------------
  class(packet_t), pointer :: packet
  integer :: pck, ray
  real(kind=8) :: ppos(3)
  !---------------------------------------------------------------------
  if (ring%BC(1)) then
    do pck =1, ring%bundle%n_packets
      packet => ring%bundle%packet(pck)
      do ray = 1, packet%n_rays
        ! for each ray index in the packet, find a point position
        ppos = self%absolute_point_position (ring, packet, ray, 1, ring%direction)
        call rt_boundary%set_obc (ppos, I_in(ray,:))
      end do
    end do
  end if

  ! inner boundary should modify the ray lengths!!
  if (ring%BC(2)) then
    error stop "inner BC for RT is not fully functional yet"
 !   packet => ring%bundle%rays
 !   do pck =1, ring%bundle%n_packets
 !     do ray = 1, packet%n_rays
        ! for each ray index in the packet, find a point position
 !       ppos = self%absolute_point_position (ring, packet, ray, 1, ring%direction)
 !       call rt_boundary%set_ibc (ppos, I_in(ray,:))
 !     end do
 !     packet => packet%next
 !   end do
  end if
  
  if (ring%BC(3)) then
      error stop "getting boundary from virtual patches or different level patches is not implemented yet"
  end if
  
END SUBROUTINE get_boundary

!> =====================================================================
!> Check if the rays in self match up in the downstream with the rays in other
!> INPUT: self and nbor are rt_patches. mine is my ring i am testing
!>   against nbor's ring other. logical connects - true if chain can be 
!>   made, logical is_upstream - if the ring is in upstream from me, but
!>   is on a different refinement level, it could still potentially
!>   give me the intensities i need.
!> RESULT: two logicals: if can be connected and if one is upstream of 
!> another, but on different refinement level.
!> SUBROUTINES CALLED FROM WITHIN: 
!> rt_mod_t::rt_t%cell2pos
!> rt_mod_t::rt_t%pos2cell
!> FLOW:
!> 1) both logicals are set to false by default. 
!> 2) both rings are checked whether they correspond to same directions and
!>   direction IDs (spatial direction).
!> 3) rt tasks, corresponding to these rings, are checked for different
!>   refinement level. If that is the case, appropriate flag is set.
!> 4) otherwise the connections are tested. Only one ray in a packet is
!>   enough to see, whether connection is good - ray_geometry module makes
!>   sure all packets are fully interconnected. 
!> 5) the last point is taken from other ring and the upstream position
!>   is taken from me. if the points coincide via cell2pos and pos2cell,
!>   there is a connection.
!> =====================================================================
SUBROUTINE ring_connects_to (self, nbor, mine, other, connects, is_upstream)
  class (rt_t) :: self
  class (rt_t), pointer :: nbor
  class (chain_ring_t) :: mine
  class (chain_ring_t) :: other
  logical :: connects, is_upstream
  !.....................................................................
  integer :: pp, ray, first_pt(3), rev(3), rid, pid
  logical :: tout(3)
  !---------------------------------------------------------------------
  connects = .false.
  is_upstream = .false.
  tout=.false.
  
  if (mine%bundle%dir_ID /= other%bundle%dir_ID)     return             ! different ray spatial direction
  if (mine%direction /= other%direction)             return             ! different direction
  
  select case(mine%direction)
  case(1)
    if (mine%bundle%b_id_fw  /= other%bundle%ID) return                 ! wrong bundle to make a chain
  case(2)
    if (mine%bundle%b_id_rev /= other%bundle%ID) return                 ! wrong bundle to make a chain
  end select
  
  if (self%level /= nbor%level) then
    if (abs(self%level-nbor%level) < 2) then
      ! test if is_upstream = .true.
      !if different refinement level, do cell2pos, pos2cell to check the connection possibility
      error stop "Different refinement level connections not done yet"
    end if
  else
    ray = 1                                                             ! testing only one ray in the bundle
    select case(mine%direction)                                           
    case(1)
      first_pt =  mine%bundle%packet(1)%coords(:,1,ray)                 ! taking the coordinates of the first point in the ray 
      rid = mine%bundle%packet(1)%r_id_rev(ray)                         ! ray id in a packet list of a ray in the previous bundle 
      pid = mine%bundle%packet(1)%p_id_rev(ray)                         ! packet id in the previous bundle 
      rev = other%bundle%packet(pid)%fwd_coords(:,rid)                  ! zeroth point in the ray (location of incoming intensity)                                     
    case(2)
      pp = mine%bundle%packet(1)%n_points
      first_pt =  mine%bundle%packet(1)%coords(:,pp,ray)                ! taking the coordinates of the last point in the ray 
      rid = mine%bundle%packet(1)%r_id_fw(ray)                          ! ray id in a packet list of a ray in the previous bundle 
      pid = mine%bundle%packet(1)%p_id_fw(ray)                          ! packet id in the previous bundle 
      rev = other%bundle%packet(pid)%rev_coords(:,rid)                  ! zeroth point in the ray (location of incoming intensity)                                                              
    end select
  
    first_pt = nbor%pos2cell(self%myposition(first_pt))
    do pp = 1,3
      if (first_pt(pp) == rev(pp)) tout(pp) = .true.
    end do
    
    if (all(tout)) connects = .true.
  end if  
END SUBROUTINE ring_connects_to

SUBROUTINE connect_to_upstream_diff_level(self, ring)
  class(chain_ring_t) :: self
  class(chain_ring_t) :: ring
  
  error stop "rt_t:: connect_to_upstream_diff_level: not implemented"
END SUBROUTINE connect_to_upstream_diff_level

!> =====================================================================
!> Connect the rings to the appropriate incoming intensities from the
!> right upstreams. 
!> INPUT: a temporary ring, which has all the master rings in its chain
!> RESULT: all the rings in all the chains know, where to get incoming
!>   intensities from.
!> FLOW: 
!> 1) first item is ommitted, since it is just a placeholder. All the 
!>   other items are oranised in "reverse", since master rings do/should
!>   not contain a %prev member.
!> 2) a loop over all master rings is made.
!> 3) a loop over all rings in this chain, belonging to a particular
!>   master is done
!> 4) master gets information about the chain length. This is utterly
!>   important later.
!> 5) ring direction is checked and based on that the right ray arrangement
!>   in upstream ring is determined
!> 6) depending on what is the periodicity of chain alterations for this
!>   orientation, one determines, where the rays are (a supplement for (5))
!> 7) the appropriate ray in the ring is connected to the outer shell of
!>   the upstream patch.
SUBROUTINE finish_chain(self)
  class(chain_ring_t) :: self
  !.....................................................................
  class(chain_ring_t), pointer :: this, prev
  class(chain_ring_t), pointer :: ring, pring
  class(task_t), pointer :: task
  class(rt_t), pointer :: rt_patch
  class(I_data_t), dimension(:,:,:), pointer:: wall
  integer :: ray, cn
  integer, dimension(:), pointer :: nb_ray_pos
  !---------------------------------------------------------------------
  this => self%prev                                                     ! first item is not needed
  do while(associated(this))                                            ! loop over master rings
    prev => this%prev                                                   ! this temporary chain is organised in reverse
    this%head => this                                                   ! link the master chain to itself (needed for later)
    this%chain_nm = 1                                                   ! i am first in the chain
   ! nullify(this%prev)                                                  ! nullify unnecessary pointer to %prev
    ring => this%next                                                   ! pick the first slave ring in this chain
    pring => this                                                       ! the %prev for that first slave ring
    cn = 2                                                              ! incremental chain number
    do while(associated(ring))                                          ! go through all the slave rings
      this%l_chain = this%l_chain + 1                                   ! increment the chain length for the master to know
      ring%chain_nm = cn
      cn = cn+1
      error stop "modify finish chain"
 !     select case(ring%direction)                                       ! select the right direction of this ring
  !    case(1)
 !       nb_ray_pos => ring%bundle%ray_chain_rev                         ! forward = upstream is in my reverse
   !   case(2)
  !      nb_ray_pos => ring%bundle%ray_chain_fw                          ! reverse = upstream is in front of me
  !    end select
   
  !    select case(ring%bundle%l_chain)                                  ! figure out the chain packet periodicity
  !    case(:1)                                                          ! if it is 0 or 1
  !      do ray=1, ring%bundle%n_rays                                      
  !        ring%packet_data(ray)%in => pring%packet_data(ray)%out        ! rays are in the same position (no alteration in periodicity)
  !      end do
  !    case default
  !      do ray=1, ring%bundle%n_rays
  !        ring%packet_data(ray)%in => pring%packet_data(nb_ray_pos(ray))%out
  !      end do
  !    end select
     
      pring => ring
      ring => ring%next
    end do
    if (this%l_chain <1) error stop "WTF?"
    allocate(this%nudged(this%l_chain))
    this%nudged = .false.
    this => prev                                                        ! point to next master ring
  end do
  if (io%verbose>1) call self%check_masters
END SUBROUTINE finish_chain

!> =====================================================================
!> to run this, disable the deallocation in finish_chain subroutine and
!> uncomment the last line there
!> =====================================================================
SUBROUTINE check_masters(self)
  class(chain_ring_t) :: self
  class(chain_ring_t), pointer :: ring
  !---------------------------------------------------------------------
  ring => self%prev
  do while(associated(ring))
  
  print *,'ring_l_chain', ring%l_chain
  print *,'associated nudged',associated(ring%nudged)
  print *,'nudged size',size(ring%nudged)
  
  ring => ring%prev
  end do
  error stop "end"
END SUBROUTINE check_masters

!> =====================================================================
!> Find the absolute position of index 'ray' in ray packet 'rays'.
!> If the optional bn parameter is present, the answer depends on
!> the direction, otherwise 'pt' chooses the point along the ray.
!> =====================================================================
FUNCTION absolute_point_position (self, ring, packet, ray, pt, bn) result (cell_pos)
  class (rt_t) :: self
  class (chain_ring_t), pointer :: ring
  class (packet_t), pointer :: packet
  integer, optional :: bn
  integer :: ray, pt, cp(3)
  real(kind=8) :: cell_pos(3)
!-----------------------------------------------------------------------
  call trace_begin ('rt_t%absolute_point_position', 4)
  if (present(bn)) then
    select case(ring%direction)
    case(1)
      ! self%my_rays
      cp = packet%fwd_coords(:,ray)
    case(2)
      cp = packet%rev_coords(:,ray)
    end select
  else
    cp = packet%coords(:,pt,ray)
  end if
  cell_pos = self%myposition(cp)
  call trace_end
END FUNCTION absolute_point_position

LOGICAL FUNCTION has_OB (self, ring)
  class (rt_t) :: self
  class (chain_ring_t) :: ring
  integer :: i(3), d(3),t
  real(8) :: p(3)
  !.....................................................................
  has_OB = .false.
  d = ring%bundle%dstep
  if (ring%direction == 1) then
    i = ring%bundle%packet(1)%rev_coords(:,1)
  else 
    i = ring%bundle%packet(1)%fwd_coords(:,1)
  end if
    p = self%myposition(i)
  
  do t = 1,3
    if (ring%direction ==1) then
      if (d(t)>0) then
        if (p(t)<rt_boundary%lbox(t)) has_OB = .true.
      else
        if (p(t)>rt_boundary%ubox(t)) has_OB = .true.
      end if
    else
      if (d(t)>0) then
        if (p(t)>rt_boundary%ubox(t)) has_OB = .true.
      else
        if (p(t)<rt_boundary%lbox(t)) has_OB = .true.
      end if
    end if
  end do

END FUNCTION has_OB

LOGICAL FUNCTION has_IB (self, ring)
  class (rt_t) :: self
  class (chain_ring_t) :: ring
  integer :: i(3), d(3),t
  real(8) :: p(3)
  !.....................................................................
  has_IB = .false.
  d = ring%bundle%dstep
  if (ring%direction == 1) then
    i = ring%bundle%packet(1)%fwd_coords(:,1)
    p = self%cell2pos(i)
  else 
    i = ring%bundle%packet(1)%rev_coords(:,1)
    p = self%cell2pos(i)
  end if
  
  do t = 1,3
    if (d(t)>0) then
      if (p(t)<rt_boundary%lbox(t)) has_IB = .true.
    else
      if (p(t)>rt_boundary%ubox(t)) has_IB = .true.
    end if
  end do
END FUNCTION has_IB
  
SUBROUTINE finish(self)
  class(rt_t) :: self  
  !---------------------------------------------------------------------
  self%updated=.false.
  self%data_ready = .false.
  call self%boundary_conditions
  call self%output
 ! print *, 'Finished rt patch'
  deallocate (self%lnd, self%lne, self%rk, self%src)                                       ! deallocate scratch
 ! print *, 'new time, dtime',self%time, self%dtime
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
SUBROUTINE solve_integral (self, ring, i_in, i_in2)
  class (rt_t) :: self
  class(chain_ring_t), pointer  :: ring
  real, dimension(:,:), pointer :: i_in
  real, dimension(:,:), pointer, optional :: i_in2
  !.....................................................................
  class(packet_t),        pointer :: packet
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  real, dimension(:,:),   pointer :: I_out, dwf, I_in0
  integer                         :: nrays, npt, nwb, ipck
  real                            :: ds
  logical :: one_p
  integer, save                   :: itimer=0
  class(task_t), pointer :: tmp
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve_integral', itimer=itimer)

  nwb   =  self%n_lambda
  do ipck = 1, ring%bundle%n_packets
    packet => ring%bundle%packet(ipck)
    nrays =  packet%n_rays

    if (packet%n_points < 3) then
      one_p = .true.
      npt = packet%n_points + 2
    else
      npt = packet%n_points
      one_p = .false.
    end if
    ds    =  packet%ds
    allocate(s   (nrays, nwb, npt))            ! TODO: boundaries modify npt
    allocate(rk  (nrays, nwb, npt))            ! TODO: boundaries modify npt
    allocate(dtau(nrays, nwb, npt))
    allocate(exp0(nrays, nwb, npt))
    allocate(q   (nrays, nwb, npt))
    allocate(dwf (nrays,      npt))
    allocate(I_out(nrays, nwb))
    allocate(I_in0(nrays, nwb))
   ! print *,'ring is master',ring%master, self%nt
   ! print *,'in solve integral, data ready, my id', self%data_ready, self%id
   ! if (associated(ring%prev)) then
   ! tmp=> ring%prev%link%task
   ! select type(tmp)
   ! class is (rt_t)
   ! print *, 'nudged by prev',tmp%time,tmp%updated, tmp%dtime, self%dtime, tmp%id
   ! end select
   ! end if
   ! if (associated(ring%next)) then
   ! tmp=> ring%next%link%task
   ! select type(tmp)
   ! class is (rt_t)
   ! print *, 'nudged by next',tmp%time,tmp%updated, tmp%dtime, self%dtime, tmp%id
   ! end select
   ! end if
   ! print *, 'time, updated, ring id',self%time, self%updated, ring%a_id
    I_in0 = I_in(packet%bnd_loc:packet%bnd_loc+packet%n_rays-1,:)
   ! print *,'self dims',self%dim(1), self%dim(2), self%dim(3)
   ! print *,'getting s'
    call self%ray_vol_map4(packet, self%src, s)                         ! map source function to rays
   ! print *,'getting rk'
    call self%ray_vol_map4(packet, self%rk, rk)                         ! map opacities to rays
    call self%optical_depth_centered(ds, rk, dtau)                      ! optical depth is centered for integral formulation
    call self%ray_vol_map3(packet, self%lnd, dwf)                       ! density weight factor
    dwf = dwf * ds  
    if (io%verbose>1) &
        print *,'dtau',self%id, minval(dtau), maxval(dtau)

#if defined __INTEL_COMPILER
      exp0 = 2**(-dtau*1.44269504) 
#else
      exp0 = exp(-dtau)                                                 ! expopential factors
#endif        
    
    ! TODO: For now only one direction per update. In no_chains mode, both
    ! directions could be done simultaneously

    if (ring%direction ==1) then
      call rt_solver%forward (dtau, exp0, s, I_in0, q, I_out, one_p)   ! solve forward
    else
      call rt_solver%reverse (dtau, exp0, s, I_in0, q, I_out, one_p)   ! solve reverse
    end if
    
   ! print *,'updating q'
    call self%update_q(packet, q, dwf, ring%a_id, one_p)
    call ring%store_I_out(packet, I_out)
    
    deallocate (s, rk, dtau, exp0, q, dwf, I_out, I_in0)   
  end do
  call trace_end (itimer)
END SUBROUTINE solve_integral

SUBROUTINE store_I_out(self, packet, I_out)
  class(chain_ring_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:), pointer :: I_out
  integer :: ray, i
  !.....................................................................
  
  !--------------------------------------------------------------------
  i=1
  !print *, associated(self%packet_data), size(self%packet_data)
  !print *, packet%bnd_loc, packet%bnd_loc+packet%n_rays-1
  do ray=packet%bnd_loc, packet%bnd_loc+packet%n_rays-1
    self%packet_data(ray)%out = I_out(i,:)
    i=i+1
  end do

END SUBROUTINE store_I_out

!> =====================================================================
!> Same as solve_integral, but for the Feautrier formulation
!> =====================================================================
SUBROUTINE solve_feautrier (self, ring, i_in, i_in2)
  class (rt_t) :: self
  class (chain_ring_t), pointer :: ring
  real, dimension(:,:), pointer :: i_in
  real, dimension(:,:), pointer, optional :: i_in2
  !.....................................................................
  class(packet_t),        pointer :: packet
  integer                         :: nrays, npt, nwb, ipck
  real :: ds
  real, dimension(:,:), pointer :: i_out1, i_out2, dwf
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  integer, save                   :: itimer=0
  logical :: one_p
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve_feautrier', itimer=itimer)
  
  if (.not.present(i_in2)) error stop "solve Feautrier: needs 2 directions to act"
  error stop "Feautrier formulation not set up properly yet"
  do ipck=1, ring%bundle%n_packets
    packet => ring%bundle%packet(ipck)
    nrays =  packet%n_rays
    npt   =  packet%n_points
    ds    =  packet%ds
    allocate(s   (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(rk  (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(dtau(nrays, nwb, npt+2))
    allocate(exp0(nrays, nwb, npt+2))
    allocate(q   (nrays, nwb, npt+2))
    allocate(dwf (nrays,      npt+2))
    allocate(I_out1(nrays, nwb))
    allocate(I_out2(nrays, nwb))
    call self%ray_vol_map4(packet, self%src, s)                         ! map source function to rays
    call self%ray_vol_map4(packet, self%rk, rk)                         ! map opacities to rays
    call self%optical_depth_staggered(ds, rk, dtau)                     ! optical depth is staggered for Feautrier formulation
    call self%ray_vol_map3(packet, self%lnd, dwf)                       ! density weight factor
    dwf = dwf * ds  
    if (io%verbose>1) &
        print *,'dtau',self%id, minval(dtau), maxval(dtau)
    
    call rt_solver%transfer (dtau, s, I_in, I_in2, q, I_out1, I_out2)        ! solve bi-directional
    call self%update_q(packet, q, dwf, ring%a_id, one_p)
    call ring%store_I_out(packet, I_out1)
    call ring%store_I_out(packet, I_out2)
    
    deallocate (s, rk, dtau, exp0, q, dwf, I_out1, I_out2)
  end do

  call trace_end (itimer)
END SUBROUTINE solve_feautrier

!> =====================================================================
!> Map cell data onto rays
!> =====================================================================
SUBROUTINE ray_vol_map4 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  integer :: p(3), ray, pt, nwb, wl
  real, dimension(:,:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:,:) , pointer, intent(in)  :: in
  integer, dimension(:,:,:), pointer :: cell
  !---------------------------------------------------------------------
  cell => packet%coords
  !print *, 'present out,in',associated(out), associated(in)
  !print *, 'shapes',shape(out), shape(in)
  nwb = size(in,4)
    do ray = 1, size(out,1)
      do pt = 2, size(out,3)-1                                          ! First and last points ommitted
        p = cell(:,pt-1,ray)
        !DIR$ SIMD
        do wl=1,nwb
          out(ray,wl,pt) = in(p(1),p(2),p(3),wl)
        end do
      end do
      !DIR$ SIMD
      do wl=1,nwb
        out(ray,wl,1) = in(packet%fwd_coords(1,ray), &
        packet%fwd_coords(2,ray), packet%fwd_coords(3,ray), wl)         ! Get the first point
        out(ray,wl,size(out,3)) = in(packet%rev_coords(1,ray), &
        packet%rev_coords(2,ray), packet%rev_coords(3,ray), wl)         ! Get the last point
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
!> Map cell data onto rays
!> =====================================================================
SUBROUTINE ray_vol_map3 (self, packet, in, out)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  integer :: p(3), ray, pt
  real, dimension(:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:) , pointer, intent(in)  :: in
  integer, dimension(:,:,:), pointer :: cell
  !---------------------------------------------------------------------
  cell => packet%coords
    do ray = 1, size(out,1)
      do pt = 2, size(out,2)-1                                          ! First and last points ommitted
        p = cell(:,pt-1,ray)
        !DIR$ SIMD
          out(ray,pt) = in(p(1),p(2),p(3))
      end do
      !DIR$ SIMD
        out(ray,1) = in(packet%fwd_coords(1,ray), &
        packet%fwd_coords(2,ray), packet%fwd_coords(3,ray))             ! Get the first point
        out(ray,size(out,2)) = in(packet%rev_coords(1,ray), &
        packet%rev_coords(2,ray), packet%rev_coords(3,ray))             ! Get the last point
    end do
    
  if (io%verbose > 1) then
   ! if (any(in==0)) then
   !   print *, in
   !   error stop "ray vol map3 in=0"
   ! end if
   ! if (any(out==0)) then
   !   do pt = 1, size(out,2)
   !   print *,"pt",pt, out(:,pt)
   !   end do
   !   error stop "ray vol map3 out=0"
   ! end if
  end if
  
END SUBROUTINE ray_vol_map3

SUBROUTINE update_q (self, packet, q, dwf, ib, one_p)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:,:), pointer :: q
  real, dimension(:,:), pointer :: dwf
  real, dimension(:,:,:), pointer :: qq
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), wl, ray, pt, ib, nr, np
  integer, save:: itimer=0
  logical :: one_p
  !---------------------------------------------------------------------
  call trace_begin ('rt_t%update_q', 4, itimer=itimer)
  qq => self%dedt(:,:,:)
  cell => packet%coords
  nr = packet%n_rays
  np = packet%n_points
  
  if (one_p) then
    do ray = 1, nr
      do pt = 1, np-1
        p = cell(:,pt,ray)
        do wl =1, size(q,2)
          qq(p(1),p(2),p(3)) = &
          qq(p(1),p(2),p(3)) + q(ray,wl,pt+1)* self%w_lambda(wl) * &
           self%w_bundle(ib) * dwf(ray,pt)
        end do
      end do
    end do
  else
    do ray = 1, nr
      do pt = 1, np
        p = cell(:,pt,ray)
        do wl =1, size(q,2)
          qq(p(1),p(2),p(3)) = &
          qq(p(1),p(2),p(3)) + q(ray,wl,pt)* self%w_lambda(wl) * &
           self%w_bundle(ib) * dwf(ray,pt)
        end do
      end do
    end do
  end if
  call trace_end (itimer)
END SUBROUTINE update_q

!> =====================================================================
!> Compute the optical depth along rays
!> =====================================================================
SUBROUTINE optical_depth_centered (self, ds, rk, dtau)
  class (rt_t) :: self
  real :: ds
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
  real :: ds
  real, dimension(:,:,:), pointer :: rk
  real, dimension(:,:,:), pointer :: dtau,dt
  !.....................................................................
  integer                         :: wl, pt
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
  real :: ds
  real, dimension(:,:,:), pointer :: dwf
  !.....................................................................
  integer                         :: iy, iz
  !---------------------------------------------------------------------
  do iz=1,size(dwf,3)
    do iy=1,size(dwf,2)
      dwf(:,iy,iz) = self%lnd(:,iy,iz) * ds
    end do
  end do
END SUBROUTINE density_weight_factor

!> =====================================================================
!> Compute the optical depth along rays
!> =====================================================================
SUBROUTINE interpolate (self, time)
  class (rt_t) :: self
  real(8):: time
  !.....................................................................
  !---------------------------------------------------------------------
END SUBROUTINE interpolate

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

!> =============================================================================
!> Apply boundary conditions
!> =============================================================================
SUBROUTINE boundary_conditions (self)
  class(rt_t):: self
  !----------------------------------------------------------------------------
  call trace%begin('experiment_t%boundary_conditions')
  associate (my => self%mesh(2)%gn)
  self%mem(:,  01:06,:,:,:,2) = 0.0
  self%mem(:,my-4:my,:,:,:,2) = 0.0
  end associate
  call trace%end
END SUBROUTINE boundary_conditions


END MODULE rt_mod
