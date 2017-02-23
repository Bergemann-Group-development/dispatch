MODULE rt_mod
  USE bits_mod
  USE eos_mod
  USE io_mod
  USE link_mod
  USE ray_geometry_v3_mod
  USE rt_boundary_mod
  !USE rt_components_mod
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
  type, public :: I_data_t
    real, dimension(:), pointer :: in  => null()                        ! (wl_bin) incoming intensity.
    real, dimension(:), pointer :: out => null()                        ! leaving intensity
    integer, pointer :: dir_ID
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
    integer :: id                                                       ! location id in the ring array of a particular patch
    logical :: master = .false.                                         ! Am I the ring to rule them all?
    logical :: BC(3) = .false.                                          ! outer/inner/virtual or level
    integer :: l_chain = 0                                              ! length of this chain, points to one number, set on the head of the chain
    integer :: n_ready                                                  ! number of ready patches in this chain. When n_ready == l_chain, the chain can be executed.
    integer :: chain_nm                                                 ! ring number in the chain
  contains
    procedure :: check_masters
    procedure :: connect_to_upstream_diff_level
    procedure :: finish_chain
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
    type(I_data_t), dimension(:,:,:,:), pointer:: I_yz => null()        ! intensity shell data for zy plane (+/-x wall), it is primaty, that is, the largest and search starts with this one.
    type(I_data_t), dimension(:,:,:,:), pointer:: I_xz => null()        ! (x dim | z dim | n_out_dirs | upper(1)/lower(2)). This is the secondary one, meaning z has -2 length. Tertiary has both dims -2.
    type(I_data_t), dimension(:,:,:,:), pointer:: I_xy => null()        ! n_out_dirs is the number of directions leaving any given surface (3 dirs = 1 n_out_dir, 13 dirs = 9 n_out_dirs). This is tertiary surface.
    type(link_t), pointer             :: mhd_link                       ! link to my mhd_patch
    integer                           :: dim(3)                         ! patch dimensions
    integer                           :: n_lambda                       ! number of wavelengths
    integer                           :: n_bundle                       ! number of bundles
    integer                           :: i_q                            ! if q is stored for printing
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
    procedure :: get_boundary_test
    procedure :: execute_chain
    procedure :: finish
    procedure :: has_OB
    procedure :: has_IB
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
  class(I_data_t), dimension(:,:,:), pointer:: wall
  integer :: ray, cn, ip, ir
  integer, dimension(:), pointer :: nb_ray_pos
  !---------------------------------------------------------------------
  this => self%prev                                                     ! first item is not needed
  
  do while(associated(this))                                            ! loop over master rings
    prev => this%prev                                                   ! this temporary chain is organised in reverse
    nullify(this%prev)                                                  ! free the %prev pointer
    ring => this
    cn=1
    
    do while(associated(ring))
      ring%head => this
      this%l_chain = this%l_chain + 1                                   ! increment the chain length (number of members)
      ring%chain_nm = cn                                                ! ring's list number in the chain
      cn = cn+1
      if(.not.ring%BC(1)) then                                          ! avoid master rings, whose incoming intensities are from the boundaries
      
        select case(ring%direction)                                     ! select the right direction of this ring
        case(1)                                                         ! forward = upstream is behind me
          do ir = 1, ring%bundle%n_rays                                 ! loop over the rays of the ring's bundle
            ip=ring%bundle%rays_rev(ir)                                 ! get the ray position in the upstream bundle
            ring%packet_data(ir)%in => pring%packet_data(ip)%out        ! connect the ray to the appropriate patch intensity shell position
          end do                       
        case(2)                                                         ! reverse = upstream is in front of me
          do ir = 1, ring%bundle%n_rays                               
            ip=ring%bundle%rays_fw(ir)                                  ! get the ray position in the upstream bundle
            ring%packet_data(ir)%in => pring%packet_data(ip)%out        ! connect the ray to the appropriate patch intensity shell position
          end do
        end select
        
      end if
      pring => ring                                                     ! this will now be the previous ring
      ring => ring%next                                                 ! move to the next ring
    end do

    allocate(this%nudged(this%l_chain))                                 ! allocate the nudging logical
    this%nudged(:) = .false.                                            ! nothing is nudged yet.
    this => prev                                                        ! point to next master ring
  end do
  !if (io%verbose>1) call self%check_masters
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

END SUBROUTINE check_masters

SUBROUTINE connect_to_upstream_diff_level(self, ring)
  class(chain_ring_t) :: self
  class(chain_ring_t) :: ring
  
  error stop "rt_t:: connect_to_upstream_diff_level: not implemented"
END SUBROUTINE connect_to_upstream_diff_level

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
  !.....................................................................
  class (ray_bndl_t), pointer :: ray_bnd
  class (ray_dirs_t), pointer :: ray_dir
  class (packet_t), pointer      :: packet
  class (chain_ring_t), pointer :: ring
  logical, save :: first_time= .true.
  logical, save :: do_output  = .false.
  real(8), save :: dtime=0.01_8, grace=0.01_8, init_time=0d0
  real,    save :: cr_heating=0.0
  character(len=16) :: solver = 'integral'
  integer :: a_id, dd, n_out(3)
  integer:: ignore, ibnd
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
  
  self%nv = 1                                                           ! time slots for dedt
  self%nt = 1                                                           ! only one time slot FIXME 
  self%nw = 1                                                           ! no time deriv  
  call self%patch_t%init                                                ! initialize patch
  call self%set (bits%static+bits%no_download)                          ! special bits
  self%dtime = dtime                                                    ! 0.1 second default
  self%time = init_time                                                 ! 0.0 default
  self%dt_fixed = dtime                                                 ! timestep is fixed for radiative transfer
  self%grace = 0.0                                                    ! grace on by default, but very small, so nboring rt_patches would nod do anything stupid
  self%kind = 'rt_patch'                                                ! set the patch kind for output post-processing purposes
  !> -------------------------------------------------------------------
  !> modify the only memory slot to hold the dedt value, which is applied
  !> used directly by the rmhd patch.
  !> -------------------------------------------------------------------
  self%i_e = 1                                                          ! dedt slot
  self%dedt => self%mem(:,:,:,self%i_e,self%it,1)                       ! map dedt to mem
  
  
  
  self%n_lambda = eos%n_lambda                                          ! number of wavelengths
  self%n_bundle = ray_geometry%head%total_n_bundles * 2                 ! double n_bundles since there are 2 directions
  self%cr_heating = cr_heating                                          ! store cr heating, if it is present
  self%do_chains = ray_geometry%head%do_chains                          ! are there chains in the system, or patches are effectively isolated
  
  !> -------------------------------------------------------------------
  !> set up rt-specific parameters : weights for frequency bins and
  !> space-angles. Also, select the appropriate solver.
  !> -------------------------------------------------------------------
  allocate (self%w_lambda(self%n_lambda))                               ! frequency bin weights
  allocate (self%w_bundle(self%n_bundle))                               ! space angle weights
  self%w_lambda = eos%w_lambda
  self%w_bundle = pi4/2./self%n_bundle                                  ! space angle weights FIXME: could be different weights for different directions
  
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
  
  
  !> -------------------------------------------------------------------
  !> Allocate the outer intensity shell. In this shell all the leaving
  !> intensities are held, therefore it also acts as a hub for communications
  !> between nboring rt patches. 
  !> -------------------------------------------------------------------
  n_out = ray_geometry%head%wall_n_dirs
  ! --- Allocate the intensity outer shell ---
  allocate(self%I_yz(self%dim(2), self%dim(3), n_out(1), 2))
  allocate(self%I_xz(self%dim(3)-2, self%dim(1), n_out(2), 2))
  allocate(self%I_xy(self%dim(1)-2, self%dim(2)-2, n_out(3), 2))
  ! FIXME: options for collapsed dimensions.
  
  !> -------------------------------------------------------------------
  !> set up all the rings, belonging to this patch. Every ring represents
  !> a separate bundle. One rt patch manages all the directions.
  !> -------------------------------------------------------------------
  a_id = 1                                                              ! ring personal ID in the ring array
  ray_dir => ray_geometry%head%ray_list%head                            ! pick a ray direction
  do while(associated(ray_dir))                                         ! loop over all ray directions
    do ibnd =1, ray_dir%n_bundles                                       ! loop over ray bundles
      ray_bnd => ray_dir%bundle(ibnd)                        
      do dd=1,2                                                         ! separate rings for forward/reverse directions
        self%my_rings(a_id)%link => self%link                           ! rings should know, what rt patch hosts them
        self%my_rings(a_id)%direction = dd                              ! forward/reverse direction rings
        self%my_rings(a_id)%id = a_id                                   ! ring personal ID in the ring array
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
  
  !> -------------------------------------------------------------------  
  !> First, we have to determine which shell walls the rays cross when 
  !> they leave the patch. 
  !> select the right wall and coordinate axis. X axis is primary, Y axis
  !> is secondary, Z axis is tertiary.
  !> -------------------------------------------------------------------  
  if      ((ring%bundle%lwall(1) /= 0).and.(dir ==2)) then
    wall => self%I_yz(:,:,:, dir)                                       ! X wall is selected
    a = 2;  b = 3                                                       ! the remaining axes, which set the coordinate on the wall
  else if ((ring%bundle%uwall(1) /= 0).and.(dir ==1)) then  
    wall => self%I_yz(:,:,:, dir)                                       ! X wall is selected
    a = 2;  b = 3                                                       ! the remaining axes, which set the coordinate on the wall
  else if ((ring%bundle%lwall(2) /= 0).and.(dir ==2)) then
    wall => self%I_xz(:,:,:, dir)                                       ! Y wall is selected
    a = 3;  b = 1                                                       ! the remaining axes, which set the coordinate on the wall
  else if ((ring%bundle%uwall(2) /= 0).and.(dir ==1)) then
    wall => self%I_xz(:,:,:, dir)                                       ! Y wall is selected
    a = 3;  b = 1                                                       ! the remaining axes, which set the coordinate on the wall
  else if ((ring%bundle%lwall(3) /= 0).and.(dir ==2)) then
    wall => self%I_xy(:,:,:, dir)                                       ! Z wall is selected
    a = 1;  b = 2                                                       ! the remaining axes, which set the coordinate on the wall
  else if ((ring%bundle%uwall(3) /= 0).and.(dir ==1)) then
    wall => self%I_xy(:,:,:, dir)                                       ! Z wall is selected
    a = 1;  b = 2                                                       ! the remaining axes, which set the coordinate on the wall
  else 
    error stop "map_intensity_shell:: the right wall was somehow not found."        
  end if
  
  rid=1
  do ipck =1, ring%bundle%n_packets
    packet => ring%bundle%packet(ipck)                                  ! loop over ray packets in the bundle                                        
    if (dir ==1) then
       pt = packet%n_points                                             ! select the last point in the ray
    else
      pt = 1                                                            ! it might be the first point for reverse direction.
    end if
    do ray = 1, packet%n_rays                                           ! loop over rays in the packet
      selected = .false.
      c = packet%coords(:,pt,ray)                                       ! get the cell position of the last point
      i=1
      do while (.not.selected)                                          ! loop over the data in the selected cell on the shell
        this => wall(c(a),c(b),i)
        if (associated(this%out)) then                                  ! this entry point is already occupied by another ray (direction)
          i = i+1                                                       ! switch to the next entry point.
        else        
          allocate(wall(c(a),c(b),i)%out(self%n_lambda))                ! allocate the outgoing intensity storage
          allocate(wall(c(a),c(b),i)%dir_id)                            ! allocate a pointer to direction id, to know, who owns the entry point
          wall(c(a),c(b),i)%dir_ID = ring%bundle%dir_ID                 ! know, which direction this intensity is from
          ring%packet_data(rid)%out => wall(c(a),c(b),i)%out            ! point the fast lookup to the shell cell position
          ring%packet_data(rid)%dir_ID =>  wall(c(a),c(b),i)%dir_ID     ! point to the direction ID to have the reference.
          selected = .true.                                             ! this ray has chosen it's entry point on the rt outer shell.
          ring%packet_data(rid)%out = 0.
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
  
  call trace%begin ('rt_t%update', itimer=itimer)
  if (io%verbose > 1) print *, 'RADIATIVE TRANSFER TIME', self%time
  
  if (self%time < 0d0) then
    self%time  = 0d0 !+ self%dtime                                      ! set current time, if it was negative (to do the first RT run before MHD updates)
  end if
  
  self%dedt = 0.0                                                       ! reset net cooling
  
  task => self%mhd_link%task                                            ! pick the rmhd task
  select type (task)
  class is (patch_t)
    mhd => task
  end select
  
  !> -------------------------------------------------------------------
  !> allocate all the scratch arrays
  !> -------------------------------------------------------------------
  allocate (self%lnd (self%dim(1), self%dim(2), self%dim(3)))           ! log density scratch
  allocate (self%lne (self%dim(1), self%dim(2), self%dim(3)))           ! log internal energy scratch
  allocate (self%rk  (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  allocate (self%src (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  
  !> -------------------------------------------------------------------
  !> interpolate in time and space to get densoties and internal energies
  !> -------------------------------------------------------------------
  call mhd%log_interpolate (self%time, mhd%i_d, self%lnd)               
  call mhd%log_interpolate (self%time, mhd%i_e, self%lne)               
  
  if (io%verbose > 1) then
    print *, self%id, 'lnd', minval(self%lnd), maxval(self%lnd)
    print *, self%id, 'lne', minval(self%lne), maxval(self%lne)
  end if
  
  self%lne = self%lne - self%lnd                                        ! internal energy per unit mass

  call eos%lookup (shape(self%lne), lne=self%lne, lnd=self%lnd, &       ! call the table lookup, ...
                   rk=self%rk, src=self%src)  
  
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
  !.....................................................................
  class(chain_ring_t), pointer :: ring
  class(task_t), pointer :: tmp
  integer :: r_id
  !---------------------------------------------------------------------
  do r_id = 1, self%n_bundle
    ring => self%my_rings(r_id)
    !$omp critical (chain_cr)
    if (io%verbose > 4) then
      print *, 'ring no.',r_id
      print *, 'before nudge, n_ready',ring%head%n_ready, ring%head%l_chain
      print *,'nudged by master', ring%master,'rt time',self%time
      print *, 'other rings updated',self%updated
      print *, 'my data is ready',self%data_ready
      print *, 'my ID is:',self%id,'ring', ring%id
    end if
    
    tmp => ring%head%link%task
    select type(tmp)
    class is (rt_t)
    end select
    
    if (io%verbose > 4) print *, 'nudging task ID:', tmp%id,'ring id', ring%head%id, &
    'nudged allocated:',associated(ring%head%nudged),'size'
    if (io%verbose > 4) print *, size(ring%head%nudged),'its content:',ring%head%nudged, 'head assoc', associated(ring%head)
    if (io%verbose > 4) print *, 'l_chain',ring%head%l_chain
    
    ring%head%n_ready = ring%head%n_ready+1                             ! nudge the master ring
    ring%head%nudged(ring%chain_nm) = .true.                            ! mark myself as done in the ready array
        
    if (all(ring%head%nudged)) call self%execute_chain(ring%head)       ! execute the chain, when all rings are ready

    !$omp end critical (chain_cr) 
  end do
END SUBROUTINE nudge

!> =====================================================================
!> Loop through all rings in the chain and solve the actual rt
!> =====================================================================
SUBROUTINE execute_chain (self, head)
  class(rt_t) :: self
  class(chain_ring_t), pointer :: head
  !.....................................................................
  class(rt_t),          pointer :: rt_patch
  class(task_t),        pointer :: task
  class(chain_ring_t),  pointer :: ring, prev
  real, dimension(:,:), pointer :: i_in, i_out
  integer :: ray
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  !call trace%begin ('rmhd_t%execute_chain', itimer=itimer)
  ring => head                                                          ! take the master ring
  nullify(prev)                                                         ! just in case - there has to be no "prev" for the master ring
  
  allocate(i_in (size(ring%packet_data),self%n_lambda))                 ! allocate the intensity scratch
  allocate(i_out(size(ring%packet_data),self%n_lambda))                 ! allocate the intensity scratch
  i_in = 0
  i_out=0
  do while (associated(ring))
    task => ring%link%task                                              ! get the patch of this chain ring
    select type (task)
    class is (rt_t)
    rt_patch => task
    end select

    if(associated(prev)) then                                           ! if there was a ring before me, take the intensities from it
      do ray=1, ring%bundle%n_rays
        I_in(ray,:) = ring%packet_data(ray)%in
      end do
    else                                                                ! if I am the master ring,
      call rt_patch%get_boundary(ring,i_in)                             ! get intensities from the boundary
    end if
        
    call rt_patch%solve(ring, i_in, i_out)                              ! do the radiative transfer for this particular ring
    do ray = 1, ring%bundle%n_rays                                      ! store the output intensities
      ring%packet_data(ray)%out = I_out(ray,:)
    end do                                                              
    rt_patch%updated(ring%id) = .true.                                  ! mark, that this ring has been updated
    if(all(rt_patch%updated)) call rt_patch%finish                      ! if all rings are done, this rt_patch is complete
    prev => ring
    ring => ring%next
  end do
  deallocate(i_in, i_out)                                               ! deallocate the scratch
  head%nudged = .false.                                                 ! reset the nudge logicals
  head%n_ready = 0                                                      ! no rings are ready in this chain
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
  !.....................................................................
  class(packet_t), pointer :: packet
  integer :: pck, ray
  real(kind=8) :: ppos(3)
  !---------------------------------------------------------------------
  !> -------------------------------------------------------------------
  !> If the ring is adjacent to the outer boundary
  !> -------------------------------------------------------------------
  if (ring%BC(1)) then
    do pck =1, ring%bundle%n_packets                                    ! loop over all packets in this bundle
      packet => ring%bundle%packet(pck)                                 ! pick a packet
      do ray = 1, packet%n_rays                                         ! loop over all rays in this packet
        ppos = self%absolute_point_position &                           ! for each ray index in the packet, find a point position
        (ring, packet, ray, 1, ring%direction)
        call rt_boundary%set_obc (ppos, I_in(ray,:))                    ! get the boundary value
      end do
    end do
  end if
 
  !> -------------------------------------------------------------------
  !> If the ring is adjacent to the inner boundary
  !> -------------------------------------------------------------------
  if (ring%BC(2)) then                                                   ! inner boundary should modify the ray lengths!!
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
  
  !> -------------------------------------------------------------------
  !> If the ring is adjacent to a virtual or different level patch
  !> -------------------------------------------------------------------
  if (ring%BC(3)) then
      error stop "getting boundary from virtual patches or different level patches is not implemented yet"
  end if
  
END SUBROUTINE get_boundary

SUBROUTINE get_boundary_test (self, ring, I_in)
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
        if (ring%direction==1) then
        I_In(ray,:) = ray*ring%bundle%axis/10.
        else
        I_In(ray,:) = -ray*ring%bundle%axis/10.
        end if
      end do
    end do
  end if
  
  ! inner boundary should modify the ray lengths!!

END SUBROUTINE get_boundary_test
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
    if (mine%bundle%b_id_rev  /= other%bundle%ID) return                 ! wrong bundle to make a chain
  case(2)
    if (mine%bundle%b_id_fw /= other%bundle%ID) return                 ! wrong bundle to make a chain
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
  if (self%do_output) call self%output

  deallocate (self%lnd, self%lne, self%rk, self%src)                                       ! deallocate scratch
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
SUBROUTINE solve_integral (self, ring, i_in, I_out, i_in2, i_out2)
  class (rt_t) :: self
  class(chain_ring_t), pointer  :: ring
  real, dimension(:,:), pointer :: i_in, i_out
  real, dimension(:,:), pointer, optional :: i_in2, i_out2
  !.....................................................................
  class(packet_t),        pointer :: packet
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  real, dimension(:,:),   pointer :: dwf, I_in0, I_out0
  integer                         :: nrays, npt, nwb, ipck, ir
  real                            :: ds
  logical :: one_p
  integer, save                   :: itimer=0
  class(task_t), pointer :: tmp
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve_integral', itimer=itimer)

  nwb   =  self%n_lambda                                                ! number of wl bins
  do ipck = 1, ring%bundle%n_packets                                    ! loop over packets
    packet => ring%bundle%packet(ipck)                                  ! pick a packet
    nrays =  packet%n_rays                                              ! number of rays in this packet
    
    !> -----------------------------------------------------------------
    !> RT solvers need at least 3 points to work with, as they use
    !> derivatives. Rays, shorter, than 3 points, have to be extended
    !> to patch boundaries to have accurate solutions.
    !> TODO: inner boundaries modify npt
    !> -----------------------------------------------------------------
    if (packet%n_points < 3) then
      one_p = .true.
      npt = packet%n_points + 2
    else
      npt = packet%n_points
      one_p = .false.
    end if
    
    ds = packet%ds                                                      ! ray patch length between the points
    
    !> -----------------------------------------------------------------
    !> allocate all the scratches
    !> -----------------------------------------------------------------
    allocate(s     (nrays, nwb, npt))                                   ! source function
    allocate(rk    (nrays, nwb, npt))                                   ! opacity
    allocate(dtau  (nrays, nwb, npt))                                   ! optical depth between the points
    allocate(exp0  (nrays, nwb, npt))                                   ! first exponential factor
    allocate(q     (nrays, nwb, npt))                                   ! heating rates
    allocate(dwf   (nrays,      npt))                                   ! density weight factor
    allocate(I_in0 (nrays, nwb))                                        ! incoming rays
    allocate(I_out0(nrays, nwb))                                        ! outgoing rays
    
    do ir=1, packet%n_rays
      I_in0(ir,:) = I_in(packet%b_id(ir),:)                             ! take the appropriate intensities
    end do
   
    !> -----------------------------------------------------------------
    !> map the quantities onto rays
    !> -----------------------------------------------------------------
    call self%ray_vol_map4(packet, self%src, s, one_p)                  ! map source function to rays
    call self%ray_vol_map4(packet, self%rk, rk, one_p)                  ! map opacities to rays
    call self%optical_depth_centered(ds, rk, dtau)                      ! optical depth is centered for integral formulation
    call self%ray_vol_map3(packet, self%lnd, dwf, one_p)                ! density weight factor
    dwf = dwf * ds  
    if (io%verbose>2) &
        print *,'dtau',self%id, minval(dtau), maxval(dtau)

#if defined __INTEL_COMPILER
      exp0 = 2**(-dtau*1.44269504) 
#else
      exp0 = exp(-dtau)                                                 ! expopential factors
#endif        
    
    if (ring%direction ==1) then
      call rt_solver%forward (dtau, exp0, s, I_in0, q, I_out0, one_p)   ! solve forward
    else
      call rt_solver%reverse (dtau, exp0, s, I_in0, q, I_out0, one_p)   ! solve reverse
    end if
    
    do ir=1, packet%n_rays
      I_out(packet%b_id(ir),:) = I_out0(ir,:)                           ! store the leaving intensities to their right spots
    end do
    
    call self%update_q(packet, q, dwf, ring%id, one_p)                  ! update the heating rates
       
    deallocate (s, rk, dtau, exp0, q, dwf, I_out0, I_in0)   
  end do
  call trace_end (itimer)
END SUBROUTINE solve_integral

!> =====================================================================
!> Same as solve_integral, but for the Feautrier formulation
!> =====================================================================
SUBROUTINE solve_feautrier (self, ring, i_in, i_out1, i_in2, i_out2)
  class (rt_t) :: self
  class (chain_ring_t), pointer :: ring
  real, dimension(:,:), pointer :: i_in, i_out1
  real, dimension(:,:), pointer, optional :: i_in2, i_out2
  !.....................................................................
  class(packet_t),        pointer :: packet
  integer                         :: nrays, npt, nwb, ipck
  real :: ds
  real, dimension(:,:), pointer :: dwf
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
    call self%ray_vol_map4(packet, self%src, s, one_p)                  ! map source function to rays
    call self%ray_vol_map4(packet, self%rk, rk, one_p)                  ! map opacities to rays
    call self%optical_depth_staggered(ds, rk, dtau)                     ! optical depth is staggered for Feautrier formulation
    call self%ray_vol_map3(packet, self%lnd, dwf, one_p)                ! density weight factor
    dwf = dwf * ds  
    if (io%verbose>1) &
        print *,'dtau',self%id, minval(dtau), maxval(dtau)
    
    call rt_solver%transfer (dtau, s, I_in, I_in2, q, I_out1, I_out2)        ! solve bi-directional
    call self%update_q(packet, q, dwf, ring%id, one_p)
    
    deallocate (s, rk, dtau, exp0, q, dwf)
  end do

  call trace_end (itimer)
END SUBROUTINE solve_feautrier

!> =====================================================================
!> Map cell data onto rays
!> =====================================================================
SUBROUTINE ray_vol_map4 (self, packet, in, out, one_p)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:,:) , pointer, intent(in)  :: in
  logical :: one_p
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), p2(3), ray, pt, nwb, wl, lp
  !---------------------------------------------------------------------
  cell => packet%coords
  nwb = size(out,2)
  lp  = size(out,3)
  !> -------------------------------------------------------------------
  !> for short rays, we need to consider the first and last points separately,
  !> as they do not belong to these rays implicitly.
  !> -------------------------------------------------------------------
  if (one_p) then
    do ray = 1, packet%n_rays                                           ! loop over all rays
      do pt = 1, packet%n_points                                        ! loop over all points
        p = cell(:,pt,ray)                                              ! get the right cell coordinate
        !DIR$ SIMD
        do wl =1, nwb                                                   ! loop over wl bins
          out(ray,wl,pt+1) = in(p(1),p(2),p(3),wl)                      ! first point is missing
        end do
        !DIR$ SIMD
      end do
      p  = packet%rev_coords(:,ray)
      p2 = packet%fwd_coords(:,ray)
      do wl=1, nwb                                                      ! fill in the first and last points
        out(ray,wl,1)  = in(p(1),p(2),p(3), wl)                         ! Get the first point
        out(ray,wl,lp) = in(p2(1),p2(2),p2(3), wl)                      ! Get the last point
      end do
    end do
  else
  !> -------------------------------------------------------------------
  !> for all other rays, no special care is needed.
  !> -------------------------------------------------------------------
    do ray = 1, packet%n_rays                                           ! loop over all rays
      do pt = 1, packet%n_points                                        ! loop over all points
        p = cell(:,pt,ray)                                              ! get the right cell coordinate
        !DIR$ SIMD
        do wl =1, size(in,4)                                            ! loop over wl bins
          out(ray,wl,pt) = in(p(1),p(2),p(3),wl)                        ! map the quantity
        end do
        !DIR$ SIMD
      end do
    end do
  end if
    
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
SUBROUTINE ray_vol_map3 (self, packet, in, out, one_p)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:) , pointer, intent(in)  :: in
  logical :: one_p
  !.....................................................................
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), ray, pt, lp
  !---------------------------------------------------------------------
  cell => packet%coords
  lp  = size(out,2)  
  if (one_p) then
    do ray = 1, packet%n_rays                                           ! loop over rays
      do pt = 1, packet%n_points                                        ! loop over points
        p = cell(:,pt,ray)                                              ! get the right cell
          out(ray,pt+1) = in(p(1),p(2),p(3))
      end do
      p  = packet%rev_coords(:,ray)
      out(ray,1) = in(p(1),p(2),p(3))                                   ! Get the first point
      p = packet%fwd_coords(:,ray)
      out(ray,lp) = in(p(1),p(2),p(3))                                  ! Get the last point
    end do
  else
    do ray = 1, packet%n_rays
      do pt = 1, packet%n_points
        p = cell(:,pt,ray)
        out(ray,pt) = in(p(1),p(2),p(3))
      end do
    end do
  end if
  
END SUBROUTINE ray_vol_map3

SUBROUTINE update_q (self, packet, q, dwf, ib, one_p)
  class(rt_t) :: self
  class(packet_t), pointer :: packet
  real, dimension(:,:,:), pointer :: q
  real, dimension(:,:), pointer :: dwf
  real, dimension(:,:,:), pointer :: qq
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), wl, ray, pt, ib
  integer, save:: itimer=0
  logical :: one_p
  !---------------------------------------------------------------------
  call trace_begin ('rt_t%update_q', 4, itimer=itimer)
  qq => self%dedt(:,:,:)                                                ! point to the memory slot
  cell => packet%coords                                                 ! point to coordinates
  
  if (one_p) then
    do ray = 1, packet%n_rays                                           ! loop over rays
      do pt = 1, packet%n_points                                        ! loop over points
        p = cell(:,pt,ray)                                              ! get the cell position
        do wl =1, size(q,2)
          qq(p(1),p(2),p(3)) = &
          qq(p(1),p(2),p(3)) + q(ray,wl,pt+1)* self%w_lambda(wl) * &
           self%w_bundle(ib) !* exp(dwf(ray,pt+1))
        end do
      end do
    end do
  else
    do ray = 1, packet%n_rays
      do pt = 1, packet%n_points
        p = cell(:,pt,ray)
        do wl =1, size(q,2)
          qq(p(1),p(2),p(3)) = &
          qq(p(1),p(2),p(3)) + q(ray,wl,pt)* self%w_lambda(wl) * &
           self%w_bundle(ib) !* exp(dwf(ray,pt))
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
