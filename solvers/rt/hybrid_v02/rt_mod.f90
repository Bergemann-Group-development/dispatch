MODULE rt_mod
  USE bits_mod
  USE eos_mod
  USE io_mod
  USE link_mod
  USE ray_geometry_v3_mod
  USE rt_boundary_mod
  USE rt_bundle_mod
  USE patch_mod
  USE task_mod
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
  type :: chain_ring_t
    class(link_t),                 pointer :: link => null()            ! link to rt_t task
    class(chain_ring_t),           pointer :: head => null()            ! pointer to the master chain ring
    class(chain_ring_t),           pointer :: tail => null()            ! pointer to the last ring
    class(chain_ring_t),           pointer :: next => null()
    class(chain_ring_t),           pointer :: prev => null()
    class(i_data_t), dimension(:), pointer :: packet_data => null()     ! quick pointing to the right locations in the data sets
    class(ray_bndl_t),             pointer :: bundle => null()          ! ray bundle of this chain ring
    integer :: direction                                                ! 1 = forward / 2 = reverse
    integer :: a_id                                                     ! location id in the ring array of a particular patch
    logical :: master = .false.
    logical :: BC(3) = .false.                                          ! outer/inner/virtual or level
    integer :: l_chain                                                  ! length of this chain, points to one number, set on the head of the chain
    integer :: n_ready                                                  ! number of ready patches in this chain. When n_ready == l_chain, the chain can be executed.
  contains
    procedure :: ring_connects_to
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
    integer                           :: dim(3)                         ! patch dimensions
    integer                           :: n_lambda                       ! number of wavelengths
    integer                           :: n_bundle                       ! number of bundles
    real, dimension(:), pointer       :: w_lambda                       ! wavelength weights
    real, dimension(:), pointer       :: w_bundle                       ! space angle weights
    real, dimension(:,:,:),   pointer :: lnd  => null()
    real, dimension(:,:,:),   pointer :: lne  => null()
    real, dimension(:,:,:,:), pointer :: rk   => null()
    real, dimension(:,:,:,:), pointer :: src  => null()
    real, dimension(:,:,:),   pointer :: dedt => null()
    character(len=16)                 :: solver ='integral'             ! choice of solver
    procedure(solve_integral), pointer:: solve => null()
  contains
    procedure :: absolute_point_position
    procedure :: boundary_conditions
    procedure :: cell2pos
    procedure :: density_weight_factor
    procedure :: dnload
    procedure :: get_boundary
    procedure :: execute_chain
    procedure :: finish
    procedure :: get_upstream
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
    procedure :: solve_integral
    procedure :: solve_feautrier
    procedure :: store_I_out
    procedure :: update => update_rt
    procedure :: update_q  
  end type
  real, parameter:: pi4 = 16.*atan(1.)
  
CONTAINS

!> =====================================================================
!> We may or may not need locally allocated arrays.  If the calling
!> procedure holds the necessary data (which it does, since it computes
!> them and is picking up the output), it just needs to set the pointers
!> in this module to point to the existing arrays.  However, the optical
!> depth increments and the corresponding exp(-dtau) factors are
!> computed here, so they need array allocation.
!> =====================================================================
SUBROUTINE init_rt (self)
  class (rt_t) :: self
  class (ray_bndl_t), pointer :: ray_bnd
  class (ray_dirs_t), pointer :: ray_dir
  class (ray_t), pointer      :: packet
  class (chain_ring_t), pointer :: ring
  logical, save:: first_time=.true.
  real(8), save:: dtime=0.01_8, grace=0.0_8
  character(len=16) :: solver = 'integral'
  integer :: a_id, ld, ud, dd, n_out
  integer:: ignore
  namelist /rt_params/ dtime, grace, solver
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%init')
  !$omp critical (init_cr)
  if (first_time) then
    rewind (io%input); read (io%input, rt_params, iostat=ignore)
    if (io%master) write (*, rt_params)
    call rt_bundle%init
  end if
  !$omp end critical (init_cr)
  self%nv = 2                                                           ! time slots for dedt and one for Q or J - TO DO
  self%nw = 1                                                           ! no time deriv
  call self%patch_t%init                                                ! initialize patch
  call self%set (bits%static+bits%no_download)                          ! special bits
  self%i_e = 1                                                          ! dedt slot
  self%n_lambda = eos%n_lambda                                          ! number of wavelengths
  self%n_bundle = ray_geometry%head%ray_list%head%total_n_bundles
  self%dedt => self%mem(:,:,:,self%i_e,self%it,1)                       ! map dedt to mem
  allocate (self%w_lambda(self%n_lambda))
  allocate (self%w_bundle(self%n_bundle))
  allocate (self%rk  (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  allocate (self%src (self%dim(1), self%dim(2), self%dim(3), self%n_lambda))
  self%w_lambda = eos%w_lambda
  self%w_bundle = pi4/2./self%n_bundle                                  ! space angle weights FIXME: could be different weights for different directions
  self%dtime = dtime                                                    ! 0.1 second default
  self%grace = grace                                                    ! grace off by default
  self%solver = solver
  select case(self%solver)
  case('integral')
    self%solve => solve_integral
  case('feautrier')
    self%solve => solve_feautrier
  case default
    print *, 'WARNING: rt solver not recognised. Using integral solver.'
    self%solve => solve_integral
  end select
  
  allocate(self%updated(self%n_bundle))
  allocate(self%my_rings(self%n_bundle))
  
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
  
  a_id = 1
  ray_dir => ray_geometry%ray_list%head
  do while(associated(ray_dir))
    ray_bnd => ray_dir%bundle%head
    do while(associated(ray_bnd))
      if (ray_geometry%head%ray_list%head%directionality(ray_bnd%axis) == 1) then
        ld = 1; ud = 2
      else if (ray_geometry%head%ray_list%head%directionality(ray_bnd%axis) == 2) then
        ld = 1; ud = 1
      else
        ld = 2; ud = 2
      end if
      do dd=ld,ud
        allocate(self%my_rings(a_id))
        ring => self%my_rings(a_id)
        ring%link => self%link
        ring%direction = dd
        ring%a_id = a_id
        ring%bundle => ray_bnd
        allocate(ring%packet_data(ray_bnd%n_rays)) 
        call self%map_intensity_shell(ring)
        a_id = a_id + 1
      end do
      ray_bnd => ray_bnd%next
    end do
    ray_dir => ray_dir%next
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
!> =====================================================================
SUBROUTINE map_intensity_shell(self, ring)
  class(rt_t) :: self
  class(chain_ring_t), pointer :: ring
  !.....................................................................
  class(ray_t), pointer :: packet
  class(I_data_t), dimension(:,:,:), pointer:: wall
  class(I_data_t), pointer :: this
  integer :: ray, rid, c(3), pt, dir, a,b, x,y, i
  logical :: selected
  !---------------------------------------------------------------------
  dir = ring%direction
  packet => ring%bundle%rays
  
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
  x=a;  y=b
  rid=1
  do while (associated(packet))                                         ! loop over ray packets in the bundle
    if (dir ==1) then
       pt = packet%n_points
    else
      pt = 1
    end if
    do ray = 1, packet%n_rays                                           ! loop over rays in the packet
      selected = .false.
      allocate(ring%packet_data(rid))                                   ! allocate fast lookup position
      c = packet%coords(:,ray,pt)
      i=1
      do while (.not.selected)
        this => wall(c(a),c(b),i)
        if (associated(this)) then                                       ! allocate wall data.
          i = i+1
        else
          allocate(wall(c(a),c(b),i))                                   ! allocate shell cell
          wall(c(a),c(b),i)%dir_ID = ring%bundle%dir_ID                 ! know, which direction this intensity is from
          allocate(wall(c(a),c(b),i)%out(self%n_lambda))                ! allocate the outgoing intensity storage
          ring%packet_data(rid)%out => wall(c(a),c(b),i)%out            ! point the fast lookup to the shell cell position
          selected = .true.
        end if
      end do
      rid=rid+1
    end do
    packet => packet%next
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
  call trace%begin ('rmhd_t%update', itimer=itimer)
  if (io%verbose > 1) print *, 'RADIATIVE TRANSFER TIME', self%time
  if (self%time < 0d0) then
    self%time  = 0.0                                                    ! set current time
  end if
  self%dedt = 0.0                                                       ! reset net cooling
  task => self%link%nbor%task                                           ! assumes the only nbor is RMHD
  select type (task)
  class is (patch_t)
    mhd => task
  end select
  allocate (self%lnd (mhd%gn(1), mhd%gn(2), mhd%gn(3)))                 ! log density scratch
  allocate (self%lne (mhd%gn(1), mhd%gn(2), mhd%gn(3)))                 ! log internal energy scratch
  call mhd%log_interpolate (self%time, mhd%i_d, self%lnd)               ! interpolate in time and space
  call mhd%log_interpolate (self%time, mhd%i_e, self%lne)               ! interpolate in time and space
  if (io%verbose > 1) then
    print *, self%id, 'lnd', minval(self%lnd), maxval(self%lnd)
    print *, self%id, 'lne', minval(self%lne), maxval(self%lne)
  end if
  self%lne = self%lne - self%lnd                                        ! internal energy per unit mass
  
  call eos%lookup (shape(self%lne), lne=self%lne, lnd=self%lnd, &       ! call the table lookup, ...
                   rk=self%rk, src=self%src)                            ! ... results go into rt_t    
  call trace%end (itimer)
  call self%nudge                                                       ! nudge the chain heads
END SUBROUTINE update_rt 

!> =====================================================================
!> Nudge the chain heads. If all patches in the chain become ready, then
!> solve for the net heating rate, storing it in a slot that rotates
!> with time for each patch in the chain.
!> =====================================================================
SUBROUTINE nudge (self)
  class(rt_t) :: self
  class(chain_ring_t), pointer :: ring
  integer :: r_id
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%nudge', itimer=itimer)
  do r_id = 1, self%n_bundle
    ring => self%my_rings(r_id)
    !$omp critical (chain_cr)
    ring%head%n_ready = ring%head%n_ready+1
    if (ring%head%n_ready >= ring%head%l_chain) call self%execute_chain(ring%head)
    !$omp end critical (chain_cr) 
  end do
  call trace%end(itimer)
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
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace%begin ('rmhd_t%execute_chain', itimer=itimer)
  ring => head
  nullify(prev)
  do while (associated(ring))
    task => ring%link%task                                              ! get the patch of this chain ring
    select type (task)
    class is (rt_t)
    rt_patch => task
    end select
    
    allocate(i_in(size(ring%packet_data),self%n_lambda))
    if(associated(prev)) then
      call rt_patch%get_upstream(ring,i_in)
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
  nullify(rt_patch, task, prev)  
  call trace%end(itimer)
END SUBROUTINE execute_chain

!> =====================================================================
!> Upstream data is on same refinement level. Direct lookup
!> =====================================================================
SUBROUTINE get_upstream (self, ring, I_in)
  class(rt_t) :: self
  class(chain_ring_t), pointer :: ring
  real, dimension(:,:), pointer :: I_in
  !---------------------------------------------------------------------
  class(chain_ring_t), pointer :: upstream
  integer :: ray
  integer, dimension(:), pointer :: nb_ray_pos
  !---------------------------------------------------------------------
  call trace%begin('rt_mod::getupstream')

  select case(ring%direction)
  case(1)
    nb_ray_pos => ring%bundle%ray_chain_rev
  case(2)
    nb_ray_pos => ring%bundle%ray_chain_fw
  end select

  upstream => ring%prev
  print *, "todo : update upstream"
  select case(ring%bundle%l_chain)
  case(:2)
 !   I_in = upstream%I_out
  case default
    do ray = 1, ring%bundle%n_rays
 !     I_in(ray,:) = upstream%I_out(nb_ray_pos(ray),:)
    end do
  end select
  
  nullify(upstream)
  call trace%end
END SUBROUTINE get_upstream
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
  class(ray_t), pointer :: packet
  integer :: pck, ray
  real(kind=8) :: ppos(3)
  !---------------------------------------------------------------------
  
  if (ring%BC(1)) then
    packet => ring%bundle%rays
    do pck =1, ring%bundle%n_packets
      do ray = 1, packet%n_rays
        ! for each ray index in the packet, find a point position
        ppos = self%absolute_point_position (ring, packet, ray, 1, ring%direction)
        call rt_boundary%set_obc (ppos, I_in(ray,:))
      end do
      packet => packet%next
    end do
  end if

  ! inner boundary should modify the ray lengths!!
  if (ring%BC(2)) then
    error stop "inner BC for RT is not fully functional yet"
    packet => ring%bundle%rays
    do pck =1, ring%bundle%n_packets
      do ray = 1, packet%n_rays
        ! for each ray index in the packet, find a point position
        ppos = self%absolute_point_position (ring, packet, ray, 1, ring%direction)
        call rt_boundary%set_ibc (ppos, I_in(ray,:))
      end do
      packet => packet%next
    end do
  end if
  
  if (ring%BC(3)) then
      error stop "getting boundary from virtual patches or different level patches is not implemented yet"
  end if
  
END SUBROUTINE get_boundary

!> =============================================================================
!> Check if the rays in self match up in the downstream with the rays in other
!> =============================================================================
FUNCTION ring_connects_to (self, other) RESULT (out)
  class (chain_ring_t)          :: self
  class (chain_ring_t), pointer :: other
  class(task_t), pointer :: task
  class (rt_t), pointer :: mine, nbor
  integer :: pp, ray, last_pt(3), rev(3)
  logical :: out, tout(3)
  !.............................................................................
  out = .false.
  tout=.false.
  if (self%bundle%dir_ID /= other%bundle%dir_ID)     return             ! different bundle
  if (self%direction /= other%direction)             return             ! different direction 
  
  ray = 1                                                               ! testing only one ray in the bundle for now... For AMR: needs to run
  select case(self%direction)                                           ! through all the rays and if ANY has a match - add to the list...
  case(1)
    last_pt =  self%bundle%rays%coords(:,self%bundle%rays%n_points,ray)! taking the coordinates of the last point in the ray 
    if (self%bundle%l_chain <2) then                                    
      pp = 1
    else
      pp = self%bundle%ray_chain_fw(1)                                  ! the ray location in the next bundle in the chain
    end if      
    rev = other%bundle%rays%fwd_coords(:,pp)                            ! zeroth point in the ray (location of incoming intensity)                             
  case(2)
    last_pt =  self%bundle%rays%coords(:,1,ray)                         ! first point is the last point in reverse   
    if (self%bundle%l_chain <2) then
      pp = 1
    else
      pp = self%bundle%ray_chain_rev(1)      
    end if 
    rev = other%bundle%rays%rev_coords(:,pp)                             
  end select
  task => self%link%task
  select type(task)
    class is (rt_t)
    mine => task
  end select
  
  task => other%link%task
  select type(task)
    class is(rt_t)
    nbor => task
  end select
  
  last_pt = nbor%pos2cell(mine%cell2pos(last_pt))
  do pp = 1,3
    if (last_pt(pp) == rev(pp)) tout(pp) = .true.
  end do
  if (all(tout)) out = .true.
  
  nullify(task)   
END function ring_connects_to

!> =====================================================================
!> Find the absolute position of index 'ray' in ray packet 'rays'.
!> If the optional bn parameter is present, the answer depends on
!> the direction, otherwise 'pt' chooses the point along the ray.
!> =====================================================================
FUNCTION absolute_point_position (self, ring, rays, ray, pt, bn) result (cell_pos)
  class (rt_t) :: self
  class (chain_ring_t), pointer :: ring
  class (ray_t), pointer :: rays
  integer, optional :: bn
  integer :: ray, pt, cp(3)
  real(kind=8) :: cell_pos(3)
!-----------------------------------------------------------------------
  call trace_begin ('rt_t%absolute_point_position', 4)
  if (present(bn)) then
    select case(ring%direction)
    case(1)
      ! self%my_rays
      cp = rays%fwd_coords(:,ray)
    case(2)
      cp = rays%rev_coords(:,ray)
    end select
  else
    cp = rays%coords(:,pt,ray)
  end if
  cell_pos = self%position-0.5*self%size+(cp-self%li)*self%ds
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
    i = ring%bundle%rays%fwd_coords(:,1)
    p = self%cell2pos(i)
  else 
    i = ring%bundle%rays%rev_coords(:,1)
    p = self%cell2pos(i)
  end if
  
  do t = 1,3
    if (d(t)>0) then
      if (p(t)<rt_boundary%lbox(t)) has_OB = .true.
    else
      if (p(t)>rt_boundary%ubox(t)) has_OB = .true.
    end if
  end do
END FUNCTION has_OB
  
SUBROUTINE finish(self)
  class(rt_t) :: self  
error stop "free mem not done"
  call self%boundary_conditions
  call self%output
  deallocate (self%lnd, self%lne)                                       ! deallocate scratch
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
  class(chain_ring_t), pointer :: ring
  real, dimension(:,:), pointer :: i_in
  real, dimension(:,:), pointer, optional :: i_in2
  !.....................................................................
  class(ray_t),           pointer :: packet
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  real, dimension(:,:),   pointer :: I_out, dwf
  integer                         :: nrays, npt, nwb
  real                            :: ds
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve_integral', itimer=itimer)

  nwb   =  self%n_lambda
  packet => ring%bundle%rays
  do while (associated(packet))
    nrays =  packet%n_rays
    npt   =  packet%n_points
    ds    =  packet%ds
    allocate(s   (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(rk  (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(dtau(nrays, nwb, npt+2))
    allocate(exp0(nrays, nwb, npt+2))
    allocate(q   (nrays, nwb, npt+2))
    allocate(dwf (nrays,      npt+2))
    allocate(I_out(nrays, nwb))
    call self%ray_vol_map4(packet, self%src, s)                         ! map source function to rays
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
      call rt_bundle%forward (dtau, exp0, s, I_in, q, I_out)            ! solve forward
    else
      call rt_bundle%reverse (dtau, exp0, s, I_in, q, I_out)            ! solve reverse
    end if
    
    call self%update_q(packet, q, dwf, ring%a_id)
    call self%store_I_out(ring,packet, I_out)
    
    deallocate (s, rk, dtau, exp0, q, dwf, I_out)   
    packet => packet%next
  end do
  call trace_end (itimer)
END SUBROUTINE solve_integral

SUBROUTINE store_I_out(self, ring, packet, I_out)
  class(rt_t) :: self
  class(ray_t), pointer :: packet
  class(chain_ring_t), pointer :: ring
  real, dimension(:,:), pointer :: I_out

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
  class(ray_t),           pointer :: packet
  integer                         :: nrays, npt, nwb
  real :: ds
  real, dimension(:,:), pointer :: i_out, i_out2, dwf
  real, dimension(:,:,:), pointer :: s, rk, dtau, exp0, q
  integer, save                   :: itimer=0
  ! --------------------------------------------------------------------
  call trace_begin ('rt_t%solve_feautrier', itimer=itimer)
  
  if (.not.present(i_in2)) error stop "solve Feautrier: needs 2 directions to act"
  
  packet => ring%bundle%rays
  do while (associated(packet))
    
    nrays =  packet%n_rays
    npt   =  packet%n_points
    ds    =  packet%ds
    allocate(s   (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(rk  (nrays, nwb, npt+2))            ! TODO: boundaries modify npt
    allocate(dtau(nrays, nwb, npt+2))
    allocate(exp0(nrays, nwb, npt+2))
    allocate(q   (nrays, nwb, npt+2))
    allocate(dwf (nrays,      npt+2))
    allocate(I_out(nrays, nwb))
    allocate(I_out2(nrays, nwb))
    call self%ray_vol_map4(packet, self%src, s)                         ! map source function to rays
    call self%ray_vol_map4(packet, self%rk, rk)                         ! map opacities to rays
    call self%optical_depth_staggered(ds, rk, dtau)                     ! optical depth is staggered for Feautrier formulation
    call self%ray_vol_map3(packet, self%lnd, dwf)                       ! density weight factor
    dwf = dwf * ds  
    if (io%verbose>1) &
        print *,'dtau',self%id, minval(dtau), maxval(dtau)
    
    call rt_bundle%transfer (dtau, s, I_in, I_out, q, I_in2, I_out2)        ! solve bi-directional
    call self%update_q(packet, q, dwf, ring%a_id)
    call self%store_I_out(ring,packet, I_out)
    call self%store_I_out(ring,packet, I_out2)
    
    deallocate (s, rk, dtau, exp0, q, dwf, I_out, I_out2)
    packet => packet%next
  end do

  call trace_end (itimer)
END SUBROUTINE solve_feautrier

!> =====================================================================
!> Map cell data onto rays
!> =====================================================================
SUBROUTINE ray_vol_map4 (self, rays, in, out)
  class(rt_t) :: self
  class(ray_t), pointer :: rays
  integer :: p(3), ray, pt, nwb, wl
  real, dimension(:,:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:,:) , pointer, intent(in)  :: in
  integer, dimension(:,:,:), pointer :: cell
  !---------------------------------------------------------------------
  cell => rays%coords
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
        out(ray,wl,1) = in(rays%fwd_coords(1,ray), &
        rays%fwd_coords(2,ray), rays%fwd_coords(3,ray), wl)             ! Get the first point
        out(ray,wl,size(out,3)) = in(rays%rev_coords(1,ray), &
        rays%rev_coords(2,ray), rays%rev_coords(3,ray), wl)             ! Get the last point
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
SUBROUTINE ray_vol_map3 (self, rays, in, out)
  class(rt_t) :: self
  class(ray_t), pointer :: rays
  integer :: p(3), ray, pt
  real, dimension(:,:)   , pointer, intent(out) :: out
  real, dimension(:,:,:) , pointer, intent(in)  :: in
  integer, dimension(:,:,:), pointer :: cell
  !---------------------------------------------------------------------
  cell => rays%coords
    do ray = 1, size(out,1)
      do pt = 2, size(out,2)-1                                          ! First and last points ommitted
        p = cell(:,pt-1,ray)
        !DIR$ SIMD
          out(ray,pt) = in(p(1),p(2),p(3))
      end do
      !DIR$ SIMD
        out(ray,1) = in(rays%fwd_coords(1,ray), &
        rays%fwd_coords(2,ray), rays%fwd_coords(3,ray))                 ! Get the first point
        out(ray,size(out,2)) = in(rays%rev_coords(1,ray), &
        rays%rev_coords(2,ray), rays%rev_coords(3,ray))                 ! Get the last point
    end do
    
  if (io%verbose > 1) then
    if (any(in==0)) then
      print *, in
      error stop "ray vol map3 in=0"
    end if
    if (any(out==0)) then
      do pt = 1, size(out,2)
      print *,"pt",pt, out(:,pt)
      end do
      error stop "ray vol map3 out=0"
    end if
  end if
  
END SUBROUTINE ray_vol_map3

SUBROUTINE update_q (self, rays, q, dwf, ib)
  class(rt_t) :: self
  class(ray_t), pointer :: rays
  real, dimension(:,:,:), pointer :: q
  real, dimension(:,:), pointer :: dwf
  real, dimension(:,:,:), pointer :: qq
  integer, dimension(:,:,:), pointer :: cell
  integer :: p(3), wl, ray, pt, ib
  integer, save:: itimer=0
  !---------------------------------------------------------------------
  call trace_begin ('rt_t%update_q', 4, itimer=itimer)
  qq => self%dedt(:,:,:)
  cell => rays%coords
  
  do ray = 1, size(q,1)
    do pt = 2, size(q,3)-1
      p = cell(:,pt-1,ray)
      do wl =1, size(q,2)
        qq(p(1),p(2),p(3)) = &
        qq(p(1),p(2),p(3)) + q(ray,wl,pt)* self%w_lambda(wl) * &
         self%w_bundle(ib) * dwf(ray,pt)
      end do
    end do
  end do
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
  p = self%position - 0.5_8*self%size + (cell-self%li)*self%ds
END FUNCTION cell2pos

!> =====================================================================
!> Convert absolute position to integer index triple 'cell'
!> =====================================================================
FUNCTION pos2cell (self, p) result (cell)
  class(rt_t)  :: self
  real(kind=8) :: p(3)
  integer      :: cell(3)
!-----------------------------------------------------------------------
  cell = nint((p-self%position+0.5_8*self%size)/self%ds)+self%li
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
