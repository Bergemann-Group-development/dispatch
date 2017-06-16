MODULE rt_components_mod
  USE link_mod
  USE ray_geometry_v3_mod
  implicit none

!> =====================================================================
!> intensity/heating rate data type is used to store intensities in the outer shell of
!> the patch as well as quickly retrieve incoming intensities from nbors
!> =====================================================================
  type, public :: q_data_t
    real, dimension(:), pointer :: in  => null()                        ! (wl_bin) incoming intensity.
    real, dimension(:), pointer :: out => null()                        ! leaving intensity
    integer, pointer :: dir_ID
    !> Subsequent pointers are used to handle upstreams with different refinement level
    type (q_data_t), pointer :: head => null()
    type (q_data_t), pointer :: next => null()
    type (q_data_t), pointer :: tail => null()
    integer :: n_upstream = 0
    real :: weight = 1.
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
    class(q_data_t), dimension(:), pointer :: packet_data => null()     ! quick pointing to the right locations in the data sets
    class(ray_bndl_t),             pointer :: bundle => null()          ! ray bundle of this chain ring
    logical, dimension(:), pointer :: nudged                            ! which rings have already nudged me
    integer :: direction                                                ! 1 = forward / 2 = reverse
    integer :: id                                                       ! location id in the ring array of a particular patch
    logical :: master = .false.                                         ! Am I the ring to rule them all?
    logical :: BC(3) = .false.                                          ! outer/inner/virtual or level
    integer :: l_chain = 0                                              ! length of this chain, points to one number, set on the head of the chain
    integer :: n_ready                                                  ! number of ready patches in this chain. When n_ready == l_chain, the chain can be executed.
    integer :: chain_nm                                                 ! ring number in the chain
    real(kind=8) :: ds                                                  ! distance between two ray points
  contains
    procedure :: add_upstream
    procedure :: check_masters
    procedure :: finish_chain
    procedure :: get_ds
    procedure :: get_dnstream_coord
    procedure :: get_upstream_coord
  end type

CONTAINS

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

!> =====================================================================
!> Get upstream cell coordinate for a particular ray
!> =====================================================================
FUNCTION get_upstream_coord(self, pck, ray) result(coords)
 class(chain_ring_t) :: self
 integer :: pck, ray, coords(3)
 !----------------------------------------------------------------------
 select case(self%direction)
 case(1)
   coords = self%bundle%packet(pck)%rev_coords(:,ray)
 case(2)
   coords = self%bundle%packet(pck)%fwd_coords(:,ray)
 end select
END FUNCTION get_upstream_coord

!> =====================================================================
!> Get downstream cell coordinate for particular ray
!> =====================================================================
FUNCTION get_dnstream_coord(self, pck, ray) result(coords)
 class(chain_ring_t) :: self
 integer :: pck, ray, coords(3)
 !----------------------------------------------------------------------
 select case(self%direction)
 case(1)
   coords = self%bundle%packet(pck)%fwd_coords(:,ray)
 case(2)
   coords = self%bundle%packet(pck)%rev_coords(:,ray)
 end select
END FUNCTION get_dnstream_coord

!> =====================================================================
!> Add a pointer to intensity data in an upstream patch. This is different
!> than just a generic upstrean chain ring, since there can be more than
!> one source of data with varying statistical weights.
!> =====================================================================
SUBROUTINE add_upstream(self, ir, ups, iw)
  class(chain_ring_t) :: self
  integer :: ir
  class(q_data_t), pointer :: ups, pd
  real :: iw
  !---------------------------------------------------------------------
  if (.not.associated(self%packet_data(ir)%head)) then
    pd => self%packet_data(ir)
    self%packet_data(ir)%head => pd
    self%packet_data(ir)%tail => pd
  else
    allocate(pd)
    self%packet_data(ir)%tail%next => pd
    self%packet_data(ir)%tail => pd
  end if  
  pd%head => self%packet_data(ir)
  self%packet_data(ir)%n_upstream = self%packet_data(ir)%n_upstream + 1
  pd%weight = iw
  pd%in => ups%out
  pd%dir_ID => ups%dir_ID

END SUBROUTINE add_upstream

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
  class(q_data_t), dimension(:,:,:), pointer:: wall
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
      if(.not.ring%BC(1).and..not.ring%BC(3)) then                      ! avoid master rings, whose incoming intensities are from the boundaries or upstreams
      
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


SUBROUTINE get_ds(self, cs, ds)
  class(chain_ring_t) :: self
  real(kind=8) :: cs(3), ds
  !.....................................................................
  real(kind=8) :: dst(3)
  !---------------------------------------------------------------------
  dst = real(self%bundle%dstep,kind=8)
  
  ds = ((cs(1)*dst(1))**2+(cs(2)*dst(2))**2+(cs(3)*dst(3))**2)**0.5
  
END SUBROUTINE get_ds

END MODULE rt_components_mod
