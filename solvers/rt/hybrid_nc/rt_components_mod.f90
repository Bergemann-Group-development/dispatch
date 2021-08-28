MODULE rt_components_mod
  USE link_mod
  USE patch_mod
  USE ray_geometry_v3_mod
  USE eos_mod
  implicit none
!=======================================================================
!> intensity/heating rate data type is used to store intensities in the outer shell of
!> the patch as well as quickly retrieve incoming intensities from nbors
!=======================================================================
  type, public :: rt_packet
    real, dimension(:,:), pointer :: in  => null()                        ! (wl_bin|time slot) incoming intensity.
  end type
!=======================================================================
!> rt_chain_t is used to make a linked list of chain link. Each rt_task
!> owns a set of links, related to particular direction. Each link in 
!> the set is then connected to a downstream nbor. 
!=======================================================================
  type, public :: chain_link_t
    class(link_t),         pointer :: link   => null()                  ! link to rt_t task
    class(chain_link_t),   pointer :: head   => null()                  ! pointer to the master chain link
    class(chain_link_t),   pointer :: tail   => null()                  ! pointer to the last link
    class(chain_link_t),   pointer :: next   => null()
    class(chain_link_t),   pointer :: prev   => null()
    class(chain_link_t),   pointer :: temp   => null()                  ! temporal chain connection to do any required modifications
    class(packet_t),       pointer :: packet => null()                  ! ray packet of this chain link
    class(rt_packet), dimension(:), pointer :: pdata => null()          ! quick pointing to the right locations in the data sets
    logical, dimension(:), pointer :: nudged => null()                  ! which links have already nudged me
    integer :: direction                                                ! 1 = forward / 2 = reverse
    integer :: dir_id
    integer :: id                                                       ! location id in the link array of a particular patch
    logical :: master = .false.                                         ! Am I the link (ring) to rule them all?
    logical :: BC(3) = .false.                                          ! outer/inner/virtual or level
    integer :: l_chain = 0                                              ! length of this chain, points to one number, set on the head of the chain
    integer :: n_ready                                                  ! number of ready patches in this chain. When n_ready == l_chain, the chain can be executed.
    integer :: chain_nm                                                 ! link number in the chain
    integer :: ivl, ivu, ivl2, ivu2
    real(kind=8) :: ds                                                  ! distance between two ray points
  contains
    procedure :: add_upstream
    procedure :: check_masters
    procedure :: finish_chain
    procedure :: get_dnstream_coord
    procedure :: get_upstream_coord_same
  end type

CONTAINS

!=======================================================================
!> to run this, disable the deallocation in finish_chain subroutine and
!> uncomment the last line there
!=======================================================================
SUBROUTINE check_masters(self)
  class(chain_link_t) :: self
  class(chain_link_t), pointer :: link
  !---------------------------------------------------------------------
  link => self%prev
  do while(associated(link))
  
  print *,'link_l_chain', link%l_chain
  print *,'associated nudged',associated(link%nudged)
  print *,'nudged size',size(link%nudged)
  
  link => link%prev
  end do
END SUBROUTINE check_masters

!=======================================================================
!> Get upstream cell coordinate for a particular ray
!> Here a modificator should come, which decides the position, depending
!> whether the no mans land is true or false.
!=======================================================================
FUNCTION get_upstream_coord_same(self, ray) result(coords)
 class(chain_link_t) :: self
 class(patch_t), pointer :: patch
 integer :: pck, ray, coords(3), i
 !----------------------------------------------------------------------
 patch => task2patch(self%link%task)
 select case(self%direction)
 case(1)
   coords = self%packet%rev_coords(:,ray)
 case(2)
   coords = self%packet%fwd_coords(:,ray)
 case default
   error stop "chain_link_t::get_upstream_coords:: direction not understood"
 end select
 do i=1,3
   if (coords(i) > patch%ui(i)) then
     coords(i) = coords(i) - patch%n(i)
   else if (coords(i) < patch%li(i)) then
     coords(i) = coords(i) + patch%n(i)
   end if
 end do
END FUNCTION get_upstream_coord_same

!=======================================================================
!> Get downstream cell coordinate for particular ray
!=======================================================================
FUNCTION get_dnstream_coord(self, ray) result(coords)
 class(chain_link_t) :: self
 integer :: pck, ray, coords(3)
 !----------------------------------------------------------------------
 select case(self%direction)
 case(1)
   coords = self%packet%fwd_coords(:,ray)
 case(2)
   coords = self%packet%rev_coords(:,ray)
 case default
   error stop "chain_link_t::get_dnstream_coords:: direction not understood"
 end select
END FUNCTION get_dnstream_coord

!=======================================================================
!> Add a pointer to intensity data in an upstream patch. This is different
!> than just a generic upstrean chain link, since there can be more than
!> one source of data with varying statistical weights.
!=======================================================================
SUBROUTINE add_upstream(self, ir, ups, iw)
  class(chain_link_t) :: self
  integer :: ir
  class(rt_packet), pointer :: ups, pd
  real :: iw
!  !---------------------------------------------------------------------
!  if (.not.associated(self%packet_data(ir)%head)) then
!    pd => self%packet_data(ir)
!    self%packet_data(ir)%head => pd
!    self%packet_data(ir)%tail => pd
!  else
!    allocate(pd)
!    self%packet_data(ir)%tail%next => pd
!    self%packet_data(ir)%tail => pd
!  end if  
!  pd%head => self%packet_data(ir)
!  self%packet_data(ir)%n_upstream = self%packet_data(ir)%n_upstream + 1
!  pd%weight = iw
!  pd%in => ups%out
!  pd%dir_ID => ups%dir_ID

END SUBROUTINE add_upstream

!=======================================================================
!> Connect the links to the appropriate incoming net heating rates from the
!> right upstreams. 
!> INPUT: a temporary link, which has all the master link in its chain
!> RESULT: all the links in all the chains know, where to get incoming
!>   intensities from.
!> FLOW: 
!> 1) first item is ommitted, since it is just a placeholder. All the 
!>   other items are oranised in "reverse", since master links do/should
!>   not contain a %prev member.
!> 2) a loop over all master links is made.
!> 3) a loop over all links in this chain, belonging to a particular
!>   master is done
!> 4) master gets information about the chain length. This is utterly
!>   important later.
!> 5) link direction is checked and based on that the right ray arrangement
!>   in upstream link is determined
!> 6) depending on what is the periodicity of chain alterations for this
!>   orientation, one determines, where the rays are (a supplement for (5))
!> 7) the appropriate ray in the link is connected to the outer shell of
!>   the upstream patch.
SUBROUTINE finish_chain(self)
  class(chain_link_t) :: self
  !.....................................................................
  class(chain_link_t), pointer :: this, prev
  class(chain_link_t), pointer :: link, plink
  class(patch_t),      pointer :: patch
  integer :: cn, ir, cc(3)
  integer, dimension(:), pointer :: nb_ray_pos
  !---------------------------------------------------------------------
  this => self%prev                                                     ! first item is not needed
  do while(associated(this))                                            ! loop over master links
    prev => this%prev                                                   ! this temporary chain is organised in reverse
    nullify(this%prev)                                                  ! free the %prev pointer
    link => this
    cn=1
    do while(associated(link))
      link%head => this
      this%l_chain = this%l_chain + 1                                   ! increment the chain length (number of members)
      link%chain_nm = cn                                                ! link's list number in the chain
      cn = cn+1
      if(.not.link%BC(1).and..not.link%BC(3)) then                      ! avoid master links, whose incoming intensities are from the boundaries or upstreams
        patch => task2patch(plink%link%task)
        do ir = 1, link%packet%n_rays
          cc = link%get_upstream_coord_same(ir)
          link%pdata(ir)%in => patch%mem(cc(1),cc(2),cc(3),link%ivl:link%ivu,:,1)
        end do
      end if
      plink => link                                                     ! this will now be the previous link
      link => link%next                                                 ! move to the next link
    end do

    allocate(this%nudged(this%l_chain))                                 ! allocate the nudging logical
    this%nudged(:) = .false.                                            ! nothing is nudged yet.
    this => prev                                                        ! point to next master link
  end do
  !if (io%verbose>1) call self%check_masters
END SUBROUTINE finish_chain

END MODULE rt_components_mod
