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
  type, public :: omega_t
    class(link_t),         pointer :: link   => null()                          ! link to rt_t task
    class(ray_bndl_t),     pointer :: bundle => null()                          ! ray packet of this omega
    class(rt_packet), dimension(:), pointer :: pdata => null()                  ! quick pointing to the right locations in the data sets
    integer :: direction                                                        ! 0 - bidirectional, 1 = forward, -1 = reverse
    integer :: id                                                               ! location id in the omega array of a particular patch
    integer :: ifl, ifu, irl, iru                                               ! lower/higher memory indices for forward/reverse ray directions
    real(kind=8) :: ds                                                          ! distance between two ray points
    real    :: weight                                                           ! statistical weight for this ray direction
  end type

CONTAINS
END MODULE rt_components_mod
