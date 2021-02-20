!===============================================================================
!> Regions are defined by their lower (fl) and upper (fu) floating point
!> coordinates in Cartesian coordinates. Important examples of regions:
!>
!> * The internal region of a patch
!> * The guard zone that covers one face of a patch
!> * The guard zone that touches an edge of a patch
!> * The guard zone that touches a corner of a patch
!> * The intersection of the internal region of a source patch and a target
!>   guard zone
!>
!> A neighbor patch needs to be in the nbor list of a patch if one of its
!> guard zones is covered by the internal region of the neighbor patch.
!> Conversely, the neighbor patch needs to have the patch in its nbor list,
!> with the nbor%needs_me flag set.
!===============================================================================
MODULE region_mod
  USE io_unit_mod
  implicit none
  private
  type, public:: region_t
    real(8), dimension(3):: l=0d0, u=0d0
  contains
    procedure:: volume
    procedure:: empty_relative_to
    procedure:: intersection
    procedure:: covers
    procedure:: overlaps
  end type region_t
  real(8), public, save:: region_margin = 1d-10
CONTAINS

!===============================================================================
!> Return the common region of self and other
!===============================================================================
FUNCTION intersection (self, other) RESULT (c)
  class(region_t) :: self, other
  type(region_t)  :: c
  !-----------------------------------------------------------------------------
  c%l = max(self%l,other%l)
  c%u = min(self%u,other%u)
END FUNCTION intersection

!===============================================================================
!> Check if two regions overlap; i.e. if the volume of the intersection is
!> positive
!===============================================================================
FUNCTION overlaps (self, region)
  class(region_t) :: self, region
  logical         :: overlaps
  type(region_t)  :: tmp
  !-----------------------------------------------------------------------------
!tmp = intersection(self,region)
!write(stdlog,'(a,3(1x,2(1x,3f8.4)),1p,2e11.2)') 'overlaps:', &
!self%l, self%u, region%l, region%u, tmp%l, tmp%u, &
!volume(intersection(self,region)), volume(region)
  overlaps =  volume(intersection(self,region)) > volume(region)*region_margin
END FUNCTION overlaps

!===============================================================================
!> Check if self covers region; i.e. if the region is completely immersed
!===============================================================================
FUNCTION covers (self, region)
  class(region_t) :: self, region
  logical         :: covers
  !-----------------------------------------------------------------------------
  covers = abs(volume(intersection(self,region))/volume(region) - 1d0) &
         < region_margin 
END FUNCTION covers

!===============================================================================
!> Return the volume of a region
!===============================================================================
FUNCTION volume (self) RESULT (v)
  class(region_t) :: self
  real(8)         :: v
  !-----------------------------------------------------------------------------
  v = product(max(self%u-self%l,0d0))
END FUNCTION volume

!===============================================================================
!> Return the volume of a region
!===============================================================================
FUNCTION empty_relative_to (self, compare) RESULT (y)
  class(region_t) :: self, compare
  logical         :: y
  !-----------------------------------------------------------------------------
  y = self%volume() < compare%volume()*region_margin
END FUNCTION empty_relative_to

END MODULE region_mod