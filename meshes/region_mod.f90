!===============================================================================
!> sources are defined by their lower (fl) and upper (fu) floating point
!> coordinates in Cartesian coordinates. Important examples of sources:
!>
!> * The internal source of a patch
!> * The guard zone that covers one face of a patch
!> * The guard zone that touches an edge of a patch
!> * The guard zone that touches a corner of a patch
!> * The intersection of the internal source of a source patch and a target
!>   guard zone
!>
!> A neighbor patch needs to be in the nbor list of a patch if one of its
!> guard zones is covered by the internal source of the neighbor patch.
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
    procedure:: contains
  end type region_t
  real(8), public, save:: region_margin = 1d-10
CONTAINS

!===============================================================================
!> Return the intersection of source of target
!===============================================================================
FUNCTION intersection (target, source, distance) RESULT (c)
  class(region_t) :: target, source
  real(8)         :: distance(3)
  type(region_t)  :: c
  !-----------------------------------------------------------------------------
  c%l = max(target%l,source%l+distance)
  c%u = min(target%u,source%u+distance)
END FUNCTION intersection

!===============================================================================
!> Check if two sources overlap; i.e. if the volume of the intersection is
!> positive
!===============================================================================
FUNCTION overlaps (target, source, distance)
  class(region_t) :: target, source
  real(8)         :: distance(3), vi, vs
  logical         :: overlaps
  type(region_t)  :: tmp
  !-----------------------------------------------------------------------------
  vi = volume (intersection (target, source, distance))
  vs = volume (source)
  overlaps =  vi  >  vs*region_margin
END FUNCTION overlaps

!===============================================================================
!> Check if target covers source; i.e. if the source is completely immersed
!===============================================================================
FUNCTION covers (target, source, distance)
  class(region_t) :: target, source
  real(8)         :: distance(3), vi, vs
  logical         :: covers
  !-----------------------------------------------------------------------------
  vi = volume (intersection (target, source, distance))
  vs = volume (source)
  covers = abs(vi-vs) < vs*region_margin 
END FUNCTION covers

!===============================================================================
!> Return the volume of a source
!===============================================================================
FUNCTION volume (target) RESULT (v)
  class(region_t) :: target
  real(8)         :: v
  !-----------------------------------------------------------------------------
  v = product (max (target%u-target%l, 0d0))
END FUNCTION volume

!===============================================================================
!> Return the volume of a source
!===============================================================================
FUNCTION empty_relative_to (target, compare) RESULT (y)
  class(region_t) :: target, compare
  logical         :: y
  !-----------------------------------------------------------------------------
  y = target%volume() < compare%volume()*region_margin
END FUNCTION empty_relative_to

!===============================================================================
!> Return true if region contains point
!===============================================================================
FUNCTION contains (self, point)
  class(region_t) :: self
  real(8)         :: point(3)
  logical         :: contains
  !-----------------------------------------------------------------------------
  contains = all(point >= self%l .and. point <= self%u)
END FUNCTION contains

END MODULE region_mod
