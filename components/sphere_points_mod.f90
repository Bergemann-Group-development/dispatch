!===============================================================================
!> Module for manipulating points on a sphere
!===============================================================================
MODULE sphere_points_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: sphere_points_t
    real(8):: x=0d0, y=0d0, z=0d0
    real(8):: theta=0d0, phi=0d0, r=0d0
  contains
    procedure:: init
    procedure:: to_cartesian
    procedure:: to_spherical
    procedure:: set_cartesian
    procedure:: set_spherical
    procedure:: go_east
    procedure:: go_north
    procedure:: permute
    procedure:: info
    procedure:: is_inside
  end type
CONTAINS

!===============================================================================
!> Initialize a point at the front of a sphere with front:x, east:y, north:z
!===============================================================================
SUBROUTINE init (self)
  class(sphere_points_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%init')
  call self%set_cartesian ([1d0, 0d0, 0d0])
  call trace%end
END SUBROUTINE init

!===============================================================================
!> Set point position using Cartesian cordinates
!===============================================================================
SUBROUTINE set_cartesian (self, p)
  class(sphere_points_t):: self
  real(8):: p(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%set_cartesian')
  self%x = p(1)
  self%y = p(2)
  self%z = p(3)
  call self%to_spherical ()
  call trace%end
END SUBROUTINE set_cartesian

!===============================================================================
!> Set point position using spherical cordinates
!===============================================================================
SUBROUTINE set_spherical (self, p)
  class(sphere_points_t):: self
  real(8):: p(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%set_spherical')
  print *, 'sphere_points_t%set_spherical:', p
  self%theta = p(1)
  self%phi   = p(2)
  self%r     = p(3)
  call self%to_cartesian
  call trace%end
END SUBROUTINE set_spherical

!===============================================================================
!> Convert spherical to Cartesian coordinates
!===============================================================================
SUBROUTINE to_cartesian (self)
  class(sphere_points_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%to_cartesian')
  self%x = self%r*cos(self%theta)*cos(self%phi)
  self%y = self%r*cos(self%theta)*sin(self%phi)
  self%z = self%r*sin(self%theta)
  call trace%end
END SUBROUTINE to_cartesian

!===============================================================================
!> Convert Cartesian to spherical coordinates
!===============================================================================
SUBROUTINE to_spherical (self)
  class(sphere_points_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%to_spherical')
  self%r     = sqrt(self%x**2 + self%y**2 + self%z**2)
  self%theta = asin (self%z/self%r)
  self%phi   = atan2 (self%y, self%x)
  call trace%end
END SUBROUTINE to_spherical

!===============================================================================
!> Take a fraction step east
!===============================================================================
SUBROUTINE go_east (self, dtheta)
  class(sphere_points_t):: self
  real(8):: dtheta
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%go_east')
  self%phi = self%phi + dtheta/cos(self%theta)
  call self%to_cartesian ()
  call trace%end
END SUBROUTINE go_east

!===============================================================================
!> Take a fraction step north
!===============================================================================
SUBROUTINE go_north (self, dtheta)
  class(sphere_points_t):: self
  real(8):: dtheta
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%go_north')
  self%theta = self%theta + dtheta
  call self%to_cartesian ()
  call trace%end
END SUBROUTINE go_north

!===============================================================================
!> Permute indices, corresponding to moving to the east face coordinate system
!===============================================================================
SUBROUTINE permute (self)
  class(sphere_points_t):: self
  real(8):: x
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%permute')
  x = self%x
  self%x = self%y
  self%y = self%z
  self%z = x
  call self%to_spherical ()
  call trace%end
END SUBROUTINE permute

!===============================================================================
!> Check if a patch is fully inside the front face
!===============================================================================
FUNCTION is_inside (self, limit) RESULT (out)
  class(sphere_points_t):: self
  real(8):: limit
  logical:: out
  type(sphere_points_t):: tmp
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%is_inside')
  out = self%theta < limit
  tmp = self
  call tmp%permute ()
  out = self%theta < limit .and. tmp%theta > limit
  call trace%end
END FUNCTION is_inside

!===============================================================================
!> Write out location info
!===============================================================================
SUBROUTINE info (self, label)
  class(sphere_points_t):: self
  character(len=*):: label
  real(8):: x
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_points_t%info')
  write (io_unit%output,'(a,1p,2(2x,3e12.2),3x,a)') 'sphere_points_t%info:', &
    self%x, self%y, self%z, self%theta, self%phi, self%r, trim(label)
  call trace%end
END SUBROUTINE info

END MODULE sphere_points_mod
