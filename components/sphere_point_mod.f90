!===============================================================================
!> Module for manipulating points on a sphere.  Each point is the center of
!> a local coordinate system, where the 1st axis is in the phi-direction, the
!> 2nd axis is in the theta direction, and the 3rd axis is in the radial
!> direction.
!===============================================================================
MODULE sphere_point_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: sphere_point_t
    integer:: id=1
    integer:: verbose=0
    integer:: status=0
    logical:: lower_boundary=.false.
    logical:: upper_boundary=.false.
    real(8):: size=1d0
    real(8):: p(3)
    real(8):: theta=0d0, phi=0d0, r=0d0
    real(8):: e(3,3)                      ! cordinate system unit vectors
    real(8):: north_proxy(3)              
  contains
    procedure:: init
    procedure:: to_cartesian
    procedure:: to_spherical
    procedure:: set_cartesian
    procedure:: set_spherical
    procedure:: go_east
    procedure:: go_north
    procedure:: east_west
    procedure:: north_south
    procedure:: permute
    procedure:: reflect
    procedure:: is_inside
    procedure:: distance_to
    procedure:: unit_vectors
    procedure:: info
    procedure:: get_corner_points
  end type
CONTAINS

!===============================================================================
!> Initialize a point at the front of a sphere with front:x, east:y, north:z
!===============================================================================
SUBROUTINE init (self)
  class(sphere_point_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%init')
  call self%set_cartesian ([1d0, 0d0, 0d0])
  self%north_proxy=[0d0,0d0,1d0]
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Set point position using Cartesian cordinates
!===============================================================================
SUBROUTINE set_cartesian (self, p)
  class(sphere_point_t):: self
  real(8):: p(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%set_cartesian')
  self%p = p
  call self%to_spherical ()
  call trace%end ()
END SUBROUTINE set_cartesian

!===============================================================================
!> Set point position using spherical cordinates
!===============================================================================
SUBROUTINE set_spherical (self, p)
  class(sphere_point_t):: self
  real(8):: p(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%set_spherical')
  print *, 'sphere_point_t%set_spherical:', p
  self%theta = p(1)
  self%phi   = p(2)
  self%r     = p(3)
  call self%to_cartesian ()
  call trace%end ()
END SUBROUTINE set_spherical

!===============================================================================
!> Convert spherical to Cartesian coordinates
!===============================================================================
SUBROUTINE to_cartesian (self)
  class(sphere_point_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%to_cartesian')
  self%p(1) = self%r*cos(self%theta)*cos(self%phi)
  self%p(2) = self%r*cos(self%theta)*sin(self%phi)
  self%p(3) = self%r*sin(self%theta)
  call trace%end ()
END SUBROUTINE to_cartesian

!===============================================================================
!> Convert Cartesian to spherical coordinates
!===============================================================================
SUBROUTINE to_spherical (self)
  class(sphere_point_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%to_spherical')
  self%r     = sqrt(self%p(1)**2 + self%p(2)**2 + self%p(3)**2)
  self%phi   = atan2(self%p(2), self%p(1))
  self%theta = asin(self%p(3)/self%r)
  call trace%end ()
END SUBROUTINE to_spherical

!===============================================================================
!> Take a fraction step east
!===============================================================================
SUBROUTINE go_east (self, dtheta)
  class(sphere_point_t):: self
  real(8):: dtheta
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%go_east')
  self%phi = self%phi + dtheta/cos(self%theta)
  call self%to_cartesian ()
  call self%info ('after go_east')
  call trace%end ()
END SUBROUTINE go_east

!===============================================================================
!> Take a fraction step north
!===============================================================================
SUBROUTINE go_north (self, dtheta)
  class(sphere_point_t):: self
  real(8):: dtheta
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%go_north')
  self%theta = self%theta + dtheta
  call self%info ('after go_north')
  call self%to_cartesian ()
  call trace%end
END SUBROUTINE go_north

!===============================================================================
!> Flip East-West
!===============================================================================
SUBROUTINE east_west (self)
  class(sphere_point_t):: self
  real(8):: x
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%permute')
  self%phi = -self%phi
  call self%to_cartesian ()
  call self%info ('after east_west')
  call trace%end ()
END SUBROUTINE east_west

!===============================================================================
!> Flip North-South
!===============================================================================
SUBROUTINE north_south (self)
  class(sphere_point_t):: self
  real(8):: x
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%permute')
  self%theta = -self%theta
  call self%to_cartesian ()
  call self%info ('after east_west')
  call trace%end ()
END SUBROUTINE north_south

!===============================================================================
!> Permute indices, corresponding to moving to the east face coordinate system
!===============================================================================
SUBROUTINE permute (self)
  class(sphere_point_t):: self
  real(8):: x
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%permute')
  x = self%p(1)
  self%p(1) = self%p(2)
  self%p(2) = self%p(3)
  self%p(3) = x
  x=self%north_proxy(1)
  self%north_proxy(1)=self%north_proxy(2)
  self%north_proxy(2)=self%north_proxy(3)
  self%north_proxy(3)=x
  call self%to_spherical ()
  call self%info ('after permute')
  call trace%end ()
END SUBROUTINE permute

!===============================================================================
!> Move to the VB system on the opposite side
!===============================================================================
SUBROUTINE reflect (self)
  class(sphere_point_t):: self
  real(8):: x
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%reflect')
  call self%info ('before reflect')
  self%p = -self%p
  call self%to_spherical ()
  call self%info ('after reflect')
  call trace%end ()
END SUBROUTINE reflect

!===============================================================================
!> Get the 8 points, representing the corner points of a patch, in the correct
!> orientation
!===============================================================================
SUBROUTINE get_corner_points(self,corners)
  class(sphere_point_t):: self
  class(sphere_point_t), pointer:: corners(:), this
  ! ----------------------------------------------------------------------------
  real(8) :: pos(3), c(2,3), np(3)
  integer:: i,j,k,ic
  ! ............................................................................
  call self%to_cartesian()   ! just in case
  pos = self%p
  call self%unit_vectors()
  c(1,:) = -0.5d0*self%size
  c(2,:) = +0.5d0*self%size
  ic=1
  do k=1,2
    do j=1,2
      do i=1,2
        this => corners(ic)
        this%p=matmul([c(i,1),c(j,2),c(k,3)],self%e)+pos
        call this%to_spherical()
        ic=ic+1
      end do
    end do
  end do

END SUBROUTINE get_corner_points

!===============================================================================
!> Check if a patch is fully inside the front face, by checking that the patch
!> is inside the latitude limit of the central VB section, and that its lower
!> south west corner is outside the east VB section
!===============================================================================
FUNCTION is_inside (self, theta0, dtheta) RESULT (in)
  class(sphere_point_t):: self
  real(8):: theta0, theta, theta1, dtheta
  logical:: in, inside1(8), inside2(8)
  integer:: i
  type(sphere_point_t):: tmp
  class(sphere_point_t), pointer:: corners(:)
  !-----------------------------------------------------------------------------
  ! Check if inside latitude theta0
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%is_inside')
  theta = self%theta
  in = theta < theta0
  tmp = self
  allocate (corners(8))
  call self%get_corner_points(corners)
  do i=1,8
    inside1(i)=corners(i)%theta<=theta0
  end do
  do i=1,8
    call corners(i)%permute()
    inside2(i)=corners(i)%theta >= theta0
  end do
  !-----------------------------------------------------------------------------
  ! Move point to south-west lower corner before converting to east VB coords
  !-----------------------------------------------------------------------------
  !tmp%theta = tmp%theta - dtheta/2d0
  !tmp%phi   = tmp%phi   - dtheta/(2d0*cos(tmp%theta))
  !tmp%r     = tmp%r*(1d0- dtheta/2d0)
  !call tmp%permute ()
  !theta1 = tmp%theta
  !-----------------------------------------------------------------------------
  ! Check if outside east VB section latitude limits
  !-----------------------------------------------------------------------------
  in = theta < theta0 .and. any(inside1) .and. any(inside2)
  deallocate (corners)
  call trace%end ()
END FUNCTION is_inside

!===============================================================================
!> Distance between position points in a sphere
!===============================================================================
FUNCTION distance_to (self, task) RESULT (out)
  class(sphere_point_t):: self, task
  real(8):: out(3)
  !-----------------------------------------------------------------------------
  out = task%p - self%p
END FUNCTION distance_to

!===============================================================================
!> Define the unit vectors of a coordinate system where x1 is East, x2 is North,
!> and x3 is Up:
!===============================================================================
SUBROUTINE unit_vectors (self)
  class(sphere_point_t):: self
  integer:: i1, i2, i3, i(1)
  real(8):: p(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sphere_point_t%unit_vectors')
  i  = maxloc(abs(self%p))                                    ! decide radial
  i3 = i(1)                                                   ! local Up
  i1 = mod(i3,3)+1                                            ! local East
  i2 = mod(i1,3)+1                                            ! local North
  i1=1;i2=2;i3=3
  !p(1) = self%p(i1)
  !p(2) = self%p(i2)
  !p(3) = self%p(i3)
  p=self%p
  self%e(:,i3) = normalize(p)                                 ! unit vector up
  self%e(:,i2) = self%north_proxy
  !self%e(:,i2) = [0d0,1d0,0d0]                                ! north proxy
  self%e(:,i1) = normalize(cross(self%e(:,i2),self%e(:,i3)))  ! unit vector east
  self%e(:,i2) = normalize(cross(self%e(:,i3),self%e(:,i1)))  ! unit vector north
  call trace%end ()
END SUBROUTINE unit_vectors

!============================i===================================================
!> Write out location info
!===============================================================================
SUBROUTINE info (self, label, verbose)
  class(sphere_point_t):: self
  character(len=*):: label
  integer, optional:: verbose
  real(8):: x
  integer:: i, verbose_l
  !-----------------------------------------------------------------------------
  verbose_l = self%verbose
  if (present(verbose)) verbose_l = verbose
  if (verbose_l > 4) then
    call trace%begin ('sphere_point_t%info')
    write (stdlog,'("patch id",i3)') self%id
    write (stdlog,'(a,1p,2(2x,3e12.2),3x,a)') &
      'sphere_point_t%info: p, theta,phi,r =', &
      self%p, self%theta, self%phi, self%r, trim(label)
    do i=1,3
      write (stdlog,'("sphere_point_t%info: e",i1,3f10.6)') i, self%e(:,i)
    end do
    call trace%end
  end if
END SUBROUTINE info

!===============================================================================
!> Utility functions for single 3-D vector operations
!===============================================================================
FUNCTION dot (a, b) RESULT (c)
  real(8):: a(3), b(3), c
  c = sum(a*b)
END FUNCTION dot
!===============================================================================
FUNCTION len (a) RESULT (c)
  real(8):: a(3), c
  c = sqrt(dot(a,a))
END FUNCTION len
!===============================================================================
FUNCTION normalize (a) RESULT (b)
  real(8):: a(3), b(3)
  b = a/sqrt(dot(a,a))
END FUNCTION normalize
!===============================================================================
FUNCTION cross (a, b) RESULT (c)
  real(8):: a(3), b(3), c(3)
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
END FUNCTION cross

END MODULE sphere_point_mod
