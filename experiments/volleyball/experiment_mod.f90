!===============================================================================
!> Minimal experiment initialization -- all other calls fall through to solver
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE link_mod
  USE task_mod
  USE solver_mod
  USE sphere_points_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    type(sphere_points_t):: place
    real(8), dimension(3):: e1, e2, e3
  contains
    procedure:: init
    procedure:: dnload
    procedure:: distance
  end type
CONTAINS

!===============================================================================
!> Generic init procedure
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i1, i2, i3, l(3), u(3)
  real(8):: r1, r2, r3, r, h
  real, save:: d0=1.0, e0=1.0, g=1.0, r0=1d0
  logical, save:: first_time=.true.
  namelist /experiment_params/ d0, e0, r0, g
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%init')
  call self%solver_t%init
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params)
    write (io%output, experiment_params)
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Setup local coordinate system, with unit vectors
  !-----------------------------------------------------------------------------
  self%position = [self%place%x, self%place%y, self%place%z]
  self%e1 = normalize(self%position)
  if (self%position(3) > 0d0) then
    self%e2 = [0d0,0d0,1d0] - self%position
  else
    self%e2 = self%position - [0d0,0d0,-1d0]
  end if
  self%e3 = normalize(cross(self%e1,self%e2))
  self%e2 = cross(self%e3,self%e1)
  !-----------------------------------------------------------------------------
  ! Set up a simple, isothermal, hydrostatic initial state.  The formulation
  ! assumes that each patch has a coodinate system initialized with the z-
  ! direction vertical, and with the self%mesh(3)%p giving the distance of the
  ! midpoint to the center of the sphere.
  !-----------------------------------------------------------------------------
  if (io%restart < 0) then
    h = (self%gamma-1d0)*e0/(d0*g)
    l = self%mesh%lb
    u = self%mesh%ub
    associate (m1=>self%mesh(1), m2=>self%mesh(2), m3=>self%mesh(3))
    m1%p = 0d0
    m2%p = 0d0
    self%position = [self%place%x, self%place%y, self%place%z]
    m3%p = sqrt(sum(self%position**2))
    do i3=l(3),u(3)
      r3 = m3%p + m3%r(i3)
      do i2=l(2),u(2)
        r2 = m2%r(i2)
        do i1=l(1),u(1)
          r1 = m1%r(i1)
          r = sqrt(r1**2 + r2**2 + r3**2) - r0
          self%mem(i1,i2,i3,self%idx%d ,self%it,1) = d0*exp(-r/h)
          self%mem(i1,i2,i3,self%idx%e ,self%it,1) = e0*exp(-r/h)
          self%mem(i1,i2,i3,self%idx%px,self%it,1) = 0.0
          self%mem(i1,i2,i3,self%idx%py,self%it,1) = 0.0
          self%mem(i1,i2,i3,self%idx%pz,self%it,1) = 0.0
if (i1==8.and.i2==8) &
print *, self%id, i1, i2, i3, m3%p, r, d0*exp(-r/h)
        end do
      end do
    end do
    end associate
  end if
  self%unsigned(self%idx%d) = .true.
  self%unsigned(self%idx%e) = .true.
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Overload the standard dnload function
!===============================================================================
SUBROUTINE dnload (self, only)
  class(experiment_t):: self
  class(link_t), pointer:: nbor
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%dnload')
  nbor => self%link%nbor
  do while (associated(nbor))
    write (stdout,*) 'experiment_t%dnload: nbor', self%id, nbor%task%id, &
      sqrt(sum((self%position-nbor%task%position)**2))
    nbor => nbor%next
  end do
  call self%solver_t%dnload (only)
  call trace%end()
END SUBROUTINE dnload

!===============================================================================
!> Overload the standard distance function
!===============================================================================
FUNCTION distance (self, task) RESULT (out)
  class(experiment_t):: self
  class(task_t):: task
  real(8):: out(3), position(3)
  !-----------------------------------------------------------------------------
  position = self%position - task%position
  out = [dot(position,self%e1), &
         dot(position,self%e2), &
         dot(position,self%e3)]
print *, out, 'distance'
END FUNCTION distance

!===============================================================================
!> Utility functions for single 3-D vector operations
!===============================================================================
FUNCTION dot (a, b) RESULT (c)
  real(8):: a(3), b(3), c
  c = sum(a*b)
END FUNCTION dot
!===============================================================================
FUNCTION normalize (a) RESULT (b)
  real(8):: a(3), b(3)
  b = a/sqrt(sum(a**2))
END FUNCTION normalize
!===============================================================================
FUNCTION cross (a, b) RESULT (c)
  real(8):: a(3), b(3), c(3)
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
END FUNCTION cross

END MODULE experiment_mod
