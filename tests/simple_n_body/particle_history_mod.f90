!===============================================================================
!> A particle history is a collection of particles snapshots at differnt times,
!> and may be used for example to compute the forces at other particles at a
!> particular time.
!===============================================================================
MODULE particle_history_mod
  USE trace_mod
  USE particle_mod
  implicit none
  private
  type, public:: particle_history_t
    integer:: id=0
    real:: mass=1.0
    integer:: n=6
    integer:: it=1, new=2
    integer:: moves=0
    real:: t=0.0
    integer, pointer:: iit(:)
    type (particle_t), pointer:: particle(:)
    type(particle_history_t), pointer:: next => null()
  contains
    procedure:: init
    procedure:: add_particle_force_to
    procedure:: rotate
    procedure:: randomize
    procedure:: print
  end type
  integer, save:: id=0
CONTAINS

!===============================================================================
!> Initialize a particle history data type, giving it a unique ID and allocated
!> space for particle snapshots
!===============================================================================
SUBROUTINE init (self)
  class (particle_history_t):: self
  !.............................................................................
  !$omp critical (particle_h_cr)
  id = id+1
  self%id = id
  !$omp end critical (particle_h_cr)
  allocate (self%particle(self%n))
  self%particle%id = id
  allocate (self%iit(self%n))
  self%iit = self%it
  self%iit(self%n) = self%new
END SUBROUTINE init

!===============================================================================
!> Add the force at time t from a particles whose history is in self
!===============================================================================
SUBROUTINE add_particle_force_to (self, particle)
  class (particle_history_t):: self
  type (particle_t), pointer:: particle
  !.............................................................................
  real, dimension(self%n):: tt, dt
  integer, dimension(self%n):: iit
  real:: pt, qt, r(3), ri3
  integer:: i
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  !call trace%begin ('particle_history_t%add_force_to', itimer=itimer)
  iit = self%iit
  tt = self%particle(iit)%t
  dt = self%particle(iit)%dt
  do i=1,self%n-1
    if (tt(i+1) > particle%t) then
      pt = (particle%t-tt(i))/dt(i)
      qt = 1.0-pt
      r = (qt*self%particle(i)%r + pt*self%particle(i+1)%r) - particle%r
      ri3 = (1./sqrt(sum(r**2)))**3
      particle%a = particle%a - self%mass*r*ri3
      return
    end if
  end do
  !call trace%end (itimer)
END SUBROUTINE add_particle_force_to

!===============================================================================
!> Rotate the time slots in the particle_history data type
!===============================================================================
SUBROUTINE rotate (self)
  class (particle_history_t):: self
  integer:: i, n, it, new
  !-----------------------------------------------------------------------------
  n = self%n
  do i=1,n-1
    self%iit(i) = self%iit(i+1)
  end do
  it = self%iit(n-1)
  new = mod(it,n)+1
  self%iit(n) = new
END SUBROUTINE rotate

!===============================================================================
!> Randomize the velocities -- for tests
!===============================================================================
SUBROUTINE randomize (self)
  class (particle_history_t):: self
  real:: v(3)
  integer:: i
  !-----------------------------------------------------------------------------
  do i=1,self%n-1
    call random_number (v)
    self%particle(i)%v = v
  end do
END SUBROUTINE randomize

!===============================================================================
!> Print the particle location and velocity at the particle time
!===============================================================================
SUBROUTINE print (self, time)
  class (particle_history_t):: self
  real, optional:: time
  real, dimension(self%n):: tt, dt
  integer, dimension(self%n):: iit
  real:: t, r(3), v(3), a(3), pt, qt
  integer:: i, n
  !-----------------------------------------------------------------------------
  if (present(time)) then
    t = time
  else
    t = self%t
  end if
  iit = self%iit
  tt  = self%particle(iit)%t
  dt  = self%particle(iit)%dt
  n   = self%n
  do i=1,n-1
    if (tt(i+1) >= t) then
      pt = (t-tt(i))/dt(i)
      qt = 1.-pt
      r = qt*self%particle(i)%r + pt*self%particle(i+1)%r
      v = qt*self%particle(i)%v + pt*self%particle(i+1)%v
      print '(i6,f10.6,1p,3(2x,3e11.3))', i, t, r, v, self%particle(n)%dt
      return
    end if
  end do
END SUBROUTINE print

END MODULE particle_history_mod
