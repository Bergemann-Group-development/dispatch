!===============================================================================
!> A particle history is a collection of particles snapshots at differnt times,
!> and may be used for example to compute the forces at other particles at a
!> particular time.  Procedures:
!>
!>  call particle_history%init (n)                : initialize n time slots
!>  call particle_history%add_force_to (particle) : force from virtual particle
!>  call particle_history%rotate                  : rotate time slots
!>  call particle_history%randomize               : for tests only
!>  call particle_history%print (time)            : print particles at time
!===============================================================================
MODULE particle_history_mod
  USE io_mod
  USE trace_mod
  implicit none
  private
  type, public:: particle_t
    integer:: id=0
    real(8):: p(3)=0.0, v(3)=0.0, a(3)=0.0, t=0.0, dt=1.0
  contains
    procedure:: kick
    procedure:: drift
  end type
  type, public:: particle_history_t
    integer:: id = 0
    integer::  n = 6
    integer:: it = 1
    integer:: new= 2
    integer:: moves  = 0
    integer:: forces = 0
    real:: t    = 0.0
    real:: mass = 1.0
    integer         , pointer:: iit(:)
    type(particle_t), pointer:: particle(:)
  contains
    procedure:: init
    procedure:: add_particle_force_to
    procedure:: rotate
    procedure:: randomize
    procedure:: info
    procedure:: delete
  end type
  integer, save:: id=0
CONTAINS

!===============================================================================
SUBROUTINE kick (self)
  class (particle_t):: self
  self%v = self%v + 0.5*self%dt*self%a
END SUBROUTINE kick

!===============================================================================
SUBROUTINE drift (self)
  class (particle_t):: self
  self%p = self%p + self%dt*self%v
  self%t = self%t + self%dt
END SUBROUTINE drift

!===============================================================================
!> Initialize a particle history data type, giving it a unique ID and allocated
!> space for particle position snapshots
!===============================================================================
SUBROUTINE init (self)
  class (particle_history_t):: self
  integer:: i
  !.............................................................................
  !$omp critical (particle_h_cr)
  id = id+1
  self%id = id
  !$omp end critical (particle_h_cr)
  allocate (self%particle(self%n))
  allocate (self%iit(self%n))
  self%iit = self%it
  self%iit(self%n) = self%new
END SUBROUTINE init

!===============================================================================
!> Add the force at time particle%t from a particles whose history is in self
!===============================================================================
SUBROUTINE add_particle_force_to (self, particle)
  class (particle_history_t):: self
  class (particle_t), pointer:: particle
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
      r = (qt*self%particle(i)%p + pt*self%particle(i+1)%p) - particle%p
      ri3 = (1./sqrt(sum(r**2)))**3
      particle%a = particle%a - self%mass*r*ri3
      self%forces = self%forces + 1
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
  self%particle(new)%p  = self%particle(it)%p
  self%particle(new)%v  = self%particle(it)%v
  self%particle(new)%t  = self%particle(it)%t
  self%particle(new)%dt = self%particle(it)%dt
  self%iit(n) = new
END SUBROUTINE rotate

!===============================================================================
!> Randomize the positions and velocities -- for tests
!===============================================================================
SUBROUTINE randomize (self)
  class (particle_history_t):: self
  real:: f(3)
  integer:: i
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_history_t%randomize')
  do i=1,self%n
    call random_number (f)
    self%particle(i)%p = (f-0.5)
    call random_number (f)
    self%particle(i)%v = f
  end do
  call trace%end()
END SUBROUTINE randomize

!===============================================================================
!> Print the particle location and velocity at the (optional) particle time
!===============================================================================
SUBROUTINE info (self, time)
  class (particle_history_t):: self
  real, optional:: time
  real, dimension(self%n):: tt, dt
  integer, dimension(self%n):: iit
  real:: t, r(3), v(3), a(3), pt, qt
  integer:: i, j, n, it
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_history_t%info')
  n   = self%n
  it  = self%it
  if (present(time)) then
    t   = time
    iit = self%iit
    tt  = self%particle(iit)%t
    dt  = self%particle(iit)%dt
    do i=1,n-1
      if (tt(i+1) >= t) then
        pt = (t-tt(i))/dt(i)
        qt = 1.-pt
        r = qt*self%particle(i)%p + pt*self%particle(i+1)%p
        v = qt*self%particle(i)%v + pt*self%particle(i+1)%v
        print '(i6,f10.6,1p,3(2x,3e11.3))', self%id, t, r, v, dt(i)
        exit
      end if
    end do
  else
    t = self%t
    do j=1,3
      r(j) = self%particle(it)%p(j)
      v(j) = self%particle(it)%v(j)
    end do
    print '(i6,f10.6,1p,3(2x,3e11.3))', self%id, t, r, v
  end if
  call trace%end()
END SUBROUTINE info

!===============================================================================
!> Deallocate the particle history array
!===============================================================================
SUBROUTINE delete (self)
  class (particle_history_t):: self
  !.............................................................................
  if (associated(self%particle)) deallocate (self%particle)
  if (associated(self%iit))      deallocate (self%iit)
END SUBROUTINE delete

END MODULE particle_history_mod
