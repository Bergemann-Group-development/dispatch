!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Demonstrate particle update method, and measure update speed
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE trace_mod
  USE random_mod
  USE task_mod
  USE star_mod
  USE particle_mod
  implicit none
  private
  type, public, extends(task_t):: particles_t
    type(particle_t), pointer:: head => null()    ! particle list head
    integer:: n
    integer:: n_step
    real:: radius
    real:: dispersion
    real:: dt_min
    real:: dt_max
    real:: courant
  contains
    procedure:: init
    procedure:: update
    procedure:: output
  end type
  type(random_t):: random
  type(star_t):: star
CONTAINS

!> =============================================================================
!> Initialize a number of particles consistent with a particle task, with a given
!> total mass, in Kepler orbits with a specified extra velocity dispersion.
!> =============================================================================
SUBROUTINE init (self)
  class(particles_t):: self
  integer:: n=10000, i
  type(particle_t), pointer:: p
  real:: radius=1.
  real:: dispersion=0.1
  real:: mass=0.01
  real:: dtime=1.0
  real:: end_time=100.
  real:: out_time=1.
  real:: courant=0.01
  logical, save:: first_time=.true.
  namelist /particle_params/ n, radius, mass, dispersion, dtime, courant, &
    end_time, out_time
  !.............................................................................
  call trace_begin ('particles_mod%init', 3)
  call star%init
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, particle_params)
    if (io%master) write (*, particle_params)
  end if
  self%n = n
  self%n_step = 0
  self%radius = radius
  self%dispersion = dispersion
  self%time = 0.0
  self%dtime = dtime
  self%end_time = end_time
  self%out_time = out_time
  self%courant = courant
  do i=1,n
    allocate (p)
    p%next => self%head
    self%head => p
    p%r = random%ran3()
    p%v = star%kepler_speed(p%r) * (1.+dispersion*random%ran3())
    p%m = mass/n
  end do
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Update particles using the Kick-Drift-Kick algorithm (Springel, 19xx).  The
!> particles are advanced individually, but for a maximum time of self%dtime,
!> which is intended to correspond to a typical patch update time.
!> =============================================================================
SUBROUTINE update (self)
  class(particles_t):: self
  real:: a(3)                                   ! force
  real:: dt                                     ! time step
  type(particle_t), pointer:: p
  !.............................................................................
  call trace_begin ('particles_mod%update', 3)
  p => self%head
  self%time = self%time + self%dtime            ! target time
  do while (associated(p%next))
    do while (p%t < self%time)
      a = star%central_force (p%r)              ! acceleration
      dt = sqrt(sum(p%v**2)/sum(a**2)) &        ! ratio of speed to acceleration
           * self%courant                       ! particle Courant number
      dt = min(self%time - p%t, dt)             ! don't overshoot
      p%v = p%v + a*dt*0.5                      ! kick
      p%r = p%r + p%v*dt                        ! drift
      a = star%central_force (p%r)              ! acceleration
      p%v = p%v + a*dt*0.5                      ! kick
      p%t = p%t + dt                            ! increment time
      self%n_step = self%n_step+1               ! count steps
    end do
    p => p%next                                 ! next particle
  end do
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Skeleton output procedure, which may also be used for -- for example --
!> measuring the total energy
!> =============================================================================
SUBROUTINE output (self)
  class(particles_t):: self
END SUBROUTINE output

END MODULE particles_mod
