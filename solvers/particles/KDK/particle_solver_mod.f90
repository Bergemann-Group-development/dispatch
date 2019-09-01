!===============================================================================
!> Kick-Drift-Kick (KDK) N-body solver.  For each particle, the steps are:
!>
!> 1. Make a prediction of the position half a step forward in time, using the
!>    particles history data in the particles_t data type.
!> 2. Decide on the time step size, using the forces and velocities estimated
!>    at that point
!> 3. Compute forces acting on the particle at the initial time
!> 4. Take a 1st kick (K) step, updating the velocity by half the time step
!> 5. Take a drift (D) step, updating the position by a full time step
!> 6. Compute the forces at the updated position, and use these to ..
!> 7. Take a 2nd kick (K) stop, updating the velocity by half the time step
!===============================================================================
MODULE particle_solver_mod
  USE io_unit_mod
  USE trace_mod
  USE particle_mod
  USE particles_mod
  USE dll_mod
  implicit none
  private
  type, public, extends(particles_t):: particle_solver_t
  contains
    procedure:: init
    procedure:: update
    procedure:: courant_time
  end type
  integer, save:: verbose=0
  real(8), save:: courant=0.5
  type(particle_solver_t), target, public:: particle_solver
CONTAINS

!===============================================================================
!> Initialize a particle list
!===============================================================================
SUBROUTINE init (self, name)
  class(particle_solver_t):: self
  character(len=*), optional:: name
  logical, save:: first_time=.true.
  integer:: iostat
  namelist /kdk_params/ verbose, courant
  !-----------------------------------------------------------------------------
  call self%particles_t%init
  if (first_time) then
    !$omp critical (particle_solver_input_cr)
    if (first_time) then
      first_time = .false.
      rewind (io_unit%input)
      read (io_unit%input, kdk_params, iostat=iostat)
      write (io_unit%output, kdk_params)
      call particles%init ('solv')
      particle_solver%verbose = verbose
      particles%verbose = verbose
    end if
    !$omp end critical (particle_solver_input_cr)
  end if
END SUBROUTINE init

!===============================================================================
!> KDK update, using the self%force() procedure provided by particles_t.
!===============================================================================
SUBROUTINE update (self, p)
  class(particle_solver_t):: self
  class(particle_t), pointer:: p
  !.............................................................................
  type(particle_t):: pp
  real(8):: dt, a(3), v(3), r(3), dv(3)
  real(8), pointer:: tnew
  integer:: prv, now, new, nt, i
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%update')
  !$omp atomic read
  now = p%it
  new = mod   (now  ,p%nt) + 1
  prv = modulo(now-2,p%nt) + 1
  dt = self%courant_time (p)                          ! estimate dt
  !$omp atomic write
  p%dtime = dt
  do i=1,3
    !$omp atomic read
    r(i) = p%r(i,now)                                 ! pick up
    !$omp atomic read
    v(i) = p%v(i,now)                                 ! pick up
  end do
  dv = 0.5_8*dt*self%force(p,now)/p%mass              ! kick 1
if (particle_solver%verbose > 1) &
print *,'K: dv =', p%id, dv
  v = v + dv           ! kick 1
if (particle_solver%verbose > 1) &
print *,'K: v1 =', p%id, v
if (particle_solver%verbose > 1) &
print *,'D: r  =', p%id, r
  r = r + v*dt                                        ! new pos
  do i=1,3
    !$omp atomic write
    p%r(i,new) = r(i)                                 ! store
  end do
  dv = 0.5_8*dt*self%force(p,new)/p%mass              ! kick 2
  v = v + dv                                          ! kick 2
if (particle_solver%verbose > 1) &
print *,'K: v2 =', p%id, v
  do i=1,3
    !$omp atomic write
    p%v(i,new) = v(i)                                 ! store
  end do
  !-----------------------------------------------------------------------------
  ! Rotate the time slice slots for the particle
  !-----------------------------------------------------------------------------
  !$omp atomic update
  p%time   = p%time+dt
  tnew => p%t(new)
  !$omp atomic write
  tnew = p%time
  !$omp atomic write
  p%it     = new
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Courant conditions:  Estimate, as inexpensively as possible, an approximately
!> reflexive time step. If an estimate of the acceleration results in it not
!> being important for the timestep size, then choose a fraction of the cell
!> travel time (whatever that means in this context)
!===============================================================================
FUNCTION courant_time (self, p) RESULT (dt)
  class(particle_solver_t):: self
  class(particle_t), pointer:: p
  real(8):: dt
  !.............................................................................
  real(8):: a(3), ac, vc, dt1, dt2
  integer:: prv, now, new
  !-----------------------------------------------------------------------------
  call p%predict()
  now = p%it
  new = mod   (now  ,p%nt) + 1
  prv = modulo(now-2,p%nt) + 1
  associate (v=>p%v(:,new))
  a = self%force(p,new)/p%mass
  ac = 1.0/sqrt(sum(a**2)+tiny(1.0_8))
  vc = 1.0/sqrt(sum(v**2)+tiny(1.0_8))
  dt1 = p%ds*vc*courant
  dt2 = p%ds*(ac/vc)*courant
  dt = min(dt1,dt2)
if (particle_solver%verbose > 0) &
print *, 'courant_time: dt1, dt2, ac, vc =', dt1, dt2, ac, vc
  end associate
END FUNCTION courant_time

END MODULE particle_solver_mod
