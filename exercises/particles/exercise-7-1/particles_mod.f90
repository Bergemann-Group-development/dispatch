!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Demonstrate particle update method, and measure update speed
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE particles_mod
  USE io_mod
  USE trace_mod
  USE omp_timer_mod
  USE random_mod
  USE task_mod
  USE star_mod
  USE particle_mod
  USE standard_disk_mod
  USE scaling_mod
  implicit none
  private
  type, public, extends(task_t):: particle_task_t
    type(particle_t), pointer:: head => null()  ! particle list head
    integer:: n                                 ! number of particles
    integer(8):: n_step                         ! number of steps taken
    real:: mass                                 ! total particle mass
    real:: radius                               ! particle cloud radius
    real:: dispersion                           ! velocity dispersion
    real:: courant                              ! "Courant" number
    real:: e_rms                                ! rms energy variation
    real:: e_max_rms                            ! rms energy variation
    real:: e_kin                                ! kinetic energy
    real:: e_aver                               ! average energy
    real:: e2_aver                              ! average energy squared
    real:: dt_min                               ! min time step
    real:: dt_aver                              ! average time step
    real:: dt_harmonic                          ! harmonic mean time step
    real:: dt_max                               ! max time step
  contains
    procedure:: init                            ! initialize
    procedure:: update_kdk                      ! update linked list
    procedure:: update_rk3                      ! update linked list
    procedure:: output                          ! output data
  end type
  type(random_t):: random
  type(standard_disk_t):: disk
CONTAINS

!> =============================================================================
!> Initialize a number of particles consistent with a typical particle task, with
!> a given total mass, in Kepler orbits with a specified extra velocity dispersion.
!> The "courant" number parameter is based on comparing the acceleration to the
!> velocity, and controls the fraction of an orbit taken.  For a circular orbit,
!> a ratio of unity corresponds to taking 2\pi steps per orbit.
!> =============================================================================
SUBROUTINE init (self)
  class(particle_task_t):: self
  type(particle_t), pointer:: p
  integer:: n=10000, i
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
  call random%init
  if (first_time) then
    call scaling%init
    call star%init
    call disk%init
    rewind (io%input)
    read (io%input, particle_params)
    if (io%master) write (*, particle_params)
  end if
  self%n = n
  self%n_step = 0
  self%mass = mass
  self%radius = radius
  self%dispersion = dispersion
  self%time = 0.0
  self%dtime = dtime
  self%end_time = end_time
  self%out_time = out_time
  self%courant = courant
  self%dt_min = 1e10
  self%dt_aver = 0.0
  self%dt_harmonic = 0.0
  self%dt_max = 0.0
  self%e_rms = 0.0
  self%e_aver = 0.0
  self%e2_aver = 0.0
  !-----------------------------------------------------------------------------
  ! Alocate a list of particles
  !-----------------------------------------------------------------------------
  if (first_time) then
    do i=1,n
      allocate (p)
      p%next => self%head
      p%r = random%ran3()
      p%v = star%kepler_speed(p%r) * (1.+dispersion*random%ran3())
      p%m = self%mass/n
      self%head => p
    end do
  else
    p => self%head
    do i=1,n
      p%r = random%ran3()
      p%v = star%kepler_speed(p%r) * (1.+dispersion*random%ran3())
      p%m = self%mass/n
      p%t = 0.0
      p => p%next
    end do
  end if
  first_time = .false.
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Update particles using the Kick-Drift-Kick algorithm (Springel, 19xx).  The
!> particles are advanced individually, for a time interval of about self%dtime,
!> which is intended to correspond to a typical patch update time. There is no
!> need to require particles to land exactly on the target time.
!> =============================================================================
SUBROUTINE update_kdk (self)
  class(particle_task_t):: self
  real:: a(3), atot(3)                          ! accelration
  real:: dt, r2, c1, c2                         ! temporary data
  type(particle_t), pointer:: p                 ! particle pointer
  real:: e, e2, e_aver, e2_aver, e_rms, e_kin
  real:: v_gas(3), c, b, am
  integer:: i, n_step
  !.............................................................................
  call trace_begin ('particles_mod%update_kdk', 3)
  c1 = star%m*star%g                            ! acceleration coefficient
  c2 = self%courant/c1                          ! time step coefficient
  i = 0                                         ! particle counter
  p => self%head                                ! pointer to 1st particle
  do while (associated(p))                      ! for each list member
    i = i+1
    e_aver = 0.0                                ! average total energy
    e2_aver = 0.0                               ! squared
    e_kin = 0.0                                 ! average kinetic energy
    n_step = 0                                  ! steps along path
    c = 0.0                                     ! no drag
    v_gas = 0.0                                 ! for testing
    do while (p%t < self%time)                  ! advance it to target time
      !------------------------------------------------------------------------
      ! First sub-step
      !------------------------------------------------------------------------
      if (n_step==0) then                       ! calculations only needed 1st time
        r2 = sum(p%r**2)                        ! radius square
        am = c1/r2                              ! magnitude of acceleration
        a = -p%r*(1./sqrt(r2))*am               ! 1./sqrt() is often fast in hardware
      end if
      n_step = n_step + 1
      dt = sqrt(sum(p%v**2))/am*self%courant    ! time step limited by curvature
      dt = min(self%dtime, dt)                  ! don't overshoot too much
     !call disk%drag (p%r+2.*dt*p%v, p%v+2.*dt*a, v_gas, c)       ! the v_gas should be at a position dt ahead of the current one
     !dt = min(dt,0.01/max(c,1e-35))
      b = 0.5                                   ! control asymptotic behavior
      atot = (a + (v_gas-p%v)*c)/(1.0+b*c*dt)   ! total acceleration
      p%v = p%v + 0.5*dt*atot                   ! update velocity
      p%r = p%r + dt*p%v                        ! drift
      p%t = p%t + dt                            ! increment time
      p%dt = dt                                 ! save for analysis
      !------------------------------------------------------------------------
      ! Second sub-step
      !------------------------------------------------------------------------
      r2 = sum(p%r**2)                          ! radius square
      am = c1/r2                                ! magnitude of acceleration
      a = -p%r*(1./sqrt(r2))*am                 ! acceleration
     !call disk%drag (p%r, p%v, v_gas, c)       ! the position for this v_gas should be at the current r
      atot = (a + (v_gas-p%v)*c)/(1.0+b*c*dt)   ! total acceleration
      p%v = p%v + 0.5*dt*atot                   ! update velocity
      !------------------------------------------------------------------------
      ! Accumulate statistics
      !------------------------------------------------------------------------
      e = 0.5*sum(p%v**2)                       ! kinetic energy
      e_kin = e_kin + e                         ! accumulate over one path
      self%e_kin = self%e_kin + e               ! accumulate over all paths
      e = e - c1/sqrt(r2)                       ! total energy
      e2 = e**2                                 ! squared
      e_aver = e_aver + e                       ! average total energy
      e2_aver = e2_aver + e2                    ! squared, for RMS
      self%e_aver = self%e_aver + e             ! accuracy measure
      self%e2_aver = self%e2_aver + e2          ! accuracy measure
      self%n_step = self%n_step+1               ! count steps
      self%dt_aver = self%dt_aver+dt            ! cost related measure
      self%dt_harmonic = self%dt_harmonic+1./dt ! cost related measure
      self%dt_min = min(self%dt_min,dt)         ! cost related measure
      self%dt_max = max(self%dt_max,dt)         ! cost related measure
    end do
    if (n_step > 0) then
      e_aver = e_aver/n_step                    ! average energy for this path
      e_kin  = e_kin /n_step                    ! average kinetic energy for this path
      e2_aver = e2_aver/n_step
      e_rms = sqrt(max(e2_aver-e_aver**2,0.0))/e_kin
      self%e_rms = self%e_rms + n_step*e_rms**2 ! Root-mean square
      self%e_max_rms = max(self%e_max_rms,e_rms)
    end if
    p => p%next                                 ! next particle
  end do
  call trace_end
END SUBROUTINE update_kdk

!> =============================================================================
!> =============================================================================
SUBROUTINE update_rk3 (self)
  class(particle_task_t):: self
  real:: a(3), atot(3)                          ! acceleration
  real, dimension(3):: a1, a2, a3               ! accelerations
  real, dimension(3):: v1, v2, v3               ! velocities
  real, dimension(3):: r1, r2, r3               ! positions
  real:: dt, rs, c1, c2                         ! temporary data
  type(particle_t), pointer:: p                 ! particle pointer
  real:: e, e2, e_aver, e2_aver, e_rms, e_kin
  real:: v_gas(3), c, b, am
  integer:: i, n_step
  !.............................................................................
  call trace_begin ('particles_mod%update_rk3', 3)
  c1 = star%m*star%g                            ! acceleration coefficient
  c2 = self%courant/c1                          ! time step coefficient
  i = 0                                         ! particle counter
  p => self%head                                ! pointer to 1st particle
  do while (associated(p))                      ! for each list member
    i = i+1
    e_aver = 0.0                                ! average total energy
    e2_aver = 0.0                               ! squared
    e_kin = 0.0                                 ! average kinetic energy
    n_step = 0                                  ! steps along path
    c = 0.0                                     ! no drag
    v_gas = 0.0                                 ! for testing
    do while (p%t < self%time)                  ! advance it to target time
      n_step = n_step + 1
      !------------------------------------------------------------------------
      ! 3rd order Runge-Kutta step
      !------------------------------------------------------------------------
      r1 = p%r                                  ! initial position
      v1 = p%v                                  ! initial speed
      rs = sum(r1**2)                           ! radius square
      am = c1/rs                                ! magnitude of acceleration
      a1 = -r1*(1./sqrt(rs))*am                 ! 1./sqrt() is often fast in hardware
      dt = sqrt(sum(p%v**2))/am*self%courant    ! time step limited by curvature
      dt = min(self%dtime, dt)                  ! don't overshoot too much
      v2 = v1 + 0.5*dt*a1                       ! update velocity
      r2 = r1 + 0.5*dt*v1                       ! drift
      rs = sum(r2**2)                           ! radius square
      am = c1/rs                                ! magnitude of acceleration
      a2 = -r2*(1./sqrt(rs))*am                 ! 1./sqrt() is often fast in hardware
      v3 = v1 + dt*(2.0*a2-a1)                  ! update velocity
      r3 = r1 + dt*(2.0*v2-v1)                  ! update position
      rs = sum(r3**2)                           ! radius square
      am = c1/rs                                ! magnitude of acceleration
      a3 = -r3*(1./sqrt(rs))*am                 ! 1./sqrt() is often fast in hardware
      p%r = r1 + (dt/6.)*(v1+4.0*v2+v3)         ! final position
      p%v = v1 + (dt/6.)*(a1+4.0*a2+a3)         ! final velocity
      p%t = p%t + dt                            ! increment time
      p%dt = dt                                 ! save for analysis
      !------------------------------------------------------------------------
      ! Accumulate statistics
      !------------------------------------------------------------------------
      e = 0.5*sum(p%v**2)                       ! kinetic energy
      e_kin = e_kin + e                         ! accumulate over one path
      self%e_kin = self%e_kin + e               ! accumulate over all paths
      e = e - c1/sqrt(rs)                       ! total energy
      e2 = e**2                                 ! squared
      e_aver = e_aver + e                       ! average total energy
      e2_aver = e2_aver + e2                    ! squared, for RMS
      self%e_aver = self%e_aver + e             ! accuracy measure
      self%e2_aver = self%e2_aver + e2          ! accuracy measure
      self%n_step = self%n_step+1               ! count steps
      self%dt_aver = self%dt_aver+dt            ! cost related measure
      self%dt_harmonic = self%dt_harmonic+1./dt ! cost related measure
      self%dt_min = min(self%dt_min,dt)         ! cost related measure
      self%dt_max = max(self%dt_max,dt)         ! cost related measure
    end do
    if (n_step > 0) then
      e_aver = e_aver/n_step                    ! average energy for this path
      e_kin  = e_kin /n_step                    ! average kinetic energy for this path
      e2_aver = e2_aver/n_step
      e_rms = sqrt(max(e2_aver-e_aver**2,0.0))/e_kin
      self%e_rms = self%e_rms + n_step*e_rms**2 ! Root-mean square
      self%e_max_rms = max(self%e_max_rms,e_rms)
    end if
    p => p%next                                 ! next particle
  end do
  call trace_end
END SUBROUTINE update_rk3
!> =============================================================================
!> Skeleton output procedure, which may also be used for -- for example --
!> measuring the total energy
!> =============================================================================
SUBROUTINE output (self)
  class(particle_task_t):: self
END SUBROUTINE output

END MODULE particles_mod
