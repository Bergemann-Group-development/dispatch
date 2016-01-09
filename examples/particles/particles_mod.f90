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
  implicit none
  private
  type, public, extends(task_t):: particle_task_t
    type(particle_t), pointer:: head => null()  ! particle list head
    integer:: n                                 ! number of particles
    integer(8):: n_step                         ! number of steps taken
    integer:: n_list                            ! switch to list when these remain
    real:: mass                                 ! total particle mass
    real:: radius                               ! particle cloud radius
    real:: dispersion                           ! velocity dispersion
    real:: e_rms                                ! rms energy variation
    real:: e_max_rms                            ! rms energy variation
    real:: e_kin                                ! kinetic energy
    real:: e_aver                               ! average energy
    real:: e2_aver                              ! average energy squared
    real:: dt_min                               ! min time step
    real:: dt_aver                              ! average time step
    real:: dt_harmonic                          ! harmonic mean time step
    real:: dt_max                               ! max time step
    real:: courant                              ! "Courant" number
  contains
    procedure:: init                            ! initialize
    procedure:: update_list                     ! update linked list
    procedure:: update_cache                    ! update via array cache
    procedure:: output                          ! output data
  end type
  type(random_t):: random
  type(star_t):: star
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
  integer:: n_list=20
  logical, save:: first_time=.true.
  namelist /particle_params/ n, radius, mass, dispersion, dtime, courant, &
    end_time, out_time, n_list
  !.............................................................................
  call trace_begin ('particles_mod%init', 3)
  call random%init
  if (first_time) then
    call star%init
    rewind (io%input)
    read (io%input, particle_params)
    if (io%master) write (*, particle_params)
  end if
  self%n = n
  self%n_step = 0
  self%n_list = n_list
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
  ! If new, allocate list of particles
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
  !-----------------------------------------------------------------------------
  ! If already initialized, reset to new, identical start
  !-----------------------------------------------------------------------------
  else
    p => self%head
    do while (associated(p))
      p%t = 0.0
      p%r = random%ran3()
      p%v = star%kepler_speed(p%r) * (1.+dispersion*random%ran3())
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
SUBROUTINE update_list (self)
  class(particle_task_t):: self
  real:: a(3)                                   ! accelration
  real:: dt, r2, c1, c2                         ! temporary data
  type(particle_t), pointer:: p                 ! particle pointer
  real:: e, e2, e_aver, e2_aver, e_rms, e_kin
  integer:: i, n_step
  !.............................................................................
  call trace_begin ('particles_mod%update', 3)
  p => self%head                                ! pointer to 1st particle
  c1 = star%m*star%g                            ! acceleration coefficient
  c2 = self%courant/c1                          ! time step coefficient
  i = 0
  do while (associated(p))                      ! for each list member
    i = i+1
    e_aver = 0.0
    e2_aver = 0.0
    n_step = 0
    do while (p%t < self%time)                  ! advance it to target time
      n_step = n_step + 1
      r2 = sum(p%r**2)                          ! radius square
      a = -p%r*(1./sqrt(r2))/r2*c1              ! 1./sqrt() is often fast
      !a = a + drag_force(p,dt)                  ! add drag force
      dt = maxval(abs(p%v))*r2*c2               ! using maxval for speed
      !dt = sqrt(sum(p%v**2)/sum(a**2))*self%courant
      dt = min(self%dtime, dt)                  ! don't overshoot too much
      p%v = p%v + a*dt*0.5                      ! kick
      p%r = p%r + p%v*dt                        ! drift
      r2 = sum(p%r**2)                          ! radius square
      a = -p%r*(1./sqrt(r2))/r2*c1              ! acceleration
      !a = a + drag_force(p,dt)                  ! add drag force
      p%v = p%v + a*dt*0.5                      ! kick
      p%t = p%t + dt                            ! increment time
      p%dt = dt                                 ! save for analysis
      e = 0.5*sum(p%v**2)                       ! kinetic energy
      e_kin = e_kin + e
      self%e_kin = self%e_kin + e
      e = e - c1/sqrt(r2)                       ! total energy
      e2 = e**2
      e_aver = e_aver + e
      e2_aver = e2_aver + e2
      self%e_aver = self%e_aver + e             ! accuracy measure
      self%e2_aver = self%e2_aver + e2          ! accuracy measure
      self%n_step = self%n_step+1               ! count steps
      self%dt_aver = self%dt_aver+dt            ! cost related measure
      self%dt_harmonic = self%dt_harmonic+1./dt ! cost related measure
      self%dt_min = min(self%dt_min,dt)         ! cost related measure
      self%dt_max = max(self%dt_max,dt)         ! cost related measure
    end do
    if (n_step > 0) then
      e_aver = e_aver/n_step
      e_kin  = e_kin /n_step
      e2_aver = e2_aver/n_step
      e_rms = sqrt(max(e2_aver-e_aver**2,0.0))/e_kin
      self%e_rms = self%e_rms + n_step*e_rms**2
      self%e_max_rms = max(self%e_max_rms,e_rms)
    end if
    p => p%next                                 ! next particle
  end do
  call trace_end
END SUBROUTINE update_list

!> =============================================================================
!> Optimized update, cashing into arrays for vectorization.  At the end of 
!> each update sweep any particles that are finished are moved to the end of
!> the array, and the number of still active particles is reduced.  When the
!> number of remaining particles drops to a small number (parameter n_list),
!> the remaining ones are updated in the linked list.
!> =============================================================================
SUBROUTINE update_arrays (self, n, m, r, v, dt, t)
  class(particle_task_t):: self
  integer:: n, m
  real, dimension(3,n):: r, v
  real, dimension(n):: dt
  real(8), dimension(n):: t
  real(8):: t1
  real:: r1(3), v1(3), dt1
  real:: c1, c2, t_max
  real, dimension(:), allocatable:: r2, r3, d, c
  real, dimension(:,:), allocatable:: v_gas
  integer:: i, j
  !...........................................................................
  allocate (r2(m), r3(m)) !, d(m), cs(m), v_gas(m,3))
  c1 = star%m*star%g                            ! acceleration coefficient
  c2 = self%courant/c1                          ! time step coefficient
  !---------------------------------------------------------------------------
  ! Kick. Optimize by using 1./sqrt(), max norm, and by keeping r2 available
  !---------------------------------------------------------------------------
  !call disk%structure (n, m, r, d, cs, v_gas)
  r2(1:m) = r(1,1:m)**2+r(2,1:m)**2+r(3,1:m)**2 
  r3(1:m) = (1./sqrt(r2(1:m)))**3
  !dt(1:m) = max(abs(v(1,1:m)),abs(v(2,1:m)),abs(v(3,1:m)))*r2(1:m)*c2
  dt(1:m) = sqrt(v(1,1:m)**2+v(2,1:m)**2+v(3,1:m)**2)*r2(1:m)*c2
  dt(1:m) = min(dt(1:m),self%dtime)
  t(1:m) = t(1:m) + dt(1:m)
  do i=1,3
    v(i,1:m) = v(i,1:m)-0.5*dt(1:m)*r(i,1:m)*r3(1:m)*c1
  end do
  !---------------------------------------------------------------------------
  ! Drift
  !---------------------------------------------------------------------------
  do i=1,3
    r(i,1:m) = r(i,1:m)+dt(1:m)*v(i,1:m)
  end do
  !---------------------------------------------------------------------------
  ! Kick
  !---------------------------------------------------------------------------
  r2(1:m) = r(1:m,1)**2+r(1:m,2)**2+r(1:m,3)**2 
  r3(1:m) = (1./sqrt(r2(1:m)))**3
  do i=1,3
    v(i,1:m) = v(i,1:m)-0.5*dt(1:m)*r(i,1:m)*r3(1:m)*c1
  end do
  self%n_step = self%n_step+m
  !print*,'b',t(1),dt(1),r2(1),0.5*sum(v(1,:)**2)-c1/sqrt(r2(1))
  !---------------------------------------------------------------------------
  ! Place any finished particles at the end.  This is NOT made faster by
  ! storing in the unused position m+1, resulting in a gap that needs to be
  ! removed when finished.
  !---------------------------------------------------------------------------
  t_max = maxval(t(1:m))
  if (t_max > self%time) then
    j = 0
    do while (j < m)
      j = j+1
      if (t(j) > self%time) then
        r1    =r(:,m); v1    =v(:,m); dt1  =dt(m); t1  =t(m)
        r(:,m)=r(:,j); v(:,m)=v(:,j); dt(m)=dt(j); t(m)=t(j) 
        r(:,j)=r1    ; v(:,j)=v1    ; dt(j)=dt1  ; t(j)=t1  ;
        m=m-1
      end if
    end do
  end if
  deallocate (r2, r3)
END SUBROUTINE update_arrays

!> =============================================================================
!> Optimized update, cashing into arrays to be able to use update_arrays(), which
!> after each update sweep moves particles that are finished to the end of the
!> arrays.  When the number of remaining particles drops to a small number
!> (parameter n_list), the remaining ones are updated directly in the linked list.
!> =============================================================================
SUBROUTINE update_cache (self)
  class(particle_task_t):: self
  type(particle_t), pointer:: p
  real, dimension(:,:), pointer:: r, v
  real, dimension(:), pointer:: dt
  real(8), dimension(:), pointer:: t
  real:: a(3)                                           ! acceleration
  real:: c1, c2                                         ! scaling constants
  integer:: i, m, n
  real(8):: t_min
  !.............................................................................
  call trace_begin ('particles_mod%update', 3)
  c1 = 1./(star%m*star%g)
  c2 = self%courant*c1
  !---------------------------------------------------------------------------
  ! Pull list date into arrays
  !---------------------------------------------------------------------------
  n = self%n
  allocate (r(3,n), v(3,n), dt(n), t(n))
  i = 0
  p => self%head
  do while (associated(p))
    i = i+1
    r(:,i) = p%r
    v(:,i) = p%v
    t(i) = p%t
    p => p%next
  end do
  !---------------------------------------------------------------------------
  ! Update, while moving the finished ones to the end, skipping them next time
  !---------------------------------------------------------------------------
  n = self%n
  m = n
  t_min = 0.0
  do while (t_min < self%time .and. n > self%n_list)
    call update_arrays (self, n, m, r, v, dt, t)
    t_min = minval(t(1:m))
  end do
  !---------------------------------------------------------------------------
  ! Push back to list
  !---------------------------------------------------------------------------
  i = 0
  p => self%head
  do while (associated(p))
    i = i+1
    p%r = r(:,i)
    p%v = v(:,i)
    p%dt = dt(i)
    p%t = t(i)
    p => p%next
  end do
  call update_list (self)
  deallocate (r, v, dt, t)
  call trace_end
END SUBROUTINE update_cache

!> =============================================================================
!> Skeleton output procedure, which may also be used for -- for example --
!> measuring the total energy
!> =============================================================================
SUBROUTINE output (self)
  class(particle_task_t):: self
END SUBROUTINE output

END MODULE particles_mod
