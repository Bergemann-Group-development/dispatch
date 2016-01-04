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
    type(particle_t), pointer:: head => null()    ! particle list head
    integer:: n
    integer:: n_step
    integer:: n_list
    real:: mass
    real:: radius
    real:: dispersion
    real:: dt_min
    real:: dt_aver
    real:: dt_harmonic
    real:: dt_max
    real:: courant
  contains
    procedure:: init
    procedure:: update
    procedure:: update_cache
    procedure:: update_bags
    procedure:: output
  end type
  type, public:: particle_bag_t
    integer:: n
    real:: dt_max
    real, dimension(:,:), pointer:: r, v
    real, dimension(:), pointer:: dt
    real(8), dimension(:), pointer:: t
  end type
  type(random_t):: random
  type(star_t):: star
CONTAINS

!> =============================================================================
!> Initialize a number of particles consistent with a particle task, with a given
!> total mass, in Kepler orbits with a specified extra velocity dispersion.
!> =============================================================================
SUBROUTINE init (self)
  class(particle_task_t):: self
  integer:: n=10000, i
  type(particle_t), pointer:: p
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
    !---------------------------------------------------------------------------
    ! If already initialized, reset to new start
    !---------------------------------------------------------------------------
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
!> particles are advanced individually, but for a maximum time of self%dtime,
!> which is intended to correspond to a typical patch update time.
!> =============================================================================
SUBROUTINE update (self)
  class(particle_task_t):: self
  real:: a(3)                                   ! force
  real:: dt, r2, c1, c2
  type(particle_t), pointer:: p
  !.............................................................................
  call trace_begin ('particles_mod%update', 3)
  p => self%head
  self%time = self%time + self%dtime            ! target time
  c1 = 1./(star%m*star%g)
  c2 = self%courant*c1
  do while (associated(p))
    do while (p%t < self%time)
      r2 = sum(p%r**2)                          ! radius square
      a = -p%r*(1./sqrt(r2))/r2*c1              ! 1./sqrt() is often fast
      dt = maxval(abs(p%v))*r2*c2               ! using maxval for speed
      dt = min(self%dtime, dt)                  ! don't overshoot too much
      p%v = p%v + a*dt*0.5                      ! kick
      p%r = p%r + p%v*dt                        ! drift
      r2 = sum(p%r**2)                          ! radius square
      a = -p%r*(1./sqrt(r2))/r2*c1              ! acceleration
      p%v = p%v + a*dt*0.5                      ! kick
      p%t = p%t + dt                            ! increment time
      p%dt = dt                                 ! save for "bagging"
      self%n_step = self%n_step+1               ! count steps
      self%dt_aver = self%dt_aver+dt            ! cost related measure
      self%dt_harmonic = self%dt_harmonic+1./dt ! cost related measure
      self%dt_min = min(self%dt_min,dt)         ! cost related measure
      self%dt_max = max(self%dt_max,dt)         ! cost related measure
    end do
    p => p%next                                 ! next particle
  end do
  call trace_end
END SUBROUTINE update

!> =============================================================================
!> Optimized update, using arrays for vectorization
!> =============================================================================
SUBROUTINE update_arrays (self, n, m, r, v, dt, t)
  class(particle_task_t):: self
  integer:: n, m
  real, dimension(n,3):: r, v
  real, dimension(n):: dt
  real(8), dimension(n):: t
  real:: r1(3), v1(3), dt1
  real(8):: t1
  real, dimension(:), allocatable:: r2, r3
  real:: c1, c2, t_max
  integer:: i
  !...........................................................................
  allocate (r2(n), r3(n))
  c1 = 1./(star%m*star%g)
  c2 = self%courant*c1
  !---------------------------------------------------------------------------
  ! Kick. Optimize by using 1./sqrt(), max norm, and by keeping r2 available
  !---------------------------------------------------------------------------
  r2(:) = r(:,1)**2+r(:,2)**2+r(:,3)**2 
  r3(:) = (1./sqrt(r2(:)))**3
  dt(:) = max(abs(v(:,1)),abs(v(:,2)),abs(v(:,3)))/r2(:)*c2
  dt(:) = min(dt(:),self%dtime)
  t(:) = t(:) + dt(:)
  do i=1,3
    v(:,i) = v(:,i)-0.5*dt(:)*r(:,i)*r3(:)*c1
  end do
  !---------------------------------------------------------------------------
  ! Drift
  !---------------------------------------------------------------------------
  do i=1,3
    r(:,i) = r(:,i)+dt(:)*v(:,i)
  end do
  !---------------------------------------------------------------------------
  ! Kick
  !---------------------------------------------------------------------------
  r2(:) = r(:,1)**2+r(:,2)**2+r(:,3)**2 
  r3(:) = (1./sqrt(r2(:)))**3
  do i=1,3
    v(:,i) = v(:,i)-0.5*dt(:)*r(:,i)*r3(:)*c1
  end do
  self%n_step = self%n_step+n
  !---------------------------------------------------------------------------
  ! Place any finished particles at the end
  !---------------------------------------------------------------------------
  t_max = maxval(t)
  m = n
  if (t_max > self%end_time) then
    i = 0
    do while (i < m)
      i = i+1
      if (t(i) > self%end_time) then
        r1    =r(m,:); v1    =v(m,:); dt1  =dt(m); t1  =t(m)
        r(m,:)=r(i,:); v(m,:)=v(i,:); dt(m)=dt(i); t(m)=t(i) 
        r(i,:)=r1    ; v(i,:)=v1    ; dt(i)=dt1  ; t(i)=t1  ;
        m=m-1
      end if
    end do
  end if
  deallocate (r2, r3)
END SUBROUTINE update_arrays

!> =============================================================================
!> Check speed of updates when caching the linked list values in an array.
!> =============================================================================
SUBROUTINE update_cache (self)
  class(particle_task_t):: self
  real:: a(3)                                   ! force
  type(particle_t), pointer:: p
  real, dimension(:,:), pointer:: r, v
  real, dimension(:), pointer:: r2, r3, dt
  real(8), dimension(:), pointer:: t
  real:: c1, c2, dt0
  integer:: i, j, m, n
  real(8):: wc(5)
  !.............................................................................
  call trace_begin ('particles_mod%update', 3)
  n = self%n
  allocate (r(n,3), v(n,3), r2(n), r3(n), dt(n), t(n))
  !-----------------------------------------------------------------------------
  ! Cache, for speed measuring purposes only, the positions and velocities
  !-----------------------------------------------------------------------------
  wc(1) = wallclock()
  call tic
  dt0 = 0.0003
  self%time = 0.0
  self%n_step = 0
  c1 = 1./(star%m*star%g)
  c2 = self%courant*c1
  !---------------------------------------------------------------------------
  ! Pull list into array
  !---------------------------------------------------------------------------
  i = 0
  p => self%head
  do while (associated(p))
    i = i+1
    r(i,:) = p%r
    v(i,:) = p%v
    t(i) = p%t
    p => p%next
  end do
  wc(2) = wallclock()
  !---------------------------------------------------------------------------
  ! Update, while moving the finished ones to the end, and skipping them
  !---------------------------------------------------------------------------
  n = self%n
  do while (self%time < self%end_time .and. n > self%n_list)
    call update_arrays (self, n, m, r, v, dt, t)
    self%time = minval(t(1:m))
    n = m
    !print *, n, m, self%time, minval(dt(1:m)), maxval(dt(1:m))
  end do
  wc(3) = wallclock()
  !---------------------------------------------------------------------------
  ! Push back to list
  !---------------------------------------------------------------------------
  i = 0
  p => self%head
  do while (associated(p))
    i = i+1
    p%r = r(i,:)
    p%v = v(i,:)
    p%dt = dt(i)
    p%t = t(i)
    p => p%next
  end do
  wc(4) = wallclock()
  self%time = self%end_time - self%dtime
  call update (self)
  wc(5) = wallclock()
  call toc ('update_cache', self%n_step)
  deallocate (r, v, r2, r3, dt, t)
  wc(1:4) = wc(2:5)-wc(1:4)
  print*,'   gather     array   scatter      list'
  print'(3p,4f10.3)', wc(1:4)
  call trace_end
END SUBROUTINE update_cache

!> =============================================================================
!> Distrubute the particles into "bags" with similar number of expected steps.
!> Each bag construction and update requires running through the particle list
!> =============================================================================
SUBROUTINE update_bags (self)
  class(particle_task_t):: self
  type(particle_bag_t), dimension(:), pointer:: bag
  type(particle_t), pointer:: p
  real, dimension(:,:), pointer:: r, v
  real, dimension(:), pointer:: r2, r3, dt
  real:: c, dt_bag, t_min, t_max
  integer:: i, j, n, n_bag, n_step, i_step
  !---------------------------------------------------------------------------
  ! Count how many bags we need, and allocate them (one extra for finished)
  !---------------------------------------------------------------------------
  print*,'============= update_bags =============='
  dt_bag = self%dtime
  n_bag = 0
  do while (dt_bag > self%dt_min)
    n_bag = n_bag+1
    dt_bag = 0.5*dt_bag
  end do
  print *, 'n_bag:', n_bag
  allocate (bag(0:n_bag))
  !---------------------------------------------------------------------------
  ! Count how many particles go into each bag, and allocate arrays
  !---------------------------------------------------------------------------
  print*,'count bag sizes'
  p => self%head
  bag(:)%n = 0
  n = 0
  p => self%head
  do while (associated(p))
    n = n+1
    if (p%t > self%end_time .or. p%t > self%time + self%dtime) then
      i = 0
    else
      i = 1.0 + log((self%time + self%dtime - p%t)/p%dt)/log(2.)
    end if
    !print*, n, i
    i = min(n_bag,max(1,i))
    bag(i)%n = bag(i)%n + 1
    p => p%next
  end do
  do i=0,n_bag
    bag(i)%dt_max = self%dt_max*0.5**(i-1)
    print *, i, bag(i)%n, dt_bag
    n = bag(i)%n
    if (n > 0) &
      allocate (bag(i)%r(n,3), bag(i)%v(n,3), bag(i)%dt(n), bag(i)%t(n))
  end do
  !---------------------------------------------------------------------------
  ! Fill the bags
  !---------------------------------------------------------------------------
  print*,'fill bags'
  bag(:)%n = 0
  n = 0
  p => self%head
  do while (associated(p))
    n = n+1
    if (p%t > self%end_time .or. p%t > self%time + self%dtime) then
      i = 0
    else
      i = 1.0 + log((self%time + self%dtime - p%t)/p%dt)/log(2.)
    end if
    i = min(n_bag,max(1,i))
    bag(i)%n = bag(i)%n+1
    !print *, n, i
    bag(i)%r(bag(i)%n,:) = p%r
    bag(i)%v(bag(i)%n,:) = p%v
    bag(i)%dt(bag(i)%n) = p%dt
    bag(i)%t(bag(i)%n) = p%t
    p => p%next
  end do
  !---------------------------------------------------------------------------
  ! Update each bag, empty it back onto the linked list, and waste it.  The
  ! timesteps in each bag are in the range [0.5,1.0]*bag(i)%dt_max
  !---------------------------------------------------------------------------
  print*,'update bags'
  p => self%head
  do i=0,n_bag
    if (i==0) then
      n_step=0
    else
      n_step = 2.0*self%dtime/bag(i)%dt_max + 0.1
    end if
    if (bag(i)%n > 0) then
      do i_step=1,n_step
        call update_arrays (self, bag(i)%n, n, bag(i)%r, bag(i)%v, bag(i)%dt, bag(i)%t)
      end do
      t_min = 1e10
      t_max = 0.0
      do j=1,bag(i)%n
        p%r  = bag(i)%r(j,:)
        p%v  = bag(i)%v(j,:)
        p%dt = bag(i)%dt(j)
        p%t  = bag(i)%t(j)
        t_min = min(t_min,p%t)
        t_max = max(t_max,p%t)
      end do
      print *, i, n_step, t_min, t_max
      deallocate (bag(i)%r, bag(i)%v, bag(i)%dt, bag(i)%t)
    end if
  end do
  deallocate (bag)
  self%time = self%time + self%dtime
  print *, self%time, t_min, t_max
END SUBROUTINE update_bags

!> =============================================================================
!> Skeleton output procedure, which may also be used for -- for example --
!> measuring the total energy
!> =============================================================================
SUBROUTINE output (self)
  class(particle_task_t):: self
END SUBROUTINE output

END MODULE particles_mod
