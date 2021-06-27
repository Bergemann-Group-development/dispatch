!===============================================================================
!> Particle solver for handling dust particles represented by array data,
!> allowing vectorization of the KDK part and using only 36 bytes per particle.
!> Updating the particles is done by looping over cells, and looping over bins
!> and particles in bins within each cells, with these steps:
!>
!> 1) compute cell gas properties needed for the update
!> 2) update particle positions with a kick-drift-kick step
!> 3) trade particles with neighboring cells
!> 4) import particles from neighboring patches
!> 5) remove holes in the arrays caused by the trading
!> 6) sumup the particle weights to get dust densities
!> 7) add or remove particles and renormalize weights
!>
!> Step 3 is the most expensive one, using about 35% of the time
!===============================================================================
MODULE particle_solver_mod
  USE omp_lib
  USE io_mod
  USE trace_mod
  USE timer_mod
  USE random_mod
  USE task_mod
  USE gpatch_mod
  USE scaling_mod
  USE units_mod
  USE list_mod
  USE link_mod
  USE omp_lock_mod
  USE omp_timer_mod
  USE particles_mod
  USE mesh_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(task_t):: particle_solver_t
    integer                    :: n=32, n_bin=5, mp=60, me, np, lb(3), ub(3)
    real, pointer, contiguous  :: s_cgs(:)            => null()
    real, pointer, contiguous  :: m_cgs(:)            => null()
    real, pointer, contiguous  :: a(:,:,:,:)          => null()
    real, pointer, contiguous  :: v(:,:,:,:)          => null()
    real, pointer, contiguous  :: d(:,:,:)            => null()
    real, pointer, contiguous  :: dd(:,:,:,:)         => null()
    type(particles_t), pointer :: particles(:,:,:,:)  => null()
    type(particles_t), pointer :: previous(:,:,:,:)   => null()
    type(particles_t), pointer :: exports(:)          => null()
    class(gpatch_t), pointer   :: patch               => null()
    class(list_t), pointer     :: task_list           => null()
    class(nbor_t), pointer     :: nbor                => null()
  contains
    procedure:: connect
    procedure:: init
    procedure:: init_task_list
    procedure:: initial_condition
    procedure:: update
    procedure, private:: setup
    procedure, private:: copy
    procedure:: pack
    procedure, private:: import
    procedure, private:: kdk
    procedure, private:: count
    procedure, private:: trade
    procedure, private:: sumup
    procedure, private:: renorm
    procedure:: output
  end type
  logical:: on = .true.
  logical:: particle_io = .false.
  logical:: count_all = .false.
  logical:: backreaction = .false.
  integer:: verbose = 0
  integer:: n_bin = 1
  integer:: particles_per_bin = 30
  integer:: output_every = 1000
  integer:: seed = 11                         ! random seed
  real:: dust_to_gas = 0.01                   ! total dust-to-gas ratio
  real:: min_fraction = 0.7                   ! particles per bin (max decrease)
  real:: max_fraction = 1.4                   ! particles per bin (max increase)
  real:: min_size_cgs = 3e-4                  ! 3 micron
  real:: max_size_cgs = 3000e-4               ! 3000 micron = 3 mm
  real:: solid_density_cgs = 2.0              ! g/cm**3
  real:: density_yr_cgs = 3e-12               ! volume density giving 1 yr stop
  real(8):: out_time = 0.1                    ! output cadence
  real(8):: start_time = 0.0                  ! output start
CONTAINS

!===============================================================================
!> Initialize patch link
!===============================================================================
SUBROUTINE connect (self, patch)
  class(particle_solver_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init')
  self%patch => patch
  call trace%end ()
END SUBROUTINE connect

!===============================================================================
!> Initialize dust parameters
!===============================================================================
SUBROUTINE init (self)
  class(particle_solver_t):: self
  integer:: iostat
  namelist /particle_solver_params/ on, verbose, particles_per_bin, n_bin, seed, &
    dust_to_gas, min_size_cgs, max_size_cgs, min_fraction, max_fraction, &
    solid_density_cgs, density_yr_cgs, particle_io, out_time, start_time, &
    output_every, backreaction
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init')
  !.............................................................................
  call io%lock%set ('particle_solver_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, particle_solver_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, particle_solver_params)
    end if
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%namelist, particle_solver_params)
    first_time = .false.
  end if
  call io%lock%unset ('particle_solver_t%init')
  !.............................................................................
  self%n_bin = n_bin
  self%mp = 2*max_fraction*particles_per_bin
  self%me = particles_per_bin*self%patch%n(1)*self%patch%n(2)
  call self%random%init (seed)
  call self%aux%init
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize dust bins
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particle_solver_t):: self
  integer:: mp, np, nmax, iostat
  integer:: i1, i2, i3, ib, j, m(3)
  real:: w
  type(particles_t), pointer:: p, o
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin ('particle_solver_t%initial_condition')
  m = self%patch%gn
  allocate (self%exports(n_bin))
  allocate (self%particles(m(1),m(2),m(3),n_bin))
  allocate (self%previous (m(1),m(2),m(3),n_bin))
  allocate (self%v(      3,m(1),m(2),m(3)))
  allocate (self%a(      3,m(1),m(2),m(3)))
  allocate (self%dd(       m(1),m(2),m(3),n_bin))
  allocate (self%s_cgs(n_bin), self%m_cgs(n_bin))
  call io%bits_mem (storage_size(self%exports), n_bin, 'mem')
  call io%bits_mem (storage_size(self%particles), &
    product(shape(self%particles)), 'mem')
  call io%bits_mem (storage_size(self%previous), &
    product(shape(self%previous)), 'mem')
  call io%bits_mem (storage_size(self%v ), product(shape(self%v )), 'v')
  call io%bits_mem (storage_size(self%a ), product(shape(self%a )), 'a')
  call io%bits_mem (storage_size(self%dd), product(shape(self%dd)), 'dd')
  self%d => self%patch%mem(:,:,:,self%patch%idx%d,self%patch%it,1)
  do ib=1,n_bin
    self%s_cgs(ib) = exp(log(min_size_cgs) &
                     + (ib-0.5)*log(max_size_cgs/min_size_cgs)/n_bin)
    self%m_cgs(ib) = solid_density_cgs*self%s_cgs(ib)**3
    write(stdlog,'(a,i4,1p,2e11.2)') &
      'dust_t%init: bin, s_cgs, m_cgs =', ib, self%s_cgs(ib), self%m_cgs(ib)
    allocate (self%exports(ib)%lock)
    call self%exports(ib)%lock%init ('exp')
    p => self%exports(ib)
    call p%lock%set ('init')
    call p%init (self%me, verbose)
    call p%lock%unset ('init')
    do i3=1,m(3)
    do i2=1,m(2)
    do i1=1,m(1)
      self%dd(i1,i2,i3,ib) = dust_to_gas*self%d(i1,i2,i3)
      self%a(:,i1,i2,i3) = 0.0              ! cell acceleration
      self%v(:,i1,i2,i3) = 0.0              ! cell velocity
      p => self%particles(i1,i2,i3,ib)
      call p%init (self%mp, verbose)
      if (verbose > 1) then
        o => self%previous(i1,i2,i3,ib)
        call o%init (self%mp, verbose)
      end if
      p%id = [self%patch%id,i1,i2,i3,ib]
      p%n = particles_per_bin
      p%lp = p%n
      w = self%dd(i1,i2,i3,ib)/(n_bin*p%n)
      do j=1,p%n
        p%r(j,:) = self%random%uniform(3)
        p%v(j,:) = 0.0
        p%q(j,:) = [i1,i2,i3]
        p%w(j) = w
        p%flag(j) = 0
      end do
    end do
    end do
    end do
  end do
  call trace%end ()
END SUBROUTINE initial_condition

!===============================================================================
!> Initialize particle solver, after nbor lists are created, making it possible
!> to create a copy of the patch nbor list, with tasks pointing to self
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(particle_solver_t):: self
  class(list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  if (.not.on) return
  call trace%begin ('particle_solver_t%init_task_list')
  self%task_list => task_list
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Update particles, using a series of operations.  Constraints: setup must be 
!> done before kdk, trade must be done after kdk, import must be done after 
!> trade and before pack, sumup is cheapest after pack, and renorm needs to be 
!> done after sumup.
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t):: self
  logical:: dbg
  real(8):: wc
  class(particles_t), pointer:: p
  if (on) then
    call trace%begin ('particle_solver_t%update')
    wc = wallclock()
    call self%output
    call self%setup               ! set up cell properties
    !call self%copy                ! copy particles
    call self%kdk                 ! move particles
    !call self%count
    call self%trade               ! trade btw cells and export
    call self%import              ! import particles from nbors
    call self%pack                ! remove array holes
    call self%sumup               ! compute densities
    call self%renorm              ! renormalize particles
    associate (t_particles => timer%t_particles)
    !$omp atomic update
    t_particles = t_particles + (wallclock()-wc)
    end associate
    call trace%end ()
  end if
END SUBROUTINE update

!===============================================================================
!> Initialize cell data
!===============================================================================
SUBROUTINE setup (self)
  class(particle_solver_t)                      :: self
  !.............................................................................
  integer                                       :: i1, i2, i3, ib, np
  real                                          :: c1, c2, c3
  real(8)                                       :: wc
  class(gpatch_t), pointer                      :: patch
  type(particles_t), pointer                    :: p
  real, pointer, contiguous, dimension(:,:,:)   :: p1, p2, p3
  integer, save                                 :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%setup', itimer=itimer)
  wc = omp_get_wtime()
  patch => self%patch
  self%dtime = patch%dtime
  self%time = patch%time + patch%dtime
  p1 => patch%mem(:,:,:,patch%idx%px,patch%it,1)
  p2 => patch%mem(:,:,:,patch%idx%py,patch%it,1)
  p3 => patch%mem(:,:,:,patch%idx%pz,patch%it,1)
  self%lb = patch%mesh%li
  self%ub = patch%mesh%ui
  np = 0
  c1 = 1.0/self%patch%ds(1)
  c2 = 1.0/self%patch%ds(2)
  c3 = 1.0/self%patch%ds(3)
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    self%v(1,i1,i2,i3) = p1(i1,i2,i3)/self%d(i1,i2,i3)*c1
    self%v(2,i1,i2,i3) = p2(i1,i2,i3)/self%d(i1,i2,i3)*c2
    self%v(3,i1,i2,i3) = p3(i1,i2,i3)/self%d(i1,i2,i3)*c3
    self%a(1,i1,i2,i3) = 0.0
    self%a(2,i1,i2,i3) = 0.0
    self%a(3,i1,i2,i3) = 0.0
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      np = np + p%n
    end do
  end do
  end do
  end do
  self%n = np
  associate (n_particles => timer%n_particles)
  !$omp atomic update
  n_particles = n_particles + np
  end associate
  wc = omp_get_wtime() - wc
  if (verbose > 1) then
    write (stdlog,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,i6)') &
      wc, wc/np, self%patch%id, np, 'setup', self%patch%id
  end if
  call trace%end (itimer)
END SUBROUTINE setup

!===============================================================================
!> Initialize cell data
!===============================================================================
SUBROUTINE copy (self)
  class(particle_solver_t)   :: self
  !.............................................................................
  integer                    :: i1, i2, i3, i4, j, k, m(4), q(3), np, np_new
  type(particles_t), pointer :: p1, p2, p3, new(:,:,:,:)
  real(8)                    :: wc
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%copy', itimer=itimer)
  wc = wallclock()
  m = shape(self%particles)
  allocate (new(m(1),m(2),m(3),m(4)))
  do i4=1,m(4)
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    p1 => self%particles(i1,i2,i3,i4)
    p2 => new(i1,i2,i3,i4)
    call p2%init (p1%np_new, self%verbose)
    p2%n = p1%np_new
  end do
  end do
  end do
  end do
  print '(6p,f10.3," mus/cell-alloc")', (wallclock()-wc)/product(m(1:3))
  wc = wallclock()
  do i4=1,m(4)
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    p1 => self%particles(i1,i2,i3,i4)
    np = p1%n
    do j=1,np
      q = p1%q(j,1:3)
      p2 => new(q(1),q(2),q(3),i4)
      k = p2%n + 1
      p2%r(k,1) = p1%r(j,1)
      p2%r(k,2) = p1%r(j,2)
      p2%r(k,3) = p1%r(j,3)
      p2%v(k,1) = p1%v(j,1)
      p2%v(k,2) = p1%v(j,2)
      p2%v(k,3) = p1%v(j,3)
      p2%q(k,1) = p1%q(j,1)
      p2%q(k,2) = p1%q(j,2)
      p2%q(k,3) = p1%q(j,3)
      p2%w(k)   = p1%w(j)
      p2%n = k
    end do
    call p1%dealloc
  end do
  end do
  end do
  end do
  deallocate (self%particles)
  self%particles => new
  print '(6p,f10.3," mus/cell-copy")', (wallclock()-wc)/product(m(1:3))
  call trace%end (itimer)
END SUBROUTINE copy

!===============================================================================
!> Kick-drift-kick update, with no hole in any cells
!===============================================================================
SUBROUTINE kdk (self)
  class(particle_solver_t)   :: self
  integer                    :: i1, i2, i3, ib, i, j, q(3), np
  real                       :: dt, drag_c, cs, c, vmin, vmax, cmin, cmax
  real                       :: d, e, pc, pwd, ds, cf
  real, dimension(3)         :: a, v, dv, pd
  real(8)                    :: wc
  class(gpatch_t), pointer   :: patch
  type(particles_t), pointer :: p, o
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%kdk', itimer=itimer)
  wc = omp_get_wtime()
  !-----------------------------------------------------------------------------
  patch => self%patch
  ds = patch%ds(1)
  dt = self%dtime*0.5d0
  drag_c = scaling%d/(density_yr_cgs*cgs%yr/scaling%t)
  !-----------------------------------------------------------------------------
  vmin = 0.0
  vmax = 0.0
  cmin = huge(1.)
  cmax = 0.0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    d = patch%mem(i1,i2,i3,patch%idx%d,patch%it,1)
    e = patch%mem(i1,i2,i3,patch%idx%e,patch%it,1)
    if (eos%gamma == 1d0) then
      cs = sqrt(e/d)/ds
    else
      cs = sqrt(e*eos%gamma*(eos%gamma-1d0)/d)/ds
    end if
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      np = p%n
      c = drag_c*d/self%s_cgs(ib)
      a = self%a(:,i1,i2,i3)
      v = self%v(:,i1,i2,i3)
      dt = self%dtime
      !-------------------------------------------------------------------------
      ! With gas backtreaction
      !-------------------------------------------------------------------------
      if (backreaction) then
        do j=1,np
          dv(1) = v(1)-p%v(j,1)                           ! speed difference
          dv(2) = v(2)-p%v(j,2)                           ! speed difference
          dv(3) = v(3)-p%v(j,3)                           ! speed difference
          pc = c*sqrt(1.0 &
              + (dv(1)**2+dv(2)**2+dv(3)**2)/cs**2)       ! supersonic Epstein
          pd(1) = pc*dv(1)                                ! drag force
          pd(2) = pc*dv(2)                                ! drag force
          pd(3) = pc*dv(3)                                ! drag force
          pwd = p%w(j)/d                                  ! density ratio
          cf = 1./(1.+pc*dt*(1.+pwd))                     ! common factor
          dv(1) = dt*(a(1)+pd(1))*cf                      ! velocity delta
          dv(2) = dt*(a(2)+pd(2))*cf                      ! velocity delta
          dv(3) = dt*(a(3)+pd(3))*cf                      ! velocity delta
          p%v(j,1) = p%v(j,1) + dv(1)                     ! particle update
          p%v(j,2) = p%v(j,2) + dv(2)                     ! particle update
          p%v(j,3) = p%v(j,3) + dv(3)                     ! particle update
          v(1) = v(1) - dv(1)*pwd                         ! gas update
          v(2) = v(2) - dv(2)*pwd                         ! gas update
          v(3) = v(3) - dv(3)*pwd                         ! gas update
          p%r(j,1) = p%r(j,1) + 2.*dt*p%v(j,1)            ! forward two dr
          p%r(j,2) = p%r(j,2) + 2.*dt*p%v(j,2)            ! forward two dr
          p%r(j,3) = p%r(j,3) + 2.*dt*p%v(j,3)            ! forward two dr
          pd(1) = pc*(v(1)-p%v(j,1))                      ! drag force
          pd(2) = pc*(v(2)-p%v(j,2))                      ! drag force
          pd(3) = pc*(v(3)-p%v(j,3))                      ! drag force
          dv(1) = dt*(a(1)+pd(1))*cf                      ! velocity delta
          dv(2) = dt*(a(2)+pd(2))*cf                      ! velocity delta
          dv(3) = dt*(a(3)+pd(3))*cf                      ! velocity delta
          p%v(j,1) = p%v(j,1) + dv(1)                     ! particle update
          p%v(j,2) = p%v(j,2) + dv(2)                     ! particle update
          p%v(j,3) = p%v(j,3) + dv(3)                     ! particle update
          v(1) = v(1) - dv(1)*pwd                         ! gas update
          v(2) = v(2) - dv(2)*pwd                         ! gas update
          v(3) = v(3) - dv(3)*pwd                         ! gas update
        end do
      !-------------------------------------------------------------------------
      ! No gas backtreaction
      !-------------------------------------------------------------------------
      else
        do j=1,np
          dv(1) = v(1)-p%v(j,1)                           ! speed difference
          dv(2) = v(2)-p%v(j,2)                           ! speed difference
          dv(3) = v(3)-p%v(j,3)                           ! speed difference
          pc = c*sqrt(1.0 &
              + (dv(1)**2+dv(2)**2+dv(3)**2)/cs**2)       ! supersonic Epstein
          dv(1) = dt*(a(1)+pc*dv(1))/(1.+pc*dt)           ! velocity delta
          dv(2) = dt*(a(2)+pc*dv(2))/(1.+pc*dt)           ! velocity delta
          dv(3) = dt*(a(3)+pc*dv(3))/(1.+pc*dt)           ! velocity delta
          p%v(j,1) = p%v(j,1) + dv(1)                     ! particle update
          p%v(j,2) = p%v(j,2) + dv(2)                     ! particle update
          p%v(j,3) = p%v(j,3) + dv(3)                     ! particle update
          p%r(j,1) = p%r(j,1) + 2.*dt*p%v(j,1)            ! forward two dr
          p%r(j,2) = p%r(j,2) + 2.*dt*p%v(j,2)            ! forward two dr
          p%r(j,3) = p%r(j,3) + 2.*dt*p%v(j,3)            ! forward two dr
          dv(1) = dt*(a(1)+pc*dv(1))/(1.+pc*dt)           ! velocity delta
          dv(2) = dt*(a(2)+pc*dv(2))/(1.+pc*dt)           ! velocity delta
          dv(3) = dt*(a(3)+pc*dv(3))/(1.+pc*dt)           ! velocity delta
          p%v(j,1) = p%v(j,1) + dv(1)                     ! particle update
          p%v(j,2) = p%v(j,2) + dv(2)                     ! particle update
          p%v(j,3) = p%v(j,3) + dv(3)                     ! particle update
        end do
      end if
      !-------------------------------------------------------------------------
      ! Update cell address and save a flag indicating motion
      !-------------------------------------------------------------------------
      do j=1,np
        q(1) = int(1.+p%r(j,1))-1
        q(2) = int(1.+p%r(j,2))-1
        q(3) = int(1.+p%r(j,3))-1
        p%r(j,1) = p%r(j,1) - q(1)
        p%r(j,2) = p%r(j,2) - q(2)
        p%r(j,3) = p%r(j,3) - q(3)
        p%q(j,1) = p%q(j,1) + q(1)
        p%q(j,2) = p%q(j,2) + q(3)
        p%q(j,3) = p%q(j,3) + q(2)
        p%flag(j) = abs(q(1)) + abs(q(2)) + abs(q(3))
      end do
    end do
  end do
  end do
  end do
  !.............................................................................
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (stdlog,'("KDK",6p,f8.3,i6,i10,1p,4e9.1)') &
      wc/self%n, self%patch%id, self%n
  end if
  call trace%end (itimer)
END SUBROUTINE kdk

!===============================================================================
!> Count new particle numbers, either directly, or differentially
!===============================================================================
SUBROUTINE count (self)
  class(particle_solver_t):: self
  integer:: i1, i2, i3, i4, j, ii(3), qq(3)
  type(particles_t), pointer :: p
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%count1', itimer=itimer)
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    if (count_all) then
      do i4=1,self%n_bin
        p => self%particles(i1,i2,i3,i4)
        do j=1,p%n
          ii = p%q(j,1:3)
          self%particles(ii(1),ii(2),ii(3),i4)%np_new = &
          self%particles(ii(1),ii(2),ii(3),i4)%np_new + 1
        end do
      end do
    else
      !-------------------------------------------------------------------------
      ! Differential counts
      !-------------------------------------------------------------------------
      ii = [i1,i2,i3]
      do i4=1,self%n_bin
        p => self%particles(i1,i2,i3,i4)
        do j=1,p%n
          qq = p%q(j,1:3)
          if (any(qq /= ii)) then
            self%particles(ii(1),ii(2),ii(3),i4)%np_new = &
            self%particles(ii(1),ii(2),ii(3),i4)%np_new - 1
            self%particles(qq(1),qq(2),qq(3),i4)%np_new = &
            self%particles(qq(1),qq(2),qq(3),i4)%np_new + 1
          end if
        end do
      end do
    end if
  end do
  end do
  end do
  call trace%end (itimer)
END SUBROUTINE count

!===============================================================================
!> Re-assign particles to other cells.  New particles will be added at the
!> end, and holes will only appear for indices up to the initial array lenght
!===============================================================================
SUBROUTINE trade (self)
  class(particle_solver_t)                    :: self
  integer                                     :: i1, i2, i3, ip, ib, i, j, ii(3)
  integer                                     :: q1, q2, q3, pq(3)
  integer                                     :: j2, ne, nm, ne1, nm1
  integer                                     :: n3, n4, n5, me, mp
  real(8)                                     :: wc
  real                                        :: dt
  type(particles_t), pointer                  :: p, pn, pe, o
  real, pointer, contiguous, dimension(:,:)   :: pv, pr
  real, pointer, contiguous, dimension(:,:,:) :: p1, p2, p3
  logical                                     :: export, move
  logical, save                               :: warn(2) = .true.
  integer, save                               :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%trade', itimer=itimer)
  wc = omp_get_wtime()
  n3 = 0
  n4 = 0
  n5 = 0
  nm = 0
  ne = 0
  do ib=1,n_bin
    pe => self%exports(ib)
    call pe%lock%set ('exp')
    if (verbose > 2) then
      write (stdlog,'(i6,2x,a,3i10)') &
      self%patch%id, 'exports before:', pe%n, pe%lp, pe%mp
    end if
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      ii = [i1,i2,i3]
      p => self%particles(i1,i2,i3,ib)
      if (verbose > 2) then
        call  p%validate ('p  before')
        call pe%validate ('pe before')
      end if
      if (p%n == 0) write(stdlog,'(a,6i5)') '0p bef trade', p%id, p%lp
      ne1 = 0
      nm1 = 0
      if (verbose > 3) then
        write (stdlog,'(a,i6,4i4,6i6)') 'trade: id,i1,i2,i3,ib,nm,ne,nlp,nle=', &
        self%patch%id,i1,i2,i3,ib,nm1,ne1,p%n,p%lp,pe%n,pe%lp
      end if
      do j=1,p%lp
        if (p%w(j) > 0.0 .and. p%flag(j) /= 0) then
          !---------------------------------------------------------------------
          ! If the position index is outside, move the particle to exports
          !---------------------------------------------------------------------
          pq(1) = p%q(j,1)
          pq(2) = p%q(j,2)
          pq(3) = p%q(j,3)
          export = .false.
          move = .false.
          do i=1,3
            if (pq(i) /= ii(i)) then
              if (pq(i) < self%lb(i) .or. pq(i) > self%ub(i)) then
                export = .true.
                exit
              end if
              move = .true.
            end if
          end do
          if (export) then
            ne1 = ne1 + 1
            if (pe%lp >= pe%mp) then
              me = 1.2*pe%mp
              call pe%resize(me)
              self%me = me
              if (verbose > 2) then
                write (stdlog,*) &
                  'export resize: mo, mp =', me, self%me
              end if
              n3 = n3+1
            end if
            pe%n = pe%n + 1           ! increment export count
            pe%lp = pe%lp + 1           ! increment export length
            j2 = pe%lp
            do i=1,3
              pe%r(j2,i) = p%r(j,i)
              pe%v(j2,i) = p%v(j,i)
              pe%q(j2,i) = p%q(j,i)
            end do
            pe%w(j2) = p%w(j)
            p%w(j) = 0.0                ! mark position empty
            p%n = p%n - 1             ! decrement actual count
            if (j == p%lp) then         ! if last slot
              p%lp = p%lp - 1           ! also decrement length
            end if
            if (verbose > 2) then
              call  p%validate ('p0 after ')
              call pe%validate ('pe after ')
              if (verbose > 3) then
                write (stdlog,'(i6,":",3i4,2i6,2x,a,i6,3i4,2i6,1p,e10.2)') &
                j, i1, i2, i3, p%n, p%lp, 'export', j2, pq, pe%n, pe%lp, pe%w(j2)
              end if
            end if
            !-------------------------------------------------------------------
            ! If the position index has changed, add to the target cell and
            ! remove from the source cell
            !-------------------------------------------------------------------
          else if (move) then
            nm1 = nm1 + 1
            !q1 = pq(1)
            !q2 = pq(2)
            !q3 = pq(3)
            pn => self%particles(pq(1),pq(2),pq(3),ib)
            pn%lp = pn%lp + 1
            j2 = pn%lp
            !--- j2 is now the index of the cell to fill ---
            if (j2 >= p%mp) then
              mp = 1.2*pn%mp
              call pn%resize(mp)
              pn%mp = mp
              if (verbose > 2) then
                write (stdlog,*) &
                  'particles resize: mo, mp =', mp, self%mp
              end if
              n4 = n4+1
            end if
            if (j2 < pn%mp) then
              pn%n = pn%n + 1             ! increment target count
              do i=1,3
                pn%r(j2,i) = p%r(j,i)
                pn%v(j2,i) = p%v(j,i)
                pn%q(j2,i) = p%q(j,i)
              end do
              pn%flag(j2) = p%flag(j)
              pn%w(j2) = p%w(j)
              p%w(j) = 0.0                ! mark position empty
              p%n = p%n - 1               ! decrement source count
              if (j == p%lp) then         ! if last slot
                p%lp = p%lp - 1           ! also decrement length
              end if
            else
              call io%abort ('should not happen')
            end if
            if (verbose > 2) then
              call  p%validate ('p1 after ')
              call pn%validate ('pn after ')
            end if
          end if
        end if
      end do
      if (p%n <= 0.5*min_fraction*particles_per_bin) then
        n5 = n5+1
      end if
      ne = ne + ne1
      nm = nm + nm1
      if (verbose > 2) then
        call pe%validate ('pe')
      end if
    end do
    end do
    end do
    if (verbose > 2) then
      write (stdlog,'(i6,2x,a,3i10)') &
      self%patch%id, 'exports  after:', pe%n, pe%lp, pe%mp
    end if
    call pe%lock%unset ('exp')
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (stdlog,'("TRADE",6p,f8.3,i6,8i10)') &
      wc/self%n, self%patch%id, self%n, nm, ne, n3, n4, n5
  end if
  call trace%end (itimer)
END SUBROUTINE trade

!===============================================================================
!> Import particles from other patch's exports
!===============================================================================
SUBROUTINE import (self)
  class(particle_solver_t)   :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: nbtask
  class(particles_t), pointer:: pe, p
  real(8)                    :: wc, r(3), dist(3)
  integer                    :: ni, ne, ib, n, j, jp, i, q(3), np
  integer                    :: i1, i2, i3
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%import' ,itimer=itimer)
  if (verbose > 2) then
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      do ib=1,n_bin
        p => self%particles(i1,i2,i3,ib)
        write(stdlog,'(a,i6,4i4,3i6)') &
        'after import: i1,i2,i3,ib,np,lp,mp0', &
        self%patch%id,i1,i2,i3,ib,p%n,p%lp,p%mp
        call pe%validate ('before import')
      end do
    end do
    end do
    end do
  end if
  ni = 0
  ne = 0
  wc = omp_get_wtime()
  nbor => self%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (particle_solver_t)
      dist = nbtask%patch%task_t%distance(self%patch)/nbtask%patch%ds
      !-------------------------------------------------------------------------
      ! Loop over particles in export list and grab those that belong here
      !-------------------------------------------------------------------------
      n = 0
      do ib=1,n_bin
        pe => nbtask%exports(ib)
        call pe%lock%set ('import')
        j = 1
        do while (j < pe%lp)
          if (pe%w(j) > 0.0) then
            r = dist + (pe%r(j,:)+pe%q(j,:))
            q = floor(r)
            if (all(q >= self%lb .and. q <= self%ub)) then
              p => self%particles(q(1),q(2),q(3),ib)
              if (p%lp >= p%mp) then
                call p%pack (np)
                if (verbose > 2) call p%validate ('p forced')
              end if
              if (p%lp < p%mp) then
                n = n  + 1
                jp = p%lp + 1
                do i=1,3
                  p%r(jp,i) = pe%r(j,i)
                  p%v(jp,i) = pe%v(j,i)
                  p%q(jp,i) = q(i)
                end do
                p%w(jp) = pe%w(j)
                p%lp = jp
                p%n = p%n+1
                pe%w(j) = 0.0
                pe%n = pe%n - 1
              else
                ne = ne + 1
              end if
            end if
          end if
          j = j + 1
        end do
        call pe%lock%unset ('import')
      end do
      ni = ni + n
    end select
    nbor => nbor%next
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 2) then
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      do ib=1,n_bin
        p => self%particles(i1,i2,i3,ib)
        write(stdlog,'(a,i6,4i4,3i6)') &
        'after import: i1,i2,i3,ib,np,lp,mp0', &
        self%patch%id,i1,i2,i3,ib,p%n,p%lp,p%mp
        call p%validate ('after import')
      end do
    end do
    end do
    end do
  end if
  if (verbose > 1) then
    write (stdlog,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,2i10)') &
      wc, wc/self%n, self%patch%id, self%n, 'import', ni, ne
  end if
  call trace%end (itimer)
END SUBROUTINE import

!===============================================================================
!> Remove the holes in the arrays
!===============================================================================
SUBROUTINE pack (self)
  class(particle_solver_t)   :: self
  integer                    :: i1, i2, i3, ib, n, lp(2), np(2), mp(2)
  real(8)                    :: wc
  type(particles_t), pointer :: p
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pack', itimer=itimer)
  !-----------------------------------------------------------------------------
  ! Cell arrays
  !-----------------------------------------------------------------------------
  wc = omp_get_wtime()
  n = 0
  lp = 0
  np = 0
  !-----------------------------------------------------------------------------
  ! Pack particles
  !-----------------------------------------------------------------------------
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      p%id = [self%patch%id,i1,i2,i3,ib]
      np(1) = max(np(1),p%n)
      lp(1) = max(lp(1),p%lp)
      if (p%n == 0) write(stdlog,'(a,7i5,3f6.1)') &
        '0p bef pack ',p%id,p%n,p%lp, self%v(:,i1,i2,i3)
      call p%pack (n)
      !if (p%n == 0) write(stdlog,'(a,7i5,3f6.1)') &
        !'0p aft pack ',p%id,p%n,p%lp, self%v(:,i1,i2,i3)
      if (verbose > 2) call p%validate ('normal p')
    end do
  end do
  end do
  end do
  mp(1) = p%mp
  !-----------------------------------------------------------------------------
  ! Pack exports
  !-----------------------------------------------------------------------------
  do ib=1,n_bin
    p => self%exports(ib)
    call p%lock%set ('pack')
    p%id = [self%patch%id,0,0,0,ib]
    np(2) = max(np(2),p%n)
    lp(2) = max(lp(2),p%lp)
    call p%pack (n)
    call p%lock%unset ('pack')
    if (verbose > 2) call p%validate ('normal pe')
  end do
  mp(2) = p%mp
  wc = omp_get_wtime() - wc
  if (verbose > 0) then
    write (stdlog,'("PACK",6p,f8.3,i6,8i10)') &
      wc/self%n, self%patch%id, self%n, n, np, lp, mp
  end if
  call trace%end (itimer)
END SUBROUTINE pack

!===============================================================================
!> Sum up density
!===============================================================================
SUBROUTINE sumup (self)
  class(particle_solver_t)    :: self
  integer                     :: i1, i2, i3, ib, j, nmax
  real                        :: dd
  real(8)                     :: wc
  type(particles_t), pointer  :: p
  integer, save               :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%sumup', itimer=itimer)
  wc = omp_get_wtime()
  nmax = 0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      nmax = max(nmax,p%lp)
      dd = 0.0
      do j=1,p%lp
        dd = dd + p%w(j)
      end do
      self%dd(i1,i2,i3,ib) = dd
      if (dd == 0.0) then
        write (stdlog,'(a,i6,4i4,2i6,1p,4e10.2)') &
        'sumup 0d: id,np,lp,v=', &
        p%id, p%n, p%lp, self%v(:,i1,i2,i3)*self%dtime
      end if
    end do
  end do
  end do
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 1) then
    write (stdlog,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,i10)') &
      wc, wc/self%n, self%patch%id, self%n, 'density', nmax
  end if
  call trace%end (itimer)
END SUBROUTINE sumup

!===============================================================================
!> Add and remove particles where needed
!===============================================================================
  SUBROUTINE renorm (self)
  class(particle_solver_t)    :: self
  integer                     :: i1, i2, i3, ib, i, j, n, na, nr, np, lp
  real                        :: n_min, n_max, w
  real(8)                     :: t, wc
  type(particles_t), pointer  :: p
  integer, save               :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%renorm', itimer=itimer)
  wc = omp_get_wtime()
  n_min = particles_per_bin*min_fraction
  n_max = particles_per_bin*max_fraction
  na = 0
  nr = 0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
  do i1=self%lb(1),self%ub(1)
    do ib=1,n_bin
      p => self%particles(i1,i2,i3,ib)
      9 format(3f8.3,i4)
      call io%assert (p%n==p%lp, 'The state must be packed before renorm')
      if (self%dd(i1,i2,i3,ib) <= 0.0) then
        write (stdlog,'(a,i6,4i3,2i4,1p,5e9.1,i10)') &
        'zero density warning: id,np,lp,v=', &
        p%id,p%n,p%lp, self%v(:,i1,i2,i3), self%d(i1,i2,i3), &
        self%patch%time, self%patch%istep
        cycle
      end if
      !-------------------------------------------------------------------------
      ! If too few particles, add new ones and set new weights to all
      !-------------------------------------------------------------------------
      if (p%n < n_min) then
        if (verbose > 1) &
          write (stdlog,'(a,i6,4i4,2x,i4,f6.1)') &
            'add', self%patch%id, i1, i2, i3, ib, p%n, n_min
        w = self%dd(i1,i2,i3,ib)/particles_per_bin
        t = self%time
        !call p%pack (n)
        do j=p%n+1,particles_per_bin
          do i=1,3
            p%r(j,i) = self%random%ran1()
            p%v(j,i) = self%v(i,i1,i2,i3)
          end do
          p%q(j,:) = [i1,i2,i3]
          p%w(j) = w
        end do
        na = na + particles_per_bin - p%n
        p%n = particles_per_bin
        p%lp = particles_per_bin
        if (verbose > 2) call p%validate ('renorm add')
      !-------------------------------------------------------------------------
      ! If too many, renormalize and remove extraneous
      !-------------------------------------------------------------------------
      else if (p%n > n_max) then
        if (verbose > 1) &
          write (stdlog,'(a,i6,4i4,2x,i4,f6.1)') &
            'rem', self%patch%id, i1, i2, i3, ib, p%n, n_max
        nr = nr + (p%n - particles_per_bin)
        j = p%lp
        np = p%n
        lp = p%lp
        do while (p%n > particles_per_bin)
          if (j <= 0) then
            write (stdlog,'(a,i6,4i4,2x,4i4)') &
              'error:', self%patch%id, i1, i2, i3, ib, np, lp, p%n, p%lp
            write (stderr,*) 'ERROR'
            stop
          end if
          if (verbose > 3) then
            write (stdlog,'("rem:",3i4,1p,e10.2 )') j, p%n, p%lp, p%w(j)
          end if
          if (p%w(j) == 0.0) then                   ! if a hole
            p%lp = p%lp - 1                         ! shorten length
          else                                      ! if a particle
            p%n = p%n - 1                         ! decrement count
            p%lp = p%lp - 1                         ! shorten length
          end if
          j = j-1
        end do
        w = self%dd(i1,i2,i3,ib)/particles_per_bin  ! new weight
        !call p%pack (n)                             ! pack and set weight
        do j=1,p%n                                 ! assumes packed state!
          p%w(j) = w
        end do
        p%n = particles_per_bin
        p%lp = particles_per_bin
        if (verbose > 2) call p%validate ('renorm rem')
      end if
    end do
  end do
  end do
  end do
  wc = omp_get_wtime() - wc
  if (verbose > 1) then
    write (stdlog,'(f6.3," s",6p,f8.3," mus/p",i6,i10,a10,2i10)') &
      wc, wc/self%n, self%patch%id, self%n, 'renorm', na, nr
  end if
  call trace%end (itimer)
END SUBROUTINE renorm

!===============================================================================
!> I/O
!===============================================================================
SUBROUTINE output (self)
  class(particle_solver_t):: self
  !.............................................................................
  class(particles_t), pointer:: p
  class(mesh_t), pointer:: mesh(:)
  real, pointer, contiguous:: tmp_r(:,:)
  integer, pointer, contiguous:: tmp_i(:)
  integer:: i, j, n, i1, i2, i3, ib
  !-----------------------------------------------------------------------------
  call trace%begin ('dust_t%output')
  if (verbose > 2) then
    write (stdlog, '(a,i6,2f8.4)') &
      'dust_t%output: id, time, aux_next =', &
      self%patch%id, self%time, self%aux%aux_next
  end if
  if (self%time >= self%aux%aux_next) then
    call self%aux%register ('dust_d', self%dd)
    !...........................................................................
    if (verbose > 2) then
      do i=1,n_bin
        write (stdlog, '(a,i6, i3,2x,1p,10e11.2)') &
          'dust_t%output: id, bin, n =', self%patch%id, i, self%dd(3,3,1:10,i)
      end do
    end if
    !...........................................................................
    if (particle_io) then
      mesh => self%patch%mesh
      do i3=self%lb(3),self%ub(3)
      do i2=self%lb(2),self%ub(1)
      do i1=self%lb(1),self%ub(1)
        do ib=1,n_bin
          n = 0
          p => self%particles(i1,i2,i3,ib)
          do j=1,p%n
            if (p%flag(j) > 0) then
             n = n+1
            end if
          end do
          allocate (tmp_r(5,n))
          i = 0
          do j=1,p%n
            if (p%flag(j) > 0) then
              i = i+1
              tmp_r(1:3,i) = mesh%p + (p%q(j,:) + p%r(j,:) - mesh%o)*mesh%d
              tmp_r(  4,i) = p%w(j)
              tmp_r(  5,i) = self%s_cgs(ib)
            end if
          end do
        end do
      end do
      end do
      end do
      call self%aux%register ('dust_p' , tmp_r)
      call self%aux%output_sequence (self%patch%id, self%time)
      if (verbose > 1) then
        write (stdlog, '(a,i6,f8.4,i3,2i8)') &
         'dust_t%output: id, particle array =', &
          self%patch%id, self%time, shape(tmp_r), i
      end if
      deallocate (tmp_r)
    else
      call self%aux%output_sequence (self%patch%id, self%time)
    end if
    if (verbose > 2) then
      write (stdlog, '(a,i6,1p,10e11.2)') &
       'dust_t%output: self%aux%aux_next =', &
        self%patch%id, self%aux%aux_next
    end if
  end if
  call trace%end ()
END SUBROUTINE output

END MODULE particle_solver_mod
