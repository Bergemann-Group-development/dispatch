!==============================================================================
!> Particle data types for dust-gas interaction:
!>
!>   particle_t         : data type for a single particle
!>   particles_t        : data type for lists of particles
!>   particle_solver_t  : data type for particle solvers
!>
!> The particle solver advances the particles past the time of the patch they
!> belong to, with timesteps that move them at most one cell forward. The time
!> step is thus generally longer than even the local Courant time for the cell,
!> and thus normally longer than the Courant time for the whole patch.  But if
!> the particle Courant factor is chosen small, it is possible that one needs
!> to advance the particle several time steps before its time exceeds the
!> partch time.
!>
!> The particle solver operates in integer cell coordinates, and thus with
!> speeds measured in cells per code time unit, instead of code lengths per
!> time unit.
!>
!> The particle solver is either called as a post_update() of a patch solver,
!> or else as an independent task, which depends on the patch task it is
!> associated with.
!===============================================================================
!define DBG
MODULE particle_solver_mod
  USE io_mod
  USE trace_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE task_mod
  USE bits_mod
  USE mesh_mod
  USE link_mod
  USE list_mod
  USE gpatch_mod
  USE particles_mod
  USE dust_mod
  USE scaling_mod
  USE units_mod
  USE timer_mod
  USE omp_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Everything needed by a particle solver
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: particle_solver_t
    type(particles_t), pointer :: particles(:,:,:,:)   ! particles in cells
    type(particles_t), pointer :: exports(:)           ! particles to export
    class(gpatch_t),   pointer :: patch => null()      ! pointer to patch
    class(nbor_t),     pointer :: nbor => null()       ! nbor list
    type(link_t)   ,   pointer :: link => null()       ! task link
    type(dust_t)               :: dust                 ! dust solver
    real, pointer, contiguous  :: dd(:,:,:,:) => null()! number density per bin
    real, pointer, contiguous  :: s_cgs(:)             ! particle size
    real, pointer, contiguous  :: m_cgs(:)             ! particle size
    logical                    :: forced = .false.     ! external force
    logical, dimension(3)      :: lower                ! has lower boundary
    logical, dimension(3)      :: upper                ! has upper boundary
    integer, dimension(3)      :: lb=0, ub=0           ! lower and upper index bounds
    integer                    :: np = 0               ! particles
    integer                    :: ni = 0               ! particles imported
    integer                    :: no = 0               ! particles outside
    real                       :: v_min, v_max         ! min/max u_phase
    class(list_t), pointer     :: task_list            ! task list
    !type(lock_t), pointer      :: lock                 ! export lock
    procedure(pforce), pointer, nopass :: pforce => null() ! point force
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: initial_condition
    procedure:: init_task_list
    procedure:: update
    procedure:: import
    procedure:: output
  end type
  type(particles_t), pointer :: pool          ! particle pool
  !$omp threadprivate (pool)                  ! for each thread
  integer:: verbose = 0                       ! noise level
  integer:: n_bin = 1                         ! number of size bins
  integer:: seed = 11                         ! random seed
  integer:: particles_per_bin = 10            ! particles per cell-bin
  integer:: output_every=1000                 ! output only that fraction
  integer:: ip_dbg = 0                        ! particle ID to debug
  integer:: id_dbg = 0                        ! patch ID to debug
  integer:: i_dbg(5) = 0                      ! cell index triplet to debug
  real:: density_yr_cgs = 3e-12               ! drag scaling constant
  real:: solid_density_cgs = 2                ! dust solid density
  real:: dust_to_gas = 0.01                   ! dust-to-gas ratio
  real:: max_size_cgs = 1e-5                  ! max particle size
  real:: min_size_cgs = 1e-5                  ! min particle size
  real:: min_fraction = 0.7                   ! lower limit
  real:: max_fraction = 1.4                   ! upper limit
  real:: drag_coeff                           ! scaled value
  logical:: particle_io = .false.             ! output particle positions?
  logical:: detailed_timer = .false.
  logical:: on = .true.
CONTAINS

!===============================================================================
!> Initialize particle solver, reading parameters and adding task link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(particle_solver_t), target  :: self
  class(gpatch_t), pointer          :: patch
  integer                           :: bin, iostat
  real                              :: solid_density, dln_m
  logical, save                     :: first_time=.true.
  namelist /particle_solver_params/ &
    on, verbose, detailed_timer, particle_io, &
    n_bin, particles_per_bin, seed, ip_dbg, id_dbg, i_dbg, &
    min_fraction, max_fraction, max_size_cgs, min_size_cgs, &
    density_yr_cgs, solid_density_cgs, dust_to_gas
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pre_init')
  self%patch => patch                           ! make the patch known
  self%id = patch%id                            ! for now, inherit
  self%dust%patch => self%patch                 ! make the patch known
  call io%lock%set ('particle_solver_t%pre_init')
  if (first_time) then
    rewind (io%input)
    read (io%input, particle_solver_params, iostat=iostat)
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%namelist, particle_solver_params)
    first_time = .false.
  end if
  call io%lock%unset ('particle_solver_t%pre_init')
  self%dust%detailed_timer = detailed_timer
  self%dust%ip_dbg = ip_dbg
  self%dust%i_dbg = i_dbg
  !-----------------------------------------------------------------------------
  ! Pre-factor for drag force.  The input parameter density_yr_cgs is the
  ! density (in CGS) for which the stopping time in 1 yr for a 1 cm particle.
  ! Since we keep size in cm, we just need to multiply with the density scaling
  ! constant, and divide with the value of one year in code units.
  !-----------------------------------------------------------------------------
  drag_coeff = scaling%d/(density_yr_cgs*cgs%yr/scaling%t)
  !-----------------------------------------------------------------------------
  ! Scaling to code units (keep dust size and mass in cgs)
  !-----------------------------------------------------------------------------
  solid_density = solid_density_cgs/(scaling%m/scaling%l**3)
  allocate (self%s_cgs(n_bin), self%m_cgs(n_bin))
  dln_m = log(max_size_cgs/min_size_cgs)/n_bin
  do bin=1,n_bin
    self%s_cgs(bin) = exp(log(min_size_cgs) + (bin-0.5)*dln_m)
    self%m_cgs(bin) = solid_density_cgs*self%s_cgs(bin)**3
    write(stdlog,'(a,i4,1p,2e11.2)') &
      'dust_t%init: bin, s_cgs, m_cgs =', bin, self%s_cgs(bin), self%m_cgs(bin)
  end do
  !-----------------------------------------------------------------------------
  ! Add a task link, in preparation for init_task_list
  !-----------------------------------------------------------------------------
  allocate (self%link)
  self%link%task => self
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE init (self)
  class(particle_solver_t) :: self
  class(gpatch_t), pointer :: pa
  integer                  :: m(3), i1, i2, i3, ib
  !-----------------------------------------------------------------------------
  ! Initialize particle lists and dust density array
  !-----------------------------------------------------------------------------
  ! if (.not.on) return
  call trace%begin ('particle_solver_t%init')
  m = self%patch%mesh%gn
  allocate (self%particles(m(1),m(2),m(3),n_bin)) 
  call io%bits_mem (storage_size(self%particles), &
                   product(shape(self%particles)), 'particles')
  allocate (pool, self%exports(n_bin))
  call io%bits_mem (storage_size(pool), omp%nthreads+n_bin, 'exports')
  call pool%init ([self%patch%id,0,0,0,0])
  do ib=1,n_bin
    call self%exports(ib)%init ([self%patch%id,0,0,0,ib])
  end do
  !-----------------------------------------------------------------------------
  allocate (self%dd(m(1),m(2),m(3),n_bin))
  call io%bits_mem (storage_size(self%dd), product(shape(self%dd)), 'dd')
  !-----------------------------------------------------------------------------
  ! Allocate force array
  !-----------------------------------------------------------------------------
  if (.not.associated(self%patch%force_per_unit_volume)) then
    allocate (self%patch%force_per_unit_volume(m(1),m(2),m(3),3))
    call io%bits_mem (storage_size(self%patch%force_per_unit_volume), &
                     product(shape(self%patch%force_per_unit_volume)), &
                     'force')
  end if
  if (.not.associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init ('part')
  end if
  call self%dust%init (n_bin, self%m_cgs, self%dd)
  call self%aux%init
  call self%random%init (seed)        ! FIXME: check seed inheritance
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particle_solver_t)    :: self
  class(particles_t), pointer :: particles
  type(particle_t), pointer   :: pn
  type(vparticle_t), pointer  :: pv
  class(particle_t), pointer  :: p
  class(gpatch_t), pointer    :: patch
  integer                     :: m(3), i1, i2, i3, bin, i, j, n, l, bytes, id(5)
  real                        :: d, w, vg(3)
  integer, save               :: itimer=0, extra_bits
  !-----------------------------------------------------------------------------
  ! Set flags indicating upper and lower boundaries.  The box size is self%box,
  ! and the patch position, size, and cell size are self%position, self%size,
  ! and self%ds, respectively.
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%initial_condition', itimer=itimer)
  self%upper = .false.
  self%lower = .false.
  patch => self%patch
  do i=1,3
    self%lb = patch%mesh%li
    self%ub = patch%mesh%ui
    if (.not.patch%periodic(i)) then
      if (patch%position(i) - patch%size(i)/2.0 - patch%ds(i) < -patch%box(i)/2.0) then
        self%lower(i) = .true.
        self%lb(i) = patch%mesh(i)%lb
      end if
      if (patch%position(i) + patch%size(i)/2.0 + patch%ds(i) > +patch%box(i)/2.0) then
        self%upper(i) = .true.
        self%ub(i) = patch%mesh(i)%ub
      end if
    end if
  end do
  self%dust%lower = self%lower
  self%dust%upper = self%upper
  self%dust%lb = self%lb
  self%dust%ub = self%ub
  !-----------------------------------------------------------------------------
  extra_bits = 200
  if (self%patch%id == 1) then
    allocate (p)
    call p%init
    bytes = (storage_size(p)+extra_bits)/8
    print '(a,i3,a,f18.3,a)', ' estimated memory for particles (', &
      bytes,' bytes):',bytes*n_bin*particles_per_bin/1024.**3 &
      *product(self%patch%mesh%b/self%patch%mesh%d), ' GB'
    call pool%append (p)
  end if
  !-----------------------------------------------------------------------------
  ! Initialize dust.  At boundaries, add particles also in guard zones
  !-----------------------------------------------------------------------------
#ifdef __GFORTRAN__
  !$omp critical (alloc_cr)
#endif
  l = 0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
    n = 0
    do i1=self%lb(1),self%ub(1)
      d = patch%mem(i1,i2,i3,patch%idx%d,patch%it,1)
      w = d*dust_to_gas/(n_bin*particles_per_bin)
      vg = patch%mem(i1,i2,i3,patch%idx%px:patch%idx%pz,patch%it,1) &
        / (patch%ds*d)
      do bin=1,n_bin
        particles => self%particles(i1,i2,i3,bin)
        id = [self%patch%id,i1,i2,i3,bin]
        call particles%init (id)
        n = n+2
        self%dd(i1,i2,i3,bin) = 0.0
        do j=1,particles_per_bin
          n = n+1
          l = l+1
          if (particle_io .and. mod(l,output_every)==0) then
            allocate(pv)
            p => pv
          else
            allocate(pn)
            p => pn
          end if
          call p%init
          p%q = [i1,i2,i3,bin]
          p%r = self%random%uniform (3)
          p%v = 0.0
          p%w = w
          p%t = 0
          self%dd(i1,i2,i3,bin) = self%dd(i1,i2,i3,bin) + w
          call particles%append (p)
        end do
      end do
      if (verbose > 0) call particles%validate ('IC')
    end do
    call io%bits_mem (storage_size(p)+extra_bits, n, 'p')
  end do
  end do
#ifdef __GFORTRAN__
  !$omp end critical (alloc_cr)
#endif
  !.............................................................................
  if (verbose > 0) then
    if (verbose > 1) then
      write (io%output,1) &
        'init: id, box, pos, size, lower, upper, n =', &
        patch%id, patch%box, patch%position, patch%size, &
        self%lower, self%upper, &
        real(self%np)
      1 format(a,i4,3(1x,3f7.3),1x,2(1x,3l1),1p,e11.2)
    end if
    write (stdlog,1) &
      'init: id, box, pos, size, lower, upper, n =', &
      patch%id, patch%box, patch%position, patch%size, &
      self%lower, self%upper, &
      real(self%np)
  else
    write(stdlog,'(a,i6,2(3f8.3,2x),1x,2(1x,3l1),1p,e11.2)') &
     'particle_solver_t%init: id, pos, box, lower, upper, n =', &
      patch%id, patch%position, patch%box, self%lower, self%upper, real(self%np)
  end if
  call self%dust%consistency ('IC')
  call trace%end (itimer)
END SUBROUTINE initial_condition

!===============================================================================
!> Initialize particle solver, after nbor lists are created, making it possible
!> to create a copy of the patch nbor list, with tasks pointing to self
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(particle_solver_t) :: self
  class(list_t), pointer   :: task_list
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init_task_list')
  self%task_list => task_list
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t)    :: self
  !.............................................................................
  real, dimension(3)          :: r, dr, v
  class(particle_t), pointer  :: p, next
  real(8)                     :: time(2), self_time, wc, start
  real                        :: ds, dt, os, n_min, n_max, ods
  integer                     :: ip, n, np, no, nk, nm, it, m(3), ii(3), dir
  integer(kind=2)             :: i, i1, i2, i3, bin, j, l(3), u(3), q(3)
  integer                     :: id, paid(5), istep
  integer(8)                  :: pid
  logical                     :: move, export
  real, pointer, contiguous   :: vg(:,:,:,:)   ! gas velocity
  class(gpatch_t), pointer    :: patch
  type(particles_t), pointer  :: exports(:), exportp
  type(particles_t), pointer  :: particles, bins(:)
  real                        :: d, up, pc, pd(3), pa(3), w
  real(8)                     :: times(9), wc1, wc2, sum
  integer, save               :: itimer=0
!-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%update', itimer=itimer)
  call self%output
  start = wallclock()
  patch => self%patch
  id = patch%id
  istep = patch%istep
  if (id == i_dbg(1)) then
    patch%dtime = max(patch%dtime,1d-10)
    self%dtime = max(self%dtime,1d-10)
    write (stdlog,'(a,2(1p,4e10.3,2x))') &
      'particle_solver_t%update begin: t, dt, t+dt =', &
      patch%time, patch%dtime, 1.0/patch%dtime, patch%time+patch%dtime, &
      self%time, self%dtime, 1.0/self%dtime, self%time+self%dtime
  end if
  self%forced = associated(self%pforce)
  !-----------------------------------------------------------------------------
  ! Only evolve particles if the last particle is behind the patch%time (which
  ! has not been updated yet, so we add the patch timestep)
  !-----------------------------------------------------------------------------
  self_time = patch%time + patch%dtime
  dt = patch%dtime/2d0
  self%v_min = patch%fminval(patch%u_phase)/patch%ds(1)
  self%v_max = patch%fmaxval(patch%u_phase)/patch%ds(1)
  if (verbose > 1) then
    write (stdlog,'(a,i6,1p,4e12.4)') &
      'particle_solver_t%update: id, time, dtime, min/max =', &
      patch%id, patch%time, patch%dtime, &
      patch%fmaxval(patch%u_phase), patch%fminval(patch%u_phase)
  end if
  !call self%particles%lock%set ('update')
  time(1) = wallclock()
  !-----------------------------------------------------------------------------
  ods = 1.0/patch%ds(1)
  self%dust%d => patch%mem(:,:,:,patch%idx%d              ,patch%it,1)
  self%dust%p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  !-----------------------------------------------------------------------------
  ! Allocate array for number of particles per cell-bin
  !-----------------------------------------------------------------------------
  m = patch%gn
  allocate (vg(3,m(1),m(2),m(3)))
  call io%bits_mem (storage_size(vg), &
                   product(shape(vg)), 'vg')
  !-----------------------------------------------------------------------------
  ! Cache the gas velocity
  !-----------------------------------------------------------------------------
  m = patch%gn
  l = self%lb
  u = self%ub
  do dir=1,3
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    vg(dir,i1,i2,i3) = patch%mem(i1,i2,i3,patch%idx%px+dir-1,patch%it,1) &
                     /(patch%mem(i1,i2,i3,patch%idx%d       ,patch%it,1) &
                      *patch%ds(dir))
  end do
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Allocate00 a new exports list, while the previous one may still be in use
  !-----------------------------------------------------------------------------
  wc = wallclock()
  allocate (exports(n_bin))
  do bin=1,n_bin
    paid = [patch%id,0,0,0,int(bin)]
    call exports(bin)%init (paid)
  end do
  !-----------------------------------------------------------------------------
  ! Import from other patches
  !-----------------------------------------------------------------------------
  call self%import
  !-----------------------------------------------------------------------------
  ! Non-vectorized pusher
  !-----------------------------------------------------------------------------
  ip = 0
  nk = 0
  no = 0
  nm = 0
  ds = minval(patch%ds, mask=patch%mesh%n > 1)
  call self%dust%consistency ('before KDK loop')
  times = 0d0
  wc1 = wc
  n_min = particles_per_bin*min_fraction
  n_max = particles_per_bin*max_fraction
  do bin=1,n_bin
    os = drag_coeff/self%s_cgs(bin)
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      ii = [i1,i2,i3]
      if (verbose > 0) then
        wc = wallclock()
      end if
      particles => self%particles(i1,i2,i3,bin)
      associate (n_particles => timer%n_particles)
      !$omp atomic update
      n_particles = n_particles + particles%n
      end associate
      paid = particles%id
      ip = ip + particles%n
      d = patch%mem(i1,i2,i3,patch%idx%d,patch%it,1)
      up = patch%u_phase(i1,i2,i3)/ds
      v(1) = vg(1,i1,i2,i3)
      v(2) = vg(2,i1,i2,i3)
      v(3) = vg(3,i1,i2,i3)
      pc = d*os
      p => particles%head%next
      j = 1
      np = particles%n
      do j=1,np
        next => p%next
        dt = 0.5*(self_time - p%t)
        if (dt > 0.0) then
          !--- kick 1 ------------------------------------------------------------
          pd = pc*(v-p%v)                                   ! drag force
          pa = 0.0                                          ! external force
          p%v = p%v + dt*(pa+pd)/(1.+pc*dt)                 ! stiff kick
          !--- drift -------------------------------------------------------------
          p%r = p%r + 2.*dt*p%v                             ! forward two dr
          p%t = self_time
          !if (self%forced) then
          !  !call self%force (p)
          !  r = mesh%p + (p%q + p%r - mesh%o)*mesh%d
          !  p%a = self%pforce (r)
          !end if
          !--- kick 2 ------------------------------------------------------------
          pd = pc*(v-p%v)                                   ! drag force
          p%v = p%v + dt*(pa+pd)/(1.+pc*dt)                 ! stiff kick
          !--- update q+r representation -----------------------------------------
          q = floor(p%r)                                    ! integer shift
          p%r = p%r - q                                     ! remove
          p%q(1:3) = p%q(1:3) + q                           ! add
          nk = nk+1                                         ! count
          if (particles%n < n_min) then
            call add_particles
          end if
          !---------------------------------------------------------------------
          move = .false.
          export = .false.
          do i=1,3
            if (p%q(i) /= ii(i)) then
              if (p%q(i) < l(i) .or. p%q(i) > u(i)) then
                export = .true.
                exit
              end if
              move = .true.
            end if
          end do
          !---------------------------------------------------------------------
          if (export) then
            call particles%remove (p)
            call exports(bin)%append (p)
            self%dd(i1,i2,i3,bin) = &
            self%dd(i1,i2,i3,bin) - p%w
            no = no+1
          !---------------------------------------------------------------------
          else if (move) then
            call particles%remove (p)
            call self%particles(p%q(1),p%q(2),p%q(3),p%q(4))%append (p)
            self%dd(i1,i2,i3,bin) = &
            self%dd(i1,i2,i3,bin) - p%w
            self%dd(p%q(1),p%q(2),p%q(3),bin) = &
            self%dd(p%q(1),p%q(2),p%q(3),bin) + p%w
            nm = nm+1
          end if
        end if
        p => next
      end do
      if (verbose > 0) call particles%validate ('after update')
      if (particles%n > n_max) then
        call remove (self, particles, self%dd(i1,i2,i3,bin))
      end if
    end do
    end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Briefly lock self%xports and copy in particles from the scratch exports
  !-----------------------------------------------------------------------------
  call self%lock%set ('copy')
  do bin=1,n_bin
    exportp => self%exports(bin)
    call exports(bin)%transfer (exportp)
  end do
  call self%lock%unset ('copy')
  deallocate (exports)
  !-----------------------------------------------------------------------------
  wc = wallclock()-wc
  associate (t_particles => timer%t_particles)
  !$omp atomic update
  t_particles = t_particles + wc
  end associate
  !-----------------------------------------------------------------------------
  np = self%np
  if (verbose > 0) then
    write (stdlog,'(2i12,f6.2,6p,2f6.2," mus/p",0p,2f6.2)') &
    np, nk, wc, wc/np, wc/nk, real(nm)/nk, real(no)/nk
    call self%dust%consistency ('after KDK + export')
    if (verbose > 1) then
      time(2) = wallclock()-time(1)
      write (stdlog, &
        '(a,i6,2x,2(a,1p,2e12.4,2x),a,2i8,2i6,2x,a,6p,f8.3,2x,a,3p,f6.1)') &
        'particle_solver_t%update: id:', patch%id, &
        't:', patch%time, self%time, &
        'dt:', patch%dtime, self%dtime, &
        'updates:', ip, nk, self%ni, self%no, &
        'mus/particle:', time(2)/max(nk,1), 'ms:', time(2)
    end if  
  end if
  !-----------------------------------------------------------------------------
  self%no = no
  self%dtime = patch%dtime
  if (nk > 0) then
    if (self%patch%id == i_dbg(1)) then
      write (stdlog,'(a,2(1p,4e10.3,2x))') &
        'particle_solver_t%update dtime: t, dt, t+dt =', &
        patch%time, patch%dtime, 1./patch%dtime, patch%time+patch%dtime, &
        self%time, self%dtime, 1./self%dtime, self%time+self%dtime
      write (stdlog,*) paid
    end if
    !---------------------------------------------------------------------------
  end if
  self%time = self_time 
  !-----------------------------------------------------------------------------
  ! Do dust evolution after move, since we need the particle timestep values
  !-----------------------------------------------------------------------------
  call self%dust%consistency ('before dust%update')
  if (self%dust%on) then
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      bins => self%particles(i1,i2,i3,:)
      call self%dust%update (bins)
    end do
    end do
    end do
  end if
  call self%dust%consistency ('after dust%update')
  !-----------------------------------------------------------------------------
  call io%bits_mem (storage_size(vg), -product(shape(vg)), 'vg')
  deallocate (vg)
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
contains

!===============================================================================
!> top up particles per bin
!===============================================================================
subroutine add_particles
  real:: d_cgs, w
  real(8):: time
  integer:: i, j, n
  !.............................................................................
  q = [i1,i2,i3]
  w = self%dd(i1,i2,i3,bin)/particles_per_bin
  n = particles%n
  p => particles%head%next
  time = 0d0
  do j=1,n
    p%w = w
    time = time + p%t
    p => p%next
  end do
  time = time/n
  do j=1,particles_per_bin-n
    p => pool%pop ()
    do i=1,3
      p%q(i) = q(i)
      p%v(i) = v(i)
      p%r(i) = self%random%ran1()
    end do
    p%q(4) = bin
    p%w = w
    p%t = time
    call particles%append (p)
  end do
end subroutine add_particles
END SUBROUTINE update

!===============================================================================
!> Bring particles per bin below max
!===============================================================================
SUBROUTINE remove (self, particles, dd)
  class(particle_solver_t)    :: self
  type(particles_t), pointer  :: particles
  class(particle_t), pointer  :: p, next
  real                        :: dd
  real                        :: w
  integer                     :: np, nf, j
  logical                     :: kill
  !.............................................................................
  nf = 0
  p => particles%head%next
  np = particles%n - particles_per_bin
  do j=1,np
    next => p%next
    select type (p)
    type is (particle_t)
      kill = .true.
    type is (vparticle_t)
      kill = .false.
      nf = nf+1
    end select
    if (kill) then
      call particles%remove (p)
      call pool%append (p)
    end if
    p => next
  end do
  if (verbose > 0 .and. nf > 0) then
    write (stdlog,*) 'remove: nf=', nf
  end if
  np = particles%n
  w = dd/np
  p => particles%head%next
  do j=1,np
    p%w = w
    p => p%next
  end do
END SUBROUTINE remove

!===============================================================================
!> Template
!===============================================================================
FUNCTION pforce (r) RESULT (f)
  real(8):: r(3), f(3)
  f = 0d0
END FUNCTION pforce

!===============================================================================
!> Import particles from nbors export lists
!===============================================================================
SUBROUTINE import (self)
  class(particle_solver_t)   :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: nbtask
  class(particles_t), pointer:: particles
  class(particle_t), pointer :: p, next
  real(8)                    :: r(3), dist(3)
  integer                    :: ni, ii(5), j, bin
  integer(kind=2)            :: q(3)
  logical                    :: inside
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  ! Loop over patch nbors
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%import',itimer=itimer)
  self%ni = 0
  nbor => self%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (particle_solver_t)
      call nbtask%lock%set ('import')
      dist = nbtask%patch%task_t%distance(self%patch)/nbtask%patch%ds
      !-------------------------------------------------------------------------
      ! Loop over particles in export list and grab those that belong here
      !-------------------------------------------------------------------------
      ni = 0
      do bin=1,n_bin
        p => nbtask%exports(bin)%head%next
        do j=1,nbtask%exports(bin)%n
          next => p%next
          r = dist + (p%r+p%q(1:3))
          q = floor(r)
          inside = all(q >= self%lb .and. q <= self%ub)
          if (inside) then
            p%r = r-q
            p%q(1:3) = q
            call nbtask%exports(bin)%remove (p)
            call self%dust%consistency2 ('before add',p%q(1),p%q(2),p%q(3),p%q(4))
            particles => self%particles(p%q(1),p%q(2),p%q(3),p%q(4))
            call particles%append (p)
            call self%dust%consistency2 ('after add',p%q(1),p%q(2),p%q(3),p%q(4))
            self%dd(q(1),q(2),q(3),p%q(4)) = &
            self%dd(q(1),q(2),q(3),p%q(4)) + p%w
            ni = ni + 1
            ii(1) = self%patch%id; ii(2:4) = p%q(1:3); ii(5) = p%q(4)
          end if
          p => next
        end do
      end do
      self%ni = self%ni + ni
      if ((ni > 0 .and. verbose > 1) .or. verbose > 3) then
        write (stdlog,'(a,4i6,2x,3i4)') 'imported: id, nbor, ni, offset =', &
          self%patch%id, nbtask%patch%id, nbtask%exports(bin)%n, ni, &
          nint(dist)
      end if
      call nbtask%lock%unset ('import')
    end select
    nbor => nbor%next
  end do
  if ((self%ni > 0 .and. verbose > 1) .or. verbose > 2) then
    write (stdlog,'(a,2i6)') 'imported: id, self%ni =', &
      self%patch%id, self%ni
    flush (stdlog)
  end if
  call trace%end (itimer)
END SUBROUTINE import

!===============================================================================
!> Output dust densities and (optionally) dust particles
!===============================================================================
SUBROUTINE output (self)
  class(particle_solver_t)      :: self
  !.............................................................................
  class(particles_t), pointer   :: particles
  class(mesh_t), pointer        :: mesh(:)
  class(particle_t), pointer    :: p
  real, pointer, contiguous     :: tmp_r(:,:)
  integer, pointer, contiguous  :: tmp_i(:)
  integer                       :: i, j, n, i1, i2, i3, bin
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particles_t%output')
  if (verbose > 1) then
    write (stdlog, '(a,i6,2f8.4)') &
      'particles_t%output: id, time, aux_next =', &
      self%patch%id, self%time, self%aux%aux_next
  end if
  if (self%time >= self%aux%aux_next) then
    call self%aux%register ('dust_s', self%s_cgs)
    call self%aux%register ('dust_d', self%dd)
    !...........................................................................
    if (verbose > 1) then
      do i=1,n_bin
        write (stdlog, '(a,i6, i3,2x,1p,10e11.2)') &
          'particles_t%output: id, bin, n =', self%patch%id, i, self%dd(3,3,1:10,i)
      end do
    end if
    !...........................................................................
    if (particle_io) then
      mesh => self%patch%mesh
      n = 0
      p => particles%head%next
      do j=1,particles%n
        select type (p)
        type is (vparticle_t)
          n = n+1
        end select
        p => p%next
      end do
      allocate (tmp_r(5,n))
      allocate (tmp_i(  n))
      i = 0
      do i3=self%lb(3),self%ub(3)
      do i2=self%lb(2),self%ub(2)
      do i1=self%lb(1),self%ub(1)
        do bin=1,n_bin
          particles => self%particles(i1,i2,i3,bin)
          p => particles%head%next
          do j=1,particles%n
            select type (p)
            type is (vparticle_t)
              i = i+1
              tmp_r(1:3,i) = mesh%p + (p%q(1:3) + p%r - mesh%o)*mesh%d
              tmp_r(  4,i) = p%w
              tmp_r(  5,i) = self%s_cgs(p%q(4))
              tmp_i(    i) = p%id
            end select
            p => p%next
          end do
        end do
        end do
      end do
      end do
      call self%aux%register ('dust_p' , tmp_r)
      call self%aux%register ('dust_id', tmp_i)
      call self%aux%output_sequence (self%patch%id, self%time)
      if (verbose > 0) then
        write (stdlog, '(a,i6,f8.4,i3,2i8)') &
         'particles_t%output: id, particle array =', &
          self%patch%id, self%time, shape(tmp_r), i
      end if
      deallocate (tmp_r, tmp_i)
    else
      call self%aux%output_sequence (self%patch%id, self%time)
    end if
    if (verbose > 1) then
      write (stdlog, '(a,i6,1p,10e11.2)') &
       'particles_t%output: self%aux%aux_next =', &
        self%patch%id, self%aux%aux_next
    end if
  end if
  call trace%end ()
END SUBROUTINE output

!===============================================================================
SUBROUTINE statistics (self)
  class(particle_solver_t):: self
  !-----------------------------------------------------------------------------
  call self%particles(1,1,1,1)%statistics
END SUBROUTINE statistics

END MODULE particle_solver_mod
