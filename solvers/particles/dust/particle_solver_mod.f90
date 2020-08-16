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
  USE mpi_mod
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
  USE scratch_mod
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
    real, pointer, contiguous  :: dg(:,:,:)            ! gas density
    real, pointer, contiguous  :: vg(:,:,:,:)          ! gas velocity
    real, pointer, contiguous  :: dd(:,:,:,:) => null()! number density per bin
    real, pointer, contiguous  :: s_cgs(:)             ! particle size
    real, pointer, contiguous  :: m_cgs(:)             ! particle size
    logical                    :: forced = .false.     ! external force
    logical, dimension(3)      :: lower                ! has lower boundary
    logical, dimension(3)      :: upper                ! has upper boundary
    integer, dimension(3)      :: lb=0, ub=0           ! lower and upper index bounds
    integer                    :: np = 0               ! particles
    integer                    :: ni = 0               ! particles imported
    integer                    :: ne = 0               ! particles outside
    integer                    :: n_bin                ! number of bins
    real                       :: dust_to_gas          ! dust:gas ratio
    real                       :: vmax                 ! max particle speed
    real                       :: v_min, v_max         ! min/max u_phase
    class(list_t), pointer     :: task_list            ! task list
    real(8)                    :: new_momentum(3,2)=0. ! new patch momentum
    real(8)                    :: old_momentum(3,2)=0. ! old patch momentum
    procedure(pforce), pointer, nopass :: pforce => null() ! point force
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: initial_condition
    procedure:: init_task_list
    procedure:: update
    procedure:: import
    procedure:: output
    procedure:: statistics
    procedure:: update_momentum
    procedure:: global_momentum
    procedure:: courant_condition
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
  real(8):: out_time = -1.0                   ! output cadence
  logical:: particle_io = .false.             ! output particle positions?
  logical:: backreaction = .false.            ! backreaction on gas
  logical:: detailed_timer = .false.
  logical:: on = .true.
  !-----------------------------------------------------------------------------
  ! Total momentum is counted using the mpi_counter data type.  The
  ! mpi_counter_t%update (increment) returns the value accumulated over ranks
  !-----------------------------------------------------------------------------
  type(mpi_counter_t), target:: momentum_counter(3,2) ! particles and gas
  real(8), target:: total_momentum(3,2)               ! their actual values
CONTAINS

!===============================================================================
!> Initialize particle solver, reading parameters and adding task link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(particle_solver_t), target  :: self
  class(gpatch_t), pointer          :: patch
  integer                           :: bin, iostat
  real                              :: solid_density, dln_m
  namelist /particle_solver_params/ &
    on, verbose, detailed_timer, particle_io, out_time, &
    n_bin, particles_per_bin, seed, ip_dbg, id_dbg, i_dbg, &
    min_fraction, max_fraction, max_size_cgs, min_size_cgs, &
    density_yr_cgs, solid_density_cgs, dust_to_gas, output_every, &
    backreaction
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pre_init')
  self%patch => patch                           ! make the patch known
  self%id = patch%id                            ! for now, inherit
  self%dust%patch => self%patch                 ! make the patch known
  !-----------------------------------------------------------------------------
  call io%lock%set ('particle_solver_t%pre_init')
  rewind (io%input)
  read (io%input, particle_solver_params, iostat=iostat)
  if (first_time) then
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%namelist, particle_solver_params)
    first_time = .false.
  end if
  call io%lock%unset ('particle_solver_t%pre_init')
  !-----------------------------------------------------------------------------
  self%n_bin = n_bin
  self%dust_to_gas = dust_to_gas
  self%out_next = 0d0
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
  self%iout = 0
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE init (self)
  class(particle_solver_t) :: self
  class(gpatch_t), pointer :: pa
  integer                  :: m(4), i, j
  !-----------------------------------------------------------------------------
  ! Initialize particle lists and dust density array
  !-----------------------------------------------------------------------------
  ! if (.not.on) return
  call trace%begin ('particle_solver_t%init')
  m(1:3) = self%patch%mesh%gn
  allocate (self%particles(m(1),m(2),m(3),n_bin)) 
  call io%bits_mem (storage_size(self%particles), &
                   product(shape(self%particles)), 'particles')
  allocate (pool, self%exports(n_bin))
  call io%bits_mem (storage_size(pool), omp%nthreads+n_bin, 'exports')
  call pool%init ([self%patch%id,0,0,0,0])
  do i=1,n_bin
    call self%exports(i)%init ([self%patch%id,0,0,0,i])
  end do
  !-----------------------------------------------------------------------------
  allocate (self%dd(m(1),m(2),m(3),n_bin))
  call io%bits_mem (storage_size(self%dd), product(shape(self%dd)), 'dd')
  !-----------------------------------------------------------------------------
  ! Allocate force array
  !-----------------------------------------------------------------------------
  m(1:3) = self%patch%mesh%gn
  m(4) = 3
  call scratch%allocate (scratch%force_per_unit_volume, m)
  if (.not.associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init ('part')
  end if
  call self%dust%init (n_bin, self%m_cgs, self%dd)
  call self%aux%init
  call self%random%init (seed)        ! FIXME: check seed inheritance
  do i=1,3
  do j=1,2
    call momentum_counter(i,j)%init (0d0)
  end do
  end do
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
  integer                     :: m(3), i1, i2, i3, i4, i, j, n, nv, l, nb, id(5)
  real                        :: d, w, vg(3), chance
  logical                     :: draw
  integer, save               :: itimer=0, extra_bits
  integer, save               :: lid = 0
  !$omp threadprivate (lid)
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
      if (patch%position(i) - patch%size(i)/2.0 - patch%ds(i) &
          < -patch%box(i)/2.0) then
        self%lower(i) = .true.
        self%lb(i) = patch%mesh(i)%lb
      end if
      if (patch%position(i) + patch%size(i)/2.0 + patch%ds(i) &
          > +patch%box(i)/2.0) then
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
    nb = (storage_size(p)+extra_bits)/8
    print '(a,i3,a,f18.3,a)', ' estimated memory for particles (', &
      nb,' bytes):',nb*n_bin*particles_per_bin/1024.**3 &
      *product(self%patch%mesh%b/self%patch%mesh%d), ' GB'
    call pool%append (p)
  end if
  !-----------------------------------------------------------------------------
  ! Initialize dust.  At boundaries, add particles also in guard zones
  !-----------------------------------------------------------------------------
#ifdef __GFORTRAN__
  !$omp critical (alloc_cr)
#endif
  if (output_every > 0) then
    chance = 1.0/output_every
  else
    chance = 0.0
  end if
  l = 0
  self%np = 0
  do i3=self%lb(3),self%ub(3)
  do i2=self%lb(2),self%ub(2)
    n = 0
    do i1=self%lb(1),self%ub(1)
      d = patch%mem(i1,i2,i3,patch%idx%d,patch%it,1)
      w = d*dust_to_gas/(n_bin*particles_per_bin)
      vg = patch%mem(i1,i2,i3,patch%idx%px:patch%idx%pz,patch%it,1) &
        / (patch%ds*d)
      do i4=1,n_bin
        particles => self%particles(i1,i2,i3,i4)
        id = [self%patch%id,i1,i2,i3,i4]
        call particles%init (id)
        n = n+2
        nv = 0
        do j=1,particles_per_bin
          n = n+1
          if (output_every > 0) then
            draw = self%random%ran1() < chance
          else
            draw = .false.
          end if
          if (draw) then  ! every output_every is for I/O
            allocate(pv)
            call pv%init
            pv%id = 1 + omp%thread_global + (lid+nv)*omp%nthreads_global
            nv = nv+1
            pv%dg = d
            pv%vg = vg
            p => pv
          else
            allocate(pn)
            call pn%init
            p => pn
          end if
          p%q = [i1,i2,i3,i4]
          p%r = self%random%uniform (3)
          p%v = 0.0
          p%w = w
          p%t = 0
          call particles%append (p)
        end do
        self%np = self%np + particles%n
      end do
      l = l+particles_per_bin
      lid = lid+nv
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
  real, dimension(3)          :: r, dr
  class(particle_t), pointer  :: p, next
  real(8)                     :: time(2), v(3), wc, start
  real                        :: ds, dt, os, n_min, n_max, ods, dd
  integer                     :: n, np, ne, nk, nm, nu, nx, m(3), ii(3), j
  integer(kind=2)             :: i, i1, i2, i3, i4, l(3), u(3), q(3)
  integer                     :: id, paid(5), istep, n_add, n_rem, np_prev
  integer                     :: dim(4)
  integer(8)                  :: pid
  logical                     :: move, export
  class(gpatch_t), pointer    :: patch
  type(particles_t), pointer  :: exports(:), exportp
  type(particles_t), pointer  :: particles, bins(:)
  real                        :: d, up, pc, pd(3), pa(3), w, pv(3), dv(3), pwd
  real(8)                     :: times(9), wc1, wc2, sum
  real(8)                     :: momentum(3,2), dvol
  real, pointer, contiguous   :: p1(:,:,:), p2(:,:,:), p3(:,:,:)
  type(mpi_counter_t), pointer:: counter
  integer, save               :: itimer=0
!-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%update', itimer=itimer)
  start = wallclock()
  patch => self%patch
  p1 => patch%mem(:,:,:,patch%idx%px,patch%new,1)
  p2 => patch%mem(:,:,:,patch%idx%py,patch%new,1)
  p3 => patch%mem(:,:,:,patch%idx%pz,patch%new,1)
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
  !write (stdlog,*) 'dtime(1)', self%dtime
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
  allocate (self%vg(3,m(1),m(2),m(3)))
  call io%bits_mem (storage_size(self%vg), &
                   product(shape(self%vg)), 'vg')
  !-----------------------------------------------------------------------------
  ! Cache the gas velocity
  !-----------------------------------------------------------------------------
  self%dg => patch%mem(:,:,:,patch%idx%d,patch%new,1)
  m = patch%gn
  l = self%lb
  u = self%ub
  do i=1,3
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    self%vg(i,i1,i2,i3) = patch%mem(i1,i2,i3,patch%idx%px+i-1,patch%it,1) &
                        / (self%dg(i1,i2,i3)*patch%ds(i))
    if (abs(self%vg(i,i1,i2,i3)) > 1e1/patch%ds(1)) then
      write(stderr,'(a,2i6,4i4,2(2x,1p,3e11.2))') &
      'warning, large vg:', patch%id,patch%istep,i1,i2,i3,i, &
      patch%mem(i1,i2,i3,patch%idx%px+i-1,patch%it,1), &
      self%dg(i1,i2,i3)
    end if
  end do
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Allocate00 a new exports list, while the previous one may still be in use
  !-----------------------------------------------------------------------------
  wc = wallclock()
  allocate (exports(n_bin))
  do i4=1,n_bin
    paid = [patch%id,0,0,0,int(i4)]
    call exports(i4)%init (paid)
  end do
  !-----------------------------------------------------------------------------
  ! Non-vectorized pusher
  !-----------------------------------------------------------------------------
  ds = minval(patch%ds, mask=patch%mesh%n > 1)
  times = 0d0
  wc1 = wc
  n_min = particles_per_bin*min_fraction
  n_max = particles_per_bin*max_fraction
  !-----------------------------------------------------------------------------
  ! Import from other patches.
  !-----------------------------------------------------------------------------
  call self%import
  !-----------------------------------------------------------------------------
  ! After the import, the interior particles should be maximally smoothly
  ! distributed, which is thus a good point to measure the dust bin density,
  ! renormalize, and do output.
  !-----------------------------------------------------------------------------
  n_add = 0
  n_rem = 0
  np_prev = self%np
  call check_count (self, self%np, 'before add/rem')
  do i4=1,n_bin
    do i3=l(3)+1,u(3)-1
    do i2=l(2)+1,u(2)-1
    do i1=l(1)+1,u(1)-1
      particles => self%particles(i1,i2,i3,i4)
      dd = 0.0
      p => particles%head%next
      do j=1,particles%n
        dd = dd + p%w
        p => p%next
      end do
      self%dd(i1,i2,i3,i4) = dd
      np = particles%n
      if (particles%n < n_min) then
        n_add = n_add+1
        if (verbose > 3) &
          write(stdlog,'(2i6,3i4,2x,a,i4)') &
            patch%id,patch%istep,i1,i2,i3,'add',particles%n
        call add (self, particles)
      else if (particles%n > n_max) then
        n_rem = n_rem+1
        if (verbose > 3) &
          write(stdlog,'(2i6,3i4,2x,a,i4)') &
            patch%id,patch%istep,i1,i2,i3,'rem',particles%n
        call remove (self, particles)
      end if
      self%np = self%np + (particles%n - np)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Copy edges from the next interior layer -- this is better than relying on
    ! particles that may be missing delayed import.  The dust densities are ONLY
    ! used for dust size evolution!
    !---------------------------------------------------------------------------
    do i3=l(3),u(3)
    do i2=l(2),u(2)
      self%dd(l(1),i2,i3,i4) = self%dd(l(1)+1,i2,i3,i4)
      self%dd(u(1),i2,i3,i4) = self%dd(u(1)-1,i2,i3,i4)
    end do
    end do
    do i3=l(3),u(3)
    do i1=l(1),u(1)
      self%dd(i1,l(2),i3,i4) = self%dd(i1,l(2)+1,i3,i4)
      self%dd(i1,u(2),i3,i4) = self%dd(i1,u(2)-1,i3,i4)
    end do
    end do
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      self%dd(i1,i2,l(3),i4) = self%dd(i1,i2,l(3)+1,i4)
      self%dd(i1,i2,u(3),i4) = self%dd(i1,i2,u(3)-1,i4)
    end do
    end do
  end do
  call self%output
  !-----------------------------------------------------------------------------
  ! Time step control
  !-----------------------------------------------------------------------------
  self%dtime = patch%dtime
  self%time = self%time + self%dtime
  !-----------------------------------------------------------------------------
  ! Now move
  !-----------------------------------------------------------------------------
  np_prev = self%np
  nu = 0
  nk = 0
  ne = 0
  nm = 0
  self%new_momentum = 0d0
  self%vmax = 0.0
  self%np = 0
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    v(1) = self%vg(1,i1,i2,i3)
    v(2) = self%vg(2,i1,i2,i3)
    v(3) = self%vg(3,i1,i2,i3)
    d = patch%mem(i1,i2,i3,patch%idx%d,patch%it,1)
    ii = [i1,i2,i3]
    momentum(:,1) = 0d0
    do i4=1,n_bin
      particles => self%particles(i1,i2,i3,i4)
      nu = nu + particles%n
      os = drag_coeff/self%s_cgs(i4)
      associate (n_particles => timer%n_particles)
      !$omp atomic update
      n_particles = n_particles + particles%n
      end associate
      paid = particles%id
      up = patch%u_phase(i1,i2,i3)/ds
      pc = d*os
      np = particles%n
      p => particles%head%next
      do j=1,np
        next => p%next
        dt = 0.5*(self%time - p%t)
        if (dt > 0.0) then
          if (abs(maxval(p%q)) > 100) then
            write (stdlog,'(i6,2x,a,2i4,2x,3i6,3(2x,1p,3e12.3))') &
            patch%id, 'move error0', j, particles%n, p%q(1:3), p%r, p%v, dt
          end if
          !- kick 1 ------------------------------------------------------------
          pa = 0.0                                     ! external force
          pd = pc*(v-p%v)                              ! drag force
          pwd = p%w/d                                  ! density ratio
          !---------------------------------------------------------------------
          ! To allow extreme dust:gas ratios, compute the largest action first
          !---------------------------------------------------------------------
          if (backreaction) then
            dv = dt*(pa+pd)/(1.+pc*dt*(1.+pwd))        ! particle velocity delta
            p%v = p%v + dv                             ! particle update
            v = v - dv*pwd                             ! gas update
          else
            p%v = p%v + dt*(pa+pd)/(1.+pc*dt)          ! particle update
          end if
          if (max(abs(p%v(1)),abs(p%v(2)),abs(p%v(3))) > 1e5) then
            write(stderr,*) 'warning, large p%v:', patch%id,i1,i2,i3,i4,j
          end if
          !- drift -------------------------------------------------------------
          p%r = p%r + 2.*dt*p%v                        ! forward two dr
          p%t = self%time
          !if (self%forced) then
          !  !call self%force (p)
          !  r = mesh%p + (p%q + p%r - mesh%o)*mesh%d
          !  p%a = self%pforce (r)
          !end if
          !- kick 2 ------------------------------------------------------------
          pd = pc*(v-p%v)                              ! drag force
          if (backreaction) then
            dv = dt*(pa+pd)/(1.+pc*dt*(1.+pwd))        ! particle velocity delta
            p%v = p%v + dv                             ! particle update
            v = v - dv*pwd                             ! gas update
          else
            p%v = p%v + dt*(pa+pd)/(1.+pc*dt)          ! particle update
          end if
          momentum(:,1) = momentum(:,1) + p%v*p%w      ! particle momentum
          !- update q+r representation -----------------------------------------
          if (abs(maxval(p%r)) > 100) then
            q = floor(p%r)
            write (stdlog,'(i6,2x,a,2i4,2(2x,3i6),3(2x,1p,3e12.3))') &
            patch%id, 'move error1', j, particles%n, p%q(1:3), q, p%r, p%v, dt
          end if
          q = floor(p%r)                               ! integer shift
          p%r = p%r - q                                ! remove
          p%q(1:3) = p%q(1:3) + q                      ! add
          if (abs(maxval(p%q)) > 100) then
            write (stdlog,'(i6,2x,a,2i4,2(2x,3i6),3(2x,1p,3e12.3))') &
            patch%id, 'move error2', j, particles%n, p%q(1:3), q, p%r, p%v, dt
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
          ! Export particle to another patch
          !---------------------------------------------------------------------
          if (export) then
            call particles%remove (p)
            call exports(i4)%append (p)
            if (abs(maxval(p%q)) > 100) &
              write (stdlog,*) 'export error', j, exports(i4)%n, p%q
            self%dd(i1,i2,i3,i4) = &
            self%dd(i1,i2,i3,i4) - p%w
            ne = ne+1
            if (verbose > 2) then
              write (stdlog,*) patch%id, 'export: q=', patch%id, p%q
            end if
          !---------------------------------------------------------------------
          ! Move particle to another cell. In case the cell has been visited,
          ! the dust bin density is updated as well.
          !---------------------------------------------------------------------
          else if (move) then
            self%vmax = max(self%vmax,abs(p%v(1)),abs(p%v(2)),abs(p%v(3)))
            call particles%remove (p)
            call self%particles(p%q(1),p%q(2),p%q(3),p%q(4))%append (p)
            self%dd(i1,i2,i3,i4) = &
            self%dd(i1,i2,i3,i4) - p%w
            self%dd(p%q(1),p%q(2),p%q(3),i4) = &
            self%dd(p%q(1),p%q(2),p%q(3),i4) + p%w
            nm = nm+1
          !---------------------------------------------------------------------
          ! Particle stays in this cell -- update vmax
          !---------------------------------------------------------------------
          else
            self%vmax = max(self%vmax,abs(p%v(1)),abs(p%v(2)),abs(p%v(3)))
          end if
        else
          nk = nk + 1
        end if
        p => next
      end do ! loop over particles
      self%np = self%np + particles%n
    end do ! loop over bins
    !---------------------------------------------------------------------------
    ! Instead of storing the new momentum directly, store it as a force that
    ! will end up updating the momentum with exactly the same amount, but will
    ! also be used in the solver predictor step.
    !---------------------------------------------------------------------------
    if (backreaction) then
      block
        real:: force(3)
        force(1) = (d*v(1)*patch%ds(1) - p1(i1,i2,i3))/self%dtime
        force(2) = (d*v(2)*patch%ds(2) - p2(i1,i2,i3))/self%dtime
        force(3) = (d*v(3)*patch%ds(3) - p3(i1,i2,i3))/self%dtime
        scratch%force_per_unit_volume(i1,i2,i3,:) = force
      end block
    end if
    !---------------------------------------------------------------------------
    ! Maintain up-to-date momentum measures of gas and particles separately,
    ! based on increments, for consistency checks
    !---------------------------------------------------------------------------
    dvol = product(patch%ds)
    momentum(:,1) = momentum(:,1)*patch%ds          ! particle momentum
    momentum(:,2) = d*v*patch%ds                    ! gas momentum
    self%new_momentum(:,1) = &
    self%new_momentum(:,1) + momentum(:,1)*dvol     ! measured gas momentum
    self%new_momentum(:,2) = &
    self%new_momentum(:,2) + momentum(:,2)*dvol     ! measured gas momentum
  end do
  end do
  end do
  self%ne = ne
  if (self%np > np_prev) then
    write(stdlog,*) patch%id, self%np,' should be <=', np_prev
    call check_count (self, self%np, 'after move')
  end if
  !-----------------------------------------------------------------------------
  ! Update the total momentum, maintained by one counter per direction and
  ! category.  The counter returns the total after the increment, and should
  ! be given the new local contribution minus the previous local contribution
  !
  ! total_momentum: momentum of entire experiment
  !   new_momentum: new momentum of this patch
  !   old_momentum: old momentum of this patch  
  !-----------------------------------------------------------------------------
  call self%update_momentum
  !-----------------------------------------------------------------------------
  ! Briefly lock self%xports and copy in particles from the scratch exports
  !-----------------------------------------------------------------------------
  call self%lock%set ('copy')
  nx = 0
  do i4=1,n_bin
    exportp => self%exports(i4)
    nx = nx + exportp%n
    call exports(i4)%transfer (exportp)
    call exports(i4)%dealloc
  end do
  call self%lock%unset ('copy')
  deallocate (exports)
  if (verbose > 1) then
    if (ne > 0) then
      write(stdlog,'(i6,2(2x,a,i8))') &
        self%patch%id, 'added', ne, 'to exports with', nx 
    end if
  end if
  !-----------------------------------------------------------------------------
  wc = wallclock()-wc
  associate (t_particles => timer%t_particles)
  !$omp atomic update
  t_particles = t_particles + wc
  end associate
  !-----------------------------------------------------------------------------
  if (verbose > 0) then
    np = self%np
    if (verbose > 1) then
      time(2) = wallclock()-time(1)
      write (stdlog, &
        '(a,i6,2x,2(a,1p,1e10.3,2x),a,6i8,2i5,3x,a,6p,f7.3,3x,a,3p,f6.1)') &
        'particle_solver_t%update id:', patch%id, &
        't:', self%time, &
        'dt:', self%dtime, &
        'n:', nu, np, self%ni, ne, nx, nk, n_add, n_rem, &
        'mus/p:', time(2)/nu, 'ms:', time(2)
    else
      write (stdlog,'(i6,2i12,f6.3,"s",6p,2f6.3," mus/p",0p,2f6.2)') &
      patch%id, np, nk, wc, wc/np, wc/nk, real(nm)/nk, real(ne)/nk
    end if
  end if
  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------
  ! Do dust evolution after move, since we need the particle timestep values
  !-----------------------------------------------------------------------------
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
  !-----------------------------------------------------------------------------
  call io%bits_mem (storage_size(self%vg), -product(shape(self%vg)), 'vg')
  deallocate (self%vg)
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> top up particles per size bin
!===============================================================================
SUBROUTINE add (self, particles)
  class(particle_solver_t)    :: self
  type(particles_t), pointer  :: particles
  class(particle_t), pointer  :: p, next
  real                        :: w, m0(3), m1(3), v(3), dd
  real(8)                     :: t
  integer                     :: n, i, j, q(4), n_add
  !.............................................................................
  m0 = 0.0                    ! initial momentum
  dd = 0.0                    ! dust bin density
  n = particles%n
  p => particles%head%next
  do j=1,n
    dd = dd + p%w
    m0 = m0 + p%v*p%w
    p => p%next
  end do
  t = 0d0                      ! time
  m1 = 0.0                     ! momentum after renormalization
  w = dd/particles_per_bin     ! new weights
  p => particles%head%next
  q = p%q
  do j=1,n
    p%w = w
    m1 = m1 + p%v*p%w
    t = t + p%t
    p => p%next
  end do
  t = t/n                       ! average time
  n_add = particles_per_bin-n   ! number to add
  if (w==0.0 .or. n_add==0) then
    v = 0.0                     ! rare case, or many zero-w particles
  else
    v = (m0-m1)/(n_add*w)       ! speed to match momentum
  end if
  do j=1,n_add
    p => pool%pop ()
    do i=1,3
      p%q(i) = q(i)
      p%v(i) = v(i)
      p%r(i) = self%random%ran1()
    end do
    p%q(4) = q(4)
    p%w = w
    p%t = t
    call particles%append (p)
  end do
  if (.false.) then
    m1 = 0.0
    p => particles%head%next
    q = p%q
    do j=1,particles%n
      m1 = m1 + p%v*p%w
      p => p%next
    end do
    print '(4i4,2(2x,1p,3e14.6))',q(1:3),n_add,m0,m1
  end if
END SUBROUTINE add

!===============================================================================
!> Bring particles per bin below max
!===============================================================================
SUBROUTINE remove (self, particles)
  class(particle_solver_t)    :: self
  type(particles_t), pointer  :: particles
  class(particle_t), pointer  :: p, next
  real                        :: dd
  real                        :: w, m0(3), m1(3), m2(3), v(3)
  integer                     :: np, nf, j, q(4)
  logical                     :: kill
  logical                     :: debug=.false.
  !-----------------------------------------------------------------------------
  ! Measure momentum
  !-----------------------------------------------------------------------------
  m0 = 0.0
  dd = 0.0
  np = particles%n
  p => particles%head%next
  do j=1,np
    m0 = m0 + p%v*p%w                   ! total momentum
    dd = dd + p%w                       ! dust bin density
    p => p%next
  end do
  !-----------------------------------------------------------------------------
  ! Remove particles
  !-----------------------------------------------------------------------------
  p => particles%head%next
  do j=1,np
    next => p%next
    if (particles%n > particles_per_bin) then
      kill = .false.
      select type (p)
      type is (particle_t)
        kill = .true.
      end select
      if (kill) then
        call particles%remove (p)
        call pool%append (p)
      end if
    else
      exit
    end if
    p => next
  end do
  !-----------------------------------------------------------------------------
  ! Measure momentum with renormalized weights
  !-----------------------------------------------------------------------------
  np = particles%n
  w = dd/np                               ! new weight factor
  m1 = 0.0
  p => particles%head%next
  do j=1,np
    p%w = w
    m1 = m1 + p%v*p%w
    p => p%next
  end do
  !-----------------------------------------------------------------------------
  ! Adjust momentum
  !-----------------------------------------------------------------------------
  if (w==0.0) then
    v = 0.0                       ! rare case, or many particle-w = 0.0
  else
    v = (m0-m1)/(w*np)            ! velocity increment to compensate
  end if
  p => particles%head%next
  do j=1,np
    p%v = p%v + v
    p => p%next
  end do
  if (debug) then
    m2 = 0.0
    p => particles%head%next
    q = p%q
    do j=1,particles%n
      m2 = m2 + p%v*p%w
      p => p%next
    end do
    print '(4i4,5(2x,1p,3f8.4))', q(1:3),np,m0,m1,m2
  end if
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
  integer                    :: ni, ii(5), j, bin, n(4)
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
      n = 0
      do bin=1,n_bin
        p => nbtask%exports(bin)%head%next
        n(1) = nbtask%exports(bin)%n
        n(3) = self%np
        do j=1,nbtask%exports(bin)%n
          if (abs(maxval(p%q)) > 100) write (stdlog,*) 'import error', j, p%q
          next => p%next
          r = dist + (p%r+p%q(1:3))
          q = floor(r)
          inside = all(q >= self%lb .and. q <= self%ub)
          if (inside) then
            p%r = r-q
            p%q(1:3) = q
            if (verbose > 2) then
              write(stdlog,*) 'import: q=', self%patch%id, p%q
            end if
            call nbtask%exports(bin)%remove (p)
            particles => self%particles(p%q(1),p%q(2),p%q(3),p%q(4))
            call particles%append (p)
            self%dd(q(1),q(2),q(3),p%q(4)) = &
            self%dd(q(1),q(2),q(3),p%q(4)) + p%w
            ni = ni + 1
            ii(1) = self%patch%id; ii(2:4) = p%q(1:3); ii(5) = p%q(4)
          else if (verbose > 2) then
            write (stdlog,'(2i6,3f10.3,2x,a,4(2x,3i6))') &
              self%patch%id,nbtask%id,dist,'rejected',p%q,q,self%lb,self%ub        
          end if
          p => next
        end do
        self%np = self%np + ni
        n(2) = nbtask%exports(bin)%n
        n(4) = self%np
        if (verbose > 1) then
          if (ni > 0) then
            write(stdlog,'(2(i6,2x,a),i6,4i8)') &
              self%patch%id, 'picked up', ni, 'from', nbtask%id, n
          end if
        end if
        if ((ni > 0 .and. verbose > 1) .or. verbose > 3) then
          write (stdlog,'(a,2i6,2i10,2x,3i4)') 'imported: id, nbor, ni, offset =', &
            self%patch%id, nbtask%patch%id, nbtask%exports(bin)%n, ni, &
            nint(dist)
        end if
      end do
      self%ni = self%ni + ni
      call nbtask%lock%unset ('import')
    end select
    nbor => nbor%next
  end do
  self%np = self%np + self%ni
  if ((self%ni > 0 .and. verbose > 1) .or. verbose > 3  ) then
    write (stdlog,'(a,i6,i10)') 'imported: id, self%ni =', &
      self%patch%id, self%ni
    flush (stdlog)
  end if
  call trace%end (itimer)
END SUBROUTINE import

!===============================================================================
!> Check that n is the correct number of solver particles
!===============================================================================
SUBROUTINE check_count (self, n, label)
  class(particle_solver_t)      :: self
  integer                       :: n
  character(len=*)              :: label
  !.............................................................................
  class(particles_t), pointer   :: particles
  class(particle_t), pointer    :: p
  integer                       :: np, i1, i2, i3, i4, j
  !-----------------------------------------------------------------------------
  if (verbose > 3) then
    np = 0
    do i4=1,n_bin
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      particles => self%particles(i1,i2,i3,i4)
      np = np + particles%n
    end do
    end do
    end do
    end do
    if (np /= n) then
      write (stdlog,*) self%id, n, 'should be', np, '  ', label
    end if
  end if
END SUBROUTINE check_count

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
  real(8), pointer              :: total_mom(:,:)
  integer, pointer, contiguous  :: tmp_i(:)
  integer                       :: i, j, n, i1, i2, i3, i4
  integer, save                 :: itimer=0
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%output', itimer=itimer)
  if (verbose > 1) then
    write (stdlog, '(a,i6,2f8.4)') &
      'particles_t%output: id, time, aux_next =', &
      self%patch%id, self%time, self%aux%aux_next
  end if
  !-----------------------------------------------------------------------------
  ! High-cadence output of v_particle_t particles
  !-----------------------------------------------------------------------------
  if (out_time > 0d0 .and. self%time >= self%out_next) then
    call voutput (self)
    self%iout = self%iout + 1
    self%out_next = self%out_next + out_time
  end if
  if (self%time >= self%aux%aux_next) then
    call self%aux%register ('dust_s', self%s_cgs)
    call self%aux%register ('dust_d', self%dd)
    total_mom => total_momentum
    call self%aux%register ('total_momentum', total_mom)
    !...........................................................................
    if (verbose > 1) then
      do i=1,n_bin
        write (stdlog, '(a,i6, i3,2x,1p,10e11.2)') &
       'particles_t%output: id, bin, n =', self%patch%id, i, self%dd(3,3,1:10,i)
      end do
    end if
    call self%aux%output_sequence (self%patch%id, self%time)
    if (verbose > 1) then
      write (stdlog, '(a,i6,1p,10e11.2)') &
       'particles_t%output: self%aux%aux_next =', &
        self%patch%id, self%aux%aux_next
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE output

!===============================================================================
!> Initialize module
!===============================================================================
SUBROUTINE voutput (self)
  class(particle_solver_t):: self
  class(particles_t), pointer:: pp
  class(particle_t), pointer:: p
  real,       pointer, contiguous:: r(:,:), v(:,:), vg(:,:), w(:), dg(:)
  real(8),    pointer, contiguous:: t(:)
  integer(2), pointer, contiguous:: q(:,:)
  integer(8), pointer, contiguous:: id(:)
  real(8):: llc
  integer:: n, m, i1, i2, i3, i4, j
  integer, save:: n_tot = 0, i_max = 0, o_max = 0, itimer(2) = 0
  integer, parameter::particle_solver_key = 3, EOF = -1
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%voutput(1)', itimer=itimer(1))
  n = 0
  m = 0
  do i4=1,n_bin
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      pp => self%particles(i1,i2,i3,i4)
      call harvest (.false.)
    end do
    end do
    end do
  end do
  !m = (product(self%ub-self%lb+1)*particles_per_bin*n_bin)/output_every + 1
  !print *,m,n
  m = n
  allocate (q(m,4), r(m,3), v(m,3), vg(m,3), dg(m), w(m), t(m), id(m))
  !-----------------------------------------------------------------------------
  n = 0
  do i4=1,n_bin
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      pp => self%particles(i1,i2,i3,i4)
      call harvest (.true.)
    end do
    end do
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Write out info on dimensions -- before data to ensure it is up-to-date
  ! during live runs
  !-----------------------------------------------------------------------------
  call io%lock%set ('particle_solver_t%voutput')
  call trace%begin ('particle_solver_t%voutput(2)', itimer=itimer(2))
  rewind (io_unit%particle_i)
  i_max = max(i_max,maxval(id))
  o_max = max(o_max,self%iout)
  write (io_unit%particle_i) n_bin, i_max+1, o_max
  flush (io_unit%particle_i)
  !-----------------------------------------------------------------------------
  ! Write out buffered data
  !-----------------------------------------------------------------------------
  write (io_unit%particle_d) particle_solver_key
  write (io_unit%particle_d) &
   self%iout, n, self%id, self%time, self%patch%position-self%patch%size/2, &
   self%patch%ds
  write (io_unit%particle_d) q(1:n,:), r(1:n,:), v(1:n,:), vg(1:n,:), &
    w(1:n), dg(1:n), t(1:n), id(1:n)
  write (io_unit%particle_d) EOF
  backspace (io_unit%particle_d)
  n_tot = n_tot + n
  if (mod(n_tot,1000000) == 0) then
    write(stdout,'(a,f8.2,a)') 'particle_io:', n_tot*1e-6, ' Mp'
  end if
  call trace%end (itimer(2))
  call io%lock%unset ('particle_solver_t%voutput')
  !-----------------------------------------------------------------------------
  deallocate (q, r, v, vg, dg, w, t, id)
  call trace%end (itimer(1))
contains

subroutine harvest (copy)
  logical:: copy
  integer:: ng(4)
  ng(1:3) = self%patch%mesh%ng
  ng(4) = 0
  p => pp%head%next
  do j=1,pp%n
    select type (p)
    type is (vparticle_t)
      n = n + 1
      if (copy) then
        q (n,:) = p%q-ng-1
        r (n,:) = p%r
        v (n,:) = p%v
        vg(n,:) = self%vg(:,p%q(1),p%q(2),p%q(3))
        dg(n)   = self%dg(  p%q(1),p%q(2),p%q(3))
        w (n)   = p%w
        t (n)   = p%t
        id(n)   = p%id-1
      end if
    end select
    p => p%next
  end do
end subroutine harvest
END SUBROUTINE voutput

!===============================================================================
SUBROUTINE statistics (self)
  class(particle_solver_t):: self
  !-----------------------------------------------------------------------------
  call self%particles(1,1,1,1)%statistics
END SUBROUTINE statistics

!===============================================================================
!> Update the global momentum, by adding the current value and subtracting
!> the previous value
!===============================================================================
SUBROUTINE update_momentum (self)
  class(particle_solver_t):: self
  integer:: i, j
  type(mpi_counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  do j=1,2
  do i=1,3
    counter => momentum_counter(i,j)
    total_momentum(i,j) = &
      counter%update (self%new_momentum(i,j) - self%old_momentum(i,j))
    self%old_momentum(i,j) = self%new_momentum(i,j)
  end do
  end do
END SUBROUTINE update_momentum

!===============================================================================
!> Return global momentum, by category
!===============================================================================
FUNCTION global_momentum (self)
  class(particle_solver_t):: self
  real(8):: global_momentum(3,3)
  !.............................................................................
  global_momentum(:,1) = total_momentum(:,1)
  global_momentum(:,2) = total_momentum(:,2)
  global_momentum(:,3) = total_momentum(:,1) + total_momentum(:,2)
END FUNCTION global_momentum

!===============================================================================
!> Allow at most mostion over one cell, given the current Courant number
!===============================================================================
SUBROUTINE courant_condition (self)
  class(particle_solver_t):: self
  class(gpatch_t), pointer:: patch
  real:: dust_velocity
  !.............................................................................
  call trace%begin ('particle_solver_t%courant_condition')
  patch => self%patch
  dust_velocity = self%vmax*minval(patch%ds)
  patch%u_max = max(patch%u_max, dust_velocity)
  self%dtime = patch%courant*minval(patch%ds)/patch%u_max       ! dt <= 1/vmax
  if (verbose > 1) then
    write(stdlog,'(a,i6,1p,7e11.2)') &
    'particle_solver_t%courant: id, u_max, vmax, dtime =', &
    patch%id, patch%u_max, self%vmax*maxval(patch%ds), self%dtime
  end if
  call trace%end ()
END SUBROUTINE courant_condition

END MODULE particle_solver_mod
