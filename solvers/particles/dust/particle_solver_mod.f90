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
!>
!> While updating, the task local particle bins are used to accumulate export
!> particles, which are then placed in visibke task%exports particle bins, from
!> where they are removed by other tasks.  Both the nbor tasks and the owner
!> task must use task locks to avoid data races; i.e. both readin and writing
!> needs to be protected.
!===============================================================================
!define DBG
MODULE particle_solver_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE omp_lib
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
    procedure:: density
    procedure:: velocities
    procedure:: dispersion
    procedure:: exports_info
  end type
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
  real:: eps=1e-5                             ! momentum error tolerance
  real(8):: out_time = -1.0                   ! output cadence
  logical:: particle_io = .false.             ! output particle positions?
  logical:: backreaction = .false.            ! backreaction on gas
  logical:: detailed_timer = .false.          ! precision timing
  logical:: do_supersonic = .true.            ! supersonic dust friction
  logical:: debug=.false.                     ! check things
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
    backreaction, do_supersonic, debug
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
  class(gpatch_t), pointer :: patch
  integer                  :: m(4), i, j
  character(len=120), save:: id= &
  '$Id$'
  !-----------------------------------------------------------------------------
  ! Initialize particle lists and dust density array
  !-----------------------------------------------------------------------------
  call trace%print_id (id)
  ! if (.not.on) return
  call trace%begin ('particle_solver_t%init')
  m(1:3) = self%patch%mesh%gn
  allocate (self%particles(m(1),m(2),m(3),n_bin)) 
  call io%bits_mem (storage_size(self%particles), &
                   product(shape(self%particles)), 'particles')
  allocate (self%exports(n_bin))
  call pool%init  ([self%patch%id,0,0,0,0], 0)
  call vpool%init ([self%patch%id,0,0,0,0], 1)
  do i=1,n_bin
    call self%exports(i)%init ([self%patch%id,0,0,0,i], 0)
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
  !-----------------------------------------------------------------------------
  self%upper = .false.
  self%lower = .false.
  patch => self%patch
  self%position = patch%position
  do i=1,3
    self%lb(i) = patch%mesh(i)%li
    self%ub(i) = patch%mesh(i)%ui
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
  class(particle_t), pointer  :: p
  type(particle_t), pointer   :: pn
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
  !-----------------------------------------------------------------------------
  patch => self%patch
  extra_bits = 200
  if (self%patch%id == 1) then
    allocate (p)
    write (stdout,'(a,4i6)') 'particle storaga_size:', &
      storage_size(p)
    write (stdout,'(a,6i6)') 'particle storaga_size:', &
      storage_size(p%t), storage_size(p%r), &
      storage_size(p%q), storage_size(p%v), storage_size(p%w), &
      storage_size(p%t) + 3*storage_size(p%r) + 4*storage_size(p%q) + &
      storage_size(p%w) + 3*storage_size(p%v) + 2*64
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
        call particles%init (id, 0)
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
            allocate(pn)
            p => pn
            call p%init_vparticle
            p%cell%id = 1 + omp%thread_global + (lid+nv)*omp%nthreads_global
            nv = nv+1
            p%cell%dg = d
            p%cell%vg = vg
          else
            allocate(p)
            call p%init
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
  real(8)                     :: time(2), v(3), wc, start, self_time
  real                        :: ds, dt, os, n_min, n_max, dd
  integer                     :: n, np, ne, nk, nm, nu, nx, m(3), ii(3), j
  integer(kind=2)             :: i, i1, i2, i3, i4, l(3), u(3), q(3)
  integer                     :: id, paid(5), istep, n_add, n_rem, n_norem
  integer                     :: dim(4), n_update
  integer(8)                  :: pid
  logical                     :: move, export
  class(gpatch_t), pointer    :: patch
  type(particles_t), pointer  :: exports(:), exportp
  type(particles_t), pointer  :: particles, bins(:)
  real                        :: d, up, pc, pd(3), pa(3), w, pv(3), dv(3), pwd
  real                        :: e, cs, pc0
  real(8)                     :: times(15), wc1, wc2, sum
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
  self%dust%d => patch%mem(:,:,:,patch%idx%d              ,patch%new,1)
  self%dust%p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%new,1)
  !-----------------------------------------------------------------------------
  ! Allocate array for number of particles per cell-bin
  !-----------------------------------------------------------------------------
  m = patch%gn
  allocate (self%vg(3,m(1),m(2),m(3)))
  call io%bits_mem (storage_size(self%vg), &
                   product(shape(self%vg)), 'vg')
  self%vg = 0.0
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
    self%vg(i,i1,i2,i3) = patch%mem(i1,i2,i3,patch%idx%px+i-1,patch%new,1) &
                        / (self%dg(i1,i2,i3)*patch%ds(i))
    if (abs(self%vg(i,i1,i2,i3)) > 1e1/patch%ds(1)) then
      write(stderr,'(a,2i6,4i4,2(2x,1p,3e11.2))') &
      'warning, large vg:', patch%id,patch%istep,i1,i2,i3,i, &
      patch%mem(i1,i2,i3,patch%idx%px+i-1,patch%new,1), &
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
    call exports(i4)%init (paid, 0)
  end do
  !-----------------------------------------------------------------------------
  ! Non-vectorized pusher
  !-----------------------------------------------------------------------------
  ds = minval(patch%ds, mask=patch%mesh%n > 1)
  n_min = particles_per_bin*min_fraction
  n_max = particles_per_bin*max_fraction
  !-----------------------------------------------------------------------------
  ! Import from other patches.
  !-----------------------------------------------------------------------------
  call self%import
  !-----------------------------------------------------------------------------
  !> Macro for performance measurement -- no impaact when TIMER is undefined
  !-----------------------------------------------------------------------------
!define TIMER
#ifdef TIMER
  wc1=omp_get_wtime()
  times = 0d0
#define TIMEIT(i) wc2=omp_get_wtime(); times(i)=times(i)+(wc2-wc1); wc1=wc2
#else
#define TIMEIT(i) 
#endif
  !-----------------------------------------------------------------------------
  ! After the import, the interior particles should be maximally smoothly
  ! distributed, which is thus a good point to measure the dust bin density,
  ! renormalize, and do output.
  !-----------------------------------------------------------------------------
  n_add = 0
  n_rem = 0
  n_norem = 0
  call check_count (self, self%np, 'before add/rem')
  do i4=1,n_bin
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      particles => self%particles(i1,i2,i3,i4)
      np = particles%n
      if (particles%n < n_min) then
        n_add = n_add+1
        if (verbose > 3) &
          write(stdlog,'(2i6,3i4,2x,a,i4)') &
            patch%id,patch%istep,i1,i2,i3,'add',particles%n
        call add (self, particles)
        TIMEIT(1)
      else if (particles%n > n_max) then
        !-----------------------------------------------------------------------
        ! Potentially too many particles, but remove vparticles from the count
        ! before deciding, since they cannot be removed in any case.
        !-----------------------------------------------------------------------
        n = particles%n
        p => particles%head%next
        do j=1,particles%n
          if (associated(p%cell)) then
            n = n-1
          end if
          p => p%next
        end do
        !-----------------------------------------------------------------------
        ! Actually remove & renormalize?
        !-----------------------------------------------------------------------
        if (n > n_max) then
          n_rem = n_rem+1
          if (verbose > 3) &
            write(stdlog,'(2i6,3i4,2x,a,i4)') &
              patch%id,patch%istep,i1,i2,i3,'rem',particles%n
          call remove (self, particles)
        else
          n_norem = n_norem+1
        end if
        TIMEIT(2)
      else
        TIMEIT(3)
      end if
      self%np = self%np + (particles%n - np)
    end do
    end do
    end do
    !---------------------------------------------------------------------------
    ! Copy edges from the next interior layer -- this is better than relying on
    ! particles that may be missing due to delayed import.  The dust densities
    ! are ONLY used for diagnostics and dust size evolution!
    !---------------------------------------------------------------------------
  end do
  !-----------------------------------------------------------------------------
  ! Do output, if arrived at self%out_time
  !-----------------------------------------------------------------------------
  self%time = patch%time
  call self%output
  TIMEIT(4)
  !-----------------------------------------------------------------------------
  ! Time step control
  !-----------------------------------------------------------------------------
  self%dtime = patch%dtime
  self%time = self%time + self%dtime
  self_time = self%time
  !-----------------------------------------------------------------------------
  ! Now move
  !-----------------------------------------------------------------------------
  nu = 0
  nk = 0
  ne = 0
  nm = 0
  self%new_momentum = 0d0
  self%vmax = 0.0
  self%np = 0
  n_update = 0
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    v(1) = self%vg(1,i1,i2,i3)
    v(2) = self%vg(2,i1,i2,i3)
    v(3) = self%vg(3,i1,i2,i3)
    d = patch%mem(i1,i2,i3,patch%idx%d,patch%new,1)
    if (do_supersonic) then
      e = patch%mem(i1,i2,i3,patch%idx%e,patch%new,1)
      if (patch%gamma==1d0) then
        cs = sqrt(e/d)
      else
        cs = sqrt(e*patch%gamma*(patch%gamma-1d0)/d)
      end if
      cs = cs/ds
    end if
    ii = [i1,i2,i3]
    momentum(:,1) = 0d0
    do i4=1,n_bin
      particles => self%particles(i1,i2,i3,i4)
      nu = nu + particles%n
      os = drag_coeff/self%s_cgs(i4)
      n_update = n_update + particles%n
      paid = particles%id
      up = patch%u_phase(i1,i2,i3)/ds
      pc0 = d*os
      np = particles%n
      p => particles%head%next
      do j=1,np
        TIMEIT(5)
        next => p%next
        dt = 0.5*(self_time - p%t)
        if (dt > 0.0) then
          TIMEIT(6)
          !- kick 1 ------------------------------------------------------------
          pa(1) = 0.0                                     ! external force
          pa(2) = 0.0                                     ! external force
          pa(3) = 0.0                                     ! external force
          if (do_supersonic) then
            dv(1) = v(1)-p%v(1)                           ! speed difference
            dv(2) = v(2)-p%v(2)                           ! speed difference
            dv(3) = v(3)-p%v(3)                           ! speed difference
            pc = pc0*sqrt(1.0 &
               + (dv(1)**2+dv(2)**2+dv(3)**2)/cs**2)      ! supersonic Epstein
            pd(1) = pc*dv(1)                              ! drag force
            pd(2) = pc*dv(2)                              ! drag force
            pd(3) = pc*dv(3)                              ! drag force
          else
            pc = pc0
            pd(1) = pc*(v(1)-p%v(1))                      ! drag force
            pd(2) = pc*(v(2)-p%v(2))                      ! drag force
            pd(3) = pc*(v(3)-p%v(3))                      ! drag force
          end if
          pwd = p%w/d                                     ! density ratio
          !---------------------------------------------------------------------
          ! To allow extreme dust:gas ratios, compute the largest action first
          !---------------------------------------------------------------------
          if (backreaction) then
            dv(1) = dt*(pa(1)+pd(1))/(1.+pc*dt*(1.+pwd))  ! velocity delta
            dv(2) = dt*(pa(2)+pd(2))/(1.+pc*dt*(1.+pwd))  ! velocity delta
            dv(3) = dt*(pa(3)+pd(3))/(1.+pc*dt*(1.+pwd))  ! velocity delta
            p%v(1) = p%v(1) + dv(1)                       ! particle update
            p%v(2) = p%v(2) + dv(2)                       ! particle update
            p%v(3) = p%v(3) + dv(3)                       ! particle update
            v(1) = v(1) - dv(1)*pwd                       ! gas update
            v(2) = v(2) - dv(2)*pwd                       ! gas update
            v(3) = v(3) - dv(3)*pwd                       ! gas update
            p%r(1) = p%r(1) + 2.*dt*p%v(1)                ! forward two dr
            p%r(2) = p%r(2) + 2.*dt*p%v(2)                ! forward two dr
            p%r(3) = p%r(3) + 2.*dt*p%v(3)                ! forward two dr
            p%t = self_time
            pd(1) = pc*(v(1)-p%v(1))                      ! drag force
            pd(2) = pc*(v(2)-p%v(2))                      ! drag force
            pd(3) = pc*(v(3)-p%v(3))                      ! drag force
            dv(1) = dt*(pa(1)+pd(1))/(1.+pc*dt*(1.+pwd))  ! velocity delta
            dv(2) = dt*(pa(2)+pd(2))/(1.+pc*dt*(1.+pwd))  ! velocity delta
            dv(3) = dt*(pa(3)+pd(3))/(1.+pc*dt*(1.+pwd))  ! velocity delta
            p%v(1) = p%v(1) + dv(1)                       ! particle update
            p%v(2) = p%v(2) + dv(2)                       ! particle update
            p%v(3) = p%v(3) + dv(3)                       ! particle update
            v(1) = v(1) - dv(1)*pwd                       ! gas update
            v(2) = v(2) - dv(2)*pwd                       ! gas update
            v(3) = v(3) - dv(3)*pwd                       ! gas update
          else
            p%v(1) = p%v(1) + dt*(pa(1)+pd(1))/(1.+pc*dt) ! particle update
            p%v(2) = p%v(2) + dt*(pa(2)+pd(2))/(1.+pc*dt) ! particle update
            p%v(3) = p%v(3) + dt*(pa(3)+pd(3))/(1.+pc*dt) ! particle update
            p%r(1) = p%r(1) + 2.*dt*p%v(1)                ! forward two dr
            p%r(2) = p%r(2) + 2.*dt*p%v(2)                ! forward two dr
            p%r(3) = p%r(3) + 2.*dt*p%v(3)                ! forward two dr
            p%t = self_time
            pd(1) = pc*(v(1)-p%v(1))                      ! drag force
            pd(2) = pc*(v(2)-p%v(2))                      ! drag force
            pd(3) = pc*(v(3)-p%v(3))                      ! drag force
            p%v(1) = p%v(1) + dt*(pa(1)+pd(1))/(1.+pc*dt) ! particle update
            p%v(2) = p%v(2) + dt*(pa(2)+pd(2))/(1.+pc*dt) ! particle update
            p%v(3) = p%v(3) + dt*(pa(3)+pd(3))/(1.+pc*dt) ! particle update
          end if
          momentum(1,1) = momentum(1,1) + p%v(1)*p%w      ! particle momentum
          momentum(2,1) = momentum(2,1) + p%v(2)*p%w      ! particle momentum
          momentum(3,1) = momentum(3,1) + p%v(3)*p%w      ! particle momentum
          TIMEIT(7)
          !- update q+r representation -----------------------------------------
          q(1) = floor(p%r(1))                            ! integer shift
          q(2) = floor(p%r(2))                            ! integer shift
          q(3) = floor(p%r(3))                            ! integer shift
          p%r(1) = p%r(1) - q(1)                          ! subtract
          p%r(2) = p%r(2) - q(2)                          ! subtract
          p%r(3) = p%r(3) - q(3)                          ! subtract
          p%q(1) = p%q(1) + q(1)                          ! add
          p%q(2) = p%q(2) + q(2)                          ! add
          p%q(3) = p%q(3) + q(3)                          ! add
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
          TIMEIT(8)
          !---------------------------------------------------------------------
          ! Export particle to another patch
          !---------------------------------------------------------------------
          if (export) then
            call particles%remove (p)
            call exports(i4)%append (p)
            ne = ne+1
            if (verbose > 2) then
              write (stdlog,*) patch%id, 'export: q=', patch%id, p%q
            end if
            TIMEIT(9)
          !---------------------------------------------------------------------
          ! Move particle to another cell. In case the cell has been visited,
          ! the dust bin density is updated as well.
          !---------------------------------------------------------------------
          else if (move) then
            self%vmax = max(self%vmax,abs(p%v(1)),abs(p%v(2)),abs(p%v(3)))
            call particles%remove (p)
            call self%particles(p%q(1),p%q(2),p%q(3),p%q(4))%append (p)
            nm = nm+1
            TIMEIT(10)
          !---------------------------------------------------------------------
          ! Particle stays in this cell -- update vmax
          !---------------------------------------------------------------------
          else
            self%vmax = max(self%vmax,abs(p%v(1)),abs(p%v(2)),abs(p%v(3)))
            TIMEIT(11)
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
    self%new_momentum(:,1) + momentum(:,1)*dvol     ! measured particle momentum
    self%new_momentum(:,2) = &
    self%new_momentum(:,2) + momentum(:,2)*dvol     ! measured gas momentum
  end do
  end do
  end do
  TIMEIT(12)
  associate (n_particles => timer%n_particles)
  !$omp atomic update
  n_particles = n_particles + n_update
  end associate
  self%ne = ne
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
  TIMEIT(13)
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
  TIMEIT(14)
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
      write (stdlog,'(i6,5i12,f6.3,"s",6p,2f8.3," mus/p",0p,2f6.2)') &
      patch%id, np, nk, n_add, n_rem, n_norem, wc, wc/np, wc/nk, &
      real(nm)/nk, real(ne)/nk
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
  TIMEIT(15)
#ifdef TIMER
  write (stdout,'("times",3p,15f8.2)') times
#endif
  !-----------------------------------------------------------------------------
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Top up particles per size bin.  Note that if there are no particles it is 
!> better to do nothing than to add particles, since they would either have to
!> have zero weight, or they would add additional dust mass.
!===============================================================================
SUBROUTINE add (self, particles)
  class(particle_solver_t)    :: self
  type(particles_t), pointer  :: particles
  class(particle_t), pointer  :: p, next
  real                        :: w, m0(3), m1(3), v(3), dd
  real(8)                     :: t
  integer                     :: n, i, j, q(4), n_add
  !-----------------------------------------------------------------------------
  n = particles%n
  if (n == 0) then
    if (debug) then
      write (stdout,'(a,i6,4i4)') &
        'particle_solver_t%add WARNING: zero particles', particles%id
    end if
    return
  end if
  !-----------------------------------------------------------------------------
  m0 = 0.0                              ! initial momentum
  dd = 0.0                              ! dust bin density
  p => particles%head%next
  do j=1,n
    dd = dd + p%w
    m0 = m0 + p%v*p%w
    p => p%next
  end do
  !-----------------------------------------------------------------------------
  t = 0d0                               ! time
  m1 = 0.0                              ! momentum after renormalization
  w = dd/particles_per_bin              ! new weights
  p => particles%head%next              ! first particle
  do j=1,n
    p%w = w
    m1 = m1 + p%v*p%w
    t = t + p%t
    p => p%next
  end do
  t = t/n                               ! average time
  !-----------------------------------------------------------------------------
  n_add = particles_per_bin-n           ! number to add
  if (w==0.0 .or. n_add==0) then
    v = 0.0                             ! rare case, or many zero-w particles
  else
    v = (m0-m1)/(n_add*w)               ! speed to match momentum
  end if
  !-----------------------------------------------------------------------------
  q(1:4) = particles%id(2:5)            ! cell indices, bin number
  do j=1,n_add
    p => pool%pop ()                    ! re-cycled particle
    do i=1,3
      p%q(i) = q(i)                     ! place it
      p%v(i) = v(i)                     ! kick it
      p%r(i) = self%random%ran1()       ! stir it
    end do
    p%q(4) = q(4)                       ! bin it
    p%w = w                             ! weigh it
    p%t = t                             ! time it
    call particles%append (p)           ! add it
  end do
  !-----------------------------------------------------------------------------
  if (debug) then
    m1 = particles%total_momentum()
    if (sum((m0-m1)**2) > eps*sum((m0+m1)**2)) then
      write (stderr,*) 'particle_solver_t%add error:', m0, m1
    end if
  end if
END SUBROUTINE add

!===============================================================================
!> Bring particles per bin below max
!===============================================================================
SUBROUTINE remove (self, particles)
  class(particle_solver_t)    :: self
  type(particles_t), pointer  :: particles
  class(particle_t), pointer  :: p, next
  real                        :: dd, dw, mm(3), dm(3), m0(3), m1(3)
  integer                     :: n0, n, j
  logical                     :: kill
  real                        :: eps=1e-5
  !-----------------------------------------------------------------------------
  if (debug) then
    m0 = particles%total_momentum()
  end if
  !-----------------------------------------------------------------------------
  ! Remove particles into the threadprivate re-use pool, while measuring the
  ! mass and momentum removed
  !-----------------------------------------------------------------------------
  dd = 0.0
  mm = 0.0
  n0 = particles%n
  p => particles%head%next
  do j=1,n0
    next => p%next
    if (particles%n > particles_per_bin) then
      if (.not.associated(p%cell)) then
        dd = dd + p%w
        mm = mm + p%w*p%v
        call particles%remove (p)
        call pool%append (p)
      end if
    else
      exit
    end if
    p => next
  end do
  !-----------------------------------------------------------------------------
  ! Adjust mass and momentum of remaining particles
  !-----------------------------------------------------------------------------
  n = particles%n
  dw = dd/n                     ! mass to add to each particle
  if (dw > 0.0) then
    dm = mm/n                   ! momentum to add to each particle
    p => particles%head%next
    do j=1,n
      p%w = p%w + dw            ! times np = dd
      p%v = p%v*(p%w-dw)/p%w    ! conserve momentum with new weight
      p%v = p%v + dm/p%w        ! p%w weighted sum = dm
      p => p%next
    end do
  else
    write (stdout,'(a,i6,4i4,2i5,1p,4e11.2)') &
      'particle_solver_t%remove WARNING:', particles%id, n0, n, dd, dm
    p => particles%head%next
    do j=1,n
      if (associated(p%cell)) then
        write (stdout,'(a,i4,2x,4i4,3f8.3,1p,2(2x,4e12.3))') &
          'remove WARNING:', j, p%q, p%r, p%w, p%v, p%cell%dg, p%cell%vg
      else
        write (stdout,'(a,i4,2x,4i4,3f8.3,1p,4e12.3)') &
          'remove WARNING:', j, p%q, p%r, p%w, p%v
      end if
      p => p%next
    end do
  end if
  !.............................................................................
  if (debug) then
    m1 = particles%total_momentum()
    if (sum((m0-m1)**2) > eps*sum((m0+m1)**2)) then
      write (stderr,'(a,5i4,3(2x,1p,3e11.2))') &
      'particle_solver_t%remove error:', particles%id, m0, m1, &
      sum((m0-m1)**2), sum((m0+m1)**2)
    end if
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
SUBROUTINE import (self, count)
  class(particle_solver_t)   :: self
  logical, optional          :: count
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: nbtask
  class(particles_t), pointer:: particles
  class(particle_t), pointer :: p, next
  real(8)                    :: r(3), dist(3), t_nbors
  integer                    :: ni, ii(5), j, bin, n(4), n_before, n_after
  integer(kind=2)            :: q(3)
  logical                    :: inside
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  ! Loop over patch nbors
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%import',itimer=itimer)
  self%ni = 0
  n_before = 0
  n_after = 0
  t_nbors = huge(1.)
  nbor => self%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (particle_solver_t)
      t_nbors = min(t_nbors,nbtask%patch%time)
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
        n_before = n_before + nbtask%exports(bin)%n
        n(3) = self%np
        do j=1,nbtask%exports(bin)%n
          if (abs(maxval(p%q)) > 100) write (stdlog,*) 'import error', j, p%q
          next => p%next
          where (p%r==1d0) p%r=0.9999d0
          r = dist + real(p%r,kind=8) + real(p%q(1:3),kind=8)
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
            ni = ni + 1
            ii(1) = self%patch%id; ii(2:4) = p%q(1:3); ii(5) = p%q(4)
          else if (verbose > 2) then
            write (stdlog,'(2i6,3f10.3,2x,a,4(2x,3i6))') &
              self%patch%id,nbtask%id,dist,'rejected',p%q,q,self%lb,self%ub        
          end if
          p => next
        end do
        n_after = n_after + nbtask%exports(bin)%n
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
  real, pointer, contiguous     :: tmp_r(:,:), vr(:,:,:,:)
  real, pointer, contiguous     :: vdx(:,:,:,:)
  real, pointer, contiguous     :: vdy(:,:,:,:)
  real, pointer, contiguous     :: vdz(:,:,:,:)
  real(8), pointer              :: total_mom(:,:)
  integer, pointer, contiguous  :: tmp_i(:)
  integer                       :: i, j, n, i1, i2, i3, i4, n_out, m(3)
  integer, save                 :: itimer=0, n_print=5
  !-----------------------------------------------------------------------------
  if (.not. on .or. .not. io%do_output) return
  call trace%begin ('particle_solver_t%output', itimer=itimer)
  if (verbose > 1) then
    write (stdlog, '(a,i6,2f8.4)') &
      'particles_t%output: id, time, aux_next =', &
      self%patch%id, self%time, self%aux%aux_next
  end if
  !-----------------------------------------------------------------------------
  ! High-cadence output of v_particle_t particles
  !-----------------------------------------------------------------------------
  if (out_time > 0d0 .and. self%time >= self%out_next &
                     .and. self%time >= io%start_time) then
    self%out_next = max(self%out_next,io%start_time)
    n_out = 0
    do while (self%time >= self%out_next)
      if (n_print > 0) then
        n_out = n_out+1
        if (n_out > 1) then
          !$omp critical (voutput_cr)
          write (stdout,*) 'WARNING: duplicated output, cadence may be too high'
          n_print = n_print-1
          !$omp end critical (voutput_cr)
        end if
      end if
      call poutput (self)
      self%iout = self%iout + 1
      self%out_next = self%out_next + out_time
    end do
  end if
  if (self%time >= self%aux%aux_next) then
    call self%import (count=debug)
    call self%density
    m = self%patch%gn
    allocate (vr(m(1),m(2),m(3),n_bin))
    allocate (vdx(m(1),m(2),m(3),n_bin))
    allocate (vdy(m(1),m(2),m(3),n_bin))
    allocate (vdz(m(1),m(2),m(3),n_bin))
    call self%dispersion (vr)
    call self%velocities (vdx, vdy, vdz)
    call self%aux%register ('dust_s', self%s_cgs)
    call self%aux%register ('dust_d', self%dd)
    call self%aux%register ('dust_vrms', vr)
    call self%aux%register ('dust_vx'  , vdx)
    call self%aux%register ('dust_vy'  , vdy)
    call self%aux%register ('dust_vz'  , vdz)
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
    deallocate (vr, vdx, vdy, vdz)
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
SUBROUTINE poutput (self)
  class(particle_solver_t):: self
  class(particles_t), pointer:: pp
  class(particle_t), pointer:: p
  real,       pointer, contiguous:: r(:,:), v(:,:), vg(:,:), w(:), dg(:), dd(:)
  real(8),    pointer, contiguous:: t(:)
  integer(2), pointer, contiguous:: q(:,:)
  integer(8), pointer, contiguous:: id(:)
  real(8):: llc
  integer:: n, m, i1, i2, i3, i4, j
  integer, save:: n_tot = 0, i_max = 0, o_max = 0, itimer(3) = 0
  integer, parameter:: particle_solver_key = 3, EOF = -1
  integer, save:: iout = -99
  character(len=120):: file
  !-----------------------------------------------------------------------------
  if (.not. particle_io) return 
  call trace%begin ('particle_solver_t%poutput(1)', itimer=itimer(1))
  call self%lock%set ('particle_solver_t%poutput')
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
    pp => self%exports(i4)
    call harvest (.false.)
  end do
  !-----------------------------------------------------------------------------
  ! Only do I/O if there are particles to write
  !-----------------------------------------------------------------------------
  if (n > 0) then
    m = n
    allocate (q(m,4), r(m,3), v(m,3), vg(m,3), dg(m), dd(m), w(m), t(m), id(m))
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
      pp => self%exports(i4)
      call harvest (.true.)
    end do
    !---------------------------------------------------------------------------
    ! Re-open a new file for each snapshot
    !---------------------------------------------------------------------------
    call trace%begin ('particle_solver_t%poutput(2)', itimer=itimer(2))
    call io%lock%set ('poutput')
    if (iout < self%iout) then
      iout = self%iout
      close (io_unit%particle_d)
      close (io_unit%particle_i)
      call os%mkdir (trim(io%rundir)//'particles/')
      write (file,'(a,i5.5,a)') trim(io%rundir)//'particles/', self%iout, '.dat'
      open (unit=io_unit%particle_d, file=file, form='unformatted', &
            status='unknown')
      write (stdout,*) trim(file)//' opened'
      write (file,'(a,i5.5,a)') trim(io%rundir)//'particles/', self%iout, '.inf'
      open (unit=io_unit%particle_i, file=file, form='unformatted', &
            status='unknown')
    end if
    call trace%end (itimer(2))
    !---------------------------------------------------------------------------
    ! Write out info on dimensions -- before data to ensure it is up-to-date
    ! during live runs
    !---------------------------------------------------------------------------
    call trace%begin ('particle_solver_t%poutput(3)', itimer=itimer(3))
    rewind (io_unit%particle_i)
    i_max = max(i_max,maxval(id))
    o_max = max(o_max,self%iout)
    write (io_unit%particle_i) n_bin, i_max+1, o_max
    flush (io_unit%particle_i)
    !---------------------------------------------------------------------------
    ! Write out buffered data
    !---------------------------------------------------------------------------
    write (io_unit%particle_d) particle_solver_key
    write (io_unit%particle_d) &
    self%iout, n, self%id, self%time, self%patch%position-self%patch%size/2, &
    self%patch%ds
    write (io_unit%particle_d) q(1:n,:), r(1:n,:), v(1:n,:), vg(1:n,:), &
      w(1:n), dg(1:n), dd(1:n), t(1:n), id(1:n)
    write (io_unit%particle_d) EOF
    backspace (io_unit%particle_d)
    n_tot = n_tot + n
    if (mod(n_tot,1000000) == 0) then
      write(stdout,'(a,f8.2,a)') 'particle_io:', n_tot*1e-6, ' Mp'
    end if
    deallocate (q, r, v, vg, dg, dd, w, t, id)
    call trace%end (itimer(3))
    call io%lock%unset ('poutput')
  end if
  !-----------------------------------------------------------------------------
  call self%lock%unset ('particle_solver_t%poutput')
  call trace%end (itimer(1))
contains

subroutine harvest (copy)
  logical:: copy
  integer:: ng(4)
  real(8):: llc(3), ds(3), pos(3)
  !-----------------------------------------------------------------------------
  llc = self%patch%position - self%patch%size/2d0 
  ds = self%patch%ds
  ng(1:3) = self%patch%mesh%ng
  ng(4) = 0
  p => pp%head%next
  do j=1,pp%n
    if (associated(p%cell)) then
      n = n + 1
      if (copy) then
        q (n,:) = p%q-ng-1
        r (n,:) = p%r
        v (n,:) = p%v
        vg(n,:) = self%vg(:,p%q(1),p%q(2),p%q(3))
        dg(n)   = self%dg(  p%q(1),p%q(2),p%q(3))
        dd(n)   = self%dd(  p%q(1),p%q(2),p%q(3),p%q(4))
        w (n)   = p%w
        t (n)   = p%t
        id(n)   = p%cell%id-1
      end if
    end if
    p => p%next
  end do
end subroutine harvest
END SUBROUTINE poutput

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
!> Allow at most mostion over one cell, given the current Courant number.
!> self%vmax is updated by self%update(), and hence there is no need to wait
!> until after self%pre_update to call this routine.  It must, on the other
!> hand, be called after the HD/MHD solver has computed its patch%u_max.
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(particle_solver_t):: self
  logical, optional:: detailed_timer
  class(gpatch_t), pointer:: patch
  real:: dust_velocity
  !.............................................................................
  if (.not. on) return
  call trace%begin ('particle_solver_t%courant_condition')
  patch => self%patch
  dust_velocity = self%vmax*minval(patch%ds)
  patch%u_max = max(patch%u_max, dust_velocity, tiny(1.))
  self%dtime = patch%courant*minval(patch%ds)/patch%u_max       ! dt <= 1/vmax
  if (verbose > 1) then
    write(stdlog,'(a,i6,1p,7e11.2)') &
    'particle_solver_t%courant: id, u_max, vmax, dtime =', &
    patch%id, patch%u_max, self%vmax*maxval(patch%ds), self%dtime
  end if
  call trace%end ()
END SUBROUTINE courant_condition

!===============================================================================
!> Compute dust density
!===============================================================================
SUBROUTINE density (self)
  class(particle_solver_t):: self
  integer:: i1, i2, i3, i4
  call self%import
  do i4=1,self%n_bin
    do i3=self%lb(3),self%ub(3)
      do i2=self%lb(2),self%ub(2)
        do i1=self%lb(1),self%ub(1)
          self%dd(i1,i2,i3,i4) = self%particles(i1,i2,i3,i4)%density()
        end do
      end do
    end do
  end do
END SUBROUTINE density

!===============================================================================
!> Compute dust density
!===============================================================================
SUBROUTINE velocities (self, vdx, vdy, vdz)
  class(particle_solver_t):: self
  real, dimension(:,:,:,:):: vdx, vdy, vdz
  real:: vd(3)
  integer:: i1, i2, i3, i4
  !-----------------------------------------------------------------
  call self%import
  do i4=1,self%n_bin
    do i3=self%lb(3),self%ub(3)
      do i2=self%lb(2),self%ub(2)
        do i1=self%lb(1),self%ub(1)
          vd = self%particles(i1,i2,i3,i4)%velocity()*self%patch%ds
          vdx(i1,i2,i3,i4) = vd(1)
          vdy(i1,i2,i3,i4) = vd(2)
          vdz(i1,i2,i3,i4) = vd(3)
        end do
      end do
    end do
  end do
END SUBROUTINE velocities

!===============================================================================
!> Compute velocity diespersion
!===============================================================================
SUBROUTINE dispersion (self, vd)
  class(particle_solver_t):: self
  real:: vd(:,:,:,:), ds
  integer:: i1, i2, i3, i4
  !-----------------------------------------------------------------------------
  ds = self%patch%ds(1)
  do i4=1,self%n_bin
    do i3=self%lb(3),self%ub(3)
      do i2=self%lb(2),self%ub(2)
        do i1=self%lb(1),self%ub(1)
          vd(i1,i2,i3,i4) = ds*self%particles(i1,i2,i3,i4)%dispersion()
        end do
      end do
    end do
  end do
END SUBROUTINE dispersion

!===============================================================================
!> Debug information about exports particle bins
!===============================================================================
SUBROUTINE exports_info (self)
  class(particle_solver_t):: self
  class(particles_t), pointer:: particles
  class(particle_t), pointer:: p
  integer:: bin, j
  !.............................................................................
  if (.not. on) return
  call trace%begin ('particle_solver_t%exports_info')
  call self%lock%set ('particle_solver_t%exports_info')
  do bin=1,n_bin
    particles => self%exports(bin)
    p => particles%head%next
    do j=1,particles%n
      write (stdlog,'(a,i6,4i4,3f12.8)') &
        'particle_solver_t%exports_info: id, j, q, p =', &
        self%id, j, p%q(1:3), p%r
      p => p%next
    end do
  end do
  call self%lock%unset ('particle_solver_t%exports_info')
  call trace%end ()
END SUBROUTINE exports_info

END MODULE particle_solver_mod
