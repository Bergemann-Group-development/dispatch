!==============================================================================
!> Particle data types for dust-gas interaction
!===============================================================================
!define DBG
MODULE particle_solver_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE omp_lib
  USE omp_mod
  USE omp_timer_mod
  USE omp_lock_mod
  USE mpi_mod
  USE task_mod
  USE bits_mod
  USE mesh_mod
  USE link_mod
  USE halo_mod
  USE list_mod
  USE gpatch_mod
  !USE dust_evolve_mod
  USE timer_mod
  USE scratch_mod
  USE math_mod
  USE eos_mod
  USE units_mod
  USE scaling_mod
  USE shared_mod
  USE particles_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Everything needed by a particle solver
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: particle_solver_t
    type(particles_t), pointer :: particles     => null()  ! particles to move
    class(gpatch_t),   pointer :: patch         => null()  ! pointer to patch
    class(link_t),     pointer :: link          => null()  ! task link
    logical                    :: on = .false.             ! on/off
    logical                    :: forced = .false.         ! external force
    logical                    :: initialized = .false.    ! ICs done
    logical, dimension(3)      :: lower                    ! has lower boundary
    logical, dimension(3)      :: upper                    ! has upper boundary
    integer, dimension(3)      :: lb=0, ub=0               ! lower and upper index bounds
    integer                    :: bytes_per_particle       ! size of one
    integer                    :: bytes_per_particles      ! size of all
    real                       :: vmax                     ! max particle speed
    real                       :: v_min, v_max             ! min/max u_phase
    real(8)                    :: new_momentum(3,2)=0.     ! new patch momentum
    real(8)                    :: old_momentum(3,2)=0.     ! old patch momentum
    !---------------------------------------------------------------------------
    ! Ad hoc added
    !---------------------------------------------------------------------------
    integer                    :: mp
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: initial_condition
    procedure:: pre_update
    procedure:: post_update
    procedure:: update
  end type
  !-----------------------------------------------------------------------------
  ! Total momentum is counted using the mpi_counter data type.  The
  ! mpi_counter_t%update (increment) returns the value accumulated over ranks
  !-----------------------------------------------------------------------------
  type(mpi_counter_t), target:: momentum_counter(3,2) ! particles and gas
  real(8), target:: total_momentum(3,2)               ! their actual values
CONTAINS

!===============================================================================
!> Initialize particle solver, reading parameters and adding task link
!> FIXME: Many of these parameters belong in particles_mod
!===============================================================================
!SUBROUTINE pre_init (self, patch)
SUBROUTINE pre_init (self)
  class(particle_solver_t), target  :: self
  class(gpatch_t), pointer          :: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pre_init')
  !---------------------------------------------------------------------------
  ! Add a task link, in preparation for init_task_list
  !-----------------------------------------------------------------------------
  allocate (halo_t:: self%link)
  self%link%task => self
  self%iout = 0
  self%out_next = 0d0
  !-----------------------------------------------------------------------------
  ! Get particle parameters
  !-----------------------------------------------------------------------------
  allocate (self%particles)
  call self%particles%pre_init (self%patch)
  self%on = self%particles%on
 call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE init (self)
  class(particle_solver_t)    :: self
  class(particles_t), pointer :: pp, pi
  class(gpatch_t),    pointer :: patch
  class(nbor_t),      pointer :: nbor, prev
  class(task_t),      pointer :: nbtask
  integer                     :: m(4), i, j, n, l(3), u(3)
  logical:: first_time = .true.
  character(len=120), save:: id= &
  '$Id$'
  !-----------------------------------------------------------------------------
  ! Initialize particle lists and dust density array
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%print_id (id)
  call trace%begin ('particle_solver_t%init')
  patch => self%patch
  self%lb = patch%li
  self%ub = patch%ui
  self%id = patch%id
  self%time = patch%time
  self%dtime = patch%dtime
  m(1:3) = patch%mesh%gn
  m(4) = 3
  !-----------------------------------------------------------------------------
  ! Allocate force array
  !-----------------------------------------------------------------------------
  call os%allocate (scratch%force_per_unit_volume, m)
  if (.not.associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init ('part')
  end if
  !-----------------------------------------------------------------------------
  ! Allocate particle data, propagating index bounds
  !-----------------------------------------------------------------------------
  l = self%lb
  u = self%ub
  pp => self%particles
  pp%lb = l
  pp%ub = u
  pp%m =  pp%per_cell*(2*product(u-l+1)-product(u-l))
  call pp%init (self%patch, pp%m)
  allocate (pp%exports)
  allocate (pp%imports)
  pi => self%particles%imports
  pi%lb = l
  pi%ub = u
  pi%m =  pp%per_cell*(product(u-l+1)-product(u-l))
  call pi%init (self%patch, pi%m)
  !-----------------------------------------------------------------------------
  ! Compute the particle memory footprint
  !-----------------------------------------------------------------------------
  n = pp%n
  self%bytes_per_particle = ( &
    storage_size(pp%p)*3 + &
    storage_size(pp%q)*3 + &
    storage_size(pp%v)*3 + &
    storage_size(pp%i) + &
    storage_size(pp%w) + &
    storage_size(pp%s))/8
  self%bytes_per_particles = n*self%bytes_per_particle
  call io%bits_mem (self%bytes_per_particle*8, n, 'particles')
  !.............................................................................
  call io%lock%set
  if (mpi%master .and. first_time) then
    write (stdout,'(i3,a,i4,a,f7.3,a)') &
      self%bytes_per_particle             , ' bytes per particle,', &
      self%particles%per_cell             , ' particles per cell,', &
      (self%bytes_per_particles + 4* &
       product(shape(patch%mem)))/1024.**2, ' MB/patch'
    first_time = .false.
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Set boundary flags
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
  !-----------------------------------------------------------------------------
  ! A few final inits
  !-----------------------------------------------------------------------------
  call self%aux%init
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
  !-----------------------------------------------------------------------------
  if (.not. self%particles%on .or. self%initialized) return
  self%initialized = .true.
  call trace%begin ('particle_solver_t%initial_condition')
  !-----------------------------------------------------------------------------
  call self%particles%initial_condition ()
  call trace%end ()
END SUBROUTINE initial_condition

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE pre_update (self)
  class(particle_solver_t), target:: self
  !.............................................................................
  call self%particles%pre_update
END SUBROUTINE pre_update

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE post_update (self)
  class(particle_solver_t):: self
  !.............................................................................
  call self%particles%post_update
END SUBROUTINE post_update

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t)    :: self
  !.............................................................................
  real, dimension(3)          :: r, dr
  real(8)                     :: time(2), v(3), wc, start, self_time
  real                        :: ds, dt, dd
  integer                     :: n, np, ne, nk, nm, nu, nx, m(3), ii(3), i, j
  integer(kind=2)             :: i1, i2, i3, l(3), u(3), q(3)
  integer                     :: id, paid(5), istep, n_add, n_rem, n_norem
  integer                     :: dim(4), n_update
  integer(8)                  :: pid
  logical                     :: move, export
  class(gpatch_t),   pointer  :: patch
  type(particles_t), pointer  :: exports
  type(particles_t), pointer  :: p
  real                        :: d, up, pc, pd(3), pa(3), w, pv(3), dv(3), pwd
  real                        :: e, cs, pc0, patch_position(3), pq(3)
  real                        :: vp(2), grad_phi(3), force(3)
  real(8)                     :: times(15), wc1, wc2, sum
  real(8)                     :: momentum(3,2), dvol
  real, pointer, contiguous   :: p1(:,:,:), p2(:,:,:), p3(:,:,:)
  type(mpi_counter_t), pointer:: counter
  integer, save               :: itimer=0, nprint=20
  !-----------------------------------------------------------------------------
  if (.not. self%on) return
  call trace%begin ('particle_solver_t%update', itimer=itimer)
  start = wallclock()
  patch => self%patch
  id = patch%id
  istep = patch%istep
  self%dtime = patch%dtime
  self%time = self%time + self%dtime
  self_time = self%time
  self%new_momentum = 0d0
  self%vmax = 0.0
  !-----------------------------------------------------------------------------
  ! Disable motion loop in test cases
  !-----------------------------------------------------------------------------
  dt = patch%dtime
  call self%particles%update (dt)
  !-----------------------------------------------------------------------------
  ! Instead of storing the new momentum directly, store it as a force that
  ! will end up updating the momentum with exactly the same amount, but will
  ! also be used in the solver predictor step.
  !-----------------------------------------------------------------------------
  if (self%particles%backreaction .and. self%dtime > 0) then
    p1 => patch%mem(:,:,:,patch%idx%px,patch%it,1)
    p2 => patch%mem(:,:,:,patch%idx%py,patch%it,1)
    p3 => patch%mem(:,:,:,patch%idx%pz,patch%it,1)
    do i3=self%lb(3),self%ub(3)
    do i2=self%lb(2),self%ub(2)
    do i1=self%lb(1),self%ub(1)
      d = self%particles%dg(  i1,i2,i3)
      v = self%particles%vg(:,i1,i2,i3)
      force(1) = (d*v(1)*patch%ds(1) - p1(i1,i2,i3))/self%dtime
      force(2) = (d*v(2)*patch%ds(2) - p2(i1,i2,i3))/self%dtime
      force(3) = (d*v(3)*patch%ds(3) - p3(i1,i2,i3))/self%dtime
      scratch%force_per_unit_volume(i1,i2,i3,:) = &
        scratch%force_per_unit_volume(i1,i2,i3,:) + force
    end do
    end do
    end do
  end if
  !---------------------------------------------------------------------------
  ! Maintain up-to-date momentum measures of gas and particles separately,
  ! based on increments, for consistency checks
  !---------------------------------------------------------------------------
  dvol = product(patch%ds)
  !momentum(:,1) = momentum(:,1)*patch%ds          ! particle momentum
  !momentum(:,2) = d*v*patch%ds                    ! gas momentum
  !self%new_momentum(:,1) = &
  !self%new_momentum(:,1) + momentum(:,1)*dvol     ! measured particle momentum
  !self%new_momentum(:,2) = &
  !self%new_momentum(:,2) + momentum(:,2)*dvol     ! measured gas momentum
  !-----------------------------------------------------------------------------
  timer%t_particles = timer%t_particles + (wallclock()-start)
  timer%n_particles = timer%n_particles + self%particles%n
  call trace%end (itimer)
END SUBROUTINE update

END MODULE particle_solver_mod
