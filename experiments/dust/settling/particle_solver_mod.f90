!===============================================================================
!> Particle data types for dust-gas interaction:
!>
!>   particle_t         : data type for a single particle
!>   particles_t        : data type for lists of particles
!>   particle_solver_t  : data type for particle solvers
!>
!===============================================================================
!define DBG
MODULE particle_solver_mod
  USE io_mod
  USE trace_mod
  USE omp_timer_mod
  USE task_mod
  USE bits_mod
  USE mesh_mod
  USE link_mod
  USE list_mod
  USE gpatch_mod
  USE particles_mod
  USE dust_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Everything needed by a particle solver
  !-----------------------------------------------------------------------------
  type, public, extends(task_t):: particle_solver_t
    type(particles_t), pointer :: particles            ! internal particles
    type(particles_t), pointer :: exports              ! particles to export
    class(gpatch_t),   pointer :: patch => null()      ! pointer to patch
    type(link_t)   ,   pointer :: link => null()       ! task link
    type(dust_t)               :: dust                 ! dust solver
    real, pointer, contiguous  :: d(:,:,:)             ! density
    logical, dimension(3)      :: lower                ! has lower boundary
    logical, dimension(3)      :: upper                ! has upper boundary
    class(list_t), pointer     :: task_list            ! task list
    procedure(pforce), pointer, nopass :: pforce => null() ! point force
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure:: update
    procedure:: force
    procedure:: import_particles
    procedure:: is_outside
  end type
  integer:: verbose = 0
  real:: courant = 1.0
  real:: v_min = 0.1
  logical:: detailed_timer = .false.
CONTAINS

!===============================================================================
!> Initialize particle solver, reading parameters and adding task link
!===============================================================================
SUBROUTINE pre_init (self)
  class(particle_solver_t), target:: self
  integer:: iostat
  namelist /particle_solver_params/ verbose, courant, v_min, detailed_timer
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pre_init')
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, particle_solver_params, iostat=iostat)
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%output, particle_solver_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  self%dust%detailed_timer = detailed_timer
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
  class(particle_solver_t):: self
  class(gpatch_t), pointer:: p
  integer:: m(3)
  !-----------------------------------------------------------------------------
  ! Initialize particle lists
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init')
  allocate (self%particles, self%exports)       ! allocate particle lists
  call self%particles%init                      ! add particles
  call self%exports%init
  !-----------------------------------------------------------------------------
  ! Set flags indicating upper and lower boundaries.  The box size is self%box,
  ! and the patch position, size, and cell size are self%position, self%size,
  ! and self%ds, respectively.
  !-----------------------------------------------------------------------------
  self%upper = .false.
  self%lower = .false.
  p => self%patch
  if (p%position(3) - p%size(3)/2.0 - p%ds(3) < -p%box(3)/2.0) then
    self%lower(3) = .true.
  end if
  if (p%position(3) + p%size(3)/2.0 + p%ds(3) > +p%box(3)/2.0) then
    self%upper(3) = .true.
  end if
  self%dust%lower = self%lower
  self%dust%upper = self%upper
  write(io_unit%log,'(a,i6,2f8.3,4(2x,2l3))') &
    'particle_solver_t%init: id, pos, box, lower, upper =', &
    p%id, p%position(3), p%box(3), self%lower, self%upper
  !-----------------------------------------------------------------------------
  ! Initialize dust
  !-----------------------------------------------------------------------------
  self%dust%patch => self%patch                 ! make the patch known to it
  call self%dust%init (self%particles)          ! initialize dust solver
  m = self%patch%mesh%gn
  allocate (self%patch%force_per_unit_volume(m(1),m(2),m(3),3))
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize particle solver, after nbor lists are created, making it possible
!> to create a copy of the patch nbor list, with tasks pointing to self
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(particle_solver_t):: self
  class(list_t), pointer:: task_list
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init_task_list')
  self%task_list => task_list
  self%dust%task_list => task_list
  call self%patch%link%copy_nbor_list (self%patch%link%nbor, self%link%nbor)
  nbor => self%link%nbor
  do while (associated(nbor))
    nbor%task => self%link%task
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t):: self
  !.............................................................................
  real, dimension(3)         :: r, dr, v
  integer, dimension(3)      :: q, ir
  real   , dimension(3)      :: ods, t, p1
  class(particle_t), pointer :: p, next
  real(8)                    :: tmin, tmax, time(2)
  real                       :: dt, v_min_cell
  logical                    :: outside
  integer                    :: ip, np, no, nk, it, m(3), i(3), bin
  class(gpatch_t), pointer   :: patch
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%update', itimer=itimer)
  call self%import_particles
  !-----------------------------------------------------------------------------
  ods = 1.0/self%patch%ds
  np = self%particles%n
  time(1) = wallclock()
  patch => self%patch
  self%dust%d => patch%mem(:,:,:,patch%idx%d              ,patch%it,1)
  self%dust%p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  !-----------------------------------------------------------------------------
  ! Allocate array for number of particles per cell-bin
  !-----------------------------------------------------------------------------
  associate (dust=>self%dust)
  m = self%patch%gn
  allocate (dust%np(m(1),m(2),m(3),dust%n_bin))  ! work array for pp-cell-bin
  dust%np = 0.0                                  ! prepare to accumulate
  !-----------------------------------------------------------------------------
  !> Non-vectorized pusher
  !-----------------------------------------------------------------------------
  ip = 0
  no = 0
  nk = 0
  dt = huge(1.0)
  v_min_cell = v_min/minval(patch%ds)
call self%dust%check_particles ('before kdk')
  p => self%particles%head
  do while (associated(p))
    !---------------------------------------------------------------------------
    ! Get pointer to next, in case the particle is removed
    !---------------------------------------------------------------------------
    next => p%next
    !---------------------------------------------------------------------------
    ! Move the particle until it is ahead of the patch time
    !---------------------------------------------------------------------------
    do while (p%t <= patch%time)
      !-------------------------------------------------------------------------
      ! Courant conditions
      !-------------------------------------------------------------------------
      p%dt = courant/max(p%v(1),p%v(2),p%v(3),v_min_cell)
      dt = min(dt,p%dt)
      call kdk (self, p, ods)
      nk = nk+1
    end do
    !---------------------------------------------------------------------------
    ! The integer+real representation does not need to be updated until here
    !---------------------------------------------------------------------------
    q = floor(p%r)
    p%r = p%r - q
    p%q = p%q + q
    !---------------------------------------------------------------------------
    ! Monitor number of cells per cell-bin
    !---------------------------------------------------------------------------
    if(verbose > 1) then
      write(io_unit%log,*) 'update', p%id, p%q, p%t <= patch%time
      flush(io_unit%log)
    endif
    i = p%q                                     ! particle space index
    bin = p%bin                                 ! particle bin
    dust%np(i(1),i(2),i(3),bin) = &             ! accumulate
    dust%np(i(1),i(2),i(3),bin) + 1             ! into self%np
    !---------------------------------------------------------------------------
    ! Move the particle to the export list if it falls outside
    !---------------------------------------------------------------------------
    outside = self%is_outside (p)
    if (outside) then
      if(verbose > 1) then
        write (io_unit%log,*) 'remove', p%id, p%q
      endif
      call self%particles%remove (p)
      call self%exports%append (p)
      no = no+1
    end if
    ip = ip+1
    if (verbose > 2) then
      write (io_unit%log,'(a,i8,l3,1p,2e9.2,0p,3i4,1x,3f7.2,2(1x,1p,3e10.2))') &
      'update: id, outside, dt, c, q, r, v, a =', &
      p%id, outside, p%dt, p%c, p%q, p%r, p%v, p%a
    end if
    p => next
  end do
!call self%dust%check_particles ('after kdk')
  if (nk > 0) then
    self%dtime = dt
    self%particles%dtime = dt
  end if
  !-----------------------------------------------------------------------------
  if (verbose > 0) then
    time(2) = wallclock()-time(1)
    write (io_unit%log,'(a,i6,2x,a,1p,2e11.2,2x,a,3i8,2x,a,6p,f8.3)') &
      'particle_solver%update: id =', patch%id, &
      'dtime =', patch%dtime, self%particles%dtime, &
      'updates:', ip, no, nk, &
      'mus/particle =', time(2)/max(nk,1)
  end if
  !-----------------------------------------------------------------------------
  ! Do dust evolution after move, since we need the particle timestep values
  !-----------------------------------------------------------------------------
  call self%dust%update (self%particles)
  deallocate (dust%np)
  end associate
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Kick-Drift-Kick (KDK) update of particles with drag.  The standard KDK steps
!> need to be modified by computing the drag force at an approximate position
!> at time t+dt, then moving the particle back for the first kick, followed by
!> drift, and a 2nd kick, now with the drag force evaluated at t+2*dt
!===============================================================================
SUBROUTINE kdk (self, p, ods)
  class(particle_solver_t):: self
  type(particle_t):: p
  integer, save:: itimer=0
  real:: ods(3)
  !-----------------------------------------------------------------------------
  ! First kick, with drag at t+dt, for v(t) -> v(t+dt)
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('particle_solver_t%kdk',itimer=itimer)
  p%r = p%r + p%dt*p%v
  p%t = p%t + p%dt
  call self%dust%drag (p)
  p%r = p%r - p%dt*p%v
  p%v = p%v + p%dt*ods*(p%a+p%d)/(1.+p%c*p%dt*ods)
  !-----------------------------------------------------------------------------
  ! Drift, r -> r+2*dt*v(t+dt*a), with v -> vg(t+dt) if drag dominates
  !-----------------------------------------------------------------------------
  p%r = p%r + 2.*p%dt*p%v
  p%t = p%t + p%dt
  !-----------------------------------------------------------------------------
  ! Get external and drag acceleration at p(t+2*dt) and do the 2nd kick
  !-----------------------------------------------------------------------------
  call self%dust%drag (p)
  call self%force (p)
  p%v = p%v + p%dt*ods*(p%a+p%d)/(1.+p%c*p%dt*ods)
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE kdk

!===============================================================================
!> Add external force (acceleration) on a particle.
!===============================================================================
SUBROUTINE force (self, p)
  class(particle_solver_t):: self
  type(particle_t):: p
  class(mesh_t), pointer:: mesh(:)
  real(8):: r(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
#ifdef DBG
  if (detailed_timer) call trace%begin ('particle_solver_t%force',itimer=itimer)
#endif
  mesh => self%patch%mesh
  r = mesh%p + (p%q + p%r - mesh%o)*mesh%d
  p%a = self%pforce (r)
#ifdef DBG
  if (verbose > 4) then
    write (io_unit%log,'(a,i6,1p,4e11.2)') 'force: id, z, a =', p%id, r(3), p%a
  end if
  if (detailed_timer) call trace%end (itimer)
#endif
END SUBROUTINE force

!===============================================================================
!> Template
!===============================================================================
FUNCTION pforce (r) RESULT (f)
  real(8):: r(3), f(3)
END FUNCTION pforce

!===============================================================================
!> Import particles from nbors export lists
!===============================================================================
SUBROUTINE import_particles (self)
  class(particle_solver_t)     :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: task, nbtask
  class(gpatch_t)  , pointer :: nbpatch
  class(particle_t), pointer :: p, next
  real(8)                    :: r(3)
  integer                    :: q(3)
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  ! Loop over patch nbors
  !-----------------------------------------------------------------------------
  if (detailed_timer) &
    call trace%begin ('particle_solver_t%import_particles',itimer=itimer)
  nbor => self%patch%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (particle_solver_t)
      nbtask => nbor%task
      select type (nbtask)
      class is (particle_solver_t)
        nbpatch => nbtask%patch
        !-----------------------------------------------------------------------
        ! Loop over particles in export list and grab those that belong here
        !-----------------------------------------------------------------------
        p => nbtask%exports%head
        do while (associated(p))
          next => p%next
          r = (nbpatch%position - self%position)/nbpatch%ds + (p%r + p%q)
          q = floor(r)
          p%r = r - q
          p%q = p%q + q
          if (is_inside (self, p)) then
            call task%exports%remove (p)
            call self%particles%append (p)
          end if
          p => next
        end do
      end select
    end select
    nbor => nbor%next
  end do
  if (detailed_timer) call trace%end (itimer)
END SUBROUTINE import_particles

!===============================================================================
!> Check if the particle is outside the field patch
!===============================================================================
LOGICAL FUNCTION is_inside (self, particle)
  class(particle_solver_t):: self
  class(particle_t):: particle
  !.............................................................................
  is_inside = .not. self%is_outside (particle)
END FUNCTION is_inside

!===============================================================================
!> Check if the particle is outside the field patch
!===============================================================================
LOGICAL FUNCTION is_outside (self, particle)
  class(particle_solver_t):: self
  !.............................................................................
  class(particle_t):: particle
  class(mesh_t), pointer:: mesh(:)
  integer:: i
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('particle_solver_t%is_outside')
  mesh => self%patch%mesh
  is_outside = .false.
  do i=1,3
    !---------------------------------------------------------------------------
    ! If the particle is below index li then it is outside, unless this is a
    ! lower boundary
    !---------------------------------------------------------------------------
    if (particle%q(i) < mesh(i)%li) then
      if (.not. self%lower(i)) then
        is_outside = .true.
        exit
      end if
    end if
    !---------------------------------------------------------------------------
    ! If the particle is aboive index ui then it is outside, unless this is an
    ! upper boundary
    !---------------------------------------------------------------------------
    if (particle%q(i) > mesh(i)%ui) then
      if (.not. self%upper(i)) then
        is_outside = .true.
        exit
      end if
    end if
  end do
  if (detailed_timer) call trace%end ()
END FUNCTION is_outside

END MODULE particle_solver_mod
