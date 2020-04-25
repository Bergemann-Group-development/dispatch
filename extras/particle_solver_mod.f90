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
  USE omp_lock_mod
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
    class(nbor_t),     pointer :: nbor => null()       ! nbor list
    type(link_t)   ,   pointer :: link => null()       ! task link
    type(dust_t)               :: dust                 ! dust solver
    real, pointer, contiguous  :: d(:,:,:)             ! density
    logical                    :: forced = .false.     ! external force
    logical, dimension(3)      :: lower                ! has lower boundary
    logical, dimension(3)      :: upper                ! has upper boundary
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
    procedure:: force
    procedure:: import_particles
    procedure:: is_outside
  end type
  integer:: verbose = 0
  integer:: id_dbg = -1
  real:: courant = 1.0
  logical:: detailed_timer = .false.
CONTAINS

!===============================================================================
!> Initialize particle solver, reading parameters and adding task link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(particle_solver_t), target:: self
  class(gpatch_t), pointer:: patch
  integer:: iostat
  namelist /particle_solver_params/ verbose, courant, detailed_timer, &
    id_dbg
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pre_init')
  self%patch => patch                           ! make the patch known
  self%id = patch%id                            ! for now, inherit
  self%dust%patch => self%patch                 ! make the patch known
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
  self%dust%id_dbg = id_dbg
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
  class(gpatch_t), pointer:: pa
  integer:: m(3), l(3), u(3), i
  !-----------------------------------------------------------------------------
  ! Initialize particle lists
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%init')
  allocate (self%particles, self%exports)       ! allocate particle lists
  call self%particles%init                      ! initialize particle list
  call self%exports%init                        ! initialize export list
  call self%dust%init (self%particles)          ! initialize dust solver
  !-----------------------------------------------------------------------------
  ! Allocate force array
  !-----------------------------------------------------------------------------
  if (.not.allocated(self%patch%force_per_unit_volume)) then
    m = self%patch%mesh%gn
    allocate (self%patch%force_per_unit_volume(m(1),m(2),m(3),3))
  end if
  if (.not.associated(self%lock)) then
    allocate (self%lock)
    call self%lock%init ('part')
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particle_solver_t):: self
  class(gpatch_t), pointer:: pa
  integer:: m(3), l(3), u(3), i
  !-----------------------------------------------------------------------------
  ! Set flags indicating upper and lower boundaries.  The box size is self%box,
  ! and the patch position, size, and cell size are self%position, self%size,
  ! and self%ds, respectively.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%initial_condition')
  self%upper = .false.
  self%lower = .false.
  pa => self%patch
  do i=1,3
    if (.not.pa%periodic(i)) then
      if (pa%position(i) - pa%size(i)/2.0 - pa%ds(i) < -pa%box(i)/2.0) then
        self%lower(i) = .true.
      end if
      if (pa%position(i) + pa%size(i)/2.0 + pa%ds(i) > +pa%box(i)/2.0) then
        self%upper(i) = .true.
      end if
    end if
  end do
  self%dust%lower = self%lower
  self%dust%upper = self%upper
  !-----------------------------------------------------------------------------
  ! Initialize dust.  At boundaries, add particles also in guard zones
  !-----------------------------------------------------------------------------
  l = pa%mesh%li
  u = pa%mesh%ui
  where (self%lower) l = pa%mesh%lb
  where (self%upper) u = pa%mesh%ub
  call self%dust%initial_condition (l, u)
  !.............................................................................
  if (verbose > 0) then
    if (verbose > 1) then
      write (io%output,1) &
        'init: id, box, pos, size, lower, upper, n =', &
        pa%id, pa%box, pa%position, pa%size, &
        self%lower, self%upper, &
        real(self%particles%n)
      1 format(a,i4,3(1x,3f7.3),1x,2(1x,3l1),1p,e11.2)
    end if
    write (io_unit%log,1) &
      'init: id, box, pos, size, lower, upper, n =', &
      pa%id, pa%box, pa%position, pa%size, &
      self%lower, self%upper, &
      real(self%particles%n)
  else
    write(io_unit%log,'(a,i6,2(3f8.3,2x),1x,2(1x,3l1),1p,e11.2)') &
     'particle_solver_t%init: id, pos, box, lower, upper, n =', &
      pa%id, pa%position, pa%box, self%lower, self%upper, real(self%particles%n)
  end if
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
  call trace%begin ('particle_solver_t%init_task_list')
  self%task_list => task_list
  self%dust%task_list => task_list
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Update particle positions
!===============================================================================
SUBROUTINE update (self)
  class(particle_solver_t):: self
  !.............................................................................
  real, dimension(3)         :: r, dr, v
  integer, dimension(3)      :: q, ir, l, u
  real                       :: ods, dtds
  class(particle_t), pointer :: p, next
  real(8)                    :: tmin, tmax, t, time(2), patch_time
  real                       :: ds, dt
  logical                    :: outside
  integer                    :: ip, np, no, nk, it, m(3), i(3), bin
  class(gpatch_t), pointer   :: patch
  class(particles_t), pointer:: exports
  integer, save              :: itimer=0, itimer1=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%update', itimer=itimer)
  call self%import_particles
  self%forced = associated(self%pforce)
  allocate (exports)
  call exports%init
  call self%particles%lock%set ('update')
  !-----------------------------------------------------------------------------
  ! Only evolve particles if the last particle is behind the patch%time (which
  ! has not been updated yet, so we add the patch timestep)
  !-----------------------------------------------------------------------------
  patch => self%patch
  patch_time = patch%time
  l = patch%mesh%li
  u = patch%mesh%ui
  self%v_min = patch%fminval(patch%u_phase)/self%patch%ds(1)
  self%v_max = patch%fmaxval(patch%u_phase)/self%patch%ds(1)
  if (verbose > 1) then
    write (io_unit%log,'(a,i6,1p,4e12.4)') &
      'particle_solver_t%update: id, time, dtime, min/max =', patch%id, patch%time, patch%dtime, &
      patch%fmaxval(patch%u_phase), patch%fminval(patch%u_phase)
  end if
  if (self%time <= patch_time) then
    time(1) = wallclock()
    !---------------------------------------------------------------------------
    ods = 1.0/patch%ds(1)
    np = self%particles%n
    self%dust%d => patch%mem(:,:,:,patch%idx%d              ,patch%it,1)
    self%dust%p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
    !---------------------------------------------------------------------------
    ! Allocate array for number of particles per cell-bin
    !---------------------------------------------------------------------------
    associate (dust=>self%dust)
    !---------------------------------------------------------------------------
    !> Non-vectorized pusher
    !---------------------------------------------------------------------------
    ip = 0
    nk = 0
    no = 0
    ds = minval(patch%ds, mask=patch%mesh%n > 1)
!call self%dust%check_particles ('before kdk')
    t = huge(1d0)
    p => self%particles%head
    do while (associated(p))
      if (detailed_timer) &
        call trace%begin ('particle_solver_t%move_particle', itimer=itimer1)
      if (verbose > 3) then
        write (io_unit%log,1) &
        'update: id, outside, dt, c, q, r, v    =', &
        patch%id, p%id, .false., p%dt, p%c, p%q, p%r, p%v
      end if
      !-------------------------------------------------------------------------
      ! Get pointer to next, in case the particle is removed
      !-------------------------------------------------------------------------
      next => p%next
      !-------------------------------------------------------------------------
      ! Move the particle until it is ahead of the patch time
      !-------------------------------------------------------------------------
      outside = .false.
      i = p%q
      if (self%is_outside(i)) then
        outside = .true.
        print '(a,i12,2x,3i4,2x,3f6.3,1p,e10.2)', &
          'particle_solver_t%update ERROR: id, q, r, u_ph =', &
          p%id, p%q, p%r, patch%u_phase(i(1),i(2),i(3))/ds
        write (io_unit%log,'(a,i12,2x,3i4,2x,3f6.3,1p,e10.2)') & 
          'particle_solver_t%update ERROR: id, q, r, u_ph =', &
          p%id, p%q, p%r, patch%u_phase(i(1),i(2),i(3))/ds
      else if (p%t <= patch_time) then
        !call self%force (p)
        !-------------------------------------------------------------------------
        ! Courant conditions; allow at most courant/patch%courant larger step
        !-------------------------------------------------------------------------
        p%dt = courant*0.5/max(patch%u_phase(i(1),i(2),i(3))/ds, &
                               abs(p%v(1)),abs(p%v(2)),abs(p%v(3)))
        dtds = p%dt*ods
        do while (p%t <= patch_time)
          call kdk (self, p, dtds)
          nk = nk+1
          !---------------------------------------------------------------------
          ! The integer+real representation may need to be updated here, since
          ! the particle might move over a cell boundary
          !---------------------------------------------------------------------
          q = floor(p%r)
          p%r = p%r - q
          p%q = p%q + q
          !---------------------------------------------------------------------
          ! Check if it has gone outside limits
          !---------------------------------------------------------------------
          i = p%q
          outside = self%is_outside (i)
          if (outside) then
            exit
          end if
        end do
      end if
      !-------------------------------------------------------------------------
      ! Move the particle to the export list if it falls outside
      !-------------------------------------------------------------------------
      if (outside) then
        call self%particles%remove (p)
        call exports%append (p)
        no = no+1
        if(verbose > 1 .or. p%id == id_dbg) then
          write (io_unit%log,'(a,3i6,i10,2x,3i4)') &
            'export: id, no, pid, q =', &
            patch%id, self%exports%n, no, p%id, p%q
        endif
      else
        !-----------------------------------------------------------------------
        ! t is the smallest particle time of the set
        !-----------------------------------------------------------------------
        t = min(t,p%t)
        ip = ip+1
        if (verbose > 2 .or. p%id == id_dbg) then
          write (io_unit%log,1) &
          'update: id,out,dt,c,q,r,v,a=', &
          patch%id, p%id, outside, p%dt, p%c, p%q, p%r, p%v, p%a, &
          patch%u_phase(i(1),i(2),i(3))
          1 format(a,i5,i8,l3,1p,2e9.2,0p,3i4,1x,3f7.2,3(1x,1p,3e9.1))
        end if
      end if
      p => next
      if (detailed_timer) call trace%end (itimer1)
    end do
    !call self%dust%check_particles ('after kdk')
    !---------------------------------------------------------------------------
    ! If updates were made to some particles, it should have been the ones with
    ! the smallest time, so the particle time should now increase
    !---------------------------------------------------------------------------
    self%no = no
    if (nk > 0) then
      self%dtime = t - self%time
      if (self%dtime <= 0d0) then
        write (stderr,'(a,2i6,1p,2(2(3e11.2,2x),2x))') &
        'particl_solver%update: ERROR', patch%id, nk, self%time, self%dtime, &
        patch_time, patch%dtime
      end if
      self%time = self%time + self%dtime
      self%dust%time = self%time
      self%dust%dtime = self%dtime
      self%particles%dtime = self%dtime
    end if
    !...........................................................................
    if (verbose > 0) then
      time(2) = wallclock()-time(1)
      write (io_unit%log, &
        '(a,i6,2x,2(a,1p,2e12.4,2x),a,2i8,2i6,2x,a,6p,f8.3,2x,a,3p,f6.1)') &
        'particle_solver_t%update: id:', patch%id, &
        't:', patch%time, self%time, &
        'dt:', patch%dtime, self%dtime, &
        'updates:', ip, nk, self%ni, self%no, &
        'mus/particle:', time(2)/max(nk,1), 'ms:', time(2)
    end if
    !---------------------------------------------------------------------------
    ! Do dust evolution after move, since we need the particle timestep values
    !---------------------------------------------------------------------------
    call self%dust%update (self%particles)
    end associate
  else if (verbose > 0) then
    write (io_unit%log,'(a,2f12.6,i10)') &
      'particle_solver_t%update: self%time, patch_time, n =', &
      self%time, patch_time, self%particles%n
  end if
  call self%lock%set ('update')
  call omp_lock%destroy (self%exports%lock)
  deallocate (self%exports)
  self%exports => exports
  call self%lock%unset ('update')
  call self%particles%lock%unset ('update')
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Kick-Drift-Kick (KDK) update of particles with drag.  The standard KDK steps
!> need to be modified by computing the drag force at an approximate position
!> at time t+dt, then moving the particle back for the first kick, followed by
!> drift, and a 2nd kick, now with the drag force evaluated at t+2*dt. Since
!> velocity is measured in cells per unit time, the acceleration is multiploied
!> by the inverse cell size, ods.
!===============================================================================
SUBROUTINE kdk (self, p, dtds)
  class(particle_solver_t):: self
  type(particle_t):: p
  real:: dtds
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! First kick, with drag at t+dt, for v(t) -> v(t+dt)
  !-----------------------------------------------------------------------------
  if (detailed_timer) call trace%begin ('particle_solver_t%kdk',itimer=itimer)
  p%r = p%r + p%dt*p%v
  p%t = p%t + p%dt
  call self%dust%drag (p)
  p%r = p%r - p%dt*p%v
  p%v = p%v + dtds*(p%a+p%d)/(1.+p%c*dtds)
  if (verbose > 3 .or. p%id == id_dbg) then
    write (io_unit%log,1) &
      '  kdk1: id,out,dt,c,q,r,v,a=', &
      self%patch%id, p%id, .false., p%dt, p%c, p%q, p%r, p%v, p%a, p%d, &
      self%patch%u_phase(p%q(1),p%q(2),p%q(3))/self%patch%ds(1), &
      self%v_min, self%v_max
    1 format(a,i5,i8,l3,1p,2e9.2,0p,3i3,3f6.2,4(1x,3f7.3))
  end if
!-----------------------------------------------------------------------------
  ! Drift, r -> r+2*dt*v(t+dt*a), with v -> vg(t+dt) if drag dominates
  !-----------------------------------------------------------------------------
  p%r = p%r + 2.*p%dt*p%v
  p%t = p%t + p%dt
  !-----------------------------------------------------------------------------
  ! Get external and drag acceleration at p(t+2*dt) and do the 2nd kick
  !-----------------------------------------------------------------------------
  call self%dust%drag (p)
  if (self%forced) then
    call self%force (p)
  end if
  p%v = p%v + dtds*(p%a+p%d)/(1.+p%c*dtds)
  if (verbose > 3 .or. p%id == id_dbg) then
    write (io_unit%log,1) &
      '  kdk2: id,out,dt,c,q,r,v,a=', &
      self%patch%id, p%id, .false., p%dt, p%c, p%q, p%r, p%v, p%a, p%d, &
      self%patch%u_phase(p%q(1),p%q(2),p%q(3))/self%patch%ds(1), &
      self%v_min, self%v_max
    flush (io_unit%log)
   end if
!-----------------------------------------------------------------------------
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
  if (self%forced) then
    mesh => self%patch%mesh
    r = mesh%p + (p%q + p%r - mesh%o)*mesh%d
    p%a = self%pforce (r)
  else
    p%a = 0.0
  end if
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
  f = 0d0
END FUNCTION pforce

!===============================================================================
!> Import particles from nbors export lists
!===============================================================================
SUBROUTINE import_particles (self)
  class(particle_solver_t)     :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: nbtask
  class(particle_t), pointer :: p, next
  real(8)                    :: r(3), dist(3)
  integer                    :: q(3), ni
  logical                    :: inside
  integer, save              :: itimer=0
  !-----------------------------------------------------------------------------
  ! Loop over patch nbors
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%import_particles',itimer=itimer)
  call self%particles%lock%set ('import')
  self%ni = 0
  nbor => self%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (particle_solver_t)
      call nbtask%lock%set ('import')
      dist = nbtask%patch%distance(self%patch)/nbtask%patch%ds
      !-------------------------------------------------------------------------
      ! Loop over particles in export list and grab those that belong here
      !-------------------------------------------------------------------------
      ni = 0
      p => nbtask%exports%head
      do while (associated(p))
        next => p%next
        r = dist + (p%r+p%q)
        q = floor(r)
        inside = is_inside (self, q)
        if (inside) then
          p%r = r-q
          p%q = q
          call nbtask%exports%remove (p)
          call self%particles%append (p)
          self%time = min(self%time, p%t)
          ni = ni + 1
          if (verbose > 1) then
            write (io_unit%log,'(a,3i6,i10,2x3i4)') &
             'import: id, nbor, ni, pid, q =', &
              self%patch%id, nbtask%patch%id, ni, p%id, q
          end if
        end if
        if (p%id == id_dbg) then
          write (io_unit%log,'(a,3i6,i10,2x3i4,l4)') &
           'import: id, nbor, ni, pid, q =', &
            self%patch%id, nbtask%patch%id, ni, p%id, q, inside
        end if
        p => next
      end do
      self%ni = self%ni + ni
      if ((ni > 0 .and. verbose > 1) .or. verbose > 3) then
        write (io_unit%log,'(a,4i6,2x,3i4)') 'imported: id, nbor, ni, offset =', &
          self%patch%id, nbtask%patch%id, nbtask%exports%n, ni, &
          nint(dist)
      end if
      call nbtask%lock%unset ('import')
    end select
    nbor => nbor%next
  end do
  if ((self%ni > 0 .and. verbose > 1) .or. verbose > 2) then
    write (io_unit%log,'(a,2i6)') 'imported: id, self%ni =', &
      self%patch%id, self%ni
    flush (io_unit%log)
  end if
  call self%particles%lock%unset ('import')
  call trace%end (itimer)
END SUBROUTINE import_particles

!===============================================================================
!> Check if the particle is outside the field patch
!===============================================================================
LOGICAL FUNCTION is_inside (self, q)
  class(particle_solver_t):: self
  integer:: q(3)
  !.............................................................................
  is_inside = .not. is_outside(self, q)
END FUNCTION is_inside

!===============================================================================
!> Check if the particle is outside the field patch
!===============================================================================
LOGICAL FUNCTION is_outside (self, q)
  class(particle_solver_t):: self
  integer:: q(3)
  !.............................................................................
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
    if (q(i) < mesh(i)%li) then
      if (.not. self%lower(i)) then
        is_outside = .true.
        exit
      end if
    end if
    !---------------------------------------------------------------------------
    ! If the particle is aboive index ui then it is outside, unless this is an
    ! upper boundary
    !---------------------------------------------------------------------------
    if (q(i) > mesh(i)%ui) then
      if (.not. self%upper(i)) then
        is_outside = .true.
        exit
      end if
    end if
  end do
  if (detailed_timer) call trace%end ()
END FUNCTION is_outside

END MODULE particle_solver_mod
