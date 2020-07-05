!===============================================================================
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
    logical                    :: forced = .false.     ! external force
    logical, dimension(3)      :: lower                ! has lower boundary
    logical, dimension(3)      :: upper                ! has upper boundary
    integer, dimension(3)      :: lb=0, ub=0           ! lower and upper index bounds
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
    procedure:: sort
    procedure:: import
  end type
  integer:: verbose = 0
  integer:: ip_dbg = 0      ! particle ID to debug
  integer:: id_dbg = 0      ! patch ID to debug
  integer:: i_dbg(5) = 0    ! cell index triplet to debug
  real:: courant = 1.0
  logical:: detailed_timer = .false.
  logical:: do_sort = .true.
  logical:: on = .true.
CONTAINS

!===============================================================================
!> Initialize particle solver, reading parameters and adding task link
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(particle_solver_t), target:: self
  class(gpatch_t), pointer:: patch
  integer:: iostat
  namelist /particle_solver_params/ on, verbose, courant, detailed_timer, &
    do_sort, ip_dbg, id_dbg, i_dbg
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%pre_init')
  self%patch => patch                           ! make the patch known
  self%id = patch%id                            ! for now, inherit
  self%dust%patch => self%patch                 ! make the patch known
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, particle_solver_params, iostat=iostat)
    call io%namelist_warning ('particle_solver_params', iostat)
    write (io%namelist, particle_solver_params)
    first_time = .false.
  end if
  call io%lock%unset
  self%dust%detailed_timer = detailed_timer
  self%dust%ip_dbg = ip_dbg
  self%dust%i_dbg = i_dbg
  self%dust%do_sort = do_sort
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
  integer:: m(3)
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
  m = self%patch%mesh%gn
  allocate (self%particles%sorted(m(1),m(2),m(3),n_bin))
  call io%bits_mem (storage_size(self%particles%sorted), &
                   product(shape(self%particles%sorted)), 'sorted')
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
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Initialize particle solver, after patch is initialized
!===============================================================================
SUBROUTINE initial_condition (self)
  class(particle_solver_t):: self
  class(gpatch_t), pointer:: pa
  integer:: m(3), i
  !-----------------------------------------------------------------------------
  ! Set flags indicating upper and lower boundaries.  The box size is self%box,
  ! and the patch position, size, and cell size are self%position, self%size,
  ! and self%ds, respectively.
  !-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%initial_condition')
  call self%lock%set ('ICs')
  self%upper = .false.
  self%lower = .false.
  pa => self%patch
  do i=1,3
    self%lb = pa%mesh%li
    self%ub = pa%mesh%ui
    if (.not.pa%periodic(i)) then
      if (pa%position(i) - pa%size(i)/2.0 - pa%ds(i) < -pa%box(i)/2.0) then
        self%lower(i) = .true.
        self%lb(i) = pa%mesh(i)%lb
      end if
      if (pa%position(i) + pa%size(i)/2.0 + pa%ds(i) > +pa%box(i)/2.0) then
        self%upper(i) = .true.
        self%ub(i) = pa%mesh(i)%ub
      end if
    end if
  end do
  self%dust%lower = self%lower
  self%dust%upper = self%upper
  self%dust%lb = self%lb
  self%dust%ub = self%ub
  !-----------------------------------------------------------------------------
  ! Initialize dust.  At boundaries, add particles also in guard zones
  !-----------------------------------------------------------------------------
#ifdef __GFORTRAN__
  !$omp critical (alloc_cr)
  call self%dust%initial_condition
  !$omp end critical (alloc_cr)
#else
  call self%dust%initial_condition
#endif
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
  call self%sort
  call self%dust%consistency ('IC')
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
  call self%dust%init_task_list (task_list)
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
  real                       :: ods, dtds
  class(particle_t), pointer :: p, next
  real(8)                    :: tmin, tmax, t, time(2), self_time, wc, start
  real                       :: ds, dt, os, n_min
  logical                    :: outside, restart
  integer                    :: ip, np, no, nk, nm, it, m(3), i(3), ii(5), dir
  integer(kind=2)            :: i1, i2, i3, bin, l(3), u(3)
  class(gpatch_t), pointer   :: patch
  class(particles_t), pointer:: exports
  class(sorted_t), pointer   :: sorted
  integer, save              :: itimer=0, itimer1=0
  real    :: d, up, pc, pd(3), pa(3)
  real(8) :: t_prv(50), dt_prv(50), times(9), wc1, wc2, sum
  integer :: j, nj
  logical :: error
!-----------------------------------------------------------------------------
  if (.not. on) return
  call trace%begin ('particle_solver_t%update', itimer=itimer)
  start = wallclock()
  patch => self%patch
  if (self%patch%id == i_dbg(1)) then
    patch%dtime = max(patch%dtime,1d-10)
    self%dtime = max(self%dtime,1d-10)
    write (io_unit%log,'(a,2(1p,4e10.3,2x))') &
      'particle_solver_t%update begin: t, dt, t+dt =', &
      patch%time, patch%dtime, 1.0/patch%dtime, patch%time+patch%dtime, &
      self%time, self%dtime, 1.0/self%dtime, self%time+self%dtime
  end if
  call self%import
  !print *,wallclock()-start,'import'
  !call self%dust%consistency ('after import')
  !call self%sort
  !call self%dust%consistency ('after sort')
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
    write (io_unit%log,'(a,i6,1p,4e12.4)') &
      'particle_solver_t%update: id, time, dtime, min/max =', patch%id, patch%time, patch%dtime, &
      patch%fmaxval(patch%u_phase), patch%fminval(patch%u_phase)
  end if
  !-----------------------------------------------------------------------------
  ! Allocate00 a new exports list, while the previous one may still be in use
  !-----------------------------------------------------------------------------
  allocate (exports)
  call exports%init
  !call self%particles%lock%set ('update')
  time(1) = wallclock()
  !-----------------------------------------------------------------------------
  ods = 1.0/patch%ds(1)
  np = self%particles%n
  self%dust%d => patch%mem(:,:,:,patch%idx%d              ,patch%it,1)
  self%dust%p => patch%mem(:,:,:,patch%idx%px:patch%idx%pz,patch%it,1)
  !-----------------------------------------------------------------------------
  ! Allocate array for number of particles per cell-bin
  !-----------------------------------------------------------------------------
  associate (dust=>self%dust)
  m = patch%gn
  allocate (dust%vg(3,m(1),m(2),m(3)))
  call io%bits_mem (storage_size(dust%vg), &
                    product(shape(dust%vg)), 'vg')
!-------------------------------------------------------------------------------
  ! Cache the gas velocity
  !-----------------------------------------------------------------------------
  m = patch%gn
  do dir=1,3
    dust%vg(dir,:,:,:) = patch%mem(:,:,:,patch%idx%px+dir-1,patch%it,1) &
                        / (patch%ds(dir) &
                          *patch%mem(:,:,:,patch%idx%d       ,patch%it,1))
    end do
  !-----------------------------------------------------------------------------
  !> Non-vectorized pusher
  !-----------------------------------------------------------------------------
  ip = 0
  nk = 0
  no = 0
  nm = 0
  ds = minval(patch%ds, mask=patch%mesh%n > 1)
  t = huge(1d0)
  l = self%lb
  u = self%ub
  call dust%consistency ('before KDK loop')
  wc = wallclock()
  times = 0d0
  wc1 = wc
  do bin=1,n_bin
    os = dust%drag_coeff/dust%s_cgs(bin)
    do i3=l(3),u(3)
    do i2=l(2),u(2)
    do i1=l(1),u(1)
      if (verbose > 0) then
        wc = wallclock()
      end if
      sorted => self%particles%sorted(i1,i2,i3,bin)
      ip = ip + sorted%n
      d = patch%mem(i1,i2,i3,patch%idx%d,patch%it,1)
      up = patch%u_phase(i1,i2,i3)/ds
      v(1) = dust%vg(1,i1,i2,i3)
      v(2) = dust%vg(2,i1,i2,i3)
      v(3) = dust%vg(3,i1,i2,i3)
      pc = d*os
      p => sorted%head
      j = 1
      do while (j <= sorted%n .and. associated(p))
        next => p%next
        p%dt = dt
        !--- kick 1 ------------------------------------------------------------
        pd = pc*(v-p%v)                                   ! drag force
        pa = 0.0                                          ! external
        p%v = p%v + dt*(pa+pd)/(1.+pc*dt)                 ! stiff kick
        !--- drift -------------------------------------------------------------
        p%r = p%r + 2.*dt*p%v                             ! forward two dr
        p%t = p%t + 2.*dt                                 ! forward two dt
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
        p%q = p%q + q                                     ! add
        nk = nk+1                                         ! count
        if (sorted%n < self%dust%n_min) then
          call add_particles
        end if
        if (any(p%q < self%lb .or. p%q > self%ub)) then
          call self%particles%remove_sorted (p,i1,i2,i3,bin)
          self%dust%dd(i1,i2,i3,bin) = &
          self%dust%dd(i1,i2,i3,bin) - p%w
          call exports%append (p)
          no = no+1
        else if (any(p%q /= [i1,i2,i3])) then
          call self%particles%remove_sorted (p,i1,i2,i3,bin)
          call self%particles%add_sorted (p,p%q(1),p%q(2),p%q(3),bin)
          self%dust%dd(i1,i2,i3,bin) = &
          self%dust%dd(i1,i2,i3,bin) - p%w
          self%dust%dd(p%q(1),p%q(2),p%q(3),bin) = &
          self%dust%dd(p%q(1),p%q(2),p%q(3),bin) + p%w
          nm = nm+1
        else if (p%id == ip_dbg) then
          write (io_unit%log,'(a,3i6,i10,2x,5i4,3f6.2,1p,2(3e11.3,2x))') &
            '   move: ', &
            patch%id, exports%n, nk, p%id, p%q, bin, &
            self%particles%sorted(i1,i2,i3,bin)%n , &
            p%r, p%v, p%t, dt, self_time
          flush (io_unit%log)
        end if
        j = j+1
        p => next
      end do
      if (sorted%n > self%dust%n_max) then
        call remove_particles
      end if
    end do
    end do
    end do
  end do
  if (verbose > 0) then
    wc = wallclock()-wc
    write (io_unit%log,'(2i12,f6.2,6p,2f6.2," mus/p",0p,2f6.2)') &
    np, nk, wc, wc/np, wc/nk, real(nm)/nk, real(no)/nk
    call self%dust%consistency ('after KDK + export')
  end if
  !-----------------------------------------------------------------------------
  ! If updates were made to some particles, it should have been the ones with
  ! the smallest time, so the particle time should now increase
  !-----------------------------------------------------------------------------
  self%no = no
  if (nk > 0) then
    self%dtime = patch%dtime
    if (self%dtime <= 0d0) then
      write (stderr,'(a,2i6,1p,2(2(3e11.2,2x),2x))') &
      'particl_solver%update: ERROR', patch%id, nk, self%time, self%dtime, &
      self_time, patch%dtime
    end if
    if (self%patch%id == i_dbg(1)) then
      write (io_unit%log,'(a,2(1p,4e10.3,2x))') &
        'particle_solver_t%update dtime: t, dt, t+dt =', &
        patch%time, patch%dtime, 1./patch%dtime, patch%time+patch%dtime, &
        self%time, self%dtime, 1./self%dtime, self%time+self%dtime
    end if
    !---------------------------------------------------------------------------
    self%time = patch%time + patch%dtime
    self%dust%time = self%time
    self%dust%dtime = self%dtime
    self%particles%dtime = self%dtime
  end if
  !.............................................................................
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
  !-----------------------------------------------------------------------------
  ! Do dust evolution after move, since we need the particle timestep values
  !-----------------------------------------------------------------------------
  call dust%consistency ('before dust%update')
  call self%dust%update (self%particles)
  call dust%consistency ('after dust%update')
  if (self%patch%id == i_dbg(1)) then
    write (io_unit%log,'(a,2(1p,4e10.3,2x))') &
      'particle_solver_t%update  dust: t, dt, t+dt =', &
      patch%time, patch%dtime, 1./patch%dtime, patch%time+patch%dtime, &
      self%time-self%dtime, self%dtime, 1./self%dtime, self%time
  end if
  !-----------------------------------------------------------------------------
  ! Briefly lock self%xports and copy in particles from the scratch exports
  !-----------------------------------------------------------------------------
  call self%exports%lock%set ('copy')
  p => exports%head
  do while (associated(p))
    next => p%next
    call exports%remove (p)
    call self%exports%append (p)
    p => next
  end do
  call self%exports%lock%unset ('copy')
  !-----------------------------------------------------------------------------
  ! Dellocate the scratch exports list
  !-----------------------------------------------------------------------------
  call exports%dealloc
  deallocate (exports)
  call io%bits_mem (storage_size(dust%vg), &
                  -product(shape(dust%vg)), '-vg')
  deallocate (dust%vg)
  end associate
  self%time = self_time 
  if (self%patch%id == i_dbg(1)) then
    write (io_unit%log,'(a,2(1p,4e10.3,2x))') &
      'particle_solver_t%update   end: t, dt, t+dt =', &
      patch%time, patch%dtime, 1./patch%dtime, patch%time+patch%dtime, &
      self%time-self%dtime, self%dtime, 1./self%dtime, self%time
  end if
  call trace%end (itimer)
contains

!===============================================================================
!> top up particles per bin
!===============================================================================
subroutine add_particles
  real:: d_cgs, w
  real(8):: time
  integer:: j, n, pp_bin
  q = [i1,i2,i3]
  pp_bin = self%dust%pp_bin
  w = self%dust%dd(i1,i2,i3,bin)/pp_bin
  n = sorted%n
  p => sorted%head
  do j=1,n
    p%w = w
    p => p%next
  end do
  time = self%patch%time
  do j=1,pp_bin-n
    allocate (p)
    call io%bits_mem (storage_size(p), 1, 'p')
    call p%init
    p%q = q
    p%bin = bin
    p%r = self%dust%random%uniform (3)
    p%t = time
    p%w = w
    p%v = v
    call self%particles%add_sorted (p,i1,i2,i3,bin)
  end do
  !call self%dust%consistency2 ('add_p',i1,i2,i3,bin)
end subroutine add_particles

!===============================================================================
!> Bring particles per bin below max
!===============================================================================
subroutine remove_particles
  real:: d_cgs, w
  integer:: np, nr
  nr = 0
  p => sorted%head
  do while (sorted%n > self%dust%n_max)
    next => p%next
    if (p%flag == 0) then
      call self%particles%remove_sorted (p, i1, i2, i3, bin)
      call io%bits_mem (storage_size(p), -1, '-p')
      deallocate (p)
      nr = nr+1
    end if
    if (associated(p,sorted%tail)) then
      exit
    end if
    p => next
  end do
  np = sorted%n
  w = self%dust%dd(i1,i2,i3,bin)/np
  p => sorted%head
  do j=1,np
    p%w = w
    p => p%next
  end do
  !call self%dust%consistency2 ('rem_p',i1,i2,i3,bin)
end subroutine remove_particles
END SUBROUTINE update

!===============================================================================
!> Sort the particle list.  This has immediate benefits, without other source
!> code changes, since subsequent references to variables indexed by p%q will
!> often be in cache.  Further speedup will be possible, by explicitly looping
!> over the indices, treating the particles cell-by-cell.
!===============================================================================
SUBROUTINE sort (self)
  class(particle_solver_t):: self
  type(particle_t), pointer:: p, next, head, tail
  type(particles_t), pointer:: prev
  !type(particles_t), allocatable, target:: sorted(:,:,:)
  type(sorted_t), pointer:: particles
  integer:: i0, i(3), l(3), u(3), n, ii(5), n0, unit
  integer(kind=2):: i1, i2, i3, bin
  real(8):: start
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('particle_solver_t%sort', itimer=itimer)
  n0 = self%particles%n
  !print *, self%patch%id, n
  !.............................................................................
  if (verbose > 0) then
    unit = io_unit%log
    start = wallclock()
    write (unit,1) (wallclock()-start)/n0, 'allocate'
    1 format(9p,f9.1,2x,a)
    start = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Initialize
  !-----------------------------------------------------------------------------
  u = self%patch%gn
  do bin=1,n_bin; do i3=1,u(3); do i2=1,u(2); do i1=1,u(1)
    nullify (self%particles%sorted(i1,i2,i3,bin)%head)
    nullify (self%particles%sorted(i1,i2,i3,bin)%tail)
    self%particles%sorted(i1,i2,i3,bin)%n = 0
  end do; end do; end do; end do
  !.............................................................................
  if (verbose > 0) then
    write (unit,1) (wallclock()-start)/n0, 'initialize'
    start = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Run through the particle list and store particles in boxes
  !-----------------------------------------------------------------------------
  p => self%particles%head
  do while (associated(p))
    next => p%next
    nullify (p%next)
    i = p%q
    if (any(i < self%lb .or. i > self%ub)) then
      print *,'error in q:', p%id, p%q, p%bin
    end if
    particles => self%particles%sorted(i(1),i(2),i(3),p%bin)
    p%prev => particles%tail
    !---------------------------------------------------------------------------
    ! If the box is empty, point head and tail to the particle and increment.
    ! If not empty, point to the particle with tail%next, and repoint tail.
    !---------------------------------------------------------------------------
    if (particles%n == 0) then
      particles%head => p
    else
      particles%tail%next => p
    end if
    particles%tail => p
    particles%n = particles%n+1
    p => next
  end do
  !.............................................................................
  if (verbose > 0) then
    write (unit,1) (wallclock()-start)/n0, 'collect'
    start = wallclock()
  end if
  !-----------------------------------------------------------------------------
  ! Chain together the self%particles%sorted particle boxes, with tail always
  ! pointing to the last known particle
  !-----------------------------------------------------------------------------
  n = 0
  nullify (head,tail)
  l = self%patch%li
  u = self%patch%ui
  do bin=1,n_bin
  do i3=l(3),u(3)
    do i2=l(2),u(2)
      do i1=l(1),u(1)
        particles => self%particles%sorted(i1,i2,i3,bin)
        !-----------------------------------------------------------------------
        ! If there are particles in a box, connect heads and tails
        !-----------------------------------------------------------------------
        if (particles%n > 0) then
          if (.not.associated(head)) then ! if head not assigned yet
            head => particles%head        ! point to the box head
          end if
          if (associated(tail)) then      ! if tail not assigned yet
            tail%next => particles%head   ! connect the tail to the head
            tail%next%prev => tail        ! make the double link back
          end if
          tail => particles%tail          ! last known particle so far
          n = n + particles%n             ! add up
        end if
      end do
    end do
  end do
  end do
  if (n /= n0) then
    print *, 'ERROR: wrong numbers (1) n0, n =', n0, n
  end if
  !-----------------------------------------------------------------------------
  ! Let the initial particle list begin with the first box
  !-----------------------------------------------------------------------------
  self%particles%head => head
  self%particles%tail => tail
  self%particles%n = n
  !-----------------------------------------------------------------------------
  ! If verbose, validate the sort
  !-----------------------------------------------------------------------------
  if (verbose > 0) then
    write (unit,1) (wallclock()-start)/n0, 'tie together'
    start = wallclock()
    p => self%particles%head
    n = 0
    do bin=1,n_bin
    do i3=l(3),u(3)
      do i2=l(2),u(2)
        do i1=l(1),u(1)
          particles => self%particles%sorted(i1,i2,i3,bin)
          do i0=1,particles%n
            if (.not.associated(p)) then
              print *,i1,i2,i3,' missing particle'
              write (io_unit%log,*) i1,i2,i3,' missing particle'
              exit
            end if
            if (any(p%q /= [i1,i2,i3])) print *, 'sort error'
            n = n+1
            p => p%next
          end do
        end do
      end do
    end do
    end do
    if (n /= n0) then
      print *, 'ERROR: wrong numbers (2) n0, n =', n0, n
    end if
    write (unit,1) (wallclock()-start)/n0, 'validate'
  end if
  !.............................................................................
  call trace%end (itimer)
END SUBROUTINE sort

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
  class(particle_solver_t)     :: self
  !.............................................................................
  class(nbor_t)    , pointer :: nbor
  class(task_t)    , pointer :: nbtask
  class(particle_t), pointer :: p, next
  real(8)                    :: r(3), dist(3)
  integer                    :: ni, ii(5)
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
      call nbtask%exports%lock%set ('import')
      dist = nbtask%patch%task_t%distance(self%patch)/nbtask%patch%ds
      !-------------------------------------------------------------------------
      ! Loop over particles in export list and grab those that belong here
      !-------------------------------------------------------------------------
      ni = 0
      p => nbtask%exports%head
      do while (associated(p))
        next => p%next
        r = dist + (p%r+p%q)
        q = floor(r)
        inside = all(q >= self%lb .and. q <= self%ub)
        if (inside) then
          p%r = r-q
          p%q = q
          call nbtask%exports%remove (p)
          call self%dust%consistency2 ('before add',p%q(1),p%q(2),p%q(3),p%bin)
          call self%particles%add_sorted (p, q(1), q(2), q(3), p%bin)
          call self%dust%consistency2 ('after add',p%q(1),p%q(2),p%q(3),p%bin)
          self%dust%dd(q(1),q(2),q(3),p%bin) = &
          self%dust%dd(q(1),q(2),q(3),p%bin) + p%w
          self%time = min(self%time, p%t)
          ni = ni + 1
          ii(1) = self%patch%id; ii(2:4) = p%q; ii(5) = p%bin
          if (verbose > 0 .or. p%id == ip_dbg .or. all(ii == i_dbg)) then
            write (io_unit%log,'(a,3i6,i10,2x3i4,l4)') &
             'import: id, nbor, ni, pid, q =', &
              self%patch%id, nbtask%patch%id, ni, p%id, q
          end if
        end if
        p => next
      end do
      self%ni = self%ni + ni
      if ((ni > 0 .and. verbose > 1) .or. verbose > 3) then
        write (io_unit%log,'(a,4i6,2x,3i4)') 'imported: id, nbor, ni, offset =', &
          self%patch%id, nbtask%patch%id, nbtask%exports%n, ni, &
          nint(dist)
      end if
      call nbtask%exports%lock%unset ('import')
    end select
    nbor => nbor%next
  end do
  if ((self%ni > 0 .and. verbose > 1) .or. verbose > 2) then
    write (io_unit%log,'(a,2i6)') 'imported: id, self%ni =', &
      self%patch%id, self%ni
    flush (io_unit%log)
  end if
  call trace%end (itimer)
END SUBROUTINE import

END MODULE particle_solver_mod
