!===============================================================================
!> A sink particle is represented as primarily a particle task, whose position
!> is updated with a particle solver. The patch is usually much smaller than an
!> MHD patch, and its values are downloaded from contributing MHD patches. 
!>
!> A particular quality level indicates to the download_mod that all cells from
!> the sink particle patch should be downloaded to overlapping patches, thus
!> communicating back the removal of accreted mass and momentum.
!>
!> The procedure is called from refine_t%extras%check_refine which handles
!> creating a task list link, adding an nbor list, and adding the task to the
!> task list.
!===============================================================================
MODULE sink_patch_mod
  USE io_mod
  USE trace_mod
  USE gpatch_mod
  USE index_mod
  USE download_mod
  USE data_hub_mod
  USE link_mod
  USE list_mod
  USE bits_mod
  USE scaling_mod
  USE kinds_mod
  USE task_mod
  USE dll_mod
  USE particle_mod
  USE particle_solver_mod
  implicit none
  private
  type, public, extends (gpatch_t):: sink_patch_t
    real(8):: mass=0d0, mom(3)
    real(8):: acceleration(3)=0d0
    real(8), dimension(:,:)  , allocatable:: p, v, a
    real(4), dimension(:,:,:), allocatable:: r
    class(particle_solver_t), pointer:: particle_solver
    class(dll_node_t), pointer:: particle
    logical:: on
  contains
    procedure:: init
    procedure:: dealloc
    procedure:: init_task_list
    procedure:: check
    procedure:: init_task
    procedure:: dnload
    procedure:: update
    procedure:: accrete
    procedure:: move
    procedure:: check_refine
    procedure, private:: refine_check
  end type
  !-----------------------------------------------------------------------------
  logical:: on=.false.                ! turn feature on/off
  logical:: do_translate=.true.       ! make a transformation to the velocity of the sink
  logical:: force_refine=.true.       ! insist on refining down to levelmax at star
  logical:: only_high_levels=.true.   ! only accrete at high levels
  integer:: id=0                      ! sink id
  integer:: center_star=0             ! index of star to zoom in to
  integer:: verbose=2                 ! verbosity level
  integer:: max_sinks=10000           ! for testing
  real:: accretion_rate=0.01          ! The rate of accretion, measured in Omega_Kepler
  real:: accretion_efficiency=0.5     ! The accretion radius, measured in cells
  real:: accretion_fraction=1e-6      ! The accretion radius, measured in cells
  real:: accretion_radius=4.          ! The accretion radius, measured in cells
  real:: exclusion_radius=32.         ! The radius, measured in cells, inside which no other sink are born
  real:: jeans_resolution=2.          ! Minimum, umber of cells per Jeans' length
  real:: rho_limit=1e6                ! The limit above which sink_particles are created
  real:: rho_limit_factor=8.          ! Express rho_limit relative to Jeans' ladder top
  real:: rho_fractionR=1e-6           ! The fration of rho_limit above which accretion happens
  real:: rho_cut=1.0                  ! ??
  real:: softening_length=1.0         ! Softening of gravitational potential
  real(8):: out_time=0.1              ! cadence of output into sink.dat
  character(len=64):: file ='sink.dat'! name of sink data file
  !-----------------------------------------------------------------------------
  integer, save:: n_sinks=0           ! counter
  type(sink_patch_t), public:: sink_patch
CONTAINS

!===============================================================================
!> Initialize an instance of the sink_patch_t data type, with parameters for
!> sink creation read from a sink_patch_params namelist.
!===============================================================================
SUBROUTINE init (self)
  class(sink_patch_t):: self
  !.............................................................................
  logical, save:: first_time=.true.
  integer:: iostat
  namelist /sink_patch_params/ &
    on, do_translate, force_refine, only_high_levels, center_star, verbose, &
    accretion_rate, accretion_efficiency, accretion_fraction, accretion_radius, &
    exclusion_radius, jeans_resolution, rho_limit, rho_limit_factor, rho_cut, &
    max_sinks, softening_length, out_time, file
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, sink_patch_params, iostat=iostat)
    if (io_unit%master) write (io_unit%output, sink_patch_params)
  end if
  !$omp end critical (input_cr)
  self%on = on
  self%particle_solver => particle_solver
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Deallocate permanent allocatables
!===============================================================================
SUBROUTINE dealloc (self)
  class(sink_patch_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%dealloc')
  if (allocated(self%r)) then
    call io%bits_mem (-storage_size(self%r), product(shape(self%r)), '-sink%r')
    deallocate (self%r)
  end if
  if (allocated(self%p)) then
    call io%bits_mem (-storage_size(self%p), product(shape(self%p)), '-sink%p')
    deallocate (self%p)
  end if
  if (allocated(self%v)) then
    call io%bits_mem (-storage_size(self%v), product(shape(self%v)), '-sink%v')
    deallocate (self%v)
  end if
  if (allocated(self%a)) then
    call io%bits_mem (-storage_size(self%a), product(shape(self%a)), '-sink%a')
    deallocate (self%a)
  end if
  call trace%end()
END SUBROUTINE dealloc

!===============================================================================
!> Allows storing a link to a task list in the sink_patch_t data type
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(sink_patch_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  self%task_list => task_list
END SUBROUTINE init_task_list

!===============================================================================
!> Test if a sink particle should be created in a given patch
!===============================================================================
FUNCTION check (self, patch)
  class(sink_patch_t):: self
  class(gpatch_t), target:: patch
  logical:: check
  !-----------------------------------------------------------------------------
  check = .false.
  if (n_sinks >= max_sinks) return
  call trace%begin ('sink_patch_t%check')
  check = .true.
  check = check .and. patch%fmaxval(patch%idx%d) > rho_limit
  if (check) then
    print *, 'sink_patch_t%check: dmax in', &
      maxloc (patch%mem(:,:,:,patch%idx%d,patch%it,1))
    !$omp atomic
    n_sinks = n_sinks+1
    print *,'n_sinks =', n_sinks
  end if
  call trace%end()
END FUNCTION check

!===============================================================================
!> Create a sink partícle 
!===============================================================================
SUBROUTINE init_task (self, patch)
  class(sink_patch_t):: self
  class(gpatch_t), target:: patch
  !.............................................................................
  class(link_t), pointer:: nbor, prev
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  real(8), dimension(:), pointer:: r
  integer:: i, l(3), u(3), m(3), ix, iy, iz
  real:: r2, dmin
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%init_task')
  allocate (self%p(3,self%nt))
  allocate (self%v(3,self%nt))
  allocate (self%a(3,self%nt))
  call io%bits_mem (storage_size(self%p), product(shape(self%p)), 'sink%p')
  call io%bits_mem (storage_size(self%v), product(shape(self%v)), 'sink%v')
  call io%bits_mem (storage_size(self%a), product(shape(self%a)), 'sink%a')
  l = patch%mesh%li
  u = patch%mesh%ui
  d => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)
  m = maxloc(d)+l-1
  do i=1,3
    r => patch%mesh(i)%r
    self%position(i) = patch%position(i) + r(m(i))
    self%ds(i) = patch%mesh(i)%d
  end do
  self%n  = 2*ceiling(accretion_radius)
  self%ng = 2
  self%nv = patch%nv
  self%nt = patch%nt
  self%nw = 1
  !-----------------------------------------------------------------------------
  ! We should use no_mans_land=.false., and an odd number for self%n
  !-----------------------------------------------------------------------------
  self%no_mans_land = .false.
  self%n  = (self%n/2)*2 + 1
  self%size = self%ds*(self%n-1)
  !-----------------------------------------------------------------------------
  ! Set a high level on the sink particle patch, to force downloading of interior
  ! data
  !-----------------------------------------------------------------------------
  call self%patch_t%setup (self%position, self%size, self%n, self%ng, self%nt, self%nv, self%nw)
  if (self%id==io%id_debug) then
    print *, 'sink_patch_t%init_task: id, patch%id, m, dmax', self%id, patch%id, m, maxval(d)
    print *, 'sink_patch_t%init_task: (pos-patch%pos)/size', (self%position-patch%position)/patch%size
    print *, 'sink_patch_t%init_task: patch%pos, d', patch%myposition(m), patch%mem(m(1),m(2),m(3),patch%idx%d,patch%it,1)
    print *, 'sink_patch_t%init_task:  self%n   , patch%n   ', self%n   , patch%n
  end if
  self%quality = 1024
  self%level   = patch%level
  self%time    = patch%time
  do i=1,3
    self%mesh(i)%h  = patch%mesh(i)%h
  end do
  self%mass = 0d0
  self%velocity = 0d0
  self%unsigned(self%idx%d) = .true.
  call self%set (bits%internal)
  !-----------------------------------------------------------------------------
  ! Pre-load radius and mask array for use in accretion calculations
  !-----------------------------------------------------------------------------
  allocate (self%mask(self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  allocate (self%r   (self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  call io%bits_mem (storage_size(self%r   ), product(shape(self%r   )), 'sink%r')
  call io%bits_mem (storage_size(self%mask), product(shape(self%mask)), 'sink%mask')
  do iz=1,self%mesh(3)%gn
  do iy=1,self%mesh(2)%gn
  do ix=1,self%mesh(1)%gn
    r2 = (ix-self%mesh(1)%o)**2 &
       + (iy-self%mesh(2)%o)**2 &
       + (iz-self%mesh(3)%o)**2
    self%r(ix,iy,iz) = max(sqrt(r2),1e-3)
    self%mask(ix,iy,iz) = merge(.true., .false., r2<accretion_radius**2)
  end do
  end do
  end do
  call trace%end()
END SUBROUTINE init_task

!===============================================================================
!> Download values to all the sinkparticle%patch cells, but only to the 
!> density and momentum slots.  This relies on that the nbor list is 
!> comprehensive, including also the parent patch
!===============================================================================
SUBROUTINE dnload (self, only)
  class (sink_patch_t):: self
  integer, optional:: only
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%dnload')
  call download%download_link (self%link, all_cells=.true., only=self%idx%d)
  call download%download_link (self%link, all_cells=.true., only=self%idx%px)
  call download%download_link (self%link, all_cells=.true., only=self%idx%py)
  call download%download_link (self%link, all_cells=.true., only=self%idx%pz)
  call trace%end()
END SUBROUTINE dnload

!===============================================================================
!===============================================================================
SUBROUTINE update (self)
  class (sink_patch_t):: self
  call trace%begin ('sink_patch_t%update')
  call self%particle_solver%force_field (self)
  call self%accrete
  call self%move
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  call trace%end()
END SUBROUTINE update

!===============================================================================
!===============================================================================
SUBROUTINE accrete (self)
  class (sink_patch_t):: self
  !.............................................................................
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, px, py, pz
  real(kind=KindScalarVar), dimension(:,:,:), allocatable:: u2
  real:: v2_kepler
  integer:: ix, iy, iz
  real:: omega, dv
  real(8):: p, q, dm
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%accrete')
  allocate (u2 (self%mesh(1)%gn,self%mesh(2)%gn,self%mesh(3)%gn))
  !-----------------------------------------------------------------------------
  ! Short hand
  !-----------------------------------------------------------------------------
  d  => self%mem(:,:,:,self%idx%d ,self%it,1)
  px => self%mem(:,:,:,self%idx%px,self%it,1)
  py => self%mem(:,:,:,self%idx%py,self%it,1)
  pz => self%mem(:,:,:,self%idx%pz,self%it,1)
  u2 = (px/d)**2
  u2 = (py/d)**2 + u2
  u2 = (px/d)**2 + u2
  self%u_max = sqrt(self%fmaxval(u2))
  deallocate (u2)
  call self%courant_condition
  dm = 0d0
  dv = product(self%ds)
  !print *, 'XX0', minval(d), maxval(d), self%fmaxval(d)
  !print *, 'XX1', dv, self%mesh%n, self%mesh%gn
  if (self%mass==0d0) then
    do iz=self%mesh(3)%li,self%mesh(3)%ui
    do iy=self%mesh(2)%li,self%mesh(2)%ui
    do ix=self%mesh(1)%li,self%mesh(1)%ui
      q = min(rho_limit/d(ix,iy,iz),1d0)
      q = merge (q, 1d0, self%mask(ix,iy,iz))
      p = 1d0-q
      if (self%mask(ix,iy,iz)) then
        dm          = dm          + d (ix,iy,iz)*dv*p
        if (verbose>1) &
          print '(3i4,3f12.6,l4)', ix, iy, iz, d(ix,iy,iz), self%r(ix,iy,iz), p, self%mask(ix,iy,iz)
        self%mom(1) = self%mom(1) + px(ix,iy,iz)*dv*p
        self%mom(2) = self%mom(2) + py(ix,iy,iz)*dv*p
        self%mom(3) = self%mom(3) + pz(ix,iy,iz)*dv*p
        d (ix,iy,iz) = d (ix,iy,iz)*q
        px(ix,iy,iz) = px(ix,iy,iz)*q
        py(ix,iy,iz) = py(ix,iy,iz)*q
        pz(ix,iy,iz) = pz(ix,iy,iz)*q
      end if
    end do
    end do
    end do
    !print *, '1 id,dm:', self%id, dm
  !-----------------------------------------------------------------------------
  ! Remove a fraction p of per-volume quantities and keep q=(1-p)
  !-----------------------------------------------------------------------------
  else
    !print *, 'XX', scaling%grav, self%mass, minval(self%r), maxval(self%r)
    do iz=1,self%mesh(3)%gn
    do iy=1,self%mesh(2)%gn
    do ix=1,self%mesh(1)%gn
      omega = sqrt(scaling%grav*self%mass/self%r(ix,iy,iz))/self%r(ix,iy,iz)
      p = merge (accretion_rate*self%dtime*omega, 0d0, self%mask(ix,iy,iz))
      q = 1d0-p
      dm          = dm          + d (ix,iy,iz)*dv*p
      self%mom(1) = self%mom(1) + px(ix,iy,iz)*dv*p
      self%mom(2) = self%mom(2) + py(ix,iy,iz)*dv*p
      self%mom(3) = self%mom(3) + pz(ix,iy,iz)*dv*p
      d (ix,iy,iz) = d (ix,iy,iz)*q
      px(ix,iy,iz) = px(ix,iy,iz)*q
      py(ix,iy,iz) = py(ix,iy,iz)*q
      pz(ix,iy,iz) = pz(ix,iy,iz)*q
    end do
    end do
    end do
    !print *, '2 id,dm:', self%id, dm
  end if
  self%mass = self%mass + dm
  self%velocity = self%mom/self%mass
  if (verbose > 0) &
    print '(4(a,1p,e12.4,2x),2(a,3e12.4,2x))', &
      'sink: time=', self%time, 'mass=', self%mass, &
      'dmin=',self%fminval(self%idx%d), 'dmax=',self%fmaxval(self%idx%d), &
      'vel=',self%velocity, 'pos=',self%position
  call trace%end()
END SUBROUTINE accrete

!===============================================================================
!> Move the sink particle one time step self%dt, using kidk-drift-kick
!===============================================================================
SUBROUTINE move (self)
  class (sink_patch_t):: self
  real(8):: acceleration(3)
  real, dimension(:,:,:), pointer:: phi
  integer:: ix, iy, iz, ii(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%move')
  !-----------------------------------------------------------------------------
  ! Kick. The sink is assumed to be in the exact grid position floor(mesh%o)
  !-----------------------------------------------------------------------------
  !phi = self%mem(:,:,:,self%idx%phi,self%it,1)
  ii = floor(self%mesh%o)
  ix = ii(1); iy=ii(2); iz=ii(3)
  self%acceleration = 0.0
  !self%acceleration = [(phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  ))/(2.*self%ds(1)), &
  !                     (phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  ))/(2.*self%ds(2)), &
  !                     (phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1))/(2.*self%ds(3))]
  self%velocity = self%velocity + 0.5d0*self%acceleration*self%dtime
  !-----------------------------------------------------------------------------
  ! Drift + update phi to the new space-time position
  !-----------------------------------------------------------------------------
  self%time     = self%time     + self%dtime
  self%position = self%position + self%dtime*self%velocity    
  !call download%download_link (self%link, all_cells=.true., only=self%idx%phi)
  !-----------------------------------------------------------------------------
  ! Kick
  !-----------------------------------------------------------------------------
  self%acceleration = 0.0
  !self%acceleration = [(phi(ix+1,iy  ,iz  )-phi(ix-1,iy  ,iz  ))/(2.*self%ds(1)), &
  !                     (phi(ix  ,iy+1,iz  )-phi(ix  ,iy-1,iz  ))/(2.*self%ds(2)), &
  !                     (phi(ix  ,iy  ,iz+1)-phi(ix  ,iy  ,iz-1))/(2.*self%ds(3))]
  self%velocity = self%velocity + 0.5d0*self%acceleration*self%dtime
  call trace%end()
END SUBROUTINE move

!===============================================================================
!===============================================================================
INTEGER FUNCTION check_refine (self, patch)
  class(sink_patch_t):: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  check_refine = -1
  call self%refine_check (patch%link, check_refine)
END FUNCTION

!===============================================================================
!> Test for the need to create a new sink particle, and handle the task list
!> aspect of that; creating an nbor list, and adding the link to the task list.
!===============================================================================
SUBROUTINE refine_check (self, link, refine)
  class(sink_patch_t):: self
  class(link_t), pointer:: link
  integer:: refine
  !.............................................................................
  logical:: check_sinkparticle
  type(list_t), pointer:: task_list
  class(task_t), pointer:: task, patch
  class(gpatch_t), pointer:: gpatch
  class(link_t), pointer:: new_link, nbor
  class(sink_patch_t), pointer:: new_sink
  class(dll_node_t), pointer:: particle
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  integer:: m(3)
  logical:: do_trace
  !-----------------------------------------------------------------------------
  check_sinkparticle = .false.
  if (.not.sink_patch%on) return
  task_list => self%task_list
  patch => link%task
  select type (patch)
  class is (sink_patch_t)
    !---------------------------------------------------------------------------
    ! Refresh nbor list
    !---------------------------------------------------------------------------
    call task_list%refresh_nbors(link)
  class is (gpatch_t)
    gpatch => patch
    check_sinkparticle = sink_patch%check (gpatch)
    if (check_sinkparticle) then
      do_trace = io%do_trace
      !io%do_trace = .true.
      !-------------------------------------------------------------------------
      ! Allocate a new sink_patch_t instance, as the task of a new link
      !-------------------------------------------------------------------------
      allocate (new_link)
      allocate (sink_patch_t:: task)
      allocate (particle_t:: particle)
      new_link%parent => link
      new_link%task => task
      !-------------------------------------------------------------------------
      ! Create the new sinkparticle, using data from the MHD patch
      !-------------------------------------------------------------------------
      select type (task)
      class is (sink_patch_t)
        task%link => new_link
        task%particle_solver => self%particle_solver
        task%particle => particle
        call self%particle_solver%append (particle)
        call task%init_task (gpatch)
        !-------------------------------------------------------------------------
        ! Initialize its neighbor list, and add it to the task list
        !-------------------------------------------------------------------------
        call task_list%init_nbors (new_link)
        call task_list%append_link (new_link)
        nbor => new_link%nbor
        do while (associated(nbor))
          print *, 'nbor:', nbor%task%id, real((nbor%task%position-task%position)/task%ds)
          nbor => nbor%next
        end do
        !------------------------------------------------------------------------
        ! Download variables and accrete
        !------------------------------------------------------------------------
        call task%dnload
        d => task%mem(:,:,:,task%idx%d,task%it,1)
        d = max(d,task%fminval(d))
        if (task%id==io%id_debug) then
          m = nint(task%offset)
          print *, 'sink_patch_t%create:  task%pos, d',  task%myposition(m), &
            task%mem(m(1),m(2),m(3),task%idx%d,task%it,1)
          print *,'MM', minval(d), task%fminval(task%idx%d), &
                        maxval(d), task%fmaxval(task%idx%d)
        end if
        call task%accrete
      end select
      !-------------------------------------------------------------------------
      ! Refresh the nbor lists of the parent and its nbors; the parent should
      ! be done after the nbors, to avoid putting the child link on its nbor
      ! list until the other ones are updated.  All nbors that happen to have
      ! overlap with the sink particle patch will have it as an nbor
      !-------------------------------------------------------------------------
      nbor => link%nbor
      do while (associated(nbor))
        call task_list%init_nbors (nbor%link)
        nbor => nbor%next
      end do
      call task_list%init_nbors (link)
      io%do_trace = do_trace
    end if
  end select
END SUBROUTINE refine_check

END MODULE sink_patch_mod
