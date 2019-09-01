!===============================================================================
!> A sink particle is represented as primarily a particle task, whose position
!> is updated with a particle solver. The patch is usually much smaller than an
!> MHD patch, and its values are downloaded from contributing MHD patches.
!>
!> A particular quality level indicates to the download_mod that all cells from
!> the sink particle patch should be downloaded to overlapping patches, thus
!> communicating back the removal of accreted mass and momentum.  As a target,
!> the sink_patch also should have values filled in from all overlapping patches,
!> arranged so the values from the highest levels are preferred.
!>
!> Often, there is only one patch that overlaps, and that patch must them be at
!> levelmax, which means it is easy to arrange for a fast and mass conserving
!> fill-in.  When a sink particle is close to a patch boundary there may be
!> contributions from nbor patches, which may have a level lower by one (unless
!> one enforces the level of all overlapping patches to be levelmax).
!>
!> Object hierarchy:
!>
!> sink_patch_t               ! Contains both a particle and a patch object
!>   gpatch_t                 ! A grid-based task holding accretion related data
!>   particle_t
!>   particle_solver_t        ! Particle solvers use particles_t%force
!>     particle_t
!>     particles_t        ! A list of particle_t objects
!>       particle_t           ! Each particle_t has a patch_forces_t list
!>         dll_node_t
!>         patch_forces_t     ! A list of patch_force_t objects
!>           dll_t
!>           patch_force_t    ! An extension of a dll_note_t list node
!>             dll_node_t
!>
!> Call sequences (relevant files: _t -> _mod.f90):
!>
!> NORMAL CREATION OF SINK PATCH:
!> ==============================
!> extras_t%check_refine                             ! refine caused by extras
!>   sink_patch_t%check_refine                       ! refine/derefine by sink
!>     sink_patch_t%refine_check                     ! new sink needed?
!>       sink_patch_t%check                          ! criteria for sink creation
!>       sink_patch_t%init_sink                      ! initialize new sink
!>         particle_t%init                           ! initialize particle
!>         particle_solver_t%append                  ! append particle to solver
!>         sink_patch_t%init_patch                   ! initialize accretion patch
!>           patch_t%setup                           ! setup patch mesh & vars
!>             task_t%init                           ! assign a task ID
!>           idx_t%init                              ! set mem slot indices
!>         sink_patch_t%init_hash                    ! initialize hash table
!>         list_t%init_nbors                         ! make nbor lists
!>         sink_patch_t%dnload                       ! fill in patch values
!>         sink_patch_t%accrete                      ! accrete to it
!>         list_t%append                             ! add to task list
!>
!> AD HOC CREATION OF SINK PATCHES:
!> ================================
!> task_list_t%init_task_list_pointers               ! pass task list to tasks
!>   extras_t%init_task_list                         ! pass task_list to extras
!>     gpatch_t%init_task_list                       ! save task_list pointer
!>     sink_patch_t%ad_hoc_add                       ! add ad hoc sink patches
!>       sink_patch_t%init_sink                      ! initialize sink
!>       ...                                         ! same as above
!>
!> UPATE OF SINK PATCH:
!> ====================
!> task_list_t%update                                 ! dispatcher selects task
!>   gpatch_t%dnload                                  ! standard call
!>     download_t%download_link                       ! special case quality
!>   sink_patch_t%update                              ! update particles & patches
!>     particle_solver_t%force                        ! fall-through
!>       particles_t%force                        ! force from particles
!>         patch_forces_t%force                       ! force from all patches
!>           patch_force_t%force                      ! force from one patches
!>     sink_patch_t%nove                              ! move particles
!>       particle_solver_t%update                     ! actual particle solver
!>         particlesolver_t%courant_time              ! determine time step
!>     sink_patch_t%accrete                           ! accrete to particles
!>       sink_patch_t%courant_condition               ! copy particle time step
!>
!> UPDATE OF MHD PATCH:
!> ====================
!> experiment_t%update                                ! update MHD patch
!>   solver_t%update                                  ! solver interface
!>     extras_t%pre_update                            ! extras coupler
!>       sink_patch%force_field                       ! interface
!>         particles_t%force_field                ! force from all sinks
!>
!> Related files:
!>   lists/dll_mod.f90                                ! doubly-linked list
!>   lists/hash_table_mod.f90                         ! hash table
!>   solvers/particles/particle_mod.f90               ! particle data type
!>   solvers/particles/patch_force_mod.f90            ! time slices of forces
!>   solvers/particles/particles_mod.f98          ! handle
!>   solvers/particles/KDK/particle_solver_mod.f90    ! KDK particle solver
!===============================================================================
MODULE sink_patch_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
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
  USE patch_force_mod
  USE particle_mod
  USE particle_solver_mod
  USE hash_table_mod
  USE omp_timer_mod
  USE shared_mod
#ifdef HDF5
  USE h5_mod
  USE hdf5, only : HID_T
#endif
  implicit none
  private
  type :: io_t
    ! name and id handles for the different datasets
    character(len=9) :: name
#ifdef HDF5
    integer(HID_T) :: gid, time_id, mass_id, pos_id, vel_id, mom_id
#endif
  end type
  type, public, extends (gpatch_t):: sink_patch_t
    real(8):: mass=0.0_8, mom(3)
    real(8):: acceleration(3)=0.0_8
    real(8), dimension(:,:)  , allocatable:: p, v, a
    real(4), dimension(:,:,:), allocatable:: r
    class(particle_solver_t), pointer:: particle_solver => null()
    class(particle_t), pointer:: particle => null()
    class(gpatch_t), pointer:: patch => null()
    logical:: on
    type(io_t):: io
  contains
    procedure:: init
    procedure:: init_patch
    procedure:: init_patch_forces
    procedure:: ad_hoc_add
    procedure:: dealloc
    procedure:: check
    procedure:: courant_condition
    procedure:: force_field
    procedure:: update
    procedure:: accrete
    procedure:: move
    procedure:: check_refine
    procedure:: output
    procedure, private:: refine_check
  end type
  !-----------------------------------------------------------------------------
#ifdef HDF5
  integer(HID_T) :: sink_file_id=-1   ! id for sink file. -1: file not open
#endif
  logical:: on=.false.                ! turn feature on/off
  logical:: do_translate=.true.       ! make a transformation to the velocity of the sink
  logical:: force_refine=.true.       ! insist on refining down to levelmax at star
  logical:: only_high_levels=.true.   ! only accrete at high levels
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
  real(8):: out_time=0.1              ! cadence of output into sink.h5
  character(len=64):: file ='sink.h5' ! name of sink data file
  !-----------------------------------------------------------------------------
  integer, save:: n_sinks=0           ! counter
  integer, save:: n_ad_hoc=0          ! number of ad_hoc sinks
  integer, parameter:: hash_size=10000
  type(sink_patch_t), public:: sink_patch
CONTAINS

!===============================================================================
!> Initialize a static instance of the sink_patch_t data type, with parameters
!> for sink creation read from a sink_patch_params namelist.
!===============================================================================
SUBROUTINE init (self)
  class(sink_patch_t)   :: self
  !.............................................................................
  logical, save:: first_time=.true.
  integer:: iostat
  namelist /sink_patch_params/ &
    on, do_translate, force_refine, only_high_levels, center_star, verbose, &
    accretion_rate, accretion_efficiency, accretion_fraction, accretion_radius, &
    exclusion_radius, jeans_resolution, rho_limit, rho_limit_factor, rho_cut, &
    max_sinks, softening_length, out_time, file, n_ad_hoc
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io_unit%input)
    read (io_unit%input, sink_patch_params, iostat=iostat)
    if (io_unit%master) write (io_unit%output, sink_patch_params)
    n_sinks = n_ad_hoc
    sink_patch%on = on
    sink_patch%particle_solver => particle_solver
    call hash_table%init (hash_size)
    call particle_solver%init
  endif
  !$omp end critical (input_cr)
  self%on = on
  self%particle_solver => particle_solver
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize an instance of the sink_patch_t data type
!===============================================================================
SUBROUTINE init_sink (self, patch, sink_position)
  type(sink_patch_t), pointer:: self
  class(gpatch_t)    , pointer:: patch
  real(kind=8), dimension(3), optional :: sink_position
  !.............................................................................
  class(link_t)      , pointer:: link, task_link
  class(task_t)      , pointer:: task
  class(particle_t)  , pointer:: particle
  type(patch_force_t), pointer:: patch_force
  class(dll_node_t)  , pointer:: dll_node
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  integer:: id, m(3)
  !-----------------------------------------------------------------------------
  ! Allocate a new sink_patch_t instance, as the task of a new link
  !-----------------------------------------------------------------------------
  call self%init
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%init_sink')
  allocate (task_link)
  task_link%parent => patch%link
  task_link%task => self
  self%link => task_link
  !-----------------------------------------------------------------------------
  !$omp atomic update
  n_sinks = n_sinks + 1
  self%task_list => patch%task_list
  self%particle_solver => particle_solver
  !-----------------------------------------------------------------------------
  ! Allocate and initialize the sink particle, and append it to the solver
  !-----------------------------------------------------------------------------
  allocate (particle)
  call particle%init
  particle%mass = tiny(1.0)
  self%particle => particle
  dll_node => particle
  call self%particle_solver%append (dll_node)
  !-----------------------------------------------------------------------------
  ! Initialize the accretion patch and its neighbor list, and download values
  !-----------------------------------------------------------------------------
  call self%init_patch (patch, sink_position)
  particle%ds = minval(patch%ds)
  call self%init_patch_forces
  !-----------------------------------------------------------------------------
  ! Allocate a task link, set points, make an nbor list, append to task
  ! list, and accrete for the first time
  !-----------------------------------------------------------------------------
  allocate (task_link)
  task_link%task => self
  self%link => task_link
  call self%task_list%init_nbors (task_link)
  link => task_link%nbor
if (verbose > 0) &
print *, self%id, self%position
  do while (associated(link))
if (verbose > 0) &
print *, link%task%id, link%task%position
    link => link%next
  end do
  call self%dnload
  call self%task_list%append (task_link)
  !-----------------------------------------------------------------------------
  ! Check guard zone values (should it all have any???) 
  !-----------------------------------------------------------------------------
  if (verbose > 1) then
    d => self%mem(:,:,:,self%idx%d,self%it,1)
    m = nint(self%offset)
if (verbose > 0) &
print *, 'sink_patch_t%create:  self%pos, d',  self%myposition(m), &
self%mem(m(1),m(2),m(3),self%idx%d,self%it,1)
if (verbose > 0) &
print *,'MM', minval(d), self%fminval(self%idx%d), &
              maxval(d), self%fmaxval(self%idx%d)
  end if
  call trace%end()
END SUBROUTINE init_sink

!===============================================================================
!> Create the accompanying cell-based patch, which records the density and
!> velocity from overlapping MHD patches, for use in the recipe for accretion.
!===============================================================================
SUBROUTINE init_patch (self, patch, sink_position)
  class(sink_patch_t):: self
  class(gpatch_t), target:: patch
  real(kind=8), dimension(3), optional :: sink_position
  !.............................................................................
  class(link_t), pointer:: nbor, prev
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  real(8), dimension(:), pointer:: r
  real(kind=8) :: pos(3)
  integer:: i, l(3), u(3), m(3), ix, iy, iz
  real:: r2, dmin
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%init_patch')
  self%patch => patch
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
  if (present(sink_position)) then
    pos = sink_position
  else
    do i=1,3
      r => patch%mesh(i)%r
      pos(i) = patch%position(i) + r(m(i))
    end do
  endif
  self%position = pos
  self%ds = patch%ds
  self%n  = 2*ceiling(accretion_radius)
  self%ng = 2                               ! FIXME: set =0 ?
  self%nv = patch%nv
  self%nt = patch%nt
  self%nw = 1
  !-----------------------------------------------------------------------------
  ! We should use no_mans_land=.true., and an even number for self%n, so the
  ! sink particle position is not at a cell center.
  !-----------------------------------------------------------------------------
  self%no_mans_land = .true.
  !-----------------------------------------------------------------------------
  ! Set a special quality, to force downloading of interior data
  !-----------------------------------------------------------------------------
  self%size = self%n*self%ds
  self%quality = 1024
  call self%patch_t%setup (self%position, self%size, self%n, self%ng, self%nt, &
                           self%nv, self%nw)
  if (self%id==io%id_debug .or. verbose > 0) then
    print *, 'sink_patch_t%init_task: id, patch%id, m, dmax', &
      self%id, patch%id, m, maxval(d)
    print *, 'sink_patch_t%init_task: (pos-patch%pos)/size', &
      (self%position-patch%position)/patch%size
    print *, 'sink_patch_t%init_task: patch%pos, d', &
      patch%myposition(m), patch%mem(m(1),m(2),m(3),patch%idx%d,patch%it,1)
    print *, 'sink_patch_t%init_task:  self%size, patch%size', &
      self%size, patch%size
  end if
  self%mem = 1.0
  self%quality = 1024
  self%level   = patch%level
  self%time    = patch%time
  do i=1,3
    self%mesh(i)%h  = patch%mesh(i)%h
  end do
  self%mass = 0d0
  self%velocity = 0d0
  call self%idx%init (4, mhd=.false.)
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
  !-----------------------------------------------------------------------------
  ! Attach a task list link, initialize nbors, and append to the task list
  !-----------------------------------------------------------------------------
  self%task_list => patch%task_list
  call self%task_list%init_nbors (self%link)
  call self%task_list%prepend_link (self%link)
  if (verbose > 1) then
    nbor => self%link%nbor
    do while (associated(nbor))
      write (stdout,*) &
        'nbor:', nbor%task%id, real((nbor%task%position-self%position)/self%ds)
      nbor => nbor%next
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Refresh the MHD patch nbor list, which adds the sink patch to it
  !-----------------------------------------------------------------------------
  call refresh_nbor_list (patch%link, self%task_list)
  call trace%end()
END SUBROUTINE init_patch

!===============================================================================
!> Add a list of patch_force_t instances to the particle%patch_forces list, 
!> indexing each patch force with the particle and patch ids in a hash table, 
!> for quick look-up.  The upshot is that each sink particle has as an attribute
!> a list of patch_force_t data types, one for each patch.  Any patch force,
!> for any combination of particle and patch, may also be located via the
!> (global) hash table, which is much faster than searching the list of patch
!> forces for a particle.
!===============================================================================
SUBROUTINE init_patch_forces (self)
  class(sink_patch_t):: self
  !.............................................................................
  class(list_t), pointer:: task_list
  class(link_t), pointer:: link
  class(task_t), pointer:: task
  class(particle_t), pointer:: particle
  class(dll_node_t), pointer:: dll_node
  class(*), pointer:: ptr
  type(patch_force_t), pointer:: patch_force
  !.............................................................................
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%init_patch_forces')
  !-----------------------------------------------------------------------------
  ! For each patch in the task list, add a patch_force_t instance to the
  ! particle%patch_forces list
  !-----------------------------------------------------------------------------
  particle => self%particle
  if (verbose > 1) &
    write (stdout,*) 'sink_patch_t: particle id =', particle%id 
  link => self%task_list%head
  do while (associated(link))
    task => link%task
    select type (task)
    !---------------------------------------------------------------------------
    ! Catch and ignore sink_patch_t tasks
    !---------------------------------------------------------------------------
    class is (sink_patch_t)
    !---------------------------------------------------------------------------
    ! Act on normal MHD tasks
    !---------------------------------------------------------------------------
    class is (gpatch_t)
      allocate (patch_force)
      ptr => patch_force
      if (verbose > 1) &
        write (stdout,*) 'hash_table%set key =', [particle%id, task%id]
      call hash_table%set ([particle%id, task%id], ptr)
      call patch_force%init
      patch_force%id = task%id
      dll_node => patch_force
      call particle%patch_forces%append (dll_node)
    end select
    link => link%next
  end do
  call trace%end()
END SUBROUTINE init_patch_forces

!===============================================================================
!> Initialize test instances of the sink_patch_t data type, with parameters for
!> sink creation read from a sink_params namelists.
!===============================================================================
SUBROUTINE ad_hoc_add (self, link)
  class(sink_patch_t)   :: self
  class(link_t), pointer:: link
  !.............................................................................
  class(task_t), pointer:: patch
  class(gpatch_t), pointer:: gpatch
  type(sink_patch_t), pointer:: task
  class(list_t), pointer:: task_list
  class(particle_t), pointer:: particle
  integer:: iostat, i_sink, m(3)
  real(kind=8):: pos(3)
  type sinks_t
    real(kind=8) :: mass, p(3), v(3)
  end type
  type(sinks_t), allocatable, save:: sinks(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  real(kind=8), save:: mass=1d0, x=0.5d0, y=0.5d0,z, vx=0d0, vy=0d0, vz=0d0
  namelist /sink_params/ mass, x, y, z, vx, vy, vz
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%ad_hoc_add')
  patch => link%task
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    !---------------------------------------------------------------------------
    ! Read in sink particles directly from namelist -- only done by master rank
    !---------------------------------------------------------------------------
    if (io_unit%master) then
      rewind (io_unit%input)
      allocate (sinks(n_ad_hoc))
      do i_sink=1,n_ad_hoc
        read (io_unit%input, sink_params, iostat=iostat)
        if (io_unit%master) write (io_unit%output, sink_params)
        sinks(i_sink)%mass = mass
        sinks(i_sink)%p = [x,y,z]
        sinks(i_sink)%v = [vx, vy, vz]
      enddo
    endif
  endif
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  ! Go through list of hand made sinks, and create them, if in patch
  !-----------------------------------------------------------------------------
  select type(patch)
  class is (gpatch_t)
    gpatch => patch
    task_list => patch%task_list
    do i_sink=1,n_ad_hoc
      pos = sinks(i_sink)%p
if (verbose > 0) &
print *, patch%id, i_sink, pos, patch%contains(pos)
      if (patch%contains(pos)) then
        !-----------------------------------------------------------------------
        ! TEST: Place the sink exactly at the patch position
        !-----------------------------------------------------------------------
        pos = patch%position
        sinks(i_sink)%p = pos
        if (verbose > 0) then
          print *, 'adding ad hoc sink to', patch%id, i_sink
          print *, 'relative position:', (pos-patch%position)/patch%size
        end if
        !-----------------------------------------------------------------------
        ! Initialize the sink_patch task (includes defining self%task_list)
        !-----------------------------------------------------------------------
        allocate (task)
        call init_sink (task, gpatch, sinks(i_sink)%p)
        task_list => task%task_list
        self%task_list => task%task_list
        task%mass = sinks(i_sink)%mass
        task%mom = sinks(i_sink)%mass * sinks(i_sink)%v
        particle => task%particle
        particle%ds = minval(patch%ds)
        particle%mass = sinks(i_sink)%mass
        particle%m(:) = particle%mass
        particle%r(1,:) = sinks(i_sink)%p(1)
        particle%r(2,:) = sinks(i_sink)%p(2)
        particle%r(3,:) = sinks(i_sink)%p(3)
        particle%v(1,:) = sinks(i_sink)%v(1)
        particle%v(2,:) = sinks(i_sink)%v(2)
        particle%v(3,:) = sinks(i_sink)%v(3)
if (verbose > 0) &
print *, 'sink_patch_t%init_sink: id,mass,velocity =', &
particle%id,particle%mass, particle%v(:,1)
        task%dtime = self%particle_solver%courant_time (particle)
        task%dt = self%dtime
        call task%accrete
      endif
    enddo
  end select
  call trace%end()
END SUBROUTINE ad_hoc_add

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
!> Output of sink data to hdf5
!===============================================================================
SUBROUTINE output (self)
  class(sink_patch_t):: self
  !.............................................................................
  logical:: new
  integer:: verbose=1
  real(kind=8),save::last_flush
  real(kind=8):: flush_frequency = 10. ! frequency with which data is flushed to disk
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%output')
#ifdef HDF5
  !-----------------------------------------------------------------------------
  ! Open file if needed
  !-----------------------------------------------------------------------------
  if (sink_file_id==-1) then
    !$omp critical
    if (sink_file_id==-1) then
      last_flush = wallclock()
      sink_file_id = h5%open (file, new, verbose=verbose)
      if (verbose > 0) &
        print *, merge('creating ','updating ',new)//trim(file)
    endif
    !$omp end critical
  endif
  !-----------------------------------------------------------------------------
  ! Open group for sink if needed
  !-----------------------------------------------------------------------------
  if (self%io%gid == -1) then
    write(self%io%name,'(a,i5.5)') 'sink', self%id
    self%io%gid = h5%group_open (sink_file_id, self%io%name)
    self%io%time_id = h5%set_open_ext_1d (self%io%gid, 'time', (/self%time/))
    self%io%mass_id = h5%set_open_ext_1d (self%io%gid, 'mass', (/self%mass/))
    self%io%pos_id  = h5%set_open_ext_1d (self%io%gid, 'position', self%position)
    self%io%vel_id  = h5%set_open_ext_1d (self%io%gid, 'velocity', self%mom / self%mass)
    self%io%mom_id  = h5%set_open_ext_1d (self%io%gid, 'momentum', self%mom)
  endif
  !-----------------------------------------------------------------------------
  ! Write sink properties
  !-----------------------------------------------------------------------------
  call h5%set_write_ext_1d (self%io%time_id, (/self%time/))
  call h5%set_write_ext_1d (self%io%mass_id, (/self%mass/))
  call h5%set_write_ext_1d (self%io%pos_id,  self%position)
  call h5%set_write_ext_1d (self%io%vel_id,  self%mom / self%mass)
  call h5%set_write_ext_1d (self%io%mom_id,  self%mom)
  !-----------------------------------------------------------------------------
  ! Flush sink file to disk at least every flush_frequency seconds
  !-----------------------------------------------------------------------------
  if (wallclock() > last_flush + flush_frequency) then
    !$omp critical
    if (wallclock() > last_flush + flush_frequency) then
      last_flush = last_flush + flush_frequency
      call h5%flush(sink_file_id)
    endif
    !$omp end critical
  endif
#endif
  call trace%end()
END SUBROUTINE output

!===============================================================================
!> Test if a sink particle should be created in a given patch
!===============================================================================
FUNCTION check (self, patch)
  class(sink_patch_t):: self
  class(gpatch_t), target:: patch
  logical:: check
  !-----------------------------------------------------------------------------
  check = .false.
  if (.not.self%on) return
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
!> Courant condition for accretion = use the sink particle timestep
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(sink_patch_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%courant_condition')
  self%dtime = self%particle%dtime
  call trace%end()
END SUBROUTINE courant_condition

!===============================================================================
!> Compute the force field from sink particles on the cells of a patch
!===============================================================================
SUBROUTINE force_field (self, patch)
  class(sink_patch_t):: self
  class(patch_t):: patch
  !.............................................................................
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%force_field')
  call self%particle_solver%force_field (patch)
  call trace%end()
END SUBROUTINE force_field

!===============================================================================
!> Update the sink particle task, by accreting mass, and by moving the sink
!===============================================================================
SUBROUTINE update (self)
  class (sink_patch_t):: self
  integer, save:: itimer=0
  !.............................................................................
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%update', itimer=itimer)
  call self%move
  call self%accrete
  call trace%end (itimer)
END SUBROUTINE update

!===============================================================================
!> Accrete mass and momentum from the sinkparticle patch
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
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%accrete', itimer=itimer)
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
  if (self%particle%mass==0d0) then
    do iz=self%mesh(3)%li,self%mesh(3)%ui
    do iy=self%mesh(2)%li,self%mesh(2)%ui
    do ix=self%mesh(1)%li,self%mesh(1)%ui
      q = min(rho_limit/d(ix,iy,iz),1d0)
      q = merge (q, 1d0, self%mask(ix,iy,iz))
      p = 1d0-q
      if (self%mask(ix,iy,iz)) then
        dm          = dm          + d (ix,iy,iz)*dv*p
if (verbose > 1) &
print '(a,3i4,1p,4e11.3,l4)', ' accrete: ii,d,r,dm,mask =', &
ix, iy, iz, d(ix,iy,iz), self%r(ix,iy,iz), p, dm
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
    if (verbose > 0) &
      print *, 'accrete(1): id,dm =', self%id, dm
  !-----------------------------------------------------------------------------
  ! Remove a fraction p of per-volume quantities and keep q=(1-p), where p is
  ! the fraction of an orbit that would be covered by in the time step, assuming
  ! Keplerian rotation. The mask variable makes sure accretion is only done
  ! inside the accretion radius.
  !-----------------------------------------------------------------------------
  else
    do iz=1,self%mesh(3)%gn
    do iy=1,self%mesh(2)%gn
    do ix=1,self%mesh(1)%gn
      if (self%mask(ix,iy,iz)) then
        omega = sqrt(scaling%grav*self%mass/self%r(ix,iy,iz))/self%r(ix,iy,iz)
        p = merge (accretion_rate*self%dtime*omega, 0d0, self%mask(ix,iy,iz))
        q = 1d0-p
if (verbose > 1) &
print '(a,3i4,1p,4e11.3,l4)', 'accrete: ii,d,r,dm,mask =', &
ix, iy, iz, d(ix,iy,iz), self%r(ix,iy,iz), p, dm
        dm          = dm          + d (ix,iy,iz)*dv*p
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
    if (verbose > 0) &
      print *, 'accrete(2): id,dm =', self%id, dm
  end if
  self%mass = self%mass + dm
  if (self%mass > 0.0) then
    self%velocity = self%mom/self%mass
  else
    self%velocity = 0.0
  end if
  if (verbose > 0) &
    print '(i4,2x,4(a,1p,e12.4,2x),2(a,3e11.3,2x))', self%id, &
      'sink: time=', self%time, 'mass=', self%mass, &
      'dmin=',self%fminval(self%idx%d), 'dmax=',self%fmaxval(self%idx%d), &
      'vel=',self%velocity, 'pos=',self%position
  call trace%end (itimer)
END SUBROUTINE accrete

!===============================================================================
!> Move the sink particle one time step, using the particle solver, which in
!> turn uses the particles_t%force() procedure to calculate the force from
!> other particles and from the patches.   Call hierarchy:
!>
!> move(self)
!>   particle_solver%update (self%particle)
!>     particles_t%force (self%particle)
!>       for other particles:
!>         force = force - other%mass*r*r3
!>       force = force + particle%patch_forces%force (particle%time)
!===============================================================================
SUBROUTINE move (self)
  class (sink_patch_t):: self
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('sink_patch_t%move', itimer=itimer)
  call self%particle_solver%update (self%particle)
  call trace%end (itimer)
END SUBROUTINE move

!===============================================================================
!> Interface to the refine_t%refine_needed procedure.
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
  class(task_t), pointer:: patch
  class(gpatch_t), pointer:: gpatch
  class(sink_patch_t), pointer:: new_sink
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  integer:: id, m(3)
  !-----------------------------------------------------------------------------
  if (.not.self%on) return
  call trace%begin ('sink_patch_tA%refine_check')
  patch => link%task
  select type (patch)
  !-----------------------------------------------------------------------------
  ! For sink_patch_t tasks: refresh nbor list (FIXME: should this be here???)
  !-----------------------------------------------------------------------------
  class is (sink_patch_t)
    call self%task_list%refresh_nbors(link)
  !-----------------------------------------------------------------------------
  ! For gpatch_t task: Check if a sink particle task should be created
  !-----------------------------------------------------------------------------
  class is (gpatch_t)
    gpatch => patch
    check_sinkparticle = sink_patch%check (patch)
    if (check_sinkparticle) then
      allocate (new_sink)
      call init_sink (new_sink, gpatch)
    end if
  end select
  call trace%end()
END SUBROUTINE refine_check

!-------------------------------------------------------------------------------
! Refresh the nbor lists of the parent and its nbors; the parent should
! be done after the nbors, to avoid putting the child link on its nbor
! list until the other ones are updated.  All nbors that happen to have
! overlap with the sink particle patch will have it as an nbor
!-------------------------------------------------------------------------------
SUBROUTINE refresh_nbor_list (link, task_list)
  implicit none
  class(link_t), pointer :: link
  class(list_t), pointer :: task_list
  class(link_t), pointer :: nbor
  !-----------------------------------------------------------------------------
  ! FIXME: This is a rather expensive way of updating the nbor lists
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%refresh_nbor_list')
  nbor => link%nbor
  do while (associated(nbor))
    call task_list%init_nbors (nbor%link)
    nbor => nbor%next
  end do
  call task_list%init_nbors (link)
  call trace%end()
END SUBROUTINE refresh_nbor_list

END MODULE sink_patch_mod
