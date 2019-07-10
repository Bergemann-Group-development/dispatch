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
  USE particle_mod
  USE particle_solver_mod
  USE omp_timer_mod
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
    class(dll_node_t), pointer:: particle => null()
    class(gpatch_t), pointer:: patch => null()
    logical:: on
    type(io_t):: io
  contains
    procedure:: init
    procedure:: init_sink
    procedure:: ad_hoc_add
    procedure:: dealloc
    procedure:: init_task_list
    procedure:: check
    procedure:: init_patch
    procedure:: courant_condition
    procedure:: dnload
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
  endif
  !$omp end critical (input_cr)
  self%on = on
  self%particle_solver => particle_solver
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Initialize an instance of the sink_patch_t data type
!===============================================================================
SUBROUTINE init_sink (self, patch)
  class(sink_patch_t), target:: self
  class(gpatch_t), pointer   :: patch
  !.............................................................................
  class(link_t), pointer:: new_link, nbor
  class(dll_node_t), pointer:: particle
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  integer:: id, m(3)
  !-----------------------------------------------------------------------------
  ! Allocate a new sink_patch_t instance, as the task of a new link
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_patch_t%init_sink')
  allocate (new_link)
  allocate (particle_t:: particle)
  new_link%parent => patch%link
  new_link%task => self
  !---------------------------------------------------------------------------
  !$omp atomic capture
  n_sinks = n_sinks + 1
  id = n_sinks
  !$omp end atomic
  self%id = id
  self%link => new_link
  self%particle_solver => particle_solver
  self%particle => particle
  call self%particle_solver%append (particle)
  call self%init_patch (patch)
  !---------------------------------------------------------------------------
  ! Initialize its neighbor list, and add it to the task list
  !---------------------------------------------------------------------------
  self%patch => patch
  self%task_list => patch%task_list
  call self%task_list%init_nbors (new_link)
  call self%task_list%append_link (new_link)
  nbor => new_link%nbor
  do while (associated(nbor))
    print *, 'nbor:', nbor%task%id, real((nbor%task%position-self%position)/self%ds)
    nbor => nbor%next
  end do
  !---------------------------------------------------------------------------
  ! Download variables and accrete
  !---------------------------------------------------------------------------
  call self%dnload
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  d = max(d,self%fminval(d))
  if (self%id==io%id_debug) then
    m = nint(self%offset)
    print *, 'sink_patch_t%create:  self%pos, d',  self%myposition(m), &
      self%mem(m(1),m(2),m(3),self%idx%d,self%it,1)
    print *,'MM', minval(d), self%fminval(self%idx%d), &
                  maxval(d), self%fmaxval(self%idx%d)
  end if
  call refresh_nbor_list (patch%link, self%task_list)
  call trace%end()
END SUBROUTINE init_sink

!===============================================================================
!> Initialize test instances of the sink_patch_t data type, with parameters for
!> sink creation read from a sink_params namelists.
!===============================================================================
SUBROUTINE ad_hoc_add (self, link)
  class(sink_patch_t)   :: self
  class(link_t), pointer:: link
  !.............................................................................
  class(task_t), pointer:: patch
  class(task_t), pointer:: task
  integer:: iostat, i_sink, m(3)
  real(kind=8):: pos(3)
  type sinks_t
    real(kind=8) :: mass,x,y,z,vx,vy,vz
  end type
  type(sinks_t), allocatable, save:: sinks(:)
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  real(kind=8), save:: mass=1d0, x=0.5d0, y=0.5d0,z, vx=0.5d0, vy=0d0, vz=0d0
  namelist /sink_params/ mass, x, y, z, vx, vy, vz
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
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
        sinks(i_sink)%x = x
        sinks(i_sink)%y = y
        sinks(i_sink)%z = z
        sinks(i_sink)%vx = vx
        sinks(i_sink)%vy = vy
        sinks(i_sink)%vz = vz
      enddo
    endif
  endif
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  ! Go through list of hand made sinks, and create them, if in patch
  !-----------------------------------------------------------------------------
  select type(patch)
  class is (gpatch_t)
    do i_sink=1,n_ad_hoc
      pos = [sinks(i_sink)%x,sinks(i_sink)%y,sinks(i_sink)%z]
      if (patch%contains(pos)) then
        print *, 'adding ad hoc sink to', patch%id, associated(self%task_list)
        allocate (sink_patch_t:: task)
        select type (task)
        class is (sink_patch_t)
          call task%init_sink (patch)
          task%mass = sinks(i_sink)%mass
          task%mom(1) = sinks(i_sink)%mass * sinks(i_sink)%vx
          task%mom(2) = sinks(i_sink)%mass * sinks(i_sink)%vy
          task%mom(3) = sinks(i_sink)%mass * sinks(i_sink)%vz
          call task%accrete
        end select
      endif
    enddo
  end select
  !-----------------------------------------------------------------------------
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
END SUBROUTINE output

!===============================================================================
!> Store a copy of the task list in the data type instance, for use when
!> initializing the nbor list
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(sink_patch_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  call trace%begin ('sink_patch_t%init_task_list')
  self%task_list => task_list
  call trace%end()
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
  self%ds = patch%mesh(:)%d
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
  call trace%end()
END SUBROUTINE init_patch

!===============================================================================
!> Courant condition for sink particles.  FIXME: For now, just copy patch dtime
!===============================================================================
SUBROUTINE courant_condition (self, detailed_timer)
  class(sink_patch_t):: self
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  self%dtime = self%patch%dtime
END SUBROUTINE courant_condition

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
!> Update the sinkparticle task, by accreting mass, and by moving the sink
!===============================================================================
SUBROUTINE update (self)
  class (sink_patch_t):: self
  !.............................................................................
  call trace%begin ('sink_patch_t%update')
  call self%particle_solver%force_field (self)
  call self%accrete
  call self%move
  self%mem(:,:,:,:,self%new,1) = self%mem(:,:,:,:,self%it,1)
  call trace%end()
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
  ! Remove a fraction p of per-volume quantities and keep q=(1-p), where p is
  ! the fraction of an orbit that would be covered by in the time step, assuming
  ! Keplerian rotation. The mask variable makes sure accretion is only done
  ! inside the accretion radius.
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
!> Move the sink particle one time step self%dt, using kick-drift-kick. We do
!> not include any potential here, since both the force from other sinks and the
!> force from the gas is evaluated directly.
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
  class(task_t), pointer:: task, patch
  class(link_t), pointer:: new_link, nbor
  class(sink_patch_t), pointer:: new_sink
  class(dll_node_t), pointer:: particle
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d
  integer:: id, m(3)
  !-----------------------------------------------------------------------------
  if (.not.sink_patch%on) return
  patch => link%task
  select type (patch)
  class is (sink_patch_t)
    !---------------------------------------------------------------------------
    ! Refresh nbor list
    !---------------------------------------------------------------------------
    call self%task_list%refresh_nbors(link)
  class is (gpatch_t)
    check_sinkparticle = sink_patch%check (patch)
    if (check_sinkparticle) then
      allocate (sink_patch_t:: task)
      select type (task)
      class is (sink_patch_t)
        call task%init_sink (patch)
        call task%accrete
      end select
    end if
  end select
END SUBROUTINE refine_check

!-------------------------------------------------------------------------
! Refresh the nbor lists of the parent and its nbors; the parent should
! be done after the nbors, to avoid putting the child link on its nbor
! list until the other ones are updated.  All nbors that happen to have
! overlap with the sink particle patch will have it as an nbor
!-------------------------------------------------------------------------
SUBROUTINE refresh_nbor_list (link, task_list)
  implicit none
  class(link_t), pointer :: link
  class(list_t), pointer :: task_list
  class(link_t), pointer :: nbor
  !
  nbor => link%nbor
  do while (associated(nbor))
    call task_list%init_nbors (nbor%link)
    nbor => nbor%next
  end do
  call task_list%init_nbors (link)
END SUBROUTINE refresh_nbor_list

END MODULE sink_patch_mod
