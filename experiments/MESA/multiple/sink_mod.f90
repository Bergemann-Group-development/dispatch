!===============================================================================
!> Simplified sink particle data type. The task makes its presence known
!> by being an nbor task to patches that it overlaps with.  When they update,
!> and encounter an nbor whose data type is sink_t, it should call the
!> sink_t%accrete procedure, which computes the mass accreted from that
!> task, and adds it to the sink_t data type. The patch_t task gets access
!> to the sink_t%accrete procedure via the pointer to the nbor task, and there
!> is thus no need for the sink_t task to compute any accretion.  Its %update
!> procedure should instead just update its position, and compute whatever
!> consequence it may have (e.g. using MESA).
!>
!> Since there are only a few scalar properties that need to be updated by the
!> patch nbors, this can be handled by atomic updates, without task locking.
!>
!> The sink_t data type should eventually have the following procedures:
!>
!>    * accrete(): called by nbor patches to compute accretion onto the sink
!>    *   force(): called by nbor patches to compute force on their gas, and
!>                 called by other sink particles to compute the force on them
!>    *  update(): update the sinkparticle position, after computing forces
!>    *    wind(): computes the wind output from the sink particle
!>    * explode(): computes the momentum input from the SN explosion of the star
!===============================================================================
MODULE sink_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE gpatch_mod
  USE link_mod
  USE list_mod
  USE mesh_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE aux_mod
  implicit none
  private
  type, public:: particle_t
    type(particle_t), pointer:: next => null()  ! next particle
    type(particle_t), pointer:: prev => null()  ! previous particle
    integer(8)      :: id = 0                   ! particle id
    real(8)         :: t = 0.0_8                ! time
    real(8)         :: dt                       ! time step
    real            :: v(3)                     ! velocity (cell units / time)
    real            :: dv(3)                    ! vgas - velocity
    real            :: a(3) = 0.0               ! acceleration
    real            :: d(3)                     ! drag
    real            :: r(3)                     ! real positions (cell units)
    integer(kind=2) :: q(3)                     ! integer position (cell units)
  end type
  type, public, extends(task_t):: sink_t
    real(8):: mass=0d0, radius=8.0
    real(8):: px=0d0, py=0d0, pz=0d0
    real(8):: gamma=1.0
    real, pointer:: mass_p=>null(), vel_p(:)=>null(), pos_p(:)=>null()
    real(8), allocatable:: pos(:,:)
    real(8):: SN_time=0.0, previous_time=0d0, dt_min=1e-7
    class(gpatch_t), pointer:: patch => null()
    class(list_t), pointer:: task_list => null()
    class(link_t), pointer:: link => null()
    type(particle_t):: particle
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: init_task_list
    procedure:: check
    procedure:: no_overlap
    procedure:: contains_point
    procedure:: create
    procedure:: add_task
    procedure:: accrete
    procedure:: pre_update
    procedure:: update
    procedure:: post_update
    procedure:: output
    procedure:: force
  end type
  logical::           on = .true.     ! on by default
  real::     jeans_limit = 4.0        ! create at jeans_limit cells per Jean's L 
  real::       per_omega = 0.001      ! fraction of mass per radian orbital
  real:: create_fraction = 0.01       ! fraction of mass to accrete initially
  integer:: accretion_cells = 8       ! accretion radius in cells
  integer::         verbose = 0       ! information level
  !type(sink_t), public:: sink
CONTAINS

!===============================================================================
!> Initialize sink task patch connection
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(sink_t):: self
  class(gpatch_t), pointer:: patch
  integer:: iostat
  namelist /sink_params/ jeans_limit, accretion_cells, create_fraction, &
    per_omega, verbose
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%pre_init')
  self%patch => patch
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, sink_params, iostat=iostat)
    call io%namelist_warning ('sink_params', iostat)
    write (io%namelist, sink_params)
    first_time = .false.
  end if
  call io%lock%unset
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize sink task patch connection
!===============================================================================
SUBROUTINE init (self)
  class(sink_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%init')
  self%size = 2d0*accretion_cells*self%patch%ds
  self%position = huge(1.0)
  self%nt      = self%patch%nt
  self%time    = self%patch%time
  self%level   = self%patch%level
  self%quality = self%patch%quality
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Save a pointer to the task list
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(sink_t), target:: self
  class(list_t), pointer:: task_list
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%init_task_list')
  self%task_list => task_list
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Check sink creation criterion, which is to have the local Jeans' length
!> resolved by at least "jeans_limit" cells.  The Jeans length is defined as
!> L = c_sound/free_fall_time, where free_fall_time = 1/sqrt(G*rho), with G
!> being the (scaled) constant of gravity.  The sound speed squared is approx
!> $c_sound^2 = \Gamma(\Gamma-1) e / d$, where $\Gamma$ is an ideal gas
!> gamma, d is density, and e is thermal energy per unit volume.
!===============================================================================
SUBROUTINE check (self, patch)
  class(sink_t):: self
  class(gpatch_t), pointer:: patch
  real, pointer, dimension(:,:,:):: d, e
  real, allocatable:: jeans_resolution(:,:,:)
  type(sink_t), pointer:: new
  real:: ds, g2, min_jeans
  integer:: m(3), ipos(3), l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%check', itimer=itimer)
  l = patch%mesh%li
  u = patch%mesh%ui
  d => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d,patch%it,1)
  e => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%e,patch%it,1)
  if (self%gamma==1d0) then
    g2 = 1d0
  else
    g2 = self%gamma*(self%gamma-1d0)
  end if
  !-----------------------------------------------------------------------------
  ! Compute the ratio of Jeans' length to cell size, where cell size is the
  ! largest of the mesh sizes in directions that are not collapsed
  !-----------------------------------------------------------------------------
  m = patch%mesh%gn
  allocate (jeans_resolution(l(1):u(1),l(2):u(2),l(3):u(3)))
  ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
  jeans_resolution = sqrt(math%pi*g2*e/scaling%grav)/(d*ds)
  min_jeans = minval(jeans_resolution)
  if (min_jeans < jeans_limit) then
    ipos = minloc(jeans_resolution)
    if (self%no_overlap (ipos+1)) then
      !...........................................................................
      if (verbose > 0) &
        write(io_unit%log,'(1x,a,i6,2x,3i4,f6.1)') &
        'sink_t%check: creating sink id, ipos, jeans =', &
          patch%id, ipos, jeans_resolution(ipos(1),ipos(2),ipos(3))
      !...........................................................................
      allocate (new, source=self)
      new%patch => patch
      new%task_list => self%task_list
      call new%create (patch, ipos+l)
    else if (verbose > 1) then
      !...........................................................................
      write(io_unit%mpi,'(1x,a,i6,f10.3,f12.6,1p,e11.2)') &
      'sink_t%check: id, min_jeans, time, max(d) =', &
        patch%id, min_jeans, patch%time, maxval(d)
      flush (io_unit%mpi)
      !...........................................................................
    end if
  end if
  deallocate (jeans_resolution)
  call trace%end (itimer)
END SUBROUTINE check

!===============================================================================
!> Check that there is no overlap of the sink accretion sphere with a point
!===============================================================================
FUNCTION no_overlap (self, ipos) RESULT (ok)
  class(sink_t):: self
  integer:: ipos(3)
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: nbtask
  logical:: ok
  !-----------------------------------------------------------------------------
  if (verbose > 1) &
    write (stdout,'(1x,a,i6,2x,3i3)') &
    'testing for overlapping sinks, ID, pos =', self%patch%id, ipos
  ok = .true.
  nbor => self%patch%link%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (sink_t)
      if (nbtask%contains_point (ipos)) then
        if (verbose > 1) &
          write (stdout,*) &
          'sink tasks ID =', nbtask%id,' contains the point'
        ok = .false.        
      else
        if (verbose > 1) &
          write (stdout,*) &
          'sink tasks ID =', nbtask%id,' does not contain the point'
      end if
    class is (patch_t)
      !write (stderr,*) 'patch nbor', nbtask%id
    class default
      write (stderr,*) &
        'sink_t%no_overlap: WARNING, other nbor', nbtask%id
    end select
    nbor => nbor%next
  end do
  if (ok .and. verbose > 1) &
    write (stdout,*) 'OK to go head with creating a sink'
END FUNCTION no_overlap

!===============================================================================
!> Check that there is no overlap of the sink accretion sphere with a point
!===============================================================================
FUNCTION contains_point (self, ipos)
  class(sink_t):: self
  integer:: ipos(3)
  logical:: contains_point
  !-----------------------------------------------------------------------------
  contains_point = sum(real(ipos-self%particle%q)**2) < accretion_cells**2
  if (verbose > 2) &
    write (stdout,'(1x,a,i6,2(2x,3i4),l4)') &
    'contains_point:', self%id, self%particle%q,  ipos, contains_point
END FUNCTION contains_point

!===============================================================================
!> create a sink particle at the give location in the patch
!===============================================================================
SUBROUTINE create (self, patch, ipos)
  class(sink_t):: self
  class(gpatch_t), pointer:: patch
  integer:: ipos(3)
  !.............................................................................
  real:: r1, r2, r3
  integer:: i1, i2, i3, m(3)
  class(mesh_t), pointer:: m1, m2, m3
  real, allocatable, dimension(:,:,:):: delta_d
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%create')
  self%patch => patch
  !$omp critical (create_cr)
  if (verbose > 1) write (stdout,*) 'initializing sink'
  call self%init
  if (verbose > 1) write (stdout,*) 'initializing sink task'
  call self%task_t%init                           ! give the task an ID
  self%particle%q = ipos
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  r1 = m1%p + m1%r(ipos(1))
  r2 = m2%p + m2%r(ipos(2))
  r3 = m3%p + m3%r(ipos(3))
  self%position = [r1,r2,r3]
  self%size = 2d0*accretion_cells*patch%ds
  self%ds = 0d0
  self%time = patch%time-patch%dtime      ! make sure it accretes
  !.............................................................................
  if (verbose > 1) then
    write (stdout,'(1x,a,i6,3i4,1p,5e11.2)') &
      'id, ipos, time, position =',self%id,ipos,self%time,self%position
    write (stdout,*) 'sink_t%create: creating sink'
  end if
  write (io_unit%log,'(1x,a,i6,2x,3i4,3f10.4)') &
     'sink_t%create: id, ipos, position =', patch%id, ipos, self%position
  !.............................................................................
  call self%accrete (self%patch)
  if (verbose > 0) write (stdout,'(1x,a,i6,1p,e11.2)') &
    'sink_t%create: ID, initial mass =', self%id, self%mass*scaling%m/cgs%m_sun
  call self%add_task
  !$omp end critical (create_cr)
  call trace%end ()
END SUBROUTINE create

!===============================================================================
!> Create a task link, with a pointer to this sink task, and append it to the
!> task list.  As a test, set a fixed dtime, and no nbor list, and check what
!> happens to the task.
!===============================================================================
SUBROUTINE add_task (self)
  class(sink_t), target:: self
  class(link_t), pointer:: link
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%add_task')
  !-----------------------------------------------------------------------------
  ! Create a task link, and initialize the task and link pointers
  !-----------------------------------------------------------------------------
  allocate (link)
  link%task => self
  self%link => link
  !-----------------------------------------------------------------------------
  ! Create an nbor node, point it to the patch task, and use it as the nbor list
  ! The sink task should not perform a dnload on the patch nbor.
  !-----------------------------------------------------------------------------
  allocate (nbor)
  !call nbor%init
  nbor%download = .false.
  nbor%task => self%patch
  nbor%link => self%patch%link
  link%nbor => nbor
  !-----------------------------------------------------------------------------
  ! Also add the sink task to the patch nbor list.  TODO: The patch can handle
  ! accretion to the sinks it contains via downld operation!
  !-----------------------------------------------------------------------------
  allocate (nbor)
  !call nbor%init
  nbor%download = .true.
  nbor%task => self
  nbor%link => self%link
  call self%patch%link%add_nbor_by_rank (self%patch%link%nbor, nbor)
  if (verbose > 1) call self%patch%link%print_nbors (' patch')
  !-----------------------------------------------------------------------------
  ! Append the task to the task list
  !-----------------------------------------------------------------------------
  call self%task_list%append (link)
  if (verbose > 1) then
    write (stdout,*) 'sink_t%add_task: appending sink task to task listø'
  end if
  call trace%end ()
END SUBROUTINE add_task

!===============================================================================
!> Accrete to sink particle from patch
!===============================================================================
SUBROUTINE accrete (self, patch)
  class(sink_t):: self
  class(gpatch_t), pointer:: patch
  integer:: pos(3)
  !.............................................................................
  real:: ds, dm, dpx, dpy, dpz, rate, omega, fraction, &
    accretion_radius, accretion_radius2, velocity(3)
  real(8):: r1, r2, r3, dr(3)
  integer:: i1, i2, i3, m(3), nm, l(3), u(3)
  real, pointer, dimension(:,:,:):: d, e, px, py, pz
  real, allocatable, dimension(:,:,:):: radius2, delta_d
  integer(2), allocatable, dimension(:,:,:):: inside
  class(mesh_t), pointer:: m1, m2, m3
  logical:: overlaps
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%accrete', itimer=itimer)
  overlaps = self%overlaps(patch)
  if (verbose > 1) write (stdout,*) &
    'accreting', self%id, patch%id, overlaps, patch%time-self%time
  if (patch%time > self%time .and. overlaps) then
    !...........................................................................
    if (verbose > 2) then
      write (io_unit%log,'(1x,a,i4,l3,4(3f8.3,2x))') &
       'sink_t%accrete: id, over, sink%size, pos, patsh%size, pos =', &
        patch%id, overlaps, self%size, self%position, patch%size, patch%position
    end if
    !...........................................................................
    m = patch%mesh%gn
    l = patch%mesh%li
    u = patch%mesh%ui
    allocate (radius2(l(1):u(1),l(2):u(2),l(3):u(3)), &
              delta_d(l(1):u(1),l(2):u(2),l(3):u(3)), &
              inside (l(1):u(1),l(2):u(2),l(3):u(3)))
    m1 => patch%mesh(1)
    m2 => patch%mesh(2)
    m3 => patch%mesh(3)
    do i3=l(3),u(3)
      r3 = m3%p + m3%r(i3)
      do i2=l(2),u(2)
        r2 = m2%p + m2%r(i2)
        do i1=l(1),u(1)
          r1 = m1%p + m1%r(i1)
          !---------------------------------------------------------------------
          ! Taking periodic boundaries into account in task_t%point_distance()
          !---------------------------------------------------------------------
          dr = self%point_distance (self%position, [r1,r2,r3])
          radius2(i1,i2,i3) = dr(1)**2 + dr(2)**2 + dr(3)**2
        end do
      end do
    end do
    !---------------------------------------------------------------------------
    ! Density, energy per unit volume, intertia
    !---------------------------------------------------------------------------
    d  => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%d ,patch%it,1)
    e  => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%e ,patch%it,1)
    px => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%px,patch%it,1)
    py => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%py,patch%it,1)
    pz => patch%mem(l(1):u(1),l(2):u(2),l(3):u(3),patch%idx%pz,patch%it,1)
    !---------------------------------------------------------------------------
    ! Accretion radius in code units, and accretion rate, scaled by rotation
    ! rate.  The accretion rate is expressed as a fraction times the rotation
    ! rate at the accretion radius, except a create_fraction is used when the
    ! mass is zero.
    !---------------------------------------------------------------------------
    ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
    accretion_radius = accretion_cells*ds
    if (self%mass == 0d0) then
      fraction = create_fraction
    else
      omega = sqrt(scaling%grav*self%mass/accretion_radius**3)
      fraction =  omega*per_omega*patch%dtime
    end if
    delta_d = 0.0
    inside = 0
    !---------------------------------------------------------------------------
    ! Subtract a fraction (tapered towards the outer edge) of all per-volume
    ! variables
    !---------------------------------------------------------------------------
    accretion_radius2 = accretion_radius**2
    where (radius2 < accretion_radius2)
      inside = 1
      delta_d = fraction*(1-radius2/accretion_radius2)**2
      d = d - delta_d*d
      e = e - delta_d*e
      px = px - delta_d*px
      py = py - delta_d*py
      pz = pz - delta_d*pz
    end where
    !---------------------------------------------------------------------------
    ! Add up the mass contributions and add to the sink mass
    !---------------------------------------------------------------------------
    nm  = sum(inside)
    dm  = product(patch%ds)*sum(delta_d*d)
    dpx = product(patch%ds)*sum(delta_d*px)
    dpy = product(patch%ds)*sum(delta_d*py)
    dpz = product(patch%ds)*sum(delta_d*pz)
    !$omp atomic update
    self%mass = self%mass + dm
    !$omp atomic update
    self%px = self%px + dpx
    !$omp atomic update
    self%py = self%py + dpy
    !$omp atomic update
    self%pz = self%pz + dpz
    if (self%mass > tiny(1.0)) then
      self%velocity = [self%px,self%py,self%pz]/self%mass
    else
      self%velocity = 0d0
    end if
    !...........................................................................
    if (verbose > 0) then
      dr = self%point_distance (self%position, patch%position)
      write (io_unit%log,'(1x,a,i4,i8,3(3f8.3,2x),1p,4e10.2)') &
       'sink_t%accrete: id, nm, pos, rpos, vel, m, dm, d_max =', &
        patch%id, nm, self%position, dr/patch%size, self%velocity, &
        self%mass, dm, maxval(d)
    end if
    !...........................................................................
    deallocate (radius2, delta_d, inside)
  else if (verbose > 2) then
    !...........................................................................
    write (io_unit%log,'(1x,a,i6,l4,5(3f10.6,1x))') 'sink_t%accrete: skip', &
      patch%id, overlaps, patch%time, self%time, self%mass, &
      self%position, self%size, patch%position, patch%size
    !...........................................................................
  end if
  call trace%end (itimer)
END SUBROUTINE accrete

!===============================================================================
!> Update patch that may contain sink tasks
!===============================================================================
SUBROUTINE pre_update (self)
  class(sink_t):: self
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: nbtask
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%pre_update')
  call self%check (self%patch)
  nbor => self%patch%link%nbor
  do while (associated(nbor))
    nbtask => nbor%task
    select type (nbtask)
    class is (sink_t)
      if (verbose > 0) write (stdout,*)  &
        self%patch%id, 'found a sink_t task: ID, mass =', &
        nbtask%id, nbtask%mass*scaling%m/cgs%m_sun
      call nbtask%accrete (self%patch)
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE update (self)
  class(sink_t):: self
  real:: velocity(3), dr(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%update')
  if (self%mass > tiny(1.0)) then
    self%dtime = self%patch%dtime
    if (self%mass > tiny(1.0)) then
      velocity = [self%px,self%py,self%pz]/self%mass
    else
      velocity = 0d0
    end if
    self%position = self%position + self%dtime*velocity
    self%particle%q = (self%position-self%patch%position)/self%patch%ds &
                    + self%patch%mesh%o
    if (verbose > 0) &
      write (stdout,'(1x,a,i6,i8,1p,e10.2,3(2x,3e10.2),2x,3i3)') &
      'updating sink: id, mass, times, velocity, pos', &
      self%id, self%istep, self%mass, self%time, self%patch%time, &
      self%patch%dtime, velocity, self%position, self%particle%q
    call self%rotate
    !...........................................................................
    if (verbose > 0) &
      write (io_unit%log,'(1x,a,i6,3i4,1p,3(2x,3e12.2))') &
       'sink_t%post_update: id, position, velocity, mass =', &
        self%patch%id, self%particle%q, self%position, velocity, self%mass
    !...........................................................................
  else if (verbose > 1) then
    dr = self%point_distance (self%position, self%patch%position)
    write (io_unit%log,'(1x,a,i6,1p,3(2x,1p,3e10.2))') &
     'sink_t%post_update: id, pos, (sink-patch)/size =', self%patch%id, &
      self%position, dr/self%patch%size
  end if
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Post-update sink task -- placing accretion here for now
!===============================================================================
SUBROUTINE post_update (self, patch)
  class(sink_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%post_update')
  call self%check (patch)
  call self%accrete (patch)
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE output (self)
  class(sink_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%output')
  if (.not.associated(self%mass_p)) then
    allocate (self%mass_p, self%pos_p(3), self%vel_p(3))
    self%mass_p = self%mass
    self%pos_p = self%position
    self%vel_p = self%velocity
  end if
  call aux%register ('sink_m', self%mass_p)
  call aux%register ('sink_p', self%pos_p)
  call aux%register ('sink_v', self%vel_p)
  call trace%end ()
END SUBROUTINE output

!===============================================================================
!> Force towards sink position
!===============================================================================
FUNCTION force (self, position)
  class(sink_t):: self
  real(8):: position(3)
  real:: force(3), dist(3), dist2
  !-----------------------------------------------------------------------------
  call trace%begin ('sink_t%force')
  dist = self%patch%point_distance (self%position, position)
  dist2 = sum(dist**2)
  force = self%mass*scaling%grav*dist/dist2/sqrt(dist2)
  call trace%end ()
END FUNCTION force

END MODULE sink_mod
