!===============================================================================
!> Simplified sink particle data type. The task makes its presence known
!> by being an nbor task to patches that it overlaps with.  When they update,
!> and encounter an nbor whose data type is sinks_t, it should call the
!> sinks_t%accrete procedure, which computes the mass accreted from that
!> task, and adds it to the sinks_t data type. The patch_t task gets access
!> to the sinks_t%accrete procedure via the pointer to the nbor task, and there
!> is thus no need for the sinks_t task to compute any accretion.  Its %update
!> procedure should instead just update its position, and compute whatever
!> consequence it may have (e.g. using MESA).
!>
!> Since there are only a few scalar properties that need to be updated by the
!> patch nbors, this can be handled by atomic updates, without task locking.
!>
!> The sink_t data type needs to have the following procedures:
!>
!>    * accrete(): called by nbor patches to compute accretion onto the sink
!>    *   force(): called by nbor patches to compute force on their gas, and
!>                 called by other sink particles to compute the force on them
!>    *  update(): update the sinkparticle position, after computing forces
!>    *    wind(): computes the wind output from the sink particle
!>    * explode(): computes the momentum input from the SN explosion of the star
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE gpatch_mod
  USE mesh_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE aux_mod
  implicit none
  private
  type, public, extends(task_t):: sinks_t
    real(8):: mass=0d0, radius=8.0
    real(8):: px=0d0, py=0d0, pz=0d0
    real(8):: gamma=1.0
    real, pointer:: mass_p=>null(), vel_p(:)=>null(), pos_p(:)=>null()
    real(8), allocatable:: pos(:,:)
    real(8):: SN_time=0.0, previous_time=0d0, dt_min=1e-7
    class(gpatch_t), pointer:: patch
    real(8):: dm_accreted=0, dm_wind=0, mdot_wind=0 ! mass change due to wind and accretion. Update when update mass.
    real(8):: v_wind=tiny(0d0)                      ! v_inf of wind
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: check_refine
    procedure:: create
    procedure:: accrete
    procedure:: pre_update
    procedure:: update
    procedure:: post_update
    procedure:: output
  end type
  real::     jeans_limit = 4.0        ! create at jeans_limit cells per Jean's L 
  real::       per_omega = 0.001      ! fraction of mass per radian orbital
  real:: create_fraction = 0.01       ! fraction of mass to accrete initially
  integer:: accretion_cells = 8       ! accretion radius in cells
  integer::         verbose = 0       ! information level
  type(sinks_t), public:: sinks
CONTAINS

!===============================================================================
!> Initialize sink task patch connection
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  integer:: iostat
  namelist /sinks_params/ jeans_limit, accretion_cells, create_fraction, &
    per_omega, verbose
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%pre_init')
  self%patch => patch
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, sinks_params, iostat=iostat)
    call io%namelist_warning ('sinks_params', iostat)
    write (io%output, sinks_params)
    write (io%namelist, sinks_params)
    first_time = .false.
  end if
  call io%lock%unset
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize sink task patch connection
!===============================================================================
SUBROUTINE init (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%init')
  self%size = 2d0*accretion_cells*self%patch%ds
  self%position = huge(1.0)
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Check sink creation criterion, which is to have the local Jeans' length
!> resolved by at least "jeans_limit" cells.  The Jeans length is defined as
!> L = c_sound/free_fall_time, where free_fall_time = 1/sqrt(G*rho), with G
!> being the (scaled) constant of gravity.  The sound speed squared is approx
!> $c_sound^2 = \Gamma(\Gamma-1) e / d$, where $\Gamma$ is an ideal gas
!> gamma, d is density, and e is thermal energy per unit volume.
!===============================================================================
SUBROUTINE check_refine (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  real, pointer, dimension(:,:,:):: d, e
  real, allocatable:: jeans_resolution(:,:,:)
  real:: ds, g2, min_jeans
  integer:: m(3), ipos(3), l(3), u(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%check_refine', itimer=itimer)
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
  if (min_jeans < jeans_limit .and. self%mass == 0d0) then
    ipos = minloc(jeans_resolution)
    !...........................................................................
    if (verbose > 0) &
      write(io_unit%log,'(1x,a,i6,2x,3i4,f6.1)') &
       'sinks_t%check_refine: creating sink id, ipos, jeans =', &
        patch%id, ipos, jeans_resolution(ipos(1),ipos(2),ipos(3))
    !...........................................................................
    call self%create (patch, ipos+l)
  else if (verbose > 1) then
    !...........................................................................
    write(io_unit%log,'(1x,a,i6,f10.3,f12.6,1p,e11.2)') &
     'sinks_t%check_refine: id, min_jeans, time, max(d) =', &
      patch%id, min_jeans, patch%time, maxval(d)
    !...........................................................................
  end if
  deallocate (jeans_resolution)
  call trace%end (itimer)
END SUBROUTINE check_refine

!===============================================================================
!> Check sink creation criterion
!===============================================================================
SUBROUTINE create (self, patch, ipos)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  integer:: ipos(3)
  !.............................................................................
  real:: r1, r2, r3
  integer:: i1, i2, i3, m(3)
  class(mesh_t), pointer:: m1, m2, m3
  real, allocatable, dimension(:,:,:):: delta_d
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%create')
  !$omp critical (create_cr)
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  r1 = m1%p + m1%r(ipos(1))
  r2 = m2%p + m2%r(ipos(2))
  r3 = m3%p + m3%r(ipos(3))
  self%position = [r1,r2,r3]
  self%size = 2d0*accretion_cells*patch%ds
  self%ds = 0d0
  self%time = patch%time-patch%dtime
  !.............................................................................
  if (verbose > 0) &
    write (io_unit%log,'(1x,a,i6,2x,3i4,3f10.4)') &
     'sinks_t%create: id, ipos, position =', patch%id, ipos, self%position
  !.............................................................................
  !$omp end critical (create_cr)
  call trace%end ()
END SUBROUTINE create

!===============================================================================
!> Accrete to sink particle from patch
!===============================================================================
SUBROUTINE accrete (self, patch)
  class(sinks_t):: self
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
  call trace%begin ('sinks_t%accrete', itimer=itimer)
  overlaps = self%overlaps(patch)
  if (patch%time > self%time .and. overlaps) then
    !...........................................................................
    if (verbose > 2) then
      write (io_unit%log,'(1x,a,i4,l3,4(3f8.3,2x))') &
       'sinks_t%accrete: id, over, sink%size, pos, patsh%size, pos =', &
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
    self%dm_accreted = self%dm_accreted + dm
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
       'sinks_t%accrete: id, nm, pos, rpos, vel, m, dm, d_max =', &
        patch%id, nm, self%position, dr/patch%size, self%velocity, &
        self%mass, dm, maxval(d)
    end if
    !...........................................................................
    deallocate (radius2, delta_d, inside)
  else if (verbose > 2) then
    !...........................................................................
    write (io_unit%log,'(1x,a,i6,l4,5(3f10.6,1x))') 'sinks_t%accrete: skip', &
      patch%id, overlaps, patch%time, self%time, self%mass, &
      self%position, self%size, patch%position, patch%size
    !...........................................................................
  end if
  call trace%end (itimer)
END SUBROUTINE accrete

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE pre_update (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%pre_update')
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE update (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%update')
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Post-update sink task -- placing accretion here for now
!===============================================================================
SUBROUTINE post_update (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  real:: velocity(3), dr(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%post_update')
  call self%check_refine (patch)
  call self%accrete (patch)
  !-----------------------------------------------------------------------------
  ! The patch that actually contains the sink is responsible for updating the
  ! time and position of the sink
  !-----------------------------------------------------------------------------
  if (self%mass > tiny(1.0) .and. patch%contains(self%position)) then
    self%dtime = patch%time - self%time
    self%time = patch%time
    if (self%mass > tiny(1.0)) then
      velocity = [self%px,self%py,self%pz]/self%mass
    else
      velocity = 0d0
    end if
    self%position = self%position + self%dtime*velocity
    !...........................................................................
    if (verbose > 0) &
      write (io_unit%log,'(1x,a,i6,1p,3(2x,3e12.2))') &
       'sinks_t%post_update: id, position, velocity, mass =', &
        patch%id, self%position, velocity, self%mass
    !...........................................................................
  else if (verbose > 1) then
    dr = self%point_distance (self%position, patch%position)
    write (io_unit%log,'(1x,a,i6,1p,3(2x,1p,3e10.2))') &
     'sinks_t%post_update: id, pos, (sink-patch)/size =', patch%id, &
      self%position, dr/patch%size
  end if
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Output sink task
!===============================================================================
SUBROUTINE output (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  !todo: fix.
  return
  call trace%begin ('sinks_t%output')
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


END MODULE sinks_mod
