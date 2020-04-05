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
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE gpatch_mod
  USE mesh_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  implicit none
  private
  type, public, extends(task_t):: sinks_t
    real(8):: mass=0.0, radius=8.0
    real(8):: gamma=1.0
    real(8), allocatable:: pos(:,:)
    real(8):: SN_time=0.0, previous_time=0d0, dt_min=1e-7
    class(gpatch_t), pointer:: patch
  contains
    procedure:: pre_init
    procedure:: init
    procedure:: check_refine
    procedure:: create
    procedure:: accrete
    procedure:: pre_update
    procedure:: post_update
  end type
  real::           v_min = 1.0
  real::     jeans_limit = 4.0
  real::       per_omega = 0.001
  real:: create_fraction = 0.10
  integer:: accretion_cells = 8
  integer::         verbose = 0
  type(sinks_t), public:: sinks
CONTAINS

!===============================================================================
!> Initialize sink task patch connection
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%pre_init')
  self%patch => patch
  self%size = 2d0*accretion_cells*patch%ds
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Initialize sink task patch connection
!===============================================================================
SUBROUTINE init (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%init')
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Check sink creation criterion
!===============================================================================
SUBROUTINE check_refine (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  real, pointer, dimension(:,:,:):: d, e
  real, allocatable:: jeans_resolution(:,:,:)
  real:: ds, g2, min_jeans
  integer:: m(3), ipos(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%check_refine')
  d = patch%mem(:,:,:,patch%idx%d,patch%it,1)
  e = patch%mem(:,:,:,patch%idx%e,patch%it,1)
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
  allocate (jeans_resolution(m(1),m(2),m(3)))
  ds = maxval(patch%mesh%d, mask=patch%mesh%n > 1)
  jeans_resolution = sqrt(math%pi*g2*e/scaling%grav)/(d*ds)
  min_jeans = minval(jeans_resolution)
  if (min_jeans < jeans_limit) then
    if (verbose > 0) &
      write(io_unit%log,*) 'sinks_t%check_refine: creating sink in patch', &
        patch%id, ipos
    ipos = minloc(jeans_resolution)
    call self%create (patch, ipos)
  else if (verbose > 1) then
    write(io_unit%log,*) 'sinks_t%check_refine: id, min_jeans =', &
      patch%id, min_jeans
  end if
  deallocate (jeans_resolution)
  call trace%end ()
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
  m1 => patch%mesh(1)
  m2 => patch%mesh(2)
  m3 => patch%mesh(3)
  r3 = m3%p + m3%r(ipos(3))
  r2 = m2%p + m2%r(ipos(2))
  r1 = m1%p + m1%r(ipos(1))
  self%position = [r1,r2,r3]
  call self%accrete (patch)
  call trace%end ()
END SUBROUTINE create

!===============================================================================
!> Check sink creation criterion
!===============================================================================
SUBROUTINE accrete (self, patch)
  class(sinks_t):: self
  class(gpatch_t), pointer:: patch
  integer:: pos(3)
  !.............................................................................
  real:: r1, r2, r3, ds, dm, rate, omega, fraction, accretion_radius, delta_t
  integer:: i1, i2, i3, m(3)
  real, pointer, dimension(:,:,:):: d, e, px, py, pz
  real, allocatable, dimension(:,:,:):: radius, delta_d
  class(mesh_t), pointer:: m1, m2, m3
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%accrete')
  if (patch%time > self%time .and. self%overlaps(patch)) then
    m = patch%mesh%gn
    allocate (radius(m(1),m(2),m(3)), delta_d(m(1),m(2),m(3)))
    m1 => patch%mesh(1)
    m2 => patch%mesh(2)
    m3 => patch%mesh(3)
    do i3=1,m(3)
      r3 = m3%p + m3%r(i3)
      do i2=1,m(2)
        r2 = m2%p + m2%r(i2)
        do i1=1,m(1)
          r1 = m1%p + m1%r(i1)
          radius = (self%position(1)-r1)**2  &
                 + (self%position(2)-r2)**2  &
                 + (self%position(3)-r3)**2
        end do
      end do
    end do
    !---------------------------------------------------------------------------
    ! Density, energy per unit volume, intertia
    !---------------------------------------------------------------------------
    d  => patch%mem(:,:,:,patch%idx%d ,patch%it,1)
    e  => patch%mem(:,:,:,patch%idx%e ,patch%it,1)
    px => patch%mem(:,:,:,patch%idx%px,patch%it,1)
    py => patch%mem(:,:,:,patch%idx%py,patch%it,1)
    pz => patch%mem(:,:,:,patch%idx%pz,patch%it,1)
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
      delta_t = patch%time - self%time
      fraction =  omega*per_omega*delta_t
    end if
    delta_d = 0.0
    !---------------------------------------------------------------------------
    ! Subtract a fraction (tapered towards the outer edge) of all per-volume
    ! variables
    !---------------------------------------------------------------------------
    where (radius < accretion_radius)
      delta_d = d*fraction*(1-radius/accretion_radius)**2
      d = d - delta_d
      e = e - delta_d*e/d
      px = px - delta_d*px/d
      py = py - delta_d*py/d
      pz = pz - delta_d*pz/d
    end where
    !---------------------------------------------------------------------------
    ! Add up the mass contributions and add to the sink mass
    !---------------------------------------------------------------------------
    dm = product(patch%ds)*sum(delta_d)
    !$omp atomic update
    self%mass = self%mass + dm
    !---------------------------------------------------------------------------
    ! The patch that actually contains the sink is responsible for updating the
    ! time
    !---------------------------------------------------------------------------
    if (patch%contains(self%position)) then
      self%time = patch%time
    end if
  end if
  call trace%end ()
END SUBROUTINE accrete

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE pre_update (self, patch)
  class(sinks_t):: self
  class(patch_t):: patch
  logical:: contains_sink
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%pre_update')
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Post-update sink task -- placing accertion here for now
!===============================================================================
SUBROUTINE post_update (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%post_update')
  call self%accrete (self%patch)
  call trace%end ()
END SUBROUTINE post_update

END MODULE sinks_mod
