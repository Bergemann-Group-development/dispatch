!===============================================================================
!> Module for handling winds from sink particles representing stars
!> As such, this module needs to know about both patches and sinks, so the
!> hierachy must be
!>
!>     experiment_t
!>         |
!>      solver_t
!>         |
!>      extras_t
!>       |  \
!>       |  winds_t
!>       |  / |
!>     sink_t |
!>         |  |
!>        patch_t
!>         |   \
!>         |  link_t
!>         |   /
!>       task_t
!>
!===============================================================================
MODULE winds_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE sinks_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  type, public:: winds_t
    real:: radius, velocity, mdot, ramp_time
    class(patch_t), pointer:: patch => null()
  contains
    procedure:: init
    procedure:: update
    procedure:: wind
  end type
  real:: velocity_kms = 2e3
  real:: ramp_time_myr = 0.1
  real:: msun_per_myr = 2e-5
  type(winds_t), public:: winds
CONTAINS

!===============================================================================
!> Pre-init, to pick up a pointer to the patch we belong to.
!===============================================================================
SUBROUTINE init (self, patch)
  class(winds_t):: self
  class(patch_t), target:: patch
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /winds_params/ velocity_kms, ramp_time_myr, msun_per_myr
  !-----------------------------------------------------------------------------
  call trace%begin ('winds_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, winds_params, iostat=iostat)
    write (io%output, winds_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)

  self%velocity = velocity_kms*1e5/scaling%u
  self%radius = patch%size(1)*0.4
  self%mdot = msun_per_myr*cgs%m_sun/scaling%m/(cgs%yr/scaling%t)
  self%ramp_time = ramp_time_myr*cgs%yr*1e6/scaling%t
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> The update procedure consists of finding the sink(s) that overlap with our
!> parent task.  The sink tasks are available in the nbor list of the patch
!===============================================================================
SUBROUTINE update (self, patch)
  class(winds_t):: self
  class(patch_t):: patch
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('winds_t%update')
  nbor => patch%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (sinks_t)
    call self%wind (patch, task)
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Compute the wind input to the patch
!===============================================================================
SUBROUTINE wind (self, patch, sink)
  class(winds_t):: self
  class(patch_t):: patch
  class(sinks_t):: sink
  class(task_t), pointer:: task
  real(8):: m_sun               ! solar mass in code units
  real(4):: dm                  ! change in mass per cell
  real(4):: dd                  ! change in mass density
  real(4):: ramp                ! ramp profile
  real(4):: mdot                ! mass loss rate
  real, pointer, dimension(:,:,:):: d, px, py, pz
  real, allocatable, dimension(:,:,:):: x, y, z, r
  logical, allocatable:: inside(:,:,:)
  integer:: i1, i2, i3, m(6)
  !-----------------------------------------------------------------------------
  call trace%begin ('winds_t%wind')
  !-----------------------------------------------------------------------------
  ! Fitting formula for the mass loss rate (Paolo)
  !-----------------------------------------------------------------------------
  m_sun = cgs%m_sun/scaling%m
  ramp = 1.0-exp(-patch%time/self%ramp_time)
  mdot = log10(sink%mass/(10.*m_sun))**3*self%mdot*ramp
  !-----------------------------------------------------------------------------
  ! Distances from sink
  !-----------------------------------------------------------------------------
  m = shape(patch%mem)
  allocate (x(m(1),m(2),m(3)))
  allocate (y(m(1),m(2),m(3)))
  allocate (z(m(1),m(2),m(3)))
  allocate (r(m(1),m(2),m(3)))
  do i3=1,m(3)
  do i2=1,m(2)
  do i1=1,m(1)
    x(i1,i2,i3) = (patch%mesh(1)%p - sink%position(1)) + patch%mesh(1)%r(i1)
    y(i1,i2,i3) = (patch%mesh(2)%p - sink%position(2)) + patch%mesh(2)%r(i2)
    z(i1,i2,i3) = (patch%mesh(3)%p - sink%position(3)) + patch%mesh(3)%r(i3)
  end do
  end do
  end do
  r = sqrt(x**2+y**2+z**2)
  !-----------------------------------------------------------------------------
  ! Compute the change in mass in the timestep, and distribute over the cells
  !-----------------------------------------------------------------------------
  allocate (inside(m(1),m(2),m(3)))
  inside = (r < self%radius)
  dm = mdot*patch%dtime
  dm = dm/count(inside)
  dd = dm/product(patch%ds)
  !-----------------------------------------------------------------------------
  ! Add mass and momentum inside the given radius, letting the velocity approach
  ! a pattern with constant divergence
  !-----------------------------------------------------------------------------
  d  => patch%mem(:,:,:,patch%idx%d ,patch%new,1)
  px => patch%mem(:,:,:,patch%idx%px,patch%new,1)
  py => patch%mem(:,:,:,patch%idx%py,patch%new,1)
  pz => patch%mem(:,:,:,patch%idx%pz,patch%new,1)
  where (inside)
    px = px + dd*(self%velocity*x/self%radius-px/d)
    py = py + dd*(self%velocity*y/self%radius-py/d)
    pz = pz + dd*(self%velocity*z/self%radius-pz/d)
    d  = d  + dd
  end where
  deallocate (inside, x, y, z, r)
  call trace%end ()
END SUBROUTINE wind

END MODULE winds_mod
