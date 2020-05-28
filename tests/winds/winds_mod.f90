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
  USE gpatch_mod
  USE patch_mod
  USE link_mod
  USE sinks_mod
  USE scaling_mod
  USE units_mod
  USE math_mod
  USE parameters_mod, only: smallc
  implicit none
  private
  type, public:: winds_t
    real:: velocity, mdot, ramp_time, energy
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure:: pre_init
    procedure:: pre_update
    procedure:: post_update
    procedure:: force
  end type
  integer::     verbose = 0
  real::   velocity_kms = 2e3
  real::  ramp_time_myr = 0.1
  real::    msun_per_yr = 1e-5
  real::         cs_kms = 0.18
  type(winds_t), public:: winds
CONTAINS

!===============================================================================
!> Pre-init, to pick up a pointer to the patch we belong to.
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(winds_t):: self
  class(gpatch_t), pointer:: patch
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /winds_params/ verbose, velocity_kms, ramp_time_myr, msun_per_yr, &
                          cs_kms
  !-----------------------------------------------------------------------------
  call trace%begin ('winds_t%pre_init')
  self%patch => patch
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, winds_params, iostat=iostat)
    call io%namelist_warning ('winds_params', iostat)
    write (io%namelist, winds_params)
    first_time = .false.
  end if
  call io%lock%unset
  self%velocity = velocity_kms*1e5/scaling%u
  self%mdot = msun_per_yr*cgs%m_sun/scaling%m/(cgs%yr/scaling%t)
  self%ramp_time = ramp_time_myr*cgs%yr*1e6/scaling%t
  self%energy = (cs_kms*cgs%kms/scaling%u)**2/(patch%gamma*(patch%gamma-1d0))
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> The update procedure consists of finding the sink(s) that overlap with our
!> parent task.  The sink tasks are available in the nbor list of the patch
!===============================================================================
SUBROUTINE pre_update (self)
  class(winds_t):: self
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('winds_t%pre_update')
  nbor => self%patch%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (sinks_t)
    call self%force (task)
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> The post_update procedure consists of finding the sink(s) that overlap with
!> our parent task.  The sink tasks are available in the nbor list of the patch
!===============================================================================
SUBROUTINE post_update (self)
  class(winds_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('winds_t%post_update')
  call trace%end ()
END SUBROUTINE post_update

!===============================================================================
!> Compute the wind input to the patch
!===============================================================================
SUBROUTINE force (self, sink)
  class(winds_t):: self
  class(sinks_t):: sink
  class(patch_t), pointer:: patch
  class(task_t), pointer:: task
  real(8):: m_sun               ! solar mass in code units
  real(4):: dt                  ! time step
  real(4):: ramp                ! ramp profile
  real(4):: mdot                ! mass loss rate
  real(4):: r0, d0, v0, e0      ! scaling constants
  real, pointer, dimension(:,:,:):: d, px, py, pz, e
  real, allocatable, dimension(:,:,:):: x, y, z, r, dd_dt, du_dt
  logical, allocatable:: inside(:,:,:)
  integer:: i, i1, i2, i3, m(6), it, n_inside
  real:: ee, ee_min, t
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  call trace%begin ('winds_t%force')
  patch => self%patch
  !-----------------------------------------------------------------------------
  ! Ramp, increasing the forcing gradually with time
  !-----------------------------------------------------------------------------
  t = patch%time/self%ramp_time
  ramp = 1.0-(1.0+t+0.5*t**2)*exp(-t)
  !-----------------------------------------------------------------------------
  ! Fitting formula for the mass loss rate (Paolo)
  !-----------------------------------------------------------------------------
  m_sun = cgs%m_sun/scaling%m
  mdot = ramp*1.19e-10*(sink%mass/m_sun-7.57)**2.60
  !-----------------------------------------------------------------------------
  ! Distances from sink
  !-----------------------------------------------------------------------------
  m = shape(patch%mem)
  allocate (x(m(1),m(2),m(3)))
  allocate (y(m(1),m(2),m(3)))
  allocate (z(m(1),m(2),m(3)))
  allocate (r(m(1),m(2),m(3)))
  allocate (dd_dt(m(1),m(2),m(3)))
  allocate (du_dt(m(1),m(2),m(3)))
  allocate (inside(m(1),m(2),m(3)))
  do i3=1,m(3)
  do i2=1,m(2)
  do i1=1,m(1)
    x(i1,i2,i3) = (patch%mesh(1)%p - sink%position(1)) + patch%mesh(1)%r(i1)
    y(i1,i2,i3) = (patch%mesh(2)%p - sink%position(2)) + patch%mesh(2)%r(i2)
    z(i1,i2,i3) = (patch%mesh(3)%p - sink%position(3)) + patch%mesh(3)%r(i3)
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Scaling constants in code units
  !-----------------------------------------------------------------------------
  r0 = sink%radius*patch%ds(1)
  v0 = velocity_kms*1e5/scaling%u
  d0 = self%mdot/(4.*math%pi*r0**2*v0)
  e0 = d0*(cs_kms*cgs%kms/scaling%u)**2/(patch%gamma*(patch%gamma-1d0))
  !-----------------------------------------------------------------------------
  ! Express the rate of adding mass and radial momentum as functions of r,
  ! leaving one power of r to the final expressions
  !-----------------------------------------------------------------------------
  r = sqrt(x**2+y**2+z**2)/r0
  ramp = 1.0
  dd_dt = ramp*d0*v0/r0*2.*r*(6.-5.*r)
  du_dt = ramp*v0**2/r0*6.*r**2*(1.-r)*(3.-2.*r)
  !print '(a,1p,10e11.2)', &
  !  'wind: radius, r0,v0,d0,dd_dt,du_dt =', &
  !  sink%radius, r0, v0, d0, maxval(dd_dt), maxval(du_dt)
  !-----------------------------------------------------------------------------
  ! If this patch has overlapping points, add mass and radial momentum
  !-----------------------------------------------------------------------------
  inside = (r < 1)
  n_inside = count(inside)
  if (n_inside > 0) then
    !---------------------------------------------------------------------------
    ! Add mass and momentum inside the given radius, letting the velocity approach
    ! a pattern with constant divergence
    !---------------------------------------------------------------------------
    it = patch%new
    d  => patch%mem(:,:,:,patch%idx%d ,it,1)
    e  => patch%mem(:,:,:,patch%idx%e ,it,1)
    px => patch%mem(:,:,:,patch%idx%px,it,1)
    py => patch%mem(:,:,:,patch%idx%py,it,1)
    pz => patch%mem(:,:,:,patch%idx%pz,it,1)
    ee_min = (smallc*cgs%kms/scaling%u)**2/(patch%gamma*(patch%gamma-1d0))
    if (verbose > 1) &
      write (io_unit%output,*) &
        'before:', patch%id, patch%fmaxval(abs(px/d)), maxval(abs(px/d))
    ee = self%energy
    dt = patch%dtime
    if (patch%time == 0d0) dt = 0d0
    if (.not.associated(patch%force_per_unit_mass)) then
      allocate (patch%force_per_unit_mass(m(1),m(2),m(3),3))
      patch%force_per_unit_mass = 0.0
      print *,patch%id, ' force allocated'
    end if
    where (inside)
      d  = d  + dt*dd_dt
      e  = e  + dt*dd_dt*e/d
      px = px + dt*(d*du_dt*x/r0 + dd_dt*px/d0)
      py = py + dt*(d*du_dt*y/r0 + dd_dt*py/d0)
      pz = pz + dt*(d*du_dt*z/r0 + dd_dt*pz/d0)
      !patch%force_per_unit_mass(:,:,:,1) = du_dt*x/r0 + dd_dt*px/(d0*d)
      !patch%force_per_unit_mass(:,:,:,2) = du_dt*y/r0 + dd_dt*py/(d0*d)
      !patch%force_per_unit_mass(:,:,:,3) = du_dt*z/r0 + dd_dt*pz/(d0*d)
    end where
    if (verbose > 1) &
      write (io_unit%output,*) &
        ' after:', patch%id, patch%fmaxval(abs(px/d)), maxval(abs(px/d))
  end if
  if (verbose > 0) then
    !print '("XXX",3(3i4,2x))',shape(e),shape(d),i
    i = patch%n(1)/2
    write (io_unit%output,'(1x,a,i5,i8,i4,1p,2e11.2)') &
      'winds_t: id, n, i, e, T =', patch%id, n_inside, i, e(i,i,i), &
      e(i,i,i)/d(i,i,i)*(patch%gamma-1.0)*scaling%temp
    write (io_unit%log,'(a,1p,11e10.2)') &
      'winds_t:  d1',     d(i:i+10,i,i)
    write (io_unit%log,'(a,1p,11e10.2)') &
      'winds_t: dd1', dd_dt(i:i+10,i,i)
  end if
  deallocate (inside, x, y, z, r, dd_dt, du_dt)
  call trace%end ()
END SUBROUTINE force

END MODULE winds_mod
