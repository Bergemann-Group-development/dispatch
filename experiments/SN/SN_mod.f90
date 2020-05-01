!===============================================================================
!> Module for handling SN from sink particles representing stars
!> As such, this module needs to know about both patches and SN, so the
!> hierachy must be
!>
!>     experiment_t
!>         |
!>      solver_t
!>         |
!>      extras_t
!>       |  \
!>       |  SN_t
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
MODULE SN_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE gpatch_mod
  USE link_mod
  USE sinks_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  real, parameter:: fourpi = 16.*atan(1.)
  type, public:: SN_t
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure:: pre_init
    procedure:: pre_update
    procedure:: courant_condition
  end type
  integer:: verbose = 0
  real::     radius = 5.0           ! initial ejecta radius, in cells
  real:: energy_foe = 1.0           ! SN energy in fiftyone ergs (foe)
  type(sinks_t):: sink
  CONTAINS

!===============================================================================
!> Pre-init, to pick up a pointer to the patch we belong to.
!===============================================================================
SUBROUTINE pre_init (self, patch)
  class(SN_t):: self
  class(gpatch_t), target:: patch
  class(sinks_t), pointer:: sink
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /SN_params/ verbose, radius
  !-----------------------------------------------------------------------------
  call trace%begin ('SN_t%pre_init')
  self%patch => patch
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, SN_params, iostat=iostat)
    call io%namelist_warning ('SN_params', iostat)
    write (io%output, SN_params)
    write (io%namelist, SN_params)
    call sink%init
    first_time = .false.
  end if
  call io%lock%unset
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> The update procedure consists of finding the sink(s) that overlap with our
!> parent task.  The sink tasks are available in the nbor list of the patch
!===============================================================================
SUBROUTINE pre_update (self)
  class(SN_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  class(nbor_t), pointer:: nbor
  class(task_t), pointer:: task
  logical:: just_happened, contains_sink
  !-----------------------------------------------------------------------------
  call trace%begin ('SN_t%pre_update')
  patch => self%patch
  contains_sink = patch%contains (sink%position)
  if (contains_sink) then
    just_happened = sink%previous_time < sink%SN_time &
              .and. sink%time          > sink%SN_time
!print *,patch%id, patch%position
!print *,patch%id, sink%position
!print *,patch%id, sink%previous_time, sink%SN_time, sink%time, just_happened
    if (just_happened) then
      call SN (self)
      sink%previous_time = sink%time
    end if
  end if
  !nbor => patch%link%nbor
  !do while (associated(nbor))
  !  task => nbor%task
  !  select type (task)
  !  class is (sinks_t)
  !  call SN (self, patch, task)
  !  end select
  !  nbor => nbor%next
  !end do
  call trace%end ()
END SUBROUTINE pre_update

!===============================================================================
!> Compute the SN input to the patch
!===============================================================================
SUBROUTINE SN (self)
  class(SN_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  class(task_t), pointer:: task
  real(8):: m_sun               ! solar mass in code units
  real(4):: dm                  ! change in mass per cell
  real(4):: dd                  ! change in mass density
  real(4):: velocity            ! change in mass density
  real(4):: r_code              ! initial SN ejecta radius in code units
  real, pointer, dimension(:,:,:):: d, px, py, pz, e
  real, allocatable, dimension(:,:,:):: x, y, z, r
  logical, allocatable:: inside(:,:,:)
  integer:: i, i1, i2, i3, m(6), n_inside
  real:: ee, ee_min, foe
  !-----------------------------------------------------------------------------
  call trace%begin ('SN_t%SN')
  patch => self%patch
  !-----------------------------------------------------------------------------
  ! Fitting formula for the mass loss rate (Paolo)
  !-----------------------------------------------------------------------------
  m_sun = cgs%m_sun/scaling%m
  !-----------------------------------------------------------------------------
  ! Distances from sink
  !-----------------------------------------------------------------------------
  m = shape(patch%mem)
  allocate (x(m(1),m(2),m(3)))
  allocate (y(m(1),m(2),m(3)))
  allocate (z(m(1),m(2),m(3)))
  allocate (r(m(1),m(2),m(3)))
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
  r = sqrt(x**2+y**2+z**2)
  r_code = radius*patch%ds(1)
  inside = (r < r_code)
  !-----------------------------------------------------------------------------
  ! Check if we have points inside and have not done this SN
  !-----------------------------------------------------------------------------
  print *, '*********************************'
  print *, '** BOOOOOOOOOOOOOOOOOOOOOOOOOM **'
  print *, '*********************************'
  dm = sink%mass/(fourpi/3.*radius**3)        ! mass per cell
  foe = 1d51/scaling%e                        ! 1 foe in code units
  ee = foe/sink%mass                          ! energy per unit mass
  velocity = sqrt(2.0*ee)                     ! ejecta velocity in code units
  dd = dm/product(patch%ds)                   ! density increase per cell
  !---------------------------------------------------------------------------
  ! Add mass and momentum inside the given radius, letting the velocity approach
  ! a pattern with constant divergence
  !---------------------------------------------------------------------------
  d  => patch%mem(:,:,:,patch%idx%d ,patch%new,1)
  e  => patch%mem(:,:,:,patch%idx%e ,patch%new,1)
  px => patch%mem(:,:,:,patch%idx%px,patch%new,1)
  py => patch%mem(:,:,:,patch%idx%py,patch%new,1)
  pz => patch%mem(:,:,:,patch%idx%pz,patch%new,1)
  i = patch%n(1)/2
  print *, 'd: now, incr =',d(i,i,i), dd
  print *, 'e: now, incr =',e(i,i,i)/d(i,i,i), ee
  where (inside)
    e  = e - 0.5*(px**2+py**2+pz**2)/d
    d  = d  + dd
    px = px + dd*(velocity*x/r_code + sink%velocity(1))
    py = py + dd*(velocity*y/r_code + sink%velocity(2))
    pz = pz + dd*(velocity*z/r_code + sink%velocity(3))
    e  = e + 0.5*(px**2+py**2+pz**2)/d
  end where
  if (verbose > 1) then
    n_inside = count(inside)
    write (stdout,'(1x,a,i5,i9,1p,8e11.2)') &
      'SN_t: id, n, m, dd, v, ee, T, vmax =', patch%id, n_inside, &
      sink%mass, dd, velocity, ee, &
      ee*(patch%gamma-1.0)*scaling%temp, &
      patch%fmaxval(abs(px/d)), maxval(abs(px/d))
  end if
  deallocate (inside, x, y, z, r)
  call trace%end ()
END SUBROUTINE SN

!===============================================================================
!> If the next SN time is within reach, make sure to overshoot with exactly
!> sink%dt_min, to land on the same time, and the right side.
!===============================================================================
SUBROUTINE courant_condition (self)
  class(SN_t):: self
  real(8):: next_time
  !-----------------------------------------------------------------------------
  call trace%begin ('SN_t%courant_condition')
  if (self%patch%time < sink%SN_time) then
    self%patch%dtime = min(self%patch%dtime, &
                           sink%SN_time - self%patch%time + sink%dt_min)
  end if
  call trace%end ()
END SUBROUTINE courant_condition

END MODULE SN_mod
