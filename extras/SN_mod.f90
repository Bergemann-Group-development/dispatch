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
  USE patch_mod
  USE link_mod
  USE sinks_mod
  USE scaling_mod
  USE units_mod
  USE parameters_mod, only: smallc
  implicit none
  private
  type, public:: SN_t
    real:: velocity, mdot, energy
    class(patch_t), pointer:: patch => null()
  contains
    procedure:: init
    procedure:: update
    procedure:: SN
  end type
  integer::     verbose = 0
  real::   velocity_kms = 2e3
  real::    msun_per_yr = 1e-5
  real::         cs_kms = 0.18
  logical:: incremental = .false.
CONTAINS

!===============================================================================
!> Pre-init, to pick up a pointer to the patch we belong to.
!===============================================================================
SUBROUTINE init (self, patch)
  class(SN_t):: self
  class(patch_t), target:: patch
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /SN_params/ verbose, velocity_kms, msun_per_yr, &
                          cs_kms, incremental
  !-----------------------------------------------------------------------------
  call trace%begin ('SN_t%init')
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, SN_params, iostat=iostat)
    call io%namelist_warning ('SN_params', iostat)
    write (io%output, SN_params)
    write (io%namelist, SN_params)
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  self%velocity = velocity_kms*1e5/scaling%u
  self%mdot = msun_per_yr*cgs%m_sun/scaling%m/(cgs%yr/scaling%t)
  self%energy = (cs_kms*cgs%kms/scaling%u)**2/(patch%gamma*(patch%gamma-1d0))
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> The update procedure consists of finding the sink(s) that overlap with our
!> parent task.  The sink tasks are available in the nbor list of the patch
!===============================================================================
SUBROUTINE update (self, patch)
  class(SN_t):: self
  class(patch_t):: patch
  class(link_t), pointer:: nbor
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('SN_t%update')
  nbor => patch%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (sinks_t)
    call self%SN (patch, task)
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Compute the SN input to the patch
!===============================================================================
SUBROUTINE SN (self, patch, sink)
  class(SN_t):: self
  class(patch_t):: patch
  class(sinks_t):: sink
  class(task_t), pointer:: task
  real(8):: m_sun               ! solar mass in code units
  real(4):: dm                  ! change in mass per cell
  real(4):: dd                  ! change in mass density
  real(4):: mdot                ! mass loss rate
  real, pointer, dimension(:,:,:):: d, px, py, pz, e
  real, allocatable, dimension(:,:,:):: x, y, z, r
  logical, allocatable:: inside(:,:,:)
  integer:: i, i1, i2, i3, m(6), n_inside
  real:: ee, ee_min, foe
  logical:: first_time = .true.
  !-----------------------------------------------------------------------------
  call trace%begin ('SN_t%SN')
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
  inside = (r < sink%radius)
  dm = sink%mass
  foe = 1d51/scaling%e
  ee = foe/dm
  self%velocity = sqrt(2.0*ee)
  n_inside = count(inside)
  !$omp critical (input_cr)
  if (n_inside > 0 .and. first_time) then
    print *, 'BOOOOOOOOOOOOOOOOOOOOOOOOOM'
    dm = dm/n_inside
    dd = dm/product(patch%ds)
    !---------------------------------------------------------------------------
    ! Add mass and momentum inside the given radius, letting the velocity approach
    ! a pattern with constant divergence
    !---------------------------------------------------------------------------
    d  => patch%mem(:,:,:,patch%idx%d ,patch%it,1)
    e  => patch%mem(:,:,:,patch%idx%e ,patch%it,1)
    px => patch%mem(:,:,:,patch%idx%px,patch%it,1)
    py => patch%mem(:,:,:,patch%idx%py,patch%it,1)
    pz => patch%mem(:,:,:,patch%idx%pz,patch%it,1)
    i = patch%n(1)/2
    print *, d(i,i,i), dd
    print *, e(i,i,i)/d(i,i,i), ee
    where (inside)
      e  = e - 0.5*(px**2+py**2+pz**2)/d
      d  = d  + dd
      px = px + dd*(self%velocity*x/sink%radius + sink%velocity(1))
      py = py + dd*(self%velocity*y/sink%radius + sink%velocity(2))
      pz = pz + dd*(self%velocity*z/sink%radius + sink%velocity(3))
      e  = e + 0.5*(px**2+py**2+pz**2)/d
    end where
    if (verbose > 1) &
      write (io_unit%output,*) &
        ' after:', patch%id, patch%fmaxval(abs(px/d)), maxval(abs(px/d))
    first_time = .false.
    write (io_unit%output,'(1x,a,i5,i8,i4,1p,2e11.2)') &
      'SN_t: id, n, i, e, T =', patch%id, n_inside, i, e(i,i,i), &
      e(i,i,i)/d(i,i,i)*(patch%gamma-1.0)*scaling%temp
  end if
  deallocate (inside, x, y, z, r)
  !$omp end critical (input_cr)
  call trace%end ()
END SUBROUTINE SN

END MODULE SN_mod
