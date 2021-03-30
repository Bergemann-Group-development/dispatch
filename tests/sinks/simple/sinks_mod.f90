!===============================================================================
!> Sink module stub -- needs only to provide mass and position
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE link_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  type, public, extends(task_t):: sinks_t
    type(link_t), pointer:: link
    real:: radius
    real(8):: mass
    real(8):: mom(3)
    real(8), allocatable:: pos(:,:)
  contains
    procedure:: init
    procedure:: update
    procedure:: move
    procedure:: accrete
  end type
  integer:: verbose=0
  real:: radius_au=4e3
  real:: accretion_rate
CONTAINS

!===============================================================================
!> Initialize sink task
!===============================================================================
SUBROUTINE init (self)
  class(sinks_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /sinks_params/ verbose, radius_au  
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%init')
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
  self%mass = 0d0
  self%radius = radius_au*cgs%au/scaling%l
  if (verbose > 0) &
    print *,'radius',self%radius
  self%nt = 5
  allocate (self%pos(3,self%nt))
  call self%task_t%init
  self%size = 2d0*self%radius
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Find the nbor patch tasks and accrete from them
!===============================================================================
SUBROUTINE update (self)
  class(sinks_t):: self
  class(link_t), pointer:: nbor, this, oldest
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%update')
  !-----------------------------------------------------------------------------
  ! Loop over nbors, sorted over time
  !-----------------------------------------------------------------------------
  this => self%link%nbor
  do while (associated(this))
    !---------------------------------------------------------------------------
    ! Find the oldest task, starting the search with this
    !---------------------------------------------------------------------------
    oldest => this
    nbor => oldest
    do while (associated(nbor))
      !------------------------------------------------------------------------
      ! If an older one is found, change the pointer
      !------------------------------------------------------------------------
      if (nbor%task%time < oldest%task%time) then
        oldest => nbor
      end if
    end do
    task => oldest%task
    select type (task)
    class is (patch_t)
      !----------------------------------------------------------------------
      ! Move to the mean patch time
      !----------------------------------------------------------------------
      self%dtime = task%time - 0.5_8*task%dtime - self%time
      call self%move
      !----------------------------------------------------------------------
      ! Accrete over the last patch time interval
      !----------------------------------------------------------------------
      call self%accrete (task, task%dtime)
    end select
    this => this%next
  end do
  call trace%end ()
END SUBROUTINE update

!===============================================================================
!> Move the sink from self%time to self%time + self%dtime.  This skeleton
!> version uses simpe 1st oder time stepping, as a demo for now
!===============================================================================
SUBROUTINE move (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  self%position = self%position + self%mom/self%mass * self%dtime
  self%pos(:,self%it) = self%position
  call self%rotate
END SUBROUTINE move

!===============================================================================
!> Accrete from the patch, at the exact patch%time
!===============================================================================
SUBROUTINE accrete (self, patch, dtime)
  class(sinks_t):: self
  class(patch_t):: patch
  real(8):: dtime
  !.............................................................................
  class(task_t), pointer:: task
  real(8):: m_sun               ! solar mass in code units
  real(4):: dm                  ! change in mass per cell
  real(4):: ramp                ! ramp profile
  real(4):: mdot                ! mass loss rate
  real, pointer, dimension(:,:,:):: d, px, py, pz, e
  real:: x, y, z, r, f, dd, dpx, dpy, dpz, profile, dv, pt(2)
  integer:: i, i1, i2, i3, m(6), n_inside, jt(2)
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%wind')
  !-----------------------------------------------------------------------------
  ! Distances from sink center
  !-----------------------------------------------------------------------------
  m = shape(patch%mem)
  !-----------------------------------------------------------------------------
  ! Compute the change in mass in the timestep, and distribute over the cells
  !-----------------------------------------------------------------------------
  if (n_inside > 0) then
    !---------------------------------------------------------------------------
    ! Pick up mass and momentum within the given radius
    !---------------------------------------------------------------------------
    d  => patch%mem(:,:,:,patch%idx%d ,patch%it,1)
    px => patch%mem(:,:,:,patch%idx%px,patch%it,1)
    py => patch%mem(:,:,:,patch%idx%py,patch%it,1)
    pz => patch%mem(:,:,:,patch%idx%pz,patch%it,1)
    dv = product(patch%ds)
    do i3=1,m(3)
    do i2=1,m(2)
    do i1=1,m(1)
      x = (patch%mesh(1)%p - self%position(1)) + patch%mesh(1)%r(i1)
      y = (patch%mesh(2)%p - self%position(2)) + patch%mesh(2)%r(i2)
      z = (patch%mesh(3)%p - self%position(3)) + patch%mesh(3)%r(i3)
      r = sqrt(x**2+y**2+z**2)/self%radius
      if (r <= 1.0) then
        profile = (1.0-3.0*r**2+2.0*r**3)
        f=self%dtime*accretion_rate*profile
        dd  = f*d (i1,i2,i3)
        dpx = f*px(i1,i2,i3)
        dpy = f*py(i1,i2,i3)
        dpz = f*pz(i1,i2,i3)
        d (i1,i2,i3)  = d (i1,i2,i3) - dm
        px(i1,i2,i3)  = px(i1,i2,i3) - dpx
        py(i1,i2,i3)  = py(i1,i2,i3) - dpy
        pz(i1,i2,i3)  = pz(i1,i2,i3) - dpz
        self%mass = self%mass + dv*dd
        self%mom(1) = self%mom(1) + dv*dpx
        self%mom(2) = self%mom(2) + dv*dpy
        self%mom(3) = self%mom(3) + dv*dpz
      end if
    end do
    end do
    end do
  end if
  call trace%end ()
END SUBROUTINE accrete

END MODULE sinks_mod
