!===============================================================================
!> Sink module stub -- needs only to provide mass and position
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE scaling_mod
  USE link_mod
  USE units_mod
  implicit none
  private
  type, public, extends(task_t):: sinks_t
    real:: mass, radius
    class(link_t), pointer:: link => null()
  contains
    procedure:: init
    procedure:: update
  end type
  integer:: verbose=0
  real:: radius=8.0
  real:: mass_msun=30
CONTAINS

!===============================================================================
!> Initialize sink task
!===============================================================================
SUBROUTINE init (self)
  class(sinks_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  namelist /sinks_params/ verbose, mass_msun, radius
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
  self%mass = mass_msun*cgs%m_sun/scaling%m
  self%radius = radius
  if (verbose > 0) &
    print *,'radius',self%radius
  self%nt = 5
  call self%task_t%init
  self%size = 2d0*self%radius
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE update (self)
  class(sinks_t):: self
  class(nbor_t), pointer:: nbor
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%update')
  self%dtime = huge(1d0)
  nbor => self%link%nbor
  do while (associated(nbor))
    self%dtime = min(self%dtime,nbor%task%dtime)
    nbor => nbor%next
  end do
  self%position = self%position + self%dtime*self%velocity
  if (verbose > 0) &
    print '(i6,2x,a,1p,2e12.4,2x,3e12.4)', &
      self%id, 'sinks_t: dt, t, moved to', self%dtime, self%time, self%position
  call trace%end ()
END SUBROUTINE update

END MODULE sinks_mod
