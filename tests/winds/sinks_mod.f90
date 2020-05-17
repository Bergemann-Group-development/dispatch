!===============================================================================
!> Sink module stub -- needs only to provide mass and position
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE patch_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  type, public, extends(task_t):: sinks_t
    real:: mass, radius
    real(8), allocatable:: pos(:,:)
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
  allocate (self%pos(3,self%nt))
  call self%task_t%init
  self%size = 2d0*self%radius
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE update (self)
  class(sinks_t):: self
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%update')
  self%dtime = 0.01
  self%position = self%position + self%dtime*self%velocity
  if (verbose > 0) &
    print *, self%id, 'sinks_t: moved to', real(self%position,kind=4)
  call trace%end ()
END SUBROUTINE update

END MODULE sinks_mod
