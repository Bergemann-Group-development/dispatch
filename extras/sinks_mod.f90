!===============================================================================
!> Sink module stub -- needs only to provide mass and position
!===============================================================================
MODULE sinks_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE gpatch_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  type, public, extends(task_t):: sinks_t
    real(8):: mass=0.0, radius=8.0
    real(8), allocatable:: pos(:,:)
    real(8):: SN_time=0.0, previous_time=0d0, dt_min=1e-7
    class(gpatch_t), pointer:: patch => null()
  contains
    procedure:: pre_init
    procedure:: connect_patch
    procedure:: pre_update
  end type
  integer:: verbose=0
  real:: radius=8.0
  real:: mass_msun=30
  real:: SN_time=0.1
CONTAINS

!===============================================================================
!> Initialize sink task
!===============================================================================
SUBROUTINE pre_init (self)
  class(sinks_t), target:: self
  class(gpatch_t), pointer:: patch
  !.............................................................................
  integer:: iostat
  real:: velocity(3)=0.0, position(3)=0.5
  logical, save:: first_time=.true.
  namelist /sinks_params/ verbose, mass_msun, radius, velocity, SN_time, &
    position
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%pre_init')
  patch => self%patch
  if (first_time) then
    !$omp critical (sinks_cr)
    if (first_time) then
      call io%lock%set
      rewind (io%input)
      read (io%input, sinks_params, iostat=iostat)
      call io%namelist_warning ('sinks_params', iostat)
      write (io%namelist, sinks_params)
      call io%lock%unset
      self%radius = radius
      self%velocity = velocity
      self%SN_time = SN_time
      if (verbose > 0) then
        write (stdout,*) 'sinks_t%init: radius (code units)', self%radius
      end if
      self%nt = 5
      allocate (self%pos(3,self%nt))
      self%position = position
      self%it = 1
      self%size = 0d0
      self%mass = mass_msun*cgs%m_sun/scaling%m
      call self%task_t%init
      self%iit = 1
      self%pos(:,1) = position
      first_time = .false.
    end if
    !$omp end critical (sinks_cr)
  end if
  call trace%end ()
END SUBROUTINE pre_init

!===============================================================================
!> Store a pointer to the patch the sink lives in
!===============================================================================
SUBROUTINE connect_patch (self, patch)
  class(sinks_t), target:: self
  class(gpatch_t), pointer:: patch
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%connect_patch')
  self%patch => patch
  call trace%end ()
END SUBROUTINE connect_patch

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE pre_update (self)
  class(sinks_t):: self
  !.............................................................................
  class(gpatch_t), pointer:: patch
  logical:: contains_sink
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%pre_update')
  call io%assert (associated(self%patch), 'sink_t%patch not associated')
  patch => self%patch
  contains_sink = patch%contains (self%position)
  if (verbose > 1) then
    print '(a,i6,2(2x,3f9.6),l5)', 'sinks_t%pre_update: patch%id, pos, patch%pos =', &
      patch%id, self%position, (patch%position-self%position)/patch%size, contains_sink
  end if
  if (contains_sink) then
    self%dtime = minval(patch%ds)/sqrt(sum(self%velocity**2))
    if (self%time < patch%time) then
      self%previous_time = self%time
      do while (self%time < patch%time)
        !-----------------------------------------------------------------------
        ! Make sure to hit the exact patch_time + a well defined extra dt_min
        !-----------------------------------------------------------------------
        self%dtime = min(self%dtime, patch%time-self%time+self%dt_min)
        self%position = self%position + self%dtime*self%velocity
        self%pos(:,self%new) = self%position
        if (verbose > 0) then
          write (stdout,*) &
            self%id, patch%id, self%time, self%it, 'sinks_t%update: moved to', real(self%position,kind=4)
        end if
        self%rotated = .false.
        call self%rotate
      end do
    end if
  end if
  call trace%end ()
END SUBROUTINE pre_update

END MODULE sinks_mod
