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
!>
!> Each patch that makes use of the force() procedure to compute the force on
!> its gas cells also removes the corresponding amount of momentum from the sink.
!> The %force() prodedure must therefor have an argument that gives the mass of
!> the object making the call.
!>
!> There is no problem computing the direct interaction between the sink and
!> all other sink particles, as long as it doesn't have to be done very often.
!> However, for gas to particle interaction, the influence of distant, low
!> mass sinks should be combined into a single mass point for each patch they
!> belong to.
!>
!> Stars should be grouped into cluster tasks, which have a total mass and 
!> a center of mass.  These should again be collected into a population task,
!> which has a list of all cluster tasks.
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
    real(8):: mass=0.0, radius=8.0
    real(8), allocatable:: pos(:,:)
    real(8):: SN_time=0.0, previous_time=0d0, dt_min=1e-7
  contains
    procedure:: init
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
SUBROUTINE init (self)
  class(sinks_t):: self
  integer:: iostat
  real:: velocity(3)=0.0, position(3)=0.5
  logical, save:: first_time=.true.
  namelist /sinks_params/ verbose, mass_msun, radius, velocity, SN_time, &
    position
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%init')
  if (first_time) then
    !$omp critical (sinks_cr)
    if (first_time) then
      !$omp critical (input_cr)
      rewind (io%input)
      read (io%input, sinks_params, iostat=iostat)
      write (io%output, sinks_params)
      !$omp end critical (input_cr)
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
END SUBROUTINE init

!===============================================================================
!> Update sink task
!===============================================================================
SUBROUTINE pre_update (self, patch)
  class(sinks_t):: self
  class(patch_t):: patch
  logical:: contains_sink
  !-----------------------------------------------------------------------------
  call trace%begin ('sinks_t%pre_update')
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
