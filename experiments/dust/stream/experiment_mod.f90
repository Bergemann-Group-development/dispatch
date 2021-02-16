!===============================================================================
!> A simple experiment setup, which defines the initial sate of each patch.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE task_mod
  USE solver_mod
  USE scaling_mod
  USE units_mod
  USE particle_solver_mod
  USE particles_mod
  USE link_mod
  USE list_mod
  USE scratch_mod
  USE eos_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
    real(8):: momentum(3,2)=0d0
  contains
    procedure:: init
    procedure:: init_bdries
    procedure:: init_task_list
    procedure:: init_stream
    procedure:: update
  end type
  real    :: T_K=10
  real(8) :: gamma=1d0
  real    :: x0=0.0, y0=0.0, v0=1.0, w0=0.2
  integer :: verbose=1
  integer :: case=1
CONTAINS

!===============================================================================
!> Setup an initial state.  Values are given in code units, defined in the
!> scaling_mod.f90 file and by the scaling_params namelist.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t) :: self
  integer             :: i, j, iostat
  real, save          :: d0, e0, csound
  logical, save       :: first_time=.true.
  real, pointer, dimension(:,:,:):: d, e, px, py, pz
  real, pointer, dimension(:,:,:,:):: p
  type(mpi_counter_t), pointer:: counter
  namelist /experiment_params/ verbose, T_K, gamma, x0, y0, v0, w0, case
  !-----------------------------------------------------------------------------
  ! Initialize the task instance, which includes allocating self%mem()
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call self%solver_t%init
  !-----------------------------------------------------------------------------
  ! Read experiment parameters and compute initial condition values
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    gamma = eos%gamma
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    eos%gamma = gamma
    csound = sqrt(cgs%k_B*T_K/(cgs%m_u*scaling%mu))/scaling%u
    d0 = 1.0
    if (gamma == 1d0) then
      e0 = self%csound**2*d0
    else
      e0 = self%csound**2*d0/(eos%gamma*(eos%gamma-1.0))
    end if
    if (verbose > 0) then
      write (io%output,'(a,f6.2)') '           gas gamma:', gamma
      write (io%output,'(a,f6.2)') 'density (code units):', d0
      write (io%output,'(a,f6.2)') ' energy (code units):', e0
      write (io%output,'(a,f6.2)') '     temperature (K):', T_K
      write (io%output,'(a,f6.2)') '  sound speed (km/s):', &
        self%csound*scaling%u/cgs%kms
    end if
    first_time = .false.
  end if
  call io%lock%unset
  !-----------------------------------------------------------------------------
  ! Initialize constant density and thermal energy, and a velocity estimated
  ! from constant acceleration for a time period of t_turb.
  !-----------------------------------------------------------------------------
  eos%gamma = gamma
  self%csound = csound
  self%mem(:,:,:,self%idx%d,:,1) = d0
  self%mem(:,:,:,self%idx%e,:,1) = e0
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  e => self%mem(:,:,:,self%idx%e,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  !-----------------------------------------------------------------------------
  ! Use pre_update to get access to the force (ifort 20.1 problem circumvent)
  !-----------------------------------------------------------------------------
  !call self%pre_update
  !p = p + d0*scratch%force_per_unit_mass*self%forces%force%t_turb
  if (self%kind(1:6) == 'ramses') then
    px => self%mem(:,:,:,self%idx%px,self%it,1)
    py => self%mem(:,:,:,self%idx%py,self%it,1)
    pz => self%mem(:,:,:,self%idx%pz,self%it,1)
    e = e + (px**2+py**2+pz**2)/d
  end if
  !-----------------------------------------------------------------------------
  ! Initialize particles
  !-----------------------------------------------------------------------------
  self%periodic(:) = .true.
  call self%particle_solver%initial_condition
  call self%init_stream
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Stub / intercept -- no boundary processing neeeded
!===============================================================================
SUBROUTINE init_bdries (self)
  class(experiment_t) :: self
  !-----------------------------------------------------------------------------
  !call self%solver_t%init_bdries
END SUBROUTINE init_bdries

!===============================================================================
!> Create an initial condition with a stream of particles in the z-direction
!===============================================================================
SUBROUTINE init_stream (self)
  class(experiment_t):: self
  integer:: i, i1, i2, i3, i4, j, lb(3), ub(3)
  real:: r(3), dp(3)
  real(8):: llc(3), momentum(3)
  type(particle_t), pointer:: p
  type(particles_t), pointer:: particles
  type(mpi_counter_t), pointer:: counter
  !-----------------------------------------------------------------------------
  call trace%begin ('init_stream')
  llc = self%position - self%size/2. - self%mesh%ng*self%mesh%d
  lb = self%mesh%li
  ub = self%mesh%ui
  momentum = 0d0
  do i3=lb(3),ub(3)
  do i2=lb(2),ub(2)
  do i1=lb(1),ub(1)
    r = ([i1,i2,i3] - 1.0)*self%mesh%d
    r = r + llc
    dp = 0.0
    do i4=1,self%particle_solver%n_bin
      if (case==2) then
        self%particle_solver%dd(i1,i2,i3,i4) = &
        self%particle_solver%dd(i1,i2,i3,i4) &
          *exp(-((r(1)-x0)**2+(r(2)-y0)**2)/w0**2)
      end if
      particles => self%particle_solver%particles(i1,i2,i3,i4)
      particles%momentum = 0d0
      p => particles%head%next
      do j=1,particles%n
        r = (p%r + p%q(1:3) - 1.0)*self%mesh%d
        r = r + llc
        p%v(1) = 0.0
        p%v(2) = 0.0
        !-----------------------------------------------------------------------
        ! Case 2 har varying dust density and constant dust velocity
        !-----------------------------------------------------------------------
        if (case==2) then
          p%v(3) = v0/self%ds(3)
          p%w = self%particle_solver%dd(i1,i2,i3,i4)/particles%n
        !-----------------------------------------------------------------------
        ! Case 1 har constant dust density and varying dust velocity
        !-----------------------------------------------------------------------
        else
          p%v(3) = v0/self%ds(3)*exp(-((r(1)-x0)**2+(r(2)-y0)**2)/w0**2)
        end if
        dp = dp + p%v*p%w
        p => p%next
      end do
    end do
    !---------------------------------------------------------------------------
    ! Accumulate particle momentum, converting to momentum per unit volume
    !---------------------------------------------------------------------------
    momentum = momentum + dp*self%ds
    !...........................................................................
    if (verbose > 0) then
      r = ([i1,i2,i3] - 1.0)*self%mesh%d
      r = r + self%position - self%size/2. - self%mesh%ng*self%mesh%d
      if (r(1)**2 + r(2)**2 < self%ds(1)**2) then
        write (stdlog,'(a,i6,3i4,2f8.3,5(2x,1p,3e12.3))') &
        'stream',self%id,i1,i2,i3,r(1:2),particles%head%next%v, dp, &
        self%particle_solver%global_momentum ()
      end if
    end if
  end do
  end do
  end do
  self%particle_solver%vmax = maxval(v0/self%ds)
  !write (stdlog,*) 'dtime(0)', self%particle_solver%dtime
  !-----------------------------------------------------------------------------
  ! Update the global particle and gas momentum counters
  !-----------------------------------------------------------------------------
  self%particle_solver%old_momentum = 0d0
  self%particle_solver%new_momentum(:,1) = momentum*product(self%ds)
  self%particle_solver%new_momentum(:,2) = 0d0
  call self%particle_solver%update_momentum
  call trace%end ()
END SUBROUTINE init_stream

!===============================================================================
!> This procedure is the place to add calls to extras 'minions' that need
!> information about the corresponding MHD patch (self), AND the nbor or task
!> list.  The procedure is called for every task on the task list, just before
!> task list updates start.
!===============================================================================
SUBROUTINE init_task_list (self, task_list)
  class(experiment_t):: self
  class(list_t), pointer:: task_list
  !.............................................................................
  class(nbor_t), pointer:: nbor, nbor2
  class(task_t), pointer:: task
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init_task_list')
  call self%particle_solver%init_task_list (task_list)                  ! part
  call self%gpatch_t%init_task_list (task_list)                         ! keep
  nbor => self%link%nbor
  do while (associated(nbor))
    task => nbor%task
    select type (task)
    class is (experiment_t)
      !-------------------------------------------------------------------------
      ! Allocate an nbor, point to the particle_solver_t part of it, and prepend
      !-------------------------------------------------------------------------
      allocate (nbor2)
      nbor2%task => task%particle_solver
      nbor2%next => self%particle_solver%nbor
      self%particle_solver%nbor => nbor2
    end select
    nbor => nbor%next
  end do
  call trace%end ()
END SUBROUTINE init_task_list

!===============================================================================
!> Update the task state, by calling the solver (this is redundant -- leaving
!> it out would automatically pass the call on to solver_t, which this data type
!> is an extension of.
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace%begin ('experiment_t%update')
  call self%solver_t%update
  if (verbose > 0) then
    write (stdlog,'(a,i6,3(2x,1p,3e12.3))') &
    'experiment_t%update: ', self%id, self%particle_solver%global_momentum ()
  end if
  call trace%end ()
END SUBROUTINE update

END MODULE
