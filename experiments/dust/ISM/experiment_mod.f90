!===============================================================================
!> A simple experiment setup, which defines the initial sate of each patch.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE task_mod
  USE solver_mod
  USE scaling_mod
  USE units_mod
  USE particle_solver_mod
  USE link_mod
  USE list_mod
  USE scratch_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
    procedure:: init_task_list
    procedure:: update
  end type
  real    :: T_K=10
  real(8) :: gamma=1d0
  integer :: verbose=1
CONTAINS

!===============================================================================
!> Setup an initial state.  Values are given in code units, defined in the
!> scaling_mod.f90 file and by the scaling_params namelist.
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t) :: self
  integer             :: iostat
  real, save          :: d0, e0, csound
  logical, save       :: first_time=.true.
  real, pointer, dimension(:,:,:):: d, e, px, py, pz
  real, pointer, dimension(:,:,:,:):: p
  namelist /experiment_params/ verbose, T_K, gamma
  character(len=120):: id= &
  '$Id$'
  !-----------------------------------------------------------------------------
  ! Initialize the task instance, which includes allocating self%mem()
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call io_unit%print_id (id)
  call self%solver_t%init
  !-----------------------------------------------------------------------------
  ! Read experiment parameters and compute initial condition values
  !-----------------------------------------------------------------------------
  call io%lock%set
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%namelist, experiment_params)
    csound = sqrt(cgs%k_B*T_K/(cgs%m_u*scaling%mu))/scaling%u
    d0 = 1.0
    if (gamma == 1d0) then
      e0 = self%csound**2*d0
    else
      e0 = self%csound**2*d0/(self%gamma*(self%gamma-1.0))
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
  self%gamma = gamma
  self%csound = csound
  self%mem(:,:,:,self%idx%d,:,1) = d0
  self%mem(:,:,:,self%idx%e,:,1) = e0
  d => self%mem(:,:,:,self%idx%d,self%it,1)
  e => self%mem(:,:,:,self%idx%e,self%it,1)
  p => self%mem(:,:,:,self%idx%px:self%idx%pz,self%it,1)
  !-----------------------------------------------------------------------------
  ! Use pre_update to get access to the force (ifort 20.1 problem circumvent)
  !-----------------------------------------------------------------------------
  call self%forces%pre_update
  p = p + d0*scratch%force_per_unit_mass*self%forces%force%t_turb
  if (self%kind(1:6) == 'ramses') then
    px => self%mem(:,:,:,self%idx%px,self%it,1)
    py => self%mem(:,:,:,self%idx%py,self%it,1)
    pz => self%mem(:,:,:,self%idx%pz,self%it,1)
    e = e + (px**2+py**2+pz**2)/d
  end if
  !call self%initial%condition (self%mesh, self%mem(:,:,:,:,1,1), self%idx)
  !-----------------------------------------------------------------------------
  ! Initialize particles
  !-----------------------------------------------------------------------------
  self%periodic(:) = .true.
  call self%particle_solver%initial_condition
  call trace%end ()
END SUBROUTINE init

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
  call trace%end ()
END SUBROUTINE update

END MODULE
