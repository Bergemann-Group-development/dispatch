!===============================================================================
!> A simple experiment setup, which defines the initial sate of each patch.
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE trace_mod
  USE solver_mod
  USE scaling_mod
  USE units_mod
  USE link_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure:: init
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
  real, save          :: d, e, csound
  logical, save       :: first_time=.true.
  namelist /experiment_params/ verbose, T_K, gamma
  !-----------------------------------------------------------------------------
  ! Initialize the task instance, which includes allocating self%mem()
  !-----------------------------------------------------------------------------
  call trace%begin ('experiment_t%init')
  call self%solver_t%init
  !-----------------------------------------------------------------------------
  ! Read experiment parameters and compute initial condition values
  !-----------------------------------------------------------------------------
  !$omp critical (input_cr)
  if (first_time) then
    rewind (io%input)
    read (io%input, experiment_params, iostat=iostat)
    call io%namelist_warning ('experiment_params', iostat)
    write (io%output, experiment_params)
    csound = sqrt(cgs%k_B*T_K/(cgs%m_u*scaling%mu))/scaling%u
    d = 1.0
    if (gamma == 1d0) then
      e = self%csound**2*d
    else
      e = self%csound**2*d/(self%gamma*(self%gamma-1.0))
    end if
    if (verbose > 0) then
      write (io%output,'(a,f6.2)') '           gas gamma:', gamma
      write (io%output,'(a,f6.2)') 'density (code units):', d
      write (io%output,'(a,f6.2)') ' energy (code units):', e
      write (io%output,'(a,f6.2)') '     temperature (K):', T_K
      write (io%output,'(a,f6.2)') '  sound speed (km/s):', &
        self%csound*scaling%u/cgs%kms
    end if
    first_time = .false.
  end if
  !$omp end critical (input_cr)
  !-----------------------------------------------------------------------------
  ! Store the saved values in each of the many task instances
  !-----------------------------------------------------------------------------
  self%gamma = gamma
  self%csound = csound
  self%mem(:,:,:,self%idx%d,1,1) = d
  self%mem(:,:,:,self%idx%e,1,1) = e
  self%periodic(:) = .true.
  call trace%end ()
END SUBROUTINE init

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
