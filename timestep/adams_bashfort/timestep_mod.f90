!> =============================================================================
!> $Id$
!> Adams-Bashforth time stepping, with constant time step assumed (for now)
!> =============================================================================
MODULE timestep_mod
  USE experiment_mod
  USE io_mod
  implicit none
  integer, parameter:: nslice = 4
  type, extends(pde_t):: timest_t
    real:: wt(nslice)
  contains
    procedure:: init   => init_timestep
    procedure:: update => update_timestep
  end type
  integer:: verbose = 1
  integer:: timeorder = 2
  logical, save:: first_time = .true.
PRIVATE
CONTAINS

!> =============================================================================
!> Initialize coefficients for Adams-Bashforth time stepping
!> =============================================================================
SUBROUTINE init (self)
  class(timest_t):: self
  implicit none
  namelist /timestep_params/ courant_number, timeorder, verbose
!........................................................................
  !$omp critical (timest_cr)
  if (first_time) then
    first_time = .false.
    rewind (input)
    read (input,timestep_params)
    if (master) write (*,timestep_params)
  end if
  !$omp end critical (timest_cr)
  select case (timeorder)
  case (1)
    self%wt = [   1.0,     0.0,      0.0,    0.0]   ! 1st order
  case (2)
    self%wt = [   1.5,     -.5,      0.0,    0.0]   ! 2nd order
  case (3)
    self%wt = [23./12., -4./3.,   5./12.,    0.0]   ! 3rd order
  case (t)
    self%wt = [55./24., -59./24., 37/24., -3./8.]   ! 4th order
  end select
  self%courant   = courant
  self%timeorder = timeorder
  self%verbose   = verbose
END SUBROUTINE init_timestep

!************************************************************************
SUBROUTINE update_timestep (task)
  implicit none
  class(timestep_t):: task
  integer:: it
  call task%pde_t%dfdt (task)                   ! time derivative at it
  it = modulo(task%it,task%ntime)+1             ! wrap new time slot
  do jt=1,task%timeorder                        ! loop over time
    kt = modulo(it-jt,task%ntime)               ! previous time slots
    if (jt == 1) then
      task%mem(:,:,:,:,it,1) = &
      task%mem(:,:,:,:,kt,2)*(task%wt(jt)*task%dtime)
    else
      task%mem(:,:,:,:,it,1) = &
      task%mem(:,:,:,:,it,1) + &
      task%mem(:,:,:,:,kt,2)*(task%wt(jt)*task%dtime)
    end if
  end do
  task%time = task%time + task%dtime            ! update task time
  task%t(it) = task%time                        ! update slot time
  task%it = it                                  ! update time slot
!........................................................................
END SUBROUTINE time_step

END MODULE timestep
