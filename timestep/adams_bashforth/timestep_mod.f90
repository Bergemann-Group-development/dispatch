!> =============================================================================
!> $Id$
!> Adams-Bashforth time stepping, with constant time step assumed (for now)
!>
!> We store results in nt memory slots, starting with number 1, and continuing
!> with numbers 2..nt, and then looping back to 1, but this is mapped by the
!> iit array so that the last nt-1 time slot are t(iit(1:nt-1)), and t(iit(nt))
!> is the next time slot.  The (1:nt-1) slots are complete, and may be used
!> freely, while the (nt) slot is the one where new values are to be stored.
!> The experiment_t%update procedure thus writes new time derivatives
!> there, and it is the task of this timestep procedure to update the values
!> there, using them time derivatives and the previous values.
!> =============================================================================
MODULE timestep_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  integer, parameter:: nslice = 4
  type, public:: timestep_t
    integer:: id                                        ! task id
    integer:: verbose=0                                 ! verbosity
    integer:: time_order=-1                             ! time order of integration
    integer, dimension(:), pointer:: iit                ! time slot order
    real:: wt(nslice)                                   ! time slot weights
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!> =============================================================================
!> Initialize coefficients for Adams-Bashforth time stepping
!> =============================================================================
SUBROUTINE init (self)
  class(timestep_t):: self
  integer, save:: verbose=0, time_order=-1
  namelist /timestep_params/ verbose, time_order
  !.............................................................................
  call trace_begin ('timestep_t%init', 1)
  !$omp critical (timestep_cr)
  if (time_order==-1) then
    time_order = 2
    rewind (io%input)
    read (io%input,timestep_params)
    if (io%master) write (*,timestep_params)
  end if
  self%verbose = verbose
  self%time_order = time_order
  !---------------------------------------------------------------------------
  select case (time_order)
  case (0)
  case (1)
    self%wt = [   1.0,     0.0,      0.0,    0.0]       ! 1st order
  case (2)
    self%wt = [   1.5,     -.5,      0.0,    0.0]       ! 2nd order
  case (3)
    self%wt = [23./12., -4./3.,   5./12.,    0.0]       ! 3rd order
  case (4)
    self%wt = [55./24., -59./24., 37/24., -3./8.]       ! 4th order
  case default
    if (io%master) print*,'ERROR: illegal timestep time_order =', time_order
    stop
  end select
  !$omp end critical (timestep_cr)
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Adams-Bashforth linear multi-step update.  We consider it to be linear
!> multi-step in time index space, which may be just as good -- or better --
!> than attempting to map into linear time space.  The time slot derivative
!> values have already been multiplied with the individual time steps dt.
!> =============================================================================
SUBROUTINE update (self, id, iit, dt, mem, mesh)
  class(timestep_t):: self
  class(mesh_t):: mesh(:)
  integer:: id, nt
  integer, dimension(:):: iit
  real(8), dimension(:):: dt
  real, dimension(:,:,:,:,:,:):: mem            ! memory buffer
  integer:: l(3), u(3)
  !......................................................................
  call trace_begin ('timestep_t%update')
  call self%init                                ! auto-initialize
  nt = size(iit)                                ! number of time slots
  l = mesh%li
  u = mesh%ui
  where (mesh%lower_boundary)
    l = mesh%lb
  end where
  where (mesh%upper_boundary)
    u = mesh%ub
  end where
  select case (self%time_order)
  case (0)
  case (1)
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt  ),1) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),1) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),2)*(self%wt(1)*dt(iit(nt-1)))
  case (2)
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt  ),1) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),1) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),2)*(self%wt(1)*dt(iit(nt-1))) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-2),2)*(self%wt(2)*dt(iit(nt-2)))
  case (3)
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt  ),1) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),1) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),2)*(self%wt(1)*dt(iit(nt-1))) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-2),2)*(self%wt(2)*dt(iit(nt-2))) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-3),2)*(self%wt(3)*dt(iit(nt-3)))
  case (4)
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt  ),1) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),1) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-1),2)*(self%wt(1)*dt(iit(nt-1))) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-2),2)*(self%wt(2)*dt(iit(nt-2))) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-3),2)*(self%wt(3)*dt(iit(nt-3))) + &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),:,iit(nt-4),2)*(self%wt(4)*dt(iit(nt-4)))
  end select
  if (self%verbose>3) then
    print'(a,i4,2x,a,1p,6e10.2,2x,6i5)', 'timestep, id:', &
      id, 'minmax:', &
      minval(mem(l(1):u(1),l(2):u(2),l(3):u(3),1,iit(nt-1),2)), &
      maxval(mem(l(1):u(1),l(2):u(2),l(3):u(3),1,iit(nt-1),2)), &
      minval(mem(l(1):u(1),l(2):u(2),l(3):u(3),1,iit(nt-1),1)), &
      maxval(mem(l(1):u(1),l(2):u(2),l(3):u(3),1,iit(nt-1),1)), &
      minval(mem(l(1):u(1),l(2):u(2),l(3):u(3),1,iit(nt  ),1)), &
      maxval(mem(l(1):u(1),l(2):u(2),l(3):u(3),1,iit(nt  ),1)), shape(mem)
  end if
  call trace_end
END SUBROUTINE update
!-------------------------------------------------------------------------------
END MODULE timestep_mod
