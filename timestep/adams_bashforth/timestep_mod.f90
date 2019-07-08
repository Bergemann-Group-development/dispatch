!===============================================================================
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
!===============================================================================
MODULE timestep_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  USE omp_lock_mod
  implicit none
  private
  integer, parameter:: nslice = 4
  type, public:: timestep_t
    integer:: verbose=0                                 ! verbosity
    integer:: time_order=-1                             ! time order of integration
    real:: wt(nslice)                                   ! time slot weights
  contains
    procedure:: init
    procedure:: update
    procedure:: update_single
  end type
  type(timestep_t), save, public:: timestep
CONTAINS

!===============================================================================
!> Initialize coefficients for Adams-Bashforth time stepping
!===============================================================================
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
  !$omp end critical (timestep_cr)
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
  call trace_end
END SUBROUTINE init

!===============================================================================
!> Adams-Bashforth linear multi-step update.  We consider it to be linear
!> multi-step in time index space, which may be just as good -- or better --
!> than attempting to map into linear time space.  The time slot derivative
!> values have already been multiplied with the individual time steps dt.
!===============================================================================
SUBROUTINE update (self, id, iit, dt, mem, mesh, lock)
  class(timestep_t):: self
  class(mesh_t):: mesh(:)
  type(lock_t):: lock
  integer:: id, nt
  integer, dimension(:):: iit
  real(8), dimension(:):: dt
  real(kind=KindScalarVar), dimension(:,:,:,:,:,:):: mem            ! memory buffer
  integer:: l(3), u(3)
  integer, save:: itimer=0
  !......................................................................
  call trace_begin ('timestep_t%update', itimer=itimer)
  if (self%time_order==-1) call io%abort('timestep_mod not initialized!')
  nt = size(iit)                                ! number of time slots
  l = mesh%lb
  u = mesh%ub
  where (mesh%lower_boundary)
    l = mesh%lb
  end where
  where (mesh%upper_boundary)
    u = mesh%ub
  end where
  if (size(mem,6)==1.and.self%time_order>0) then
    if (io%master) print *,id,'WARNING: setting time_order=0 in timestep%update'
    self%time_order=0
  end if
  call lock%set ('timestep')
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
  call lock%unset ('timestep')
  call trace_end (itimer)
END SUBROUTINE update

!===============================================================================
SUBROUTINE update_single (self, id, iit, dt, mem, dmem, mesh, lock, time_order)
  class(timestep_t):: self
  class(mesh_t):: mesh(:)
  type(lock_t):: lock
  integer:: id, nt
  integer, dimension(:):: iit
  real(8), dimension(:):: dt
  real(kind=KindScalarVar), dimension(:,:,:,:):: mem, dmem       ! memory buffers
  integer, intent(in):: time_order
  integer:: l(3), u(3)
  real:: wt(nslice)
  integer, save:: itimer=0
  !......................................................................
  call trace_begin ('timestep_t%update_single', itimer=itimer)
  if (time_order <= 0) return

  nt = size(iit)                                ! number of time slots
  l = mesh%lb
  u = mesh%ub
  where (mesh%lower_boundary)
    l = mesh%lb
  end where
  where (mesh%upper_boundary)
    u = mesh%ub
  end where

  select case (time_order)
  case (0)
  case (1)
    wt = [   1.0,     0.0,      0.0,    0.0]       ! 1st order
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt  )) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1)) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))*(wt(1)*dt(iit(nt-1)))
  case (2)
    wt = [   1.5,     -.5,      0.0,    0.0]       ! 2nd order
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt  )) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1)) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))*(wt(1)*dt(iit(nt-1))) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-2))*(wt(2)*dt(iit(nt-2)))
  case (3)
    wt = [23./12., -4./3.,   5./12.,    0.0]       ! 3rd order
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt  )) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1)) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))*(wt(1)*dt(iit(nt-1))) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-2))*(wt(2)*dt(iit(nt-2))) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-3))*(wt(3)*dt(iit(nt-3)))
  case (4)
    wt = [55./24., -59./24., 37/24., -3./8.]       ! 4th order
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt  )) = &
    mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1)) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))*(wt(1)*dt(iit(nt-1))) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-2))*(wt(2)*dt(iit(nt-2))) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-3))*(wt(3)*dt(iit(nt-3))) + &
    dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-4))*(wt(4)*dt(iit(nt-4)))
  end select
  if (self%verbose>3) then
    print'(a,i4,2x,a,1p,6e10.2,2x,6i5)', 'timestep, id:', &
      id, 'minmax:', &
      minval(dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))), &
      maxval(dmem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))), &
      minval(mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))), &
      maxval(mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt-1))), &
      minval(mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt  ))), &
      maxval(mem(l(1):u(1),l(2):u(2),l(3):u(3),iit(nt  ))), shape(mem)
  end if
  call trace_end (itimer)
END SUBROUTINE update_single

END MODULE timestep_mod
