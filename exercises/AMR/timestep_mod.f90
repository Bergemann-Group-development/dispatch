!> =============================================================================
!> $Id$
!> Adams-Bashforth time stepping, with constant time step assumed (for now)
!> =============================================================================
MODULE timestep_mod
  USE io_mod
  USE trace_mod
  USE pde_mod
  implicit none
  private
  integer, parameter:: nslice = 4
  type, public, extends(pde_t):: timestep_t
    real:: wt(nslice)
    integer:: time_order
  contains
    procedure:: init   => init_timestep
    procedure:: update => update_timestep
  end type
  real:: courant = 0.4
  integer:: verbose = 1
  integer:: time_order = 2
  logical, save:: first_time = .true.
CONTAINS

!> =============================================================================
!> Initialize coefficients for Adams-Bashforth time stepping
!> =============================================================================
SUBROUTINE init_timestep (self, position, size)
  class(timestep_t):: self
  real:: position(3), size(3)
  namelist /timestep_params/ courant, time_order, verbose
!........................................................................
  call trace_begin ('timestep_t%init', 1)
  call self%pde_t%init (position, size)
  print*,'ts',self%id
  print*,'ts',self%mem(1:6,1,1,1,1,1)
  !$omp critical (timest_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input,timestep_params)
    if (io%master) write (*,timestep_params)
  end if
  !$omp end critical (timest_cr)
  select case (time_order)
  case (1)
    self%wt = [   1.0,     0.0,      0.0,    0.0]   ! 1st order
  case (2)
    self%wt = [   1.5,     -.5,      0.0,    0.0]   ! 2nd order
  case (3)
    self%wt = [23./12., -4./3.,   5./12.,    0.0]   ! 3rd order
  case (4)
    self%wt = [55./24., -59./24., 37/24., -3./8.]   ! 4th order
  end select
  self%courant   = courant
  self%time_order = time_order
  self%verbose   = verbose
  call trace_end
END SUBROUTINE init_timestep

!************************************************************************
SUBROUTINE update_timestep (self)
  class(timestep_t):: self
  integer:: it, jt, kt
  !......................................................................
  call trace_begin ('timestep_t%update', 3)
  call self%pde_t%update                        ! time derivative at it
  it = modulo(self%it,self%nt)+1                ! wrap to new time slot
  if (io%verbose>3) then
    print*, 'timestep_t%update, id =', self%id, self%it, it
    print*, 'timestep_t%update, d =', self%mem(1:5,1,1,1,1,1)
  end if
  do jt=1,self%time_order                       ! loop over time
    kt = modulo(it-jt,self%nt)                  ! previous time slots
    if (jt == 1) then
      self%mem(:,:,:,:,it,1) = &
      self%mem(:,:,:,:,kt,2)*(self%wt(jt)*self%dtime)
    else
      self%mem(:,:,:,:,it,1) = &
      self%mem(:,:,:,:,it,1) + &
      self%mem(:,:,:,:,kt,2)*(self%wt(jt)*self%dtime)
    end if
  end do
  if (io%verbose>3) then
    print*, 'timestep_t%update, t =', self%time, self%dtime
    print*, 'timestep_t%update, d =', self%mem(1:5,1,1,1,1,1)
  end if
  self%time = self%time + self%dtime            ! update self time
  self%t(it) = self%time                        ! update slot time
  self%it = it                                  ! update time slot
  call trace_end
!........................................................................
END SUBROUTINE update_timestep

END MODULE timestep_mod
