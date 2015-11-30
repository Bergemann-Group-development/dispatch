!> =============================================================================
!> $Id$
!> Adams-Bashforth time stepping, with constant time step assumed (for now)
!> =============================================================================
MODULE timestep_mod
  USE io_mod
  USE pde_mod
  implicit none
  private
  integer, parameter:: nslice = 4
  type, public, extends(pde_t):: timestep_t
    real:: wt(nslice)
    integer:: timeorder
  contains
    procedure:: init   => init_timestep
    procedure:: update => update_timestep
  end type
  integer:: verbose = 1
  integer:: timeorder = 2
  logical, save:: first_time = .true.
CONTAINS

!> =============================================================================
!> Initialize coefficients for Adams-Bashforth time stepping
!> =============================================================================
SUBROUTINE init_timestep (self, dims, position, size)
  class(timestep_t):: self
  integer:: dims(5)
  real:: position(3), size(3)
  real:: courant = 0.4
  namelist /timestep_params/ courant, timeorder, verbose
!........................................................................
  !$omp critical (timest_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input,timestep_params)
    if (io%master) write (*,timestep_params)
  end if
  !$omp end critical (timest_cr)
  select case (timeorder)
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
  self%timeorder = timeorder
  self%verbose   = verbose
END SUBROUTINE init_timestep

!************************************************************************
SUBROUTINE update_timestep (self)
  class(timestep_t):: self
  integer:: it, jt, kt
  call self%pde_t%update                          ! time derivative at it
  it = modulo(self%it,self%nt)+1             ! wrap new time slot
  do jt=1,self%timeorder                        ! loop over time
    kt = modulo(it-jt,self%nt)               ! previous time slots
    if (jt == 1) then
      self%mem(:,:,:,:,it,1) = &
      self%mem(:,:,:,:,kt,2)*(self%wt(jt)*self%dtime)
    else
      self%mem(:,:,:,:,it,1) = &
      self%mem(:,:,:,:,it,1) + &
      self%mem(:,:,:,:,kt,2)*(self%wt(jt)*self%dtime)
    end if
  end do
  self%time = self%time + self%dtime            ! update self time
  self%t(it) = self%time                        ! update slot time
  self%it = it                                  ! update time slot
!........................................................................
END SUBROUTINE update_timestep

END MODULE timestep_mod
