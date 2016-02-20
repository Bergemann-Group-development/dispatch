!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE experiment_mod
  USE io_mod
  USE omp_timer_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  USE laplace_mod
  USE ghd_mod
  USE boundary_mod
  USE timestep_mod
  implicit none
  private
  type, public, extends(ghd_t):: experiment_t
    real:: w_BE
    real:: p_BE
    real:: a_BE
    real:: r_BE
    real:: amplitude
    type(timestep_t):: timestep
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!> =============================================================================
!> Setup a patch with density and gravitational potential for the Poisson solver
!> =============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer:: i, j, k
  integer, save:: levelmax=0
  real:: x, y, z
  logical:: first_time = .true.
  real:: w_BE=5.0, p_BE=1.0, a_BE=0.10, r_BE=0.45, amplitude=10., r
  namelist /experiment_params/ w_BE, p_BE, a_BE, r_BE, amplitude, levelmax
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  call io%init
  call trace_begin ('experiment_t%init')
  call self%ghd_t%init
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, experiment_params)
    if (io%master) write (*,experiment_params)
  end if
  self%w_BE      = w_BE
  self%p_BE      = p_BE
  self%a_BE      = a_BE
  self%r_BE      = r_BE
  boundary%r_BE  = r_BE
  self%amplitude = amplitude
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Approximate a Bonnor-Ebert sphere
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do k=1,self%gn(3)
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5)*self%ds(1) + self%position(1)
      y = (j-(self%gn(2)+1)*0.5)*self%ds(2) + self%position(2)
      z = (k-(self%gn(3)+1)*0.5)*self%ds(3) + self%position(3)
      r = sqrt(x**2+y**2+z**2)/self%w_BE
      r = min(r,r_BE/w_BE)
      self%d%f(i,j,k) = &
        self%amplitude*exp(-r**self%p_BE/(1.+r**self%p_BE/log(r**2/self%a_BE)))
    end do
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Initial condition for phi
  !-----------------------------------------------------------------------------
  self%phi%f = 0.0
  !
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Enforce constant values outside a radius r=r_BE (not equalt to the BE-radius)
!> =============================================================================
SUBROUTINE update(self)
  class(experiment_t):: self
  integer:: i, j, k, ii(3), iv
  real:: x, y, z, r, ntot, nout
  !.............................................................................
  call trace_begin('experiment_t%update')
associate (li=>self%li,ui=>self%ui,gn=>self%gn,ds=>self%ds,pos=>self%position)
  !-----------------------------------------------------------------------------
  ! Set phi=0 at and outside the boundary
  !-----------------------------------------------------------------------------
  do k=1,gn(3)
  do j=1,gn(2)
  do i=1,gn(1)
    x = (i-(gn(1)+1)*0.5)*ds(1) + pos(1)
    y = (j-(gn(2)+1)*0.5)*ds(2) + pos(2)
    z = (k-(gn(3)+1)*0.5)*ds(3) + pos(3)
    r = sqrt(x**2+y**2+z**2)
    if (r > self%r_BE) then
      self%mem(i,j,k,self%i_phi,self%it,1) = 0.0
    end if
  end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! This where d(f)/dt is computed, so zap to zero outside the boundary, treating
  ! staggered components separately
  !-----------------------------------------------------------------------------
  call self%ghd_t%update
  ntot = 0.0
  nout = 0.0
  do k=1,gn(3)
  do j=1,gn(2)
  do i=1,gn(1)
    ntot = ntot+1.0
    x = (i-(gn(1)+1)*0.5)*ds(1) + pos(1)
    y = (j-(gn(2)+1)*0.5)*ds(2) + pos(2)
    z = (k-(gn(3)+1)*0.5)*ds(3) + pos(3)
    r = sqrt(x**2+y**2+z**2)
    if (r > self%r_BE) then
      nout = nout+1.0
      self%mem(i,j,k,self%i_d  ,self%it,2) = 0.0
      self%mem(i,j,k,self%i_s  ,self%it,2) = 0.0
    end if
    x = (i-(gn(1)+2)*0.5)*ds(1) + pos(1)
    r = sqrt(x**2+y**2+z**2)
    if (r > self%r_BE) then
      self%mem(i,j,k,self%i_px,self%it,2) = 0.0
    end if
    x = (i-(gn(1)+1)*0.5)*ds(1) + pos(1)
    y = (j-(gn(2)+2)*0.5)*ds(2) + pos(2)
    r = sqrt(x**2+y**2+z**2)
    if (r > self%r_BE) then
      self%mem(i,j,k,self%i_py,self%it,2) = 0.0
    end if
    y = (j-(gn(2)+1)*0.5)*ds(2) + pos(2)
    z = (k-(gn(3)+2)*0.5)*ds(3) + pos(3)
    r = sqrt(x**2+y**2+z**2)
    if (r > self%r_BE) then
      self%mem(i,j,k,self%i_pz,self%it,2) = 0.0
    end if
  end do
  end do
  end do
  !print*,'fraction outside:',nout/ntot, pos, ds
end associate
  !-----------------------------------------------------------------------------
  ! Call timestep update procedure to fill the next time slot with values
  !-----------------------------------------------------------------------------
  if (io%verbose > 1) then
    print'(a,7i10)',   'experiment_t%update: iit =', self%iit
    print'(a,7f10.6)', 'experiment_t%update:   t =', self%t(self%iit)
  end if
  call self%courant_condition
  call self%timestep%update (self%id, self%iit, self%dt, self%mem)
  !-----------------------------------------------------------------------------
  ! Periodic wrapping, time extrapolation, time slot rotation
  !-----------------------------------------------------------------------------
  if (self%id==1) then
    call self%periodic_grid                     ! periodic root grid
  else
    call self%extrapolate_guards                ! time-extrapolate guard cells
  end if
  call self%patch_t%update                      ! rotate time slots
  !
  call trace_end
END SUBROUTINE update
!> =============================================================================
END MODULE experiment_mod
