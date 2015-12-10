!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE io_mod
  USE omp_timer_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  USE laplace_mod
  USE ghd_mod
  implicit none
  private
  type, public, extends(ghd_t):: pde_t
    real:: w_BE
    real:: p_BE
    real:: a_BE
    real:: r_BE
    real:: amplitude
  contains
    procedure:: init
    procedure:: update
  end type
CONTAINS

!> =============================================================================
!> Setup a patch with density and gravitational potential for the Poisson solver
!> =============================================================================
SUBROUTINE init (self, position, size)
  class(pde_t):: self
  real:: position(3), size(3)
  integer:: i, j, k
  real:: x, y, z
  logical:: first_time = .true.
  real:: w_BE=5.0, p_BE=1.0, a_BE=0.10, r_BE=0.45, amplitude=10., r
  namelist /BE_params/ w_BE, p_BE, a_BE, r_BE, amplitude
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  !experimentname = 'BE'
  call io%init
  call trace_begin ('BE_mod::pde_t%init')
  call self%ghd_t%init (position, size)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, BE_params)
    if (io%master) write (*,BE_params)
  end if
  self%size      = 1.
  self%position  = 0.
  self%w_BE      = w_BE
  self%p_BE      = p_BE
  self%a_BE      = a_BE
  self%r_BE      = r_BE
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
  ! Cheat refinement into making Rubiks Cubes, by insisting that the central
  ! density is huge, thus triggering Jeans' length refinement
  !-----------------------------------------------------------------------------
  if (self%level<io%levelmax .and. all(self%position==0.0)) then
    self%d%f(self%gn(1)/2,self%gn(2)/2,self%gn(3)/2) = 1e8
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
  class(pde_t):: self
  integer:: i, j, k, ii(3), iv
  real:: x, y, z, r, ntot, nout
  !.............................................................................
  call trace_begin('BE_mod::pde_t%update')
associate (li=>self%li,ui=>self%ui,gn=>self%gn,ds=>self%ds,pos=>self%position)
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
  ! This where d(f)/dt is computed, so zap to zero outside r=r_BE
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
  call trace_end
END SUBROUTINE update
!> =============================================================================
END MODULE pde_mod
