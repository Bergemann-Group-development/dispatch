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
  integer, parameter:: i_d=1, i_phi=2, i_px=3, i_py=4, i_pz=5
  type, public, extends(ghd_t):: pde_t
    real:: w_BE
    real:: p_BE
    real:: a_BE
    real:: amplitude
  contains
    procedure:: init
  end type
CONTAINS

!> =============================================================================
!> Setup a patch with density and gravitational potential for the Poisson solver
!> =============================================================================
SUBROUTINE init (self)
  class(pde_t):: self
  integer:: i, j, k
  real:: x, y, z
  logical:: first_time = .true.
  real:: w_BE=5.0, p_BE=1.0, a_BE=0.10, amplitude=10., r
  namelist /BE_params/ w_BE, p_BE, a_BE, amplitude
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  experimentname = 'BE'
  call io%init (experimentname)
  call self%ghd_t%init
  call trace_begin ('hd_mod:pde_t%init', 2)
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
      self%d%f(i,j,k) = &
        1.0 + self%amplitude*exp(-r**self%p_BE/(1.+r**self%p_BE/log(r**2/self%a_BE)))
    end do
    end do
    end do
  end if
  !-----------------------------------------------------------------------------
  ! Initial and boundary conditions for phi
  !-----------------------------------------------------------------------------
  self%phi%f = 0.0
  !
  call trace_end
END SUBROUTINE init
!> =============================================================================
END MODULE pde_mod
