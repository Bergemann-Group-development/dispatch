!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the Poisson problem in 3-D
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE io_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_phi=2
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, phi
    real:: gamma=5./3.
    real:: w_gauss
    real:: p_gauss
    real:: amplitude
  contains
    procedure:: init
    procedure:: poisson
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
  real:: w_gauss=5.0, p_gauss=1.0, amplitude=1.0
  namelist /poisson_params/ w_gauss, p_gauss, amplitude
  !-----------------------------------------------------------------------------
  ! Read standard parameters and copy to structure
  !-----------------------------------------------------------------------------
  experimentname = 'poisson'
  call io%init (experimentname)
  call trace_begin ('hd_mod:pde_t%init', 2)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, poisson_params)
    if (io%master) write (*,poisson_params)
  end if
  self%size      = 1.
  self%position  = 0.
  self%w_gauss   = w_gauss
  self%p_gauss   = p_gauss
  self%amplitude = amplitude
  self%nv        = 2
  call self%patch_t%init
  !-----------------------------------------------------------------------------
  ! Memory mapping og density and potentiall field
  !-----------------------------------------------------------------------------
  call    self%d%init (self%gn, self%size/self%n, 'd')
  call  self%phi%init (self%gn, self%size/self%n, 'phi')
  self%d%f   => self%mem(:,:,:,i_d   ,self%it,1)
  self%phi%f => self%mem(:,:,:,i_phi ,self%it,1)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: gaussian top hat
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do k=1,self%gn(3)
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5) + self%position(1)/self%ds(1)
      y = (j-(self%gn(2)+1)*0.5) + self%position(2)/self%ds(2)
      y = (j-(self%gn(3)+1)*0.5) + self%position(3)/self%ds(3)
      self%d%f(i,j,k) = &
        exp(self%amplitude*exp(-sqrt(x**2+y**2)/self%w_gauss**self%p_gauss))
    end do
    end do
    end do
  end if
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Compute the time derivatives, using the symbolic names mapped into task%mem,
!> and an associate block to make the code maximally readable
!> =============================================================================
SUBROUTINE poisson (self)
  class(pde_t):: self
  !.............................................................................
  call trace_begin ('poisson_mod:pde_t%update', 2)
  !...
  call trace_end
END SUBROUTINE poisson

!> =============================================================================
END MODULE pde_mod
