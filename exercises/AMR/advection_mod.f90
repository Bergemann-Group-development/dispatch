!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the induction problem in 3-D, for now with the veloocity field constant
!> in space
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE trace_mod
  USE patches_mod
  USE scalar_mod
  USE vector_mod
  USE tensor_mod
!  USE upwind_mod, only: d
  implicit none
  private
  integer, parameter:: i_d=1
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, dddt
    real:: u(3)=(/1.0,0.,0.1/)
    real:: gamma=5./3.
    real:: w_gauss
    real:: p_gauss
    real:: u_max
    real:: amplitude
  contains
    procedure:: init
    procedure:: update => update_pde
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  type(pde_t):: mhd
  class(central_t), pointer:: deriv
  public experiment_name
CONTAINS

!> =============================================================================
FUNCTION experiment_name()
  character(len=64):: experiment_name
  experimentname = 'advection'                  ! shared from patches_mod
  experiment_name = experimentname
END FUNCTION
!> =============================================================================
!> Map variables to the current time slot. The 'it' slot is the one we are to
!> update, with new values going into the 'new' time slot
!> =============================================================================
SUBROUTINE init (self, position, size)
  class(pde_t):: self
  real:: position(3), size(3)
  type(tensor_t):: tensor
  integer:: i, j
  real:: x, y
  logical:: first_time = .true.
  real:: w_gauss=5.0, p_gauss=1.0, amplitude=1.0
  namelist /advection_params/ w_gauss, p_gauss, amplitude
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init', 2)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, advection_params)
    if (io%master) write (*,advection_params)
  end if
  self%w_gauss = w_gauss
  self%p_gauss = p_gauss
  self%amplitude = amplitude
  !-----------------------------------------------------------------------------
  ! Physical initial conditions: Gaussian top hat.  Only at time=0.0
  !-----------------------------------------------------------------------------
  if (self%time==0.0) then
    do j=1,self%gn(2)
    do i=1,self%gn(1)
      x = (i-(self%gn(1)+1)*0.5) + self%position(1)/self%ds(1)
      y = (j-(self%gn(2)+1)*0.5) + self%position(2)/self%ds(2)
      self%mem(i,j,:,1,1,1) = &
        1.0 + self%amplitude*exp(-(x**2+y**2)/self%w_gauss**2)**self%p_gauss
    end do
    end do
    self%nv = 1
    if (io%verbose>1) then
      print'(i4,2x,3i4,1p,6e10.2)', &
        self%id,self%n,self%w_gauss,self%p_gauss,self%position
    else
      print'(1x,a,i4,2(3f8.3,2x),a,1p,2e11.3)', &
        'initial:', self%id, self%position, self%size, &
        '  min, max:', &
        minval(self%mem(:,:,:,1,1,1)), &
        maxval(self%mem(:,:,:,1,1,1))
    end if
    self%u_max = sqrt(sum(self%u**2))
  end if
  !-----------------------------------------------------------------------------
  ! Memory mapping; this is cheap, and may be done every time
  !-----------------------------------------------------------------------------
  allocate (deriv)
  call     deriv%init (self%gn, self%size/self%n)
  call    self%d%init (self%gn, self%size/self%n, 'd')
  call self%dddt%init (self%gn, self%size/self%n, 'dddt')
  self%d%f    => self%mem(:,:,:,i_d ,self%it,1)
  self%dddt%f => self%mem(:,:,:,i_d ,self%it,2)
  call check_dims (self%dddt%f, 'init')
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Compute the time derivatives, using the symbolic names mapped into task%mem,
!> and an associate block to make the code maximally readable
!> =============================================================================
SUBROUTINE update_pde (self)
  USE scalar_mod
  class(pde_t):: self
  type(vector_t), save:: p
  !.............................................................................
  call trace_begin ('advection_mod:pde_t%update', 2)
  call p%init (self%gn, self%size/self%n, 'p')
  self%dddt%f => self%mem(:,:,:,1,self%it,2)
  associate (d=>self%d,u=>self%u,dddt=>self%dddt)
    p = u
    p = p*d
    dddt%f = -deriv%dx(p%x,self%gn)-deriv%dy(p%y,self%gn)
    !call minus_div (dddt, p)
    if (io%verbose>2) &
    print'(1x,a,1p,i4,4(2x,2e10.2))', 'update_pde', &
      self%id,minval(self%dddt),maxval(self%dddt), &
              minval(self%dddt%f),maxval(self%dddt%f), &
              minval(self%mem(:,:,:,1,self%it,2)), &
              maxval(self%mem(:,:,:,1,self%it,2)), &
              minval(self%mem(:,:,:,1,self%new,2)), &
              maxval(self%mem(:,:,:,1,self%new,2))
  end associate
  call trace_end
END SUBROUTINE update_pde

!> =============================================================================
!> Compute the log density -- for use in refinement
!> =============================================================================
SUBROUTINE log_density (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  call trace_begin ('hd_mod:pde_t%log_density')
  if (io%verbose>3) then
    print'(1x,a,4i4,3(2x,3i4),3x,a,1p,4e10.2)', &
      'log_d', self%id, self%it, self%new, i_d, self%li, self%ui, self%n, &
      'minmax:', minval(self%d%f(self%li(1):self%ui(1), &
                                 self%li(2):self%ui(2), &
                                 self%li(3):self%ui(3))),&
                 maxval(self%d%f(self%li(1):self%ui(1), &
                                 self%li(2):self%ui(2), &
                                 self%li(3):self%ui(3))), &
      minval(self%mem(self%li(1):self%ui(1), &
                      self%li(2):self%ui(2), &
                      self%li(3):self%ui(3), i_d, self%it, 1)), &
      maxval(self%mem(self%li(1):self%ui(1), &
                      self%li(2):self%ui(2), &
                      self%li(3):self%ui(3), i_d, self%it, 1))
  end if
  out = log(self%mem(self%li(1):self%ui(1), &
                     self%li(2):self%ui(2), &
                     self%li(3):self%ui(3), i_d, self%it, 1))
  call trace_end
END SUBROUTINE log_density

!> =============================================================================
SUBROUTINE log_pressure (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE log_pressure

!> =============================================================================
SUBROUTINE velocity_magnitude (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE velocity_magnitude

!> =============================================================================
SUBROUTINE magnetic_field_magnitude (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = 0.0
END SUBROUTINE magnetic_field_magnitude

!> =============================================================================
END MODULE pde_mod
