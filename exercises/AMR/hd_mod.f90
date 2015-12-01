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
  integer, parameter:: i_d=1, i_s=2, i_px=3, i_py=4, i_pz=5
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, dddt
    type(scalar_t):: s, dsdt
    type(vector_t):: p, dpdt
    real:: gamma=5./3.
  contains
    procedure:: init
    procedure:: init_values
    procedure:: update => update_pde
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  type(pde_t):: mhd
CONTAINS

!> =============================================================================
!> Map variables to the current time slot. The 'it' slot is the one we are to
!> update, with new values going into the 'new' time slot
!> =============================================================================
SUBROUTINE init (self, position, size)
  class(pde_t):: self
  real:: position(3), size(3)
  type(tensor_t):: tensor
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init', 2)
  call self%patch_t%init (position, size)
  !-----------------------------------------------------------------------------
  ! Memory mapping
  !-----------------------------------------------------------------------------
  call    self%d%init (self%gn, self%size/self%n, self%derivs, 'd')
  call    self%s%init (self%gn, self%size/self%n, self%derivs, 's')
  call    self%p%init (self%gn, self%size/self%n, self%derivs, 'p')
  call self%dddt%init (self%gn, self%size/self%n, self%derivs, 'dddt')
  call self%dsdt%init (self%gn, self%size/self%n, self%derivs, 'dsdt')
  call self%dpdt%init (self%gn, self%size/self%n, self%derivs, 'dpdt')
  self%d%f    => self%mem(:,:,:,i_d ,self%it,1)
  self%s%f    => self%mem(:,:,:,i_s ,self%it,1)
  self%p%x    => self%mem(:,:,:,i_px,self%it,1)
  self%p%y    => self%mem(:,:,:,i_py,self%it,1)
  self%p%z    => self%mem(:,:,:,i_pz,self%it,1)
  self%dddt%f => self%mem(:,:,:,i_d ,self%it,2)
  self%dsdt%f => self%mem(:,:,:,i_s ,self%it,2)
  self%dpdt%x => self%mem(:,:,:,i_px,self%it,2)
  self%dpdt%y => self%mem(:,:,:,i_py,self%it,2)
  self%dpdt%z => self%mem(:,:,:,i_pz,self%it,2)
  call check_dims (self%dddt%f, 'init')
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Make sure the memory is allocated, then map variables to the current time slot
!> The 'it' slot must be the one we are to update, with new values going into the
!> 'new' time slot
!> =============================================================================
SUBROUTINE init_values (self, position, size)
  class(pde_t):: self
  real:: position(3), size(3)
  type(tensor_t):: tensor
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init_values', 2)
  call self%init (position, size)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions
  !-----------------------------------------------------------------------------
  self%s = 0.01
  self%d = 1.0
  self%p = 0.0
  print*,'patch:', self%id, '  min, max:', minval(self%d), maxval(self%d)
  call trace_end
END SUBROUTINE init_values

!> =============================================================================
!> Compute the time derivatives, using the symbolic names mapped into task%mem,
!> and an associate block to make the code maximally readable
!> =============================================================================
SUBROUTINE update_pde (self)
  USE scalar_mod
  class(pde_t):: self
  type(scalar_t), save:: pg
  type(vector_t), save:: u
  type(tensor_t), save:: pu
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%update', 2)
  call pg%init (self%gn, self%size/self%n, self%derivs, 'pg')
  call  u%init (self%gn, self%size/self%n, self%derivs, 'u' )
  call pu%init (self%gn, self%size/self%n, self%derivs, 'pu')
  associate (d=>self%d, p=>self%p, s=>self%s, gamma=>self%gamma)
    call trace('self%d**gamma')  ! OK
    pg = self%d**self%gamma
    call trace('s/d')       ! NOK
    pg = s/d
    call trace('s/d*(gamma-1)')
    pg = (s/d)*(gamma-1.0)
    call trace('pg')
    pg = d**gamma*exp((s/d)*(gamma-1.0))
    call trace('u=p/d')
    u = p/d
    pu = outer(p,u)
    call trace('dddt')
    self%dddt = -div(p)
    call trace('dpdt')
    self%dpdt = -div(pu)
    call trace('dpdt')
    self%dpdt = -div(pu) - grad(pg)
    call trace('dsdt')
    self%dsdt = -div(p*s)
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
    print*, 'log_d', self%id, self%it, self%new, i_d
    print*, 'log_d', self%li, self%ui, self%n
    print *, self%mem(1:5,1,1,1,1,1)
    print *, 'min,max:', minval(self%d%f), maxval(self%d%f)
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
  out = log(self%mem(self%li(1):self%ui(1), &
                     self%li(2):self%ui(2), &
                     self%li(3):self%ui(3), i_d, self%it, 1))
END SUBROUTINE log_pressure

!> =============================================================================
SUBROUTINE velocity_magnitude (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  out = log(self%mem(self%li(1):self%ui(1), &
                     self%li(2):self%ui(2), &
                     self%li(3):self%ui(3), i_d, self%it, 1))
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
