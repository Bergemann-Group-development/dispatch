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
  implicit none
  private
  integer, parameter:: i_d=1, i_s=2, i_px=3, i_py=4, i_pz=5
  type, public, extends(patch_t):: pde_t
    type(scalar_t):: d, dddt
    type(scalar_t):: s, dsdt
    type(vector_t):: p, dpdt
    real:: gamma=5./3.
  contains
    procedure:: init => init_hd
    procedure:: pde
    procedure:: log_density
    procedure:: log_pressure
    procedure:: velocity_magnitude
    procedure:: magnetic_field_magnitude
  end type
  type(pde_t):: mhd
  type(scalar_t):: pg
  type(vector_t):: u, ld, uup
  type(tensor_t):: uu
CONTAINS

!> =============================================================================
!> Make sure the memory is allocated, then map variables to the current time slot
!> The 'it' slot must be the one we are to update, with new values going into the
!> 'new' time slot
!> =============================================================================
SUBROUTINE init_hd (self, position, size)
  class(pde_t):: self
  real:: position(3), size(3)
  type(tensor_t):: tensor
  !.............................................................................
  call trace_begin ('hd_mod:pde_t%init', 2)
  call self%patch_t%init (position, size)
  !-----------------------------------------------------------------------------
  ! Memory mapping
  !-----------------------------------------------------------------------------
  call self%d%init (self%gn, size/self%n)
  call self%p%init (self%gn, size/self%n)
  call tensor%init (self%gn, size/self%n)
  self%d%f    => self%mem(:,:,:,i_d ,self%it,1)
  call trace ('1')
  self%s%f    => self%mem(:,:,:,i_s ,self%it,1)
  self%p%x    => self%mem(:,:,:,i_px,self%it,1)
  self%p%y    => self%mem(:,:,:,i_py,self%it,1)
  self%p%z    => self%mem(:,:,:,i_pz,self%it,1)
  self%dddt%f => self%mem(:,:,:,i_d ,self%it,2)
  self%dsdt%f => self%mem(:,:,:,i_s ,self%it,2)
  self%dpdt%x => self%mem(:,:,:,i_px,self%it,2)
  self%dpdt%y => self%mem(:,:,:,i_py,self%it,2)
  self%dpdt%z => self%mem(:,:,:,i_pz,self%it,2)
  !-----------------------------------------------------------------------------
  ! Physical initial conditions
  !-----------------------------------------------------------------------------
  call trace ('2')
  self%d%id = 1
  self%s%id = 2
  self%p%id = 3
  self%d = 1.0
  self%s = 0.01
  self%p = 0.0
  print*,'hd',self%id, self%it, i_d
  print*,'hd',self%li,self%ui
  print*,'hd',self%mem(1:6,1,1,1,1,1)
  call trace_end
END SUBROUTINE init_hd

!> =============================================================================
!> Compute the time derivatives, using the symbolic names mapped into task%mem
!> =============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(pde_t):: self
  integer:: it
  type(vector_t):: u
  type(tensor_t):: pu
  !.............................................................................
  call  u%allocate
  call pu%allocate
  call trace_begin ('hd_t%pde')
associate (d=>self%d, p=>self%p, s=>self%s, gamma=>self%gamma)
  pg = d**gamma*exp(s/d*(gamma-1.0))
  u  = p/d
  pu = d*outer(u,u)
  self%dddt = -div(p)
  self%dpdt = -div(pu) - grad(pg)
  self%dsdt = -div(p*s)
  call trace_end
end associate
END SUBROUTINE pde

!> =============================================================================
!> Compute the log density -- for use in refinement
!> =============================================================================
SUBROUTINE log_density (self, out)
  class(pde_t):: self
  real, dimension(self%n(1),self%n(2),self%n(3)):: out
  !...............................................................................
  call trace_begin ('hd_mod:pde_t%log_density')
  if (io%verbose>3) then
    print*,'log_d',self%id, self%it, i_d
    print*,'log_d',self%li,self%ui
    print *,self%mem(1:5,1,1,1,1,1)
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
