!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the induction problem in 3-D, for now with the veloocity field constant
!> in space
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE pde_mod
  USE task_mod
  USE scalar_mod
  USE vector_mod
  USE tensor_mod
  USE stagger_mod
  implicit none
  private
  integer, parameter:: i_d=1, i_s=2, i_px=3, i_py=4, i_pz=5
  type, public, extends(task_t):: pde_t
    type(scalar_t):: d, dddt
    type(scalar_t):: s, dsdt
    type(vector_t):: p, dpdt
    real:: gamma=5./3.
  contains
    procedure:: init => init_hd
    procedure:: pde
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
SUBROUTINE init_hd (self, dims, position, size)
  class(pde_t):: self
  integer:: dims(5)
  real:: position(3), size(3)
  !.............................................................................
  call self%task_t%init (dims, position, size)
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
END SUBROUTINE init_hd
!> =============================================================================
!> Compute the time derivatives, using the symbolic names mapped into task%mem
!> above
!> =============================================================================
SUBROUTINE pde (self)
  USE scalar_mod
  class(pde_t):: self
  integer:: it
  type(scalar_t):: tmp
  type(vector_t):: vec
  !.............................................................................
associate (d=>self%d, p=>self%p, s=>self%s, gamma=>self%gamma)
  tmp = self%s/self%d*(gamma-1.0)
  pg = d**gamma*exp(tmp)
  ld = dnv(log(d))
  u  = p/exp(ld)
  self%dddt = -div(p)
  uup = upv(u)
  uu = outer(uup,uup)
  uu = d*uu
  self%dpdt = -div(uu) - grad(pg)
  vec = dnv(s)
  vec = vec*p
  self%dsdt = -div(vec)
end associate
END SUBROUTINE pde
!> =============================================================================
END MODULE pde_mod
