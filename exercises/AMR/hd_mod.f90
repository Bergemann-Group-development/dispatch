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
  implicit none
  private
  type, public, extends(task_t):: mhd_t
    type(scalar_t):: d, s
    type(vector_t):: p
    real:: gamma=5./3.
  contains
    procedure:: init
    procedure:: pde
  end type
  type(mhd_t):: mhd
  type(scalar_t):: pg, ld
  type(vector_t):: u
CONTAINS
!> =============================================================================
!> Make sure the memory is allocated, then map variables to the current time slot
!> The 'it' slot must be the one we are to update, with new values going into the
!> 'new' time slot
!> =============================================================================
SUBROUTINE init (self)
  type(mhd_t):: self
  integer:: it
  !.............................................................................
  call self%task_t%init
  d%f    => task%mem(:,:,:,i_d ,task%it,1)
  s%f    => task%mem(:,:,:,i_s ,task%it,1)
  p%x    => task%mem(:,:,:,i_px,task%it,1)
  p%y    => task%mem(:,:,:,i_py,task%it,1)
  p%z    => task%mem(:,:,:,i_pz,task%it,1)
  dddt%f => task%mem(:,:,:,i_d ,task%it,2)
  dsdt%f => task%mem(:,:,:,i_s ,task%it,2)
  dpdt%x => task%mem(:,:,:,i_px,task%it,2)
  dpdt%y => task%mem(:,:,:,i_py,task%it,2)
  dpdt%z => task%mem(:,:,:,i_pz,task%it,2)
END SUBROUTINE init
!> =============================================================================
!> Compute the time derivatives, using the symbolic names mapped into task%mem
!> above
!> =============================================================================
SUBROUTINE pde (self)
  type(mhd_t):: self
  integer:: it
  !.............................................................................
associate (d=>self%d, p=>self%p, s=>self%s, gsmma=>self%gamma}
  pg = d**gamma*exp((s/d)*(gamma-1.0))
  ld = dn(log(d))
  u  = p/exp(ld)
  self%dddt = -div(p)
  self%dpdt = -div(d*up(u)**2) - grad(pg)
  self%dsdt = -div(dn(s)*p)
end associate
END SUBROUTINE pde
!> =============================================================================
END MODULE pde_mod
