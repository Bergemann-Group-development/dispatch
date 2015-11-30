!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> Solve the induction problem in 3-D, for now with the veloocity field constant
!> in space
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE mhd_mod
  USE scalar_mod
  USE vector_mod
  USE tensor_mod
  implicit none
  private
  type, public, extends(task_t):: mhd_t
    type(scalar_t):: d, s
    type(vector_t):: p
  contains
    procedure:: init
    procedure:: update
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
  b%x    => task%mem(:,:,:,i_bx,task%it,1)
  b%y    => task%mem(:,:,:,i_by,task%it,1)
  b%z    => task%mem(:,:,:,i_bz,task%it,1)
  dddt%f => task%mem(:,:,:,i_d ,task%it,2)
  dsdt%f => task%mem(:,:,:,i_s ,task%it,2)
  dpdt%x => task%mem(:,:,:,i_px,task%it,2)
  dpdt%y => task%mem(:,:,:,i_py,task%it,2)
  dpdt%z => task%mem(:,:,:,i_pz,task%it,2)
  dbdt%x => task%mem(:,:,:,i_bx,task%it,2)
  dbdt%y => task%mem(:,:,:,i_by,task%it,2)
  dbdt%z => task%mem(:,:,:,i_bz,task%it,2)
END SUBROUTINE init
!> =============================================================================
!> Make sure the memory is allocated, then map variables to the current time slot
!> The 'it' slot must be the one we are to update, with new values going into the
!> 'new' time slot
!> =============================================================================
SUBROUTINE pde (self)
  type(mhd_t):: self
  integer:: it
  !.............................................................................
  pg = self%pressure ()
  ld = dn(log(d))
  u  = p/exp(ld)
  dddt = -div(p)
  dpdt = -div(d*up(u)**2) - grad(pg)
  dsdt = -div(dn(s)*p)
END SUBROUTINE pde
!> =============================================================================
END MODULE mhd_mod
