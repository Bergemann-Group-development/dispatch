!> =============================================================================
!> $Id$
!> PDE for simple advection
!> =============================================================================
MODULE pde_mod
  USE patch_mod
  USE mesh_mod
  USE io_mod
  implicit none
  type, extends(patch_t):: pde_t
    real:: v(3) = 0.0
  contains
    procedure:: init     => init_advection
    procedure:: allocate => allocate_advection
    procedure:: pde
  end type
PRIVATE
PUBLIC pde_t
CONTAINS

!> =============================================================================
!> Set the number of variables and use the generic allocating routine
!> =============================================================================
SUBROUTINE init_advection (self)
  class(pde_t):: self
  real:: velocity(3)
  namelist /advection_params/ velocity
  !.............................................................................
  call io%init
  rewind (io%input)
  read (io%input,advection_params)
  if (io%master) write(*,advection_params)
  self%v = velocity
  self%nvar = 1
  call self%patch_t%allocate
END SUBROUTINE

!> =============================================================================
!> Set the number of variables and use the generic allocating routine
!> =============================================================================
SUBROUTINE allocate_advection (self)
  class(pde_t):: self
  self%nvar = 1
  call self%patch_t%allocate
END SUBROUTINE

!> =============================================================================
!> The PDE computes the time derivative of the variables.
!> =============================================================================
SUBROUTINE pde (self)
  class(pde_t):: self
  real, dimension(:,:,:), pointer:: rho, drhodt
  !.............................................................................
  rho    => self%mem(:,:,:,1,self%it,1)
  drhodt => self%mem(:,:,:,1,self%it,2)
  drhodt = - self%m(1)%ddx(rho*self%v(1)) &
           - self%m(2)%ddy(rho*self%v(2)) &
           - self%m(3)%ddz(rho*self%v(3))
END SUBROUTINE

END MODULE pde_mod
