!> =============================================================================
!> Stub for self-contained radiative transfer module
!> =============================================================================
MODULE radiation_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE link_mod
  implicit none
  private
  type, public:: radiation_t
    real, dimension(:,:,:), pointer:: rho       ! pointer to density array
    real, dimension(:,:,:), pointer:: temp      ! pointer to temperature array
    real, dimension(:,:,:), pointer:: qq        ! pointer to cooling / heating
    type(mesh_t), dimension(:), pointer:: mesh  ! pointer to mesh
  contains
    procedure:: init
    procedure:: solve
  end type
CONTAINS

!> =============================================================================
!> Do one-time work for each instance of this solver, such as allocating Q, 
!> =============================================================================
SUBROUTINE init (self)
  class(radiation_t):: self
  associate (mx=>self%mesh(1)%gn, my=>self%mesh(2)%gn, mz=>self%mesh(3)%gn)
  !.............................................................................
  call trace%begin ('radiation_mod::init')
  allocate (self%qq(mx,my,mz))
  call trace%end
  end associate
END SUBROUTINE init

!> =============================================================================
!> Solve for the cooling / heating, assuming self has pointers to whatever is
!> relevant (mesh, density, temperature)
!> =============================================================================
SUBROUTINE solve (self)
  class(radiation_t):: self
  !.............................................................................
  call trace%begin ('radiation_mod::solve')
  !-----------------------------------------------------------------------------
  ! Here we could compute cooling via an optically thin cooling function, or
  ! solve the full problem
  !-----------------------------------------------------------------------------
  self%qq = 0.0
  call trace%end
END SUBROUTINE solve

END MODULE radiation_mod
