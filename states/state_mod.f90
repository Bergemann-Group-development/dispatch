!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for states (of variables) of physical systems, e.g., for MHD.
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE state_mod
  USE variable_mod
  implicit none

  type, abstract :: state_t
    ! a state should be composed of at least one scalar or vector variable
  contains
    procedure(initInterface), deferred :: init
  end type
  abstract interface
    subroutine initInterface (self, mem)
      import state_t
      class(state_t) :: self
      real, dimension(:,:,:,:), pointer, intent(in) :: mem
    end subroutine initInterface
  end interface

  type, extends(state_t) :: mhd_state
    real, dimension(:,:,:), pointer :: d, s, phi
    type(vector_t), allocatable :: p, b
  contains
    procedure :: init => init_mhd_state
    procedure :: deallocate => deallocate_mhd_state
  end type

CONTAINS

SUBROUTINE init_mhd_state (self, mem)
  class(mhd_state) :: self
  real, dimension(:,:,:,:), pointer, intent(in) :: mem
  !.............................................................................
  call ScalarFactory (mem, 1, self%d)       ! density
  call ScalarFactory (mem, 2, self%s)       ! entropy
  call ScalarFactory (mem, 3, self%phi)     ! grav. potential
  call VectorFactory (mem, [4,5,6], self%p) ! momentum
  call VectorFactory (mem, [7,8,9], self%b) ! magnetic field
END SUBROUTINE init_mhd_state

SUBROUTINE deallocate_mhd_state (self)
  class(mhd_state) :: self
  !.............................................................................
  call ScalarDestroy (self%d)
  call ScalarDestroy (self%s)
  call ScalarDestroy (self%phi)
  call VectorDestroy (self%p)
  call VectorDestroy (self%b)
END SUBROUTINE deallocate_mhd_state

END MODULE state_mod
