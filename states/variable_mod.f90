!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> Template module for variables and states, e.g., for MHD.
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE variable_mod
  USE mpi_mod
  USE trace_mod
  USE kinds_mod
  implicit none
  private
  
  type, public :: scalar_t
    real(kind=KindScalarVar), dimension(:,:,:), pointer :: q => null()
    real(kind=KindScalarVar), dimension(:,:,:), pointer :: d => null()
    character(len=32) :: name = ''
  contains
    procedure :: init => ScalarFactory
    procedure :: deallocate => ScalarRecycler
  end type
  
  type, public :: vector_t
    type(scalar_t) :: x, y, z  ! components of the vector
    real, dimension(3) :: h = 0.0  ! staggering (if any) of each vector component
    character(len=32) :: name = ''
  contains
    procedure :: init => VectorFactory
    procedure :: deallocate => VectorRecycler
  end type

CONTAINS

!> =============================================================================
!> A modified "factory" design pattern used to create scalar variables.
!> =============================================================================
SUBROUTINE ScalarFactory (self, mem, qidx, tidx)
  class(scalar_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:,:), pointer, intent(in) :: mem
  integer, intent(in) :: qidx, tidx
  !.............................................................................
  if (associated(mem)) then
    self%q => mem(:,:,:,qidx,tidx,1)
    self%d => mem(:,:,:,qidx,tidx,2)
  else
    call mpi%abort('Variable memory has not been allocated! Abort!')
  end if
  
END SUBROUTINE ScalarFactory

!> =============================================================================
!> A means to deallocate a scalar variable
!> =============================================================================
SUBROUTINE ScalarRecycler (self)
  class(scalar_t) :: self
  !.............................................................................
  nullify(self%q)
  nullify(self%d)

END SUBROUTINE ScalarRecycler

!> =============================================================================
!> A modified "factory" design pattern used to create vectors.
!> =============================================================================
SUBROUTINE VectorFactory (self, mem, vidx, tidx, h)
  class(vector_t) :: self
  real(kind=KindScalarVar), dimension(:,:,:,:,:,:), pointer, intent(in) :: mem
  integer, intent(in) :: tidx
  integer, dimension(3), intent(in) :: vidx
  real, dimension(3), optional :: h
  !.............................................................................
  if (present(h)) then
    self%h = h
  endif
  call self%x%init(mem, vidx(1), tidx)
  call self%y%init(mem, vidx(2), tidx)
  call self%z%init(mem, vidx(3), tidx)

END SUBROUTINE VectorFactory

!> =============================================================================
!> A means to deallocate a vector type
!> =============================================================================
SUBROUTINE VectorRecycler (self)
  class(vector_t) :: self
  !.............................................................................
  call self%x%deallocate()
  call self%y%deallocate()
  call self%z%deallocate()

END SUBROUTINE VectorRecycler

END MODULE variable_mod
