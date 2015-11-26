!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE operators_mod
  type:: scalar_t
    real, dimension(:,:,:), pointer:: f
    logical:: allocated=.false.
    integer:: n(3)
  contains
    procedure:: allocate => allocate_scalar
    procedure, nopass:: dims => set_dims
  end type
  interface assignment(=)
    module procedure scalar_assign
  end interface
  interface operator(+)
    module procedure scalar_real_add, scalar_int_add, scalar_scalar_add
  end interface
  integer, private:: n(3)=1
CONTAINS
!> =============================================================================
SUBROUTINE scalar_assign (out, in)
  type(scalar_t), intent(out):: out
  real, dimension(:,:,:), intent(in), pointer:: in
  call out%allocate
  out%f = in
  print*,out%n
END SUBROUTINE scalar_assign
!> =============================================================================
FUNCTION scalar_real_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  call out%allocate
  out%f = in1%f + in2
  print*,out%n
END FUNCTION scalar_real_add
!> =============================================================================
FUNCTION scalar_int_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%allocate
  out%f = in1%f + in2
  print*,out%n
END FUNCTION scalar_int_add
!> =============================================================================
FUNCTION scalar_scalar_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call out%allocate
  out%f = in1%f + in2%f
END FUNCTION scalar_scalar_add
!> =============================================================================
SUBROUTINE allocate_scalar (self)
  class(scalar_t):: self
  if (self%allocated) return
  self%n = n
  allocate (self%f(self%n(1), self%n(2), self%n(3)))
  self%allocated = .true.
END SUBROUTINE allocate_scalar
!> =============================================================================
SUBROUTINE set_dims (dims)
  integer:: dims(3)
  n = dims
  print*,n
END SUBROUTINE set_dims
!> =============================================================================
END MODULE operators_mod
