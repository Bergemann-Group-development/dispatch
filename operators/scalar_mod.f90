!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE scalar_mod
  USE deriv_mod
  implicit none
  type:: scalar_t
    real, dimension(:,:,:), pointer:: f
    logical:: allocated=.false.
    integer:: n(3)
  contains
    procedure:: allocate => allocate_scalar
    procedure, nopass:: dims => set_dims
  end type
  interface assignment(=)
    module procedure scalar_assign, scalar_assign_real, scalar_assign_int
  end interface
  interface operator(+)
    module procedure scalar_real_add, scalar_int_add, scalar_scalar_add
  end interface
  interface operator(*)
    module procedure scalar_real_mul, scalar_int_mul, scalar_scalar_mul
  end interface
  interface operator(/)
    module procedure scalar_scalar_div
  end interface
  interface operator(-)
    module procedure scalar_real_sub, scalar_int_sub, scalar_scalar_sub
  end interface
  interface exps
    module procedure scalar_exp
  end interface
  interface maxval
    module procedure scalar_maxval
  end interface
  integer:: n(3)=1
CONTAINS
!> =============================================================================
SUBROUTINE scalar_assign_real (out, in)
  type(scalar_t), intent(out):: out
  real, intent(in):: in
  call out%allocate
  out%f = in
  print*,out%n
END SUBROUTINE scalar_assign_real
!> =============================================================================
SUBROUTINE scalar_assign_int (out, in)
  type(scalar_t), intent(out):: out
  integer, intent(in):: in
  call out%allocate
  out%f = in
END SUBROUTINE scalar_assign_int
!> =============================================================================
SUBROUTINE scalar_assign (out, in)
  type(scalar_t), intent(out):: out
  real, dimension(:,:,:), intent(in), pointer:: in
  n = [size(in,dim=1), size(in,dim=2), size(in,dim=3)]
  call out%allocate
  out%f = in
END SUBROUTINE scalar_assign
!> =============================================================================
FUNCTION scalar_maxval (in) RESULT (out)
  real:: out
  type(scalar_t), intent(in):: in
  n = [size(in%f,dim=1), size(in%f,dim=2), size(in%f,dim=3)]
  out = maxval(in%f)
END FUNCTION scalar_maxval
!> =============================================================================
FUNCTION scalar_exp (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  n = [size(in%f,dim=1), size(in%f,dim=2), size(in%f,dim=3)]
  out%f = exp(in%f)
END FUNCTION scalar_exp

!> =============================================================================
FUNCTION scalar_real_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f + in2
  print*,out%n
END FUNCTION scalar_real_add
!> =============================================================================
FUNCTION scalar_int_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f + in2
  print*,out%n
END FUNCTION scalar_int_add
!> =============================================================================
FUNCTION scalar_scalar_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f + in2%f
END FUNCTION scalar_scalar_add

!> =============================================================================
FUNCTION scalar_real_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  real, intent(in):: in1
  type(scalar_t), intent(in):: in2
  n = [size(in2%f,dim=1), size(in2%f,dim=2), size(in2%f,dim=3)]
  call out%allocate
  out%f = in1 * in2%f
  print*,out%n
END FUNCTION scalar_real_mul
!> =============================================================================
FUNCTION scalar_int_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f * in2
  print*,out%n
END FUNCTION scalar_int_mul
!> =============================================================================
FUNCTION scalar_scalar_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f * in2%f
END FUNCTION scalar_scalar_mul
!> =============================================================================
FUNCTION scalar_scalar_div (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f / in2%f
END FUNCTION scalar_scalar_div

!> =============================================================================
FUNCTION scalar_real_sub (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f + in2
  print*,out%n
END FUNCTION scalar_real_sub
!> =============================================================================
FUNCTION scalar_int_sub (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f + in2
  print*,out%n
END FUNCTION scalar_int_sub
!> =============================================================================
FUNCTION scalar_scalar_sub (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  n = [size(in1%f,dim=1), size(in1%f,dim=2), size(in1%f,dim=3)]
  call out%allocate
  out%f = in1%f + in2%f
END FUNCTION scalar_scalar_sub

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
END MODULE scalar_mod
