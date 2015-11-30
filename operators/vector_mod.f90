!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE vector_mod
  USE io_mod
  USE trace_mod
  USE scalar_mod
  USE derivatives_mod
  implicit none
  public
  type, public:: vector_t
    real, dimension(:,:,:), pointer:: x, y, z
    logical:: allocated=.false.
    integer:: n(3)
    integer:: id=0
  contains
    procedure:: init
    procedure:: allocate => allocate_vector
  end type
  interface assignment(=)
    module procedure vector_assign_array, vector_assign_real, vector_assign_int, vector_assign_real_3
  end interface
  interface operator(*)
    module procedure real4_vector_mul, real8_vector_mul, vector_vector_mul, vector_scalar_mul
  end interface
  interface dot
    module procedure vector_vector_dot
  end interface
  interface operator(**)
    module procedure vector_power_int
  end interface
  interface operator(/)
    module procedure  vector_scalar_divide, vector_vector_divide
  end interface
  interface operator(+)
    module procedure vector_real_add, vector_int_add, vector_vector_add
  end interface
  interface operator(-)
    module procedure vector_vector_sub, vector_minus
  end interface
  interface grad
    module procedure vector_real_grad, vector_int_grad, vector_scalar_grad
  end interface
  interface div
    module procedure vector_scalar_div
  end interface
  interface cross
    module procedure real4_vector_cross, real8_vector_cross, vector_cross
  end interface
  interface curl
    module procedure vector_curl
  end interface
  interface exp
    module procedure vector_exp
  end interface
  integer, private:: n(3)
  class(derivatives_t), private, pointer:: d
CONTAINS
!> =============================================================================
SUBROUTINE init (self, dims, ds)
  class(vector_t):: self
  integer:: dims(3)
  real:: ds(3)
  !.............................................................................
  call trace_begin ('vector_t%init', 2)

  n = dims
  self%n = dims
  call trace_end
END SUBROUTINE init
!> =============================================================================
FUNCTION vector_power_int (in1, in2) RESULT(out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  integer, intent(in):: in2
  !-----------------------------------------------------------------------------
  out%x = in1%x**in2
  out%y = in1%y**in2
  out%z = in1%z**in2
END FUNCTION vector_power_int
!> =============================================================================
FUNCTION vector_exp (in) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in
  n = [size(in%x,dim=1), size(in%y,dim=2), size(in%z,dim=3)]
  out%x = exp(in%x)
  out%y = exp(in%y)
  out%z = exp(in%z)
END FUNCTION vector_exp
!> =============================================================================
SUBROUTINE vector_assign_array (out, in)
  type(vector_t):: out
  real, dimension(:,:,:,:), intent(in), pointer:: in
  out%x = in(:,:,:,1)
  out%y = in(:,:,:,2)
  out%z = in(:,:,:,3)
END SUBROUTINE vector_assign_array
!> =============================================================================
SUBROUTINE vector_assign_real_3 (out, in)
  type(vector_t), intent(inout):: out
  real, intent(in):: in(3)
  print *,'v_a_e',in
  out%x = in(1)
  out%y = in(2)
  out%z = in(3)
END SUBROUTINE vector_assign_real_3
!> =============================================================================
SUBROUTINE vector_assign_real (out, in)
  type(vector_t):: out
  real, intent(in):: in
  out%x = in
  out%y = in
  out%z = in
END SUBROUTINE vector_assign_real
!> =============================================================================
SUBROUTINE vector_assign_int (out, in)
  type(vector_t):: out
  integer, intent(in):: in
  out%x = in
  out%y = in
  out%z = in
END SUBROUTINE vector_assign_int
!> =============================================================================
FUNCTION vector_real_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  real, intent(in):: in2
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in1%x + in2
  out%y = in1%y + in2
  out%z = in1%z + in2
END FUNCTION vector_real_add
!> =============================================================================
FUNCTION vector_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1,in2
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in2%x * in1%x
  out%y = in2%y * in1%y
  out%z = in2%z * in1%z
END FUNCTION vector_vector_mul
!> =============================================================================
FUNCTION vector_vector_dot (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(vector_t), intent(in):: in1,in2
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%f = in2%x * in1%x &
        + in2%y * in1%y &
        + in2%z * in1%z
END FUNCTION vector_vector_dot
!> =============================================================================
FUNCTION vector_scalar_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(scalar_t), intent(in):: in2
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in1%x * in2%f
  out%y = in1%y * in2%f
  out%y = in1%z * in2%f
END FUNCTION vector_scalar_mul
!> =============================================================================
FUNCTION real4_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in1
  type(vector_t), intent(in):: in2
  n = [size(in2%x,dim=1), size(in2%x,dim=2), size(in2%x,dim=3)]
  out%x = in2%x * in1
  out%y = in2%y * in1
  out%z = in2%z * in1
END FUNCTION real4_vector_mul
!> =============================================================================
FUNCTION real8_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real(8), intent(in):: in1
  type(vector_t), intent(in):: in2
  n = [size(in2%x,dim=1), size(in2%x,dim=2), size(in2%x,dim=3)]
  out%x = in2%x * in1
  out%y = in2%y * in1
  out%z = in2%z * in1
END FUNCTION real8_vector_mul
!> =============================================================================
FUNCTION vector_int_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  integer, intent(in):: in2
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in1%x + in2
  out%y = in1%y + in2
  out%z = in1%z + in2
END FUNCTION vector_int_add
!> =============================================================================
FUNCTION vector_vector_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in1%x + in2%x
  out%y = in1%y + in2%y
  out%z = in1%z + in2%z
END FUNCTION vector_vector_add
!> =============================================================================
FUNCTION vector_vector_sub (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in1%x - in2%x
  out%y = in1%y - in2%y
  out%z = in1%z - in2%z
END FUNCTION vector_vector_sub
!> =============================================================================
FUNCTION vector_square (in1) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in1%x * in1%x
  out%y = in1%y * in1%y
  out%z = in1%z * in1%z
END FUNCTION vector_square
!> =============================================================================
FUNCTION vector_minus (in1) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = -in1%x
  out%y = -in1%y
  out%z = -in1%z
END FUNCTION vector_minus
!> =============================================================================

FUNCTION vector_real_grad (in) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in
  out%x = 0.0
  out%y = 0.0
  out%z = 0.0
  print*,out%n
END FUNCTION vector_real_grad
!> =============================================================================
FUNCTION vector_int_grad (in) RESULT (out)
  type(vector_t):: out
  integer, intent(in):: in
  out%x = 0.0
  out%y = 0.0
  out%z = 0.0
END FUNCTION vector_int_grad
!> =============================================================================
FUNCTION vector_scalar_grad (in) RESULT (out)
  type(vector_t):: out
  type(scalar_t), intent(in):: in
  !-----------------------------------------------------------------------------
  n = [size(in%f,dim=1), size(in%f,dim=2), size(in%f,dim=3)]
  out%x = d%dx (in%f)
  out%y = d%dy (in%f)
  out%z = d%dz (in%f)
END FUNCTION vector_scalar_grad
!> =============================================================================
FUNCTION vector_vector_divide (in1,in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  type(vector_t), intent(in):: in1
  !-----------------------------------------------------------------------------
  n = [size(in2%x,dim=1), size(in2%y,dim=2), size(in2%z,dim=3)]
  out%x = in1%x/in2%x
  out%y = in1%y/in2%y
  out%z = in1%z/in2%z
END FUNCTION vector_vector_divide
!> =============================================================================
FUNCTION vector_scalar_divide (in1,in2) RESULT (out)
  type(vector_t):: out
  type(scalar_t), intent(in):: in2
  type(vector_t), intent(in):: in1
  !-----------------------------------------------------------------------------
  n = [size(in2%f,dim=1), size(in2%f,dim=2), size(in2%f,dim=3)]
  out%x = in1%x/in2%f
  out%y = in1%y/in2%f
  out%z = in1%z/in2%f
END FUNCTION vector_scalar_divide
!> =============================================================================
FUNCTION vector_scalar_div (in) RESULT (out)
  type(scalar_t):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  n = [size(in%x,dim=1), size(in%x,dim=2), size(in%x,dim=3)]
  out%f = d%dx (in%x) + d%dy (in%y) + d%dz (in%z)
END FUNCTION vector_scalar_div
!> =============================================================================
FUNCTION real8_vector_cross (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  real(8):: in1(3)
  !-----------------------------------------------------------------------------
  n = [size(in2%x,dim=1), size(in2%x,dim=2), size(in2%x,dim=3)]
  out%x = in1(2)*in2%z - in1(3)*in2%y
  out%y = in1(3)*in2%x - in1(1)*in2%z
  out%z = in1(1)*in2%y - in1(2)*in2%x
  print*,'xx'
END FUNCTION real8_vector_cross
!> =============================================================================
FUNCTION real4_vector_cross (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  real(4):: in1(3)
  !-----------------------------------------------------------------------------
  n = [size(in2%x,dim=1), size(in2%x,dim=2), size(in2%x,dim=3)]
  out%x = in1(2)*in2%z - in1(3)*in2%y
  out%y = in1(3)*in2%x - in1(1)*in2%z
  out%z = in1(1)*in2%y - in1(2)*in2%x
END FUNCTION real4_vector_cross
!> =============================================================================
FUNCTION vector_cross (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  !-----------------------------------------------------------------------------
  n = [size(in1%x,dim=1), size(in1%x,dim=2), size(in1%x,dim=3)]
  out%x = in1%y*in2%z - in1%z*in2%y
  out%y = in1%z*in2%x - in1%x*in2%z
  out%z = in1%x*in2%y - in1%y*in2%x
END FUNCTION vector_cross
!> =============================================================================
FUNCTION vector_curl (in) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  n = [size(in%x,dim=1), size(in%x,dim=2), size(in%x,dim=3)]
  out%x = d%dy(in%z) - d%dz(in%y)
  out%y = d%dz(in%x) - d%dx(in%z)
  out%z = d%dx(in%y) - d%dy(in%x)
END FUNCTION vector_curl
!> =============================================================================
SUBROUTINE allocate_vector (self)
  class(vector_t):: self
  if (self%allocated) return
  if (io%verbose>3) print*,'allocate_vector',n,self%n
  call flush()
  self%n = n
  allocate (self%x(self%n(1), self%n(2), self%n(3)))
  allocate (self%y(self%n(1), self%n(2), self%n(3)))
  allocate (self%z(self%n(1), self%n(2), self%n(3)))
  self%allocated = .true.
END SUBROUTINE allocate_vector
!> =============================================================================
END MODULE vector_mod
