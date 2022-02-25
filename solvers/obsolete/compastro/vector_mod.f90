!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE vector_mod
  USE io_mod
  USE trace_mod
  USE scalar_mod
  !USE central_mod
  implicit none
  public
  type, public, extends(scalar_t):: vector_t
    real, dimension(:,:,:), pointer:: x=>null(), y=>null(), z=>null()
  contains
    procedure:: init => init_vector
    procedure:: allocate => allocate_vector
  end type
  interface assignment(=)
    module procedure vector_assign_array, vector_assign_real, vector_assign_int, &
                     vector_assign_scalar, vector_assign_real_3
  end interface
  interface operator(*)
    module procedure real4_vector_mul, real8_vector_mul, vector_vector_mul, vector_scalar_mul, &
                     const3_scalar_mul, const3_vector_mul
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
  interface minusdiv
    module procedure vector_scalar_minusdiv
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
  private:: allocate_vector
CONTAINS

!===============================================================================
SUBROUTINE minus_div (out, in)
  type(vector_t), intent(in):: in
  type(scalar_t), intent(inout):: out
  !.............................................................................
  !print*,in%d%n, in%d%ds, associated(out%f)
  !out%f = - in%d%dx(in%x) - in%d%dy(in%y) - in%d%dz(in%z)
END SUBROUTINE

!===============================================================================
SUBROUTINE init_vector (self, dims, ds, name)
  class(vector_t):: self
  integer:: dims(3)
  real(8):: ds(3)
  character(len=*):: name
  !.............................................................................
  call trace_begin ('vector_t%init '//trim(self%name), 2)
  self%n = dims
  self%ds = ds
  self%name = name
  call self%allocate (self)
  call trace_end
END SUBROUTINE init_vector
!===============================================================================
FUNCTION vector_power_int (in1, in2) RESULT(out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  integer, intent(in):: in2
  !-----------------------------------------------------------------------------
  call out%allocate (in1)
  out%x = in1%x**in2
  out%y = in1%y**in2
  out%z = in1%z**in2
END FUNCTION vector_power_int
!===============================================================================
FUNCTION vector_exp (in) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in
  call out%allocate (in)
  out%x = exp(in%x)
  out%y = exp(in%y)
  out%z = exp(in%z)
END FUNCTION vector_exp
!===============================================================================
SUBROUTINE vector_assign_array (out, in)
  type(vector_t), intent(inout):: out
  real, dimension(:,:,:,:), intent(in), pointer:: in
  call out%allocate (out)
  out%x = in(:,:,:,1)
  out%y = in(:,:,:,2)
  out%z = in(:,:,:,3)
END SUBROUTINE vector_assign_array
!===============================================================================
SUBROUTINE vector_assign_scalar (out, in)
  type(vector_t), intent(inout):: out
  type(scalar_t):: in
  call out%allocate (out)
  out%x = in%f
  out%y = in%f
  out%z = in%f
END SUBROUTINE vector_assign_scalar
!===============================================================================
SUBROUTINE vector_assign_real_3 (out, in)
  type(vector_t), intent(inout):: out
  real, intent(in):: in(3)
  if (io%verbose>3) print *,'v_a_e',in
  call out%allocate (out)
  out%x = in(1)
  out%y = in(2)
  out%z = in(3)
END SUBROUTINE vector_assign_real_3
!===============================================================================
SUBROUTINE vector_assign_real (out, in)
  type(vector_t), intent(inout):: out
  real, intent(in):: in
  call out%allocate (out)
  out%x = in
  out%y = in
  out%z = in
END SUBROUTINE vector_assign_real
!===============================================================================
SUBROUTINE vector_assign_int (out, in)
  type(vector_t), intent(inout):: out
  integer, intent(in):: in
  call out%allocate (out)
  out%x = in
  out%y = in
  out%z = in
END SUBROUTINE vector_assign_int
!===============================================================================
FUNCTION vector_real_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  real, intent(in):: in2
  call out%allocate (in1)
  out%x = in1%x + in2
  out%y = in1%y + in2
  out%z = in1%z + in2
END FUNCTION vector_real_add
!===============================================================================
FUNCTION vector_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1,in2
  call out%allocate (in1)
  out%x = in2%x * in1%x
  out%y = in2%y * in1%y
  out%z = in2%z * in1%z
END FUNCTION vector_vector_mul
!===============================================================================
FUNCTION vector_vector_dot (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(vector_t), intent(in):: in1,in2
  call out%allocate (in1)
  out%f = in2%x * in1%x &
        + in2%y * in1%y &
        + in2%z * in1%z
END FUNCTION vector_vector_dot
!===============================================================================
FUNCTION vector_scalar_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(scalar_t), intent(in):: in2
  call out%allocate (in1)
  out%x = in1%x * in2%f
  out%y = in1%y * in2%f
  out%y = in1%z * in2%f
  !out%d => in2%d
END FUNCTION vector_scalar_mul
!===============================================================================
FUNCTION const3_scalar_mul (in2, in1) RESULT (out)
  type(scalar_t):: out
  real, intent(in):: in1(3)
  type(vector_t), intent(in):: in2
!  call out%allocate (in2)
  out%f = in2%x * in1(1)
  out%f = in2%y * in1(2)
  out%f = in2%z * in1(3)
END FUNCTION const3_scalar_mul
!===============================================================================
FUNCTION const3_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in1(3)
  type(vector_t), intent(in):: in2
!  call out%allocate (in2)
  out%x = in2%x * in1(1)
  out%y = in2%y * in1(2)
  out%z = in2%z * in1(3)
END FUNCTION const3_vector_mul
!===============================================================================
FUNCTION real4_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in1
  type(vector_t), intent(in):: in2
  call out%allocate (in2)
  out%x = in2%x * in1
  out%y = in2%y * in1
  out%z = in2%z * in1
END FUNCTION real4_vector_mul
!===============================================================================
FUNCTION real8_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real(8), intent(in):: in1
  type(vector_t), intent(in):: in2
  call out%allocate (in2)
  out%x = in2%x * in1
  out%y = in2%y * in1
  out%z = in2%z * in1
END FUNCTION real8_vector_mul
!===============================================================================
FUNCTION vector_int_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  integer, intent(in):: in2
  call out%allocate (in1)
  out%x = in1%x + in2
  out%y = in1%y + in2
  out%z = in1%z + in2
END FUNCTION vector_int_add
!===============================================================================
FUNCTION vector_vector_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  call out%allocate (in1)
  out%x = in1%x + in2%x
  out%y = in1%y + in2%y
  out%z = in1%z + in2%z
END FUNCTION vector_vector_add
!===============================================================================
FUNCTION vector_vector_sub (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  call out%allocate (in1)
  out%x = in1%x - in2%x
  out%y = in1%y - in2%y
  out%z = in1%z - in2%z
END FUNCTION vector_vector_sub
!===============================================================================
FUNCTION vector_square (in1) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  call out%allocate (in1)
  out%x = in1%x * in1%x
  out%y = in1%y * in1%y
  out%z = in1%z * in1%z
END FUNCTION vector_square
!===============================================================================
FUNCTION vector_minus (in1) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  call out%allocate (in1)
  out%x = -in1%x
  out%y = -in1%y
  out%z = -in1%z
END FUNCTION vector_minus
!===============================================================================

FUNCTION vector_real_grad (in) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in
  call out%allocate (out)
  out%x = 0.0
  out%y = 0.0
  out%z = 0.0
  if (io%verbose>3) print*,'v_r_g',out%n
END FUNCTION vector_real_grad
!===============================================================================
FUNCTION vector_int_grad (in) RESULT (out)
  type(vector_t):: out
  integer, intent(in):: in
  call out%allocate (out)
  out%x = 0.0
  out%y = 0.0
  out%z = 0.0
END FUNCTION vector_int_grad
SUBROUTINE check_dims (f, label)
  real, dimension(:,:,:), intent(inout), pointer:: f
  character(len=*):: label
  if (io%verbose>2) &
  print *, 'check_dims: '//label, size(f,dim=1), size(f,dim=2), size(f,dim=3)
END
!===============================================================================
FUNCTION vector_scalar_grad (in) RESULT (out)
  type(vector_t):: out
  type(scalar_t), intent(inout):: in
  !-----------------------------------------------------------------------------
  call trace_begin ('vector_scalar_grad')
  call check_dims(in%f, 'in%f')
  print*, associated(out%x), out%allocated
  call check_dims(out%x, 'out%x')
  call check_dims(out%y, 'out%y')
  call check_dims(out%z, 'out%z')
  call out%allocate (in)
  print*, associated(out%x), out%allocated
  !out%x = out%d%dx (in%f)
  !out%y = out%d%dy (in%f)
  !out%z = out%d%dz (in%f)
  call trace_end
END FUNCTION vector_scalar_grad
!===============================================================================
FUNCTION grads (in) RESULT(out)
  type(vector_t):: out
  type(scalar_t), intent(in):: in
  !-----------------------------------------------------------------------------
  call trace_begin ('grad')
  print *, storage_size(out%x)
  print *, storage_size(out%y)
  print *, storage_size(out%z)
  call out%allocate (in)
  !out%x = out%d%dx (in%f)
  !out%y = out%d%dy (in%f)
  !out%z = out%d%dz (in%f)
  call trace_end
END FUNCTION grads
!===============================================================================
FUNCTION vector_vector_divide (in1,in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  type(vector_t), intent(in):: in1
  !-----------------------------------------------------------------------------
  call out%allocate (in1)
  out%x = in1%x/in2%x
  out%y = in1%y/in2%y
  out%z = in1%z/in2%z
END FUNCTION vector_vector_divide
!===============================================================================
FUNCTION vector_scalar_divide (in1,in2) RESULT (out)
  type(vector_t):: out
  type(scalar_t), intent(in):: in2
  type(vector_t), intent(in):: in1
  !-----------------------------------------------------------------------------
  call trace_begin ('vector_scalar_divide '//trim(in1%name)//'/'//trim(in2%name), 3)
  call out%allocate (in1)
  out%x = in1%x/in2%f
  out%y = in1%y/in2%f
  out%z = in1%z/in2%f
  call trace_end
END FUNCTION vector_scalar_divide
!===============================================================================
SUBROUTINE vector_scalar_div (out, in)
  type(scalar_t), intent(inout):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  call trace_begin ('vector_scalar_div', 3)
  !  out%f = in%d%dx (in%x) + in%d%dy (in%y) + in%d%dz (in%z)
!  out%d => in%d
  call trace_end
END SUBROUTINE vector_scalar_div
!===============================================================================
SUBROUTINE vector_scalar_minusdiv (out, in)
  type(scalar_t):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  call trace_begin ('vector_scalar_minus_div', 3)
  !out%f = -in%d%dx (in%x) - in%d%dy (in%y) - in%d%dz (in%z)
!  out%d => in%d
  call trace_end
END SUBROUTINE vector_scalar_minusdiv

!===============================================================================
FUNCTION real8_vector_cross (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  real(8):: in1(3)
  !-----------------------------------------------------------------------------
  call out%allocate (in2)
  out%x = in1(2)*in2%z - in1(3)*in2%y
  out%y = in1(3)*in2%x - in1(1)*in2%z
  out%z = in1(1)*in2%y - in1(2)*in2%x
END FUNCTION real8_vector_cross
!===============================================================================
FUNCTION real4_vector_cross (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  real(4):: in1(3)
  !-----------------------------------------------------------------------------
  call out%allocate (in2)
  out%x = in1(2)*in2%z - in1(3)*in2%y
  out%y = in1(3)*in2%x - in1(1)*in2%z
  out%z = in1(1)*in2%y - in1(2)*in2%x
END FUNCTION real4_vector_cross
!===============================================================================
FUNCTION vector_cross (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  !-----------------------------------------------------------------------------
  call out%allocate (in1)
  out%x = in1%y*in2%z - in1%z*in2%y
  out%y = in1%z*in2%x - in1%x*in2%z
  out%z = in1%x*in2%y - in1%y*in2%x
END FUNCTION vector_cross
!===============================================================================
FUNCTION vector_curl (in) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  call out%allocate (in)
  !out%x = out%d%dy(in%z) - out%d%dz(in%y)
  !out%y = out%d%dz(in%x) - out%d%dx(in%z)
  !out%z = out%d%dx(in%y) - out%d%dy(in%x)
END FUNCTION vector_curl
!===============================================================================
SUBROUTINE allocate_vector (self, in)
  class(vector_t):: self
  class(scalar_t):: in
  integer:: id=0
  !.............................................................................
  if (associated(self%x)) return
  call trace_begin ('allocate vector '//trim(self%name))
  !$omp atomic
  id = id+1
  self%id = id
  if (io%verbose>3) print*,'allocate_vector',self%n, in%n
  self%n = in%n
  select type (self)
  type is (vector_t)
    allocate (self%x(self%n(1), self%n(2), self%n(3)))
    allocate (self%y(self%n(1), self%n(2), self%n(3)))
    allocate (self%z(self%n(1), self%n(2), self%n(3)))
  end select
  self%allocated = .true.
  call trace_end
END SUBROUTINE allocate_vector
!===============================================================================
END MODULE vector_mod
