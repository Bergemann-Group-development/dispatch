!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE vector_mod
  USE scalar_mod
  implicit none
  public
  type:: vector_t
    type(scalar_t):: x, y, z
    logical:: allocated=.false.
    integer:: n(3)
    integer:: id
    real:: ds(3)
    character(len=32):: name = ''
  contains
    procedure:: init => init_vector
    procedure:: allocate => allocate_vector
    procedure:: tempdealloc => tempdealloc_vector
  end type
  interface assignment(=)
    module procedure vector_assign_array, vector_assign_real, vector_assign_int, &
                     vector_assign_vector, vector_assign_scalar, vector_assign_real_3
  end interface
  interface operator(*)
    module procedure real4_vector_mul, real8_vector_mul, vector_vector_mul, vector_scalar_mul, &
                     scalar_vector_mul, const3_scalar_mul, const3_vector_mul
  end interface
  interface operator(**)
    module procedure vector_power_int
  end interface
  interface operator(/)
    module procedure  vector_scalar_divide, vector_vector_divide
  end interface
  interface operator(+)
    module procedure vector_real_add, vector_int_add, vector_vector_add, vector_scalar_add
  end interface
  interface operator(-)
    module procedure vector_vector_sub, vector_scalar_sub, vector_minus
  end interface
  interface dot
    module procedure vector_vector_dot
  end interface
  interface exp
    module procedure vector_exp
  end interface
  interface log
    module procedure vector_log
  end interface
  private:: allocate_vector
CONTAINS

!> =============================================================================
SUBROUTINE init_vector (self, dims, ds, name)
  class(vector_t):: self
  type(scalar_t):: a, b
  integer:: dims(3)
  real:: ds(3)
  integer, save:: id=0
  character(len=*), optional:: name
  !.............................................................................
  !$omp atomic
  id = id+1
  self%id = id
  self%n = dims
  self%ds = ds
  self%name=''; if (present(name)) self%name = name
  !print*,'init',self%n,self%ds
  call self%x%init (dims, ds, 'x')
  call self%y%init (dims, ds, 'y')
  call self%z%init (dims, ds, 'z')
  self%allocated = .true.
END SUBROUTINE init_vector
!> =============================================================================
SUBROUTINE init_vectors (dims, ds, v1, v2, v3, v4, v5, v6, v7, v8, v9)
  integer:: dims(3)
  real:: ds(3)
  type(vector_t), optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9
  !.............................................................................
  if (present(v1)) call init_vector (v1, dims, ds)
  if (present(v2)) call init_vector (v2, dims, ds)
  if (present(v3)) call init_vector (v3, dims, ds)
  if (present(v4)) call init_vector (v4, dims, ds)
  if (present(v5)) call init_vector (v5, dims, ds)
  if (present(v6)) call init_vector (v6, dims, ds)
  if (present(v7)) call init_vector (v7, dims, ds)
  if (present(v8)) call init_vector (v8, dims, ds)
  if (present(v9)) call init_vector (v9, dims, ds)
END SUBROUTINE init_vectors
!> =============================================================================
SUBROUTINE allocate_vector (self)
  class(vector_t):: self
  !.............................................................................
  if (self%allocated) return
  self%allocated = .true.
END SUBROUTINE allocate_vector
!> =============================================================================
SUBROUTINE tempdealloc_vector (self)
  class(vector_t):: self
  call self%x%tempdealloc
  call self%y%tempdealloc
  call self%z%tempdealloc
END SUBROUTINE tempdealloc_vector

!> =============================================================================
SUBROUTINE minus_div (out, in)
  type(vector_t), intent(in):: in
  type(scalar_t), intent(inout):: out
  !.............................................................................
  call out%tempalloc (in%x)
  !out = - ddxup(in%x) - ddyup(in%y) - ddzup(in%z)
  call in%tempdealloc
END SUBROUTINE

!> =============================================================================
FUNCTION vector_power_int (in1, in2) RESULT(out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  integer, intent(in):: in2
  !-----------------------------------------------------------------------------
  out%x = in1%x**in2
  out%y = in1%y**in2
  out%z = in1%z**in2
  call in1%tempdealloc
END FUNCTION vector_power_int
!> =============================================================================
FUNCTION vector_exp (in) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in
  out%x = exp(in%x)
  out%y = exp(in%y)
  out%z = exp(in%z)
  call in%tempdealloc
END FUNCTION vector_exp
!> =============================================================================
FUNCTION vector_log (in) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in
  out%x = log(in%x)
  out%y = log(in%y)
  out%z = log(in%z)
  call in%tempdealloc
END FUNCTION vector_log

!> =============================================================================
SUBROUTINE vector_assign_vector (out, in)
  type(vector_t), intent(inout):: out
  type(vector_t), intent(in):: in
  out%x = in%x
  out%y = in%y
  out%z = in%z
  call in%tempdealloc
END SUBROUTINE vector_assign_vector
!> =============================================================================
SUBROUTINE vector_assign_array (out, in)
  type(vector_t), intent(inout):: out
  real, dimension(:,:,:,:), intent(in), pointer:: in
  out%x = in(:,:,:,1)
  out%y = in(:,:,:,2)
  out%z = in(:,:,:,3)
END SUBROUTINE vector_assign_array
!> =============================================================================
SUBROUTINE vector_assign_scalar (out, in)
  type(vector_t), intent(inout):: out
  type(scalar_t):: in
  out%x = in
  out%y = in
  out%z = in
  call in%tempdealloc
END SUBROUTINE vector_assign_scalar
!> =============================================================================
SUBROUTINE vector_assign_real_3 (out, in)
  type(vector_t), intent(inout):: out
  real, intent(in):: in(3)
  out%x = in(1)
  out%y = in(2)
  out%z = in(3)
END SUBROUTINE vector_assign_real_3
!> =============================================================================
SUBROUTINE vector_assign_real (out, in)
  type(vector_t), intent(inout):: out
  real, intent(in):: in
  out%x = in
  out%y = in
  out%z = in
END SUBROUTINE vector_assign_real
!> =============================================================================
SUBROUTINE vector_assign_int (out, in)
  type(vector_t), intent(inout):: out
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
  out%x = in1%x + in2
  out%y = in1%y + in2
  out%z = in1%z + in2
  call in1%tempdealloc
END FUNCTION vector_real_add
!> =============================================================================
FUNCTION vector_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1,in2
  out%x = in2%x * in1%x
  out%y = in2%y * in1%y
  out%z = in2%z * in1%z
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_vector_mul
!> =============================================================================
FUNCTION vector_vector_dot (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(vector_t), intent(in):: in1,in2
  out = in2%x * in1%x &
      + in2%y * in1%y &
      + in2%z * in1%z
END FUNCTION vector_vector_dot
!> =============================================================================
FUNCTION vector_scalar_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(scalar_t), intent(in):: in2
  out%x = in1%x * in2
  out%y = in1%y * in2
  out%y = in1%z * in2
  call in1%tempdealloc
  call in2%tempdealloc
  !out%d => in2%d
END FUNCTION vector_scalar_mul
!> =============================================================================
FUNCTION scalar_vector_mul (in2, in1) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(scalar_t), intent(in):: in2
  out%x = in1%x * in2
  out%y = in1%y * in2
  out%y = in1%z * in2
  call in1%tempdealloc
  call in2%tempdealloc
  !out%d => in2%d
END FUNCTION scalar_vector_mul
!> =============================================================================
FUNCTION const3_scalar_mul (in2, in1) RESULT (out)
  type(scalar_t):: out
  real, intent(in):: in1(3)
  type(vector_t), intent(in):: in2
  out = in2%x * in1(1)
  out = in2%y * in1(2)
  out = in2%z * in1(3)
  call in2%tempdealloc
END FUNCTION const3_scalar_mul
!> =============================================================================
FUNCTION const3_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in1(3)
  type(vector_t), intent(in):: in2
  out%x = in2%x * in1(1)
  out%y = in2%y * in1(2)
  out%z = in2%z * in1(3)
  call in2%tempdealloc
END FUNCTION const3_vector_mul
!> =============================================================================
FUNCTION real4_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real, intent(in):: in1
  type(vector_t), intent(in):: in2
  out%x = in2%x * in1
  out%y = in2%y * in1
  out%z = in2%z * in1
  call in2%tempdealloc
END FUNCTION real4_vector_mul
!> =============================================================================
FUNCTION real8_vector_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  real(8), intent(in):: in1
  type(vector_t), intent(in):: in2
  out%x = in1 * in2%x
  out%y = in1 * in2%y
  out%z = in1 * in2%z
  call in2%tempdealloc
END FUNCTION real8_vector_mul
!> =============================================================================
FUNCTION vector_int_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  integer, intent(in):: in2
  out%x = in1%x + in2
  out%y = in1%y + in2
  out%z = in1%z + in2
  call in1%tempdealloc
END FUNCTION vector_int_add
!> =============================================================================
FUNCTION vector_vector_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  out%x = in1%x + in2%x
  out%y = in1%y + in2%y
  out%z = in1%z + in2%z
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_vector_add
!> =============================================================================
FUNCTION vector_scalar_add (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(scalar_t), intent(in):: in2
  out%x = in1%x + in2
  out%y = in1%y + in2
  out%z = in1%z + in2
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_scalar_add
!> =============================================================================
FUNCTION vector_vector_sub (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  out%x = in1%x - in2%x
  out%y = in1%y - in2%y
  out%z = in1%z - in2%z
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_vector_sub
!> =============================================================================
FUNCTION vector_scalar_sub (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(scalar_t), intent(in):: in2
  out%x = in1%x - in2
  out%y = in1%y - in2
  out%z = in1%z - in2
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_scalar_sub
!> =============================================================================
FUNCTION vector_square (in1) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  out%x = in1%x * in1%x
  out%y = in1%y * in1%y
  out%z = in1%z * in1%z
  call in1%tempdealloc
END FUNCTION vector_square
!> =============================================================================
FUNCTION vector_minus (in1) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  out%x = -in1%x
  out%y = -in1%y
  out%z = -in1%z
  call in1%tempdealloc
END FUNCTION vector_minus

!> =============================================================================
FUNCTION vector_vector_divide (in1,in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  type(vector_t), intent(in):: in1
  !-----------------------------------------------------------------------------
  out%x = in1%x/in2%x
  out%y = in1%y/in2%y
  out%z = in1%z/in2%z
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_vector_divide
!> =============================================================================
FUNCTION vector_scalar_divide (in1,in2) RESULT (out)
  type(vector_t):: out
  type(scalar_t), intent(in):: in2
  type(vector_t), intent(in):: in1
  !-----------------------------------------------------------------------------
  out%x = in1%x/in2
  out%y = in1%y/in2
  out%z = in1%z/in2
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_scalar_divide

!> =============================================================================
FUNCTION real8_vector_cross (in1, in2) result (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  real(8):: in1(3)
  !-----------------------------------------------------------------------------
  out%x = in1(2)*in2%z - in1(3)*in2%y
  out%y = in1(3)*in2%x - in1(1)*in2%z
  out%z = in1(1)*in2%y - in1(2)*in2%x
  call in2%tempdealloc
END FUNCTION real8_vector_cross
!> =============================================================================
FUNCTION real4_vector_cross (in1, in2) result (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in2
  real(4):: in1(3)
  !-----------------------------------------------------------------------------
  out%x = in1(2)*in2%z - in1(3)*in2%y
  out%y = in1(3)*in2%x - in1(1)*in2%z
  out%z = in1(1)*in2%y - in1(2)*in2%x
  call in2%tempdealloc
END FUNCTION real4_vector_cross
!> =============================================================================
FUNCTION vector_cross (in1, in2) result (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1, in2
  !-----------------------------------------------------------------------------
  out%x = in1%y*in2%z - in1%z*in2%y
  out%y = in1%z*in2%x - in1%x*in2%z
  out%z = in1%x*in2%y - in1%y*in2%x
  call in1%tempdealloc
  call in2%tempdealloc
END FUNCTION vector_cross
!> =============================================================================
FUNCTION real4_array_cross (f, g) result (h)
  real, dimension(:,:,:,:), intent(in):: f, g
  real, dimension(size(f,1),size(f,2),size(f,3),size(f,4)):: h
  !-----------------------------------------------------------------------------
  h(:,:,:,1) = f(:,:,:,2)*g(:,:,:,3)-f(:,:,:,3)*g(:,:,:,2)
  h(:,:,:,2) = f(:,:,:,3)*g(:,:,:,1)-f(:,:,:,1)*g(:,:,:,3)
  h(:,:,:,3) = f(:,:,:,1)*g(:,:,:,2)-f(:,:,:,2)*g(:,:,:,1)
END FUNCTION real4_array_cross

!> =============================================================================
END MODULE vector_mod
