!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE scalar_mod
  USE io_mod
  USE trace_mod
  !USE central_mod
  !USE derivatives_mod
  implicit none
  public
  type, public:: scalar_t
    real, dimension(:,:,:), pointer:: f=>null()
    logical:: allocated=.false.
    integer:: n(3)
    integer:: id
    real:: ds(3)
    character(len=32):: name = ''
    !class(central_t), pointer:: d => null()
  contains
    procedure:: init => init_scalar
    procedure:: allocate => allocate_scalar
    procedure, nopass:: dims => set_dims
  end type
  interface assignment(=)
    module procedure scalar_assign, scalar_assign_real, scalar_assign_int
  end interface
  interface operator(+)
    module procedure scalar_real_add, scalar_int_add, scalar_scalar_add
  end interface
  interface operator(**)
    module procedure scalar_power_int, scalar_power_real
  end interface
  interface operator(*)
    module procedure scalar_real_mul, real_scalar_mul, scalar_int_mul, scalar_scalar_mul
  end interface
  interface operator(/)
    module procedure scalar_scalar_divide
  end interface
  interface operator(-)
    module procedure scalar_real_sub, scalar_int_sub, scalar_scalar_sub, scalar_minus
  end interface
  interface exp
    module procedure scalar_exp
  end interface
  interface log
    module procedure scalar_log
  end interface
  interface sqrt
    module procedure scalar_sqrt
  end interface
  interface minval
    module procedure scalar_minval
  end interface
  interface maxval
    module procedure scalar_maxval
  end interface
  private init_scalar, allocate_scalar
CONTAINS
!===============================================================================
SUBROUTINE init_scalar (self, dims, ds, name)
  class(scalar_t):: self
  integer:: dims(3)
  real(8):: ds(3)
  character(len=*):: name
  !.............................................................................
  call trace_begin ('scalar_t%init '//trim(self%name), 2)
  self%n = dims
  self%ds = ds
  self%name = name
  call self%allocate (self)
  call trace_end
END SUBROUTINE init_scalar
!===============================================================================
SUBROUTINE scalar_assign_real (out, in)
  type(scalar_t), intent(inout):: out
  real, intent(in):: in
  if (io%verbose>2) &
  print*,'sc', out%id, out%allocated, associated(out%f), out%name
  if (io%verbose>2) &
  print*,'sc', in, size(out%f)
  call out%allocate (out)
  out%f = in
  if (io%verbose>2) &
  print*,'sc', out%f(1:6,1,1)
END SUBROUTINE scalar_assign_real
!===============================================================================
SUBROUTINE scalar_assign_int (out, in)
  type(scalar_t), intent(inout):: out
  integer, intent(in):: in
  call out%allocate (out)
  out%f = in
END SUBROUTINE scalar_assign_int
!===============================================================================
SUBROUTINE scalar_assign (out, in)
  type(scalar_t), intent(inout):: out
  real, dimension(:,:,:), intent(in), pointer:: in
  call out%allocate (out)
  out%f = in
END SUBROUTINE scalar_assign
!===============================================================================
FUNCTION scalar_minus (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  call out%allocate (in)
  out%f = - in%f
END FUNCTION scalar_minus
!===============================================================================
FUNCTION scalar_maxval (in) RESULT (out)
  real:: out
  type(scalar_t), intent(in):: in
  out = maxval(in%f)
END FUNCTION scalar_maxval
!===============================================================================
FUNCTION scalar_minval (in) RESULT (out)
  real:: out
  type(scalar_t), intent(in):: in
  out = minval(in%f)
END FUNCTION scalar_minval
!===============================================================================
FUNCTION scalar_exp (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  call out%allocate (in)
  out%f = exp(in%f)
END FUNCTION scalar_exp
!===============================================================================
FUNCTION scalar_log (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
!  call out%allocate (in)
  out%f = log(in%f)
END FUNCTION scalar_log
!===============================================================================
FUNCTION scalar_sqrt (in) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in
  call out%allocate (in)
  out%f = sqrt(in%f)
END FUNCTION scalar_sqrt

!===============================================================================
FUNCTION scalar_real_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  call out%allocate (in1)
  out%f = in1%f + in2
  if (io%verbose>3) print*,'scalar_mod:s_r-a', out%n
END FUNCTION scalar_real_add
!===============================================================================
FUNCTION scalar_int_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%allocate (in1)
  out%f = in1%f + in2
  print*,out%n
END FUNCTION scalar_int_add
!===============================================================================
FUNCTION scalar_scalar_add (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call out%allocate (in1)
  out%f = in1%f + in2%f
END FUNCTION scalar_scalar_add

!===============================================================================
FUNCTION real_scalar_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  real, intent(in):: in1
  type(scalar_t), intent(in):: in2
  call out%allocate (in2)
  out%f = in1 * in2%f
  print*,out%n
END FUNCTION real_scalar_mul
!===============================================================================
FUNCTION scalar_real_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  real, intent(in):: in2
  type(scalar_t), intent(in):: in1
  call out%allocate (in1)
  out%f = in2 * in1%f
  print*,out%n
END FUNCTION scalar_real_mul
!===============================================================================
FUNCTION scalar_int_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%allocate (in1)
  out%f = in1%f * in2
  print*,out%n
END FUNCTION scalar_int_mul
!===============================================================================
FUNCTION scalar_scalar_mul (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call out%allocate (in1)
  out%f = in1%f * in2%f
END FUNCTION scalar_scalar_mul
!===============================================================================
FUNCTION scalar_scalar_divide (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call trace_begin ('scalar_scalar_divide '//trim(in1%name)//'/'//trim(in2%name), 3)
  call out%allocate (in2)
  out%f = in1%f / in2%f
  call trace_end
END FUNCTION scalar_scalar_divide

!===============================================================================
FUNCTION scalar_power_int (in1, in2) RESULT(out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%allocate (in1)
  out%f = in1%f**in2
  print*,out%n
END FUNCTION scalar_power_int
!===============================================================================
FUNCTION scalar_power_real (in1, in2) RESULT(out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  if (io%verbose>3) print*,'scalar:', out%id
  call out%allocate (in1)
  out%f = in1%f**in2
  print*,out%n
END FUNCTION scalar_power_real
!===============================================================================
FUNCTION scalar_real_sub (in1, in2) RESULT(out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  real, intent(in):: in2
  call out%allocate (in1)
  out%f = in1%f - in2
  print*,out%n
END FUNCTION scalar_real_sub
!===============================================================================
FUNCTION scalar_int_sub (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1
  integer, intent(in):: in2
  call out%allocate (in1)
  out%f = in1%f - in2
  print*,out%n
END FUNCTION scalar_int_sub
!===============================================================================
FUNCTION scalar_scalar_sub (in1, in2) RESULT (out)
  type(scalar_t):: out
  type(scalar_t), intent(in):: in1, in2
  call out%allocate (in1)
  out%f = in1%f - in2%f
END FUNCTION scalar_scalar_sub

!===============================================================================
SUBROUTINE allocate_scalar (self, in)
  class(scalar_t):: self
  class(scalar_t):: in
  integer:: id=0
  !.............................................................................
  if (associated(self%f)) return
  call trace_begin ('allocate_scalar '//trim(self%name))
  !$omp atomic
  id = id+1
  self%id = id
  if (io%verbose>3) print*,'sc', self%id, in%id, self%allocated, associated(self%f)
  if (io%verbose>3) print*,'allocate_scalar', self%n, in%n
  self%n = in%n
  allocate (self%f(self%n(1), self%n(2), self%n(3)))
  self%allocated = .true.
  call trace_end
END SUBROUTINE allocate_scalar
!===============================================================================
SUBROUTINE set_dims (dims)
  integer:: dims(3)
END SUBROUTINE set_dims
!===============================================================================
END MODULE scalar_mod
