!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE tensor_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  !USE central_mod
  implicit none
  type, public, extends(vector_t):: tensor_t
    real, dimension(:,:,:), pointer:: xx=>null(),xy=>null(),xz=>null(), &
                                      yx=>null(),yy=>null(),yz=>null(), &
                                      zx=>null(),zy=>null(),zz=>null()
  contains
    procedure:: init => init_tensor
    procedure:: allocate => allocate_tensor
  end type
  interface assignment(=)
    module procedure tensor_tensor_assign
  end interface
  interface operator(+)
    module procedure tensor_tensor_add
  end interface
  interface operator(*)
    module procedure vector_tensor_mul, scalar_tensor_mul
  end interface
  interface operator(-)
    module procedure tensor_tensor_sub
  end interface
  interface grad
    module procedure vector_tensor_grad
  end interface
  interface outers
    module procedure vector_vector_outer
  end interface
  interface outer
    module procedure vector_vector_outerf
  end interface
  interface div
    module procedure vector_tensor_div
  end interface
CONTAINS
!===============================================================================
SUBROUTINE init_tensor (self, dims, ds, name)
  class(tensor_t):: self
  integer:: dims(3)
  real(8):: ds(3)
  character(len=*):: name
  !.............................................................................
  call trace_begin ('tensor_t%init '//trim(self%name), 2)
  self%n = dims
!  self%ds = ds
  self%name = name
  call self%allocate (self)
  call trace_end
END SUBROUTINE init_tensor
!===============================================================================
SUBROUTINE tensor_tensor_assign (out, in)
  type(tensor_t), intent(out):: out
  type(tensor_t), intent(in):: in
  !.............................................................................
  call out%allocate (in)
  out%xx = in%xx
  out%yx = in%yx
  out%zx = in%zx
  out%xy = in%xy
  out%yy = in%yy
  out%zy = in%zy
  out%xz = in%xz
  out%yz = in%yz
  out%zz = in%zz
END SUBROUTINE tensor_tensor_assign
!===============================================================================
FUNCTION tensor_vector_grad (in) RESULT (out)
  type(tensor_t):: out
  type(vector_t), intent(in):: in
  !.............................................................................
  call out%allocate (in)
!  out%xx = out%d%dx (in%x)
!  out%yx = out%d%dy (in%x)
!  out%zx = out%d%dz (in%x)
!  out%xy = out%d%dx (in%y)
!  out%yy = out%d%dy (in%y)
!  out%zy = out%d%dz (in%y)
!  out%xz = out%d%dx (in%z)
!  out%yz = out%d%dy (in%z)
!  out%zz = out%d%dz (in%z)
END FUNCTION tensor_vector_grad
!===============================================================================
FUNCTION vector_tensor_grad (in) RESULT (out)
  type(tensor_t):: out
  type(vector_t), intent(in):: in
  !.............................................................................
  call out%allocate (in)
!  out%xx = out%d%dx (in%x)
!  out%yx = out%d%dy (in%x)
!  out%zx = out%d%dz (in%x)
!  out%xy = out%d%dx (in%y)
!  out%yy = out%d%dy (in%y)
!  out%zy = out%d%dz (in%y)
!  out%xz = out%d%dx (in%z)
!  out%yz = out%d%dy (in%z)
!  out%zz = out%d%dz (in%z)
END FUNCTION vector_tensor_grad
!===============================================================================
SUBROUTINE vector_tensor_div (out, in)
  type(vector_t), intent(inout):: out
  class(tensor_t), intent(in):: in
  !.............................................................................
  call trace_begin ('vector_tensor_div', 3)
  print *,trim(out%name)//' = div('//trim(in%name)
  !call out%allocate (in)
  select type (in)
  type is (tensor_t)
    print*,'v_t_div'
!    print*,'v_t_div', in%d%ds
!    print*,'v_t_div', in%d%n
!    out%x = in%d%dx (in%xx) + in%d%dy(in%yx) + in%d%dz(in%zx)
!    out%y = in%d%dx (in%xy) + in%d%dy(in%yy) + in%d%dz(in%zy)
!    out%z = in%d%dx (in%xz) + in%d%dy(in%yz) + in%d%dz(in%zz)
  class default
    print*,'WARNING'
  end select
  call trace_end
END SUBROUTINE vector_tensor_div
!===============================================================================
FUNCTION tensor_tensor_add (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(tensor_t), intent(in):: in1, in2
  !.............................................................................
  call out%allocate (in1)
  out%xx = in1%xx + in2%xx
  out%yx = in1%yx + in2%yx
  out%zx = in1%zx + in2%zx
  out%xy = in1%xy + in2%xy
  out%yy = in1%yy + in2%yy
  out%zy = in1%zy + in2%zy
  out%xz = in1%xz + in2%xz
  out%yz = in1%yz + in2%yz
  out%zz = in1%zz + in2%zz
END FUNCTION tensor_tensor_add
!===============================================================================
FUNCTION vector_tensor_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(tensor_t), intent(in):: in2
  !.............................................................................
  call out%allocate (in1)
  out%x = in1%x * in2%xx + in1%y * in2%xy + in1%z * in2%xz
  out%y = in1%y * in2%yx + in1%y * in2%yy + in1%z * in2%yz
  out%z = in1%z * in2%zx + in1%y * in2%zy + in1%z * in2%zz
END FUNCTION vector_tensor_mul
!===============================================================================
FUNCTION scalar_tensor_mul (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(scalar_t), intent(in):: in1
  type(tensor_t), intent(in):: in2
  !.............................................................................
  call out%allocate (in1)
  out%xx = in1%f - in2%xx
  out%yx = in1%f - in2%yx
  out%zx = in1%f - in2%zx
  out%xy = in1%f - in2%xy
  out%yy = in1%f - in2%yy
  out%zy = in1%f - in2%zy
  out%xz = in1%f - in2%xz
  out%yz = in1%f - in2%yz
  out%zz = in1%f - in2%zz
END FUNCTION scalar_tensor_mul
!===============================================================================
FUNCTION tensor_tensor_sub (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(tensor_t), intent(in):: in1, in2
  !.............................................................................
  call out%allocate (in1)
  out%xx = in1%xx - in2%xx
  out%yx = in1%yx - in2%yx
  out%zx = in1%zx - in2%zx
  out%xy = in1%xy - in2%xy
  out%yy = in1%yy - in2%yy
  out%zy = in1%zy - in2%zy
  out%xz = in1%xz - in2%xz
  out%yz = in1%yz - in2%yz
  out%zz = in1%zz - in2%zz
END FUNCTION tensor_tensor_sub

!===============================================================================
SUBROUTINE vector_vector_outer (out, in1, in2)
  type(tensor_t), intent(inout):: out
  type(scalar_t), intent(in):: in1
  type(vector_t), intent(in):: in2
  !.............................................................................
  call out%allocate (in1)
  out%xx = in1%f * in2%x * in2%x
  out%yx = in1%f * in2%y * in2%x
  out%zx = in1%f * in2%z * in2%x
  out%xy = in1%f * in2%x * in2%y
  out%yy = in1%f * in2%y * in2%y
  out%zy = in1%f * in2%z * in2%y
  out%xz = in1%f * in2%x * in2%z
  out%yz = in1%f * in2%y * in2%z
  out%zz = in1%f * in2%z * in2%z
END SUBROUTINE vector_vector_outer
!===============================================================================
FUNCTION vector_vector_outerf (in1, in2) RESULT(out)
  type(tensor_t):: out
  type(vector_t), intent(in):: in1, in2
  !.............................................................................
  call out%allocate (in1)
  out%xx = in1%x - in2%x
  out%yx = in1%y - in2%x
  out%zx = in1%z - in2%x
  out%xy = in1%x - in2%y
  out%yy = in1%y - in2%y
  out%zy = in1%z - in2%y
  out%xz = in1%x - in2%z
  out%yz = in1%y - in2%z
  out%zz = in1%z - in2%z
END FUNCTION vector_vector_outerf

!===============================================================================
SUBROUTINE allocate_tensor (self, in)
  class(tensor_t):: self
  class(scalar_t):: in
  !.............................................................................
  if (associated(self%xx)) return
  call trace_begin('allocate_tensor '//trim(self%name), 3)
  if (io%verbose>3) print*,'allocate_tensor',self%n, in%n
  self%n = in%n
  allocate (self%xx (self%n(1), self%n(2), self%n(3)))
  allocate (self%yx (self%n(1), self%n(2), self%n(3)))
  allocate (self%zx (self%n(1), self%n(2), self%n(3)))
  allocate (self%xy (self%n(1), self%n(2), self%n(3)))
  allocate (self%yy (self%n(1), self%n(2), self%n(3)))
  allocate (self%zy (self%n(1), self%n(2), self%n(3)))
  allocate (self%xz (self%n(1), self%n(2), self%n(3)))
  allocate (self%yz (self%n(1), self%n(2), self%n(3)))
  allocate (self%zz (self%n(1), self%n(2), self%n(3)))
  self%allocated = .true.
  call trace_end
END SUBROUTINE allocate_tensor
!===============================================================================
END MODULE tensor_mod
