!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE tensor_mod
  USE trace_mod
  USE vector_mod
  USE derivatives_mod
  implicit none
  type:: tensor_t
    real, dimension(:,:,:), pointer:: xx,xy,xz,yx,yy,yz,zx,zy,zz
    logical:: allocated=.false.
    integer:: n(3)
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
  interface outer
    module procedure vector_vector_outer
  end interface
  interface div
    module procedure vector_tensor_div
  end interface
  integer, private:: n(3)
  class(derivatives_t), private, pointer:: d
CONTAINS
!> =============================================================================
SUBROUTINE init_tensor (self, dims, ds)
  class(tensor_t):: self
  integer:: dims(3)
  real:: ds(3)
  !.............................................................................
  call trace_begin ('tensor_t%init', 2)
  n = dims
  self%n = dims
  call trace_end
END SUBROUTINE init_tensor
!> =============================================================================
SUBROUTINE tensor_tensor_assign (out, in)
  type(tensor_t), intent(out):: out
  type(tensor_t), intent(in):: in
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
!> =============================================================================
FUNCTION tensor_vector_grad (in) RESULT (out)
  type(tensor_t):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  n = [size(in%x,dim=1), size(in%x,dim=2), size(in%x,dim=3)]
  out%xx = d%dx (in%x)
  out%yx = d%dy (in%x)
  out%zx = d%dz (in%x)
  out%xy = d%dx (in%y)
  out%yy = d%dy (in%y)
  out%zy = d%dz (in%y)
  out%xz = d%dx (in%z)
  out%yz = d%dy (in%z)
  out%zz = d%dz (in%z)
END FUNCTION tensor_vector_grad
!> =============================================================================
FUNCTION vector_tensor_grad (in) RESULT (out)
  type(tensor_t):: out
  type(vector_t), intent(in):: in
  !-----------------------------------------------------------------------------
  n = [size(in%x,dim=1), size(in%x,dim=2), size(in%x,dim=3)]
  out%xx = d%dx (in%x)
  out%yx = d%dy (in%x)
  out%zx = d%dz (in%x)
  out%xy = d%dx (in%y)
  out%yy = d%dy (in%y)
  out%zy = d%dz (in%y)
  out%xz = d%dx (in%z)
  out%yz = d%dy (in%z)
  out%zz = d%dz (in%z)
END FUNCTION vector_tensor_grad
!> =============================================================================
FUNCTION vector_tensor_div (in) RESULT (out)
  type(vector_t):: out
  type(tensor_t), intent(in):: in
  !-----------------------------------------------------------------------------
  n = [size(in%xx,dim=1), size(in%xx,dim=2), size(in%xx,dim=3)]
  out%x = d%dx (in%xx) + d%dy(in%yx) + d%dz(in%zx)
  out%y = d%dx (in%xy) + d%dy(in%yy) + d%dz(in%zy)
  out%z = d%dx (in%xz) + d%dy(in%yz) + d%dz(in%zz)
END FUNCTION vector_tensor_div
!> =============================================================================
FUNCTION tensor_tensor_add (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(tensor_t), intent(in):: in1, in2
  !-----------------------------------------------------------------------------
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
!> =============================================================================
FUNCTION vector_tensor_mul (in1, in2) RESULT (out)
  type(vector_t):: out
  type(vector_t), intent(in):: in1
  type(tensor_t), intent(in):: in2
  !-----------------------------------------------------------------------------
  out%x = in1%x * in2%xx + in1%y * in2%xy + in1%z * in2%xz
  out%y = in1%y * in2%yx + in1%y * in2%yy + in1%z * in2%yz
  out%z = in1%z * in2%zx + in1%y * in2%zy + in1%z * in2%zz
END FUNCTION vector_tensor_mul
!> =============================================================================
FUNCTION scalar_tensor_mul (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(scalar_t), intent(in):: in1
  type(tensor_t), intent(in):: in2
  !-----------------------------------------------------------------------------
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
!> =============================================================================
FUNCTION tensor_tensor_sub (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(tensor_t), intent(in):: in1, in2
  !-----------------------------------------------------------------------------
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
!> =============================================================================
FUNCTION vector_vector_outer (in1, in2) RESULT (out)
  type(tensor_t):: out
  type(vector_t), intent(in):: in1, in2
  !-----------------------------------------------------------------------------
  out%xx = in1%x - in2%x
  out%yx = in1%y - in2%x
  out%zx = in1%z - in2%x
  out%xy = in1%x - in2%y
  out%yy = in1%y - in2%y
  out%zy = in1%z - in2%y
  out%xz = in1%x - in2%z
  out%yz = in1%y - in2%z
  out%zz = in1%z - in2%z
END FUNCTION vector_vector_outer
!> =============================================================================
SUBROUTINE allocate_tensor (self)
  class(tensor_t):: self
  if (self%allocated) return
  self%n = n
  print*,'allocate_tensor',n
  call flush()
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
END SUBROUTINE allocate_tensor
!> =============================================================================
END MODULE tensor_mod
