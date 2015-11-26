!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE tensor_mod
  USE deriv_mod
  USE vector_mod
  implicit none
  type:: tensor_t
    real, dimension(:,:,:), pointer:: xx,xy,xz,yx,yy,yz,zx,zy,zz
    logical:: allocated=.false.
    integer:: n(3)
  contains
    procedure:: allocate => allocate_tensor
  end type
  interface assignment(=)
    module procedure tensor_tensor_assign
  end interface
  interface operator(+)
    module procedure tensor_tensor_add
  end interface
  interface operator(-)
    module procedure tensor_tensor_sub
  end interface
  interface grad
    module procedure tensor_vector_grad
  end interface
  private n
CONTAINS
!> =============================================================================
SUBROUTINE tensor_tensor_assign (out, in)
  type(tensor_t), intent(out):: out
  type(tensor_t), intent(in):: in
  call out%allocate
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
  call deriv%init (n=n)
  call out%allocate
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
SUBROUTINE allocate_tensor (self)
  class(tensor_t):: self
  if (self%allocated) return
  self%n = n
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
