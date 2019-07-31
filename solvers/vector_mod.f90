!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE vector_mod
  USE io_mod
  USE trace_mod
  implicit none
  public
  integer, private:: verbose=0
  interface assignment(=)
    module procedure vector_assign_scalar
  end interface
  interface dot
    module procedure vector_vector_dot
  end interface
  interface operator(*)
    module procedure vector_scalar_mul, scalar_vector_mul, const3_vector_mul
  end interface
  interface operator(+)
    module procedure vector_scalar_add
  end interface
  interface operator(-)
    module procedure vector_scalar_sub
  end interface
  interface operator(/)
    module procedure vector_scalar_divide
  end interface
  interface cross
    module procedure vector_cross_vector, vector_cross_real4, real4_cross_vector
  end interface
CONTAINS

!===============================================================================
SUBROUTINE allocate_vectors_a (n, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:,:), allocatable, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('allocate_vectors_a')
  if (present(v1)) allocate(v1(n(1),n(2),n(3),3))
  if (present(v2)) allocate(v2(n(1),n(2),n(3),3))
  if (present(v3)) allocate(v3(n(1),n(2),n(3),3))
  if (present(v4)) allocate(v4(n(1),n(2),n(3),3))
  if (present(v5)) allocate(v5(n(1),n(2),n(3),3))
  if (present(v6)) allocate(v6(n(1),n(2),n(3),3))
  if (present(v7)) allocate(v7(n(1),n(2),n(3),3))
  if (present(v8)) allocate(v8(n(1),n(2),n(3),3))
  if (present(v9)) allocate(v9(n(1),n(2),n(3),3))
  if (present(v10)) allocate(v10(n(1),n(2),n(3),3))
  if (present(v1)) call io%bits_mem(storage_size(v1),product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),product(shape(v10)),'v10')
  call trace_end
END SUBROUTINE allocate_vectors_a

!===============================================================================
SUBROUTINE allocate_vectors (n, v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:,:), pointer, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('allocate_vectors')
  if (present(v1)) allocate(v1(n(1),n(2),n(3),3))
  if (present(v2)) allocate(v2(n(1),n(2),n(3),3))
  if (present(v3)) allocate(v3(n(1),n(2),n(3),3))
  if (present(v4)) allocate(v4(n(1),n(2),n(3),3))
  if (present(v5)) allocate(v5(n(1),n(2),n(3),3))
  if (present(v6)) allocate(v6(n(1),n(2),n(3),3))
  if (present(v7)) allocate(v7(n(1),n(2),n(3),3))
  if (present(v8)) allocate(v8(n(1),n(2),n(3),3))
  if (present(v9)) allocate(v9(n(1),n(2),n(3),3))
  if (present(v10)) allocate(v10(n(1),n(2),n(3),3))
  if (present(v1)) call io%bits_mem(storage_size(v1),product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),product(shape(v10)),'v10')
  call trace_end
END SUBROUTINE allocate_vectors

!===============================================================================
SUBROUTINE deallocate_vectors_a (v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:,:), allocatable, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('deallocate_vectors_a')
  if (present(v1)) call io%bits_mem(storage_size(v1),-product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),-product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),-product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),-product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),-product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),-product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),-product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),-product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),-product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),-product(shape(v10)),'v10')
  if (present(v1)) deallocate(v1)
  if (present(v2)) deallocate(v2)
  if (present(v3)) deallocate(v3)
  if (present(v4)) deallocate(v4)
  if (present(v5)) deallocate(v5)
  if (present(v6)) deallocate(v6)
  if (present(v7)) deallocate(v7)
  if (present(v8)) deallocate(v8)
  if (present(v9)) deallocate(v9)
  if (present(v10)) deallocate(v10)
  call trace_end
END SUBROUTINE

!===============================================================================
SUBROUTINE deallocate_vectors (v1, v2, v3, v4, v5, v6, v7, v8, v9, v10)
  integer:: n(3)
  real, dimension(:,:,:,:), pointer, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9, v10
  !.............................................................................
  call trace_begin('deallocate_vectors')
  if (present(v1)) call io%bits_mem(storage_size(v1),-product(shape(v1)),'v1')
  if (present(v2)) call io%bits_mem(storage_size(v2),-product(shape(v2)),'v2')
  if (present(v3)) call io%bits_mem(storage_size(v3),-product(shape(v3)),'v3')
  if (present(v4)) call io%bits_mem(storage_size(v4),-product(shape(v4)),'v4')
  if (present(v5)) call io%bits_mem(storage_size(v5),-product(shape(v5)),'v5')
  if (present(v6)) call io%bits_mem(storage_size(v6),-product(shape(v6)),'v6')
  if (present(v7)) call io%bits_mem(storage_size(v7),-product(shape(v7)),'v7')
  if (present(v8)) call io%bits_mem(storage_size(v8),-product(shape(v8)),'v8')
  if (present(v9)) call io%bits_mem(storage_size(v9),-product(shape(v9)),'v9')
  if (present(v10)) call io%bits_mem(storage_size(v10),-product(shape(v10)),'v10')
  if (present(v1)) deallocate(v1)
  if (present(v2)) deallocate(v2)
  if (present(v3)) deallocate(v3)
  if (present(v4)) deallocate(v4)
  if (present(v5)) deallocate(v5)
  if (present(v6)) deallocate(v6)
  if (present(v7)) deallocate(v7)
  if (present(v8)) deallocate(v8)
  if (present(v9)) deallocate(v9)
  if (present(v10)) deallocate(v10)
  call trace_end
END SUBROUTINE

!===============================================================================
SUBROUTINE vector_assign_scalar (out, in)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3), intent(inout):: out
  out(:,:,:,1) = in
  out(:,:,:,2) = in
  out(:,:,:,3) = in
END SUBROUTINE vector_assign_scalar

!===============================================================================
FUNCTION vector_vector_dot (in, in2) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in,in2
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  out = in(:,:,:,1)*in2(:,:,:,1) &
      + in(:,:,:,2)*in2(:,:,:,2) &
      + in(:,:,:,3)*in2(:,:,:,3)
END FUNCTION vector_vector_dot

!===============================================================================
FUNCTION vector_scalar_mul (in, in2) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(:,:,:), intent(in):: in2
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  out(:,:,:,1) = in(:,:,:,1)*in2
  out(:,:,:,2) = in(:,:,:,2)*in2
  out(:,:,:,3) = in(:,:,:,3)*in2
END FUNCTION
!===============================================================================
FUNCTION scalar_vector_mul (in2, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(:,:,:), intent(in):: in2
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  out(:,:,:,1) = in(:,:,:,1)*in2
  out(:,:,:,2) = in(:,:,:,2)*in2
  out(:,:,:,3) = in(:,:,:,3)*in2
END FUNCTION
!===============================================================================
FUNCTION const3_vector_mul (in1, in) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, dimension(:), intent(in):: in1
  out(:,:,:,1) = in(:,:,:,1)*in1(1)
  out(:,:,:,2) = in(:,:,:,2)*in1(2)
  out(:,:,:,3) = in(:,:,:,3)*in1(3)
END FUNCTION

!===============================================================================
FUNCTION vector_scalar_add (in, in2) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(:,:,:), intent(in):: in2
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  out(:,:,:,1) = in(:,:,:,1)+in2
  out(:,:,:,2) = in(:,:,:,2)+in2
  out(:,:,:,3) = in(:,:,:,3)+in2
END FUNCTION
!===============================================================================
FUNCTION vector_scalar_sub (in, in2) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(:,:,:), intent(in):: in2
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  out(:,:,:,1) = in(:,:,:,1)-in2
  out(:,:,:,2) = in(:,:,:,2)-in2
  out(:,:,:,3) = in(:,:,:,3)-in2
END FUNCTION
!===============================================================================
FUNCTION vector_scalar_divide (in,in2) RESULT (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(:,:,:), intent(in):: in2
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = in(:,:,:,1)/in2
  out(:,:,:,2) = in(:,:,:,2)/in2
  out(:,:,:,3) = in(:,:,:,3)/in2
END FUNCTION

!===============================================================================
FUNCTION real4_cross_vector (in1, in) result (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, intent(in):: in1(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = in1(2)*in(:,:,:,3) - in1(3)*in(:,:,:,2)
  out(:,:,:,2) = in1(3)*in(:,:,:,1) - in1(1)*in(:,:,:,3)
  out(:,:,:,3) = in1(1)*in(:,:,:,2) - in1(2)*in(:,:,:,1)
END FUNCTION real4_cross_vector
!===============================================================================
FUNCTION vector_cross_real4 (in, in2) result (out)
  real, dimension(:,:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  real, intent(in):: in2(3)
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = in2(3)*in(:,:,:,2) - in2(2)*in(:,:,:,3)
  out(:,:,:,2) = in2(1)*in(:,:,:,3) - in2(3)*in(:,:,:,1)
  out(:,:,:,3) = in2(2)*in(:,:,:,1) - in2(1)*in(:,:,:,2)
END FUNCTION vector_cross_real4
!===============================================================================
FUNCTION vector_cross_vector (in, in2) result (out)
  real, dimension(:,:,:,:), intent(in):: in, in2
  real, dimension(size(in,1),size(in,2),size(in,3),3):: out
  !-----------------------------------------------------------------------------
  out(:,:,:,1) = in(:,:,:,2)*in2(:,:,:,3) - in(:,:,:,3)*in2(:,:,:,2)
  out(:,:,:,2) = in(:,:,:,3)*in2(:,:,:,1) - in(:,:,:,1)*in2(:,:,:,3)
  out(:,:,:,3) = in(:,:,:,1)*in2(:,:,:,2) - in(:,:,:,2)*in2(:,:,:,1)
END FUNCTION vector_cross_vector

!===============================================================================
END MODULE vector_mod
