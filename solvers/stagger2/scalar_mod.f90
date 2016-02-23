!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
MODULE scalar_mod
  USE trace_mod
  implicit none
  public
  !
  integer, private:: verbose=0
  type, private:: node_t
    type(node_t), pointer:: next=>null()
  end type
  type(node_t), private, pointer:: head=>null(), tail=>null()
  integer, private, save:: id=0
  !
  integer, private:: nalloc=0
  interface log
    module procedure scalar_log_scalar
  end interface
  interface exp
    module procedure scalar_exp_scalar
  end interface
CONTAINS

FUNCTION scalar_log_scalar(in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  integer:: i
  do i=1,size(in,3)
    out(:,:,i) = log(in(:,:,i))
  end do
END FUNCTION

FUNCTION scalar_exp_scalar(in) RESULT (out)
  real, dimension(:,:,:), intent(in):: in
  real, dimension(size(in,1),size(in,2),size(in,3)):: out
  integer:: i
  do i=1,size(in,3)
    out(:,:,i) = exp(in(:,:,i))
  end do
END FUNCTION

!> =============================================================================
SUBROUTINE allocate_scalars (n, v1, v2, v3, v4, v5, v6, v7, v8, v9)
  integer:: n(3)
  real, dimension(:,:,:), pointer, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9
  !.............................................................................
  call trace_begin('allocate_scalars')
  if (present(v1)) allocate(v1(n(1),n(2),n(3)))
  if (present(v2)) allocate(v2(n(1),n(2),n(3)))
  if (present(v3)) allocate(v3(n(1),n(2),n(3)))
  if (present(v4)) allocate(v4(n(1),n(2),n(3)))
  if (present(v5)) allocate(v5(n(1),n(2),n(3)))
  if (present(v6)) allocate(v6(n(1),n(2),n(3)))
  if (present(v7)) allocate(v7(n(1),n(2),n(3)))
  if (present(v8)) allocate(v8(n(1),n(2),n(3)))
  if (present(v9)) allocate(v9(n(1),n(2),n(3)))
  call trace_end
END SUBROUTINE

!> =============================================================================
SUBROUTINE deallocate_scalars (v1, v2, v3, v4, v5, v6, v7, v8, v9)
  integer:: n(3)
  real, dimension(:,:,:), pointer, optional:: v1, v2, v3, v4, v5, v6, v7, v8, v9
  !.............................................................................
  call trace_begin('init_vectors')
  if (present(v1)) deallocate(v1)
  if (present(v2)) deallocate(v2)
  if (present(v3)) deallocate(v3)
  if (present(v4)) deallocate(v4)
  if (present(v5)) deallocate(v5)
  if (present(v6)) deallocate(v6)
  if (present(v7)) deallocate(v7)
  if (present(v8)) deallocate(v8)
  if (present(v9)) deallocate(v9)
  call trace_end
END SUBROUTINE

!> =============================================================================
END MODULE scalar_mod
