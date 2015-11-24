!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!> $Id$
!> |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
PROGRAM experiment
  USE io_mod
  USE mpi_mod
  USE patch_mod
  implicit none
  type(patch_t), dimension(:,:,:), pointer:: patch
  type(patch_t), pointer:: p
  real:: size(3) = 1.0
  integer:: i, j, k, n(3)=[4,4,1]
  !.............................................................................
  call mpi%init
  call io%init
  print*, io%master
  do i=1,n(1)
  do j=1,n(1)
  do k=1,n(3)
    allocate (patch(i,j,k))
    p => patch(i,j,k)
    call p%init(size,n,i,j,k)
    print *,i,j,k,p%mesh%p
  end do
  end do
  end do
  call mpi%end
END PROGRAM experiment
