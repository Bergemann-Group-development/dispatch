MODULE patch_mod
  USE omp_lib
  implicit none
  type patch_t
    real, dimension(:,:,:,:,:), pointer:: mem
  contains
    procedure:: allocate
    procedure:: reallocate
  end type
CONTAINS

SUBROUTINE allocate (self, n, time)
  class(patch_t):: self
  integer:: n(5)
  real(8):: time
  !-----------------------------------------------------------------------------
  time = omp_get_wtime()
  allocate (self%mem(n(1),n(2),n(3),n(4),n(5)))
  self%mem = 0.0
  time = omp_get_wtime() - time
END SUBROUTINE allocate 

SUBROUTINE reallocate (self, n, time)
  class(patch_t):: self
  integer:: n(5)
  real(8):: time
  integer:: m(5)
  real, dimension(:,:,:,:,:), pointer:: mem
  !-----------------------------------------------------------------------------
  time = omp_get_wtime()
  m = shape(self%mem)
  allocate (mem(n(1),n(2),n(3),n(4),n(5)))
  mem(1:m(1),1:m(2),1:m(3),1:m(4),1:m(5)) = self%mem
  deallocate (self%mem)
  self%mem => mem
  time = omp_get_wtime() - time
END SUBROUTINE reallocate 

END MODULE patch_mod

PROGRAM TEST
  USE patch_mod
  integer:: n(5)
  type(patch_t):: patch
  real(8):: total, time
  integer(8):: ntotal, ni
  n = [20,20,20,5,5]
  ni = 100
  total = 0d0
  do i=1,ni
    call patch%allocate (n, time)
    total = total+time
    deallocate (patch%mem)
  end do
  ntotal = product(n(1:3))*ni
  print '(f6.3,"s",i12,f10.1,2x,a)', total, ntotal, total*1e9/ntotal, 'ns/pt allocate'
  total = 0d0
  do i=1,ni
    n = [20,20,20,5,5]
    call patch%allocate (n, time)
    n = n + [2,3,0,0,0]
    call patch%reallocate (n, time)
    total = total+time
    deallocate (patch%mem)
  end do
  ntotal = product(n(1:3))*ni
  print '(f6.3,"s",i12,f10.1,2x,a)', total, ntotal, total*1e9/ntotal, 'ns/pt deallocate'
END PROGRAM
