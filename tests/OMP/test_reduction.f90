USE mpi_timer
implicit none
real, dimension(:), allocatable:: a, b, c
integer:: i, n
real:: s, sp

n = 1e8
allocate (a(n), b(n), c(n))

!$omp parallel do
do i=1,n
  a(i) = i
  b(i) = 0.1*(n-i)
  c(i) = ((i-1.0)/n)**2
end do

call tic
s = 0.0
!$omp parallel do reduction(+:s)
do i=1,n
  a(i) = b(i) + c(i)
  s = s + a(i)
end do
call toc ('omp parallel do           ', n)
print *,s

call tic
s = 0.0
!$omp parallel default(none) private(i,sp) shared(a,b,c,n,s)
sp = 0.0
!$omp do
do i=1,n
  a(i) = b(i) + c(i)
  sp = sp + a(i)
end do
!$omp end do
!$omp atomic
s = s + sp
!$omp end parallel
call toc ('omp parallel, omp do      ', n)
print *,s

END
