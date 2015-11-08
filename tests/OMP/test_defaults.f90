USE mpi_timer
implicit none
real, dimension(:), allocatable:: a, b, c
integer:: i, n

n = 1e8
allocate (a(n), b(n), c(n))

!$omp parallel do
do i=1,n
  a(i) = i
  b(i) = 0.1*(n-i)
  c(i) = ((i-1.0)/n)**2
end do

call tic
!$omp parallel do
do i=1,n
  a(i) = b(i) + c(i)
end do
call toc ('omp parallel do           ', n)

call tic
!$omp parallel do private(i)
do i=1,n
  a(i) = b(i) + c(i)
end do
call toc ('omp parallel do private(i)', n)

call tic
!$omp parallel default(none) private(i) shared(a,b,c,n)
!$omp do
do i=1,n
  a(i) = b(i) + c(i)
end do
!$omp end do
!$omp end parallel
call toc ('omp parallel, omp do      ', n)

END
