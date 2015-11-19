!*******************************************************************************
!
! Example OpenMP program.  Try compiling this with OpenMP option:
!
!   ifort -openmp openmp1.f90 -o openmp1.x
!   gfortran -fopenmp openmp1.f90 -o openmp1.x
!
! Then run it with one core, and with -- for openmp -- 4 cores::
!
!   export OMP_NUM_THREADS=1; ./openmp1.x      # with bash shell
!   export OMP_NUM_THREADS=4; ./openmp1.x      # with bash shell
!
!   setenv OMP_NUM_THREADS 1; ./openmp1.x      # with tcsh shell
!   setenv OMP_NUM_THREADS 4; ./openmp1.x      # with tcsh shell
!
! You probably don't get the same answer.   Correct the code!
!
!*******************************************************************************
PROGRAM openmp1
  implicit none
  integer, parameter:: n=100
  real, dimension(n):: a, b, c
  real:: f, s
  integer:: i

  !$omp parallel default(none) &
  !$omp   private(i,f) shared(a,b,c,s)

  !$omp do
  do i=1,n
    b(i) = i-1
    c(i) = 0.01*i**2
  end do
  !$omp enddo

  f = 10.
  s = 0.
  
  !$omp do
  do i=1,n
     a(i) = b(i) + c(i) * f
     s = s + a(i)
  enddo
  !$omp enddo
  !$omp end parallel

  print *,'s =',s

END PROGRAM openmp1
