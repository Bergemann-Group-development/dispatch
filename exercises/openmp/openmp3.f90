!*******************************************************************************
!
! Example OpenMP program.  Try compiling this with OpenMP option:
!
!   ifort -openmp openmp3.f90 -o openmp3.x
!   gfortran -fopenmp openmp3.f90 -o openmp3.x
!
! Then run it with one core, and with -- for openmp -- 4 cores::
!
!   export OMP_NUM_THREADS=1; ./openmp3.x      # with bash shell
!   export OMP_NUM_THREADS=4; ./openmp3.x      # with bash shell
!
!   setenv OMP_NUM_THREADS 1; ./openmp3.x      # with tcsh shell
!   setenv OMP_NUM_THREADS 4; ./openmp3.x      # with tcsh shell
!
! You still probably don't get the same answer in the two cases (but you 
! may have to the second case several times -- it may work sometimes).
!
! Correct the code!
!
!*******************************************************************************
PROGRAM openmp3
  implicit none
  integer, parameter:: n=100
  real, dimension(n):: a, b, c
  real:: f, s
  integer:: i

  !$omp parallel default(none) &
  !$omp   private(i) shared(a,b,c,f,s)

  !$omp do
  do i=1,n
    b(i) = i-4
    c(i) = 0.01*i**2
  end do
  !$omp enddo

  f = 10.
  s = 0.
  
  !$omp do
  do i=1,n
     a(i) = b(i) + c(i) * f
     !$omp critical
     s = s + a(i)
     !$omp end critical
  enddo
  !$omp enddo
  !$omp end parallel

  print *,'s =',s

END PROGRAM openmp3
