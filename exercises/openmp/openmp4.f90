!*******************************************************************************
!
! Example OpenMP program.  Try compiling this with OpenMP option:
!
!   ifort -openmp openmp4.f90 -o openmp4.x
!   gfortran -fopenmp openmp4.f90 -o openmp4.x
!
! Then run it with one core, and with -- for openmp -- 4 cores::
!
!   export OMP_NUM_THREADS=1; ./openmp4.x      # with bash shell
!   export OMP_NUM_THREADS=4; ./openmp4.x      # with bash shell
!
!   setenv OMP_NUM_THREADS 1; ./openmp4.x      # with tcsh shell
!   setenv OMP_NUM_THREADS 4; ./openmp4.x      # with tcsh shell
!
! You still probably don't get the same answer in the two cases (but you 
! may have to the second case several times -- it may work sometimes; and
! it may work more often with ifort than with gfortran, so try both).
!
! Correct the code!
!
!*******************************************************************************
PROGRAM openmp4
  implicit none
  integer, parameter:: n=100
  real, dimension(n):: a, b, c
  real:: f, s, s1
  integer:: i

  s = 0.
  f = 10.
 
  !$omp  parallel default(none) &
  !$omp&   private(i,s1) shared(a,b,c,f,s)

  !$omp do
  do i=1,n
    b(i) = i-4
    c(i) = 0.01*i**2
  end do
  !$omp enddo

  
  s1 = 0.
  !$omp do
  do i=1,n
     a(i) = b(i) + c(i) * f
     s1 = s1 + a(i)
  enddo
  !$omp enddo
  !$atomic
  s = s + s1
  !$omp end parallel

  print *,'s =',s

END PROGRAM openmp4
