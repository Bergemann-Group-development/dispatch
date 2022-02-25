!*******************************************************************************
!> Test the speed of exp() and log(), for single and double precision
!*******************************************************************************
SUBROUTINE cases
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer(8):: ntot, i
  real(8):: wc, sum8, c8
  real:: sum, c, arg, ex, b

  call omp%init

  sum = 0.0
  ntot=1e9
  c = 1./ntot
  call tic
  do i=1,ntot
    sum=sum+c*real(i)
  end do
  call timer%toc ('sum(real)', ntot)
  print*,sum

  sum8 = 0.0
  ntot=1e9
  c8 = 1./ntot
  call tic
  do i=1,ntot
    sum8=sum8+c8*real(i,kind=8)
  end do
  call timer%toc ('sum(double)', ntot)
  print*,sum8

  sum = 0.0
  ntot=5e7
  call tic
  do i=1,ntot
    sum=sum+exp(-c*real(i))
  end do
  call timer%toc ('sum(exp(real))', ntot)
  print*,sum

  sum = 0.0
  ntot=5e7
  b = 1.4426950*c
  call tic
  do i=1,ntot
#if defined __INTEL_COMPILER
    sum=sum+2.**(-b*real(i))
#elif defined __GFORTRAN__
    sum=sum+exp(-c*real(i))
#endif
  end do
  call timer%toc ('sum(exp2(real))', ntot)
  print*,sum

  sum = 0.0
  ntot=5e7
  c = 1./ntot
  call tic
  do i=1,ntot
    arg = c*real(i)
    if (arg < 1e-3) then
      ex = 1.-arg
    else
      ex = exp(-arg)
    end if
    sum=sum+ex
  end do
  call timer%toc ('sum(exp_series(real))', ntot)
  print*,sum

  sum = 0.0
  ntot=5e7
  c = 1./ntot
  call tic
  do i=1,ntot
    arg = c*real(i)
    if (arg < 0.0) then
      ex = 1.-arg
    else
      ex = exp(-arg)
    end if
    sum=sum+ex
  end do
  call timer%toc ('non_vector_sum(exp(real))', ntot)
  print*,sum

  sum8 = 0.0
  ntot=5e7
  c8 = 1./ntot
  call tic
  do i=1,ntot
    sum8=sum8+exp(-c8*real(i,kind=8))
  end do
  call timer%toc ('sum(exp(double))', ntot)
  print*,sum8

  sum = 0.0
  ntot=5e7
  call tic
  do i=1,ntot
    sum=sum+log(real(i))
  end do
  call timer%toc ('sum(log(real))', ntot)
  print*,sum

  sum8 = 0.0
  ntot=5e7
  call tic
  do i=1,ntot
    sum8=sum8+log(real(i,kind=8))
  end do
  call timer%toc ('sum(log(double))', ntot)
  print*,sum8

END SUBROUTINE
