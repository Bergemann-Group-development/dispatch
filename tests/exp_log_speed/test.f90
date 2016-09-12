!> *****************************************************************************
!> Test the speed of exp() and log(), for single and double precision
!> *****************************************************************************
PROGRAM test
  USE omp_mod
  USE omp_timer_mod
  implicit none
  integer(8):: ntot, i
  real(8):: wc, sum8, c8
  real:: sum, c, arg, ex

  call omp%init

  sum = 0.0
  ntot=1e9
  c = 1./ntot
  call tic
  do i=1,ntot
    sum=sum+c*real(i)
  end do
  call toc ('sum(real)', ntot)
  print*,sum

  sum8 = 0.0
  ntot=1e9
  c8 = 1./ntot
  call tic
  do i=1,ntot
    sum8=sum8+c8*real(i,kind=8)
  end do
  call toc ('sum(double)', ntot)
  print*,sum8

  sum = 0.0
  ntot=5e7
  call tic
  do i=1,ntot
    sum=sum+exp(-c*real(i))
  end do
  call toc ('sum(exp(real))', ntot)
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
  call toc ('sum(exp_series(real))', ntot)
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
  call toc ('non_vector_sum(exp(real))', ntot)
  print*,sum

  sum8 = 0.0
  ntot=5e7
  c8 = 1./ntot
  call tic
  do i=1,ntot
    sum8=sum8+exp(-c8*real(i,kind=8))
  end do
  call toc ('sum(exp(double))', ntot)
  print*,sum8

  sum = 0.0
  ntot=5e7
  call tic
  do i=1,ntot
    sum=sum+log(real(i))
  end do
  call toc ('sum(log(real))', ntot)
  print*,sum

  sum8 = 0.0
  ntot=5e7
  call tic
  do i=1,ntot
    sum8=sum8+log(real(i,kind=8))
  end do
  call toc ('sum(log(double))', ntot)
  print*,sum8

END PROGRAM
