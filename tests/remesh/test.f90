PROGRAM test
  USE io_mod
  USE omp_timer_mod
  USE patch_mod
  USE setup_mod
  USE remesh_mod
  implicit none
  class(patch_t), pointer:: coarse, fine
  integer:: ll(3), i, u(3)
  real:: w(8,8)
  real(8):: time
  !-----------------------------------------------------------------------------
  call setup%init
  allocate (coarse, fine)
  call coarse%init
  call fine%init
  coarse%mem = 1.0
  fine%mem = 0.0
  time = wallclock()
  call remesh%prolong (coarse, fine, coarse%li)
  time = wallclock()-time
  w = remesh%weights(1)
  do i=1,8
    write(stdout,*) 'weight sums:', i, sum(w(i,:))
  end do
  write (stdout,*) 'internal min/max:', fine%fminval(1), fine%fmaxval(1)
  write (stdout,*) 'guard cells:', maxval(fine%mem(1:2,:,:,1,1,1))
  write (stdout,*) 'guard cells:', maxval(fine%mem(:,1:2,:,1,1,1))
  write (stdout,*) 'guard cells:', maxval(fine%mem(:,:,1:2,1,1,1))
  u = fine%mesh%uo
  write (stdout,*) 'guard cells:', maxval(fine%mem(u(1):u(1)+1,:,:,1,1,1))
  write (stdout,*) 'guard cells:', maxval(fine%mem(:,u(2):u(2)+1,:,1,1,1))
  write (stdout,*) 'guard cells:', maxval(fine%mem(:,:,u(3):u(3)+1,1,1,1))
  write (stdout,'(1x,a,f6.3)') 'mus per cell:', time*1e6/product(coarse%n)
END PROGRAM test

