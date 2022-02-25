PROGRAM test
  USE omp_timer_mod
  USE io_mod
  USE patch_mod
  USE poisson_mod
  implicit none
  type (patch_t), pointer:: patch
  type (poisson_t) :: poisson
  real, dimension(:,:,:), pointer :: rho, phi
  integer :: mid(3), iters
  !-----------------------------------------------------------------------------
  ! Allocate and initialize a patch (reads size and dimension from namelist)
  !-----------------------------------------------------------------------------
  call io%init
  allocate (patch)
  call patch%init
  rho => patch%mem(:,:,:,1,1,1)
  phi => patch%mem(:,:,:,2,1,1)
  print*,shape(patch%mem)
  !-----------------------------------------------------------------------------
  ! Set the density to a constant, except for a value 1.0 times larger at the
  ! center
  !-----------------------------------------------------------------------------
  rho = 0.1
  mid = (1+patch%mesh%gn)/2
  rho(mid(1),mid(2),mid(3)) = rho(mid(1),mid(2),mid(3)) + 1.0
  phi = 0.0
  call poisson%init (patch%mesh)
  call tic
  call poisson%solve (patch%mesh, rho, phi, iters)
  call timer%toc ('poisson', product(patch%mesh%n)*iters)
  !call timer%toc ('poisson', product(patch%mesh%n))
END PROGRAM test
