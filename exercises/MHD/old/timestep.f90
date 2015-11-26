! $Id$
!************************************************************************
MODULE timestep
  USE mhd_mod
  USE pde_mod
  implicit none
  real, dimension(3) ::  alpha, beta
  integer:: timeorder=3
  type(mhd_t):: d
PUBLIC
  real(kind=8):: courant_number=0.2
  real(kind=8):: t, dt
CONTAINS

!************************************************************************
SUBROUTINE init_timestep
  implicit none
  namelist /time/ courant_number, timeorder
!........................................................................
  call trace_begin ('init_timestep')
  rewind(input); read(input,time); if(master) write(*,time)
  if      (timeorder==1) then
                              alpha=(/ 0.00000,  0.000000,  0.000000 /)
                              beta =(/ 1.00000,  0.000000,  0.000000 /)
  else if (timeorder==2) then
                              alpha=(/ 0.00000, -0.500000,  0.000000 /)
                              beta =(/ 0.50000,  1.000000,  0.000000 /)
  else if (timeorder==3) then
                              alpha=(/ 0.00000, -0.641874, -1.310210 /)
                              beta =(/ 0.46173,  0.924087,  0.390614 /)
  end if
  call allocate_mhd (d)
  t = 0d0
  call trace_end
END SUBROUTINE init_timestep

!************************************************************************
SUBROUTINE time_step (f)
  implicit none
  type(mhd_t):: f
  integer:: iz, substep
  character(len=mch):: id="timestep $Id$"
!........................................................................
  call print_id (id)
  call trace_begin ('timestep')
!------------------------------------------------------------------------
!  Loop over Runge-Kutta stages
!------------------------------------------------------------------------
  do substep=1,timeorder
    if (substep==1) then
      !$omp parallel do
      do iz=m(3)%lb,m(3)%ub
          d%d(:,:,iz) = 0.
          d%s(:,:,iz) = 0.
        d%p%x(:,:,iz) = 0.
        d%p%y(:,:,iz) = 0.
        d%p%z(:,:,iz) = 0.
        d%B%x(:,:,iz) = 0.
        d%B%y(:,:,iz) = 0.
        d%B%z(:,:,iz) = 0.
      end do
    else
      !$omp parallel do
      do iz=m(3)%lb,m(3)%ub
          d%d(:,:,iz) = alpha(substep)*  d%d(:,:,iz)
          d%s(:,:,iz) = alpha(substep)*  d%s(:,:,iz)
        d%p%x(:,:,iz) = alpha(substep)*d%p%x(:,:,iz)
        d%p%y(:,:,iz) = alpha(substep)*d%p%y(:,:,iz)
        d%p%z(:,:,iz) = alpha(substep)*d%p%z(:,:,iz)
        d%B%x(:,:,iz) = alpha(substep)*d%B%x(:,:,iz)
        d%B%y(:,:,iz) = alpha(substep)*d%B%y(:,:,iz)
        d%B%z(:,:,iz) = alpha(substep)*d%B%z(:,:,iz)
      end do
    end if
    call pde (f, d)
    if (substep==1) call courant (f, courant_number, dt)
      !$omp parallel do
    do iz=m(3)%lb,m(3)%ub
        f%d(:,:,iz) =   f%d(:,:,iz) + (dt*beta(substep))*  d%d(:,:,iz)
        f%s(:,:,iz) =   f%s(:,:,iz) + (dt*beta(substep))*  d%s(:,:,iz)
      f%p%x(:,:,iz) = f%p%x(:,:,iz) + (dt*beta(substep))*d%p%x(:,:,iz)                      
      f%p%y(:,:,iz) = f%p%y(:,:,iz) + (dt*beta(substep))*d%p%y(:,:,iz)
      f%p%z(:,:,iz) = f%p%z(:,:,iz) + (dt*beta(substep))*d%p%z(:,:,iz)
      f%B%x(:,:,iz) = f%B%x(:,:,iz) + (dt*beta(substep))*d%B%x(:,:,iz)
      f%B%y(:,:,iz) = f%B%y(:,:,iz) + (dt*beta(substep))*d%B%y(:,:,iz)
      f%B%z(:,:,iz) = f%B%z(:,:,iz) + (dt*beta(substep))*d%B%z(:,:,iz)
    end do
  end do
  t = t + dt
  call trace_end
END SUBROUTINE time_step
END MODULE timestep
