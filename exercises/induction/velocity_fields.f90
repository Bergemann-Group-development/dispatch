MODULE velocity_fields
  character(len=16):: type
CONTAINS

!*******************************************************************************
SUBROUTINE init_velocity (UU,u)
  USE mpi_base, only: master
  USE mesh,     only: m
  USE vectors,  only: vector, allocate_vector, deallocate_vector
  USE io,       only: input
  implicit none
  type(vector):: UU, k
  integer ix,iy,iz
  real:: px, py, pz, u(3)
  real:: pi2 = 8.0*atan(1.0)
  namelist /velocity/ type
!-------------------------------------------------------------------------------
! Select type => coefficients in expressions
!-------------------------------------------------------------------------------
  rewind (input); read(input,velocity); if (master) write(*,velocity)
!-------------------------------------------------------------------------------
! (x,y,z)-cordinates normalized to (-pi,pi)
!-------------------------------------------------------------------------------
  call allocate_vector(k)
  do iz=m(3)%lb,m(3)%ub
  do iy=m(2)%lb,m(2)%ub
    k%x(:,iy,iz) = pi2*m(1)%r/m(1)%s
  end do
  end do
  do iz=m(3)%lb,m(3)%ub
  do ix=m(1)%lb,m(1)%ub
    k%y(ix,:,iz) = pi2*m(2)%r/m(2)%s
  end do
  end do
  do iy=m(2)%lb,m(2)%ub
  do ix=m(1)%lb,m(1)%ub
    k%z(ix,iy,:) = pi2*m(3)%r/m(3)%s
  end do
  end do
!-------------------------------------------------------------------------------
! Model expressions for velocity field UU
!-------------------------------------------------------------------------------
  select case (trim(type))
  case('advect')
    UU%x =  u(1)                        ! constant advection speed
    UU%y =  u(2)                        ! constant advection speed
    UU%z =  u(3)                        ! constant speed along B
  case('converge')
    UU%x = -u(1)*sin(k%x)               ! convergence towards x=y=0
    UU%y = -u(2)*sin(k%y)               ! convergence towards x=y=0
    UU%z =  u(3)                        ! constant speed along B
  case('diverge')
    UU%x =  u(1)*sin(k%x)               ! divergence away from x=y=0
    UU%y =  u(2)*sin(k%y)               ! divergence away from x=y=0
    UU%z =  u(3)                        ! constant speed along B
  case('bend')
    UU%x =  u(1)*sin(k%z)               ! shear in the z-direction
    UU%y =  u(2)*sin(k%z)               ! shear in the z-direction
    UU%z =  u(3)                        ! constant speed along B
  case('twist')
    UU%x = +u(1)*sin(k%y)*sin(k%z)      ! rotation increasing with z
    UU%y = -u(2)*sin(k%x)*sin(k%z)      ! rotation increasing with z
    UU%z =  u(3)                        ! constant speed along B
  case default
    print*,'unknown velocity type',type
  end select
  call deallocate_vector(k)
END SUBROUTINE init_velocity

END MODULE velocity_fields
