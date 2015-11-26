!> =============================================================================
!> To use this module to generate a velocity field, do
!> 
!> USE velocity_fields
!> type(velocity_t):: vel
!> integer:: dims(3)                    ! dimensions
!> real:: u(3)                          ! velocity amplitudes
!> ...
!> call vel%init(self%n)
!>   ... = ... vel%ux(:,:,:) ...
!>   ... = ... vel%uy(:,:,:) ...
!>   ... = ... vel%ux(:,:,:) ...
!>
!> The input.nml file should contain entries such as
!>
!> &velocity_params
!>  type='advect' u=0.5,0.8,0.4
!> /
!> =============================================================================
MODULE velocity_fields
  USE io_mod
  USE mpi_mod
  implicit none
  private
  type, public:: velocity_t
    real, dimension(:,:,:), pointer:: Ux, Uy, Uz
    character(len=16):: type
    integer:: n(3)
  contains
    procedure:: init
  end type
  character(len=16):: type
  real, dimension(:,:,:), pointer:: kx, ky, kz
CONTAINS

!*******************************************************************************
SUBROUTINE init (s, n)
  class(velocity_t):: s
  integer n(3), ix, iy, iz
  real:: u(3)
  real:: pi2 = 8.0*atan(1.0)
  namelist /velocity_params/ type, u
!-------------------------------------------------------------------------------
! Select type => coefficients in expressions
!-------------------------------------------------------------------------------
  rewind (io%input); read(io%input,velocity_params)
  if (io%master) write(*,velocity_params)
!-------------------------------------------------------------------------------
! (x,y,z)-cordinates normalized to (-pi,pi)
!-------------------------------------------------------------------------------
  allocate (  kx(n(1),n(2),n(3)),   ky(n(1),n(2),n(3)),   kz(n(1),n(2),n(3)))
  allocate (s%ux(n(1),n(2),n(3)), s%uy(n(1),n(2),n(3)), s%uz(n(1),n(2),n(3)))
  s%n = n
  s%type = type
  do iz=1,n(1)
  do iy=1,n(2)
  do ix=1,n(3)
    kx(ix,iy,iz) = ix*pi2/n(1)
    ky(ix,iy,iz) = iy*pi2/n(2)
    kz(ix,iy,iz) = iz*pi2/n(3)
  end do
  end do
  end do
!-------------------------------------------------------------------------------
! Model expressions for velocity field UU
!-------------------------------------------------------------------------------
  select case (trim(type))
  case('advect')
    s%Ux =  u(1)                      ! constant advection speed
    s%Uy =  u(2)                      ! constant advection speed
    s%Uz =  u(3)                      ! constant speed along B
  case('converge')
    s%Ux = -u(1)*sin(kx)              ! convergence towards x=y=0
    s%Uy = -u(2)*sin(ky)              ! convergence towards x=y=0
    s%Uz =  u(3)                      ! constant speed along B
  case('diverge')
    s%Ux =  u(1)*sin(kx)              ! divergence away from x=y=0
    s%Uy =  u(2)*sin(ky)              ! divergence away from x=y=0
    s%Uz =  u(3)                      ! constant speed along B
  case('bend')
    s%Ux =  u(1)*sin(kz)              ! shear in the z-direction
    s%Uy =  u(2)*sin(kz)              ! shear in the z-direction
    s%Uz =  u(3)                      ! constant speed along B
  case('twist')
    s%Ux = +u(1)*sin(ky)*sin(kz)      ! rotation increasing with z
    s%Uy = -u(2)*sin(kx)*sin(kz)      ! rotation increasing with z
    s%Uz =  u(3)                      ! constant speed along B
  case default
    print*,'unknown velocity type',type
  end select
  deallocate (kx, ky, kz)
END SUBROUTINE init

END MODULE velocity_fields
