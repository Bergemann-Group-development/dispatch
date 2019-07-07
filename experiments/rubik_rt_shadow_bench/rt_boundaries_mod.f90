!===============================================================================
!> Boundary conditions for radiative heat transfer
!===============================================================================
MODULE rt_boundaries_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mpi_mod
  USE bits_mod
  USE mesh_mod
  USE eos_mod
  USE units_mod
  USE scaling_mod
  USE math_mod
  implicit none
  private
  !-----------------------------------------------------------------------
  type, public :: rt_boundaries_t
    class(mesh_t), pointer:: mesh(:)
    real:: cool, hot
  contains
    procedure:: init
    procedure:: condition                                               ! modify outer BCs if necessary
  end type
!-----------------------------------------------------------------------
  type(rt_boundaries_t), public :: rt_boundaries
CONTAINS

!=======================================================================
!> Initialize the boundary data type
!=======================================================================
SUBROUTINE init (self, mesh)
  class(rt_boundaries_t):: self
  class(mesh_t), pointer:: mesh(:)
  !.....................................................................
  real, dimension(:,:,:), pointer :: d, e
  real, dimension(:,:,:,:), pointer :: src, rk
  real :: t, d0
  integer:: m(4)
  logical, save:: first_time
  real, save :: hot_temp = 1000, cool_temp = 10
  real, save :: cool = 0.0, hot = 0.0
  namelist /rt_bc_params/ cool_temp, hot_temp
  !-----------------------------------------------------------------------------
  self%mesh => mesh
  !
  !$omp critical (input_cr)
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, rt_bc_params)
    if (io%master) write (*, rt_bc_params)
  end if
  !$omp end critical (input_cr)
  allocate(d(1,1,2))
    allocate(e(1,1,2))
    allocate(src(1,1,2, eos%n_lambda))
    allocate(rk(1,1,2, eos%n_lambda))
    d = 1.0/scaling%d
    e(:,:,1) = (cool_temp/scaling%temp)/(eos%gamma-1.0)
    e(:,:,2) = (hot_temp/scaling%temp)/(eos%gamma-1.0)
    m = shape(rk)
    call eos%lookup_table (m, d=d, e=e, rk=rk, src=src)
    cool = 0.0      
    hot  = src(1,1,2,1) - src(1,1,1,1)          
    deallocate(d,e,src, rk)
  self%cool = cool
  self%hot  = hot
END SUBROUTINE init

!===============================================================================
!===============================================================================
SUBROUTINE condition (self, mem, src, rk,tt, mu, w)
  class(rt_boundaries_t) :: self
  real, dimension(:,:,:,:):: mem, src, rk
  real, dimension(:,:,:), pointer:: tt
  real   :: mu, w, tau
  !.....................................................................
  class(mesh_t), pointer:: m1, m2, m3
  integer:: i3, i3b, i2, i1, n(4), iv
  !---------------------------------------------------------------------
  n = shape(rk)
  
  if (self%mesh(3)%lower_boundary) then
    i3 =  self%mesh(3)%lo
    i3b = self%mesh(3)%lb
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i2,i3b:i3,iv) = self%cool*rk(i1,i2,i3b:i3,iv)!/abs(mu)/w*math%pi4
      end do
      end do
      end do
  end if
  !print *, 'lower bnd', self%mesh(1)%lower_boundary,self%mesh(2)%lower_boundary,self%mesh(3)%lower_boundary
  !print *, 'upper bnd', self%mesh(1)%upper_boundary,self%mesh(2)%upper_boundary,self%mesh(3)%upper_boundary
  
  if (self%mesh(3)%upper_boundary) then
    i3  = self%mesh(3)%uo
    i3b = self%mesh(3)%ub
      do i2=self%mesh(2)%lb,self%mesh(2)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i2,i3:i3b,iv) = self%cool*rk(i1,i2,i3:i3b,iv)!/abs(mu)/w*math%pi4
      end do
      end do
      end do
  end if
!
if (self%mesh(2)%lower_boundary) then
    i3  = self%mesh(2)%lo
    i3b = self%mesh(2)%lb
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3b:i3,i2,iv) = self%cool*rk(i1,i3b:i3,i2,iv)!/abs(mu)/w*math%pi4
      end do
      end do
      end do
  end if
  
  if (self%mesh(2)%upper_boundary) then
    i3  = self%mesh(2)%uo
    i3b = self%mesh(2)%ub
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(1)%lb,self%mesh(1)%ub
      do iv=1, n(4)
        mem(i1,i3:i3b,i2,iv) = self%cool*rk(i1,i3:i3b,i2,iv)!/abs(mu)/w*math%pi4
      end do
      end do
      end do
  end if
  !
  if (self%mesh(1)%lower_boundary) then
    i3  = self%mesh(1)%lo
    i3b = self%mesh(1)%lb
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3b:i3,i1,i2,iv) = self%hot*rk(i3b:i3,i1,i2,iv)!/abs(mu)/w*math%pi4
        !print *, 'boundary set',i3b, i1, i2, mem(i3b:i3,i1,i2,iv)
      end do
      end do
      end do
  end if
  
  if (self%mesh(1)%upper_boundary) then
    i3  = self%mesh(1)%uo
    i3b = self%mesh(1)%ub
      do i2=self%mesh(3)%lb,self%mesh(3)%ub
      do i1=self%mesh(2)%lb,self%mesh(2)%ub
      do iv=1, n(4)
        mem(i3:i3b,i1,i2,iv) = self%cool*rk(i3:i3b,i1,i2,iv)!/abs(mu)/w*math%pi4
      end do
      end do
      end do
  end if  
END SUBROUTINE condition

END MODULE rt_boundaries_mod

