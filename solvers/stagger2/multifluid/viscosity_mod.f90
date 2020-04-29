!===============================================================================
!> This is a start on a module that could become a general viscosity module
!===============================================================================
MODULE viscosity_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE force_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: viscosity_t
    integer:: id
    integer:: i_gas, i_dust, nc, i=1
    character(len=64):: solver
    real(4):: nu
    procedure(constant_viscosity) , pointer:: selected =>null()
    procedure(diffusivity), pointer:: selected2=>null()
    procedure(dispersion) , pointer:: selected3=>null()
  contains
    procedure init
  end type
  logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Initialize viscosity_mod
!===============================================================================
SUBROUTINE init (self, solver)
   class(viscosity_t):: self
   character(len=64):: solver
   character(len=32), save:: type='zero'
   integer, save:: i_gas=1, i_dust=2, nc=1
   real(4), save:: nu=1.0e-5
   namelist /viscosity_params/ type, nu, i_gas, i_dust, nc
   !----------------------------------------------------------------------------
   call trace_begin('viscosity_t%init')
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, viscosity_params)
     if (io%master) write (*, viscosity_params)
   end if
   call io%lock%unset
   self%solver = trim(solver)
   self%nu     = nu
   self%nc     = nc
   self%i_gas  = i_gas
   self%i_dust = i_dust
   select case (type)
   case ('constant_viscosity')
     self%selected => constant_viscosity
   case ('zero')
     nullify(self%selected)
   case default
     print*,'UNKNOWN VISCOSITY TYPE'
   end select
   nullify(self%selected2)
   nullify(self%selected3)
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Apply constant viscosity parameter
!===============================================================================
FUNCTION constant_viscosity (self, time, mem, U, cs, m) RESULT (nu)
  class(viscosity_t)    :: self
  real(8)               :: time
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: nu
  real, dimension(:,:,:):: cs
  real, dimension(:,:,:,:):: mem, U
  !-----------------------------------------------------------------------------
  call trace_begin('viscosity_t:constant_viscosity')
  if (self%i==self%i_gas) then
    nu = self%nu
  else if (self%i==self%i_dust) then
    nu = 0.0
  end if
  self%i = self%i + 1
  if (self%i > self%nc) self%i=1
  call trace_end
END FUNCTION constant_viscosity

FUNCTION diffusivity (self, time, mem, U, m, force) RESULT (dddt)
  class(viscosity_t)    :: self
  real(8)               :: time
  class(mesh_t), pointer, dimension(:) :: m
  type(force_t):: force
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: dddt
  real, dimension(:,:,:,:):: mem, U
  !-----------------------------------------------------------------------------
  dddt = 0.0
END FUNCTION diffusivity

FUNCTION dispersion (self, time, mem, U, m, force) RESULT (dpdt)
  class(viscosity_t)    :: self
  real(8)               :: time
  class(mesh_t), pointer, dimension(:) :: m
  type(force_t):: force
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3) :: dpdt
  real, dimension(:,:,:,:):: mem, U
  !-----------------------------------------------------------------------------
  dpdt = 0.0
END FUNCTION dispersion

END MODULE
