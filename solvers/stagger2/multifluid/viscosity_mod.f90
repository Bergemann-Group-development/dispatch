!> =============================================================================
!> This is a start on a module that could become a general viscosity module
!> =============================================================================
MODULE viscosity_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: viscosity_t
    integer:: id
    character(len=64):: solver
    real(4):: nu
    procedure(constant_viscosity), pointer:: selected=>null()
  contains
    procedure init
  end type
  logical, save:: first_time=.true.
CONTAINS

!> =============================================================================
!> Initialize viscosity_mod
!> =============================================================================
SUBROUTINE init (self, solver)
   class(viscosity_t):: self
   character(len=64):: solver
   character(len=32), save:: type='zero'
   real(4), save:: nu=0.1
   namelist /viscosity_params/ type, nu
   !----------------------------------------------------------------------------
   call trace_begin('viscosity_t%init')
   !$omp critical (input_cr)
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, viscosity_params)
     if (io%master) write (*, viscosity_params)
   end if
   !$omp end critical (input_cr)
   self%solver = trim(solver)
   self%nu = nu
   select case (type)
   case ('constant_viscosity')
     self%selected => constant_viscosity
   case ('zero')
     nullify(self%selected)
   case default
     print*,'UNKNOWN FORCING TYPE'
   end select
   call trace_end
END SUBROUTINE init

!> =============================================================================
!> Apply constant viscosity parameter
!> =============================================================================
FUNCTION constant_viscosity (self, time, mem, U, cs, m) RESULT (nu)
  class(viscosity_t)    :: self
  real(8)               :: time
  class(mesh_t), pointer, dimension(:) :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn) :: nu
  real, dimension(:,:,:):: cs
  real, dimension(:,:,:,:):: mem, U
  !-----------------------------------------------------------------------------
  call trace_begin('viscosity_t:constant_viscosity')
  nu = self%nu
  call trace_end
END FUNCTION constant_viscosity

END MODULE
