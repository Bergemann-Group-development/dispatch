!===============================================================================
!> This is a start on a module that could become a general turbulence forcing
!> module, where the standard forcing we have been using in RAMSES is implemented
!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE vector_mod
  USE mesh_mod
  USE random_mod
  implicit none
  private
  real(8), parameter:: pi2=8.0*atan(1.0)
  type, public:: force_t
    real(8)                          :: grav
    procedure(const_gravity), pointer:: selected=>null()
  contains
    procedure init
  end type
 logical, save:: first_time=.true.
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
   class(force_t):: self
   real(8), save               :: grav
   character(len=32), save:: type='zero'
   namelist /force_params/ type, grav
   !----------------------------------------------------------------------------
   call trace_begin('force_t%init')
   call io%lock%set
   if (first_time) then
     first_time = .false.
     rewind (io%input)
     read (io%input, force_params)
     if (io%master) write (*, force_params)
   end if
   call io%lock%unset
   self%grav = grav
   select case (type)
   case ('const_gravity')
     self%selected => const_gravity
   case ('zero')
     nullify(self%selected)
   case default
     print*,'UNKNOWN FORCING TYPE'
   end select
   call trace_end
END SUBROUTINE init

!===============================================================================
!> Forcing a constant gravitational acceleration in x
!===============================================================================
SUBROUTINE const_gravity (self, time, m, dd, Ux, Uy, Uz, dpdt)
  class(force_t)           :: self
  real(8)                  :: time
  class(mesh_t), pointer   :: m(:)
  real, dimension(:,:,:,:) :: dd, dpdt
  real, dimension(:,:,:)   :: Ux, Uy, Uz
  !-----------------------------------------------------------------------------
  call trace_begin('force_t:const_gravity')
  dpdt(:,:,:,1) = dpdt(:,:,:,1) + self%grav*dd(:,:,:,1)
  call trace_end
END SUBROUTINE const_gravity

END MODULE
