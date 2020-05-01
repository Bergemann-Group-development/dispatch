!===============================================================================
MODULE force_mod
  USE io_mod
  USE trace_mod
  USE kinds_mod
  USE mesh_mod
  implicit none
  private

  type, public:: force_t
    integer:: id
    character(len=64):: solver
    real(4):: y_grav_acc
    procedure(vertical_acceleration), pointer:: selected=>null()
  contains
    procedure init
  end type
 logical, save:: first_time=.true.
CONTAINS

!===============================================================================
SUBROUTINE init (self, solver)
  class(force_t):: self
  integer:: iostat
  character(len=64):: solver
  real(4):: y_grav_acc = -0.5
  namelist /force_params/ y_grav_acc
  !----------------------------------------------------------------------------
  call trace%begin('force_t%init')
  call io%lock%set
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, force_params, iostat=iostat)
    if (io%master) write (*, force_params)
  end if
  call io%lock%unset

  self%y_grav_acc = y_grav_acc
  self%selected => vertical_acceleration
  self%solver = trim(solver)
  call trace%end()
END SUBROUTINE init

!===============================================================================
FUNCTION vertical_acceleration (self, time, d, p, Ux, Uy, Uz, m) RESULT (ff)
  class(force_t)                                :: self
  real(8)                                       :: time
  real, dimension(:,:,:)  , pointer             :: Ux, Uy, Uz
  real, dimension(:,:,:)  , pointer, contiguous :: d
  real, dimension(:,:,:,:), pointer, contiguous :: p
  class(mesh_t), dimension(:), pointer          :: m
  real, dimension(m(1)%gn,m(2)%gn,m(3)%gn,3)    :: ff
  !.............................................................................
  call trace%begin('force_t:vertical_acceleration')

  ! apply uniform gravitational force in the y-direction
  ff(:,:,:,1) = 0.0
  ff(:,:,:,2) = self%y_grav_acc
  ff(:,:,:,3) = 0.0

  call trace%end()
END FUNCTION vertical_acceleration

END MODULE force_mod
