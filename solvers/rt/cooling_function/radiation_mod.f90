!===============================================================================
!> Stub for self-contained radiative transfer module
!===============================================================================
MODULE radiation_mod
  USE io_mod
  USE trace_mod
  USE mesh_mod
  USE link_mod
  implicit none
  private
  type, public:: radiation_t
    real, dimension(:,:,:), pointer:: rho      => null() ! pointer to density array
    real, dimension(:,:,:), pointer:: temp     => null() ! pointer to temperature array
    real, dimension(:,:,:), pointer:: qq       => null() ! pointer to cooling / heating
    type(mesh_t), dimension(:), pointer:: mesh => null() ! pointer to mesh
    logical:: do_radiation = .false.
    integer:: type
  contains
    procedure:: init
    procedure:: solve
  end type
  real, save:: T_MC = 5.
  integer, parameter:: dp=4
CONTAINS
  include 'cooling_Dalgarno.f90'
  include 'cooling_Rosen_etal.f90'
  include 'cooling_Shure.f90'
  include 'cooling_Casiana.f90'
  include 'cooling_Gnedin.f90'

!===============================================================================
!> Do one-time work for each instance of this solver, such as allocating Q.
!> Make this a zero-filled array in any case, to avoid having to use if's in
!> the solvers, but allocate a temperature array only if we are going to do
!> actual radiation
!===============================================================================
SUBROUTINE init (self)
  class(radiation_t):: self
  logical, save:: first_time=.true., do_radiation=.false.
  character(len=8), save:: type = 'void'
  integer:: void
  namelist /radiation_params/ do_radiation, type
  associate (mx=>self%mesh(1)%gn, my=>self%mesh(2)%gn, mz=>self%mesh(3)%gn)
  !.............................................................................
  call trace%begin ('radiation_mod::init')
  call io%lock%set ('radiation_t%init')
  if (first_time) then
    rewind (io_unit%input); read (io_unit%input, radiation_params, iostat=void)
    !if (io_unit%master) &
      write (*, radiation_params)
  end if
  call io%lock%unset ('radiation_t%init')
  self%do_radiation = do_radiation
  select case (type)
  case ("Gnedin")
    self%type = 5
  case ("Casiana")
    self%type = 4
  case ("Shure")
    self%type = 3
  case ("Rosen")
    self%type = 2
  case ("Dalgarno")
    self%type = 1
  case default
    self%type = 0
    if (first_time .and. io_unit%master) &
      print*,'WARNING, unknown cooling function type: '//trim(type)
  end select
  if (do_radiation) then
    allocate (self%temp(mx,my,mz))
  end if
  allocate (self%qq(mx,my,mz))
  self%qq = 0.0
  first_time = .false.
  call trace%end()
  end associate
END SUBROUTINE init

!===============================================================================
!> Compute the cooling, assuming self has valid pointers
!===============================================================================
SUBROUTINE solve (self)
  class(radiation_t):: self
  integer:: ix, iy, iz
  real:: qq, scaling
  !.............................................................................
  call trace%begin ('radiation_mod::solve')
  !-----------------------------------------------------------------------------
  ! The scaling parameter needs to be specified
  !-----------------------------------------------------------------------------
  scaling = 1.0
  !-----------------------------------------------------------------------------
  ! This would be faster if the functions accpted an array argument
  !-----------------------------------------------------------------------------
  do iz=1,self%mesh(3)%gn 
  do iy=1,self%mesh(2)%gn 
  do ix=1,self%mesh(1)%gn 
    associate (temp=>self%temp(ix,iy,iz))
    select case (self%type)
    case (0)
      qq = 0.0
    case (1)
      qq = cooling_Dalgarno (temp)
    case (2)
      qq = cooling_Rosen_etal (temp)
    case (3)
      qq = cooling_Casiana (temp)
    case (4)
      qq = cooling_Shure (temp)
    case (5)
      qq = cooling_Gnedin (temp)
    end select
    end associate
    self%qq = -scaling*qq*self%rho**2
  end do
  end do
  end do
  call trace%end()
END SUBROUTINE solve

END MODULE radiation_mod
