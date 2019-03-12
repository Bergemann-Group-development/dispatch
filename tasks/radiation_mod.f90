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
  contains
    procedure:: init
    procedure:: solve
  end type
CONTAINS

!===============================================================================
!> Do one-time work for each instance of this solver, such as allocating Q.
!> Make this a zero-filled array in any case, to avoid having to use if's in
!> the solvers, but allocate a temperature array only if we are going to do
!> actual radiation
!===============================================================================
SUBROUTINE init (self)
  class(radiation_t):: self
  logical, save:: first_call=.true., do_radiation=.false.
  integer:: void
  namelist /radiation_params/ do_radiation
  associate (mx=>self%mesh(1)%gn, my=>self%mesh(2)%gn, mz=>self%mesh(3)%gn)
  !.............................................................................
  call trace%begin ('radiation_mod::init')
  !$omp critical (radiation_cr)
  if (first_call) then
    first_call = .false.
    rewind (io_unit%input); read (io_unit%input, radiation_params, iostat=void)
    if (io_unit%master) write (io_unit%output, radiation_params)
  end if
  !$omp end critical (radiation_cr)
  self%do_radiation = do_radiation
  if (do_radiation) then
    allocate (self%temp(mx,my,mz))
  end if
  allocate (self%qq(mx,my,mz))
  self%qq = 0.0
  call trace%end
  end associate
END SUBROUTINE init

!===============================================================================
!> Solve for the cooling / heating, assuming self has pointers to whatever is
!> relevant (mesh, density, temperature)
!===============================================================================
SUBROUTINE solve (self)
  class(radiation_t):: self
  !.............................................................................
  call trace%begin ('radiation_mod::solve')
  !-----------------------------------------------------------------------------
  ! Here we could compute cooling via an optically thin cooling function, or
  ! solve the full problem
  !-----------------------------------------------------------------------------
  self%qq = 0.0
  call trace%end
END SUBROUTINE solve

END MODULE radiation_mod
