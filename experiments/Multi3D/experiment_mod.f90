!===============================================================================
!> Top level MHD task definition -- can either rely entirely on lower levels,
!> including extras_t%initial_conditions, or can optionally do such things here
!===============================================================================
MODULE experiment_mod
  USE io_mod
  USE bits_mod
  USE trace_mod
  USE solver_mod
  USE multi3d_mod
  implicit none
  private
  type, public, extends(solver_t):: experiment_t
  contains
    procedure init
    procedure update
    procedure read_txt_atmos
  end type
CONTAINS

!===============================================================================
!> Setup a uniform initial state, with density=d0 and B=B0
!===============================================================================
SUBROUTINE init (self)
  class(experiment_t):: self
  integer, save:: itimer=0
  !----------------------------------------------------------------------------
  call trace_begin('experiment_mod::init', itimer=itimer)
  call self%read_txt_atmos()
  !call self%solver_t%init
  call trace_end (itimer=itimer)
END SUBROUTINE init

!===============================================================================
!> This update only concerns the MHD task; the RT tasks are called directly
!===============================================================================
SUBROUTINE update (self)
  class(experiment_t):: self
  !----------------------------------------------------------------------------
  call trace_begin('experiment_t%update')
  !call self%solver_t%update
  call trace_end
END SUBROUTINE update

!===============================================================================
!> Reads a tabulated model atmosphere from a text file
!===============================================================================
SUBROUTINE read_txt_atmos(self)
  use m3io_mod, only: NextEntry

  class(experiment_t) :: self
  character(len=126) :: atmos_file, atmos_name
  integer :: lunit, nx, ny, nz, i
  real(8), dimension(:), allocatable :: zi, temp_1d, pe_1d, rho_1d, vturb_1d
  real(8), dimension(:,:,:), allocatable :: temp, pe, rho, vturb
  namelist /input_params/ atmos_file
  !----------------------------------------------------------------------------
  open(file='atmos_input.nml', newunit=lunit)
  read(unit=lunit, nml=input_params) ! reads atmos_file
  close(lunit)

  open(file=atmos_file, newunit=lunit)

  nx = 5
  ny = 5

  call nextentry(lunit)
  read(lunit, *) atmos_name
  call nextentry(lunit)
  read(lunit, *) nz

  print *, 'Reading atmosphere ' // atmos_name

  allocate(zi(nz))
  allocate(temp_1d, pe_1d, rho_1d, vturb_1d, mold=zi)
  allocate(temp(nx,ny,nz))
  allocate(pe, rho, vturb, mold=temp)

  do i=1, nz
    call nextentry(lunit)
    read(lunit, *) zi(i), temp_1d(i), pe_1d(i), rho_1d(i), vturb_1d(i)
    
    temp(:,:,i)  = temp_1d(i)
    pe(:,:,i)    = pe_1d(i)
    rho(:,:,i)   = rho_1d(i)
    vturb(:,:,i) = vturb_1d(i)
  end do
END SUBROUTINE

END MODULE
