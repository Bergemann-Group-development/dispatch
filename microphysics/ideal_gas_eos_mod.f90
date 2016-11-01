!> =====================================================================
!> EOS module, assuming an ideal gas law
!> =====================================================================
MODULE eos_mod
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  implicit none
  private
  logical, save:: first_time=.true.
  !---------------------------------------------------------------------
  ! Public part
  !---------------------------------------------------------------------
  type, public:: eos_t
    logical:: on = .false.
    integer :: eos_type = 1
    procedure(lookup_stagger), pointer :: lookup => null()
  contains
    procedure:: init
  end type
  type (eos_t), public:: eos
CONTAINS

!> =====================================================================
!> Initialize EOS module (only once)
!> =====================================================================
SUBROUTINE init (self, solver)
  class (eos_t):: self
  character(len=8) :: solver 
  !---------------------------------------------------------------------
  !$omp critical (eos_cr)
  if (first_time) then
    first_time = .false.
    self%on = .true.
    select case(solver)
    case('stagger ')
      self%lookup => lookup_stagger
    case('stagger2')
      self%lookup => lookup_stagger
    case('zeus    ')
      self%lookup => lookup_zeus
    case('ramses  ')
      self%lookup => lookup_ramses
    case default
      error stop "ideal_gas_eos::init: unknown solver or solver not set up"
    end select
  end if
  !$omp end critical (eos_cr)
  !$omp barrier
END SUBROUTINE init

!> =====================================================================
!> General EOS lookup routine, returning values in optional arguments
!> rho - density; s - entropy; pg - gas pressure; tt - gas temperature;
!> ne - electron pressure; rk - rho*kappa; src - source function;
!> gam - adiabatic index, csound - sound speed.
!> =====================================================================
SUBROUTINE lookup_stagger (self, rho, ee, pg, tt, ne, rk, src, gam, csound)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in) :: rho, ee
  real, dimension(:,:,:),   optional :: pg, tt, ne, rk
  real, dimension(:,:,:,:), optional :: src
  real, optional :: gam, csound
  integer:: x, y, z
  !.....................................................................

  !---------------------------------------------------------------------
  ! Gas pressure
  !---------------------------------------------------------------------
  if (present(pg)) then
    if (gam==1.0) then
      pg = rho*csound**2
    else
      pg = rho**gam*exp(ee/rho*(gam-1))
    end if
  end if
  !---------------------------------------------------------------------
  ! Gas temperature
  !---------------------------------------------------------------------
  if (present(tt)) then
    if (present(pg)) then
      tt = pg/rho
    else
      if (gam==1.0) then
        tt = csound**2
      else
        do x = 1,size(rho,1)
        do y = 1,size(rho,2)
        do z = 1,size(rho,3)
        !print *,rho(x,y,z), gam, s(x,y,z), 'xyz', x, y, z
        tt(x,y,z) = (rho(x,y,z)**gam*exp(ee(x,y,z)/rho(x,y,z)*(gam-1)))/rho(x,y,z)
        end do
        end do
        end do
      end if
    end if
  end if
  !---------------------------------------------------------------------
  ! Electron pressure
  !---------------------------------------------------------------------
  if (present(ne)) then
    error stop "ideal_gas_eos_mod: ne requested, but not coded"
  end if
  !---------------------------------------------------------------------
  ! Source function
  !---------------------------------------------------------------------
  if (present(src)) then
    error stop "ideal_gas_eos_mod: src requested, but not coded"
  end if
  !---------------------------------------------------------------------
  ! rho*kappa
  !---------------------------------------------------------------------
  if (present(rk)) then
    error stop "ideal_gas_eos_mod: rk requested, but not coded"
  end if
END SUBROUTINE lookup_stagger

!> =====================================================================
!> General EOS lookup routine, returning values in optional arguments
!> rho - density; ie - internal energy; pg - gas pressure; tt - gas temperature;
!> ne - electron pressure; rk - rho*kappa; src - source function;
!> gam - adiabatic index, csound - sound speed.
!> =====================================================================
SUBROUTINE lookup_zeus (self, rho, ie, pg, tt, ne, rk, src, gam, csound)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in) :: rho, ie
  real, dimension(:,:,:),   optional :: pg, tt, ne, rk
  real, dimension(:,:,:,:), optional :: src
  real, optional :: gam, csound
  !.....................................................................

  !---------------------------------------------------------------------
  ! Gas pressure
  !---------------------------------------------------------------------
  if (present(pg)) then
    if (gam==1.0) then
      pg = rho
    else
      pg = ie * (gam - 1.0)
    end if
  end if
  !---------------------------------------------------------------------
  ! Gas temperature
  !---------------------------------------------------------------------
  if (present(tt)) then
    if (present(pg)) then
      tt = pg/rho
    else
      if (gam==1.0) then
        tt = 1
      else
        tt = (ie * (gam - 1.0))/rho
      end if
    end if
  end if
  !---------------------------------------------------------------------
  ! Electron pressure
  !---------------------------------------------------------------------
  if (present(ne)) then
    error stop "ideal_gas_eos_mod: ne requested, but not coded"
  end if
  !---------------------------------------------------------------------
  ! Source function
  !---------------------------------------------------------------------
  if (present(src)) then
    error stop "ideal_gas_eos_mod: src requested, but not coded"
  end if
  !---------------------------------------------------------------------
  ! rho*kappa
  !---------------------------------------------------------------------
  if (present(rk)) then
    error stop "ideal_gas_eos_mod: rk requested, but not coded"
  end if
END SUBROUTINE lookup_zeus

!> =====================================================================
!> General EOS lookup routine, returning values in optional arguments
!> rho - density; s - entropy; pg - gas pressure; tt - gas temperature;
!> ne - electron pressure; rk - rho*kappa; src - source function;
!> gam - adiabatic index, csound - sound speed.
!> =====================================================================
SUBROUTINE lookup_ramses (self, rho, s, pg, tt, ne, rk, src, gam, csound)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in) :: rho, s
  real, dimension(:,:,:),   optional :: pg, tt, ne, rk
  real, dimension(:,:,:,:), optional :: src
  real, optional :: gam, csound
  !.....................................................................

  !---------------------------------------------------------------------
  ! Gas pressure
  !---------------------------------------------------------------------
  if (present(pg)) then
    if (gam==1.0) then
      pg = rho*csound**2
    else
      pg = rho**gam*exp(s/rho*(gam-1))
    end if
  end if
  !---------------------------------------------------------------------
  ! Gas temperature
  !---------------------------------------------------------------------
  if (present(tt)) then
    if (present(pg)) then
      tt = pg/rho
    else
      if (gam==1.0) then
        tt = csound**2
      else
        tt = (rho**gam*exp(s/rho*(gam-1)))/rho
      end if
    end if
  end if
  !---------------------------------------------------------------------
  ! Electron pressure
  !---------------------------------------------------------------------
  if (present(ne)) then
    error stop "ideal_gas_eos_mod: ne requested, but not coded"
  end if
  !---------------------------------------------------------------------
  ! Source function
  !---------------------------------------------------------------------
  if (present(src)) then
    error stop "ideal_gas_eos_mod: src requested, but not coded"
  end if
  !---------------------------------------------------------------------
  ! rho*kappa
  !---------------------------------------------------------------------
  if (present(rk)) then
    error stop "ideal_gas_eos_mod: rk requested, but not coded"
  end if
END SUBROUTINE lookup_ramses

END MODULE eos_mod
