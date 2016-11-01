!> =====================================================================
!> EOS module, reading tabular EOS compatible with STAGGER code
!> =====================================================================
MODULE eos_mod
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  implicit none
  private
  !---------------------------------------------------------------------
  ! Table related parameters
  !---------------------------------------------------------------------
  ! real(kind=4):: eosxmin,ur,ul,ut,eps,tff,grv,abnd,dbox  
  !mr=101, ms=101, mv=8, mbox=4
  !real, parameter:: lnrmin=-6.0, lnrmax=6.0, eemin=1.0, eemax=10.
  !real, parameter:: dlnr=(lnrmax-lnrmin)/(mr-1), dee=(smax-smin)/(ms-1)
  !real, parameter:: clnr=1./dlnr, cs=1./ds
  
  logical, save:: first_time=.true.
  !---------------------------------------------------------------------
  ! Public part
  !---------------------------------------------------------------------
  type, public:: eos_t
    integer :: mbox = 4                                                 ! # of opacity bins
    integer :: nvar = 8                                                 ! number of variables
    real(4) :: dbox = 10.                                               ! opacity bin spacing
    real(4) :: ul, ut, ur, eps, tff, grv, abnd                          ! some other basic parameters (??)
    integer :: nlnrho, nee                                              ! $ of ln(rho), ee entries
    real(4) :: lnrho_min, dlnrho                                        ! ln(rho) scale
    real(4) :: ee_min, dee                                              ! ee scale (energy or entropy)
    real(4) :: clnrho, cee                                              ! 1/dlnrho, 1/dee
    real(4), allocatable, dimension(:,:,:):: table                      ! table
    integer :: eos_type = 2
    logical:: on = .false.
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
  real(4) :: xmin !, dbox, ul,ut,ur,eos,tff,grv,abnd
  !real(4) :: lnrho_max, lnrho_min, dlnrho
  !real(4) :: ee_max, ee_min, dee
  !integer :: mbox, nlnrho, nee
  character(len=64):: tablefile
  logical:: do_eos
  namelist /eos_params/ do_eos, tablefile
  !---------------------------------------------------------------------
  !$omp critical (eos_cr)
  if (first_time) then
    first_time = .false.
    self%on = .true.
    tablefile = '../../data/eos/stagger_eos_square_table.dat'
    rewind (io%input); read (io%input,eos_params)
    if (io%master) write (*,eos_params)
    !-----------------------------------------------------------------------------
    ! Read the table file
    !-----------------------------------------------------------------------------
    open (io_unit%data, file=tablefile, form='unformatted', status='old', action='read')
    read (io_unit%data) self%nvar, xmin, self%dbox, self%ul, &
    self%ut, self%ur, self%eps, self%tff, self%grv, self%abnd
    print *,'mbox, xmin, dbox, ul, ut, ur, eps, tff, grv, abnd'
    print *, self%mbox, xmin, self%dbox, self%ul, &
    self%ut, self%ur, self%eps, self%tff, self%grv, self%abnd
    read (io_unit%data) self%nlnrho, self%lnrho_min, self%dlnrho
    read (io_unit%data) self%nee, self%ee_min, self%dee
    read (io_unit%data) self%nlnrho, self%lnrho_min, self%dlnrho

    read (io_unit%data) self%nee, self%ee_min, self%dee
    allocate (self%table(self%nlnrho,self%nee,self%nvar))
    read (io_unit%data) self%table
    self%clnrho=1./self%dlnrho
    self%cee=1./self%dee
    select case(solver)
    case('stagger ')
      self%lookup => lookup_stagger
    case('stagger2')
      self%lookup => lookup_stagger
    case('zeus    ')
      print *, "WARNING: stagger EOS might not be compatible with",solver," solver"
      self%lookup => lookup_zeus
    case('ramses  ')
      print *, "WARNING: stagger EOS might not be compatible with",solver," solver"
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
  !.....................................................................
  real, dimension(size(rho,1),size(rho,2),size(rho,3)):: lnr
  real, dimension(size(rho,1)):: pr, pe, qr, qe
  real, dimension(size(rho,1),4,self%nvar):: f
  real, dimension(size(rho,1),4):: g
  integer, dimension(size(rho,1)):: ir, ie
  integer:: ix, iy, iz, iv
  integer:: mx, my, mz
  logical:: skip
!-----------------------------------------------------------------------
  ! Allocate temporary arrays
  !---------------------------------------------------------------------
  lnr = log(rho)
  mx=size(rho,1)
  my=size(rho,2)
  mz=size(rho,3)
  if (present(pg).or.present(tt).or.present(ne).or.present(rk).or.present(src)) then
    skip = .false.
  else
    skip = .true.
  end if
  !-----------------------------------------------------------------------------
  ! Loop over points
  !-----------------------------------------------------------------------------
  do iz=1,mz
  do iy=1,my
    !-------------------------------------------------------------------------
    ! Addresses and weights, vectorized
    !-------------------------------------------------------------------------
    pr = 1.+(lnr(:,iy,iz)-self%lnrho_min)*self%clnrho
    pe = 1.+( ee(:,iy,iz)-self%ee_min )*self%cee
    ir = pr
    ie = pe
    ir = max(1,min(ir,self%nlnrho))
    ie = max(1,min(ie,self%nee))
    pr = pr-ir
    pe = pe-ie
    qr = 1.0-pr
    qe = 1.0-pe
    g(:,1) = qr*qe
    g(:,2) = pr*qe
    g(:,3) = qr*pe
    g(:,4) = pr*pe
    !---------------------------------------------------------------------------
    !  Table values, indirect addressing
    !---------------------------------------------------------------------------
    do ix=1,mx
    do iv=1,self%nvar
      f(ix,1,iv) = self%table(ir(ix)  ,ie(ix)  ,iv)
      f(ix,2,iv) = self%table(ir(ix)+1,ie(ix)  ,iv)
      f(ix,3,iv) = self%table(ir(ix)  ,ie(ix)+1,iv)
      f(ix,4,iv) = self%table(ir(ix)+1,ie(ix)+1,iv)
    end do
    end do
    if (skip) cycle
    !---------------------------------------------------------------------------
    !  Loop over table entries
    !---------------------------------------------------------------------------
    if (present(pg)) then
      pg(:,iy,iz) = &
         f(:,1,1)*g(:,1) &
       + f(:,2,1)*g(:,2) &
       + f(:,3,1)*g(:,3) &
       + f(:,4,1)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(rk)) then
      rk(:,iy,iz) = &
         f(:,1,2)*g(:,1) &
       + f(:,2,2)*g(:,2) &
       + f(:,3,2)*g(:,3) &
       + f(:,4,2)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(tt)) then
      tt(:,iy,iz) = &
         f(:,1,3)*g(:,1) &
       + f(:,2,3)*g(:,2) &
       + f(:,3,3)*g(:,3) &
       + f(:,4,3)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(ne)) then
      ne(:,iy,iz) = &
         f(:,1,4)*g(:,1) &
       + f(:,2,4)*g(:,2) &
       + f(:,3,4)*g(:,3) &
       + f(:,4,4)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(src)) then
     do iv=1,size(src,4)
      src(:,iy,iz,iv) = &
         f(:,1,4+iv)*g(:,1) &
       + f(:,2,4+iv)*g(:,2) &
       + f(:,3,4+iv)*g(:,3) &
       + f(:,4,4+iv)*g(:,4)
     end do
    end if
  end do
  end do
END SUBROUTINE lookup_stagger

!> =====================================================================
!> General EOS lookup routine, returning values in optional arguments
!> rho - density; s - internal energy; pg - gas pressure; tt - gas temperature;
!> ne - electron pressure; rk - rho*kappa; src - source function;
!> gam - adiabatic index, csound - sound speed.
!> =====================================================================
SUBROUTINE lookup_zeus (self, rho, s, pg, tt, ne, rk, src, gam, csound)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in) :: rho, s
  real, dimension(:,:,:),   optional :: pg, tt, ne, rk
  real, dimension(:,:,:,:), optional :: src
  real, optional :: gam, csound
  !.....................................................................
  logical:: skip
  integer:: mx,my,mz,mbox

  !---------------------------------------------------------------------
  ! Gas pressure
  !---------------------------------------------------------------------
  if (present(pg)) then
    if (gam==1.0) then
      pg = rho
    else
      pg = s * (gam - 1.0)
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
        tt = (s * (gam - 1.0))/rho
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
  logical:: skip
  integer:: mx,my,mz,mbox

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
