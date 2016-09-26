!> =============================================================================
!> EOS module, reading tabular EOS compatible with STAGGER code
!> =============================================================================
MODULE eos_mod
  USE io_mod
  USE io_unit_mod
  USE mpi_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Table related parameters
  !-----------------------------------------------------------------------------
  integer, parameter:: mr=101, me=101, mv=8, mbox=4
  real, parameter:: lnrmin=-6.0, lnrmax=6.0, eemin=1.0, eemax=10.
  real, parameter:: dlnr=(lnrmax-lnrmin)/(mr-1), dee=(eemax-eemin)/(me-1)
  real, parameter:: clnr=1./dlnr, cee=1./dee
  integer, allocatable, dimension(:,:,:):: table
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Public part
  !-----------------------------------------------------------------------------
  type, public:: eos_t
    logical:: on = .false.
  contains
    procedure:: init
    procedure:: lookup
  end type
  type (eos_t), public:: eos
CONTAINS

!> =============================================================================
!> Initialize EOS module (only once)
!> =============================================================================
SUBROUTINE init (self, mx ,my ,mz)
  class (eos_t):: self
  integer:: mx, my, mz, ie, ir, iv
  !-----------------------------------------------------------------------------
  !$omp critical (eos_cr)
  if (first_time) then
    first_time = .false.
    allocate (table(mr,me,mv))
    do ir=1,mr
    do ie=1,me
    do iv=1,mv
      table(ir,ie,iv) = lnrmin+dlnr*(ir-1) + eemin+dee*(ie-1) + iv
    end do
    end do
    end do
  end if
  !$omp end critical (eos_cr)
  !$omp barrier
END SUBROUTINE init

!> =============================================================================
!> General EOS lookup routine, returning values in optional arguments
!> =============================================================================
SUBROUTINE lookup (self, lnr, ee, pg, tt, ne, rk, src)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in) :: lnr, ee
  real, dimension(:,:,:),   optional :: pg, tt, ne, rk
  real, dimension(:,:,:,:), optional :: src
  !.............................................................................
  real, dimension(size(lnr,1)):: pr, pe, qr, qe
  real, dimension(size(lnr,1),4,mv):: f
  real, dimension(size(lnr,1),4):: g
  integer, dimension(size(lnr,1)):: ir, ie
  integer:: ix, iy, iz, iv
  integer:: mx, my, mz
  logical:: skip
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays
  !-----------------------------------------------------------------------------
  mx=size(lnr,1)
  my=size(lnr,2)
  mz=size(lnr,3)
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
    pr = 1.+(lnr(:,iy,iz)-lnrmin)*clnr
    pe = 1.+( ee(:,iy,iz)-eemin )*cee
    ir = pr
    ie = pe
    ir = max(1,min(ir,mr))
    ie = max(1,min(ie,me))
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
    do iv=1,mv
      f(ix,1,iv) = table(ir(ix)  ,ie(ix)  ,iv)
      f(ix,2,iv) = table(ir(ix)+1,ie(ix)  ,iv)
      f(ix,3,iv) = table(ir(ix)  ,ie(ix)+1,iv)
      f(ix,4,iv) = table(ir(ix)+1,ie(ix)+1,iv)
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
END SUBROUTINE lookup

END MODULE eos_mod
