!> =============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!> =============================================================================
MODULE eos_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE eos_reader_mod
  USE opacity_reader_mod
  USE lagrange_mod
  implicit none
  private
  type, public:: eos_t
    integer:: nd, ne, nv
    real(8):: lnd_min, lnd_max, lnd_stp
    real(8):: lne_min, lne_max, lne_stp
    real(8), dimension(:), pointer:: lnd, lne
    real(8), dimension(:,:,:), pointer:: table
  contains
    procedure:: init
    procedure:: pressure
    procedure:: lookup
  end type
  type (eos_t), public:: eos
CONTAINS

!> =============================================================================
!> The eos reader is supposed to deliver table date in the eos_reader data
!> type, with pressure and internal energy in an (nT,nd,2) table, with NT
!> temperature value T(nT), and nd density values d(nd).
!> =============================================================================
SUBROUTINE init (self)
  class (eos_t):: self
  real(8), dimension(:), pointer:: x, y
  real(8):: Emin, Emax, E, T, P
  integer:: i, j, nd, ne
  associate (reader=>eos_reader)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  call reader%init
  nd = reader%nd
  ne = reader%nT
  self%nd = nd
  self%ne = ne
  allocate (self%lnd(nd), self%lne(ne), self%table(ne,nd,2))
  !-----------------------------------------------------------------------------
  ! Fint extremum values of internal energy
  !-----------------------------------------------------------------------------
  Emin = reader%table(1,1,2)
  Emax = reader%table(1,1,2)
  do j=1,nd
    do i=1,ne
      Emin = min(Emin,reader%table(i,j,2))
      Emax = max(Emax,reader%table(i,j,2))
    end do
    !print '(a,1p,2e15.6)', 'Emin, Emax =', minval(reader%table(:,j,2)), maxval(reader%table(:,j,2))
  end do
  print '(a,1p,2e15.6)', 'Emin, Emax =', Emin, Emax
  !-----------------------------------------------------------------------------
  ! Construct table x- and y-values
  !-----------------------------------------------------------------------------
  self%lne_min = log(Emin)
  self%lne_max = log(Emax)
  self%lne_stp = (self%lne_max-self%lne_min)/(ne-1.0)
  do i=1,ne
    self%lne(i) = self%lne_min + (i-1)*self%lne_stp
  end do
  self%lnd = log(reader%d)
  self%lnd_min = self%lnd(1)
  self%lnd_max = self%lnd(nd)
  self%lnd_stp = (self%lnd_max-self%lnd_min)/(nd-1.0)
  do j=1,nd
    self%lnd(j) = self%lnd_min + (j-1)*self%lnd_stp
  end do
  !-----------------------------------------------------------------------------
  ! For each density, interpolate in the table to find the T corrsponding to
  ! the new table E, and then again, to find the pressure at that T
  !-----------------------------------------------------------------------------
  do i=1,ne
    do j=1,nd
      E = exp(self%lne(i))              ! new table E-value
      x => reader%table(:,j,2)          ! x-values = known E-values
      y => reader%T                     ! y-values = known T-values
      T = interpolate (x, y, E)         ! find the corresponding T
      x => reader%T                     ! x-values = table T-values
      y => reader%table(:,j,1)          ! y-values = table P-values
      P = interpolate (x, y, T)         ! find the pressure at T
      self%table(i,j,1) = P             ! store the pressure
      self%table(i,j,2) = T             ! store the temperature
      !print '(2i5,1p,3e15.6)', i, j, T, P, E
    end do
  end do
  open (io_unit%data, file='tmp.dat', form='formatted', status='unknown')
  write (io_unit%data,*) ne
  write (io_unit%data,*) self%lne
  write (io_unit%data,*) nd
  write (io_unit%data,*) self%lnd
  do j=1,nd
    write (io_unit%data,*) j, exp(self%lnd(j))
    do i=1,ne
      write (io_unit%data,'(1p,2e15.6)') self%table(i,j,:)
    end do
  end do
  close (io_unit%data)
  call trace%end
  end associate
END SUBROUTINE init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION interpolate (x, y, x0) RESULT (out)
  !USE lagrange_mod
  real(8), dimension(:):: x, y
  real(8):: x0, out, p
  integer:: i
  !-----------------------------------------------------------------------------
  out = lagrange%sequence (x0, x, y)
END FUNCTION interpolate

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION pressure (self, lnd, lne) RESULT (pg)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: lnd, lne
  real, dimension(size(lnd,1),size(lnd,2),size(lnd,3)):: pg
  !-----------------------------------------------------------------------------
  call self%lookup ([size(lnd,1),size(lnd,2),size(lnd,3)], lnd, lne, pg=pg)
END FUNCTION pressure

!> =============================================================================
!> General EOS lookup routine, returning values in optional arguments
!> =============================================================================
SUBROUTINE lookup (self, dim, lnd, lne, d, e, pg, tt, ne, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnd, lne
  real, dimension(:,:,:),   optional :: d, e, pg, tt, ne, rk
  real, dimension(:,:,:,:), optional :: src
  real, optional:: gamma
  !.............................................................................
  real, dimension(:,:,:), pointer:: lnd_loc, lne_loc
  real, dimension(dim(1)):: pr, pe, qr, qe
  real, dimension(dim(1),4,size(self%table,3)):: f
  real, dimension(dim(1),4):: g
  integer, dimension(dim(1)):: ir, ie
  integer:: ix, iy, iz, iv
  integer:: md, me, mv
  real:: clnd, clne
  logical:: skip
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays
  !-----------------------------------------------------------------------------
  mv=  self%nv
  md = self%nd
  me = self%ne
  if (present(pg).or.present(tt).or.present(ne).or.present(rk).or.present(src)) then
    skip = .false.
  else
    skip = .true.
  end if
  if (present(lnd)) then
    lnd_loc => lnd
  else if (present(d)) then
    allocate (lnd_loc(dim(1),dim(2),dim(3)))
    lnd = log(d)
  else
    call mpi%abort ('eos_t: missing density')
  end if
  if (present(lne)) then
    lne_loc => lne
  else if (present(e)) then
    allocate (lne_loc(dim(1),dim(2),dim(3)))
    lne = log(e)
  else
    call mpi%abort ('eos_t: missing internal energy')
  end if
  !-----------------------------------------------------------------------------
  ! Loop over points
  !-----------------------------------------------------------------------------
  clnd = 1./self%lnd_stp
  clne = 1./self%lne_stp
  do iz=1,dim(3)
  do iy=1,dim(2)
    !-------------------------------------------------------------------------
    ! Addresses and weights, vectorized
    !-------------------------------------------------------------------------
    pr = 1.+(lnd_loc(:,iy,iz)-self%lnd_min)*clnd
    pe = 1.+(lne_loc(:,iy,iz)-self%lne_min)*clne
    ir = pr
    ie = pe
    ir = max(1,min(ir,md))
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
    do ix=1,dim(1)
    do iv=1,mv
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
  if (present(d)) deallocate (lnd_loc)
  if (present(e)) deallocate (lne_loc)
END SUBROUTINE lookup

END MODULE eos_mod

