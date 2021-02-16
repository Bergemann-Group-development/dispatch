!===============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE eos_reader_mod
  USE lagrange_mod
  USE scaling_mod
  USE timer_mod
  implicit none
  private
  type, public:: eos_t
    integer:: nd, ne
    integer:: nv      = 2
    integer:: i_pg    = 1
    integer:: i_tt    = 2
    integer:: i_kappa = 3
    integer:: i_src   = 4
    integer:: n_lambda = 1
    integer:: w_lambda = 1
    real:: lnd_min, lnd_max, lnd_stp
    real:: lne_min, lne_max, lne_stp
    real(8):: gamma=1.4d0    ! required for syntax in some solvers
    real, dimension(:), pointer:: lnd, lne
    real, dimension(:,:,:), pointer:: table
    logical :: do_rt
  contains
    procedure:: init
    procedure:: pressure
    procedure:: lookup
  end type
  type (eos_t), public:: eos
CONTAINS

!===============================================================================
!> The eos reader is supposed to deliver table date in the eos_reader data
!> type, with pressure and internal energy in an (nT,nd,2) table, with NT
!> temperature value T(nT), and nd density values d(nd).
!===============================================================================
SUBROUTINE init (self,nv)
  class (eos_t):: self
  real, dimension(:), pointer:: x, y
  real:: log_Emin, log_Emax, Pmin, Pmax, E, T, P, log_T, log_P, log_E
  integer:: i, j, nd, ne
  integer,optional :: nv
  associate (reader=>eos_reader)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  call reader%init
  call timer%tic
  if (io%master) print *, &
 '------------------------ constructing fast EOS table -------------------------'
  if (present(nv)) self%nv = nv
  nd = reader%nd
  ne = reader%nT
  self%nd = nd
  self%ne = ne
  allocate (self%lnd(nd), self%lne(ne), self%table(ne,nd,self%nv))
  !-----------------------------------------------------------------------------
  ! Fint extremum values of internal energy
  !-----------------------------------------------------------------------------
  log_Emin = maxval(reader%table( 1,:,2))
  log_Emax = minval(reader%table(ne,:,2))
  !-----------------------------------------------------------------------------
  ! Construct table x- and y-values
  !-----------------------------------------------------------------------------
  self%lne_min = log_Emin*log(10.)
  self%lne_max = log_Emax*log(10.)
  self%lne_stp = (self%lne_max-self%lne_min)/(ne-1.0)
  do i=1,ne
    self%lne(i) = self%lne_min + (i-1)*self%lne_stp
  end do
  self%lnd = reader%log_d*log(10.)
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
      E = exp(self%lne(i))              ! new table lnE-value
      log_E = log10(E)
      x => reader%table(:,j,2)          ! x-values = known log_E-values
      y => reader%log_T                 ! y-values = known T-values
      log_T = interpolate (x, y, log_E) ! find the corresponding log_T
      x => reader%log_T                 ! x-values = table log_T-values
      y => reader%table(:,j,1)          ! y-values = table log_P-values
      log_P = interpolate (x, y, log_T) ! find the log pressure at log_T
      self%table(i,j,1) = 10.**log_P    ! store the linear pressure
      self%table(i,j,2) = 10.**log_T    ! store the linear temperature
      if (io%verbose > 0 .and. io%master) then
        if ((i==1.or.i==ne).and.(j==1.or.j==nd))then
          print '(2i5,3(2x,a,1p,2e12.4))', i, j, &
            ' log_E, E:', log_E, 10.**log_E, &
            ' log_T, T:', log_T, 10.**log_T, &
            ' log_P, P:', log_P, 10.**log_P
        end if
      end if
      !print '(2i5,1p,3e15.6)', i, j, T, P, E
    end do
  end do
  !-----------------------------------------------------------------------------
  ! Rescaling table to code units
  !-----------------------------------------------------------------------------
  if (io%verbose > 0 .and. io%master) then
    print 1, 'Emin, Emax [cgs]  =', 10.**log_Emin, 10.**log_Emax
    1 format(1x,a,1p,2e15.6)
    log_Emin = log_Emin - log10(scaling%u**2)
    log_Emax = log_Emax - log10(scaling%u**2)
    print 1, 'Emin, Emax [code] =', 10.**log_Emin, 10.**log_Emax
    Pmin = minval(self%table(:,:,1))
    Pmax = maxval(self%table(:,:,1))
    print 1, 'Pmin, Pmax [cgs]  =', Pmin, Pmax
    Pmin = Pmin/(scaling%d*scaling%u**2)
    Pmax = Pmax/(scaling%d*scaling%u**2)
    print 1, 'Pmin, Pmax [code] =', Pmin, Pmax
    print 1, 'Tmin, Tmax [K]    =', minval(self%table(:,:,2)), &
                                    maxval(self%table(:,:,2))
  end if
  self%lnd = self%lnd - log(scaling%d)
  self%lne = self%lne - log(scaling%u**2)
  self%table(:,:,1) = self%table(:,:,1)/(scaling%d*scaling%u**2)
  !-----------------------------------------------------------------------------
  ! Deallocate the reader table, which we don't need anymore
  !-----------------------------------------------------------------------------
  deallocate (eos_reader%table)
  call timer%toc ('eos_t%init')
  call trace%end()
  end associate
END SUBROUTINE init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION interpolate (x, y, x0) RESULT (out)
  !USE lagrange_mod
  real, dimension(:):: x, y
  real:: x0, out, p
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

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup (self, dim, lnd, lne, d, ee, pg, tt, ne, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnd, lne
  real, dimension(:,:,:),   optional :: d, ee, pg, tt, ne
  real, dimension(:,:,:,:), optional :: src, rk
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
    lnd_loc = log(d)
  else
    call mpi%abort ('eos_t: missing density')
  end if
  if (present(lne)) then
    lne_loc => lne
  else if (present(ee)) then
    allocate (lne_loc(dim(1),dim(2),dim(3)))
    lne_loc = log(ee)
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
         f(:,1,self%i_pg)*g(:,1) &
       + f(:,2,self%i_pg)*g(:,2) &
       + f(:,3,self%i_pg)*g(:,3) &
       + f(:,4,self%i_pg)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(tt)) then
      tt(:,iy,iz) = &
         f(:,1,self%i_tt)*g(:,1) &
       + f(:,2,self%i_tt)*g(:,2) &
       + f(:,3,self%i_tt)*g(:,3) &
       + f(:,4,self%i_tt)*g(:,4)
    end if
    !---------------------------------------------------------------------------
    if (present(rk)) then
      do iv=0,size(rk,4)-1
        rk(:,iy,iz, iv+1) = &
           f(:,1,self%i_kappa+iv)*g(:,1) &
         + f(:,2,self%i_kappa+iv)*g(:,2) &
         + f(:,3,self%i_kappa+iv)*g(:,3) &
         + f(:,4,self%i_kappa+iv)*g(:,4)
      end do
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
     do iv=0,size(src,4)-1
      src(:,iy,iz,iv+1) = &
         f(:,1,self%i_src+iv)*g(:,1) &
       + f(:,2,self%i_src+iv)*g(:,2) &
       + f(:,3,self%i_src+iv)*g(:,3) &
       + f(:,4,self%i_src+iv)*g(:,4)
     end do
    end if
  end do
  end do
  if (present(d))  deallocate (lnd_loc)
  if (present(ee)) deallocate (lne_loc)
END SUBROUTINE lookup

END MODULE eos_mod
