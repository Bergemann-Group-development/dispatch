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
  integer, parameter:: kind1=4, kind2=8, kind3=4, nf=10
  type, public:: eos_t
    integer:: nd, ne
    integer:: nv      = 2
    integer:: i_pg    = 1
    integer:: i_tt    = 2
    integer:: i_ee    = 2
    real(8):: gamma   = 1.4d0  ! required for syntax in some solvers
    real(kind=kind2):: lnd_min, lnd_max, lnd_stp
    real(kind=kind2):: lne_min, lne_max, lne_stp
    real(kind=kind2), dimension(:), pointer:: lne
    real(kind=kind2), dimension(:), pointer:: lnd
    real(kind=kind1), dimension(:,:,:), allocatable:: table
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
  integer,optional :: nv
  !.............................................................................
  real(8), dimension(:), pointer:: x, y
  real(kind=kind1), dimension(:), pointer:: x1, y1
  real(kind=kind2), dimension(:), pointer:: lnd
  real(kind=kind2), dimension(:,:,:), allocatable, target:: table
  real(8):: Pmin, Pmax, E, T, P, log_T, log_P, log_E
  integer:: i, j, iv, nd, ne
  associate (reader=>eos_reader)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  call scaling%init
  call reader%init
  call timer%tic
  if (io%master) print *, &
 '------------------------ constructing fast EOS table -------------------------'
  if (present(nv)) then
    self%nv = nv
  else
    self%nv = reader%nv
  end if
  nd = reader%nd
  ne = reader%nT*nf
  self%nd = nd
  self%ne = ne
  allocate (self%lne(ne))
  allocate (lnd(nd))
  allocate (table(ne,nd,self%nv))
  print *, 'table size:', storage_size(table)*(product(shape(table))/(8*1024.**2)),' MB'
  !-----------------------------------------------------------------------------
  ! Construct table x- and y-values
  !-----------------------------------------------------------------------------
  self%lne_min = maxval(reader%table(        1,:,2))
  self%lne_max = minval(reader%table(reader%nT,:,2))
  self%lne_stp = (self%lne_max-self%lne_min)/(ne-1.0)
  do i=1,ne
    self%lne(i) = self%lne_min + (i-1)*self%lne_stp
  end do
  self%lnd_min = reader%log_d(1)
  self%lnd_max = reader%log_d(reader%nd)
  self%lnd_stp = (self%lnd_max-self%lnd_min)/(nd-1.0)
  do i=1,nd
    lnd(i) = self%lnd_min + (i-1)*self%lnd_stp
  end do
  !-----------------------------------------------------------------------------
  ! For each density, interpolate in the table to find the T corrsponding to
  ! the new table E, and then again, to find the pressure at that T
  !-----------------------------------------------------------------------------
  do i=1,ne
    do j=1,nd
      log_E = self%lne_min + (i-1)*self%lne_stp ! table lnE value
      x => reader%table(:,j,2)                  ! x-values = known log_E-values
      y => reader%log_T                         ! y-values = known T-values
      log_T = interpolate (x, y, log_E)         ! find the corresponding log_T
      x => reader%log_T                         ! x-values = table log_T-values
      y => reader%table(:,j,1)                  ! y-values = table log_P-values
      log_P = interpolate (x, y, log_T)         ! find the log pressure at log_T
      table(i,j,1) = log_P                      ! store the log pressure
      table(i,j,2) = log_T                      ! store the log temperature
      if (io%verbose > 1 .and. io%master) then
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
  ! Interpolate in reader%table, to recover known values
  !-----------------------------------------------------------------------------
  print *, 'Interpolation in reader%table:'
  j = 61
  x => reader%log_T
  y => reader%table(:,j,self%i_ee)
  log_T = log(1e3)
  log_E = interpolate (x, y, log_T)
  y => reader%table(:,j,self%i_pg)
  log_P = interpolate (x, y, log_T)
  print *, 'rho:', exp(self%lnd_min+(j-1)*self%lnd_stp)/scaling%d, exp(reader%log_d(j))/scaling%d
  print *, '  T:', exp(log_T)
  print *, ' ee:', exp(log_E)/scaling%u**2
  print *, ' Pg:', exp(log_P)/scaling%p
  !-----------------------------------------------------------------------------
  ! Interpolate in table, to recover the values above
  !-----------------------------------------------------------------------------
  print *, 'Interpolation in table:'
  x => self%lne
  y => table(:,j,self%i_tt)
  !print *, 'x1'
  !print '(1p,10e12.3)', x1
  !print *, 'y1'
  !print '(1p,10e12.3)', y1
  T = exp(interpolate (x, y, log_E))
  y => table(:,j,self%i_pg)
  P = exp(interpolate (x, y, log_E))
  print *, 'lnE:', log_E
  print *, '  T:', T
  print *, ' ee:', exp(log_E)/scaling%u**2
  print *, ' Pg:', P/scaling%p
  !-----------------------------------------------------------------------------
  ! Make the table 10 times better resolved in log(d)
  !-----------------------------------------------------------------------------
  self%nd = (nd-1)*nf+1
  self%lnd_stp = self%lnd_stp/real(nf,kind=8)
  allocate (self%table(self%ne,self%nd,self%nv))
  allocate (self%lnd(self%nd))
  do j=1,self%nd
    self%lnd(j) = self%lnd_min + (j-1)*self%lnd_stp
  end do
  do iv=1,self%nv
  do i=1,self%ne
    x => lnd
    y => table(i,:,iv)
    do j=1,self%nd
      self%table(i,j,iv) = interpolate (x, y, self%lnd(j))
      self%table(i,j,iv) = exp(self%table(i,j,iv))
    end do
  end do
  end do
  !-----------------------------------------------------------------------------
  ! Rescaling table to code units
  !-----------------------------------------------------------------------------
  if (io%verbose > 0 .and. io%master) then
    print 1, 'Tmin, Tmax [K]    =', minval(self%table(:,:,self%i_tt)), &
                                    maxval(self%table(:,:,self%i_tt))
    1 format(1x,a,1p,2e15.6)
    Pmin = minval(self%table(:,:,self%i_pg))
    Pmax = maxval(self%table(:,:,self%i_pg))
    print 1, 'Pmin, Pmax [cgs]  =', Pmin, Pmax
    print 1, 'Emin, Emax [cgs]  =', exp(self%lne_min), exp(self%lne_max)
  end if
  self%lnd_min = self%lnd_min - log(scaling%d)
  self%lnd_max = self%lnd_max - log(scaling%d)
  self%lnd     = self%lnd     - log(scaling%d)
  self%lne_min = self%lne_min - log(scaling%u**2)
  self%lne_max = self%lne_max - log(scaling%u**2)
  self%lne     = self%lne     - log(scaling%u**2)
  self%table(:,:,self%i_pg) = self%table(:,:,self%i_pg)/scaling%p
  if (io%verbose > 0 .and. io%master) then
    Pmin = Pmin/scaling%p
    Pmax = Pmax/scaling%p
    print 1, 'Pmin, Pmax [code] =', Pmin, Pmax
    print 1, 'Emin, Emax [code] =', exp(self%lne_min), exp(self%lne_max)
  end if
  !-----------------------------------------------------------------------------
  ! Deallocate the reader table, which we don't need anymore
  !-----------------------------------------------------------------------------
  deallocate (eos_reader%table, table)
  call timer%toc ('eos_t%init')
  call trace%end()
  end associate
END SUBROUTINE init

!===============================================================================
!> Accurate interpolation, used only to generate the working table
!===============================================================================
FUNCTION interpolate (x, y, x0) RESULT (out)
  real(8), dimension(:):: x, y
  real(8):: x0, out, p
  integer:: i
  !-----------------------------------------------------------------------------
  out = lagrange%sequence (x0, x, y)
END FUNCTION interpolate

!===============================================================================
!> Accurate interpolation, used only to generate the working table
!===============================================================================
FUNCTION interpolate4 (x, y, x0) RESULT (out)
  real, dimension(:):: x, y
  real:: x0, out, p
  integer:: i
  !-----------------------------------------------------------------------------
  out = lagrange%sequence (x0, x, y)
END FUNCTION interpolate4

!===============================================================================
!> By convention, there should be a pressure (lnd, lne) function
!===============================================================================
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
SUBROUTINE lookup (self, dim, lnd, lne, d, ee, pg, tt, iskip)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: lnd, lne
  real, dimension(:,:,:), intent(in), optional :: d, ee
  real, dimension(:,:,:), optional :: pg, tt
  integer, optional:: iskip
  !.............................................................................
  real(kind=kind1), dimension(dim(1)):: pr, pe, qr, qe
  real(kind=kind2), dimension(dim(1)):: lnd1, lne1
  real(kind=kind1), dimension(dim(1),4,size(self%table,3)):: f
  real(kind=kind1), dimension(dim(1),4):: g
  integer, dimension(dim(1)):: ir, ie
  integer:: ix, iy, iz, iv
  integer:: md, me, mv, lskip
  real(kind=kind2):: clnd, clne
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays
  !-----------------------------------------------------------------------------
  mv = self%nv
  md = self%nd
  me = self%ne
  if (present(iskip)) then
    lskip = iskip
  else
    lskip = 0
  end if
  !-----------------------------------------------------------------------------
  ! Loop over points
  !-----------------------------------------------------------------------------
  clnd = 1./self%lnd_stp
  clne = 1./self%lne_stp
  do iz=1,dim(3)
  do iy=1,dim(2)
    !-------------------------------------------------------------------------
    ! Pick up one-dimensional values via vectorized assigments
    !-------------------------------------------------------------------------
    if (present(lne)) then
      lne1 = lne(:,iy,iz)
    else if (present(ee)) then
      lne1 = log(real(ee(:,iy,iz),kind=kind3))
    else
      call mpi%abort ('eos_t: missing internal energy')
    end if
    if (lskip==1) cycle
    if (present(lnd)) then
      lnd1 = lnd(:,iy,iz)
    else if (present(d)) then
      lnd1 = log(real(d(:,iy,iz),kind=kind3))
    else
      call mpi%abort ('eos_t: missing density')
    end if
    if (lskip==2) cycle
    !-------------------------------------------------------------------------
    ! Addresses and weights, vectorized
    !-------------------------------------------------------------------------
    pr = 1.+(lnd1(:)-self%lnd_min)*clnd
    pe = 1.+(lne1(:)-self%lne_min)*clne
    ir = pr
    ie = pe
    ir = max(1,min(ir,md))
    ie = max(1,min(ie,me))
    pr = pr-ir
    pe = pe-ie
    qr = 1.0-pr
    qe = 1.0-pe
    g(:,1) = qe*qr
    g(:,2) = pe*qr
    g(:,3) = qe*pr
    g(:,4) = pe*pr
    if (lskip==3) cycle
    !---------------------------------------------------------------------------
    !  Table values, indirect addressing
    !---------------------------------------------------------------------------
    do ix=1,dim(1)
    do iv=1,mv
      f(ix,1,iv) = self%table(ie(ix)  ,ir(ix)  ,iv)
      f(ix,2,iv) = self%table(ie(ix)+1,ir(ix)  ,iv)
      f(ix,3,iv) = self%table(ie(ix)  ,ir(ix)+1,iv)
      f(ix,4,iv) = self%table(ie(ix)+1,ir(ix)+1,iv)
    end do
    end do
    if (lskip==4) cycle
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
    if (lskip==5) cycle
    if (present(tt)) then
      tt(:,iy,iz) = &
         f(:,1,self%i_tt)*g(:,1) &
       + f(:,2,self%i_tt)*g(:,2) &
       + f(:,3,self%i_tt)*g(:,3) &
       + f(:,4,self%i_tt)*g(:,4)
    end if
    !---------------------------------------------------------------------------
  end do
  end do
END SUBROUTINE lookup

END MODULE eos_mod
