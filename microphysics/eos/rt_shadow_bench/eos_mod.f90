!===============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!=======================================================================
!> Simple opacity, expressed as
!> kappa = (T/T0)^-3.5*(rho/rho0)^2
!> T0 = 290K, rho0 = 1 g cm-3
!=======================================================================
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE scaling_mod
  USE units_mod
  USE kinds_mod
  implicit none
  private
  type, public:: eos_t
    real(8):: gamma = 1.6666667                  ! in case eos%init is never called
    integer :: eos_type = 0             ! 0 - ideal, x - any other
    real :: t0    = 290.
    real :: rho0 = 1.
    real :: alfa0 = 1.
    real :: ee_min = 0.
    integer :: n_lambda = 1
    real :: w_lambda = 1.
    logical :: do_rt = .false.
  contains
    procedure:: init
    procedure:: pressure
    procedure:: temperature
    procedure:: lookup
    procedure:: lookup_table
  end type
  type (eos_t), public:: eos
CONTAINS

!===============================================================================
!> The eos reader is supposed to deliver table date in the eos_reader data
!> type, with pressure and internal energy in an (nT,nd,2) table, with NT
!> temperature value T(nT), and nd density values d(nd).
!===============================================================================
SUBROUTINE init (self)
  class (eos_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  real(8), save:: gamma=1.66666667
  logical, save:: do_rt = .false.
  integer, save:: itimer=0
  namelist /eos_params/ gamma, do_rt
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init', itimer=itimer)
  call io%header ('microphysics/eos/ideal/eos_mod:')
  call io%lock%set ('eos_t%init')
  if (first_time) then
    gamma = self%gamma
    do_rt = self%do_rt
    rewind (io%input)
    read (io%input, eos_params, iostat=iostat)
    call io%namelist_warning ('eos_params', iostat)
    if (io%master) write (io%namelist, eos_params)
    self%gamma = gamma
    self%t0   = 1./(290./scaling%temp)
    self%rho0 = 1./(1.  /scaling%d)
    self%alfa0 = 1. * scaling%m/scaling%l**2
    first_time = .false.
    self%do_rt = do_rt
  end if
  call io%lock%unset ('eos_t%init')
  call trace%end(itimer)
END SUBROUTINE init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION pressure (self, d, e) RESULT (pg)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: d, e
  real, dimension(size(d,1),size(d,2),size(d,3)):: pg
  !-----------------------------------------------------------------------------
  pg = (self%gamma-1.0)*d*e
END FUNCTION

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION temperature (self, d, e) RESULT (tt)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: d, e
  real, dimension(size(d,1),size(d,2),size(d,3)):: tt
  !-----------------------------------------------------------------------------
  tt = (self%gamma-1.0)*e
END FUNCTION

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup (self, dim, lnx, x, lny, y, pg, tt, ss, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, optional:: gamma
  real, dimension(:,:,:), intent(in), pointer, optional :: lnx, x, lny, y
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, dimension(:,:,:), pointer :: x_loc, y_loc
  real:: stefan
  integer:: i
  integer, save:: itimer = 0
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  !-----------------------------------------------------------------------------
  if (present(y)) then
    y_loc => y
  else if (present(lny)) then
    allocate (y_loc(dim(1),dim(2),dim(3)))
    y_loc = exp(lny)
  else
    call mpi%abort ('eos_t: missing 2nd argument')
  end if
  !-----------------------------------------------------------------------------
  if (present(x)) then
    x_loc => x
  else if (present(lnx)) then
    allocate (x_loc(dim(1),dim(2),dim(3)))
    x_loc = exp(lnx)
  else
    call mpi%abort ('eos_t: missing 1st argument')
  end if
  !-----------------------------------------------------------------------------
  if (present(tt)) then
    tt = (self%gamma-1.0)*x_loc
  end if
  if (present(src)) then
    stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
    do i=1,size(src,4)
      src(:,:,:,i) = stefan*real(((self%gamma-1.0)*x_loc),kind=8)**4
    end do
  end if
  if (present(pg)) then
    pg = (self%gamma-1.0)*y_loc*x_loc
  end if
  if (present(rk)) then
    do i=1,size(rk,4)
      rk(:,:,:,i) = self%alfa0 *((self%gamma-1.0)*x_loc * self%t0)**(-3.5)*(y_loc*self%rho0)**2
    end do
  end if
  if (present(lnx)) deallocate(x_loc)
  if (present(lny)) deallocate(y_loc)
  call trace%end(itimer)
END SUBROUTINE lookup

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table (self, dim, e, ee, d, lnd, lne,lnee, pg, tt, ne, rk, &
                         ss, s, lns, src, gamma,e_below, label, id)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:)  , pointer, optional :: lnd,d,lne,lnee,e,ee,ne
  real, dimension(:,:,:)  , pointer, optional :: ss, s, lns
  real, dimension(:,:,:,:)         , optional :: src, rk
  real, dimension(:,:,:)           , optional :: pg, tt
  integer, optional:: id
  real, optional:: gamma
  logical, optional:: e_below
  character(len=*),optional::label
  ! --------------------------------------------------------------------
  real, dimension(:,:,:), pointer    :: d_loc, e_loc
  integer,         dimension(dim(1)) :: il, jl
  real(kind=4),    dimension(dim(1)) :: f1, f2, f3, f4
  real    :: dmin,dmax,dd, emin,emax,de,fr, m1t, m1u, t, u, stefan
  integer :: i
  integer :: mx,my,mz
  integer, save:: itimer=0
  ! ....................................................................
  call trace%begin ('eos_t%lookup_table', itimer=itimer)
  if (io%verbose>2) then
    print *, 'present(lnd) =', present(lnd)
    print *, 'present(ee)  =', present(ee)
    print *, 'present(pg)  =', present(pg)
    print *, 'present(src) =', present(src)
    print *, 'present(rk)  =', present(rk)
    print *, 'present(tt)  =', present(tt)
    print *, 'present(ne)  =', present(ne)
    print *, 'present(gamma)  =', present(gamma)
    print *, 'dim:', dim
  end if
  mx=dim(1)
  my=dim(2)
  mz=dim(3)
  ! --------------------------------------------------------------------
  ! Allocate temporary arrays, if needed
  ! --------------------------------------------------------------------
  ! For log(density)
  if (present(d)) then
    d_loc => d
  else if (present(lnd)) then
    allocate (d_loc(dim(1),dim(2),dim(3)))
    d_loc = exp(lnd)
  else
    call io%abort ('eos_t: missing density')
  end if
  ! For log(internal energy)
  if (present(lnee)) then
    allocate (e_loc(dim(1),dim(2),dim(3)))
    e_loc = exp(lnee)*d_loc
  else if (present(e)) then
    e_loc => e
  else if (present(lne)) then
    allocate (e_loc(dim(1),dim(2),dim(3)))
    e_loc = exp(lne)
  else if (present(ee)) then
    allocate (e_loc(dim(1),dim(2),dim(3)))
    e_loc = ee*d_loc
  else
    call io%abort ('eos_t: missing internal energy')
  end if
  !
  if (present(tt)) then
    tt = (self%gamma-1.0)*e_loc
  end if
  if (present(src)) then
    stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
    do i=1,size(src,4)
      src(:,:,:,i) = stefan*real(((self%gamma-1.0)*e_loc),kind=8)**4
    end do
  end if
  if (present(pg)) then
    pg = (self%gamma-1.0)*d_loc*e_loc
  end if
  if (present(rk)) then
    do i=1,size(rk,4)
      rk(:,:,:,i) = self%alfa0 *((self%gamma-1.0)*e_loc * self%t0)**(-3.5)*(d_loc*self%rho0)**2
    end do
  end if
  !
  if (present(lnd)) then
    deallocate (d_loc)
  end if
  if (present(lnee).or.present(lne).or.present(ee)) then
    deallocate (e_loc)
  end if
  call trace%end (itimer)
END SUBROUTINE lookup_table

END MODULE eos_mod

