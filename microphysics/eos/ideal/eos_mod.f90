!===============================================================================
!> Equation of state module for any sort of tables, provided by a reader
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE kinds_mod
  USE scaling_mod
  USE units_mod
  implicit none
  private
  type, public:: eos_t
    real:: gamma = 1.4                  ! in case eos%init is never called
    integer :: eos_type = 0             ! 0 - ideal, x - any other
    integer :: n_lambda = 1
    real :: w_lambda = 1
    logical :: do_rt = .false.
    real :: kappa
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
  real, save:: gamma=1.4
  real, save:: kappa=0.1
  logical, save:: do_rt = .false.
  integer, save:: itimer=0
  namelist /eos_params/ gamma, do_rt, kappa
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init', itimer=itimer)
  if (io%master) print *, &
 '------------------------ ideal equation-of-state -----------------------------'
  !$omp critical (init_cr)
  if (first_time) then
    gamma = self%gamma
    do_rt = self%do_rt
    kappa = self%kappa
    rewind (io%input); read (io%input, eos_params, iostat=iostat)
    if (io%master) write (io%output, eos_params)
  end if
  !$omp end critical (init_cr)
  self%gamma = gamma
  self%do_rt = do_rt
  self%kappa = kappa * scaling%m/scaling%l**2
  call trace%end(itimer)
END SUBROUTINE init

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
FUNCTION pressure (self, d, e) RESULT (pg)
  class (eos_t):: self
  real(kind=KindScalarVar), dimension(:,:,:), pointer:: d, e
  real(kind=KindScalarVar), dimension(size(d,1),size(d,2),size(d,3)):: pg
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
SUBROUTINE lookup (self, dim, lnd, ee, lnx, x, lny, y, pg, tt, ss, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, optional:: gamma
  real, dimension(:,:,:), intent(in), pointer, optional :: lnx, x, lny, y, lnd, ee
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, dimension(:,:,:), pointer :: x_loc, y_loc
  integer:: i
  integer, save:: itimer=0
  real:: stefan
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  if (present(lnd)) then
    call mpi%abort('trying to use lnd in eos/ideal/eos_mod.f90')
  else if (present(ee)) then
    call mpi%abort('trying to use ee in eos/ideal/eos_mod.f90')
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
      rk(:,:,:,i) = y_loc*self%kappa
    end do
  end if
  if (present(lnx)) deallocate(x_loc)
  if (present(lny)) deallocate(y_loc)
  call trace%end(itimer)
END SUBROUTINE lookup

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table (self, dim, e, ee, d, lnd, lne,lnee, pg, tt, ne, rk, src, gamma)

  class (eos_t):: self
  integer:: dim(3)
  real, optional:: gamma
  real, dimension(:,:,:), intent(in), pointer, optional :: e, lne, d, lnee, lnd, ee
  real, dimension(:,:,:),                      optional :: pg, tt, ne
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, dimension(:,:,:), pointer :: d_loc, ee_loc
  integer:: i
  integer, save:: itimer=0
  real:: stefan
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  if (present(lnd)) then
    call mpi%abort('trying to use lnd in eos/ideal/eos_mod.f90')
  else if (present(ee)) then
    call mpi%abort('trying to use ee in eos/ideal/eos_mod.f90')
  end if
  !-----------------------------------------------------------------------------
  if (present(d)) then
    d_loc => d
  else if (present(lnd)) then
    allocate (d_loc(dim(1),dim(2),dim(3)))
    d_loc = exp(lnd)
  else
    call mpi%abort ('eos_t: missing 1st argument')
  end if
  !-----------------------------------------------------------------------------
  if (present(ee)) then
    ee_loc => ee
  else if (present(lne)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    ee_loc = exp(lne)/d_loc
  else if (present(e)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    ee_loc = e/d_loc
  else
    call mpi%abort ('eos_t: missing 2nd argument')
  end if
  !-----------------------------------------------------------------------------
  if (present(tt)) then
    tt = (self%gamma-1.0)*ee_loc
  end if
  if (present(src)) then
    stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
    do i=1,size(src,4)
      src(:,:,:,i) = stefan*real(((self%gamma-1.0)*ee_loc),kind=8)**4
    end do
  end if
  if (present(pg)) then
    pg = (self%gamma-1.0)*ee_loc*d_loc
  end if
  if (present(rk)) then
    do i=1,size(rk,4)
      rk(:,:,:,i) = d_loc*self%kappa
    end do
  end if
  if (present(lne).or.present(e).or.present(lnee)) deallocate(ee_loc)
  if (present(lnd)) deallocate(d_loc)
  call trace%end(itimer)
END SUBROUTINE lookup_table

END MODULE eos_mod

