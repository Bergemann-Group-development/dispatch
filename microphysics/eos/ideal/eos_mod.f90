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
    real(8):: gamma = 1d0               ! in case eos%init is never called
    integer :: eos_type = 0             ! 0 - ideal, x - any other
    integer :: n_lambda = 1
    real :: w_lambda = 1
    real :: csound=1.0
    real :: mu=2.33
    real:: ee_min = 0.0
    logical :: do_rt = .false.
    real :: kappa
  contains
    procedure:: init
    procedure:: pressure
    procedure:: temperature
    procedure:: lookup
    procedure:: lookup_table
    procedure:: lookup_table_e
    procedure:: lookup_table_s
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
  real(8), save:: gamma=1d0
  real, save:: kappa=0.1, mu=2.33
  logical, save:: do_rt = .false.
  integer, save:: itimer=0
  namelist /eos_params/ gamma, mu, do_rt, kappa
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init', itimer=itimer)
  call io%header ('microphysics/eos/ideal/eos_mod:')
  call io%lock%set ('eos_t%init')
  if (first_time) then
    mu = self%mu
    gamma = self%gamma
    do_rt = self%do_rt
    kappa = self%kappa
    rewind (io%input)
    read (io%input, eos_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, eos_params)
    end if
    call io%namelist_warning ('eos_params', iostat)
    write (io%namelist, eos_params)
  end if
  self%mu = mu
  self%gamma = gamma
  self%do_rt = do_rt
  self%kappa = kappa * scaling%m/scaling%l**2
  call io%lock%unset ('eos_t%init')
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
SUBROUTINE lookup (self, dim, lnd, ee, lnx, x, lny, y, pg, tt, ss, rk, src, &
                   gamma, e_below)
  class (eos_t):: self
  integer:: dim(3)
  real, optional:: gamma
  logical, optional:: e_below
  real, dimension(:,:,:), intent(in), pointer, contiguous,optional :: lnx, x, lny, y, lnd, ee
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, dimension(:,:,:), pointer, contiguous :: x_loc, y_loc
  integer:: i
  integer, save:: itimer=0
  real:: stefan
  logical:: x_alloc, y_alloc
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  e_below = .false.
  if (present(gamma)) then
    self%gamma = gamma
  end if
  x_alloc = .false.
  y_alloc = .false.
  if (present(lnd)) then
    allocate (x_loc(dim(1),dim(2),dim(3)))
    x_alloc = .true.
    x_loc = exp(lnd)
  end if
  if (present(ee)) then
    y_loc => ee
  end if
  !-----------------------------------------------------------------------------
  if (present(y)) then
    y_loc => y
  else if (present(lny)) then
    allocate (y_loc(dim(1),dim(2),dim(3)))
    y_alloc = .true.
    y_loc = exp(lny)
  else
    call mpi%abort ('eos_t: missing 2nd argument')
  end if
  !-----------------------------------------------------------------------------
  if (present(x)) then
    x_loc => x
  else if (present(lnx)) then
    allocate (x_loc(dim(1),dim(2),dim(3)))
    x_alloc = .true.
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
  if (x_alloc) deallocate(x_loc)
  if (y_alloc) deallocate(y_loc)
  call trace%end(itimer)
END SUBROUTINE lookup

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table (self, dim, e, ee, d, lnd, lne,lnee, pg, tt, ne, rk, &
                         ss, s, lns, src, gamma,e_below, label, id)
  class (eos_t):: self
  integer:: dim(3)
  integer, optional:: id
  real, optional:: gamma
  real, dimension(:,:,:), intent(in), pointer, optional :: e, lne, d, lnee, lnd, ee
  real, dimension(:,:,:), intent(in), pointer, optional :: ss, s, lns
  real, dimension(:,:,:),                      optional :: pg, tt, ne
  real, dimension(:,:,:,:),                    optional :: src, rk
  logical, optional:: e_below
  character(len=*),optional::label
  real, dimension(:,:,:), pointer :: d_loc, ee_loc
  integer:: i
  integer, save:: itimer=0
  real:: stefan
  logical:: d_alloc, e_alloc
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  nullify (d_loc, ee_loc)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  d_alloc = .false.
  e_alloc = .false.
   !-----------------------------------------------------------------------------
  if (present(d)) then
    d_loc => d
  else if (present(lnd)) then
    allocate (d_loc(dim(1),dim(2),dim(3)))
    d_alloc = .true.
    d_loc = exp(lnd)
  else
    call mpi%abort ('eos_t: missing 1st argument')
  end if
  !-----------------------------------------------------------------------------
  if (present(ee) .or. present(e) .or. present(lne)) then
    if (present(ee)) then
      ee_loc => ee
    else if (present(lne)) then
      allocate (ee_loc(dim(1),dim(2),dim(3)))
      e_alloc = .true.
      ee_loc = exp(lne)/d_loc
    else if (present(e)) then
      allocate (ee_loc(dim(1),dim(2),dim(3)))
      e_alloc = .true.
      ee_loc = e/d_loc
    else
      call mpi%abort ('eos_t: missing 2nd argument')
    end if
    call self%lookup_table_e(dim, ee=ee_loc, d=d_loc, pg=pg, tt=tt, ne=ne, rk=rk, &
                           src=src, gamma=gamma,e_below=e_below, label=label, id=id)
  else if (present(ss) .or. present(s) .or. present(lns)) then
    if (present(ss)) then
      ee_loc => ss
    else if (present(lns)) then
      allocate (ee_loc(dim(1),dim(2),dim(3)))
      e_alloc = .true.
      ee_loc = exp(lns)/d_loc
    else if (present(s)) then
      allocate (ee_loc(dim(1),dim(2),dim(3)))
      e_alloc = .true.
      ee_loc = s/d_loc
    else
      call mpi%abort ('eos_t: missing 2nd argument')
    end if
    call self%lookup_table_s(dim, ss=ee_loc, d=d_loc, pg=pg, tt=tt, ne=ne, rk=rk, &
                           src=src, gamma=gamma,e_below=e_below, label=label, id=id)
  end if
  !-----------------------------------------------------------------------------
  if (present(e_below)) e_below=.false.
  if (d_alloc) deallocate(d_loc)
  if (e_alloc) deallocate(ee_loc)
  call trace%end(itimer)
END SUBROUTINE lookup_table

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table_e (self, dim, ee, d, pg, tt, ne, rk, &
                         src, gamma,e_below, label, id)
  class (eos_t):: self
  integer:: dim(3)
  integer, optional:: id
  real, optional:: gamma
  real, dimension(:,:,:), intent(in), pointer, optional :: d, ee
  real, dimension(:,:,:),                      optional :: pg, tt, ne
  real, dimension(:,:,:,:),                    optional :: src, rk
  logical, optional:: e_below
  character(len=*),optional::label
  real, dimension(:,:,:), pointer :: d_loc, ee_loc
  integer:: i
  integer, save:: itimer=0
  real:: stefan
  logical:: d_alloc, e_alloc
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup_e', itimer=itimer)
  nullify (d_loc, ee_loc)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  d_alloc = .false.
  e_alloc = .false.
  !-----------------------------------------------------------------------------
  d_loc => d
  ee_loc => ee
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
  if (present(e_below)) e_below=.false.
  call trace%end(itimer)
END SUBROUTINE lookup_table_e

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table_s (self, dim, ss, d, pg, tt, ne, rk, &
                           src, gamma,e_below, label, id)
  class (eos_t):: self
  integer:: dim(3)
  integer, optional:: id
  real, optional:: gamma
  real, dimension(:,:,:), intent(in), pointer, optional :: ss, d
  real, dimension(:,:,:),                      optional :: pg, tt, ne
  real, dimension(:,:,:,:),                    optional :: src, rk
  logical, optional:: e_below
  character(len=*),optional::label
  real, dimension(:,:,:), pointer :: d_loc, ss_loc
  integer:: i
  integer, save:: itimer=0
  real:: stefan, g1
  logical:: d_alloc, e_alloc
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup_s', itimer=itimer)
  nullify (d_loc, ss_loc)
  if (present(gamma)) then
    self%gamma = gamma
  end if
  !-----------------------------------------------------------------------------
  d_loc => d
  ss_loc => ss
  !-----------------------------------------------------------------------------
  if (present(tt)) then
    g1 = self%gamma-1.0
    tt = d_loc**self%gamma*exp(ss_loc*g1)/d_loc
  end if

  if (present(src)) then
    stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
    g1 = self%gamma-1.0
    do i=1,size(src,4)
      src(:,:,:,i) = stefan*real((d_loc**self%gamma*exp(ss_loc*g1)/d_loc),kind=8)**4
    end do
  end if
  if (present(pg)) then
    if (self%gamma==1d0) then
      pg = d_loc*self%csound**2
    else
      pg = d_loc**self%gamma*exp(ss_loc*(self%gamma-1d0))
    end if
  end if
  if (present(rk)) then
    do i=1,size(rk,4)
      rk(:,:,:,i) = d_loc*self%kappa
    end do
  end if
  if (present(e_below)) e_below=.false.
  call trace%end(itimer)
END SUBROUTINE lookup_table_s

END MODULE eos_mod

