!===============================================================================
!> Ideal equation of state module
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
    integer :: eos_type = 0
    integer :: n_lambda = 1
    real(8) :: gamma = 1.4d0
    real    :: w_lambda = 1.0
    real    :: s0 = 0.0
    real    :: mu = 1.0
    real    :: kappa = 0.1
    real    :: csound =1.0
    logical :: do_rt = .false.
    logical :: initialized = .false.
  contains
    procedure:: init
    procedure:: lookup
    procedure:: inverse_lookup
  end type
  real(8), save :: gamma
  real(8), save :: s0
  real,    save :: mu
  real,    save :: kappa
  logical, save :: do_rt
  real,    save :: csound = 1.0
  type (eos_t), public:: eos
CONTAINS

!===============================================================================
!> Initialize parameters
!===============================================================================
SUBROUTINE init (self)
  class (eos_t):: self
  integer:: iostat
  logical, save:: first_time=.true.
  integer, save:: itimer=0
  namelist /eos_params/ gamma, mu, csound, s0, do_rt, kappa
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init', itimer=itimer)
  if (first_time) then
    call io%lock%set ('eos_t%init')
    if (first_time) then
      call io%header ('microphysics/eos/ideal/eos_mod:')
      gamma = self%gamma
      s0    = self%s0
      mu    = self%mu
      kappa = self%kappa
      do_rt = self%do_rt
      rewind (io%input)
      read (io%input, eos_params, iostat=iostat)
      if (iostat > 0) then
        rewind (io%input)
        read (io%input, eos_params)
      end if
      call io%namelist_warning ('eos_params', iostat)
      write (io%namelist, eos_params)
      first_time=.false.
    end if
    call io%lock%unset ('eos_t%init')
  end if
  self%initialized=.true.
  self%do_rt = do_rt
  self%gamma = gamma
  self%kappa = kappa * scaling%m/scaling%l**2
  self%mu    = mu
  self%s0    = s0
  call trace%end (itimer)
END SUBROUTINE init

!===============================================================================
!> Generic EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup (self, d, lnd, e, ee, lne, lnee, s, ss, pg, tt, src, rk, &
                   gamma, label, e_below, id)
  class (eos_t):: self
  real, dimension(:,:,:)   :: lnd, lne, lnee, ee, e, d, ss, s
  pointer                  :: lnd, lne, lnee, ee, e, d, ss, s
  intent(in)               :: lnd, lne, lnee, ee, e, d, ss, s
  optional                 :: lnd, lne, lnee, ee, e, d, ss, s
  real, dimension(:,:,:)   :: pg, tt
  optional                 :: pg, tt
  real, dimension(:,:,:,:) :: src, rk
  optional                 :: src, rk
  integer, optional:: id
  real(8),      optional :: gamma
  character(*), optional :: label
  logical,      optional :: e_below
  !.............................................................................
  real, dimension(:,:,:), pointer, contiguous :: d_loc, e_loc, lnd_loc
  integer :: i, dim(3)
  real(8) :: g1
  logical :: d_alloc, e_alloc
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  if (.not. self%initialized) call self%init
  if (present(gamma)) then
    self%gamma = gamma
  end if
  g1 = merge(1d0, self%gamma-1d0, self%gamma==1d0)
  d_alloc = .false.
  !-----------------------------------------------------------------------------
  ! Establish d_loc as the density
  !-----------------------------------------------------------------------------
  nullify (lnd_loc, d_loc, e_loc)
  if (present(d)) then
    dim = shape(d)
    d_loc => d
  end if
  if (present(lnd)) then
    dim = shape(lnd)
    lnd_loc => lnd
  end if
  if (.not.(associated(d_loc) .or. associated(lnd_loc))) then
    call mpi%abort ('eos_t: missing d or lnd argument')
  end if
  !-----------------------------------------------------------------------------
  ! Establish e_loc as the internal energy per unit volume
  !-----------------------------------------------------------------------------
  e_alloc = .true.
  if (present(ee)) then
    !--- energy per unit mass given ---
    allocate (e_loc(dim(1),dim(2),dim(3)))
    if (associated(lnd_loc)) then
      d_loc = exp(lnd)
    end if
    e_loc = ee*d_loc
  else if (present(lnee)) then
    !--- log energy per unit mass given ---
    allocate (e_loc(dim(1),dim(2),dim(3)))
    if (associated(lnd_loc)) then
      d_loc = exp(lnd)
    end if
    e_loc = exp(lnee)*d_loc
  else if (present(lne)) then
    !--- log energy per unit volume given ---
    allocate (e_loc(dim(1),dim(2),dim(3)))
    e_loc = exp(lne)
  else if (present(ss)) then
    !--- entropy per unit mass given ---
    allocate (e_loc(dim(1),dim(2),dim(3)))
    if (associated(lnd_loc)) then
      e_loc = (exp(eos%gamma*lnd_loc + ss*g1 - self%s0))/g1
    else
      e_loc = d_loc**self%gamma*(exp(ss*g1 - self%s0))/g1
    end if
  else if (present(s)) then
    !--- entropy per unit volume given ---
    allocate (e_loc(dim(1),dim(2),dim(3)))
    if (associated(lnd_loc)) then
      if (.not.associated(d_loc)) then
        d_loc = exp(lnd_loc)
      end if
      e_loc = (exp(eos%gamma*lnd_loc + s/d_loc*g1 - self%s0))/g1
    else
      e_loc = d_loc**self%gamma*(exp(s/d_loc*g1 - self%s0))/g1
    end if
  else if (present(e)) then
    !--- energy per unit volume given ---
    e_alloc = .false.
    e_loc => e
  else
    call mpi%abort ('eos_t: missing energy/entropy argument')
  end if
  !-----------------------------------------------------------------------------
  ! Express answers in density (d) and energy per unit volume (e)
  !-----------------------------------------------------------------------------
  if (present(pg)) then
    !--- Pg = E*(gamma-1) ---
    pg = e_loc*g1
  end if
  if (present(tt)) then
    !--- T = Pg/rho ---
    tt = e_loc*g1*self%mu/d_loc
  end if
  if (present(src)) then
    !--- T**4 source function ---
    do i=1,size(src,4)
      src(:,:,:,i) = ((e_loc*g1/d_loc)**2)**2
    end do
  end if
  if (present(rk)) then
    !--- scaled from constant kappa ---
    do i=1,size(src,4)
      rk(:,:,:,i) = d_loc*self%kappa
    end do
  end if
  !-----------------------------------------------------------------------------
  if (present(e_below)) e_below = .false.
  if (d_alloc) deallocate(d_loc)
  if (e_alloc) deallocate(e_loc)
  call trace%end(itimer)
END SUBROUTINE lookup

!===============================================================================
!> Inverse lookup, where the 2nd independent variable is log temperature, and
!> the single output variable is the energy variable
!===============================================================================
SUBROUTINE inverse_lookup (self, d, v, lnd, name, s, ss, ee, label)
  class (eos_t):: self
  real, dimension(:,:,:), pointer :: d, lnd, v
  intent(in)                      :: d, lnd, v
  optional                        :: d, lnd, v
  real, dimension(:,:,:), pointer :: s, ss, ee
  optional                        :: s, ss, ee
  character(*), optional:: name
  character(*), optional:: label
  !.............................................................................
  real, dimension(:,:,:), pointer :: d_loc, tt
  integer:: m(3)
  real(8):: g1
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%inverse_lookup', itimer=itimer)
  if (.not. self%initialized) call self%init
  !-----------------------------------------------------------------------------
  ! Temperature argument must be present
  !-----------------------------------------------------------------------------
  if (.not.present(v).or. trim(name)/='T') then
    call io%abort ('eos_t%inverse_lookup: temperature argument missing or wrong name'//trim(name))
  end if
  tt => v
  g1 = merge (self%gamma-1d0, 1d0, self%gamma>1d0)
  !-----------------------------------------------------------------------------
  ! Select density argument
  !-----------------------------------------------------------------------------
  if (present(d)) then
    d_loc => d
  else if (present(lnd)) then
    m = shape(lnd)
    allocate (d_loc(m(1),m(2),m(3)))
    d_loc = exp(lnd)
  else
    call io%abort ('eos_t%inverse_lookup: density argument missing')
  end if
  !-----------------------------------------------------------------------------
  ! Select output argument
  !-----------------------------------------------------------------------------
  if (present(ee)) then
    ee = tt*self%mu/g1
  else if (present(s)) then
    s = (self%s0+(log(tt*self%mu*d)-log(d)*eos%gamma))  &
                                  * d / (merge(self%gamma-1d0,1d0,self%gamma>1))
  else if (present(ss)) then
    ss = (self%s0+(log(tt*self%mu*d)-log(d)*eos%gamma))  &
                                 / (merge(self%gamma-1d0,1d0,self%gamma>1))
  end if
  if (present(lnd)) deallocate(d_loc)
  call trace%end (itimer)
END SUBROUTINE inverse_lookup

END MODULE eos_mod
