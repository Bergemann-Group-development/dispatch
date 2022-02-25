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
    real:: s0 = 0.0
    integer :: n_lambda = 1
    real :: w_lambda = 1.
    logical :: do_rt = .false.
    logical :: initialized = .false.
  contains
    procedure:: init
    procedure:: inverse_lookup
    procedure:: lookup
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
  if (self%initialized) return
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
  self%initialized = .true.
  call trace%end(itimer)
END SUBROUTINE init

!===============================================================================
!> Inverse lookup, where the 2nd independent variable is log temperature, and
!> the single output variable is the energy variable
!===============================================================================
SUBROUTINE inverse_lookup (self, d, lnd, v, name, ss,ee, label)
  class (eos_t):: self
  real, dimension(:,:,:),             pointer, optional :: d, lnd, v
  real, dimension(:,:,:),             pointer :: d_loc
  real, dimension(:,:,:)            , pointer, optional :: ss, ee
  character(*), optional:: name
  integer:: n(4), m(3)
  real(8):: g1
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%inverse_lookup')
  call io%abort("legacy rt_shadow_bench EOS table module does not support inverse lookup for now.")
  call trace%end ()
END SUBROUTINE inverse_lookup

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup (self, e, ee, d, lnd, lne,lnee, pg, tt, ne, rk, &
                         ss, s, lns, src, gamma,e_below, label, id)
  class (eos_t):: self
  real, dimension(:,:,:)  , pointer, optional :: lnd,d,lne,lnee,e,ee,ne
  real, dimension(:,:,:)  , pointer, optional :: ss, s, lns
  real, dimension(:,:,:,:)         , optional :: src, rk
  real, dimension(:,:,:)           , optional :: pg, tt
  integer, optional:: id
  real(8), optional:: gamma
  logical, optional:: e_below
  character(len=*),optional::label
  ! --------------------------------------------------------------------
  integer:: dim(3)
  real, dimension(:,:,:), pointer    :: d_loc, e_loc, ss_loc

  real    :: dmin,dmax,dd, emin,emax,de,fr, m1t, m1u, t, u, stefan, g1
  integer :: i
  integer :: mx,my,mz
  integer, save:: itimer=0
  ! ....................................................................
  if(present(d)) then
    dim=shape(d)
  else 
    dim=shape(lnd)
  end if
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
  else if (present(ss)) then
    allocate (ss_loc(dim(1),dim(2),dim(3)))
    ss_loc = ss*(self%gamma-1d0)
  else
    call io%abort ('eos_t: missing internal energy')
  end if
  !
  if (associated(e_loc)) then
    if (present(tt)) then
      tt = (self%gamma-1.0)*e_loc/d_loc
    end if
    if (present(src)) then
      stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
      do i=1,size(src,4)
        src(:,:,:,i) = stefan*real(((self%gamma-1.0)*e_loc/d_loc),kind=8)**4
      end do
    end if
    if (present(pg)) then
      pg = (self%gamma-1.0)*e_loc
    end if
    if (present(rk)) then
      do i=1,size(rk,4)
        rk(:,:,:,i) = self%alfa0 *((self%gamma-1.0)*e_loc/d_loc * self%t0)**(-3.5)*(d_loc*self%rho0)**2
      end do
    end if

  else if (associated(ss_loc)) then
    if (present(tt)) then
      g1 = self%gamma-1.0
      tt = d_loc**self%gamma*exp(ss_loc)/d_loc
    end if
    if (present(src)) then
      stefan = cgs%stefan/(scaling%p*scaling%u)*scaling%temp**4
      do i=1,size(src,4)
        src(:,:,:,i) = stefan*real((d_loc**self%gamma*exp(ss_loc)/d_loc),kind=8)**4
      end do
    end if
    if (present(pg)) then
      g1 = self%gamma-1.0
      pg = d_loc**self%gamma*exp(ss_loc)
    end if
   ! kappa = (T/T0)^-3.5*(rho/rho0)^2
    if (present(rk)) then
      do i=1,size(rk,4)
        rk(:,:,:,i) = self%alfa0 *(d_loc**self%gamma*exp(ss_loc)/d_loc * self%t0)**(-3.5)*(d_loc*self%rho0)**2
      end do
    end if
  end if
  !
  if (present(lnd)) then
    deallocate (d_loc)
  end if
  if (present(lnee).or.present(lne).or.present(ee)) then
    deallocate (e_loc)
  end if
  if (present(ss)) deallocate(ss_loc)
  call trace%end (itimer)
END SUBROUTINE lookup

END MODULE eos_mod

