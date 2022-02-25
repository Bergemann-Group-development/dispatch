!===============================================================================
!> Module using a table prepared by a Python procedure.  This module is
!> initialized while still in single-thread mode, so no OMP locks are not
!> needed. The input data table is in a custom format, which may be written by
!> a custom reader, written in Fortran or Python.
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE scaling_mod
  USE trace_mod
  USE table_mod
  private
  type, public:: eos_t
    character(len=16) :: name='Tomida+Hori(e)'
    type(table_t)     :: table
    integer :: ne, nd, nv
    integer :: w_lambda = 1
    integer :: n_lambda = 1
    real, pointer:: lnd(:) => null()
    real(8) :: gamma = 1.4
  contains
    procedure :: init
    procedure :: lookup
    procedure :: inverse_lookup
  end type
  logical, save:: first_time=.true.
  type(eos_t), public:: eos
CONTAINS

!===============================================================================
! Initialize the eos table data, table axes and dimensions
!===============================================================================
SUBROUTINE init (self)
  class(eos_t):: self
  integer:: i, iostat
  real(8):: gamma
  real, pointer:: lnd(:,:,:), e(:,:,:), tg(:,:,:), pg(:,:,:)
  integer, save:: itimer=0
  namelist /eos_params/ gamma
  !-----------------------------------------------------------------------------
  call io%lock%set ('eos_t%init')
  if (first_time) then
    first_time=.false.
    rewind (io%input)
    read (io%input, eos_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, eos_params)
    end if
    call io%namelist_warning ('eos_params', iostat)
    write (io%namelist, eos_params)
    call self%table%init ('data/eos/tomida/e/eos.dat')
    self%ne = self%table%n1
    self%nd = self%table%n2
    self%nv = self%table%n3
    self%table%x_min = self%table%x_min - log(scaling%d)
    self%table%y_min = self%table%y_min - log(scaling%e)
    allocate (self%lnd(self%nd))
    do i=1,self%nd
      self%lnd(i) = self%table%x_min + (i-1)*self%table%x_step
    end do
    ! test
    ! d=2.025e-07 T=5.800e+03 E=2.058e+12
    allocate(lnd(1,1,1))
    allocate(e(1,1,1))
    allocate(tg(1,1,1))
    allocate(pg(1,1,1))
    lnd=log(2.0253139e-07)-log(scaling%d)
    e= 2057648493081.844/scaling%e!/exp(lnd)
    print *, 'table test'
    print *, 'input d, e', exp(lnd*scaling%d), e*exp(lnd)*scaling%e
    print *, 'code input d, e', exp(lnd), e*exp(lnd)
    print *, 'input lnd, ee, logee', lnd, e, log(e)
    call self%lookup(lnd=lnd,ee=e,tt=tg,pg=pg)
    print *, 'result t, pg(code), pg', tg*scaling%temp,pg, pg*scaling%p
    !minmax d   160522.5       381970.3    
    !minmax e  1.6495171E+09  6.5163566E+09
    !lnd=log(160522.5)
    !lnd=log(424414.3)
    !e=6.5163566E+16
    !e=7.7174094E+16
    !print *, 'dd', exp(lnd)*scaling%d
    !call self%lookup(lnd=lnd,e=e,tt=tg,pg=pg)
    !print *, 'result t, pg(code), pg', tg*scaling%temp,pg, pg*scaling%p
    deallocate(lnd,e,tg, pg)
  end if
  call io%lock%unset ('eos_t%init')
END SUBROUTINE init


!===============================================================================
!> Forward table lookup, using the generic energy variable y (which typically
!> is lnE (log of internal energy) or S (entropy per unit mass)
!===============================================================================
SUBROUTINE lookup (self, d, lnd, e, ee, s, pg, tt, ss, rk, src, gamma, label, id)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in), pointer, optional :: d, lnd, e, ee, s
  real, dimension(:,:,:),                      optional :: pg, tt, ss
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, optional:: gamma(:,:,:)
  integer, optional::id
  character(*), optional:: label
  ! ----------------------------------------------------------------------------
  real, dimension(:,:,:),             pointer           :: lnd_l, lnee_l
  integer:: m(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup_table')
  if (present(d)) then
    m = shape(d)
  else
    m = shape(lnd)
  end if
  if (present(d)) then
    allocate (lnd_l(m(1),m(2),m(3)))
    lnd_l = log(d)
  else
    lnd_l => lnd
  end if
  if (present(e)) then
    allocate (lnee_l(m(1),m(2),m(3)))
    lnee_l = log(e)-lnd_l
  else if (present(ee)) then
    allocate (lnee_l(m(1),m(2),m(3)))
    lnee_l=log(ee)
  end if
  call self%table%lookup (m, lnd_l, lnee_l, pg, tt, ss, rk, src, gamma)
  if (present(d)) then
    deallocate (lnd_l)
  end if
  deallocate (lnee_l)
  call trace%end ()
END SUBROUTINE lookup

!===============================================================================
!> Inverse lookup, where the 2nd independent variable is log temperature, and
!> the single output variable is the energy variable
!===============================================================================
SUBROUTINE inverse_lookup (self, d, lnd, tt, e, ee, ss, label)
  class (eos_t):: self
  real, dimension(:,:,:),             pointer, optional :: d, lnd, tt
  real, dimension(:,:,:),             pointer, contiguous:: lnd_loc
  real, dimension(:,:,:)            , pointer, optional :: e, ee, ss
  integer:: n(4), m(3)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%inverse_lookup')
  call io%abort ('tomida/e eos inverse_lookup called - need to set up!')
  call trace%end ()
END SUBROUTINE inverse_lookup



END MODULE eos_mod
