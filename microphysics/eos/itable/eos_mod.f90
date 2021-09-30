!===============================================================================
!> Module using a table prepared by a Python procedure such as
!> 
!> utilities/python/dispatch/EOS/stagger.py
!> utilities/python/dispatch/EOS/bifrost.py
!> utilities/python/dispatch/EOS/tomida.py
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE trace_mod
  USE itable_mod
  implicit none
  private
  type, public:: eos_t
    character(len=16) :: name='itable'
    type(itable_t)    :: table
    integer           :: nd
    integer           :: ne
    integer           :: w_lambda = 1
    integer           :: n_lambda = 1
    integer           :: itt=0, ipg=0, irk=0, ik1=0, ik2=0, is1=0, is2=0
    real              :: ee_min = 0.0, ee_max = huge(1.)
    real(8)           :: gamma = 1.4
    real, pointer, contiguous :: results(:,:,:,:) => null()
    real, pointer, contiguous :: lnd(:,:,:) => null()
  contains
    procedure         :: init
    procedure         :: lookup
    procedure         :: inverse_lookup
    procedure         :: lookup_table
    procedure         :: ddeos
  end type
  integer:: verbose=0
  character(64):: file = 'input_data/dispatch.eos'
  type(eos_t), public:: eos
  !$omp threadprivate (eos)
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(eos_t):: self
  !.............................................................................
  integer:: iostat, iv
  logical, save:: first_time = .true.
  character(32):: name
  namelist /eos_params/ verbose, file
  !-----------------------------------------------------------------------------
  call trace%begin ('itable/eos_t%init')
  if (first_time) then
    rewind (io%input)
    read (io%input, eos_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, eos_params)
    end if
    write (io%namelist, eos_params)
    first_time = .false.
  end if
  call self%table%init (trim(file))
  !-----------------------------------------------------------------------------
  ! Find the indices of the variables from their names
  !-----------------------------------------------------------------------------
  self%itt = 0
  do iv=1,self%table%nv
    name = self%table%names(iv)
    write (stdout,*) iv, name
    if (name(1:6) == 'log(T)' .or. name(1:5) == 'ln(T)' .or. name(1:1) == 'T') then
      self%itt = iv
      self%table%igen = iv
    else if (name(1:1) == 'P') then
      self%ipg = iv
    else if (name(1:2) == 'rk') then
      self%irk = iv
    else if (name(1:2) == 'K1') then
      self%ik1 = iv
    else if (name(1:2) == 'K4') then
      self%ik2 = iv
    else if (name(1:2) == 'S1') then
      self%is1 = iv
    else if (name(1:2) == 'S4') then
      self%is2 = iv
    end if 
  end do
  write (stdout,'(a,6i3)') 'itt,ipg,irk,is1,is2=', &
    self%itt, self%ipg, self%irk, self%is1, self%is2
  if (self%itt==0) then
    call io%warning ('no inverse lookup available')
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Forward table lookup, using the generic energy variable y (which typically
!> is lnE (log of internal energy) or S (entropy per unit mass)
!===============================================================================
SUBROUTINE lookup (self, d, lnd, e, ee, pg, tt, ss, rk, src, gamma, label)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in), pointer, optional :: d, lnd, e, ee
  real, dimension(:,:,:),             pointer, optional :: pg, tt, ss, rk
  real, dimension(:,:,:,:),           pointer, optional :: src
  real, dimension(:,:,:),             pointer           :: lnd_l, ee_l
  real, optional:: gamma
  character(*), optional:: label
  !.............................................................................
  integer:: nbin, m(3), n(4)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup')
  nbin = self%n_lambda
  if (present(d)) then
    m = shape(d)
  else
    m = shape(lnd)
  end if
  if (associated(self%results)) then
    n = shape(self%results)
    if (any(m /= n(1:3))) then
      deallocate (self%results)
    end if
  end if
  if (.not.associated(self%results)) then
    allocate (self%results(m(1),m(2),m(3),self%is2))
  end if
  !-----------------------------------------------------------------------------
  ! Make sure to use lnd and ee in the call to table%lookup
  !-----------------------------------------------------------------------------
  if (present(d)) then
    allocate (lnd_l(m(1),m(2),m(3)))
    lnd_l = log(d)
  else
    lnd_l => lnd
  end if
  if (present(e)) then
    allocate (ee_l(m(1),m(2),m(3)))
    if (present(d)) then
      ee_l = e/d
    else
      ee_l = e/exp(lnd_l)
    end if
  else
    ee_l => ee
  end if
 !-----------------------------------------------------------------------------
  call self%table%lookup (lnd_l, ee_l, self%results)
  !-----------------------------------------------------------------------------
  ! Clean up
  !-----------------------------------------------------------------------------
  if (present(d)) then
    deallocate (lnd_l)
  end if
  if (present(e)) then
    deallocate (ee_l)
  end if
  if (present(tt )) tt  => self%results(:,:,:,self%itt)
  if (present(pg )) pg  => self%results(:,:,:,self%ipg)
  if (present(rk )) rk  => self%results(:,:,:,self%irk)
  if (present(src)) src => self%results(:,:,:,self%is1:self%is2)
  call trace%end ()
END SUBROUTINE lookup

!===============================================================================
!> Inverse lookup, where the 2nd independent variable is log temperature, and
!> the single output variable is the energy variable
!===============================================================================
SUBROUTINE inverse_lookup (self, d, lnd, tt, y, label)
  class (eos_t):: self
  real, dimension(:,:,:),             pointer, optional :: d, lnd, tt
  real, dimension(:,:,:)            , pointer, optional :: y
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%inverse_lookup')
  call self%table%inverse_lookup (lnd, tt, self%results)
  y => self%results(:,:,:,1)
  call trace%end ()
END SUBROUTINE inverse_lookup

!===============================================================================
!> Forward table lookup, using the generic energy variable y (which typically
!> is lnE (log of internal energy) or S (entropy per unit mass)
!===============================================================================
SUBROUTINE lookup_table (self, dim, d, lnd, e, ee, pg, tt, ss, rk, src, gamma, &
                         id, e_below, label)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:), intent(in), pointer, optional :: d, lnd, e, ee
  real, dimension(:,:,:),             pointer, optional :: pg, tt, ss
  real, dimension(:,:,:,:),           pointer, optional :: src, rk
  real, optional:: gamma
  logical, optional:: e_below
  integer, optional:: id
  character(*), optional:: label
  !.............................................................................
  integer:: nbin, m(4), l(4), u(4), ng, iv
  real:: pg_min, pg_max
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup')
  nbin = self%n_lambda
  !-----------------------------------------------------------------------------
  ! Determine if we can reuse self%results
  !-----------------------------------------------------------------------------
  if (associated(self%results)) then
    m = shape(self%results)
    if (any(m(1:3) /= dim)) then
      deallocate (self%results)
    end if
  end if
  if (.not.associated(self%results)) then
    allocate (self%results(dim(1),dim(2),dim(3),self%is2))
  end if
  !---------------------.-------------------------------------------------------
  ! Detect if lnd or d is present (don't even bother to check for stupitidy)
  !-----------------------------------------------------------------------------
  call self%lookup (d=d, lnd=lnd, e=e, ee=ee)
  !-----------------------------------------------------------------------------
  ! Handle, as a non-preferred fall-back, also the case where there are
  ! pre-allocated arrays to hold the results
  !-----------------------------------------------------------------------------
  if (present(tt)) then
    if (associated(tt)) then
      tt =  self%results(:,:,:,self%itt)
    else
      tt => self%results(:,:,:,self%itt)
    end if
  end if
  if (present(pg)) then
    if (associated(pg)) then
      pg =  self%results(:,:,:,self%ipg)
    else
      pg => self%results(:,:,:,self%ipg)
    end if
    pg_min = minval(pg)
    pg_max = maxval(pg)
    if (verbose > 1) then
      write (stdlog,*) 'pg_min,pg_max=', pg_min, pg_max
    end if
  end if
  if (present(rk)) then
    if (associated(rk)) then
      rk =  self%results(l(1):u(1),l(2):u(2),l(3):u(3),self%ik1:self%ik2)
    else
      rk => self%results(:,:,:,self%ik1:self%ik2)
    end if
  end if
  if (present(src)) then
    if (associated(src)) then
      src =  self%results(l(1):u(1),l(2):u(2),l(3):u(3),self%is1:self%is2)
    else
      src => self%results(:,:,:,self%is1:self%is2)
    end if
  end if
  call trace%end ()
END SUBROUTINE lookup_table

!===============================================================================
SUBROUTINE ddeos (self, d, ee, dlnPdE_r, dlnPdlnr_E)
  class (eos_t):: self
  real, dimension(:,:), intent(in), pointer :: d, ee
  real, dimension(:,:), optional,   pointer :: dlnPdE_r
  real, dimension(:,:), optional,   pointer :: dlnPdlnr_E
  !.............................................................................
  integer:: nbin, m(3)
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%ddeos')
  if (present(dlnPdE_r  )) dlnPdE_r   = 1.0/ee
  if (present(dlnPdlnr_E)) dlnPdlnr_E = 1.0
  call trace%end ()
END SUBROUTINE ddeos

END MODULE eos_mod
