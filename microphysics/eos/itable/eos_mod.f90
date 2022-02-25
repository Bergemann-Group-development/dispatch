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
  USE scaling_mod
  USE units_mod
  implicit none
  private
  real, pointer, contiguous:: results(:,:,:,:)
  type, public:: eos_t
    character(len=16) :: name='itable'
    type(itable_t)    :: table
    integer           :: nd
    integer           :: ne
    integer           :: w_lambda = 1
    integer           :: n_lambda = 1
    integer           :: itt=0, ipg=0, irk=0, ik1=0, ik2=0, is1=0, is2=0
    real              :: ee_min = 0.0, ee_max = huge(1.)
    real(8)           :: gamma = 1.2
    real(8)           :: s0 = 0.0
    !real, pointer, contiguous :: results(:,:,:,:) => null()
    real, pointer, contiguous :: lnd(:,:,:) => null()
    logical           :: initialized=.false.
  contains
    procedure         :: init
    procedure         :: init_warning
    procedure         :: lookup
    procedure         :: inverse_lookup
    procedure         :: lookup_table
    procedure         :: ddeos
    procedure         :: rescale
  end type
  integer:: verbose=0
  character(64):: file = 'input_data/dispatch.eos'
  character(64):: dir = ''
  character(64):: top = '../../'
  logical, save:: convert_to_code=.false.         ! Convert to code units?
  type(eos_t), public:: eos
  !$omp threadprivate (results)
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (self)
  class(eos_t):: self
  !.............................................................................
  integer:: iostat, iv
  real(8),save:: gamma=1.4
  logical, save:: first_time = .true.
  logical, save:: first_time_print = .true.
  character(32):: name
  real, pointer:: lnd(:,:,:), e(:,:,:),e2(:,:,:), tg(:,:,:), pg(:,:,:)
  namelist /eos_params/ verbose, top, dir, file, gamma, convert_to_code
  !-----------------------------------------------------------------------------
  call trace%begin ('itable/eos_t%init')
  if (self%initialized) return
  if (first_time) then
    rewind (io%input)
    read (io%input, eos_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, eos_params)
    end if
    self%gamma=gamma
    write (io%namelist, eos_params)
    first_time = .false.
  end if
  call self%table%init (trim(top)//trim(dir)//trim(file))
  !-----------------------------------------------------------------------------
  ! Find the indices of the variables from their names
  !-----------------------------------------------------------------------------
  self%itt = 0
  if (first_time_print) then
      write (stdout,*) 'scale1 ', trim(self%table%scale1%name), &
                       ',scale2 ', trim(self%table%scale2%name)
  end if
  do iv=1,self%table%nv
    name = self%table%names(iv)
    if (first_time_print) then
      write (stdout,*) iv, name
    end if
    if (name(1:6) == 'log(T)' .or. name(1:5) == 'ln(T)' .or. name(1:1) == 'T') then
      self%itt = iv
      self%table%igen = iv
    else if (name(1:1) == 'P') then
      self%ipg = iv
    else if (name(1:4) == 'rkap') then
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
  if (convert_to_code) call self%rescale(first_time_print)
  !!!
  self%initialized=.true.
  if (first_time_print) then
    write (stdout,'(a,6i3)') 'itt,ipg,irk,is1,is2=', &
    self%itt, self%ipg, self%irk, self%is1, self%is2
    first_time_print = .false.
    allocate(lnd(1,1,1))
    allocate(e(1,1,1))
    allocate(e2(1,1,1))
    allocate(tg(1,1,1))
    allocate(pg(1,1,1))
    lnd=log(2.0251139e-7)-log(scaling%d)
    e= 2057648493081.844/scaling%e*scaling%d!exp(lnd)
    print *, 'minmax ee', self%table%scale2%min, self%table%scale2%max
    print *, 'ee', e, e*exp(lnd)
    print *, 'table test'
    print *, 'input d, e', lnd+log(scaling%d), e*scaling%e/scaling%d*exp(lnd)*scaling%d
    print *, 'code input d, e', exp(lnd), e*exp(lnd)
    print *, 'input lnd, ee, logee', lnd, e, log(e)
    call self%lookup(lnd=lnd,ee=e,tt=tg,pg=pg)
    print *, 'result t, pg(code), pg', tg*scaling%temp,pg, pg*scaling%p
    call self%inverse_lookup(lnd=lnd,v=tg,name='T',ee=e2)
    print *, 'inverse lookup ee', e2, 'diff', e2-e, (e2-e)/e
    lnd=exp(lnd)
    call self%inverse_lookup(d=lnd,v=tg,name='T',ee=e2)
    print *, 'inverse lookup ee', e, e2-e, (e2-e)/e
    deallocate(lnd,e,e2,tg, pg)
    if (self%table%get_index('T')==self%itt) then
      print *, 'itt index OK', self%itt, self%table%get_index('T')
    end if
    if (self%table%get_index('P')==self%ipg) then
      print *, 'ipg index OK', self%ipg, self%table%get_index('P')
    end if
  end if
    
  if (self%itt==0) then
    call io%warning ('no inverse lookup available')
  end if
  call trace%end ()
END SUBROUTINE init

!===============================================================================
!> Optionally rescale the table to code units
!===============================================================================
SUBROUTINE rescale(self, first_time_print)
  class(eos_t):: self
  logical:: first_time_print
  ! ----------------------------------------------------------------------------
  class(independent_t), pointer:: scale
  class(dependent_t), pointer:: val
  character(32):: name
  real(8):: sc
  integer:: i, iv
  ! ............................................................................
  ! First, rescale the scales
  name = self%table%scale1%name
  if (verbose>0.and.first_time_print) then
    write (stdout,*) 'scale 1 '//name//' min, max, step', &
    self%table%scale1%min, self%table%scale1%max, self%table%scale1%step
  end if
  if (name(1:3)=='rho' .or. name(1:1)=='d') then
    sc = 1d0/scaling%d
    self%table%scale1%min  = scale_mult(self%table%scale1%min ,sc)
    self%table%scale1%max  = scale_mult(self%table%scale1%max ,sc)
    self%table%scale1%step = scale_mult(self%table%scale1%step,sc)
  else if (name(1:7)=='ln(rho)' .or. name(1:5)=='ln(d)') then
    sc = -log(scaling%d)
    self%table%scale1%min  = scale_add(self%table%scale1%min , sc)
    self%table%scale1%max  = scale_add(self%table%scale1%max , sc)
  end if
  if (verbose>0.and.first_time_print) then
    write (stdout,*) 'scale 1 (CODE) '//name//' min, max, step', &
    self%table%scale1%min, self%table%scale1%max, self%table%scale1%step
  end if
  !
  name = self%table%scale2%name
  if (verbose>0.and.first_time_print) then
    write (stdout,*) 'scale 2 '//name//' min, max, step', &
    self%table%scale2%min, self%table%scale2%max, self%table%scale2%step
  end if
  if (name(1:2)=='ee') then
    sc = scaling%d / scaling%e
    self%table%scale2%min  = scale_mult(self%table%scale2%min , sc) 
    self%table%scale2%max  = scale_mult(self%table%scale2%max , sc) 
    self%table%scale2%step = scale_mult(self%table%scale2%step, sc) 
    do i=1,self%table%scale1%n
      scale => self%table%scales(i)
      scale%min  = scale_mult(scale%min ,sc)
      scale%max  = scale_mult(scale%max ,sc)
      scale%step = scale_mult(scale%step,sc)
    end do
  else if (name(1:6)=='ln(ee)') then
    sc = - log(scaling%e) + log(scaling%d)
    self%table%scale2%min  = scale_add(self%table%scale2%min ,sc)
    self%table%scale2%max  = scale_add(self%table%scale2%max ,sc)
    do i=1,self%table%scale1%n
      scale => self%table%scales(i)
      scale%min  = scale_add(scale%min ,sc)
      scale%max  = scale_add(scale%max ,sc)
      scale%step = scale_add(scale%step,sc)
    end do
  else if (name(1:5)=='ln(E)') then
    sc = - log(scaling%e)
    self%table%scale2%min  = scale_add(self%table%scale2%min ,sc)
    self%table%scale2%max  = scale_add(self%table%scale2%max ,sc)
    do i=1,self%table%scale1%n
      scale => self%table%scales(i)
      scale%min  = scale_add(scale%min ,sc)
      scale%max  = scale_add(scale%max ,sc)
      scale%step = scale_add(scale%step,sc)
    end do
  else if (name(1:1)=='E') then
    sc = 1d0 / scaling%e
    self%table%scale2%min  = scale_mult(self%table%scale2%min , sc) 
    self%table%scale2%max  = scale_mult(self%table%scale2%max , sc) 
    self%table%scale2%step = scale_mult(self%table%scale2%step, sc) 
    do i=1,self%table%scale1%n
      scale => self%table%scales(i)
      scale%min  = scale_mult(scale%min ,sc)
      scale%max  = scale_mult(scale%max ,sc)
      scale%step = scale_mult(scale%step,sc)
    end do
  else
    call io%abort('eos%rescale:: '//trim(name)//' rescaling not set up')
  end if
  if (verbose>0.and.first_time_print) then
    write (stdout,*) 'scale 2 (CODE) '//name//' min, max, step', &
    self%table%scale2%min, self%table%scale2%max, self%table%scale2%step
  end if
  ! ----------------------------------------------------------------------------
  ! now rescale all dependent values
  ! ----------------------------------------------------------------------------
  do i=1,self%table%scale1%n
    val => self%table%values(i)
    do iv=1, self%table%nv
    name = self%table%names(iv)
      if (name(1:6) == 'log(T)') then
        val%v(iv,:)=val%v(iv,:) - log10(scaling%temp)
      else if (name(1:5) == 'ln(T)') then
        val%v(iv,:) = val%v(iv,:) - log(scaling%temp) 
      else if(name(1:1) == 'T') then
        val%v(iv,:) = val%v(iv,:) / scaling%temp
      else if (name(1:1) == 'P') then
        val%v(iv,:) = val%v(iv,:) / scaling%p
      else if (name(1:4) == 'rkap' .or. &
               name(1:2) == 'kr'   .or. &
               name(1:2) == 'kp') then
        val%v(iv,:) = val%v(iv,:) /(scaling%l**2/scaling%m)
      else if (name(1:2) == 'cs'.or. name(1:3)=='C_T')  then
        val%v(iv,:) = val%v(iv,:) / scaling%u
      else if (name(1:3) == 'Src')  then
        val%v(iv,:) = val%v(iv,:) / scaling%e
      end if 
    end do
  end do

  ! ----------------------------------------------------------------------------
  contains
    function scale_mult (i, s) result (v)
      real(8):: v,i,s
      v=i*s
    end function scale_mult
    function scale_add (i, s) result (v)
      real(8):: v,i,s
      v=i+s
    end function scale_add
END SUBROUTINE rescale

SUBROUTINE init_warning(self)
  class(eos_t):: self
  ! ----------------------------------------------------------------------------
  if (self%initialized) return
  print *, 'WARNING: EOS was not initialized!! -- remove this if better solution found.'
  call self%init
END SUBROUTINE init_warning

!===============================================================================
!> Forward table lookup, using the generic energy variable y (which typically
!> is lnE (log of internal energy) or S (entropy per unit mass)
!===============================================================================
SUBROUTINE lookup (self, d, lnd, e, ee, pg, tt, ss, s, rk, src, gamma, var, label, name, id)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in), pointer, optional :: d, lnd, e, ee
  real, dimension(:,:,:),             pointer, optional :: pg, tt, ss, s, rk, gamma, var
  real, dimension(:,:,:,:),           pointer, optional :: src
  real, dimension(:,:,:),             pointer           :: lnd_l, ee_l
  integer, optional:: id
  character(*), optional:: label
  character(*), optional:: name
  !.............................................................................
  integer:: nbin, m(3), n(4), idx, ivar
  !-----------------------------------------------------------------------------
  call self%init_warning
  call trace%begin ('eos_t%lookup')
  nbin = self%n_lambda
  if (present(d)) then
    m = shape(d)
  else
    m = shape(lnd)
  end if
  if (associated(results)) then
    n = shape(results)
    if (any(m /= n(1:3))) then
      deallocate (results)
    end if
  end if
  if (.not.associated(results)) then
    allocate (results(m(1),m(2),m(3),self%table%nv))
  end if
  ! ----------------------------------------------------------------------------
  ! Look up custom item in the table. Note: the EOS talbes do not guarantee,
  ! that this item will be present. If it is not in the EOS tables, an abort
  ! signal will be triggered.
  ! ----------------------------------------------------------------------------
  if (present(name).and.present(var)) then
    ivar=self%table%get_index(name)
  else
    ivar=-1
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
  if (self%table%scale2%name(1:2)=='ee') then
    if (present(e)) then
      allocate (ee_l(m(1),m(2),m(3)))
      if (present(d)) then
        ee_l = e/d
      else
        ee_l = e/exp(lnd_l)
      end if
    else if (present(ee)) then
      ee_l => ee
    end if
    if (present(ss)) then
      allocate (ee_l(m(1),m(2),m(3)))
      ee_l = ss*(self%gamma-1.0)!*exp(lnd_l)
    end if
  else if (self%table%scale2%name(1:6)=='ln(ee)') then
    allocate (ee_l(m(1),m(2),m(3)))
    if (present(e)) then
      if (present(d)) then
        ee_l = log(e/d)
      else
        ee_l = log(e)-lnd_l
      end if
    else if (present(ee)) then
      ee_l = log(ee)
    end if
    if (present(ss)) then
      allocate (ee_l(m(1),m(2),m(3)))
      ee_l = log(ss*(self%gamma-1.0))!*exp(lnd_l)
    end if
  end if
    

 !-----------------------------------------------------------------------------
  !call self%table%lookup (lnd_l, ee_l, self%results)
  call self%table%lookup (lnd_l, ee_l, results)
  !-----------------------------------------------------------------------------
  ! Clean up
  !-----------------------------------------------------------------------------
  if (present(d)) then
    deallocate (lnd_l)
  end if
  if (present(e).or.present(ss).or.self%table%scale2%name(1:6)=='ln(ee)') then
    deallocate (ee_l)
  end if
  if (present(tt )) then
    if (associated(tt )) then
      tt  =  results(:,:,:,self%itt)
    else
      tt  => results(:,:,:,self%itt)
    end if
  end if
  if (present(pg )) then
    if (associated(pg )) then
      pg  = results(:,:,:,self%ipg)
    else
      pg  => results(:,:,:,self%ipg)
    end if
  end if
  if (present(rk )) then
    if (associated(rk )) then
      rk  = results(:,:,:,self%irk)
    else
      rk  => results(:,:,:,self%irk)
    end if
  end if
  if (present(gamma )) then
    idx = self%table%get_index('gamma')
    if (associated(gamma )) then
      gamma  = results(:,:,:,idx)
    else
      gamma  => results(:,:,:,idx)
    end if
  end if
  if (present(src)) then
    if (associated(src)) then
      src = results(:,:,:,self%is1:self%is2)
    else
      src => results(:,:,:,self%is1:self%is2)
    end if
  end if
  if (ivar>0) then
    var = results(:,:,:,ivar)
  end if
  call trace%end ()
END SUBROUTINE lookup

!===============================================================================
!> Inverse lookup, where the 2nd independent variable is any variable in EOS, and
!> the single output variable is the energy variable
!===============================================================================
SUBROUTINE inverse_lookup (self, d, lnd, v, name, ee, ss, label)
  class (eos_t):: self
  real, dimension(:,:,:),             pointer, optional :: d, lnd, v
  character(*), optional:: name
  real, dimension(:,:,:),             pointer, contiguous:: lnd_loc
  real, dimension(:,:,:)            , pointer, optional :: ee, ss
  integer:: n(4), m(3)
  character(*), optional:: label
  !-----------------------------------------------------------------------------
  call self%init_warning
  call trace%begin ('eos_t%inverse_lookup')
  if (present(d)) then
    m=shape(d)
    allocate(lnd_loc(m(1),m(2),m(3)))
    lnd_loc=log(d)
  else if (present(lnd)) then
    lnd_loc=> lnd
    m=shape(lnd)
  end if
  if (associated(results)) then
    n = shape(results)
    if (any(m /= n(1:3).or.n(4)/=1)) then
      deallocate (results)
    end if
  end if
  if (.not.associated(results)) then
    allocate (results(m(1),m(2),m(3),1))
  end if
  if (present(name)) then
    call self%table%inverse_lookup (lnd_loc, v, self%table%get_index(name),results)
  else
    call self%table%inverse_lookup (lnd_loc, v, self%itt,results)
  end if
  if (present(ee)) then
    if (self%table%scale2%name(1:6)=='ln(ee)') then
      ee = exp(results(:,:,:,1))
    else
      ee = results(:,:,:,1)
    end if
  end if
  if (present(ss)) ss = results(:,:,:,1)
  if (present(d)) deallocate(lnd_loc)
  call trace%end ()
END SUBROUTINE inverse_lookup

!===============================================================================
!> Forward table lookup, using the generic energy variable y (which typically
!> is lnE (log of internal energy) or S (entropy per unit mass)
!===============================================================================
SUBROUTINE lookup_table (self, dim, d, lnd, e, ee, pg, tt, ss, rk, src, gamma, &
                         id, e_below, label)
  class (eos_t):: self
  integer:: dim(3), irk
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
  call self%init_warning
  call trace%begin ('eos_t%lookup')
  nbin = self%n_lambda
  !-----------------------------------------------------------------------------
  ! Determine if we can reuse self%results
  !-----------------------------------------------------------------------------
  if (associated(results)) then
    m = shape(results)
    if (any(m(1:3) /= dim)) then
      deallocate (results)
    end if
  end if
  if (.not.associated(results)) then
    allocate (results(dim(1),dim(2),dim(3),self%table%nv))
  end if
  !---------------------.-------------------------------------------------------
  ! Detect if lnd or d is present (don't even bother to check for stupitidy)
  !-----------------------------------------------------------------------------
  call self%lookup (d=d, lnd=lnd, e=e, ee=ee,ss=ss)
  !-----------------------------------------------------------------------------
  ! Handle, as a non-preferred fall-back, also the case where there are
  ! pre-allocated arrays to hold the results
  !-----------------------------------------------------------------------------
  if (present(tt)) then
    if (associated(tt)) then
      tt =  results(:,:,:,self%itt)
    else
      tt => results(:,:,:,self%itt)
    end if
  end if
  if (present(pg)) then
    if (associated(pg)) then
      pg =  results(:,:,:,self%ipg)
    else
      pg => results(:,:,:,self%ipg)
    end if
    pg_min = minval(pg)
    pg_max = maxval(pg)
    if (verbose > 1) then
      write (stdlog,*) 'pg_min,pg_max=', pg_min, pg_max
    end if
  end if
  if (present(rk)) then
    if (associated(rk)) then
      !rk(:,:,:,1)=self%results(:,:,:,self%irk)
      !do irk=2,4
      !  rk(:,:,:,irk) =  rk(:,:,:,irk-1)*10.0
      !end do
      rk = results(:,:,:,self%ik1:self%ik2)
    else
      rk => results(:,:,:,self%ik1:self%ik2)
    end if
  end if
  if (present(src)) then
    if (associated(src)) then
      src =  results(:,:,:,self%is1:self%is2)
    else
      src => results(:,:,:,self%is1:self%is2)
    end if
    src = exp(src)
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
  call self%init_warning
  call trace%begin ('eos_t%ddeos')
  if (present(dlnPdE_r  )) dlnPdE_r   = 1.0/ee
  if (present(dlnPdlnr_E)) dlnPdlnr_E = 1.0
  call trace%end ()
END SUBROUTINE ddeos

END MODULE eos_mod
