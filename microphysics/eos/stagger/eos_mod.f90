!===============================================================================
!> EOS module, reading tabular EOS compatible with STAGGER code
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE mpi_mod
  USE units_mod
  USE scaling_mod
  USE kinds_mod
  implicit none
  private
  !-----------------------------------------------------------------------------
  ! Table related parameters
  !-----------------------------------------------------------------------------
  integer iupdte,nvar,md,mtable,njon
  real(kind=4):: eosxmin,ur,ul,ut,eps,tff,grv,abnd,dbox
  real(kind=4), allocatable, dimension(:):: tmean,tamp,rhm,xcorr,thmin,thmax, &
    dth,eemin,eemax,deetab,tab
  integer, allocatable, dimension(:):: itab,mtab
  integer, save:: ncall=0, mbox
  character(len=64):: top, tablefile, file
  logical, save:: do_eos
  logical, save:: first_time=.true.
  !-----------------------------------------------------------------------------
  ! Public part
  !-----------------------------------------------------------------------------
  integer, parameter:: n_lambda=4
  type, public:: eos_t
    integer :: n_lambda = n_lambda
    real    :: w_lambda(n_lambda) = 1.0
    real    :: gamma
    real    :: kappa
    logical :: on = .false.
  contains
    procedure:: init
    procedure:: lookup
    procedure:: lookup_table
    procedure:: pressure
    procedure:: temperature
  end type
  integer, save:: verbose=0
  type (eos_t), public:: eos
CONTAINS

!===============================================================================
!> Initialize EOS module (only once)
!===============================================================================
SUBROUTINE init (self)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: d, ee, pg
  real, save :: kappa = 1.0
  integer:: ir
  namelist /eos_params/ do_eos, top, tablefile, kappa, verbose
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  !$omp critical (eos_cr)
  if (first_time) then
    first_time = .false.
    do_eos = .true.
    top       = '../../'
    tablefile = 'data/eos/stagger/table.dat'
    rewind (io%input); read (io%input,eos_params)
    if (io%master) write (*,eos_params)
    !-----------------------------------------------------------------------------
    ! Read the table file
    !-----------------------------------------------------------------------------
    file = trim(top)//'/'//trim(tablefile)
    open (io_unit%data, file=trim(file), form='unformatted', status='old', action='read')
    read (io_unit%data) md,iupdte,nvar,mbox,eosxmin,dbox,ul,ut,ur,eps,tff,grv,abnd
    njon = nvar-mbox
    allocate (tmean(md),tamp(md),rhm(md),xcorr(md),thmin(md),thmax(md))
    allocate (dth(md),eemin(md),eemax(md),deetab(md),itab(md),mtab(md))
    if (io%do_trace.and.io%master) print *,mpi%rank, 'reading 12 times', md, ' values'
    read (io_unit%data) tmean,tamp,rhm,xcorr,thmin,thmax,dth,eemin,eemax,deetab,itab,mtab
    read (io_unit%data) mtable
    allocate (tab(mtable))
    read (io_unit%data) tab
    close (io_unit%data)
    if (mpi%master) then
      print *, trim(tablefile)//': md, mtable =', md, mtable
      allocate (d(1,1,1), ee(1,1,1), pg(1,1,1))
      d = 1.0
      ee  = 5.0
      !call self%lookup (shape(ee), ee=ee, d=d, pg=pg)
      print *,'d, ee, pg =', d, ee, pg
      deallocate (d, ee, pg)
      if (verbose > 0) then
        print *, '      density       min(E)      max(E)        min(T)      max(T)'
        do ir=1,md
          print '(i4,1p,e12.4,2(2x,2e12.4))', ir, rhm(ir), eemin(ir), eemax(ir), &
            tab(itab(ir)+2*3*mtab(ir)), tab(itab(ir)+(2*3+1)*mtab(ir)-1)
        end do
      end if
    end if
  end if
  self%on = do_eos
  self%kappa = kappa
  !$omp end critical (eos_cr)
  !$omp barrier
  call trace%end()
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
  real, dimension(:,:,:), pointer :: x_loc, y_loc, lnd_l, ee_l
  integer:: i
  integer, save:: itimer=0
  real:: stefan
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
  if (present(lnd).and.present(ee)) then
    if (present(src)) then
      call lookup_table (self, dim, ee=ee, lnd=lnd, rk=rk, src=src)
    else
      call lookup_table (self, dim, ee=ee, lnd=lnd, rk=rk)
    end if
  else if (present(src).and.present(rk)) then
    allocate (lnd_l(dim(1),dim(2),dim(3)))
    allocate (ee_l (dim(1),dim(2),dim(3)))
    lnd_l = log(x_loc)
    ee_l = y_loc/x_loc
    call lookup_table (self, dim, ee=ee_l, lnd=lnd_l, src=src, rk=rk)
    deallocate (lnd_l, ee_l)
  else
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
  end if
  if (present(lnx)) deallocate(x_loc)
  if (present(lny)) deallocate(y_loc)
  call trace%end(itimer)
END SUBROUTINE lookup

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table (self, dim, e, ee, d, lnd, lne, pg, tt, ne, rk, src, &
                         ss, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:)  , pointer, optional :: lnd,d,lne,e,ee,tt,ne,ss
  real, dimension(:,:,:,:)         , optional :: src,rk
  real, dimension(:,:,:)           , optional :: pg
  real, optional:: gamma
  !.............................................................................
  real, dimension(:,:,:), pointer    :: lnd_loc, ee_loc
  integer,         dimension(dim(1)) :: np,np1,ik,ntab,ik1,ntab1
  real(kind=4),    dimension(dim(1)) :: px,py,f00,f01,f10,f11, &
                                        fx00,fx01,fx10,fx11,fy00,fy01,fy10,fy11
  real    :: rhm1,rhm2,drhm,algrk,eek,eek1,py1
  real    :: qx,qy,pxqx,pxpx,qxqx,pypy,pyqy,qyqy,pxpy,qxqy,pxqy,qxpy
  logical :: newtable
  integer :: ix,iy,iz,kee,kee1,j
  integer :: mx,my,mz,mbox,nbelow_d, nbelow_e
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup_table', itimer=itimer)
  if (io%verbose>2) then
    print *, 'present(lnd) =', present(lnd)
    print *, 'present(ee)  =', present(ee)
    print *, 'present(pg)  =', present(pg)
    print *, 'present(src) =', present(src)
    print *, 'present(rk)  =', present(rk)
    print *, 'dim:', dim
  end if
  mx=dim(1)
  my=dim(2)
  mz=dim(3)
  newtable = .false.
  if (present(lnd)) then
    lnd_loc => lnd
  else if (present(d)) then
    allocate (lnd_loc(dim(1),dim(2),dim(3)))
    lnd_loc = log(d)
  else
    call mpi%abort ('eos_t: missing density')
  end if
  if (present(ee)) then
    ee_loc => ee
  else if (present(e)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    if (present(d)) then
      ee_loc = e/d
    else
      ee_loc = e/exp(lnd_loc)
    end if
  else if (present(lne)) then
    allocate (ee_loc(dim(1),dim(2),dim(3)))
    ee_loc = exp(lne)
  else
    call mpi%abort ('eos_t: missing internal energy')
  end if
  !
  rhm1=log(rhm(1))
  rhm2=log(rhm(md))
  drhm=(rhm2-rhm1)/(md-1)
  !-----------------------------------------------------------------------------
  ! Loop over points
  !-----------------------------------------------------------------------------
  nbelow_d = 0
  nbelow_e = 0
  do iz=1,mz
  do iy=1,my
    !---------------------------------------------------------------------------
    ! Density index (6 lops)
    !---------------------------------------------------------------------------
    do ix=1,mx
      algrk=1.+(lnd_loc(ix,iy,iz)-rhm1)/drhm
      nbelow_d = nbelow_d + merge(1,0,lnd_loc(ix,iy,iz)<rhm1)
      np(ix)=max0(1,min0(md-1,int(algrk)))
      px(ix)=algrk-np(ix)
    end do
    do ix=1,mx
      !-------------------------------------------------------------------------
      ! Internal energy index -- this loop has indirect referencing (9 flops)
      !-------------------------------------------------------------------------
      ntab(ix)  = mtab(np(ix))
      eek       = 1.+(ee_loc(ix,iy,iz)-eemin(np(ix)))/deetab(np(ix))
      nbelow_e  = nbelow_e + merge(1,0,eek<1.0)
      eek       = max(eek,1.0)
      kee       = min(ntab(ix)-1,max(1,int(eek)))
      py(ix)    = eek-kee
      ik(ix)    = itab(np(ix))+kee-1
      ntab1(ix) = mtab(np(ix)+1)
      ! -- the integer offset below adjusts for a possible larger eemin at the
      ! -- next table density.  Zero for square tables (e.g. default table.dat)
      kee1      = kee - nint((eemin(np(ix)+1)-eemin(np(ix))) &
                                            /deetab(np(ix)))
      kee1      = min0(ntab1(ix)-1,max0(1,kee1))
      ik1(ix)   = itab(np(ix)+1)+kee1-1
    end do
    !---------------------------------------------------------------------------
    ! Separate loop for weights, vectorizes well
    !---------------------------------------------------------------------------
    do ix=1,mx
! +5 = 20 flops
      qx   = 1. - px(ix)
      pxqx = px(ix) * qx
      pxpx = px(ix) * px(ix)
      qxqx = qx * qx
! +4 = 24 flops
      qy   = 1. - py(ix)
      pyqy = py(ix) * qy
      pypy = py(ix) * py(ix)
      qyqy = qy     * qy
! +4 = 28 flops
      pxqy = px(ix) * qy
      pxpy = px(ix) * py(ix)
      qxqy = qx     * qy
      qxpy = qx     * py(ix)
! +5 = 33 flops
      f00(ix) = qxqy * (1. + pxqx - pxpx + pyqy - pypy)
      f01(ix) = qxpy * (1. + pxqx - pxpx + pyqy - qyqy)
      f10(ix) = pxqy * (1. - qxqx + pxqx + pyqy - pypy)
      f11(ix) = pxpy * (1. - qxqx + pxqx + pyqy - qyqy)
! +4 = 37 flops
      fx00(ix) =    qxqy * pxqx
      fx01(ix) =    qxpy * pxqx
      fx10(ix) =  - pxqy * pxqx
      fx11(ix) =  - pxpy * pxqx
! +4 = 41 flops
      fy00(ix) =    qxqy * pyqy
      fy01(ix) =  - qxpy * pyqy
      fy10(ix) =    pxqy * pyqy
      fy11(ix) =  - pxpy * pyqy
    end do
    !---------------------------------------------------------------------------
    !  Loop over table entries, nvar = 4.  Would this be more efficent
    !  if split into one gather loop and one calculation loop?
    !---------------------------------------------------------------------------
    if (present(pg)) then
     do ix=1,mx
      pg(ix,iy,iz) = exp( &
         f00(ix) * tab(ik (ix))+  f01(ix) * tab(ik (ix) + 1) &
      +  f10(ix) * tab(ik1(ix))+  f11(ix) * tab(ik1(ix) + 1) &
      + fx00(ix) * tab(ik (ix) + 2 * ntab (ix)    ) &
      + fx01(ix) * tab(ik (ix) + 2 * ntab (ix) + 1) &
      + fx10(ix) * tab(ik1(ix) + 2 * ntab1(ix)    ) &
      + fx11(ix) * tab(ik1(ix) + 2 * ntab1(ix) + 1) &
      + fy00(ix) * tab(ik (ix) +     ntab (ix)    ) &
      + fy01(ix) * tab(ik (ix) +     ntab (ix) + 1) &
      + fy10(ix) * tab(ik1(ix) +     ntab1(ix)    ) &
      + fy11(ix) * tab(ik1(ix) +     ntab1(ix) + 1) &
        )
     end do
    end if
    ik  = ik  + 3*ntab
    ik1 = ik1 + 3*ntab1
    !---------------------------------------------------------------------------
    if (present(rk)) then
      do ix=1,mx
       rk(ix,iy,iz,1) = exp( &
          f00(ix) * tab(ik (ix))+  f01(ix) * tab(ik (ix) + 1) &
       +  f10(ix) * tab(ik1(ix))+  f11(ix) * tab(ik1(ix) + 1) &
       + fx00(ix) * tab(ik (ix) + 2 * ntab (ix)    ) &
       + fx01(ix) * tab(ik (ix) + 2 * ntab (ix) + 1) &
       + fx10(ix) * tab(ik1(ix) + 2 * ntab1(ix)    ) &
       + fx11(ix) * tab(ik1(ix) + 2 * ntab1(ix) + 1) &
       + fy00(ix) * tab(ik (ix) +     ntab (ix)    ) &
       + fy01(ix) * tab(ik (ix) +     ntab (ix) + 1) &
       + fy10(ix) * tab(ik1(ix) +     ntab1(ix)    ) &
       + fy11(ix) * tab(ik1(ix) +     ntab1(ix) + 1) &
         )
      end do
      if (size(rk,4)>1) then
        do j=2,size(rk,4)
          rk(:,iy,iz,j) = rk(:,iy,iz,j-1)*10.0
        end do
      end if
    end if
    ik  = ik  + 3*ntab
    ik1 = ik1 + 3*ntab1
    !---------------------------------------------------------------------------
    if (present(tt)) then
     do ix=1,mx
      tt(ix,iy,iz) = ( &
         f00(ix) * tab(ik (ix))+  f01(ix) * tab(ik (ix) + 1) &
      +  f10(ix) * tab(ik1(ix))+  f11(ix) * tab(ik1(ix) + 1) &
      + fx00(ix) * tab(ik (ix) + 2 * ntab (ix)    ) &
      + fx01(ix) * tab(ik (ix) + 2 * ntab (ix) + 1) &
      + fx10(ix) * tab(ik1(ix) + 2 * ntab1(ix)    ) &
      + fx11(ix) * tab(ik1(ix) + 2 * ntab1(ix) + 1) &
      + fy00(ix) * tab(ik (ix) +     ntab (ix)    ) &
      + fy01(ix) * tab(ik (ix) +     ntab (ix) + 1) &
      + fy10(ix) * tab(ik1(ix) +     ntab1(ix)    ) &
      + fy11(ix) * tab(ik1(ix) +     ntab1(ix) + 1) &
        )
     end do
    end if
    ik  = ik  + 3*ntab
    ik1 = ik1 + 3*ntab1
    !---------------------------------------------------------------------------
    if (present(ne)) then
      do ix=1,mx
        ne(ix,iy,iz) = exp( &
         f00(ix) * tab(ik (ix))+  f01(ix) * tab(ik (ix) + 1) &
      +  f10(ix) * tab(ik1(ix))+  f11(ix) * tab(ik1(ix) + 1) &
      + fx00(ix) * tab(ik (ix) + 2 * ntab (ix)    ) &
      + fx01(ix) * tab(ik (ix) + 2 * ntab (ix) + 1) &
      + fx10(ix) * tab(ik1(ix) + 2 * ntab1(ix)    ) &
      + fx11(ix) * tab(ik1(ix) + 2 * ntab1(ix) + 1) &
      + fy00(ix) * tab(ik (ix) +     ntab (ix)    ) &
      + fy01(ix) * tab(ik (ix) +     ntab (ix) + 1) &
      + fy10(ix) * tab(ik1(ix) +     ntab1(ix)    ) &
      + fy11(ix) * tab(ik1(ix) +     ntab1(ix) + 1) &
        )
      end do
    end if
    !---------------------------------------------------------------------------
    if (present(src)) then
      do j=1,size(src,4)
      do ix=1,mx
      ik (ix) = ik (ix) + 3*ntab (ix)
      ik1(ix) = ik1(ix) + 3*ntab1(ix)
      src(ix,iy,iz,j) = exp( &
         f00(ix) * tab(ik (ix))+  f01(ix) * tab(ik (ix) + 1) &
      +  f10(ix) * tab(ik1(ix))+  f11(ix) * tab(ik1(ix) + 1) &
      + fx00(ix) * tab(ik (ix) + 2 * ntab (ix)    ) &
      + fx01(ix) * tab(ik (ix) + 2 * ntab (ix) + 1) &
      + fx10(ix) * tab(ik1(ix) + 2 * ntab1(ix)    ) &
      + fx11(ix) * tab(ik1(ix) + 2 * ntab1(ix) + 1) &
      + fy00(ix) * tab(ik (ix) +     ntab (ix)    ) &
      + fy01(ix) * tab(ik (ix) +     ntab (ix) + 1) &
      + fy10(ix) * tab(ik1(ix) +     ntab1(ix)    ) &
      + fy11(ix) * tab(ik1(ix) +     ntab1(ix) + 1) &
        )
      end do
      end do
    end if
  end do
  end do
  if (nbelow_d+nbelow_e > 0) then
    write(io%output,*) &
      'WARNING: number of points below table d & e:', nbelow_d, nbelow_e
  end if
  if (present(d)) then
    deallocate (lnd_loc)
    nullify (lnd_loc)
  end if
  if (present(e).or.present(lne)) then
    deallocate (ee_loc)
    nullify (ee_loc)
  end if
  call trace%end (itimer)
END SUBROUTINE lookup_table

END MODULE eos_mod
