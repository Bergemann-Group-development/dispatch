!===============================================================================
!> EOS module, reading tabular EOS compatible with STAGGER code
!===============================================================================
MODULE eos_mod
  USE io_mod
  USE io_unit_mod
  USE trace_mod
  USE mpi_mod
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
  integer, save:: nbelow=0, ncall=0, mbox
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
    real(8) :: gamma=1.2d0  ! required for syntax in some solvers
    logical :: on = .false.
    logical :: do_rt = .false.
  contains
    procedure:: init
    procedure:: lookup
  end type
  type (eos_t), public:: eos
CONTAINS

!===============================================================================
!> Initialize EOS module (only once)
!===============================================================================
SUBROUTINE init (self)
  class (eos_t):: self
  real, dimension(:,:,:), pointer:: d, ee, pg
  logical, save :: do_rt = .false.
  namelist /eos_params/ do_eos, top, tablefile, do_rt
  !-----------------------------------------------------------------------------
  do_rt = self%do_rt
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
    end if
  end if
  self%on = do_eos
  self%do_rt = do_rt
  !$omp end critical (eos_cr)
  !$omp barrier
END SUBROUTINE init

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup (self, dim, ee, d, lnd, lne, pg, tt, ne, rk, src, gamma)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:)  , pointer, optional :: ee,d,lne,lnd,pg,tt,ne
  real, dimension(:,:,:,:), pointer, optional :: src,rk
  real, optional:: gamma
  !.............................................................................
  real, dimension(:,:,:),   pointer  :: lnd_loc, ee_loc
  integer,         dimension(dim(1)) :: np,np1,ik,ntab,ik1,ntab1
  real(kind=4),    dimension(dim(1)) :: px,py,f00,f01,f10,f11, &
                                        fx00,fx01,fx10,fx11,fy00,fy01,fy10,fy11
  real    :: rhm1,rhm2,drhm,algrk,eek,eek1,py1
  real    :: qx,qy,pxqx,pxpx,qxqx,pypy,pyqy,qyqy,pxpy,qxqy,pxqy,qxpy
  logical :: newtable
  integer :: ix,iy,iz,kee,kee1,j
  integer :: mx,my,mz,mbox
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
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
  do iz=1,mz
  do iy=1,my
    !---------------------------------------------------------------------------
    ! Density index (6 lops)
    !---------------------------------------------------------------------------
    do ix=1,mx
      algrk=1.+(lnd_loc(ix,iy,iz)-rhm1)/drhm
      np(ix)=max0(1,min0(md-1,int(algrk)))
      px(ix)=algrk-np(ix)
    end do
    do ix=1,mx
      !-------------------------------------------------------------------------
      ! Internal energy index -- this loop has indirect referencing (9 flops)
      !-------------------------------------------------------------------------
      ntab(ix)   = mtab(np(ix))
      eek       = 1.+(ee_loc(ix,iy,iz)-eemin(np(ix)))/deetab(np(ix))
      nbelow    = nbelow + merge(1,0,eek<1.0)           ! FIXME
      eek       = max(eek,1.0)                          ! FIXME
      kee       = min(ntab(ix)-1,max(1,floor(eek)))
      py(ix)     = eek-kee
      ik(ix)     = itab(np(ix))+kee-1
      ntab1(ix)  = mtab(np(ix)+1)
      kee1      = kee - nint((eemin(np(ix)+1)-eemin(np(ix))) &
                                            /deetab(np(ix)))
      kee1      = min0(ntab1(ix)-1,max0(1,kee1))
      ik1(ix)    = itab(np(ix)+1)+kee1-1
      !-------------------------------------------------------------------------
      ! The three first lines below may be commented out when we are satisfied
      ! it never panics
      !-------------------------------------------------------------------------
      eek1      = 1.+(ee_loc(ix,iy,iz)-eemin(np(ix)+1))/deetab(np(ix)+1)
      nbelow    = nbelow + merge(1,0,eek1<kee1)         ! FIXME
      eek1      = max(eek1,real(kee1))                  ! FIXME
      py1       = eek1-kee1
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
  if (present(d))   deallocate (lnd_loc)
  if (present(ee)) deallocate (ee_loc)
  call trace%end (itimer)
END SUBROUTINE lookup

END MODULE eos_mod
