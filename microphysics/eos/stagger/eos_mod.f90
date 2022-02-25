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
    real    :: lnd_min, lnd_max, lnd_stp
    real    ::  ee_min,  ee_max,  ee_stp
    real(8) :: gamma=1.4d0
    real    :: kappa
    logical :: on = .false.
    logical :: initialized = .false.
  contains
    procedure:: init
    procedure:: lookup
    procedure:: inverse_lookup
    procedure, private:: lookup_table
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
  real, save :: kappa = 1.0, gamma=1.2d0
  integer:: ir, iostat
  namelist /eos_params/ do_eos, top, tablefile, kappa, verbose, gamma
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%init')
  !----------------------------------------------------------------------------
  ! Note that this is called from setup_mod.f90, while the code is still
  ! single-threaded, so there should be no critical region or io%lock here
  !----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    do_eos = .true.
    top       = '../../'
    tablefile = 'data/eos/stagger/table.dat'
    rewind (io%input)
    read (io%input,eos_params, iostat=iostat)
    call io%namelist_warning ('eos_params', iostat)
    write (io%namelist, eos_params)
    !-----------------------------------------------------------------------------
    ! Read the table file
    !-----------------------------------------------------------------------------
    file = trim(top)//'/'//trim(tablefile)
    open (io_unit%data, file=trim(file), form='unformatted', status='old', &
      action='read')
    read (io_unit%data) md, iupdte, nvar, mbox, eosxmin, dbox, ul, ut, ur, &
      eps, tff, grv, abnd
    njon = nvar-mbox
    allocate (tmean(md), tamp(md), rhm(md), xcorr(md), thmin(md), thmax(md))
    allocate (dth(md), eemin(md), eemax(md), deetab(md), itab(md), mtab(md))
    if (io%do_trace .and. io%master) &
      write (stdout,*) 'eos_t%init:', 'reading 12 times', md, ' values'
    read (io_unit%data) tmean, tamp, rhm, xcorr, thmin, thmax, dth, eemin, &
      eemax, deetab, itab, mtab
    read (io_unit%data) mtable
    allocate (tab(mtable))
    read (io_unit%data) tab
    close (io_unit%data)
      write (stdout,*) trim(tablefile)//': md, mtable =', md, mtable
      allocate (d(1,1,1), ee(1,1,1), pg(1,1,1))
      d = 1.0
      ee  = 5.0
      !call self%lookup (shape(ee), ee=ee, d=d, pg=pg)
      write (stdout,*) 'd, ee, pg =', d, ee, pg
      deallocate (d, ee, pg)
      if (verbose > 0) then
        write (stdout,*) &
        '      density       min(E)      max(E)        min(T)      max(T)'
        do ir=1,md
          write (stdout,'(i4,1p,e12.4,2(2x,2e12.4))') &
            ir, rhm(ir), eemin(ir), eemax(ir), &
            tab(itab(ir)+2*3*mtab(ir)), tab(itab(ir)+(2*3+1)*mtab(ir)-1)
        end do
      end if
  end if
  self%gamma   = gamma
  self%ee_min  = minval(eemin)
  self%ee_max  = maxval(eemax)
  self%ee_stp  = deetab(1)
  self%lnd_min = log(rhm(1))
  self%lnd_max = log(rhm(md))
  self%lnd_stp = (self%lnd_max-self%lnd_min)/(md-1)
  self%on      = do_eos
  self%kappa   = kappa
  self%Initialized=.true.
  call trace%end()
END SUBROUTINE init

!===============================================================================
!> Interface to call with ideal equation of state fall-back
!> FIXME: Check if this is actually used, and replace with generic call
!===============================================================================
SUBROUTINE lookup (self, lnd, d, lne, lnee, e, ee, ss, s, pg, tt, ne, rk, src, &
                   dlnpdE_r, dlnpdlnr_E, gamma, id, e_below, label)
  class (eos_t):: self
  real, dimension(:,:,:), intent(in), pointer, optional :: lnd, d, lne, lnee, e, ee
  real, dimension(:,:,:), intent(in), pointer, optional :: ss, s
  real, dimension(:,:,:),                      optional :: pg, tt, ne
  real, dimension(:,:,:),                      optional :: dlnpdE_r, dlnpdlnr_E
  real, dimension(:,:,:,:),                    optional :: src, rk
  real, optional:: gamma
  integer, optional:: id
  logical, optional:: e_below
  character(len=*), optional::label
  !.............................................................................
  real, dimension(:,:,:), pointer :: lnd_L, ee_L
  real(8):: g1
  integer:: dim(3)
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup', itimer=itimer)
  !-----------------------------------------------------------------------------
  if (.not. self%initialized) call self%init
  if (present(gamma)) then
    self%gamma = gamma
  end if
  g1 = merge(1d0, self%gamma-1d0, self%gamma==1d0)
  !-----------------------------------------------------------------------------
  ! Establish d_L as the density
  !-----------------------------------------------------------------------------
  nullify (lnd_L, ee_L)
  if (present(d)) then
    dim = shape(d)
    allocate (lnd_L(dim(1),dim(2),dim(3)))
    lnd_L = log(d)
  else if (present(lnd)) then
    dim = shape(lnd)
    lnd_L => lnd
  else
    call mpi%abort ('eos_t: missing 1st argument')
  end if
  !-----------------------------------------------------------------------------
  ! Establish ee_L as the internal energy per unit volume
  !-----------------------------------------------------------------------------
  if (present(ee)) then
    !--- energy per unit mass given ---
    ee_L => ee
  else if (present(lnee)) then
    !--- log energy per unit mass given ---
    allocate (ee_L(dim(1),dim(2),dim(3)))
    ee_L = exp(lnee)
  else if (present(lne)) then
    !--- log energy per unit volume given ---
    allocate (ee_L(dim(1),dim(2),dim(3)))
    ee_L = exp(lne - lnd_L)
  else if (present(e)) then
    !--- energy per unit volume given ---
    allocate (ee_L(dim(1),dim(2),dim(3)))
    if (.not. present(d)) then
      ee_L = e/exp(lnd_L)
    else
      ee_L = e/d
    end if
  else if (present(ss).or.present(s)) then
    call mpi%abort('stagger EOS table not set up for entropy')
  else
    call mpi%abort ('eos_t: missing 2nd argument')
  end if
  !-----------------------------------------------------------------------------
  ! Do the table lookup and deallocate locals
  !-----------------------------------------------------------------------------  
  call lookup_table (self, dim, lnd_L, ee_L, pg=pg, tt=tt, ne=ne, rk=rk, src=src, &
                     dlnpdE_r=dlnpdE_r, dlnpdlnr_E=dlnpdlnr_E, gamma=gamma, &
                     id=id, e_below=e_below, label=label)
  if (.not.present(ee)) deallocate( ee_L)
  if (present(d)) deallocate(lnd_L)
  !-----------------------------------------------------------------------------
  call trace%end(itimer)
END SUBROUTINE lookup

!===============================================================================
!> General EOS lookup routine, returning values in optional arguments
!===============================================================================
SUBROUTINE lookup_table (self, dim, lnd, ee, pg, tt, ne, rk, src, &
                         dlnpdE_r, dlnpdlnr_E, gamma, id, e_below, label)
  class (eos_t):: self
  integer:: dim(3)
  real, dimension(:,:,:)            :: lnd, ee
  real, dimension(:,:,:)  , optional:: pg, tt, ne
  real, dimension(:,:,:,:), optional:: src, rk
  real, dimension(:,:,:)  , optional:: dlnpdE_r, dlnpdlnr_E
  real, optional:: gamma
  integer, optional:: id
  logical, optional:: e_below
  character(len=*), optional::label
  !.............................................................................
  integer,         dimension(dim(1)) :: np,np1,ik,ntab,ik1,ntab1
  real(kind=4),    dimension(dim(1)) :: px,py,f00,f01,f10,f11, &
                                        fx00,fx01,fx10,fx11,fy00,fy01,fy10,fy11
  real    :: rhm1,rhm2,drhm,algrk,eek,eek1,py1
  real    :: qx,qy,pxqx,pxpx,qxqx,pypy,pyqy,qyqy,pxpy,qxqy,pxqy,qxpy
  logical :: newtable
  integer :: ix, iy, iz, kee, kee1, j, id_L
  integer :: mx, my, mz, mbox, nbelow_d, nbelow_e
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Allocate temporary arrays
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%lookup_table', itimer=itimer)
  mx=dim(1)
  my=dim(2)
  mz=dim(3)
  newtable = .false.
  !
  rhm1 = log(rhm(1))
  rhm2 = log(rhm(md))
  drhm = (rhm2-rhm1)/(md-1)
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
    !dir$ ivdep
    !$omp simd private (algrk)
    do ix=1,mx
      algrk=1.+(lnd(ix,iy,iz)-rhm1)/drhm
      nbelow_d = nbelow_d + merge(1,0,lnd(ix,iy,iz)<rhm1)
      np(ix)=max0(1,min0(md-1,int(algrk)))
      px(ix)=algrk-np(ix)
    end do
    !dir$ ivdep
    !$omp simd private (eek, kee, kee1)
    do ix=1,mx
      !-------------------------------------------------------------------------
      ! Internal energy index -- this loop has indirect referencing (9 flops)
      !-------------------------------------------------------------------------
      ntab(ix)  = mtab(np(ix))
      eek       = 1.+(ee(ix,iy,iz)-eemin(np(ix)))/deetab(np(ix))
      !nbelow_e  = nbelow_e + merge(1,0,eek<1.0)
      if (lnd(ix,iy,iz)< rhm1) &
        write (stdout,*) 'd below in eos..', &
        ix,iy, iz, lnd(ix,iy,iz), rhm1, np(ix), 'from ',label, id
      if (ee(ix,iy,iz)< eemin(np(ix))) &
        write (stdout,*) 'e below in eos..', &
        ix,iy, iz, ee(ix,iy,iz), np(ix), shape(eemin)
      if (ee(ix,iy,iz)< eemin(np(ix))) &
       write (stdout,*) 'eemin', eemin(np(ix)), 'from ',label, id
      ee(ix,iy,iz) = max (ee(ix,iy,iz), eemin(np(ix)))
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
    !dir$ ivdep
    !$omp simd private (qx, pxqx, pxpx, qxqx, qy, pyqy, pypy, qyqy, &
    !$omp    & pxqy, pxpy, qxqy, qxpy)
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
    ik  = ik  + ntab
    ik1 = ik1 + ntab1
    !---------------------------------------------------------------------------
    ! Pressure derivatives -- just bi-linear interpolation here
    !---------------------------------------------------------------------------
    if (present(dlnpdE_r)) then
      do ix=1,mx
        dlnpdE_r(ix,iy,iz) = ( &
          f00(ix)*tab(ik (ix)) +  f01(ix)*tab(ik (ix)+1) &
        + f10(ix)*tab(ik1(ix)) +  f11(ix)*tab(ik1(ix)+1) &
                           ) / deetab(np(ix))
      end do
    end if
    ik  = ik  + ntab
    ik1 = ik1 + ntab1
    if (present(dlnpdlnr_E)) then
      do ix=1,mx
        dlnpdlnr_E(ix,iy,iz) = ( &
          f00(ix)*tab(ik (ix)) +  f01(ix)*tab(ik (ix)+1) &
        + f10(ix)*tab(ik1(ix)) +  f11(ix)*tab(ik1(ix)+1) &
                          ) / drhm
      end do
    end if
    ik  = ik  + ntab
    ik1 = ik1 + ntab1
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
    if (present(id)) id_L=id
    write (stdout,*) &
      'WARNING: number of points below table d & e:', nbelow_d, nbelow_e, &
      'for ',trim(label), id_L, minval(lnd), minval(ee)
  end if
  if (present(e_below)) then
    if (nbelow_e > 0) then
      e_below = .true.
    else
      e_below = .false.
    end if
  end if
  call trace%end (itimer)
END SUBROUTINE lookup_table

!===============================================================================
!> Interface to call with ideal equation of state fall-back
!> FIXME: Check if this is actually used, and replace with generic call
!===============================================================================
SUBROUTINE inverse_lookup (self, d, lnd, v, name, ss, ee, label)
  class (eos_t):: self
  real, dimension(:,:,:), pointer :: d, lnd, v
  optional                        :: d, lnd, v
  real, dimension(:,:,:), pointer :: ss, ee
  optional                        :: ss, ee
  character(*), optional:: name
  character(*), optional:: label
  !.............................................................................
  real, dimension(:,:,:), pointer :: x_L, y_L, lnd_l, ee_l
  integer:: dim(3)
  integer:: i
  integer, save:: itimer=0
  real:: stefan
  !-----------------------------------------------------------------------------
  call trace%begin ('eos_t%inverse_lookup', itimer=itimer)
  call io%abort ('eos_t%inverse_lookup not ready yet')
  call trace%end(itimer)
END SUBROUTINE inverse_lookup

END MODULE eos_mod
