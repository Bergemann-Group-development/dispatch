! $Id$
!*******************************************************************************
MODULE transfer_m
  integer verbose
  real, save, dimension(20):: time=0.
  integer, save:: ntime=0,itime,nit,ma

CONTAINS
!-------------------------------------------------------------------------------
SUBROUTINE transfer (mx,my,mz,ny,dxdy,dzdy,xm,ym,zm,dx,dz,sx,sz,lnss0,lnrk0,qq0)
  implicit none
  integer             :: mx,my,mz,ny
  real                :: xm(mx),ym(my),ds(my),zm(mz)
  real                :: dxdy,dzdy,dx,dz,sx,sz
  real, dimension(mx,my,mz) :: lnss0,lnrk0,qq0
  real, dimension(mx)       :: rk1,rk2,ss0,ss1,ss2,d2ss1,d2ss2, &
                               dii1,dii2,d1ss1,d1ss2,dtaus,edtau1
  integer, dimension(my)    :: ixma,ixpa,izma,izpa
  real                      :: px0,px1,px2,pz0,pz1,pz2,ds1,ds2,dtau1,dtau2, &
                               dtauc,d1ss,d2ss,px,qx,pz,qz,ex1,ex2
  integer                   :: ix,iy,iz,nitx,nitz,i,ixm,ixp,izm,izp
  real                      :: cpu(2),dtime
  character(len=10)         :: name(20)
  real, allocatable, dimension(:,:,:) :: rk,ss,ex0,dqqp,dqqm,lnrk,lnss, &
                                         qqp,qqm,etp,etm,dtau,dssdtau
!...............................................................................
!
! Compute the maximum ghost zone size needed
!
  if (verbose > 1) print*,'margins:'
  ntime=ntime+1
  itime=1; time(itime) = time(itime)+dtime(cpu); name(itime)='start'
  ma = 0
  ixma=0; ixpa=0; izma=0; izpa=0
  do iy=2,my-1
    px0 = (xm(1)+dxdy*ym(iy-1))/dx 
    px1 = (xm(1)+dxdy*ym(iy  ))/dx
    px2 = (xm(1)+dxdy*ym(iy+1))/dx
    pz0 = (zm(1)+dzdy*ym(iy-1))/dz
    pz1 = (zm(1)+dzdy*ym(iy  ))/dz
    pz2 = (zm(1)+dzdy*ym(iy+1))/dz
    ds(iy) = sqrt((px2-px1)**2+(ym(iy+1)-ym(iy))**2+(pz2-pz1)**2)
    ixma(iy) = floor(px0)-floor(px1)
    ixpa(iy) = floor(px2)-floor(px1)
    izma(iy) = floor(pz0)-floor(pz1)
    izpa(iy) = floor(pz2)-floor(pz1)
    ma = max(ma,max(abs(ixpa(iy)),abs(izpa(iy))))
    if (verbose > 2) print'(5i4)',iy,ixma(iy),ixpa(my),izma(iy),izpa(iy)
  end do
  px1 = (xm(1)+dxdy*ym(1))/dx
  px2 = (xm(1)+dxdy*ym(2))/dx
  pz1 = (zm(1)+dzdy*ym(1))/dz
  pz2 = (zm(1)+dzdy*ym(2))/dz
  ixpa( 1) = floor(px2)-floor(px1)
  izpa( 1) = floor(pz2)-floor(pz1)
  ds(1) = sqrt((px2-px1)**2+(ym(2)-ym(1))**2+(pz2-pz1)**2)
  if (verbose > 2) print'(5i4)',1,ixpa(1),izpa(1)
  px0 = (xm(1)+dxdy*ym(my-1))/dx
  px1 = (xm(1)+dxdy*ym(my  ))/dx
  pz0 = (zm(1)+dzdy*ym(my-1))/dz
  pz1 = (zm(1)+dzdy*ym(my  ))/dz
  ixma(my) = floor(px0)-floor(px1)
  izma(my) = floor(pz0)-floor(pz1)
  if (verbose > 2) print'(5i4)',1,ixma(my),izma(my)
  itime=itime+1; time(itime) = time(itime)+dtime(cpu); name(itime)='ghost'
  if (verbose > 1) print*,'ma =',ma
!
! Allocate scratch variables with ghost zones
!
  if (verbose > 1) print*,'allocate:'
  if (verbose > 1) print'(3(i3,1h:,i2,5x))',1-ma,mx+ma,my,my,1-mz,mz+ma
  allocate (  rk(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! rho*kappa
              ss(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! source
             qqp(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! I+
             qqm(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! I-
             ! etp(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! exp(-tau[+])
             ! etm(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! exp(-tau[-])
            lnss(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! ln(S)
            lnrk(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! ln(rho*kappa)
            dqqp(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! dI+
            dqqm(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! dI-
             ex0(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! exp(-dtau)
            dtau(1-ma:mx+ma,0:my+1,1-ma:mz+ma), &                               ! exp(-dtau)
         dssdtau(1-ma:mx+ma,0:my+1,1-ma:mz+ma))                                 ! dtau
!
! Compute source functions and opacities
!
  lnrk(1:mx,1:my,1:mz) = lnrk0                                                  ! add ghost zones
  lnss(1:mx,1:my,1:mz) = lnss0
  call ghost_fill (mx,my,mz,ma,lnrk)                                            ! fill in ghost zones
  call ghost_fill (mx,my,mz,ma,lnss)
  itime=itime+1; time(itime) = time(itime)+dtime(cpu); name(itime)='cp+fill'
  
  do iz=1,mz
  do iy=1,my
    px1 = (xm(1)+dxdy*ym(iy))/dx
    pz1 = (zm(1)+dzdy*ym(iy))/dz
    px = px1-floor(px1)
    pz = pz1-floor(pz1)
    qx = 1.-px
    qz = 1.-pz
    ss(1:mx,iy,iz) = (lnss(1:mx,iy,iz  )*qx+lnss(2:mx+1,iy,iz  )*px)*qz &       ! bi-linear
                   + (lnss(1:mx,iy,iz+1)*qx+lnss(2:mx+1,iy,iz+1)*px)*pz         ! interpolation
    rk(1:mx,iy,iz) = (lnrk(1:mx,iy,iz  )*qx+lnrk(2:mx+1,iy,iz  )*px)*qz &       ! bi-linear
                   + (lnrk(1:mx,iy,iz+1)*qx+lnrk(2:mx+1,iy,iz+1)*px)*pz         ! interpolation
    ss(1:mx,iy,iz) = exp(ss(1:mx,iy,iz))
    rk(1:mx,iy,iz) = exp(rk(1:mx,iy,iz))
  end do
  end do
  call ghost_fill (mx,my,mz,ma,ss)
  call ghost_fill (mx,my,mz,ma,rk)
  itime=itime+1; time(itime) = time(itime)+dtime(cpu); name(itime)='interp'
!
! Optical depth increments and atennuation factors
!
  nitx = ceiling(dxdy*(ym(my)-ym(1))/sx)
  nitz = ceiling(dzdy*(ym(my)-ym(1))/sz)
  nit = 1+max(nitx,nitz)
  if (verbose > 2) print*,'nit =',nit
  qqp = 0.
  qqp(:,0:1,:) = -ss(:,0:1,:)
  ! etp(:,0:1,:) = 1.
  qqm = 0.
  qqm(:,my:my+1,:) = 0.
  ! etm(:,my:my+1,:) = 1.
  do iz=1,mz
    do iy=1,my-1
      ixp = ixpa(iy)
      izp = izpa(iy)
      do ix=1,mx
        dtau(ix,iy,iz)    = (rk(ix+ixp,iy+1,iz+izp) + rk(ix,iy,iz))*(0.5*ds(iy))
	dssdtau(ix,iy,iz) = (ss(ix+ixp,iy+1,iz+izp) - ss(ix,iy,iz))/dtau(ix,iy,iz)
      end do
      ex0(:,iy,iz) = exp(-dtau(:,iy,iz))
    end do
    ex0(:,0,iz) = 1.
    dtau(:,0,iz) = dtau(:,1,iz)
    dssdtau(:,0,iz) = dssdtau(:,1,iz)
    ex0(:,my:my+1,iz) = 0.
    dtau(:,my,iz) = dtau(:,my-1,iz)
    dssdtau(:,my,iz) = dssdtau(:,my-1,iz)
  end do
  call ghost_fill (mx,my,mz,ma,ex0)
  call ghost_fill (mx,my,mz,ma,dtau)
  call ghost_fill (mx,my,mz,ma,dssdtau)
  do iz=1,mz
!
! Source function 1st and 2nd derivatives
!
    do iy=2,my
      ixm = ixma(iy)
      izm = izma(iy)
      do ix=1,mx
	dtauc = 1./(dtau(ix,iy,iz)+dtau(ix+ixm,iy-1,iz+izm))
        d1ss = (dssdtau(ix,iy,iz)*dtau(ix+ixm,iy-1,iz+izm)+dssdtau(ix+ixm,iy-1,iz+izm)*dtau(ix,iy,iz))
        d2ss = (dssdtau(ix,iy,iz)-dssdtau(ix+ixm,iy-1,iz+izm))*2.*dtauc
        ex1 = 1.-ex0(ix,iy,iz)
        ex2 = ex1 - dtau(ix,iy,iz)*ex0(ix,iy,iz)
	dqqm(ix,iy,iz) = + d1ss*ex1 + d2ss*ex2
        ex1 = 1.-ex0(ix+ixm,iy-1,iz+izm)
        ex2 = ex1 - dtau(ix+ixm,iy-1,iz+izm)*ex0(ix+ixm,iy-1,iz+izm)
	dqqp(ix,iy,iz) = - d1ss*ex1 + d2ss*ex2
      end do
    end do
  end do
  itime=itime+1; time(itime) = time(itime)+dtime(cpu); name(itime)='dqq'
  ! call ghost_fill (mx,my,mz,ma,dqqp)
  ! call ghost_fill (mx,my,mz,ma,dqqm)
!
! Sum contributions and multiply attenuation factors
!
  do i=1,nit
    do iz=1,mz
      do iy=2,my
        ixm = ixma(iy)
        izm = izma(iy)
        do ix=1,mx
          qqp(ix,iy,iz) = qqp(ix+ixm,iy-1,iz+izm)*ex0(ix+ixm,iy-1,iz+izm) + dqqp(ix,iy,iz)
          ! etp(ix+ixm,iy+1,iz+izp) = etp(ix,iy,iz)*ex0(ix,iy,iz)
        end do
      end do
      do iy=my-1,1,-1
        ixp = ixpa(iy)
        izp = izpa(iy)
        do ix=1,mx
	  qqm(ix,iy,iz) = qqm(ix+ixp,iy+1,iz+izp)*ex0(ix,iy,iz) + dqqm(ix,iy,iz)
	  ! etm(ix,iy,iz) = etm(ix+ixp,iy+1,iz+izp)*ex0(ix,iy,iz)
        end do
      end do
    end do
    call ghost_fill (mx,my,mz,ma,qqp)
    call ghost_fill (mx,my,mz,ma,qqm)
  end do
  if (verbose > 2) then
    if (mz == 1) then
      print*,'source(x):'
      do iy=1,my; print'(i4,12f9.3)',iy,ss(1-ma:mx+ma,iy,1); end do
      print*,'opac:'
      do iy=1,my; print'(i4,12f9.3)',iy,rk(1-ma:mx+ma,iy,1); end do
      print*,'dtau:'
      do iy=1,my; print'(i4,12f9.3)',iy,dtau(1-ma:mx+ma,iy,1); end do
      print*,'ex0:'
      do iy=1,my; print'(i4,12f9.3)',iy,ex0(1-ma:mx+ma,iy,1); end do
      print*,'dssdtau:'
      do iy=1,my; print'(i4,12f9.3)',iy,dssdtau(1-ma:mx+ma,iy,1); end do
      print*,'dqqm:'
      do iy=1,my; print'(i4,12f9.3)',iy,dqqm(1-ma:mx+ma,iy,1); end do
      print*,'dqqp:'
      do iy=1,my; print'(i4,12f9.3)',iy,dqqp(1-ma:mx+ma,iy,1); end do
      print*,'qqm:'
      do iy=1,my; print'(i4,12f9.3)',iy,qqm(1-ma:mx+ma,iy,1); end do
      print*,'qqp:'
      do iy=1,my; print'(i4,12f9.3)',iy,qqp(1-ma:mx+ma,iy,1); end do
    else
      print*,'source:'
      do iy=1,my; print'(i4,12f9.3)',iy,ss(1,iy,1-ma:mz+ma); end do
      print*,'opac:'
      do iy=1,my; print'(i4,12f9.3)',iy,rk(1,iy,1-ma:mz+ma); end do
      print*,'dtau:'
      do iy=1,my; print'(i4,12f9.3)',iy,dtau(1,iy,1-ma:mz+ma); end do
      print*,'ex0:'
      do iy=1,my; print'(i4,12f9.3)',iy,ex0(1,iy,1-ma:mz+ma); end do
      print*,'dssdtau:'
      do iy=1,my; print'(i4,12f9.3)',iy,dssdtau(1,iy,1-ma:mz+ma); end do
      print*,'dqqm:'
      do iy=1,my; print'(i4,12f9.3)',iy,dqqm(1,iy,1-ma:mz+ma); end do
      print*,'dqqp:'
      do iy=1,my; print'(i4,12f9.3)',iy,dqqp(1,iy,1-ma:mz+ma); end do
      print*,'qqm:'
      do iy=1,my; print'(i4,12f9.3)',iy,qqm(1,iy,1-ma:mz+ma); end do
      print*,'qqp:'
      do iy=1,my; print'(i4,12f9.3)',iy,qqp(1,iy,1-ma:mz+ma); end do
    end if
  end if
  itime=itime+1; time(itime) = time(itime)+dtime(cpu); name(itime)='qq'
!
! Interpolate back to cartesian mesh
!
  qqp = qqm+qqp
  itime=itime+1; time(itime) = time(itime)+dtime(cpu); name(itime)='sum'
  do iz=1,mz
  do iy=1,my
    px1 = -(xm(1)+dxdy*ym(iy))/dx
    pz1 = -(zm(1)+dzdy*ym(iy))/dz
    px = px1-floor(px1)
    pz = pz1-floor(pz1)
    qx = 1.-px
    qz = 1.-pz
    qq0(1:mx,iy,iz) = (qqp(1:mx,iy,iz  )*qx+qqp(2:mx+1,iy,iz  )*px)*qz &        ! bi-linear
                    + (qqp(1:mx,iy,iz+1)*qx+qqp(2:mx+1,iy,iz+1)*px)*pz          ! interpolation
  end do
  end do
  itime=itime+1; time(itime) = time(itime)+dtime(cpu); name(itime)='interp'
  
  if (verbose > 1) print'(a,i4,5x,20a10)','nit,times:',nit,name(1:itime)
  if (verbose > 1) print'(a,i4,20f10.3)','nit,times:',nit,time(1:itime),sum(time(2:itime))
  if (verbose > 1) print'(a,i4,20f10.3)','nit,times:',nit,time(1:itime)/sum(time(2:itime)+1e-30)
  if (verbose > 1) print'(a,i4,20f10.1)','nit,times:',nit,time(1:itime)*1e9/(2.*ntime*real(mx)*real(my)*real(mz))
  if (verbose > 1) print*,sum(time(2:itime))*1e9/(2.*ntime*real(mx)*real(my)*real(mz)),' ns/pt'
  deallocate (rk,ss,qqp,qqm,lnss,lnrk,dtau,dqqp,dqqm)  
END SUBROUTINE

!-------------------------------------------------------------------------------
SUBROUTINE ghost_fill (mx,my,mz,ma,ss)
  implicit none
  integer mx,my,mz,ma
  real ss(1-ma:mx+ma,0:my+1,1-ma:mz+ma)
!...............................................................................
  ss(1-ma:0,:,:) = ss(mx-ma+1:mx,:,:)
  ss(:,:,1-ma:0) = ss(:,:,mz-ma+1:mz)
  ss(mx+1:mx+ma,:,:) = ss(1:ma,:,:)
  ss(:,:,mz+1:mz+ma) = ss(:,:,1:ma)
  if (verbose > 2) print'(12f10.5)',ss(:,1,1)
END SUBROUTINE
END MODULE
!*******************************************************************************

!-------------------------------------------------------------------------------
PROGRAM test_transfer                                                           ! sets up properties ..
  USE transfer_m
  implicit none                                                                 ! .. similar to full code
  integer mx,my,mz,ny
  real dxdy,dzdy,dx,dy,dz,sx,sy,sz
  real, parameter:: pi=3.14159265
  real, allocatable, dimension(:):: xm,ym,zm
  real, allocatable, dimension(:,:,:):: lnss,lnrk,qq
  integer i,ix,iy,iz
!...............................................................................
  mx=100
  my=100
  mz=100
  !mx=1
  !my=20
  !mz=10
  ny=my
  dx=0.1
  dy=0.05
  dz=0.1

  allocate(xm(mx),ym(my),zm(mz))
  xm = dx*(/((i-mx/2),i=1,mx)/)
  ym = dy*(/((i-my/2),i=1,my)/)
  zm = dz*(/((i-mz/2),i=1,mz)/)
  sx = dx*mx
  sy = dy*my
  sz = dz*mz

  dxdy = 2.1
  dzdy = 2.1
  if (mx==1) dxdy=0.
  if (mz==1) dzdy=0.
  verbose = 2

  allocate (lnss(mx,my,mz),lnrk(mx,my,mz),qq(mx,my,mz))
  do i=0,min(12,max(mx+1,mz+1))
    qq = 0.
    do iy=1,my
      do iz=1,mz
      do ix=1,mx
        lnss(ix,iy,iz) = sin(2.*pi*(xm(ix)+(i+2)*dx+zm(iz))/max(sx,sz))
        lnrk(ix,iy,iz) = 1.+ym(iy)/0.2
      end do
      end do
      if (verbose > 2) print'(i4,5g12.4)',iy,ym(iy), &
        exp(lnss(1,iy,1)),exp(lnss(mx/2,iy,1)),exp(lnrk(1,iy,1))
    end do
    call transfer (mx,my,mz,ny,dxdy,dzdy,xm,ym,zm,dx,dz,sx,sz,lnss,lnrk,qq)
    if (verbose > 0 .and. mz==1) print'(i3,10f9.5)',i,qq(:,2,1)
    if (verbose > 0 .and. mx==1) print'(i3,10f9.5)',i,qq(1,2,:)
    if (verbose > 0) print*,'END'
  end do
  print*,ma,nit,sum(time(2:itime))*1e9/(2.*ntime*real(mx)*real(my)*real(mz)),' ns/pt'

  deallocate(xm,ym,zm,lnss,lnrk,qq)
END
