MODULE stagger_mod2
CONTAINS
!*******************************************************************************
FUNCTION xdn2 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  n = shape(a)
  c = 0.5
  b(1,:,:) = a(1,:,:)
  do iz=1,n(3)
  do iy=1,n(2)
  do ix=2,n(1)
    b(ix,iy,iz) = c*(a(ix  ,iy,iz) + a(ix-1,iy,iz))
  end do
  end do
  end do
END FUNCTION

!*******************************************************************************
FUNCTION ydn2 (a) RESULT (b)
  real, dimension(:,:,:):: a
  real, dimension(size(a,1),size(a,2),size(a,3)):: b
  real    :: c
  integer :: ix, iy, iz, n(3)
  integer, save:: itimer=0
!...............................................................................
  n = shape(a)
  c = 0.5
  b(:,1,:) = a(:,1,:)
  do iz=1,n(3)
  do iy=2,n(2)
  do ix=1,n(1)
    b(ix,iy,iz) = c*(a(ix,iy,iz) + a(ix,iy-1,iz))
  end do
  end do
  end do
END FUNCTION
END MODULE

PROGRAM test
  USE stagger_mod
  USE stagger_mod2
  USE omp_timer_mod
  implicit none
  real, dimension(:,:,:), pointer:: fxy, fyz
  real, dimension(:,:,:), pointer:: Ux, Uy, Uz
  real, dimension(:,:,:), pointer:: fdx, fdy, fdz
  real, dimension(:,:,:), pointer:: w1, w2, w3, w4
  integer, parameter:: m=64, iter=10
  integer, parameter:: mx=m, my=m, mz=m
  integer:: i, ix, iy, iz
  real(8):: wc, mtot
  real(8):: nstagger=4, fl_per_stagger=2.
  real:: xd1, xd2, yd1, yd2

  mtot=(m-1)**3*(iter-1)

  !-----------------------------------------------------------------------------
  ! Function calls
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    if (i==2) wc=wallclock()
    fxy = fxy + (xdn(fdy)*ydn(Ux)+ydn(fdx)*xdn(Uy))
  end do
  wc=wallclock()-wc
  print 1,'stagger operations via function calls      :', &
    wc*1e9/mtot,' ns/upd',mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  1 format(a,2(f7.1,2x,a))
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Local function calls, encouraging compiler to in-line
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    if (i==2) wc=wallclock()
    fxy = fxy + (xdn2(fdy)*ydn2(Ux)+ydn2(fdx)*xdn2(Uy))
  end do
  wc=wallclock()-wc
  print 1,'stagger operations via local function calls:', &
   wc*1e9/mtot,' ns/upd',mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Temporary arrays
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  allocate (w1(mx,my,mz), w2(mx,my,mz), w3(mx,my,mz), w4(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    if (i==2) wc=wallclock()
    w1 = xdn(fdy)
    w2 = ydn(Ux)
    w3 = ydn(fdx)
    w4 = xdn(Uy)
    fxy = fxy + (w1*w2+w3*w4)
  end do
  wc=wallclock()-wc
  print 1,'stagger operations via temporary arrays    :', &
    wc*1e9/mtot,' ns/upd',mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  deallocate (w1, w2, w3, w4)
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Merge the four stagger operations into one loop
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    if (i==2) wc=wallclock()
    fxy(1,:,:)=0.0
    fxy(:,1,:)=0.0
    fxy(:,:,1)=0.0
    do iz=2,mz; do iy=2,my; do ix=2,mx
      xd1 = (fdy(ix,iy,iz)-fdy(ix-1,iy,iz))*0.5
      yd1 = (Ux(ix,iy,iz)-Ux(ix,iy-1,iz))*0.5
      xd2 = (Uy(ix,iy,iz)-Uy(ix-1,iy,iz))*0.5
      yd2 = (fdx(ix,iy,iz)-fdx(ix,iy-1,iz))*0.5
      fxy(ix,iy,iz) = fxy(ix,iy,iz) + xd1*yd1+xd2*yd2
    end do; end do; end do
  end do
  wc=wallclock()-wc
  print 1,'stagger operations merged                  :', &
    wc*1e9/mtot,' ns/upd',mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Eliminate temporary scalars
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    fxy(:,1,:)=0.0
    fxy(:,:,1)=0.0
    if (i==2) wc=wallclock()
    do iz=2,mz
      do iy=2,my
        fxy(1,iy,iz)=0.0
        do ix=2,mx
          fxy(ix,iy,iz) = fxy(ix,iy,iz) &
            + (fdy(ix,iy,iz)-fdy(ix-1,iy,iz))*0.5*( Ux(ix,iy,iz)- Ux(ix,iy-1,iz))*0.5 &
            + ( Uy(ix,iy,iz)- Uy(ix-1,iy,iz))*0.5*(fdx(ix,iy,iz)-fdx(ix,iy-1,iz))*0.5
        end do
      end do
    end do
  end do
  wc=wallclock()-wc
  print 1,'temporary scalars eliminated               :', &
    wc*1e9/mtot,' ns/upd',mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Eliminate temporary scalars
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    fxy(:,1,:)=0.0
    fxy(:,:,1)=0.0
    if (i==2) wc=wallclock()
    do iz=1,mz
      do iy=1,my
        fxy(1,iy,iz)=0.0
        do ix=1,mx
          fxy(ix,iy,iz) = fxy(ix,iy,iz) &
            + (fdy(ix,iy,iz)-fdy(max(1,ix-1),iy,iz))*0.5*( Ux(ix,iy,iz)- Ux(ix,max(1,iy-1),iz))*0.5 &
            + ( Uy(ix,iy,iz)- Uy(max(1,ix-1),iy,iz))*0.5*(fdx(ix,iy,iz)-fdx(ix,max(1,iy-1),iz))*0.5
        end do
      end do
    end do
  end do
  wc=wallclock()-wc
  print 1,'temporary scalars eliminated, min/max index:', &
    wc*1e9/mtot,' ns/upd',mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Eliminate temporary scalars
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    fxy(:,1,:)=0.0
    fxy(:,:,1)=0.0
    if (i==2) wc=wallclock()
    do iz=3,mz-1
      do iy=3,my-1
        fxy(1,iy,iz)=0.0
        do ix=3,mx-1
          fxy(ix,iy,iz) = fxy(ix,iy,iz) &
            + ((fdy(ix,iy,iz)-fdy(ix-1,iy,iz))*0.4 + (fdy(ix+1,iy,iz)-fdy(ix-2,iy,iz))*0.1) &
            * (( Ux(ix,iy,iz)- Ux(ix,iy-1,iz))*0.4 + ( Ux(ix,iy+1,iz)- Ux(ix,iy-2,iz))*0.1) &
            + (( Uy(ix,iy,iz)- Uy(ix-1,iy,iz))*0.4 + ( Uy(ix+1,iy,iz)- Uy(ix-2,iy,iz))*0.1) &
            * ((fdx(ix,iy,iz)-fdx(ix,iy-1,iz))*0.4 + (fdx(ix,iy+1,iz)-fdx(ix,iy-2,iz))*0.1)
        end do
      end do
    end do
  end do
  wc=wallclock()-wc
  print 1,'temporary scalars eliminated; 4th order    :', &
    wc*1e9/mtot,' ns/upd',mtot*(4*6)/(wc*1e9),' Gfl/s'
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Eliminate temporary scalars
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), fyz(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), Uz(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz), fdz(mx,my,mz))
  fxy=0.0; fyz=0.0; Ux=1.1; Uy=Ux-sin(Ux); Uz=Ux+sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy); fdz=fdy-cos(Ux)
  do i=1,iter
    fxy(:,1,:)=0.0
    fxy(:,:,1)=0.0
    fyz(:,1,:)=0.0
    fyz(:,:,1)=0.0
    if (i==2) wc=wallclock()
    do iz=2,mz
      do iy=2,my
        fxy(1,iy,iz)=0.0
        fyz(1,iy,iz)=0.0
        do ix=2,mx
          fxy(ix,iy,iz) = fxy(ix,iy,iz) &
            + (fdy(ix,iy,iz)-fdy(ix-1,iy,iz))*0.5*( Ux(ix,iy,iz)- Ux(ix,iy-1,iz))*0.5 &
            + ( Uy(ix,iy,iz)- Uy(ix-1,iy,iz))*0.5*(fdx(ix,iy,iz)-fdx(ix,iy-1,iz))*0.5
          fyz(ix,iy,iz) = fyz(ix,iy,iz) &
            + (fdz(ix,iy,iz)-fdz(ix,iy-1,iz))*0.5*( Uy(ix,iy,iz)- Uy(ix,iy,iz-1))*0.5 &
            + ( Uz(ix,iy,iz)- Uz(ix,iy-1,iz))*0.5*(fdy(ix,iy,iz)-fdy(ix,iy,iz-1))*0.5
        end do
      end do
    end do
  end do
  wc=wallclock()-wc
  print 1,'two expressions merged                     :', &
    wc*1e9/mtot,' ns/upd',2*mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  deallocate (fxy, Ux, Uy, fdx, fdy)

  !-----------------------------------------------------------------------------
  ! Local functions
  !-----------------------------------------------------------------------------
  allocate (fxy(mx,my,mz), Ux(mx,my,mz), Uy(mx,my,mz), fdx(mx,my,mz), fdy(mx,my,mz))
  fxy = 0.0; Ux=1.1; Uy=Ux-sin(Ux); fdx=Uy+5.; fdy=fdx+cos(Uy)
  do i=1,iter
    if (i==2) wc=wallclock()
    fxy(1,:,:)=0.0
    fxy(:,1,:)=0.0
    fxy(:,:,1)=0.0
    do iz=2,mz; do iy=2,my; do ix=2,mx
      fxy(ix,iy,iz) = fxy(ix,iy,iz) &
        + xdn3(fdy,ix,iy,iz)*ydn3(Ux,ix,iy,iz) + xdn3(Uy,ix,iy,iz)*ydn3(fdx,ix,iy,iz)
    end do; end do; end do
  end do
  wc=wallclock()-wc
  print 1,'local scalar functions                     :', &
    wc*1e9/mtot,' ns/upd',mtot*(1+nstagger*fl_per_stagger)/(wc*1e9),' Gfl/s'
  deallocate (fxy, Ux, Uy, fdx, fdy)
contains
  real function xdn3(f,ix,iy,iz)
  real:: f(:,:,:)
  integer:: ix,iy,iz
  xdn3=0.5*(f(ix,iy,iz)+f(ix-1,iy,iz))
  end function
  real function ydn3(f,ix,iy,iz)
  real:: f(:,:,:)
  integer:: ix,iy,iz
  ydn3=0.5*(f(ix,iy,iz)+f(ix,iy-1,iz))
  end function

END PROGRAM
