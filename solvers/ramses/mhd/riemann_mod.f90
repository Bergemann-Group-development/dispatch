module riemann_mod
  USE amr_parameters
  USE hydro_parameters
  USE const
  USE io_mod
  USE mpi_mod
  USE trace_mod
  implicit none
  private
  type riemann_t
    integer(8):: ndiff=0, nsolu=0
    real:: max_dlogd=10.
  contains
    procedure:: hlld
    procedure, nopass:: hlld_2d
    procedure, nopass:: llf
    procedure, nopass:: llf_2d
  end type
  type(riemann_t), public:: rieman
contains

!===============================================================================
!===============================================================================
subroutine hlld (self, qleft, qright, fgdnv, ngrid, detailed_timer)
  class(riemann_t):: self
  integer:: ngrid
  real(dp),dimension(:,:)::qleft,qright,fgdnv
  LOGICAL, OPTIONAL:: detailed_timer
  REAL(dp)::SL,SR
  REAL(dp)::entho,lim
  REAL(dp)::rl,pl,ul,ecinl,etotl,ptotl,Bl,Cl,wstarl,vtotBl,emagl,SAL,sqrrstarl,vl,wl
  REAL(dp)::rr,pr,ur,ecinr,etotr,ptotr,Br,Cr,wstarr,vtotBr,emagr,SAR,sqrrstarr,vr,wr
  REAL(dp)::cfastl,rcl,rstarl,el,vdotBstarl,vstarl,etotstarstarl,etotstarl,vdotBl
  REAL(dp)::cfastr,rcr,rstarr,er,vdotBstarr,vstarr,etotstarstarr,etotstarr,vdotBr
  REAL(dp)::ustar,ptotstar,estar,vstarstar,wstarstar
  REAL(dp)::ro,uo,ptoto,etoto,vo,wo
  REAL(dp)::A,sgnm
  REAL(dp)::Bstarl,Bstarr,Cstarl,Cstarr,calfvenl,calfvenr
  REAL(dp)::Bo,Co,Bstarstar,Cstarstar,vdotBstarstar,vdotBo,emago
#define INLINE
#ifdef INLINE
  REAL(dp)::d,P,u,v,w,B,C,B2,c2,d2,cf
#endif
  INTEGER ::ivar,i,icase, ndiff
  logical:: flag
  integer, save:: itimer=0
  !.............................................................................
  !call trace%begin('riemann_t%hlld',  itimer=itimer, detailed_timer=detailed_timer)
  !
  ! constants
  if (isothermal) then
    entho = one
  else
    entho = one/max(gamma-one,1e-6)
  end if
  !
  fgdnv      = 0d0
  ndiff = 0
  do i=1,ngrid
    !print '(i3,1p,8e10.2,2x,8e10.2)',i,qleft(i,:),qright(i,:)
    !
    ! Enforce continuity of normal component
    A=half*(qleft(i,4)+qright(i,4))
    sgnm=sign(one,A)
    qleft(i,4)=A; qright(i,4)=A
    !
    ! Left variables
    rl=qleft(i,1); Pl=qleft(i,2); ul=qleft(i,3)
    vl=qleft(i,5); Bl=qleft(i,6); wl=qleft(i,7); Cl=qleft(i,8)
    ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
    emagl = half*(A*A+Bl*Bl+Cl*Cl)
    etotl = Pl*entho+ecinl+emagl
    Ptotl = Pl + emagl
    vdotBl= ul*A+vl*Bl+wl*cl
    !
    ! Right variables
    rr=qright(i,1); Pr=qright(i,2); ur=qright(i,3)
    vr=qright(i,5); Br=qright(i,6); wr=qright(i,7); Cr=qright(i,8)
    ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
    emagr = half*(A*A+Br*Br+Cr*Cr)
    etotr = Pr*entho+ecinr+emagr
    Ptotr = Pr + emagr
    vdotBr= ur*A+vr*Br+wr*Cr
    !
    ! Find the largest eigenvalues in the normal direction to the interface
#ifdef INLINE
    d=qleft(i,1); P=qleft(i,2); u=qleft(i,3); A=qleft(i,4)
    v=qleft(i,5); B=qleft(i,6); w=qleft(i,7); C=qleft(i,8)
    B2 = A*A+B*B+C*C
    c2 = gamma*P/d
    d2 = half*(B2/d+c2)
    cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
    cfastl=cf
    d=qright(i,1); P=qright(i,2); u=qright(i,3); A=qright(i,4)
    v=qright(i,5); B=qright(i,6); w=qright(i,7); C=qright(i,8)
    B2 = A*A+B*B+C*C
    c2 = gamma*P/d
    d2 = half*(B2/d+c2)
    cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
    cfastr=cf
#else
    CALL find_speed_fast(qleft (i,:),cfastl)
    CALL find_speed_fast(qright(i,:),cfastr)
#endif
    !
    ! Compute HLL wave speed
    SL=min(ul,ur)-max(cfastl,cfastr)
    SR=max(ul,ur)+max(cfastl,cfastr)
    !
    lim = 1.0_dp
    !
    ! Compute lagrangian sound speed
    rcl=rl*(ul-SL)
    rcr=rr*(SR-ur)
    !
    ! Compute acoustic star state
    ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
    Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)
    !
    ! Left star region variables
    rstarl=rl*(SL-ul)/(SL-ustar)
    estar =rl*(SL-ul)*(SL-ustar)-A**2
    el    =rl*(SL-ul)*(SL-ul   )-A**2
    !if (dbg) print 1,'hlld4', estar,rl,SL,ul,ustar,A
    if(abs(estar) < 1d-4*A**2)then
       vstarl=vl
       Bstarl=Bl
       wstarl=wl
       Cstarl=Cl
    else
       vstarl=vl-A*Bl*(ustar-ul)/estar
       Bstarl=Bl*el/estar
       wstarl=wl-A*Cl*(ustar-ul)/estar
       Cstarl=Cl*el/estar
    endif
    vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
    etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
    sqrrstarl=sqrt(rstarl)
    calfvenl=abs(A)/sqrrstarl
    SAL=ustar-calfvenl
    !
    ! Right star region variables
    rstarr=rr*(SR-ur)/(SR-ustar)
    estar =rr*(SR-ur)*(SR-ustar)-A**2
    er    =rr*(SR-ur)*(SR-ur   )-A**2
    if(abs(estar) < 1d-4*A**2)then
       vstarr=vr
       Bstarr=Br
       wstarr=wr
       Cstarr=Cr
    else
       vstarr=vr-A*Br*(ustar-ur)/estar
       Bstarr=Br*er/estar
       wstarr=wr-A*Cr*(ustar-ur)/estar
       Cstarr=Cr*er/estar
    endif
    !
    vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
    etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
    sqrrstarr=sqrt(rstarr)
    calfvenr=abs(A)/sqrrstarr
    SAR=ustar+calfvenr
    !
    ! Double star region variables
    vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
    wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
    Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
    Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
    vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
    etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
    etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)

    ! Sample the solution at x/t=0
    if(SL>0d0)then
       icase=1
       ro=rl
       uo=ul
       vo=vl
       wo=wl
       Bo=Bl
       Co=Cl
       Ptoto=Ptotl
       etoto=etotl
       vdotBo=vdotBl
       emago=emagl
    else if(SAL>0d0)then
       icase=2
       ro=rstarl
       uo=ustar
       vo=vstarl
       wo=wstarl
       Bo=Bstarl
       Co=Cstarl
       Ptoto=Ptotstar
       etoto=etotstarl
       vdotBo=vdotBstarl
       emago=half*(A*A+Bo*Bo+Co*Co)
    else if(ustar>0d0)then
       icase=3
       ro=rstarl
       uo=ustar
       vo=vstarstar
       wo=wstarstar
       Bo=Bstarstar
       Co=Cstarstar
       Ptoto=Ptotstar
       etoto=etotstarstarl
       vdotBo=vdotBstarstar
       emago=half*(A*A+Bo*Bo+Co*Co)
    else if(SAR>0d0)then
       icase=4
       ro=rstarr
       uo=ustar
       vo=vstarstar
       wo=wstarstar
       Bo=Bstarstar
       Co=Cstarstar
       Ptoto=Ptotstar
       etoto=etotstarstarr
       vdotBo=vdotBstarstar
       emago=half*(A*A+Bo*Bo+Co*Co)
    else if (SR>0d0)then
       icase=5
       ro=rstarr
       uo=ustar
       vo=vstarr
       wo=wstarr
       Bo=Bstarr
       Co=Cstarr
       Ptoto=Ptotstar
       etoto=etotstarr
       vdotBo=vdotBstarr
       emago=half*(A*A+Bo*Bo+Co*Co)
    else
       icase=6
       ro=rr
       uo=ur
       vo=vr
       wo=wr
       Bo=Br
       Co=Cr
       Ptoto=Ptotr
       etoto=etotr
       vdotBo=vdotBr
       emago=emagr
    end if
    ! Add a diffusive flux in the down grad(ln(rho)) direction, which is only nonzero when the
    ! smallest of the left/right density is within a factor 10 from smallr, and which has at most
    ! magnitude of 0.25 times the largest of the left/right fast mode speeds.
    ustar = alog(rl/rr)*max(0.0, alog10((10.*smallr)/min(rl,rr)))
    ustar = 0.25*max(cfastr,cfastl)*min(1.0, max(-1.0, ustar))
    uo = uo + ustar
    ndiff = ndiff + merge(1,0,abs(ustar)>0.0)
    !
    ! Compute the Godunov flux
    fgdnv(i,1) = ro*uo
    fgdnv(i,2) = (etoto+Ptoto)*uo-A*vdotBo
    fgdnv(i,3) = ro*uo*uo+Ptoto+emago*(lim-1.0_dp)-A*A*lim
    fgdnv(i,4) = zero
    fgdnv(i,5) = ro*uo*vo-A*Bo*lim
    fgdnv(i,6) = Bo*uo-A*vo
    fgdnv(i,7) = ro*uo*wo-A*Co*lim
    fgdnv(i,8) = Co*uo-A*wo
  end do
  !$omp atomic
  self%ndiff = self%ndiff + ndiff
  !$omp atomic
  self%nsolu = self%nsolu + ngrid
  !call trace%end (itimer, detailed_timer=detailed_timer)
end subroutine hlld

!===============================================================================
!> HLLD 2D fluxes
!===============================================================================
SUBROUTINE hlld_2d (qLL,qRL,qLR,qRR,emf,j,k,l0,l1,nvar,xdim,detailed_timer)
  INTEGER, INTENT(IN) :: j,k,l0,l1,nvar,xdim
  REAL(dp), DIMENSION(:,:),INTENT(IN)::qLL,qRL,qLR,qRR
  REAL(dp), DIMENSION(:,:,:):: emf
  LOGICAL, OPTIONAL:: detailed_timer
  ! local variables
  INTEGER:: l
  REAL(dp) :: E,ELL,ERL,ELR,ERR,S,SL,SR,SB,ST,SAL,SAR,SAT,SAB
  REAL(dp) :: cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  REAL(dp) :: calfvenR,calfvenL,calfvenT,calfvenB
  REAL(dp) :: vLLx,vRLx,vLRx,vRRx,vLLy,vRLy,vLRy,vRRy
  REAL(dp) :: PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  REAL(dp) :: ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  REAL(dp) :: rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  REAL(dp) :: EstarLLx,EstarLRx,EstarRLx,EstarRRx,EstarLLy,EstarLRy,EstarRLy,EstarRRy,EstarLL,EstarLR,EstarRL,EstarRR
  REAL(dp) :: AstarT,AstarB,BstarR,BstarL

  REAL(dp) :: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR, &
                       ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  REAL(dp), DIMENSION(1:nvar) :: qtmp
  REAL(dp) :: B2, c2, d2, d1
  !integer, save:: itimer(10)=0
  integer, save:: itimer=0
  !.............................................................................
  !call trace%begin ('riemann_t%hlld_2d', itimer=itimer, detailed_timer=detailed_timer)
  do l=l0,l1
    !call trace%begin ('hlld_2d find_fast_speed', itimer=itimer(1))
    rLL=qLL(l,1); pLL=qLL(l,2); uLL=qLL(l,3); vLL=qLL(l,4); ALL=qLL(l,6); BLL=qLL(l,7) ; CLL=qLL(l,8)
    rLR=qLR(l,1); pLR=qLR(l,2); uLR=qLR(l,3); vLR=qLR(l,4); ALR=qLR(l,6); BLR=qLR(l,7) ; CLR=qLR(l,8)
    rRL=qRL(l,1); pRL=qRL(l,2); uRL=qRL(l,3); vRL=qRL(l,4); ARL=qRL(l,6); BRL=qRL(l,7) ; CRL=qRL(l,8)
    rRR=qRR(l,1); pRR=qRR(l,2); uRR=qRR(l,3); vRR=qRR(l,4); ARR=qRR(l,6); BRR=qRR(l,7) ; CRR=qRR(l,8)

    ! Compute 4 fast magnetosonic velocity relative to x direction
    qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
    qtmp(3)=qLL(l,3); qtmp(4)=qLL(l,6); qtmp(5)=qLL(l,4); qtmp(6)=qLL(l,7)
    !call find_speed_fast(qtmp,cfastLLx)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
    cfastLLx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
    qtmp(3)=qLR(l,3); qtmp(4)=qLR(l,6); qtmp(5)=qLR(l,4); qtmp(6)=qLR(l,7)
    !call find_speed_fast(qtmp,cfastLRx)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
    cfastLRx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
    qtmp(3)=qRL(l,3); qtmp(4)=qRL(l,6); qtmp(5)=qRL(l,4); qtmp(6)=qRL(l,7)
    !call find_speed_fast(qtmp,cfastRLx)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
!    print*,'qtmp(4),qtmp(6),qtmp(8)',qtmp(4),qtmp(6),qtmp(8)
    cfastRLx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
    qtmp(3)=qRR(l,3); qtmp(4)=qRR(l,6); qtmp(5)=qRR(l,4); qtmp(6)=qRR(l,7)
    qtmp(1)=max(qtmp(1),smallr)
    !call find_speed_fast(qtmp,cfastRRx)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
    cfastRRx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    ! Compute 4 fast magnetosonic velocity relative to y direction
    qtmp(1)=qLL(l,1); qtmp(2)=qLL(l,2); qtmp(7)=qLL(l,5); qtmp(8)=qLL(l,8)
    qtmp(3)=qLL(l,4); qtmp(4)=qLL(l,7); qtmp(5)=qLL(l,3); qtmp(6)=qLL(l,6)
    !call find_speed_fast(qtmp,cfastLLy)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
    cfastLLy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    qtmp(1)=qLR(l,1); qtmp(2)=qLR(l,2); qtmp(7)=qLR(l,5); qtmp(8)=qLR(l,8)
    qtmp(3)=qLR(l,4); qtmp(4)=qLR(l,7); qtmp(5)=qLR(l,3); qtmp(6)=qLR(l,6)
    !call find_speed_fast(qtmp,cfastLRy)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
    cfastLRy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    qtmp(1)=qRL(l,1); qtmp(2)=qRL(l,2); qtmp(7)=qRL(l,5); qtmp(8)=qRL(l,8)
    qtmp(3)=qRL(l,4); qtmp(4)=qRL(l,7); qtmp(5)=qRL(l,3); qtmp(6)=qRL(l,6)
    !call find_speed_fast(qtmp,cfastRLy)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
    cfastRLy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    qtmp(1)=qRR(l,1); qtmp(2)=qRR(l,2); qtmp(7)=qRR(l,5); qtmp(8)=qRR(l,8)
    qtmp(3)=qRR(l,4); qtmp(4)=qRR(l,7); qtmp(5)=qRR(l,3); qtmp(6)=qRR(l,6)
    !call find_speed_fast(qtmp,cfastRRy)

    B2 = qtmp(4)**2+qtmp(6)**2+qtmp(8)**2
    d1 = 1./qtmp(1)
    c2 = gamma*qtmp(2)*d1
    d2 = half*(B2*d1+c2)
    cfastRRy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp(4)**2*d1,zero)))

    SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
    SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
    SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
    ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
    !call trace%end (itimer(1))

    if(SB>0d0)then
       if(SL>0d0)then
          !call trace%begin ('hlld_2d(2)', itimer=itimer(2))
          ELL=uLL*BLL-vLL*ALL

          E=ELL
          !call trace%end (itimer(2))
       else if(SR<0d0)then
          !call trace%begin ('hlld_2d(3)', itimer=itimer(3))
          ERL=uRL*BRL-vRL*ARL

          E=ERL
          !call trace%end (itimer(3))
       else
          !call trace%begin ('hlld_2d(4)', itimer=itimer(4))
          PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
          PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
          PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
          PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)

          rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
          rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
          rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
          rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)

          ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
          vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)

          rstarLLx=rLL*(SL-uLL)/(SL-ustar)
          rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
          AstarLL=ALL*(SB-vLL)/(SB-vstar)

          rstarLRx=rLR*(SL-uLR)/(SL-ustar)
          rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
          AstarLR=ALR*(ST-vLR)/(ST-vstar)

          rstarRLx=rRL*(SR-uRL)/(SR-ustar)
          rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
          AstarRL=ARL*(SB-vRL)/(SB-vstar)

          rstarRRx=rRR*(SR-uRR)/(SR-ustar)
          rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
          AstarRR=ARR*(ST-vRR)/(ST-vstar)

          calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
          calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)

          SAL=min(ustar-calfvenL,zero)
          SAR=max(ustar+calfvenR,zero)

          BstarLL=BLL*(SL-uLL)/(SL-ustar)
          BstarRL=BRL*(SR-uRL)/(SR-ustar)

          EstarLLx=ustar*BstarLL-vLL*ALL
          EstarRLx=ustar*BstarRL-vRL*ARL

          E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
          !call trace%end (itimer(4))
       endif
    else if (ST<0d0)then
       if(SL>0d0)then
          !call trace%begin ('hlld_2d(5)', itimer=itimer(5))
          ELR=uLR*BLR-vLR*ALR

          E=ELR
          !call trace%end (itimer(5))
       else if(SR<0d0)then
          !call trace%begin ('hlld_2d(6)', itimer=itimer(6))
          ERR=uRR*BRR-vRR*ARR

          E=ERR
          !call trace%end (itimer(6))
       else
          !call trace%begin ('hlld_2d(7)', itimer=itimer(7))
          PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
          PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
          PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
          PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)

          rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
          rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
          rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
          rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)

          ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
          vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)

          rstarLLx=rLL*(SL-uLL)/(SL-ustar)
          rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
          AstarLL=ALL*(SB-vLL)/(SB-vstar)

          rstarLRx=rLR*(SL-uLR)/(SL-ustar)
          rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
          AstarLR=ALR*(ST-vLR)/(ST-vstar)

          rstarRLx=rRL*(SR-uRL)/(SR-ustar)
          rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
          AstarRL=ARL*(SB-vRL)/(SB-vstar)

          rstarRRx=rRR*(SR-uRR)/(SR-ustar)
          rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
          AstarRR=ARR*(ST-vRR)/(ST-vstar)

          calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
          calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)

          SAL=min(ustar-calfvenL,zero)
          SAR=max(ustar+calfvenR,zero)

          BstarLR=BLR*(SL-uLR)/(SL-ustar)
          BstarRR=BRR*(SR-uRR)/(SR-ustar)

          EstarLRx=ustar*BstarLR-vLR*ALR
          EstarRRx=ustar*BstarRR-vRR*ARR

          E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
          !call trace%end (itimer(7))
       endif
    else if(SL>0d0)then
       !call trace%begin ('hlld_2d(8)', itimer=itimer(8))
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)

       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)

       ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
       vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)

       rstarLLy=rLL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
       BstarLL=BLL*(SL-uLL)/(SL-ustar)

       rstarLRy=rLR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
       BstarLR=BLR*(SL-uLR)/(SL-ustar)

       rstarRLy=rRL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
       BstarRL=BRL*(SR-uRL)/(SR-ustar)

       rstarRRy=rRR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
       BstarRR=BRR*(SR-uRR)/(SR-ustar)

       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)

       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)

       AstarLL=ALL*(SB-vLL)/(SB-vstar)
       AstarLR=ALR*(ST-vLR)/(ST-vstar)

       EstarLLy=uLL*BLL-vstar*AstarLL
       EstarLRy=uLR*BLR-vstar*AstarLR

       E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
       !call trace%end (itimer(8))
    else if(SR<0d0)then
       !call trace%begin ('hlld_2d(9)', itimer=itimer(9))
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)

       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)

       ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
       vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)

       rstarLLy=rLL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
       BstarLL=BLL*(SL-uLL)/(SL-ustar)

       rstarLRy=rLR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
       BstarLR=BLR*(SL-uLR)/(SL-ustar)

       rstarRLy=rRL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
       BstarRL=BRL*(SR-uRL)/(SR-ustar)

       rstarRRy=rRR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
       BstarRR=BRR*(SR-uRR)/(SR-ustar)

       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)

       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)

       AstarRL=ARL*(SB-vRL)/(SB-vstar)
       AstarRR=ARR*(ST-vRR)/(ST-vstar)

       EstarRLy=uRL*BRL-vstar*AstarRL
       EstarRRy=uRR*BRR-vstar*AstarRR

       E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
       !call trace%end (itimer(9))
    else
       !call trace%begin ('hlld_2d(10)', itimer=itimer(10))
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)

       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)

       ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
       vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)

       rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
       rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)

       rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
       rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)

       rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
       rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)

       rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
       rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)

       calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                    abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
       calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                    abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
       SAL=min(ustar-calfvenL,zero)
       SAR=max(ustar+calfvenR,zero)
       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)
       EstarLL =ustar*BstarLL-vstar*AstarLL
       EstarRL =ustar*BstarRL-vstar*AstarRL
       EstarLR =ustar*BstarLR-vstar*AstarLR
       EstarRR =ustar*BstarRR-vstar*AstarRR

       AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL); AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
       BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB); BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)

       E=(SAL*SAB*EstarRR-SAL*SAT*EstarRL-SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
            & -SAT*SAB/(SAT-SAB)*(AstarT-AstarB)+SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
       !call trace%end (itimer(10))
    endif
    if (xdim == 1) then
       emf(l,j,k) = E
    else if (xdim == 2) then
       emf(k,l,j) = E
    else if (xdim == 3) then
       emf(j,k,l) = E
    endif
  end do
  !call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE hlld_2d

!===============================================================================
!> LLF 2D fluxes
!===============================================================================
SUBROUTINE llf_2d (qLL,qRL,qLR,qRR,emf,j,k,l0,l1,nvar,xdim,detailed_timer)
  INTEGER, INTENT(IN) :: j,k,l0,l1,nvar,xdim
  REAL(dp), DIMENSION(:,:),INTENT(IN)::qLL,qRL,qLR,qRR
  REAL(dp), DIMENSION(:,:,:):: emf
  LOGICAL, OPTIONAL:: detailed_timer
  ! local variables
  INTEGER:: l
  REAL(dp) :: E, EV, ELL, ERL, ELR, ERR, qleft(8), qright(8), fmean_x(8), fmean_y(8)
  integer, save:: itimer=0
  !.............................................................................
  !call trace%begin ('riemann_t%hlld_2d', itimer=itimer, detailed_timer=detailed_timer)
  !
  do l=l0,l1
    ! vx*by - vy*bx at the four edge centers
    ELL = qLL(l,3)*qLL(l,7) - qLL(l,4)*qLL(l,6)
    ERL = qRL(l,3)*qRL(l,7) - qRL(l,4)*qRL(l,6)
    ELR = qLR(l,3)*qLR(l,7) - qLR(l,4)*qLR(l,6)
    ERR = qRR(l,3)*qRR(l,7) - qRR(l,4)*qRR(l,6)
    !
    ! find the average value of E
    Ev = forth*(ELL+ERL+ELR+ERR)
    !
    ! call the first solver in the x direction
    ! density
    qleft (1) = half*(qLL(l,1)+qLR(l,1))
    qright(1) = half*(qRR(l,1)+qRL(l,1))
    !
    ! pressure
    qleft (2) = half*(qLL(l,2)+qLR(l,2))
    qright(2) = half*(qRR(l,2)+qRL(l,2))
    !
    ! vt1 becomes normal velocity
    qleft (3) = half*(qLL(l,3)+qLR(l,3))
    qright(3) = half*(qRR(l,3)+qRL(l,3))
    !
    ! bt1 becomes normal magnetic field
    qleft (4) = half*(qLL(l,6)+qLR(l,6))
    qright(4) = half*(qRR(l,6)+qRL(l,6))
    !
    ! vt2 becomes transverse velocity field
    qleft (5) = half*(qLL(l,4)+qLR(l,4))
    qright(5) = half*(qRR(l,4)+qRL(l,4))
    !
    ! bt2 becomes transverse magnetic field
    qleft (6) = half*(qLL(l,7)+qLR(l,7))
    qright(6) = half*(qRR(l,7)+qRL(l,7))
    !
    ! velocity component perp. to the plane is now transverse
    qleft (7) = half*(qLL(l,5)+qLR(l,5))
    qright(7) = half*(qRR(l,5)+qRL(l,5))
    !
    ! magnetic field component perp. to the plane is now transverse
    qleft (8) = half*(qLL(l,8)+qLR(l,8))
    qright(8) = half*(qRR(l,8)+qRL(l,8))
    !
    CALL lax_friedrich(qleft,qright,fmean_x)
    !
    ! call the second solver in the y direction
    ! density
    qleft (1) = half*(qLL(l,1)+qRL(l,1))
    qright(1) = half*(qRR(l,1)+qLR(l,1))
    !
    ! pressure
    qleft (2) = half*(qLL(l,2)+qRL(l,2))
    qright(2) = half*(qRR(l,2)+qLR(l,2))
    !
    ! vt2 becomes normal velocity
    qleft (3) = half*(qLL(l,4)+qRL(l,4))
    qright(3) = half*(qRR(l,4)+qLR(l,4))
    !
    ! bt2 becomes normal magnetic field
    qleft (4) = half*(qLL(l,7)+qRL(l,7))
    qright(4) = half*(qRR(l,7)+qLR(l,7))
    !
    ! vt1 becomes transverse velocity field
    qleft (5) = half*(qLL(l,3)+qRL(l,3))
    qright(5) = half*(qRR(l,3)+qLR(l,3))
    !
    ! bt1 becomes transverse magnetic field
    qleft (6) = half*(qLL(l,6)+qRL(l,6))
    qright(6) = half*(qRR(l,6)+qLR(l,6))
    !
    ! velocity component perp. to the plane is now transverse
    qleft (7) = half*(qLL(l,5)+qRL(l,5))
    qright(7) = half*(qRR(l,5)+qLR(l,5))
    !
    ! magnetic field component perp. to the plane is now transverse
    qleft (8) = half*(qLL(l,8)+qRL(l,8))
    qright(8) = half*(qRR(l,8)+qLR(l,8))
    !
    CALL lax_friedrich(qleft,qright,fmean_y)
    !
    E = Ev + (fmean_x(6) - fmean_y(6))
    if (xdim == 1) then
      emf(l,j,k) = E
    else if (xdim == 2) then
      emf(k,l,j) = E
    else if (xdim == 3) then
      emf(j,k,l) = E
    endif
  end do
  !call trace%end (itimer, detailed_timer=detailed_timer)
END SUBROUTINE llf_2d

!===============================================================================
!> LLF
!===============================================================================
SUBROUTINE llf (qleft, qright, fgdnv, ngrid,  detailed_timer)
  real(DP), dimension(:,:):: qleft, qright, fgdnv
  integer:: ngrid, l
  logical, optional:: detailed_timer
  !-----------------------------------------------------------------------------
  do l=1,ngrid
    call lax_friedrich (qleft(l,:), qright(l,:), fgdnv(l,:))
  end do
END SUBROUTINE llf

!===============================================================================
!> 1D local Lax-Friedrich Riemann solver
!===============================================================================
SUBROUTINE lax_friedrich (qleft, qright, fgdnv)
  REAL(dp), DIMENSION(:):: qleft, qright, fgdnv
  REAL(dp), DIMENSION(size(qleft)):: fleft, fright, fmean
  REAL(dp), DIMENSION(size(qleft)):: uleft ,uright, udiff
  REAL(dp):: vleft, vright, bx_mean
  !-----------------------------------------------------------------------------
  ! Enforce continuity of normal component
  bx_mean = half*(qleft(4)+qright(4))
  qleft (4) = bx_mean
  qright(4) = bx_mean
  !
  CALL find_mhd_flux (qleft ,uleft ,fleft )
  CALL find_mhd_flux (qright,uright,fright)
  !
  ! find the mean flux
  fmean =  half * ( fright + fleft )
  !
  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_info (qleft ,vleft )
  CALL find_speed_info (qright,vright)
  !
  ! difference between the 2 states
  udiff = half * ( uright - uleft )
  !
  ! the local Lax-Friedrich flux
  fgdnv = fmean - MAX(vleft,vright) * udiff
END SUBROUTINE lax_friedrich

!===============================================================================
!===============================================================================
SUBROUTINE find_mhd_flux (qvar, cvar, ff)
  !! compute the 1D MHD fluxes from the conservative variables
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1, Btransverse1, Vtransverse2, Btransverse2
  IMPLICIT NONE
  INTEGER :: i , ivar, ib
  REAL(dp), DIMENSION(:):: qvar, cvar, ff
  REAL(dp):: ecin,emag,etot,d,u,v,w,A,B,C,P,Ptot,entho
  !-----------------------------------------------------------------------------
  ! Local variables
  entho = one/(gamma-one)
  if (isothermal) then
    entho = one
  else
    entho = one/max(gamma-one,1e-6)
  end if
  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  ecin = half*(u*u+v*v+w*w)*d
  emag = half*(A*A+B*B+C*C)
  etot = P*entho+ecin+emag
  Ptot = P + emag
  !
  ! Compute conservative variables
  cvar(1) = d
  cvar(2) = etot
  cvar(3) = d*u
  cvar(4) = A
  cvar(5) = d*v
  cvar(6) = B
  cvar(7) = d*w
  cvar(8) = C
  !
  ! Compute fluxes
  ff(1) = d*u
  ff(2) = (etot+Ptot)*u-A*(A*u+B*v+C*w)
  ff(3) = d*u*u+Ptot-A*A
  ff(4) = zero
  ff(5) = d*u*v-A*B
  ff(6) = B*u-A*v
  ff(7) = d*u*w-A*C
  ff(8) = C*u-A*w
END SUBROUTINE find_mhd_flux

!===============================================================================
!===============================================================================
SUBROUTINE find_speed_fast(qvar,vel_info)
  !! calculate the fast magnetosonic velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  REAL(dp), DIMENSION(:):: qvar
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  d=max(d,smallr)
  P=max(P,smallr*smallc*smallc)
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
  vel_info = cf
END SUBROUTINE find_speed_fast

!===============================================================================
!===============================================================================
SUBROUTINE find_speed_info (qvar, vel_info)
  !! calculate the fastest velocity at which information is exchanged
  !! at the interface
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
  REAL(dp), DIMENSION(:):: qvar
  REAL(dp):: vel_info
  REAL(dp):: d, P, u, v, w, A, B, C, B2, c2, d2, cf
  !-----------------------------------------------------------------------------
  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(d2**2-c2*A*A/d) )
  vel_info = cf+abs(u)
END SUBROUTINE find_speed_info

END MODULE riemann_mod
