
module riemann_mod
  USE parameters_mod
  USE io_mod
  USE mpi_mod
  USE trace_mod
  USE timer_mod
  USE eos_mod
  implicit none
  private
  type riemann_t
    integer(8):: ndiff=0, nsolu=0
    real:: max_dlogd=10.
  contains
    procedure:: hlld
    procedure:: hllc
    procedure:: roe_ec
    procedure, nopass:: hlld_2d
    procedure, nopass:: llf
    procedure, nopass:: llf_2d
  end type
  type(riemann_t), public:: rieman
contains

!===============================================================================
!===============================================================================
subroutine hlld (self, qleft, qright, fgdnv, ngrid, detailed_timer)
  IMPLICIT none
  class(riemann_t):: self
  integer:: ngrid, scase
  real(dp),dimension(:,:)::qleft,qright,fgdnv
  LOGICAL:: detailed_timer
  REAL(dp)::SSL,SSR, c_v
  REAL(dp)::re,le,lim, tgl, tgr, tgstarl, tgstarr
  REAL(dp)::rl,pl,sl,ul,ecinl,stotl,ptotl,Bl,Cl,wstarl,vtotBl,emagl,SAL,sqrrstarl,vl,wl
  REAL(dp)::rr,pr,sr,ur,ecinr,stotr,ptotr,Br,Cr,wstarr,vtotBr,emagr,SAR,sqrrstarr,vr,wr
  REAL(dp)::cfastl,rcl,rstarl,el,vdotBstarl,vstarl,stotstarstarl,stotstarl,vdotBl
  REAL(dp)::cfastr,rcr,rstarr,er,vdotBstarr,vstarr,stotstarstarr,stotstarr,vdotBr
  REAL(dp)::ustar,ptotstar,estar,vstarstar,wstarstar
  REAL(dp)::ro,uo,ptoto,so,vo,wo,visc, scond, smax, temp, M
  REAL(dp)::A,sgnm, somax, csl2, csr2, csl, csr, tgo,aa2
  REAL(dp)::alfa_s, alfa_f, alfa_m,tt1,tt2,tt3, tt4, tt5, tt6, bp1, bp2
  REAL(dp)::Bstarl,Bstarr,Cstarl,Cstarr,calfvenl,calfvenr,cslowl, cslowr
  REAL(dp)::Bo,Co,Bstarstar,Cstarstar,vdotBstarstar,vdotBo,emago,msl,msr, mso,mss, mmax
#define INLINE
#ifdef INLINE
  REAL(dp)::d,P,u,v,w,B,C,B2,d2,cf, gamma
#endif
  INTEGER ::ivar,i,icase, ndiff
  logical:: flag, illconditioned
  integer, save:: itimer=0
  !.............................................................................
  if (detailed_timer) &
    call timer%begin('riemann_t%hlld',  itimer=itimer)
  !
  ndiff = 0
!define PRIVATE
#ifdef PRIVATE
  !$omp simd private(sgnm,d,P,u,v,w,A,B,C,B2,d2,cf,Ptotstar,estar, &
  !$omp &  sl,rl,vl,ecinl,emagl,stotl,Ptotl,vdotBl,Bl,Cl,cfastl,SSL, &
  !$omp &  sr,rr,vr,ecinr,emagr,stotr,Ptotr,vdotBr,Br,Cr,cfastr,SSR, &
  !$omp &  rcl,rstarl,el,Bstarl,Cstarl,vstarl,wstarl, tgl, tgr, tgo &
  !$omp &  rcr,rstarr,er,Bstarr,Cstarr,vstarr,wstarr, M,&
  !$omp &  vdotBstarl,stotstarl,sqrrstarl,calfvenl,SAR, c2l, c2r, cl, cr&
  !$omp &  vdotBstarr,stotstarr,sqrrstarr,calfvenr,SAL, &
  !$omp &  vstarstar,wstarstar,Bstarstar,Cstarstar,vdotBstarstar,visc, &
  !$omp &  scond, tt1, tt2, tt3, tt4, tt5, tt6, aa2, alfa_s, alfa_f, bp1,bp2, &
  !$omp &  icase,ro,uo,vo,wo,Bo,Co,Ptoto,so,vdotBo,emago,ustar)
#endif
gamma = eos%gamma
smax=0.0
somax=0.0
scase=0
msl=0.
msr=0.
mso=0.
mss=0.
mmax=0.
fgdnv(:,9) =0.0
  do i=1,ngrid
    !print '(i3,1p,8e10.2,2x,8e10.2)',i,qleft(i,:),qright(i,:)
    !
    ! Enforce continuity of normal component
    A=half*(qleft(i,4)+qright(i,4))
    sgnm=sign(one,A)
    scond=0.0
    qleft(i,4)=A; qright(i,4)=A
    !
    ! Left variables
    rl=qleft(i,1); Pl=qleft(i,2); ul=qleft(i,3); sl=qleft(i,10)
    vl=qleft(i,5); Bl=qleft(i,6); wl=qleft(i,7); Cl=qleft(i,8); tgl=qleft(i,9)
    ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
    emagl = half*(A*A+Bl*Bl+Cl*Cl)
    !etotl = El*rl+ecinl+emagl
    Ptotl = Pl + emagl
  !  ptotl=rl**gamma*exp(sl/rl*(gamma-one))
    vdotBl= ul*A+vl*Bl+wl*cl
    !
    ! Right variables
    rr=qright(i,1); Pr=qright(i,2); ur=qright(i,3); sr=qright(i,10)
    vr=qright(i,5); Br=qright(i,6); wr=qright(i,7); Cr=qright(i,8); tgr=qright(i,9)
    ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
    emagr = half*(A*A+Br*Br+Cr*Cr)
    !etotr = Er*rr+ecinr+emagr
    Ptotr = Pr + emagr
   ! ptotr=rr**gamma*exp(sr/rr*(gamma-one))
    vdotBr= ur*A+vr*Br+wr*Cr
    !
    ! Find the largest eigenvalues in the normal direction to the interface
!#ifdef INLINE
    d=qleft(i,1); P=qleft(i,2); u=qleft(i,3); A=qleft(i,4)
    v=qleft(i,5); B=qleft(i,6); w=qleft(i,7); C=qleft(i,8)
    B2 = A*A+B*B+C*C
    csl2 = gamma*P/d
    csl = sqrt(csl2)
    d2 = half*(B2/d+csl2)
    cf = sqrt( d2 + sqrt(max(d2*d2-csl2*A*A/d,zero)) )
    cslowl = sqrt(max(d2 - sqrt(max(d2*d2-csl2*A*A/d,zero)),zero))
    cfastl=cf
    d=qright(i,1); P=qright(i,2); u=qright(i,3); A=qright(i,4)
    v=qright(i,5); B=qright(i,6); w=qright(i,7); C=qright(i,8)
    B2 = A*A+B*B+C*C
    csr2 = gamma*P/d
    csr=sqrt(csr2)
    d2 = half*(B2/d+csr2)
    cf = sqrt( d2 + sqrt(max(d2*d2-csr2*A*A/d,zero)) )
    cslowr = sqrt(max(d2 - sqrt(max(d2*d2-csr2*A*A/d,zero)),zero))
    cfastr=cf
!#else
!    CALL find_speed_fast(qleft (i,:),cfastl)
!    CALL find_speed_fast(qright(i,:),cfastr)
!#endif
    !
    ! Compute HLL wave speed
    SSL=min(ul,ur)-max(cfastl,cfastr)!*sqrt((gamma-one)/(two*gamma))
    SSR=max(ul,ur)+max(cfastl,cfastr)!*sqrt((gamma-one)/(two*gamma))
    !
    lim = 1.0_dp
    !
    ! Compute lagrangian sound speed
    rcl=rl*(ul-SSL)
    rcr=rr*(SSR-ur)
    !
    ! Compute acoustic star state
    ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
    Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)
    !
    ! Left star region variables
    rstarl=rl*(SSL-ul)/(SSL-ustar)
    estar =rl*(SSL-ul)*(SSL-ustar)-A**2
    el    =rl*(SSL-ul)*(SSL-ul   )-A**2
    !if (dbg) print 1,'hlld4', estar,rl,SL,ul,ustar,A'
!define MERGE
#ifdef MERGE
    illconditioned = abs(estar) < 1d-4*A**2
    vstarl = merge (vl, vl-A*Bl*(ustar-ul)/estar, illconditioned)
    Bstarl = merge (Bl, Bl*el/estar             , illconditioned)
    wstarl = merge (wl, wl-A*Cl*(ustar-ul)/estar, illconditioned)
    cstarl = merge (cl, Cl*el/estar             , illconditioned)
#else
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
#endif
    vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
    tgstarl = tgl*(SSL-ul)/(SSL-ustar)
    stotstarl=sl
    c_v=one/(eos%gamma-one)
    !stotstarl=((sl-eos%s0)/rl*rstarl)*(SSL-ul)/(SSL-ustar)+eos%s0
    !sl = sl + c_v*log(ptotl) + c_v*eos%gamma*log(rl)
    sqrrstarl=sqrt(rstarl)
    calfvenl=abs(A)/sqrrstarl
    SAL=ustar-calfvenl
    !
    ! Right star region variables
    rstarr=rr*(SSR-ur)/(SSR-ustar)
    estar =rr*(SSR-ur)*(SSR-ustar)-A**2
    er    =rr*(SSR-ur)*(SSR-ur   )-A**2
#ifdef MERGE
    illconditioned = abs(estar) < 1d-4*A**2
    vstarr = merge (vr, vr-A*Br*(ustar-ur)/estar, illconditioned)
    Bstarr = merge (Br, Br*er/estar             , illconditioned)
    wstarr = merge (wr, wr-A*Cr*(ustar-ur)/estar, illconditioned)
    cstarr = merge (cr, Cr*er/estar             , illconditioned)
#else
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
#endif
    !
    vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
    tgstarr = tgr*(SSR-ur)/(SSR-ustar)
    stotstarr=sr
    !sr = sr + c_v*log(ptotr) + c_v*eos%gamma*log(rr)
    !stotstarr = ((sr-eos%s0)/rr*rstarr)*(SSR-ur)/(SSR-ustar)+eos%s0
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
    !etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
    !etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)
    stotstarstarl=stotstarl!-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
    stotstarstarr=stotstarr!+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)
    !
    ! Sample the solution at x/t=0
    tt3=0
    tt4=0
    tt2=0
    if(SSL>=0d0)then
       icase=1
       ro=rl
       uo=ul
       vo=vl
       wo=wl
       Bo=Bl
       Co=Cl
       Ptoto=Ptotl
       tgo=tgl
       so=sl
       vdotBo=vdotBl
       emago=emagl
       M = (ul/csl)-(ur/csr)
       ! This works (no temperature multiplier below):
       ! tt1 =(ul-ssl)*log(pl/ro)
       !scond =(ul-ssl)*log(ecinl)!(log(ro*uo*uo*half))
       !works shu-osher
       !scond=ssl*log(pl/ro)*uo
       tt1=ssl*ssl*half*ro*uo!(ecinl-ecinr)
       !scond = ptoto/ro
       scond=0.
       !print *, ptotl, ptotr, rl, rr
       if (M>=1.0) then
      !  call io%abort('case 1 sonic boom')
       ! estimate P2
       tt1=ptoto*(one+ two*eos%gamma/(eos%gamma+one)*(M*M-one))
       ! estimate rho2
       tt2=ro*((eos%gamma+one)*M*M/((eos%gamma-one)*M*M+two))
       tt3=c_v*log(tt1) - c_v*log(ptotl) + c_v*eos%gamma*log(tt2) - c_v*eos%gamma*log(rl)
       tt3=tt3
       !tt3=(M*gamma)/(gamma+one)
       ! Other option
       tt1=(two*gamma*(M*M-one)+(gamma+one)) &
          /(gamma+one)
       tt2=(two*gamma*(M*M-one)+(gamma+one))*((gamma-one)*(M*M-one)+(gamma+one)) &
          /((gamma+one)*(gamma+one)*(M*M-one)+(gamma+one)*(gamma+one))
       !   tt2=(tt2*tgo-tgo)**2/3
       
       !tt3=log(ptotl)-gamma*log(rl)
       !tt4=(gamma+1)*csl*csl/12.*((ul-ur)/ul)**3/tgl
       end if
       tt3=gamma/(gamma-one)*log(tgl/tgr)+log(ptotl/ptotr)
       tt3=tt3*abs(ul-ur)
       if (M>=1.0) then
         tt1=abs(ul-ur)
         tt2=(ur-tt1)
         tt5=(uo-tt1)
         tt4=(gamma+1)*csr*csr/12.*((tt5-tt2)/tt1)**3/tgl
         
       !  print *, 'sc1 tt4', tt3, tt4, tt1, tt2, tt5
       end if
       !tt3=log(ptotl/ptotr)-gamma*log(rl/rr)
       
       !print *, 'stress',ul, ssl, pl+rl*csl*ul
       !print *,'RIEMANN: left state'
       !tt1 =ssl*log(pl/rl)!+(pl+pr)*half*(1/rl-1/rr)
       !scond =(ul-ssl)*(sqrt(uo*uo*ro*half)*(pl/(pl-one))-A*vdotBo)
       !scond =(ul-ssl)*(((ptoto)/(ptoto-one)))*(log(ro*uo*uo))
       !print *, tt1, scond
    else if(SAL>0d0)then
       icase=2
       ro=rstarl
       uo=ustar
       vo=vstarl
       wo=wstarl
       Bo=Bstarl
       Co=Cstarl
       Ptoto=Ptotstar
       so=stotstarl 
       tgo=tgstarl
      M = (uo/csl)-(ur/csr)
       tt1=0
       vdotBo=vdotBstarl
       emago=half*(A*A+Bo*Bo+Co*Co)
       scond=0
       !if (M>=1.0) then
       tt3=c_v*log(ptoto) - c_v*log(ptotl) + c_v*eos%gamma*log(ro) - c_v*eos%gamma*log(rl)
       ! estimate P2
       tt1=ptoto*(one+ two*eos%gamma/(eos%gamma+one)*(M*M-one))
       ! estimate rho2
       tt2=ro*((eos%gamma+one)*M*M/((eos%gamma-one)*M*M+two))
       tt3=c_v*log(tt1) - c_v*log(ptotl) + c_v*eos%gamma*log(tt2) - c_v*eos%gamma*log(rl)
       !tt3=(M*gamma)/(gamma+one)
       !tt2=(tt2*tgo-tgl)**2/3
       tt1=(two*gamma*(M*M-one)+(gamma+one)) &
          /(gamma+one)
       tt2=(two*gamma*(M*M-one)+(gamma+one))*((gamma-one)*(M*M-one)+(gamma+one)) &
          /((gamma+one)*(gamma+one)*(M*M-one)+(gamma+one)*(gamma+one))
       
      ! print *, ptotl, ptotr, rr, rl
      ! tt3=gamma/(gamma-one)*log(tgo/tgl)-log(ptoto/ptotl)
        !print *, 'newcond2', uo,csl,M, tt1, tt2, tt3
       !tt3=tt3
       !tt3=tt3+ecinl-(uo*uo+vo*vo+wo*wo)*ro*half/tgl
     !  tt4=(gamma+1)*csl*csl/12.*((ul-ustar)/ul)**3/tgl
       tt3=0
      ! print *, '2nd', tt3, tt4, tt2
       if (M>=1.0) then
         tt1=abs(ul-ur)
         tt2=(ur-tt1)
         tt5=(uo-tt1)
         tt4=(gamma+1)*csr*csr/12.*((tt5-tt2)/tt1)**3/tgl
         
       !  tt3=log(ptotl/ptotr)-gamma*log(rl/rr)
       end if
       tt3=gamma/(gamma-one)*log(tgo/tgr)+log(ptoto/ptotr)
       tt3=tt3*abs(uo-ur)
       !tt3=log(ptotl/ptotr)-gamma*log(rl/rr)
       !end if
       
    else if(ustar>0d0)then
       icase=3
       ro=rstarl
       uo=ustar
       vo=vstarstar
       wo=wstarstar
       Bo=Bstarstar
       Co=Cstarstar
       Ptoto=Ptotstar
       so=stotstarstarl
       tgo=tgstarl
       M = (uo-ur)/csl
       tt1=0
       scond=0
       if (M>=1.0) then
       tt3=c_v*log(ptoto) - c_v*log(ptotl) + c_v*eos%gamma*log(ro) - c_v*eos%gamma*log(rl)
       tt3=tt3
       end if
       tt3=log(ptotl/ptotr)-gamma*log(rl/rr)
       tt3=0
       vdotBo=vdotBstarstar
       emago=half*(A*A+Bo*Bo+Co*Co)
       tt3=gamma/(gamma-one)*log(tgo/tgr)+log((ptoto+emago)/pr)
       tt3=tt3*abs(uo-ur)!-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)/tgo

      !print *, 'ic3', so, tt3, -sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)/tgo
      ! print *, 'ic3',so, tt3, , (gamma/(gamma-one)*log(tgo/tgr)+log((ptoto+emago)/pr))*abs(uo-ur)
       !print *, 'ic3'
    else if(SAR>0d0)then
       icase=4
       !print *, 'c4'
       ro=rstarr
       uo=ustar
       vo=vstarstar
       wo=wstarstar
       Bo=Bstarstar
       Co=Cstarstar
       Ptoto=Ptotstar
       so=stotstarstarr
       tgo=tgstarr
       vdotBo=vdotBstarstar
       emago=half*(A*A+Bo*Bo+Co*Co)
       M = (ul/csl)-(ur/csr)
       tt1=0
       scond=0
       if (M>=1.0) then
       tt3=c_v*log(ptoto) - c_v*log(ptotr) + c_v*eos%gamma*log(ro) - c_v*eos%gamma*log(rr)
       tt3=tt3
       end if
       tt3=log(ptotl/ptotr)-gamma*log(rl/rr)
       tt3=0
       tt3=gamma/(gamma-one)*log(tgl/tgo)+log(pl/(ptoto+emago))
       tt3=tt3*abs(ul-uo)!+sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)/tgo
    else if (SSR>0d0)then
       icase=5
       ro=rstarr
       uo=ustar
       vo=vstarr
       wo=wstarr
       Bo=Bstarr
       Co=Cstarr
       Ptoto=Ptotstar
       so = stotstarr
       tgo=tgstarr
       vdotBo=vdotBstarr
       emago=half*(A*A+Bo*Bo+Co*Co)
       M=(ul/csl)-(uo/csr)
       scond=0
       tt1=0
       !print *, 'fifth', M, uo, ul, ur,ssr,csr, ssr-ur
       !if (M>=1.0) then
       tt3=c_v*log(ptoto) - c_v*log(ptotr) + c_v*eos%gamma*log(ro) - c_v*eos%gamma*log(rr)
       tt3=uo*tt3
       ! estimate P2
       tt1=ptoto*(one+ two*eos%gamma/(eos%gamma+one)*(M*M-one))
       ! estimate rho2
       tt2=ro*((eos%gamma+one)*M*M/((eos%gamma-one)*M*M+two))
       tt3=tt3
       tt3=(c_v*log(tt1) - c_v*log(ptoto) + c_v*eos%gamma*log(tt2) - c_v*eos%gamma*log(rr))
      ! tt3=(2*gamma)/(gamma+one)
     ! tt4=(gamma+1)*csr*csr/12.*((ur-ustar)/ur)**3/tgr
     ! tt1=(two*gamma*(M*M-one)+(gamma+one)) &
     !    /(gamma+one)
     ! tt2=(two*gamma*(M*M-one)+(gamma+one))*((gamma-one)*(M*M-one)+(gamma+one)) &
     !    /((gamma+one)*(gamma+one)*(M*M-one)+(gamma+one)*(gamma+one))
     ! tt3=gamma/(gamma-one)*log(tt2)+log(tt1)
 !     print *, ptotl, ptotr, rr, rl
      tt3=0
      if (M>=1.0) then
         tt1=abs(ur-ul)
         tt2=(ul-tt1)
         tt5=(uo-tt1)
         tt4=(gamma+1)*csr*csr/12.*((tt5-tt2)/tt1)**3/tgl
         !print *, 'ic5 tt4', tt3, tt4, tt1, tt2, tt5
         tt3=(log(ptoto)-gamma*log(ro))
         tt3=log(ptotl/ptotr)-gamma*log(rl/rr)
       
       end if
       tt3=gamma/(gamma-one)*log(tgl/tgo)+log(ptotl/ptoto)
       tt3=tt3*abs(ul-uo)
       !   tt3=(M*M*gamma)/(gamma+one)
      ! tt2=(tt2*tgo-tgl)**2/3
      ! end if
    else
       icase=6
       ro=rr
       uo=ur
       vo=vr
       wo=wr
       Bo=Br
       Co=Cr
       so=sr 
       tgo=tgr
       Ptoto=Ptotr
       vdotBo=vdotBr
       emago=emagr
       M = (ul/csl)-(ur/csr)
       !Works
       !tt1 =-(ssr-ur)*log(pr/ro)
       !scond =-(ssr-ur)*log(ecinr)!(log(ro*uo*uo*half))
       !scond=ssr*ptoto/ro
       !print *,'RIEMANN: right state'
       scond=0.
       tt1=ssr*ssr*half*ro*uo
       if (M>=1.0 .or. M<=-1.0) then
       ! estimate P2
       tt1=ptoto*(one+ two*eos%gamma/(eos%gamma+one)*(M*M-one))
       ! estimate rho2
       tt2=ro*((eos%gamma+one)*M*M/((eos%gamma-one)*M*M+two))
       tt3=(c_v*log(tt1) - c_v*log(ptoto) + c_v*eos%gamma*log(tt2) - c_v*eos%gamma*log(rr))
       tt3=tt3
       tt3=(log(ptoto)-gamma*log(ro))
       tt4=(gamma+1)*csr*csr/12.*((ur-ul)/ur)**3/tgr
       end if
       tt3=log(ptotl/ptotr)-gamma*log(rl/rr)
       tt3=gamma/(gamma-one)*log(tgl/tgr)+log(ptotl/ptotr)
       tt3=tt3*abs(ul-ur)
       !tt3=(M*gamma)/(gamma+one)
       
       !tt2=(tt2*tgo-tgr)**2/3
       !print *, tt1, scond
       !tt1 =ssr*log(pr/ro)
       !scond =-(ssr-ur)*(sqrt(uo*uo*ro*half)*(pr/(pr-one))-A*vdotBo)
       !scond =-(ssr-ur)*(ptoto/(ptoto-one))*(log(ro*uo*uo))
    end if
   !if (i<25)  print *, i,icase, uo*ro*so, uo, tt3, tt4, ptoto
    !---------------------------------------------------------------------------
    ! Add a diffusive flux in the down grad(ln(rho)) direction, which is only
    ! nonzero when the smallest of the left/right density is within a factor 10
    ! from smallr, and which has at most magnitude of 0.25 times the largest of
    ! the left/right fast mode speeds.
    !---------------------------------------------------------------------------
    ustar = alog(rl/rr)*max(0.0, alog10((10.*smallr)/min(rl,rr)))
    ustar = 0.25*max(cfastr,cfastl)*min(1.0, max(-1.0, ustar))
    uo = uo + ustar

    msl=max(msl,sl)
    msr=max(msr,sr)
    mso=max(mso,so)
    mss=max(mss,M*ro/(eos%gamma-1d0)*(tt3+tt4))
    mmax=max(mmax,M)
   ! print *, icase, uo, ssl, ul-ur, tt3, M, tt4
    if (abs(so*uo)>abs(somax)) somax=so*uo
    !
    ndiff = ndiff + merge(1,0,abs(ustar)>0.0)
    !
    ! Compute the Godunov flux
    fgdnv(i,1) = ro*uo
    !fgdnv(i,2) = so*uo + log(ptoto+ro)+one/tgo *((uo*uo+vo*vo+wo+wo)*half)!scond
    fgdnv(i,2) = uo*ro*(so)!+tt1!ro*log(ptoto/ro**gamma)/(gamma-one)
    fgdnv(i,3) = ro*uo*uo+Ptoto+emago*(lim-one)-A*A*lim
    fgdnv(i,4) = zero
    fgdnv(i,5) = ro*uo*vo-A*Bo*lim
    fgdnv(i,6) = Bo*uo-A*vo
    fgdnv(i,7) = ro*uo*wo-A*Co*lim
    fgdnv(i,8) = Co*uo-A*wo
    fgdnv(i,9) = (tt3)!(abs(tt3)+abs(tt4))!(tt4+tt3)!+tt4!ro*tt3!tt3!tt3!(tt1)!*tgo
   ! print *, i, icase, so, tt3,tt4, uo, ul-ur, M, M*M!,tt4, ro/(eos%gamma-1d0)*(tt3+tt4), ptotl, ptotr, rl, rr, log(ptotl/ptotr),log(rl/rr)
   
    if (abs(scond)>abs(smax)) then
      scase=icase
      smax = scond
    end if
  !  if (M*ro/(eos%gamma-1d0)*(tt3+tt4)<0) then
  !    print *, 'tt3 leq 0', M*ro/(eos%gamma-1d0)*(tt3+tt4)
  !    error stop "end"
  !  end if
  !  print *, icase, so, uo, tt3+tt4, uo*ro/(eos%gamma-1d0)*(tt3+tt4),ro/(eos%gamma-1d0)*(tt3+tt4)
  end do
  !print *, msl, msr, mso, mss, mmax
 ! error stop "end"
!#define ATOMIC
#ifdef ATOMIC
  !$omp atomic
  self%ndiff = self%ndiff + ndiff
  !$omp atomic
  self%nsolu = self%nsolu + ngrid
#endif
  if (detailed_timer) &
    call timer%end (itimer)

end subroutine hlld

!===============================================================================
!===============================================================================
subroutine roe_ec (self, qleft, qright, fgdnv, ngrid, detailed_timer)
  IMPLICIT none
  class(riemann_t):: self
  integer:: ngrid
  real(dp),dimension(:,:)::qleft,qright,fgdnv
  LOGICAL:: detailed_timer
  !
  REAL(dp)::A,sgnm, arbit
  REAL(dp)::z1l, z1r, z5l, z5r, lnz5, lnz1, avz5, avz1, ecin
  REAL(dp)::d, P, u, v, w, B, C, sqd
  REAL(dp)::r_hat, p_hat1, u_hat1, v_hat1, w_hat1
  REAL(dp)::       p_hat2, u_hat2, v_hat2, w_hat2
  REAL(dp)::       A_hat , B_hat , C_hat , AB_hat 
  REAL(dp)::       A_dot , B_dot , C_dot , BC_hat 
  REAL(dp):: alfa_f, alfa_s, beta1, beta2, beta3, aa, aa2, bb2, cc, cs, cf, cf2, cs2
  REAL(dp):: b_perp
  REAL(dp):: matR(8,8), matRT(8,8), Tdiag(8), diss(8), jump(8), inter(8,8),stab(8)
  REAL(dp)::rl,pl,sl,ul,Bl,Cl,wstarl,vl,wl
  REAL(dp)::rr,pr,sr,ur,Br,Cr,wstarr,vr,wr
  REAL(dp)::ro,uo,so,vo,wo, Bo, Co, gamma
  INTEGER ::ivar,i,icase, ndiff, im
  logical:: flag, illconditioned
  integer, save:: itimer=0
  !.............................................................................
  if (detailed_timer) &
    call timer%begin('riemann_t%hlld',  itimer=itimer)
  !
  ndiff = 0
!define PRIVATE
#ifdef PRIVATE
  !$omp simd private(sgnm,d,P,u,v,w,A,B,C,B2,c2,d2,cf,Ptotstar,estar, &
  !$omp &  sl,rl,vl,ecinl,emagl,stotl,Ptotl,vdotBl,Bl,Cl,cfastl,SSL, &
  !$omp &  sr,rr,vr,ecinr,emagr,stotr,Ptotr,vdotBr,Br,Cr,cfastr,SSR, &
  !$omp &  rcl,rstarl,el,Bstarl,Cstarl,vstarl,wstarl, &
  !$omp &  rcr,rstarr,er,Bstarr,Cstarr,vstarr,wstarr, &
  !$omp &  vdotBstarl,stotstarl,sqrrstarl,calfvenl,SAR, &
  !$omp &  vdotBstarr,stotstarr,sqrrstarr,calfvenr,SAL, &
  !$omp &  vstarstar,wstarstar,Bstarstar,Cstarstar,vdotBstarstar,visc, &
  !$omp &  icase,ro,uo,vo,wo,Bo,Co,Ptoto,so,vdotBo,emago,ustar)
#endif
  gamma = eos%gamma
  do i=1,ngrid
    !print '(i3,1p,8e10.2,2x,8e10.2)',i,qleft(i,:),qright(i,:)
    !
    ! Enforce continuity of normal component
    A=half*(qleft(i,4)+qright(i,4))
    sgnm=sign(one,A)
    arbit=1.0_dp/sqrt(2.0_dp)
    qleft(i,4)=A; qright(i,4)=A
    !
    ! Left variables
    rl=qleft(i,1); Pl=qleft(i,2); ul=qleft(i,3); sl=qleft(i,9)/qleft(i,1)*(gamma-one)
    vl=qleft(i,5); Bl=qleft(i,6); wl=qleft(i,7); Cl=qleft(i,8)
    !
    ! Right variables
    rr=qright(i,1); Pr=qright(i,2); ur=qright(i,3); sr=qright(i,9)/qright(i,1)*(gamma-one)
    vr=qright(i,5); Br=qright(i,6); wr=qright(i,7); Cr=qright(i,8)
    ! Compute Z-values
    z1l  = sqrt(rl/Pl)
    z1r  = sqrt(rr/Pr)
    z5l  = sqrt(rl*Pl)
    z5r  = sqrt(rr*Pr)
    lnz1 = log_mean(z1l,z1r)
    avz1 = aver(z1l,z1r)
    lnz5 = log_mean(z5l,z5r)
    avz5 = aver(z5l,z5r)
    ! Compute averaged quantities and products
    u = aver(ul,ur)
    v = aver(vl,vr)
    w = aver(wl,wr)
    d = aver(rl,rr)
    sqd = sqrt(d)
    P = aver(Pl,Pr)
    B = aver(Bl,Br)
    C = aver(Cl,Cr)
    r_hat  = avz1*lnz5
    p_hat1 = avz5/avz1
    p_hat2 = (gamma+1.0_dp)/(2.0_dp*gamma)*lnz5/lnz1+(gamma-1.0_dp)/(2.0_dp*gamma)*avz5/avz1
    u_hat1 = aver(z1l*ul,z1r*ul)/avz1
    v_hat1 = aver(z1l*vl,z1r*vl)/avz1
    w_hat1 = aver(z1l*wl,z1r*wl)/avz1
    z1l = z1l*z1l; z1r=z1r*z1r
    avz1 = aver(z1l,z1r)
    u_hat2 = aver(z1l*ul,z1r*ul)/avz1
    v_hat2 = aver(z1l*vl,z1r*vl)/avz1
    w_hat2 = aver(z1l*wl,z1r*wl)/avz1
    A_hat  = A
    B_hat  = aver(Bl,Br)
    C_hat  = aver(Cl,Cr)
    A_dot  = A * A
    B_dot  = aver(Bl*Bl,Br*Br)
    C_dot  = aver(Cl*Cl,Cr*Cr)
    AB_hat = aver(A*Bl,A*Br)
    BC_hat = aver(A*Cl,A*Cr)
    ! Compute entropy conserving flux
    ro = r_hat*u_hat1
    uo = p_hat1 + r_hat*u_hat1*u_hat1 + 0.5_dp*(A_dot+B_dot+C_dot) - A_dot
    vo = r_hat*u_hat1*v_hat1 - AB_hat
    wo = r_hat*u_hat1*w_hat1 - BC_hat
    so = gamma/(gamma-1.0_dp)*u_hat1*p_hat2 + 0.5_dp*r_hat*u_hat1*&
         (u_hat1*u_hat1 + v_hat1*v_hat1 + w_hat1*w_hat1) &
         + u_hat2*(B_hat*B_hat + C_hat*C_hat) &
         - A_hat *(v_hat2*B_hat + w_hat2*C_hat)
    Bo = u_hat2*B_hat - v_hat2*A
    Co = u_hat2*C_hat - w_hat2*A
    ! Stabilize based on propagating plane-wave solutions
    ! Compute the convenience functions
   ! print *, 'src', ro, uo, vo, wo, so, bo, co
    ecin = (u*u + v*v + w*w)* 0.5_dp
    beta1 = merge(A/(B*C),A/arbit,B*C>0.0_dp)
    beta2 = merge(B/(A*C),B/arbit,A*C>0.0_dp)
    beta3 = merge(C/(A*B),C/arbit,A*B>0.0_dp)
    aa2 = gamma*p/d
    aa  = sqrt(aa2)
    bb2 = A*A+B*B+C*C
    b_perp = sqrt(B*B+C*C)
    cc = sqrt((aa2+bb2)*(aa2+bb2)-4.0_dp*aa2*A*A)
    cf2 = 0.5_dp*((aa2+bb2)+cc)              ! Fast mode wave
    cs2 = 0.5_dp*((aa2+bb2)-cc)              ! Slow mode wave
    cf  = sqrt(cf2)
    cs  = sqrt(cs2)
    alfa_f = (aa2-cs2)/(cf2-cs2)
    alfa_s = (cf2-aa2)/(cf2-cs2)
    ! ---------------------------------------------
    ! Fill in the matrix of eigenvectors
    ! ---------------------------------------------
    ! Magnetoacoustic wave +f
    matR(1,1)=alfa_f*d               
    matR(2,1)=alfa_f*d*(u+cf) 
    matR(3,1)=d*(alfa_f*v - alfa_s*cs*beta2*sgnm)
    matR(4,1)=d*(alfa_f*w - alfa_s*cs*beta3*sgnm)
    matR(5,1)=alfa_f*d*ecin+aa*alfa_s*d*b_perp+(alfa_f*d*aa2/(gamma-1)) &
            + alfa_f*cf*d*u - alfa_s*cs*d*sgnm*(v*beta2+w*beta3)
    matR(6,1)=0.0_dp    
    matR(7,1)=alfa_s*aa*beta2*sqd      
    matR(8,1)=alfa_s*aa*beta3*sqd
    ! Alfven wave +a
    matR(1,2)=0_dp
    matR(2,2)=0_dp
    matR(3,2)=+d**1.5*beta3
    matR(4,2)=-d**1.5*beta2
    matR(5,2)=-d**1.5*(beta2*w-beta3*v)
    matR(6,2)=0_dp
    matR(7,2)=-d*beta3
    matR(8,2)=+d*beta2
    ! Magnetoacoustic wave +s
    matR(1,3)=alfa_s*d
    matR(2,3)=alfa_s*d*(u + cs)
    matR(3,3)=d*(alfa_s*v + alfa_f*cf*beta2*sgnm)
    matR(4,3)=d*(alfa_s*w + alfa_f*cf*beta3*sgnm)
    matR(5,3)=alfa_s*d*ecin-aa*alfa_f*d*b_perp+(alfa_s*d*aa2/(gamma-1)) &
            + alfa_s*cs*d*u + alfa_f*cf*d*sgnm*(v*beta2+w*beta3)
    matR(6,3)=0.0_dp    
    matR(7,3)=-alfa_f*aa*beta2*sqd 
    matR(8,3)=-alfa_f*aa*beta3*sqd 
    ! Entropy wave
    matR(1,4)=1.0_dp
    matR(2,4)=u
    matR(3,4)=v
    matR(4,4)=w
    matR(5,4)=ecin
    matR(6,4)=0.0_dp
    matR(7,4)=0.0_dp
    matR(8,4)=0.0_dp
    ! Divergence wave
    matR(:,5)=0.0_dp
    matR(5,5)=A
    matR(6,5)=1.0_dp
    ! Magnetoacoustic wave -s
    matR(1,6)=matR(1,3)
    matR(2,6)=alfa_s*d*(u - cs)
    matR(3,6)=d*(alfa_s*v - alfa_f*cf*beta2*sgnm)
    matR(4,6)=d*(alfa_s*w - alfa_f*cf*beta3*sgnm)
    matR(5,6)=alfa_s*d*ecin-aa*alfa_f*d*b_perp+(alfa_s*d*aa2/(gamma-1)) &
            - alfa_s*cs*d*u - alfa_f*cf*d*sgnm*(v*beta2+w*beta3)
    matR(6,6)=0.0_dp  
    matR(7,6)=matR(7,3) 
    matR(8,6)=matR(8,3)   
    ! Alfven wave -a
    matR(1,7)= 0_dp
    matR(2,7)= 0_dp
    matR(3,7)= -d**1.5*beta3
    matR(4,7)= +d**1.5*beta2
    matR(5,7)= +d**1.5*(beta2*w-beta3*v)
    matR(6,7)= 0_dp
    matR(7,7)= matR(7,2)
    matR(8,7)= matR(8,2)
    ! Magnetoacoustic wave -f
    matR(1,8)=matR(1,1)
    matR(2,8)=alfa_f*d*(u-cf)
    matR(3,8)=d*(alfa_f*v + alfa_s*cs*beta2*sgnm)
    matR(4,8)=d*(alfa_f*w + alfa_s*cs*beta3*sgnm)
    matR(5,8)=(alfa_f*d*ecin)+aa*alfa_s*d*b_perp+(alfa_f*d*aa2/(gamma-1)) &
            - alfa_f*cf*d*u + alfa_s*cs*d*sgnm*(v*beta2+w*beta3)
    matR(6,8)=0.0_dp
    matR(7,8)=matR(7,1)
    matR(8,8)=matR(8,1)
    ! -----------------------------
    ! Diagonal scaling matrix
    ! -----------------------------
    Tdiag(1) = 1.0_dp/(2.0_dp*d*gamma)
    Tdiag(2) = p/(2.0_dp*d*d*d)
    Tdiag(3) = Tdiag(1)
    Tdiag(4) = d*(gamma-1.0_dp)/gamma
    Tdiag(5) = p/d
    Tdiag(6) = Tdiag(1)
    Tdiag(7) = Tdiag(2)
    Tdiag(8) = Tdiag(1)

    do im=1,8
      matR(:,im)=matR(:,im)*sqrt(Tdiag(im))
    end do
    ! ------------------------------
    ! Transpose the matrix
    ! ------------------------------
    matRT = transpose(matR)
    ! ------------------------------
    ! Select scalar dissipation type
    ! ------------------------------
    !Roe
    diss(1)=abs(u+cf)
    diss(2)=abs(u+A)
    diss(3)=abs(u+cs)
    diss(4)=abs(u)
    diss(5)=abs(u)
    diss(6)=abs(u-cs)
    diss(7)=abs(u-A)
    diss(8)=abs(u-cf)
    if (sqrt(max(0.0,(Pl-Pr)/(Pl+Pr)))>0.5) then
      !LLF
      diss(:) = maxval(diss)
    end if
    do im=1,8
      matR(:,im)=matR(:,im)*diss(im)
    end do
    ! ------------------------------
    ! Compute jump in the entropy vector, if s scaled in negative direction, remove minus in front,
    ! ------------------------------
    jump(1) = ((gamma-sr)/(gamma-1.0_dp) - rr*(ur*ur+vr*vr+wr*wr)/(2.0_dp*pr)) &
            - ((gamma-sl)/(gamma-1.0_dp) - rl*(ul*ul+vl*vl+wl*wl)/(2.0_dp*pl))
    jump(2) = rr*ur/pr - rl*ul/pl
    jump(3) = rr*vr/pr - rl*vl/pl
    jump(4) = rr*wr/pr - rl*wl/pl
    jump(5) = (-rr/pr) - (-rl/pl)
    jump(6) = rr/pr*A - rl/pl*A
    jump(7) = rr/pr*B - rl/pl*B
    jump(8) = rr/pr*C - rl/pl*C
    jump = -jump
    ! -----------------------------
    ! Multiply the matrixes
    ! ----------------------------
    inter = matmul(matR,matRT)
    stab  = matmul(inter,jump)
    !print *, 'flx',fgdnv(i,:)
    !print *, 'stab',0.5*stab
    stab = 0.5_dp*stab
    !
    ! Compute the Godunov flux
    !visc = ro*uo/2*(uo*uo+vo*vo+wo*wo)+uo*Bo*Bo+uo*Co*Co-vo*A*Bo-w*A*Co
    fgdnv(i,1) = ro+stab(1) 
    fgdnv(i,2) = so+stab(5)
    !fgdnv(i,2) = ((gamma*uo*Ptoto)/(gamma-one) + ro*uo/2d0*(uo*uo + vo*vo + wo*wo) + &
    !             uo*(Bo*Bo + Co*Co) - vo*A*Bo - wo*A*Co)
    fgdnv(i,3) = uo+stab(2)
    fgdnv(i,4) = zero
    fgdnv(i,5) = vo+stab(3)
    fgdnv(i,6) = Bo+stab(7)
    fgdnv(i,7) = wo+stab(4)
    fgdnv(i,8) = Co+stab(8)
  end do
!#define ATOMIC
#ifdef ATOMIC
  !$omp atomic
  self%ndiff = self%ndiff + ndiff
  !$omp atomic
  self%nsolu = self%nsolu + ngrid
#endif
  if (detailed_timer) &
    call timer%end (itimer)
  contains
  function log_mean(aL, aR) result (F)
    real(dp):: aL, aR, F
    real(dp):: psi, ff, u, F1
    real(dp):: eps=1e-2
    ! --------------------------------------------------------------------------
    psi = aL/aR
    ff = (psi-1.0_dp)/(psi+1.0_dp)
    u  = ff * ff
    F1 = merge(1.0+u/3.0_dp+u*u/5.0_dp+u*u*u/7.0_dp,log(psi)/2.0_dp/ff,u<eps)
    F  = (aL+aR)/(2.0_dp*F1)
  end function log_mean
  function aver(aL, aR) result (F)
    real(dp):: aL, aR, F
    ! --------------------------------------------------------------------------
    F  = (aL+aR)*0.5_dp
  end function aver
END SUBROUTINE roe_ec

!===============================================================================
!> Solve for the HLLC flux give two HD state vectors (left and right of the
!> interface) stored in fgdnv
!===============================================================================
SUBROUTINE hllc (self, qleft, qrght, fgdnv, ngrid, detailed_timer)
  implicit none
  class(riemann_t):: self
  !-----------------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------------
  real, parameter:: half=0.5, zero=0.0, one=1.0
  LOGICAL:: detailed_timer
  !--- single precision
  REAL(dp),dimension(:,:):: qleft,qrght,fgdnv
  REAL(dp):: ll,rl,sl,pl,ul,vl,wl,ecinl,stotl,ptotl
  REAL(dp):: lr,rr,sr,pr,ur,vr,wr,ecinr,stotr,ptotr
  REAL(dp):: rstarl,stotstarl
  REAL(dp):: rstarr,stotstarr
  REAL(dp):: ro,uo,so,ptoto,stoto
  REAL(dp):: ptotstar
  !--- double precision
  REAL(dp):: cfastl,SSL,rcl
  REAL(dp):: cfastr,SSR,rcr
  REAL(dp):: ustar, gamma
  !---
  INTEGER:: ngrid,i2,i3,dir,n_cheap,n_expensive
  INTEGER:: i,j1,j2,j3,v1,v2,v3
  INTEGER, SAVE:: itimer=0
  n_cheap = 0
  n_expensive = 0
  gamma = eos%gamma
    !---------------------------------------------------------------------------
    ! Loop over 1st index
    !---------------------------------------------------------------------------
    do i=1,ngrid
      !
      ! Left variables
      ll = qleft(i,1)
      pl = qleft(i,2)
      sl = qleft(i,9)
      ul = qleft(i,3)
      vl = qleft(i,5)
      wl = qleft(i,7)
      rl = exp(ll)
      !
      ! Right variables
      lr = qrght(i,1)
      pr = qrght(i,2)
      sr = qrght(i,9)
      ur = qrght(i,3)
      vr = qrght(i,5)
      wr = qrght(i,5)
      rr = exp(lr)
      !
      ! Find the largest eigenvalues in the normal direction to the interface
      cfastl = sqrt(gamma*pl/rl)
      cfastr = sqrt(gamma*pr/rr)
      !
      ! Compute HLL wave speed; SL positive when positive supersonic, SR negative when negative ss
      SSL = min(ul,ur) - max(cfastl,cfastr)
      SSR = max(ul,ur) + max(cfastl,cfastr)
      if (SSL > 0.0) then
        !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'SL case',i
        ro = rl
        uo = ul
        so =sl
        Ptoto = pl
        n_cheap = n_cheap+1
      else if (SSR < 0.0) then
        !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'SR case',i
        ro = rr
        uo = ur
        so = sr
        Ptoto = pr
        n_cheap = n_cheap+1
      else
        n_expensive = n_expensive+1
        !
        ! Compute acoustic star state and lagrangian sound speed
        rcl = rl*(ul-SSL)
        rcr = rr*(SSR-ur)
        Ptotl = pl
        Ptotr = pr
        ustar =    (rcr*ur    + rcl*ul    +   (Ptotl-Ptotr))/(rcr+rcl)
        Ptotstar = (rcr*Ptotl + rcl*Ptotr + rcl*rcr*(ul-ur))/(rcr+rcl)
        !
        ! Sample the solution at x/t=0
        if (ustar > 0.0)then
          !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'rstarl case',i
          ! Left star region variables
          ro = rl*(SSL-ul)/(SSL-ustar)
          so = sl*(SSL-ul)/(SSL-ustar)
          uo = ustar
          Ptoto = Ptotstar
          stotl = sl
          stoto = ((SSL-ul)*stotl)/(SSL-ustar)
        else
          ! Right star region variables
          !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'rstarr case',i
          ro = rr*(SSR-ur)/(SSR-ustar)
          so = sr*(SSR-ur)/(SSR-ustar)
          uo = ustar
          Ptoto = Ptotstar
          stotr = sr
          stoto = ((SSR-ur)*stotr)/(SSR-ustar)
        end if
      end if
      !
      ! Compute the Godunov flux
      fgdnv(i,1) = ro*uo
      fgdnv(i,2) = uo*so
      fgdnv(i,3) = ro*uo*uo + Ptoto
      fgdnv(i,5) = fgdnv(i,1)*merge(vl,vr,uo>0.0)
      fgdnv(i,7) = fgdnv(i,1)*merge(wl,wr,uo>0.0)
    end do
#ifdef DBGC
  stop
#endif
  !-----------------------------------------------------------------------------
  ! Update solver statistics
  !-----------------------------------------------------------------------------
  !$omp atomic update
  timer%n_cheap = timer%n_cheap + n_cheap
  !$omp atomic update
  timer%n_expensive = timer%n_expensive + n_expensive
END SUBROUTINE hllc

!===============================================================================
!> HLLD 2D fluxes
!===============================================================================
SUBROUTINE hlld_2d (qLL,qRL,qLR,qRR,emf,j,k,l0,l1,nvar,xdim,detailed_timer)
  IMPLICIT none
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
  !
  REAL(dp) :: rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR, &
                       ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  REAL(dp) :: qtmp1,qtmp2,qtmp3,qtmp4,qtmp5,qtmp6,qtmp7,qtmp8
  REAL(dp) :: B2, c2, d2, d1, gamma
  integer, save:: itimer=0
  !.............................................................................
  if (detailed_timer) &
    call timer%begin ('riemann_t%hlld_2d', itimer=itimer)
  gamma = eos%gamma
  do l=l0,l1
    rLL=qLL(l,1); pLL=qLL(l,2); uLL=qLL(l,3); vLL=qLL(l,4); ALL=qLL(l,6); BLL=qLL(l,7) ; CLL=qLL(l,8)
    rLR=qLR(l,1); pLR=qLR(l,2); uLR=qLR(l,3); vLR=qLR(l,4); ALR=qLR(l,6); BLR=qLR(l,7) ; CLR=qLR(l,8)
    rRL=qRL(l,1); pRL=qRL(l,2); uRL=qRL(l,3); vRL=qRL(l,4); ARL=qRL(l,6); BRL=qRL(l,7) ; CRL=qRL(l,8)
    rRR=qRR(l,1); pRR=qRR(l,2); uRR=qRR(l,3); vRR=qRR(l,4); ARR=qRR(l,6); BRR=qRR(l,7) ; CRR=qRR(l,8)
    !
    ! Compute 4 fast magnetosonic velocity relative to x direction
    qtmp1=qLL(l,1); qtmp2=qLL(l,2); qtmp7=qLL(l,5); qtmp8=qLL(l,8)
    qtmp3=qLL(l,3); qtmp4=qLL(l,6); qtmp5=qLL(l,4); qtmp6=qLL(l,7)
    !call find_speed_fast(qtmp,cfastLLx)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
    cfastLLx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    qtmp1=qLR(l,1); qtmp2=qLR(l,2); qtmp7=qLR(l,5); qtmp8=qLR(l,8)
    qtmp3=qLR(l,3); qtmp4=qLR(l,6); qtmp5=qLR(l,4); qtmp6=qLR(l,7)
    !call find_speed_fast(qtmp,cfastLRx)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
    cfastLRx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    qtmp1=qRL(l,1); qtmp2=qRL(l,2); qtmp7=qRL(l,5); qtmp8=qRL(l,8)
    qtmp3=qRL(l,3); qtmp4=qRL(l,6); qtmp5=qRL(l,4); qtmp6=qRL(l,7)
    !call find_speed_fast(qtmp,cfastRLx)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
!    print*,'qtmp4,qtmp6,qtmp8',qtmp4,qtmp6,qtmp8
    cfastRLx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    qtmp1=qRR(l,1); qtmp2=qRR(l,2); qtmp7=qRR(l,5); qtmp8=qRR(l,8)
    qtmp3=qRR(l,3); qtmp4=qRR(l,6); qtmp5=qRR(l,4); qtmp6=qRR(l,7)
    qtmp1=max(qtmp1,smallr)
    !call find_speed_fast(qtmp,cfastRRx)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
    cfastRRx=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    ! Compute 4 fast magnetosonic velocity relative to y direction
    qtmp1=qLL(l,1); qtmp2=qLL(l,2); qtmp7=qLL(l,5); qtmp8=qLL(l,8)
    qtmp3=qLL(l,4); qtmp4=qLL(l,7); qtmp5=qLL(l,3); qtmp6=qLL(l,6)
    !call find_speed_fast(qtmp,cfastLLy)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
    cfastLLy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    qtmp1=qLR(l,1); qtmp2=qLR(l,2); qtmp7=qLR(l,5); qtmp8=qLR(l,8)
    qtmp3=qLR(l,4); qtmp4=qLR(l,7); qtmp5=qLR(l,3); qtmp6=qLR(l,6)
    !call find_speed_fast(qtmp,cfastLRy)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
    cfastLRy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    qtmp1=qRL(l,1); qtmp2=qRL(l,2); qtmp7=qRL(l,5); qtmp8=qRL(l,8)
    qtmp3=qRL(l,4); qtmp4=qRL(l,7); qtmp5=qRL(l,3); qtmp6=qRL(l,6)
    !call find_speed_fast(qtmp,cfastRLy)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
    cfastRLy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    qtmp1=qRR(l,1); qtmp2=qRR(l,2); qtmp7=qRR(l,5); qtmp8=qRR(l,8)
    qtmp3=qRR(l,4); qtmp4=qRR(l,7); qtmp5=qRR(l,3); qtmp6=qRR(l,6)
    !call find_speed_fast(qtmp,cfastRRy)
    !
    B2 = qtmp4**2+qtmp6**2+qtmp8**2
    d1 = 1./qtmp1
    c2 = gamma*qtmp2*d1
    d2 = half*(B2*d1+c2)
    cfastRRy=sqrt(d2 + sqrt(max(d2*d2-c2*qtmp4**2*d1,zero)))
    !
    SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
    SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
    SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
    ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
    !
    if(SB>0d0)then
       if(SL>0d0)then
          ELL=uLL*BLL-vLL*ALL
    !
          E=ELL
       else if(SR<0d0)then
          ERL=uRL*BRL-vRL*ARL
    !
          E=ERL
       else
          PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
          PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
          PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
          PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
          rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
          rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
          rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
          rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
          ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
          vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
    !
          rstarLLx=rLL*(SL-uLL)/(SL-ustar)
          rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
          AstarLL=ALL*(SB-vLL)/(SB-vstar)
    !
          rstarLRx=rLR*(SL-uLR)/(SL-ustar)
          rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
          AstarLR=ALR*(ST-vLR)/(ST-vstar)
    !
          rstarRLx=rRL*(SR-uRL)/(SR-ustar)
          rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
          AstarRL=ARL*(SB-vRL)/(SB-vstar)
    !
          rstarRRx=rRR*(SR-uRR)/(SR-ustar)
          rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
          AstarRR=ARR*(ST-vRR)/(ST-vstar)
    !
          calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
          calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
    !
          SAL=min(ustar-calfvenL,zero)
          SAR=max(ustar+calfvenR,zero)
    !
          BstarLL=BLL*(SL-uLL)/(SL-ustar)
          BstarRL=BRL*(SR-uRL)/(SR-ustar)
    !
          EstarLLx=ustar*BstarLL-vLL*ALL
          EstarRLx=ustar*BstarRL-vRL*ARL
    !
          E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
       endif
    else if (ST<0d0)then
       if(SL>0d0)then
          ELR=uLR*BLR-vLR*ALR
    !
          E=ELR
       else if(SR<0d0)then
          ERR=uRR*BRR-vRR*ARR
    !
          E=ERR
       else
          PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
          PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
          PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
          PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
          rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
          rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
          rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
          rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
          ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
          vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
    !
          rstarLLx=rLL*(SL-uLL)/(SL-ustar)
          rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
          AstarLL=ALL*(SB-vLL)/(SB-vstar)
    !
          rstarLRx=rLR*(SL-uLR)/(SL-ustar)
          rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
          AstarLR=ALR*(ST-vLR)/(ST-vstar)
    !
          rstarRLx=rRL*(SR-uRL)/(SR-ustar)
          rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
          AstarRL=ARL*(SB-vRL)/(SB-vstar)
    !
          rstarRRx=rRR*(SR-uRR)/(SR-ustar)
          rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
          AstarRR=ARR*(ST-vRR)/(ST-vstar)
    !
          calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
                       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
          calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
                       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
    !
          SAL=min(ustar-calfvenL,zero)
          SAR=max(ustar+calfvenR,zero)
    !
          BstarLR=BLR*(SL-uLR)/(SL-ustar)
          BstarRR=BRR*(SR-uRR)/(SR-ustar)
    !
          EstarLRx=ustar*BstarLR-vLR*ALR
          EstarRRx=ustar*BstarRR-vRR*ARR
    !
          E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
       endif
    else if(SL>0d0)then
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
       ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
       vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
    !
       rstarLLy=rLL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
       BstarLL=BLL*(SL-uLL)/(SL-ustar)
    !
       rstarLRy=rLR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
       BstarLR=BLR*(SL-uLR)/(SL-ustar)
    !
       rstarRLy=rRL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
       BstarRL=BRL*(SR-uRL)/(SR-ustar)
    !
       rstarRRy=rRR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
       BstarRR=BRR*(SR-uRR)/(SR-ustar)
    !
       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
    !
       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)
    !
       AstarLL=ALL*(SB-vLL)/(SB-vstar)
       AstarLR=ALR*(ST-vLR)/(ST-vstar)
    !
       EstarLLy=uLL*BLL-vstar*AstarLL
       EstarLRy=uLR*BLR-vstar*AstarLR
    !
       E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
    else if(SR<0d0)then
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
       ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
       vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
    !
       rstarLLy=rLL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
       BstarLL=BLL*(SL-uLL)/(SL-ustar)
    !
       rstarLRy=rLR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
       BstarLR=BLR*(SL-uLR)/(SL-ustar)
    !
       rstarRLy=rRL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
       BstarRL=BRL*(SR-uRL)/(SR-ustar)
    !
       rstarRRy=rRR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
       BstarRR=BRR*(SR-uRR)/(SR-ustar)
    !
       calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
                    abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
       calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
                    abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
    !
       SAB=min(vstar-calfvenB,zero)
       SAT=max(vstar+calfvenT,zero)
    !
       AstarRL=ARL*(SB-vRL)/(SB-vstar)
       AstarRR=ARR*(ST-vRR)/(ST-vstar)
    !
       EstarRLy=uRL*BRL-vstar*AstarRL
       EstarRRy=uRR*BRR-vstar*AstarRR
    !
       E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
    else
       PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
       PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
       PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
       PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)
    !
       rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
       rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
       rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
       rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)
    !
       ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
       vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)
    !
       rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
       rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
       rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
    !
       rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
       rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
       rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
    !
       rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
       rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
       rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
    !
       rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
       rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
       rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
    !
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
    !
       AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL); AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
       BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB); BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)
    !
       E=(SAL*SAB*EstarRR-SAL*SAT*EstarRL-SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
            & -SAT*SAB/(SAT-SAB)*(AstarT-AstarB)+SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
    endif
    emf(l,j,k) = E
  end do
  if (detailed_timer) &
    call timer%end (itimer)
END SUBROUTINE hlld_2d

!===============================================================================
!> LLF 2D fluxes
!===============================================================================
SUBROUTINE llf_2d (qLL,qRL,qLR,qRR,emf,j,k,l0,l1,nvar,xdim,detailed_timer)
  IMPLICIT none
  INTEGER, INTENT(IN) :: j,k,l0,l1,nvar,xdim
  REAL(dp), DIMENSION(:,:),INTENT(IN)::qLL,qRL,qLR,qRR
  REAL(dp), DIMENSION(:,:,:):: emf
  LOGICAL, OPTIONAL:: detailed_timer
  ! local variables
  INTEGER:: l
  REAL(dp) :: E, EV, ELL, ERL, ELR, ERR, qleft(8), qright(8), fmean_x(8), fmean_y(8)
  integer, save:: itimer=0
  !.............................................................................
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
END SUBROUTINE llf_2d

!===============================================================================
!> LLF
!===============================================================================
SUBROUTINE llf (qleft, qright, fgdnv, ngrid,  detailed_timer)
  IMPLICIT none
  real(DP), dimension(:,:):: qleft, qright, fgdnv
  integer:: ngrid, l
  logical:: detailed_timer
  !-----------------------------------------------------------------------------
  do l=1,ngrid
    call lax_friedrich (qleft(l,:), qright(l,:), fgdnv(l,:))
  end do
END SUBROUTINE llf

!===============================================================================
!> 1D local Lax-Friedrich Riemann solver
!===============================================================================
SUBROUTINE lax_friedrich (qleft, qright, fgdnv)
  IMPLICIT none
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
  REAL(dp):: ecin,emag,etot,d,u,v,w,A,B,C,P,Ptot,e
  !-----------------------------------------------------------------------------
  ! Local variables
  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4); e=qvar(9)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  ecin = half*(u*u+v*v+w*w)*d
  emag = half*(A*A+B*B+C*C)
  etot = e+ecin+emag
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
  IMPLICIT none
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
  c2 = eos%gamma*P/d
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
  IMPLICIT none
  REAL(dp), DIMENSION(:):: qvar
  REAL(dp):: vel_info
  REAL(dp):: d, P, u, v, w, A, B, C, B2, c2, d2, cf
  !-----------------------------------------------------------------------------
  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
  c2 = eos%gamma*P/d
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(d2**2-c2*A*A/d) )
  vel_info = cf+abs(u)
END SUBROUTINE find_speed_info

!===============================================================================
!===============================================================================
SUBROUTINE find_speed_fast_vec (qvar, vel_info, n)
  !! calculate the fast magnetosonic velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT none
  INTEGER :: n, i ,ib, iv
  REAL(dp),DIMENSION(:,:),target:: qvar
  REAL(dp),DIMENSION(:), pointer:: d, P, u, v, w, A, B, C
  REAL(dp),DIMENSION(n):: vel_info
  REAL(dp),DIMENSION(n):: B2, c2, d2, cf
  !-----------------------------------------------------------------------------
  d=>qvar(:,1); P=>qvar(:,2); u=>qvar(:,3); A=>qvar(:,4)
  v=>qvar(:,5); B=>qvar(:,6); w=>qvar(:,7); C=>qvar(:,8)
  B2 = A*A+B*B+C*C
#ifdef KROME
  c2 = merge(qvar(:,ichem),eos%gamma,chemistry)*P/d
#else
  c2 = eos%gamma*P/d
#endif
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
  vel_info = cf
END SUBROUTINE find_speed_fast_vec

END MODULE riemann_mod
