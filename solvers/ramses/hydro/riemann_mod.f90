module riemann_mod
  use amr_parameters
  use hydro_parameters
  use const
  USE io_mod
  USE trace_mod
  USE mpi_mod
  implicit none
  private
  type riemann_t
  contains
    procedure, nopass:: approx
    procedure, nopass:: acoustic
    procedure, nopass:: llf
    procedure, nopass:: hll
    procedure, nopass:: hllc
    procedure, nopass:: hllc_a
    procedure, nopass:: hllc_v
  end type
  type(riemann_t), public:: rieman
contains

!===============================================================================
!===============================================================================
subroutine approx(qleft,qright,fgdnv,ngrid,nvar)
  integer::ngrid,nvar
  real(dp),dimension(1:ngrid,1:nvar)::qleft,qright,qgdnv,fgdnv

  ! local arrays
  real(dp),dimension(1:ngrid)::rl   ,ul   ,pl   ,cl
  real(dp),dimension(1:ngrid)::rr   ,ur   ,pr   ,cr
  real(dp),dimension(1:ngrid)::ro   ,uo   ,po   ,co
  real(dp),dimension(1:ngrid)::rstar,ustar,pstar,cstar
  real(dp),dimension(1:ngrid)::wl   ,wr   ,wo
  real(dp),dimension(1:ngrid)::sgnm ,spin ,spout,ushock
  real(dp),dimension(1:ngrid)::frac ,delp ,pold
  integer ,dimension(1:ngrid)::ind  ,ind2

  ! local variables
  real(dp)::smallp, gamma6, ql, qr, usr, usl, wwl, wwr, smallpp, entho, etot
  integer ::i, j, n, iter, n_new

  ! constants
  smallp = smallc**2/gamma
  smallpp = smallr*smallp
  gamma6 = (gamma+one)/(two*gamma)
  entho = one/(gamma-one)

  ! Pressure, density and velocity
  do i=3,ngrid-1
     rl(i)=MAX(qleft (i,1),smallr)
     ul(i)=    qleft (i,2)
     pl(i)=MAX(qleft (i,3),rl(i)*smallp)
     rr(i)=MAX(qright(i,1),smallr)
     ur(i)=    qright(i,2)
     pr(i)=MAX(qright(i,3),rr(i)*smallp)
  end do

  ! Lagrangian sound speed
  do i=3,ngrid-1
     cl(i) = gamma*pl(i)*rl(i)
     cr(i) = gamma*pr(i)*rr(i)
  end do

  ! First guess
  do i=3,ngrid-1
     wl(i)= sqrt(cl(i)); wr(i) = sqrt(cr(i))
     pstar(i) = ((wr(i)*pl(i)+wl(i)*pr(i))+wl(i)*wr(i)*(ul(i)-ur(i)))/(wl(i)+wr(i))
     pstar(i) = MAX(pstar(i),0.0_dp)
     pold (i)= pstar(i)
  end do
  n = ngrid
  do i = 1,n
     ind(i)=i
  end do

  ! Newton-Raphson iterations to find pstar at the required accuracy
  ! for a two-shock Riemann problem
  do iter = 1,niter_riemann
     do i=1,n
        wwl=sqrt(cl(ind(i))*(one+gamma6*(pold(i)-pl(ind(i)))/pl(ind(i))))
        wwr=sqrt(cr(ind(i))*(one+gamma6*(pold(i)-pr(ind(i)))/pr(ind(i))))
        ql=two*wwl**3/(wwl**2+cl(ind(i)))
        qr=two*wwr**3/(wwr**2+cr(ind(i)))
        usl=ul(ind(i))-(pold(i)-pl(ind(i)))/wwl
        usr=ur(ind(i))+(pold(i)-pr(ind(i)))/wwr
        delp(i)=MAX(qr*ql/(qr+ql)*(usl-usr),-pold(i))
     end do
     do i=1,n
        pold(i)=pold(i)+delp(i)
     end do
     ! Convergence indicator
     do i=1,n
        uo(i)=ABS(delp(i)/(pold(i)+smallpp))
     end do
     n_new=0
     do i=1,n
        if(uo(i)>1.d-06)then
           n_new=n_new+1
           ind2(n_new)=ind (i)
           po  (n_new)=pold(i)
        end if
     end do
     j=n_new
     do i=1,n
        if(uo(i)<=1.d-06)then
           n_new=n_new+1
           ind2(n_new)=ind (i)
           po  (n_new)=pold(i)
        end if
     end do
     ind (1:n)=ind2(1:n)
     pold(1:n)=po  (1:n)
     n=j
  end do

  ! Star region pressure
  ! for a two-shock Riemann problem
  do i=3,ngrid-1
     pstar(ind(i))=pold(i)
  end do
  do i=3,ngrid-1
     wl(i) = sqrt(cl(i)*(one+gamma6*(pstar(i)-pl(i))/pl(i)))
     wr(i) = sqrt(cr(i)*(one+gamma6*(pstar(i)-pr(i))/pr(i)))
  end do

  ! Star region velocity
  ! for a two shock Riemann problem
  do i=3,ngrid-1
     ustar(i) = half*(ul(i) + (pl(i)-pstar(i))/wl(i) + &
          &           ur(i) - (pr(i)-pstar(i))/wr(i) )
  end do

  ! Left going or right going contact wave
  do i=3,ngrid-1
     sgnm(i) = sign(one,ustar(i))
  end do

  ! Left or right unperturbed state
  do i=3,ngrid-1
     if(sgnm(i)==one)then
        ro(i) = rl(i)
        uo(i) = ul(i)
        po(i) = pl(i)
        wo(i) = wl(i)
     else
        ro(i) = rr(i)
        uo(i) = ur(i)
        po(i) = pr(i)
        wo(i) = wr(i)
     end if
  end do
  do i=3,ngrid-1
     co(i) = max(smallc,sqrt(abs(gamma*po(i)/ro(i))))
  end do

  ! Star region density
  do i=3,ngrid-1
     if(pstar(i)>= po(i))then
        ! Shock
        rstar(i) = ro(i)/(one+ro(i)*(po(i)-pstar(i))/wo(i)**2)
     else
        ! Rarefaction
        rstar(i) = ro(i)*(pstar(i)/po(i))**(one/gamma)
     end if
  end do

  do i=3,ngrid-1
     ! Prevent vacuum formation in star region
     rstar(i) = max(rstar(i),smallr)
     ! Star region sound speed
     cstar(i) = sqrt(abs(gamma*pstar(i)/rstar(i)))
     cstar(i) = max(cstar(i),smallc)
     ! Compute rarefaction head and tail speed
     spout(i) = co   (i)-sgnm(i)*uo   (i)
     spin (i) = cstar(i)-sgnm(i)*ustar(i)
     ! Compute shock speed
     ushock(i) = wo(i)/ro(i)-sgnm(i)*uo(i)
  end do

  do i=3,ngrid-1
     if(pstar(i)>=po(i))then
        spout(i)=ushock(i)
        spin (i)=spout (i)
     end if
  end do

  ! Sample the solution at x/t=0
  do i=3,ngrid-1
     if(spout(i)<=zero)then
        qgdnv(i,1) = ro(i)
        qgdnv(i,2) = uo(i)
        qgdnv(i,3) = po(i)
     else if(spin(i)>=zero)then
        qgdnv(i,1) = rstar(i)
        qgdnv(i,2) = ustar(i)
        qgdnv(i,3) = pstar(i)
     else
        frac(i)=spout(i)/(spout(i)-spin(i))
        qgdnv(i,2) = frac(i)*ustar(i) + (one - frac(i))*uo(i)
        qgdnv(i,3) = frac(i)*pstar(i) + (one - frac(i))*po(i)
        qgdnv(i,1) = ro(i)*(qgdnv(i,3)/po(i))**(one/gamma)
     end if
  end do

  ! Passive scalars
  do n = 4,nvar
     do i=3,ngrid-1
        if(sgnm(i)==one)then
           qgdnv(i,n) = qleft(i,n)
        else
           qgdnv(i,n) = qright(i,n)
        end if
     end do
  end do

  ! Compute fluxes
  do i = 3,ngrid-1
     fgdnv(i,1) = qgdnv(i,1)*qgdnv(i,2)  ! Mass density
     fgdnv(i,2) = qgdnv(i,3)+qgdnv(i,1)*qgdnv(i,2)**2  ! Normal momentum
     etot = qgdnv(i,3)*entho + half*qgdnv(i,1)*qgdnv(i,2)**2
     etot = etot + half*qgdnv(i,1)*qgdnv(i,4)**2
     etot = etot + half*qgdnv(i,1)*qgdnv(i,5)**2
     fgdnv(i,3) = qgdnv(i,2)*(etot+qgdnv(i,3))     ! Total energy
  end do
  ! Other advected quantities
  do n = 4, nvar
     do i = 3,ngrid-1
        fgdnv(i,n) = fgdnv(i,1)*qgdnv(i,n)
     end do
  end do
end subroutine approx

!===============================================================================
!===============================================================================
subroutine acoustic(qleft,qright,fgdnv,ngrid,nvar)
  integer::ngrid,nvar
  real(dp),dimension(1:ngrid,1:nvar)::qleft,qright,qgdnv,fgdnv

  ! local variables
  integer::i,n
  real(dp)::smallp, entho, etot

  ! local arrays
  real(dp),dimension(1:ngrid)::rl   ,ul   ,pl   ,cl
  real(dp),dimension(1:ngrid)::rr   ,ur   ,pr   ,cr
  real(dp),dimension(1:ngrid)::ro   ,uo   ,po   ,co
  real(dp),dimension(1:ngrid)::rstar,ustar,pstar,cstar
  real(dp),dimension(1:ngrid)::wl   ,wr   ,wo
  real(dp),dimension(1:ngrid)::sgnm ,spin ,spout,ushock
  real(dp),dimension(1:ngrid)::frac

  ! constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  ! Initial states pressure, density and velocity
  do i=3,ngrid-1
     rl(i)=max(qleft (i,1),smallr)
     ul(i)=    qleft (i,2)
     pl(i)=max(qleft (i,3),rl(i)*smallp)
     rr(i)=max(qright(i,1),smallr)
     ur(i)=    qright(i,2)
     pr(i)=max(qright(i,3),rr(i)*smallp)
  end do

  ! Acoustic star state
  do i=3,ngrid-1
     cl(i) = sqrt(gamma*pl(i)/rl(i))
     cr(i) = sqrt(gamma*pr(i)/rr(i))
     wl(i) = cl(i)*rl(i)
     wr(i) = cr(i)*rr(i)
     pstar(i) = ((wr(i)*pl(i)+wl(i)*pr(i))+wl(i)*wr(i)*(ul(i)-ur(i))) &
          &   / (wl(i)+wr(i))
     ustar(i) = ((wr(i)*ur(i)+wl(i)*ul(i))+(pl(i)-pr(i))) &
          &   / (wl(i)+wr(i))
!!$     pstar(i) = MAX(pstar(i),zero)
  end do

  ! Left going or right going contact wave
  do i=3,ngrid-1
     sgnm(i) = sign(one,ustar(i))
  end do

  ! Left or right unperturbed state
  do i=3,ngrid-1
     if(sgnm(i)==one)then
        ro(i) = rl(i)
        uo(i) = ul(i)
        po(i) = pl(i)
        wo(i) = wl(i)
        co(i) = cl(i)
     else
        ro(i) = rr(i)
        uo(i) = ur(i)
        po(i) = pr(i)
        wo(i) = wr(i)
        co(i) = cr(i)
     end if
  end do

  ! Star region density and sound speed
  do i=3,ngrid-1
     rstar(i) = ro(i)+(pstar(i)-po(i))/co(i)**2
     rstar(i) = max(rstar(i),smallr)
     cstar(i) = sqrt(abs(gamma*pstar(i)/rstar(i)))
     cstar(i) = max(cstar(i),smallc)
  end do

  ! Head and tail speed of rarefaction
  do i=3,ngrid-1
     spout(i) = co   (i)-sgnm(i)*uo   (i)
     spin (i) = cstar(i)-sgnm(i)*ustar(i)
  end do

  ! Shock speed
  do i=3,ngrid-1
     ushock(i) = half*(spin(i)+spout(i))
     ushock(i) = max(ushock(i),-sgnm(i)*ustar(i))
  end do
  do i=3,ngrid-1
     if(pstar(i)>=po(i))then
        spout(i)=ushock(i)
        spin (i)=spout (i)
     end if
  end do

  ! Sample the solution at x/t=0
  do i=3,ngrid-1
     if(spout(i)<zero)then      ! Initial state
        qgdnv(i,1) = ro(i)
        qgdnv(i,2) = uo(i)
        qgdnv(i,3) = po(i)
     else if(spin(i)>=zero)then  ! Star region
        qgdnv(i,1) = rstar(i)
        qgdnv(i,2) = ustar(i)
        qgdnv(i,3) = pstar(i)
     else                        ! Rarefaction
        frac(i) = spout(i)/(spout(i)-spin(i))
        qgdnv(i,1) = frac(i)*rstar(i) + (one - frac(i))*ro(i)
        qgdnv(i,2) = frac(i)*ustar(i) + (one - frac(i))*uo(i)
        qgdnv(i,3) = frac(i)*pstar(i) + (one - frac(i))*po(i)
     end if
  end do

  ! Passive scalars
  do n = 4,nvar
     do i=3,ngrid-1
        if(sgnm(i)==one)then
           qgdnv(i,n) = qleft (i,n)
        else
           qgdnv(i,n) = qright(i,n)
        end if
     end do
  end do

  ! Compute fluxes
  do i = 3,ngrid-1
     fgdnv(i,1) = qgdnv(i,1)*qgdnv(i,2)  ! Mass density
     fgdnv(i,2) = qgdnv(i,3)+qgdnv(i,1)*qgdnv(i,2)**2  ! Normal momentum
     etot = qgdnv(i,3)*entho + half*qgdnv(i,1)*qgdnv(i,2)**2
     etot = etot             + half*qgdnv(i,1)*qgdnv(i,4)**2
     etot = etot             + half*qgdnv(i,1)*qgdnv(i,5)**2
     fgdnv(i,3) = qgdnv(i,2)*(etot+qgdnv(i,3))     ! Total energy
  end do
  ! Other advected quantities
  do n = 4, nvar
     do i = 3,ngrid-1
        fgdnv(i,n) = fgdnv(i,1)*qgdnv(i,n)
     end do
  end do
end subroutine acoustic

!===============================================================================
!===============================================================================
subroutine llf(qleft,qright,fgdnv,ngrid,nvar)
  integer::ngrid,nvar
  real(dp),dimension(1:ngrid,1:nvar)::qleft,qright,fgdnv

  ! local arrays
  real(dp),dimension(1:ngrid,1:nvar)::fleft,fright,uleft,uright
  real(dp),dimension(1:ngrid)::cmax

  ! local variables
  integer::i,n
  real(dp)::smallp, entho
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr

  ! constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  ! Maximum wave speed
  do i=3,ngrid-1
     rl=max(qleft (i,1),smallr)
     ul=    qleft (i,2)
     pl=max(qleft (i,3),rl*smallp)
     rr=max(qright(i,1),smallr)
     ur=    qright(i,2)
     pr=max(qright(i,3),rr*smallp)
     cl= sqrt(gamma*pl/rl)
     cr= sqrt(gamma*pr/rr)
     cmax(i)=max(abs(ul)+cl,abs(ur)+cr)
  end do

  ! Compute conservative variables
  do i = 3,ngrid-1
     ! Mass density
     uleft (i,1) = qleft (i,1)
     uright(i,1) = qright(i,1)
     ! Normal momentum
     uleft (i,2) = qleft (i,1)*qleft (i,2)
     uright(i,2) = qright(i,1)*qright(i,2)
     ! Total energy
     uleft (i,3) = qleft (i,3)*entho + half*qleft (i,1)*qleft (i,2)**2
     uright(i,3) = qright(i,3)*entho + half*qright(i,1)*qright(i,2)**2
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,4)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,4)**2
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,5)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,5)**2
  end do
  ! Other advected quantities
  do n = 4, nvar
     do i = 3,ngrid-1
        uleft (i,n) = qleft (i,1)*qleft (i,n)
        uright(i,n) = qright(i,1)*qright(i,n)
     end do
  end do

  ! Compute left and right fluxes
  do i = 3,ngrid-1
     ! Mass density
     fleft (i,1) = uleft (i,2)
     fright(i,1) = uright(i,2)
     ! Normal momentum
     fleft (i,2) = qleft (i,3)+uleft (i,2)*qleft (i,2)
     fright(i,2) = qright(i,3)+uright(i,2)*qright(i,2)
     ! Total energy
     fleft (i,3) = qleft (i,2)*(uleft (i,3)+qleft (i,3))
     fright(i,3) = qright(i,2)*(uright(i,3)+qright(i,3))
  end do
  ! Other advected quantities
  do n = 4, nvar
     do i = 3,ngrid-1
        fleft (i,n) = fleft (i,1)*qleft (i,n)
        fright(i,n) = fright(i,1)*qright(i,n)
     end do
  end do

  ! Compute Lax-Friedrich fluxes
  do n = 1, nvar
     do i = 3,ngrid-1
        fgdnv(i,n) = half*(fleft(i,n)+fright(i,n)-cmax(i)*(uright(i,n)-uleft(i,n)))
     end do
  end do
end subroutine llf

!===============================================================================
!===============================================================================
subroutine hllc (qleft,qright,fgdnv,ngrid,nvar,isothermal,error,detailed_timer)
  integer::ngrid,nvar
  real(dp),dimension(1:ngrid,1:nvar)::qleft,qright,fgdnv
  REAL(dp)::SL,SR
  REAL(dp)::entho
  REAL(dp)::rl,pl,ul,ecinl,etotl,ptotl
  REAL(dp)::rr,pr,ur,ecinr,etotr,ptotr
  REAL(dp)::cfastl,rcl,rstarl
  REAL(dp)::cfastr,rcr,rstarr
  REAL(dp)::etotstarl,etotstarr
  REAL(dp)::ustar,ptotstar
  REAL(dp)::ro,uo,ptoto,etoto
  REAL(dp)::smallp
  INTEGER ::ivar, i
  LOGICAL :: isothermal,error,detailed_timer
  INTEGER, SAVE:: itimer=0
  !
  if (detailed_timer) call trace_begin ('riemann::hllc', itimer=itimer)
  ! constants
  smallp = smallc**2/gamma
  if (isothermal .or. gamma==1d0) then
    entho=one
  else
    entho=one/(gamma-one)
  end if
  error = .false.
  !
  if (io%verbose>1) then
    do i=3,ngrid-1
      !
      ! Left variables
      rl=max(qleft (i,1),smallr)
      Pl=max(qleft (i,3),rl*smallp)
      ul=    qleft (i,2)
      !
      ecinl = half*rl*(ul*ul+qleft(i,4)**2+qleft(i,5)**2)
      etotl = Pl*entho+ecinl
      Ptotl = Pl
      !
      ! Right variables
      rr=max(qright(i,1),smallr)
      Pr=max(qright(i,3),rr*smallp)
      ur=    qright(i,2)
      !
      ecinr = half*rr*(ur*ur+qright(i,4)**2+qright(i,5)**2)
      etotr = Pr*entho+ecinr
      Ptotr = Pr
      !
      ! Find the largest eigenvalues in the normal direction to the interface
      cfastl=sqrt(max(gamma*Pl/rl,smallc**2))
      cfastr=sqrt(max(gamma*Pr/rr,smallc**2))
      !
      ! Compute HLL wave speed
      SL=min(ul,ur)-max(cfastl,cfastr)
      SR=max(ul,ur)+max(cfastl,cfastr)
      !
      ! Compute lagrangian sound speed
      rcl=rl*(ul-SL)
      rcr=rr*(SR-ur)
      if ((rcl+rcl)==0.0) then
        print 1,'HLLC: rl,Pl,ul,ecinl,etotl,Ptotl,cfastl,SL =', &
          i,rl,Pl,ul,ecinl,etotl,Ptotl,cfastl,SL,SL-ul
        print 1,'HLLC: rr,Pr,ur,ecinr,etotr,Ptotr,cfastr,SR =', &
          i,rr,Pr,ur,ecinr,etotr,Ptotr,cfastr,SR,SR-ur
        1 format(a,i4,1p,10e12.3)
        error = .true.
        return
      end if
      !
      ! Compute acoustic star state
      ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
      Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)
      if ((SL-ustar)==0.0 .or. (SR-ustar)==0.0) then
        print 1,'HLLC: rl,Pl,ul,ecinl,etotl,Ptotl,cfastl,SL =', &
          i,rl,Pl,ul,ecinl,etotl,Ptotl,cfastl,SL,SL-ul
        print 1,'HLLC: rr,Pr,ur,ecinr,etotr,Ptotr,cfastr,SR =', &
          i,rr,Pr,ur,ecinr,etotr,Ptotr,cfastr,SR,SR-ur
        print 1, 'HLLC: ustar,Ptotstar', i,ustar,Ptotstar
        call mpi%abort ('HLLC')
        error = .true.
        return
      end if
      !
      ! Left star region variables
      rstarl=rl*(SL-ul)/(SL-ustar)
      etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
      !
      ! Right star region variables
      rstarr=rr*(SR-ur)/(SR-ustar)
      etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
      !
      ! Sample the solution at x/t=0
      if(SL>0d0)then
         ro=rl
         uo=ul
         Ptoto=Ptotl
         etoto=etotl
      else if(ustar>0d0)then
         ro=rstarl
         uo=ustar
         Ptoto=Ptotstar
         etoto=etotstarl
      else if (SR>0d0)then
         ro=rstarr
         uo=ustar
         Ptoto=Ptotstar
         etoto=etotstarr
      else
         ro=rr
         uo=ur
         Ptoto=Ptotr
         etoto=etotr
      end if
      !
      ! Compute the Godunov flux
      fgdnv(i,1) = ro*uo
      fgdnv(i,2) = ro*uo*uo+Ptoto
      fgdnv(i,3) = (etoto+Ptoto)*uo
    end do
  else
    do i=3,ngrid-1
      !
      ! Left variables
      rl=max(qleft (i,1),smallr)
      Pl=max(qleft (i,3),rl*smallp)
      ul=    qleft (i,2)
      !
      ecinl = half*rl*(ul*ul+qleft(i,4)**2+qleft(i,5)**2)
      etotl = Pl*entho+ecinl
      Ptotl = Pl
      !
      ! Right variables
      rr=max(qright(i,1),smallr)
      Pr=max(qright(i,3),rr*smallp)
      ur=    qright(i,2)
      !
      ecinr = half*rr*(ur*ur+qright(i,4)**2+qright(i,5)**2)
      etotr = Pr*entho+ecinr
      Ptotr = Pr
      !
      ! Find the largest eigenvalues in the normal direction to the interface
      cfastl=sqrt(max(gamma*Pl/rl,smallc**2))
      cfastr=sqrt(max(gamma*Pr/rr,smallc**2))
      !
      ! Compute HLL wave speed
      SL=min(ul,ur)-max(cfastl,cfastr)
      SR=max(ul,ur)+max(cfastl,cfastr)
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
      etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
      !
      ! Right star region variables
      rstarr=rr*(SR-ur)/(SR-ustar)
      etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
      !
      ! Sample the solution at x/t=0
!define MERGE
#ifdef MERGE
      ro=rr
      uo=ur
      Ptoto=Ptotr
      etoto=etotr
      ro   =merge(rstarr   ,ro     ,SR>0d0)
      uo   =merge(ustar    ,uo     ,SR>0d0)
      Ptoto=merge(Ptotstar ,Ptoto  ,SR>0d0)
      etoto=merge(etotstarr,etoto  ,SR>0d0)
      ro   =merge(rstarr   ,ro     ,ustar>0d0)
      uo   =merge(ustar    ,uo     ,ustar>0d0)
      Ptoto=merge(Ptotstar ,Ptoto  ,ustar>0d0)
      etoto=merge(etotstarl,etoto  ,ustar>0d0)
      ro   =merge(rl       ,ro     ,SL>0d0)
      uo   =merge(ul       ,uo     ,SL>0d0)
      Ptoto=merge(Ptotl    ,Ptoto  ,SL>0d0)
      etoto=merge(etotl    ,etoto  ,SL>0d0)
#else
      if(SL>0d0)then
         ro=rl
         uo=ul
         Ptoto=Ptotl
         etoto=etotl
      else if(ustar>0d0)then
         ro=rstarl
         uo=ustar
         Ptoto=Ptotstar
         etoto=etotstarl
      else if (SR>0d0)then
         ro=rstarr
         uo=ustar
         Ptoto=Ptotstar
         etoto=etotstarr
      else
         ro=rr
         uo=ur
         Ptoto=Ptotr
         etoto=etotr
      end if
#endif
      !
      ! Compute the Godunov flux
      fgdnv(i,1) = ro*uo
      fgdnv(i,2) = ro*uo*uo+Ptoto
      fgdnv(i,3) = (etoto+Ptoto)*uo
    end do
  end if
  !
  do ivar=4,nvar
    do i=3,ngrid-1
      fgdnv(i,ivar) = fgdnv(i,1)*merge(qleft(i,ivar),qright(i,ivar),fgdnv(i,1)>0)
    end do
  end do
  if (detailed_timer) call trace_end (itimer)
end subroutine hllc

!===============================================================================
!===============================================================================
subroutine hllc_a (qleft,qright,fgdnv,ngrid,nvar)
  integer::ngrid,nvar
  real(dp),dimension(1:ngrid,1:nvar)::qleft,qright,fgdnv
  REAL(dp)::entho
  REAL(dp)::pl,ecinl
  REAL(dp)::pr,ecinr
  REAL(dp)::cfastl,rcl
  REAL(dp)::cfastr,rcr
  REAL(dp)::smallp
  REAL(dp),dimension(ngrid)::SL,SR
  REAL(dp),dimension(ngrid)::ro,uo,Ptoto,Etoto
  REAL(dp),dimension(ngrid)::rl,ul,Ptotl,Etotl
  REAL(dp),dimension(ngrid)::rr,ur,Ptotr,Etotr
  REAL(dp),dimension(ngrid)::rstarl,Etotstarl
  REAL(dp),dimension(ngrid)::rstarr,Etotstarr
  REAL(dp),dimension(ngrid)::ustar,Ptotstar
  INTEGER ::ivar, i
  !
  ! constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)
  !
  do i=3,ngrid-1
     !
     ! Left variables
     rl(i)=max(qleft (i,1),smallr)
     Pl=max(qleft (i,3),rl(i)*smallp)
     ul(i)=    qleft (i,2)
     !
     ecinl = half*rl(i)*(ul(i)*ul(i)+qleft(i,4)**2+qleft(i,5)**2)
     etotl(i) = Pl*entho+ecinl
     Ptotl(i) = Pl
     !
     ! Right variables
     rr(i)=max(qright(i,1),smallr)
     Pr=max(qright(i,3),rr(i)*smallp)
     ur(i)=    qright(i,2)
     !
     ecinr = half*rr(i)*(ur(i)*ur(i)+qright(i,4)**2+qright(i,5)**2)
     etotr(i) = Pr*entho+ecinr
     Ptotr(i) = Pr
     !
     ! Find the largest eigenvalues in the normal direction to the interface
     cfastl=sqrt(max(gamma*Pl/rl(i),smallc**2))
     cfastr=sqrt(max(gamma*Pr/rr(i),smallc**2))
     !
     ! Compute HLL wave speed
     SL(i)=min(ul(i),ur(i))-max(cfastl,cfastr)
     SR(i)=max(ul(i),ur(i))+max(cfastl,cfastr)
     !
     ! Compute lagrangian sound speed
     rcl=rl(i)*(ul(i)-SL(i))
     rcr=rr(i)*(SR(i)-ur(i))
     !
     ! Compute acoustic star state
     ustar(i)   =(rcr*ur(i)   +rcl*ul(i)   +  (Ptotl(i)-Ptotr(i)))/(rcr+rcl)
     Ptotstar(i)=(rcr*Ptotl(i)+rcl*Ptotr(i)+rcl*rcr*(ul(i)-ur(i)))/(rcr+rcl)
  end do
  do i=3,ngrid-1
     !
     ! Left star region variables
     rstarl(i)=rl(i)*(SL(i)-ul(i))/(SL(i)-ustar(i))
     etotstarl(i)=((SL(i)-ul(i))*etotl(i)-Ptotl(i)*ul(i)+Ptotstar(i)*ustar(i))/(SL(i)-ustar(i))
     !
     ! Right star region variables
     rstarr(i)=rr(i)*(SR(i)-ur(i))/(SR(i)-ustar(i))
     etotstarr(i)=((SR(i)-ur(i))*etotr(i)-Ptotr(i)*ur(i)+Ptotstar(i)*ustar(i))/(SR(i)-ustar(i))
  end do
     !
     ! Sample the solution at x/t=0
  do i=3,ngrid-1
     if(SL(i)>0d0)then
        ro(i)=rl(i)
        uo(i)=ul(i)
        Ptoto(i)=Ptotl(i)
        etoto(i)=etotl(i)
     else if(ustar(i)>0d0)then
        ro(i)=rstarl(i)
        uo(i)=ustar(i)
        Ptoto(i)=Ptotstar(i)
        etoto(i)=etotstarl(i)
     else if (SR(i)>0d0)then
        ro(i)=rstarr(i)
        uo(i)=ustar(i)
        Ptoto(i)=Ptotstar(i)
        etoto(i)=etotstarr(i)
     else
        ro(i)=rr(i)
        uo(i)=ur(i)
        Ptoto(i)=Ptotr(i)
        etoto(i)=etotr(i)
     end if
  end do
     !
     ! Compute the Godunov flux
  do i=3,ngrid-1
     fgdnv(i,1) = ro(i)*uo(i)
     fgdnv(i,2) = ro(i)*uo(i)*uo(i)+Ptoto(i)
     fgdnv(i,3) = (etoto(i)+Ptoto(i))*uo(i)
  end do
     !
  do ivar=4,nvar
  do i=3,ngrid-1
     fgdnv(i,ivar) = merge(fgdnv(i,1)*qleft (i,ivar), &
                           fgdnv(i,1)*qright(i,ivar), fgdnv(i,1)>0)
  end do
  end do
end subroutine hllc_a

!===============================================================================
!===============================================================================
subroutine hllc_v (qleft,qright,fgdnv,ngrid,nvar)
  integer::ngrid,nvar
  real(dp),dimension(1:ngrid,1:nvar)::qleft,qright,fgdnv
  REAL(dp),dimension(ngrid)::SL,SR
  REAL(dp),dimension(ngrid)::entho
  REAL(dp),dimension(ngrid)::rl,pl,ul,ecinl,etotl,ptotl
  REAL(dp),dimension(ngrid)::rr,pr,ur,ecinr,etotr,ptotr
  REAL(dp),dimension(ngrid)::cfastl,rcl,rstarl
  REAL(dp),dimension(ngrid)::cfastr,rcr,rstarr
  REAL(dp),dimension(ngrid)::etotstarl,etotstarr
  REAL(dp),dimension(ngrid)::ustar,ptotstar
  REAL(dp),dimension(ngrid)::ro,uo,ptoto,etoto
  REAL(dp),dimension(ngrid)::smallp
  INTEGER ::ivar, i
  !
  ! constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)
  !
     !
     ! Left variables
     rl=max(qleft (:,1),smallr)
     Pl=max(qleft (:,3),rl*smallp)
     ul=    qleft (:,2)
     !
     ecinl = half*rl*(ul*ul+qleft(:,4)**2+qleft(:,5)**2)
     etotl = Pl*entho+ecinl
     Ptotl = Pl
     !
     ! Right variables
     rr=max(qright(:,1),smallr)
     Pr=max(qright(:,3),rr*smallp)
     ur=    qright(:,2)
     !
     ecinr = half*rr*(ur*ur+qright(:,4)**2+qright(:,5)**2)
     etotr = Pr*entho+ecinr
     Ptotr = Pr
     !
     ! Find the largest eigenvalues in the normal direction to the interface
     cfastl=sqrt(max(gamma*Pl/rl,smallc**2))
     cfastr=sqrt(max(gamma*Pr/rr,smallc**2))
     !
     ! Compute HLL wave speed
     SL=min(ul,ur)-max(cfastl,cfastr)
     SR=max(ul,ur)+max(cfastl,cfastr)
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
     etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
     !
     ! Right star region variables
     rstarr=rr*(SR-ur)/(SR-ustar)
     etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
     !
     ! Sample the solution at x/t=0
     where (SL>0d0)
        ro=rl
        uo=ul
        Ptoto=Ptotl
        etoto=etotl
     else where (ustar>0d0)
        ro=rstarl
        uo=ustar
        Ptoto=Ptotstar
        etoto=etotstarl
     else where  (SR>0d0)
        ro=rstarr
        uo=ustar
        Ptoto=Ptotstar
        etoto=etotstarr
     else where
        ro=rr
        uo=ur
        Ptoto=Ptotr
        etoto=etotr
     end where
     !
     ! Compute the Godunov flux
     fgdnv(:,1) = ro*uo
     fgdnv(:,2) = ro*uo*uo+Ptoto
     fgdnv(:,3) = (etoto+Ptoto)*uo
     do ivar = 4,nvar
        where (fgdnv(:,1)>0)
           fgdnv(:,ivar) = fgdnv(:,1)*qleft (:,ivar)
        else where
           fgdnv(:,ivar) = fgdnv(:,1)*qright(:,ivar)
        end where
     end do
end subroutine hllc_v

!===============================================================================
!===============================================================================
subroutine hll(qleft,qright,fgdnv,ngrid,nvar)
  integer::ngrid,nvar
  real(dp),dimension(1:ngrid,1:nvar)::qleft,qright,fgdnv
  real(dp),dimension(1:ngrid,1:nvar)::fleft,fright,uleft,uright
  real(dp),dimension(1:ngrid)::SL,SR
  integer::i,n
  real(dp)::smallp, entho
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr

  ! constants
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  ! Maximum wave speed
  do i=3,ngrid-1
     rl=max(qleft (i,1),smallr)
     ul=    qleft (i,2)
     pl=max(qleft (i,3),rl*smallp)
     rr=max(qright(i,1),smallr)
     ur=    qright(i,2)
     pr=max(qright(i,3),rr*smallp)
     cl= sqrt(gamma*pl/rl)
     cr= sqrt(gamma*pr/rr)
     SL(i)=min(min(ul,ur)-max(cl,cr),zero)
     SR(i)=max(max(ul,ur)+max(cl,cr),zero)
  end do

  ! Compute conservative variables
  do i = 3,ngrid-1
     uleft (i,1) = qleft (i,1)
     uright(i,1) = qright(i,1)
     uleft (i,2) = qleft (i,1)*qleft (i,2)
     uright(i,2) = qright(i,1)*qright(i,2)
     uleft (i,3) = qleft (i,3)*entho + half*qleft (i,1)*qleft (i,2)**2
     uright(i,3) = qright(i,3)*entho + half*qright(i,1)*qright(i,2)**2
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,4)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,4)**2
     uleft (i,3) = uleft (i,3)       + half*qleft (i,1)*qleft (i,5)**2
     uright(i,3) = uright(i,3)       + half*qright(i,1)*qright(i,5)**2
  end do
  ! Other advected quantities
  do n = 4, nvar
     do i = 3,ngrid-1
        uleft (i,n) = qleft (i,1)*qleft (i,n)
        uright(i,n) = qright(i,1)*qright(i,n)
     end do
  end do

  ! Compute left and right fluxes
  do i = 3,ngrid-1
     fleft (i,1) = uleft (i,2)
     fright(i,1) = uright(i,2)
     fleft (i,2) = qleft (i,3)+uleft (i,2)*qleft (i,2)
     fright(i,2) = qright(i,3)+uright(i,2)*qright(i,2)
     fleft (i,3) = qleft (i,2)*(uleft (i,3)+qleft (i,3))
     fright(i,3) = qright(i,2)*(uright(i,3)+qright(i,3))
  end do
  ! Other advected quantities
  do n = 4, nvar
     do i = 3,ngrid-1
        fleft (i,n) = fleft (i,1)*qleft (i,n)
        fright(i,n) = fright(i,1)*qright(i,n)
     end do
  end do

  ! Compute HLL fluxes
  do n = 1, nvar
     do i = 3,ngrid-1
        fgdnv(i,n) = (SR(i)*fleft(i,n)-SL(i)*fright(i,n) &
             & + SR(i)*SL(i)*(uright(i,n)-uleft(i,n)))/(SR(i)-SL(i))
     end do
  end do
end subroutine hll

end module riemann_mod
