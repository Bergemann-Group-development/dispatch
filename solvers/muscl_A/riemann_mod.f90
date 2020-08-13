!define DBG
!===============================================================================
!> Riemann solvers for HD
!===============================================================================
MODULE riemann_mod
  USE io_mod
  USE timer_mod
  implicit none
  private
  integer, parameter:: sp=4, dp=8
  real, parameter:: one=1.0, half=0.5, zero=0.0
  type, public:: riemann_t
    real(8)          :: gamma,g1, g2
    character(len=4) :: method='hll'
  contains
    procedure        :: solve
  end type
  type(riemann_t), public:: riemann
CONTAINS

!===============================================================================
!> Interface to choose method
!===============================================================================
SUBROUTINE solve (self, dir, l, u, left, rght, flux)
  class(riemann_t):: self
  integer:: dir, l(3), u(3)
  real, dimension(:,:,:,:), contiguous, intent(in) :: left
  real, dimension(:,:,:,:),             intent(in) :: rght
  real, dimension(:,:,:,:), contiguous, intent(out):: flux
  !.............................................................................
  select case (trim(self%method))
  case ('hll')
    call hll  (self, dir, l, u, left, rght, flux)
  case ('hllc')
    call hllc (self, dir, l, u, left, rght, flux)
  case ('hlld')
    call hlld (self, dir, l, u, left, rght, flux)
  end select
END SUBROUTINE solve

!===============================================================================
!> Solve for the HLL flux give two HD state vectors (left and right of the
!> interface) stored in s
!===============================================================================
SUBROUTINE hll (self, dir, l, u, left, rght, flux)
  class(riemann_t):: self
  integer:: dir, l(3), u(3)
  real, dimension(:,:,:,:), contiguous, intent(in) :: left
  real, dimension(:,:,:,:),             intent(in) :: rght
  real, dimension(:,:,:,:), contiguous, intent(out):: flux
  !.............................................................................
  real(sp):: &
    Ll, Dl, Ul, Pl, Sr, Vl, Wl, Tl, cl, fDl, fSl, fUl, fVl, fWl, SSl, &
    Lr, Dr, Ur, Pr, Sl, Vr, Wr, Tr, cr, fDr, fSr, fUr, fVr, fWr, SSr, &
    cmax
  integer:: i1, i2, i3, j1, j2, j3, v1, v2, v3
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  !  Parallel and perpendicular velocity components
  !-----------------------------------------------------------------------------
  v1 = 2 + mod(dir-1,3)
  v2 = 2 + mod(dir  ,3)
  v3 = 2 + mod(dir+1,3)
  !-----------------------------------------------------------------------------
  ! For the left states, we need to refer to index-1 in the flux direction
  !-----------------------------------------------------------------------------
  if (dir==1) then
    j1=1; j2=0; j3=0
  else if (dir==2) then
    j1=0; j2=1; j3=0
  else
    j1=0; j2=0; j3=1
  end if
  !-----------------------------------------------------------------------------
  do i3=l(3),u(3)
  do i2=l(2),u(2)
  do i1=l(1),u(1)
    ll = left(i1-j1,i2-j2,i3-j3,1)
    sl = left(i1-j1,i2-j2,i3-j3,5)
    ul = left(i1-j1,i2-j2,i3-j3,v1)
    vl = left(i1-j1,i2-j2,i3-j3,v2)
    wl = left(i1-j1,i2-j2,i3-j3,v3)
    lr = left(i1-j1,i2-j2,i3-j3,1)
    sr = left(i1-j1,i2-j2,i3-j3,5)
    ur = left(i1-j1,i2-j2,i3-j3,v1)
    vr = left(i1-j1,i2-j2,i3-j3,v2)
    wr = left(i1-j1,i2-j2,i3-j3,v3)
    dl = exp(ll)
    dr = exp(lr)
    pl = exp(self%g1*sl + self%gamma*ll)
    pr = exp(self%g1*sr + self%gamma*lr)
    !---------------------------------------------------------------------------
    ! Pressure and sound speed for each side of interface (l==left, r==right)
    !---------------------------------------------------------------------------
    cl = sqrt(self%gamma*pl/dl)
    cr = sqrt(self%gamma*pr/dr)
    cmax = max(cl,cr)
    !---------------------------------------------------------------------------
    ! Hydro fluxes
    !---------------------------------------------------------------------------
    fDl = Dl*Ul
    fDr = Dr*Ur
    fUl = Dl*Ul*Ul + Pl
    fUr = Dr*Ur*Ur + Pr
    fSl = Dl*Ul*sl
    fSr = Dl*Ur*sr
    fVl = Dl*Ul*Vl
    fVr = Dr*Ur*Vr
    fWl = Dl*Ul*Wl
    fWr = Dr*Ur*Wr
    !---------------------------------------------------------------------------
    ! maximum wave speeds to the l and r (guaranteed to have correct sign)
    !---------------------------------------------------------------------------
    SSL = min(min(Ul,Ur)-cmax,0.0)
    SSR = max(max(Ul,Ur)+cmax,0.0)
    !---------------------------------------------------------------------------
    ! HLL flux based on wavespeeds.  The general form is 
    !    (SR*F_l - SL*F_r + SL*SR *(u_r - u_l))/(SR-SL)
    ! where u = (rho, rho*v, E_tot) are the conserved variables
    !---------------------------------------------------------------------------
    flux(i1,i2,i3, 1) = (SSR*fDl - SSL*fDr + SSL*SSR*(Dr    - Dl   )) / (SSR-SSL)
    flux(i1,i2,i3,v1) = (SSR*fUl - SSL*fUr + SSL*SSR*(Dr*Ur - Dl*Ul)) / (SSR-SSL)
    flux(i1,i2,i3,v2) = (SSR*fVl - SSL*fVr + SSL*SSR*(Dr*Vr - Dl*Vl)) / (SSR-SSL)
    flux(i1,i2,i3,v3) = (SSR*fWl - SSL*fWr + SSL*SSR*(Dr*Wr - Dl*Wl)) / (SSR-SSL)
    flux(i1,i2,i3, 5) = (SSR*fSl - SSL*fSr + SSL*SSR*(Sr    - Sl   )) / (SSR-SSL)
#ifdef DBG
    if (dir==1 .and. i2==l(2) .and. i3==l(3)) then
      print '(i3,3(2x,5f10.4))', &
        i1, Dl,ul,sl,Pl,SSL, Dr,ur,sr,Pr,SSR, flux(i1,i2,i3,1:5)
      end if
#endif
  end do
  end do
  end do
END SUBROUTINE hll

!===============================================================================
!> Solve for the HLLC flux give two HD state vectors (left and right of the
!> interface) stored in fgdnv
!===============================================================================
SUBROUTINE hllc (self, dir, l, u, qleft, qrght, fgdnv)
  implicit none
  class(riemann_t):: self
  !-----------------------------------------------------------------------------
  ! Parameters
  !-----------------------------------------------------------------------------
  real, parameter:: half=0.5, zero=0.0, one=1.0
  !--- single precision
  REAL(sp),dimension(:,:,:,:):: qleft,qrght,fgdnv
  REAL(dp):: ll,rl,sl,pl,ul,vl,wl,ecinl,stotl,ptotl
  REAL(dp):: lr,rr,sr,pr,ur,vr,wr,ecinr,stotr,ptotr
  REAL(dp):: rstarl,stotstarl
  REAL(dp):: rstarr,stotstarr
  REAL(dp):: ro,uo,so,ptoto,stoto
  REAL(dp):: ptotstar
  !--- double precision
  REAL(dp):: cfastl,SSL,rcl
  REAL(dp):: cfastr,SSR,rcr
  REAL(dp):: ustar
  !---
  INTEGER:: l(3),u(3),i2,i3,dir,n_cheap,n_expensive
  INTEGER:: i,j1,j2,j3,v1,v2,v3
  INTEGER, SAVE:: itimer=0
  !-----------------------------------------------------------------------------
  ! Choose variable indices; v1=parallel, v2-v3: perpendicular
  !-----------------------------------------------------------------------------
  v1 = 2 + mod(dir-1,3)
  v2 = 2 + mod(dir  ,3)
  v3 = 2 + mod(dir+1,3)
  n_cheap = 0
  n_expensive = 0
  !-----------------------------------------------------------------------------
  ! For the left states, we need to refer to index-1 in the flux direction
  !-----------------------------------------------------------------------------
  if (dir==1) then
    j1=1; j2=0; j3=0
  else if (dir==2) then
    j1=0; j2=1; j3=0
  else
    j1=0; j2=0; j3=1
  end if
  !---------------------------------------------------------------------------
  ! Loops over 2nd and 3rd index
  !---------------------------------------------------------------------------
  do i3=l(3),u(3)
  do i2=l(2),u(2)
    !---------------------------------------------------------------------------
    ! Loop over 1st index
    !---------------------------------------------------------------------------
    do i=l(1),u(1)
      !
      ! Left variables
      ll = qleft(i-j1,i2-j2,i3-j3,1)
      sl = qleft(i-j1,i2-j2,i3-j3,5)
      ul = qleft(i-j1,i2-j2,i3-j3,v1)
      vl = qleft(i-j1,i2-j2,i3-j3,v2)
      wl = qleft(i-j1,i2-j2,i3-j3,v3)
      rl = exp(ll)
      sl = sl*rl
      pl = exp(self%g1*sl + self%gamma*ll)
      !
      ! Right variables
      lr = qrght(i,i2,i3,1)
      sr = qrght(i,i2,i3,5)
      ur = qrght(i,i2,i3,v1)
      vr = qrght(i,i2,i3,v2)
      wr = qrght(i,i2,i3,v3)
      rr = exp(lr)
      sr = sr*rr
      pr = exp(self%g1*sr + self%gamma*lr)
      !
      ! Find the largest eigenvalues in the normal direction to the interface
      cfastl = sqrt(self%gamma*pl/rl)
      cfastr = sqrt(self%gamma*pr/rr)
      !
      ! Compute HLL wave speed; SL positive when positive supersonic, SR negative when negative ss
      SSL = min(ul,ur) - max(cfastl,cfastr)
      SSR = max(ul,ur) + max(cfastl,cfastr)
      if (SSL > 0.0) then
        !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'SL case',i
        ro = rl
        uo = ul
        Ptoto = pl
        stoto = sl
        n_cheap = n_cheap+1
      else if (SSR < 0.0) then
        !if (dir==1 .and. i2==l(2) .and. i3==l(3)) print *,'SR case',i
        ro = rr
        uo = ur
        Ptoto = pr
        stoto = sr
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
      fgdnv(i,i2,i3, 1) = ro*uo
      fgdnv(i,i2,i3, 5) = so*uo
      fgdnv(i,i2,i3,v1) = ro*uo*uo + Ptoto
      fgdnv(i,i2,i3,v2) = fgdnv(i,i2,i3,1)*merge(vl,vr,uo>0.0)
      fgdnv(i,i2,i3,v3) = fgdnv(i,i2,i3,1)*merge(wl,wr,uo>0.0)
!define DBGC
#ifdef DBGC
      if (dir==1 .and. i2==l(2) .and. i3==l(3)) then
        print '(i3,4(2x,5f10.4))', &
          i, rl,ul,sl,Pl,SSL, rr,ur,sr,Pr,SSR, &
             ro,uo,so,Ptoto,0.0, fgdnv(i,i2,i3,1:5)
      end if
#endif
    end do
  end do
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
!> HLLD method based on log density, velocity, entropy, and magnetic vector
!> potential
!===============================================================================
!subroutine hlld (self, qleft, qrght, fgdnv, ngrid, detailed_timer)
SUBROUTINE hlld (self, dir, li, ui, qleft, qrght, fgdnv)
  IMPLICIT none
  class(riemann_t):: self
  integer:: dir,li(3),ui(3)
  real(sp),dimension(:,:,:,:)::qleft,qrght,fgdnv
  !.............................................................................
  REAL(dp)::SSL,SSR
  REAL(dp)::re,le,lim
  REAL(dp)::ecinl,stotl,ptotl,wstarl,vtotBl,emagl,SAL,sqrrstarl
  REAL(dp)::ecinr,stotr,ptotr,wstarr,vtotBr,emagr,SAR,sqrrstarr
  REAL(dp)::cfastl,rcl,rstarl,sl,vdotBstarl,vstarl,stotstarstarl,stotstarl,vdotBl
  REAL(dp)::cfastr,rcr,rstarr,sr,vdotBstarr,vstarr,stotstarstarr,stotstarr,vdotBr
  REAL(dp)::ustar,ptotstar,sstar,vstarstar,wstarstar
  REAL(dp)::ro,uo,ptoto,stoto,vo,wo
  REAL(sp)::d,P,u,v,w,A,B,C,sgnm
  REAL(sp)::ll,rl,el,pl,ul,vl,wl,Al,Bl,Cl
  REAL(sp)::lr,rr,er,pr,ur,vr,wr,Ar,Br,Cr
  REAL(dp)::Bstarl,Bstarr,Cstarl,Cstarr,calfvenl,calfvenr
  REAL(dp)::Bo,Co,Bstarstar,Cstarstar,vdotBstarstar,vdotBo,emago,estar,so
  REAL(dp)::B2,c2,d2,cf
  INTEGER ::ivar,i,i2,i3,j1,j2,j3,v1,v2,v3,icase
  logical:: flag
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  ! Choose variable indices; v1=parallel, v2-v3: perpendicular
  !-----------------------------------------------------------------------------
  v1 = 2 + mod(dir-1,3)
  v2 = 2 + mod(dir  ,3)
  v3 = 2 + mod(dir+1,3)
  !-----------------------------------------------------------------------------
  ! For the left states, we need to refer to index-1 in the flux direction
  !-----------------------------------------------------------------------------
  if (dir==1) then
    j1=1; j2=0; j3=0
  else if (dir==2) then
    j1=0; j2=1; j3=0
  else
    j1=0; j2=0; j3=1
  end if
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  do i3=li(3),ui(3)
  do i2=li(2),ui(2)
!define PRIVATE
#ifdef PRIVATE
  !$omp simd private(sgnm,d,P,u,v,w,A,B,C,B2,c2,d2,cf,Ptotstar,sstar, &
  !$omp &  sl,rl,vl,ecinl,emagl,stotl,Ptotl,vdotBl,Bl,Cl,cfastl,SSL, &
  !$omp &  sr,rr,vr,ecinr,emagr,stotr,Ptotr,vdotBr,Br,Cr,cfastr,SSR, &
  !$omp &  rcl,rstarl,sl,Bstarl,Cstarl,vstarl,wstarl, &
  !$omp &  rcr,rstarr,sr,Bstarr,Cstarr,vstarr,wstarr, &
  !$omp &  vdotBstarl,stotstarl,sqrrstarl,calfvenl,SAR, &
  !$omp &  vdotBstarr,stotstarr,sqrrstarr,calfvenr,SAL, &
  !$omp &  vstarstar,wstarstar,Bstarstar,Cstarstar,vdotBstarstar, &
  !$omp &  icase,ro,so,uo,vo,wo,Bo,Co,Ptoto,etoto,vdotBo,emago,ustar)
#endif
  do i=li(1),ui(1)
    !print '(i3,1p,8e10.2,2x,8e10.2)',i,qleft(i,:),qrght(i,:)
    !
    ! Enforce continuity of normal component
    A = half*(qleft(i-j1,i2-j2,i3-j3,v1+4)+&
              qrght(i   ,i2   ,i3   ,v1+4))
    sgnm=sign(one,A)
    qleft(i-j1,i2-j2,i3-j3,v1+4) = A
    qrght(i   ,i2   ,i3   ,v1+4) = A
    !
    ! Left variables
    ll = qleft(i-j1,i2-j2,i3-j3,1)
    ul = qleft(i-j1,i2-j2,i3-j3,v1)
    vl = qleft(i-j1,i2-j2,i3-j3,v2)
    wl = qleft(i-j1,i2-j2,i3-j3,v3)
    sl = qleft(i-j1,i2-j2,i3-j3,5)
    Bl = qleft(i-j1,i2-j2,i3-j3,v2+4)
    Cl = qleft(i-j1,i2-j2,i3-j3,v3+4)
    !
    rl = exp(ll)
    pl = exp(self%g1*sl + self%gamma*ll)
    ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
    emagl = half*(A*A+Bl*Bl+Cl*Cl)
    Ptotl = Pl + emagl
    vdotBl= ul*A+vl*Bl+wl*cl
    !
    ! Right variables
    lr = qrght(i,i2,i3,1)
    ur = qrght(i,i2,i3,v1)
    vr = qrght(i,i2,i3,v2)
    wr = qrght(i,i2,i3,v3)
    sr = qrght(i,i2,i3,5)
    Br = qrght(i,i2,i3,v2+4)
    Cr = qrght(i,i2,i3,v3+4)
    rr = exp(lr)
    pr = exp(self%g1*sr + self%gamma*lr)
    er = pr/self%g1
    ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
    emagr = half*(A*A+Br*Br+Cr*Cr)
    Ptotr = Pr + emagr
    vdotBr= ur*A+vr*Br+wr*Cr
    !
    ! Find the largest eigenvalues in the normal direction to the interface
    B2 = A*A+Br*Br+Cr*Cr
    c2 = self%gamma*Pr/rr
    d2 = half*(B2/rr+c2)
    cfastl = sqrt(d2 + sqrt(max(d2*d2-c2*A*A/d,zero)))
    B2 = A*A+Bl*Bl+Cl*Cl
    c2 = self%gamma*Pl/rl
    d2 = half*(B2/rl+c2)
    cfastr = sqrt(d2 + sqrt(max(d2*d2-c2*A*A/d,zero)))
    !
    ! Compute HLL wave speeds
    SSL=min(ul,ur)-max(cfastl,cfastr)
    SSR=max(ul,ur)+max(cfastl,cfastr)
    !
    lim = 1.0_dp
    !
    ! Compute lagrangian wave speed
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
!define PROTECT
#ifdef PROTECT
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
#else
    vstarl=vl-A*Bl*(ustar-ul)/estar
    Bstarl=Bl*el/estar
    wstarl=wl-A*Cl*(ustar-ul)/estar
    Cstarl=Cl*el/estar
#endif
    vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
    stotstarl=((SL-ul)*stotl)/(SL-ustar)
    sqrrstarl=sqrt(rstarl)
    calfvenl=abs(A)/sqrrstarl
    SAL=ustar-calfvenl
    !
    ! Right star region variables
    rstarr=rr*(SR-ur)/(SR-ustar)
    estar =rr*(SR-ur)*(SR-ustar)-A**2
    er    =rr*(SR-ur)*(SR-ur   )-A**2
#ifdef PROTECT
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
#else
    vstarr=vr-A*Br*(ustar-ur)/estar
    Bstarr=Br*er/estar
    wstarr=wr-A*Cr*(ustar-ur)/estar
    Cstarr=Cr*er/estar
#endif
    !
    vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
    stotstarr=((SR-ur)*stotr)/(SR-ustar)
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
    stotstarstarl=stotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
    stotstarstarr=stotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)
    !
    ! Sample the solution at x/t=0
    if(SSL>0d0)then
       icase=1
       ro=rl
       uo=ul
       vo=vl
       wo=wl
       Bo=Bl
       Co=Cl
       so=sl
       Ptoto=Ptotl
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
       stoto=stotstarl
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
       so=stotstarstarl
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
       so=stotstarstarr
       vdotBo=vdotBstarstar
       emago=half*(A*A+Bo*Bo+Co*Co)
    else if (SSR>0d0)then
       icase=5
       ro=rstarr
       uo=ustar
       vo=vstarr
       wo=wstarr
       Bo=Bstarr
       Co=Cstarr
       Ptoto=Ptotstar
       so=stotstarr
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
       so=sr
       Ptoto=Ptotr
       stoto=stotr
       vdotBo=vdotBr
       emago=emagr
    end if
!define DIFFUSE
#ifdef DIFFUSE
    ! Add a diffusive flux in the down grad(ln(rho)) direction, which is only nonzero when the
    ! smallest of the left/right density is within a factor 10 from smallr, and which has at most
    ! magnitude of 0.25 times the largest of the left/right fast mode speeds.
    ustar = alog(rl/rr)*max(0.0, alog10((10.*smallr)/min(rl,rr)))
    ustar = 0.25*max(cfastr,cfastl)*min(1.0, max(-1.0, ustar))
    uo = uo + ustar
    !
    ndiff = ndiff + merge(1,0,abs(ustar)>0.0)
#endif
    !
    ! Compute the Godunov flux
    fgdnv(i,i2,i3,1 ) = ro*uo
    fgdnv(i,i2,i3,v1) = ro*uo*uo+Ptoto+emago*(lim-1.0_dp)-A*A*lim
    fgdnv(i,i2,i3,v2) = ro*uo*vo-A*Bo*lim
    fgdnv(i,i2,i3,v3) = ro*uo*wo-A*Co*lim
    fgdnv(i,i2,i3,5 ) = so*uo
!print'(1p,i3,12e11.2)',i,fgdnv(i,i2,i3,1:5)
  end do
  end do
  end do
END SUBROUTINE hlld

END MODULE riemann_mod
