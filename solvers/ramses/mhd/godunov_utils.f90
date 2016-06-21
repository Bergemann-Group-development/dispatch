!###########################################################
!###########################################################
!###########################################################
!###########################################################
MODULE find_for_solvers
CONTAINS
SUBROUTINE find_mhd_flux_vec(qvar,cvar,ff,vel)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! compute the 1D MHD fluxes from the conservative variables
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1, Btransverse1, Vtransverse2, Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i , ivar, ib
  REAL(dp),DIMENSION(1:nvector,1:nvar), target, intent(in) :: qvar
  REAL(dp),DIMENSION(1:nvector),                intent(in) :: vel
  REAL(dp),DIMENSION(1:nvector,1:nvar)                     :: cvar, ff
  REAL(dp),DIMENSION(:), pointer:: d,u,v,w,A,B,C,P
  REAL(dp),DIMENSION(1:nvector):: ecin,emag,etot,Ptot,du,lim
  REAL(dp) :: entho
  integer  :: l

  call trace_mpi('find_mhd_flux_vec(qvar,cvar,ff,vel)','enter',3)

  ! Local variables
  entho = one/(gamma-one)
  d=>qvar(:,1); P=>qvar(:,2); u=>qvar(:,3); A=>qvar(:,4)
  v=>qvar(:,5); B=>qvar(:,6); w=>qvar(:,7); C=>qvar(:,8)
  ecin = half*(u*u+v*v+w*w)*d
  emag = half*(A*A+B*B+C*C)
#ifdef KROME
  etot = P / (merge(qvar(:,ichem),gamma,chemistry) - one) + ecin + emag
#else
  etot = P*entho+ecin+emag
#endif
  Ptot = P + emag
  du   = d*u
  
  ! Compute conservative variables
  cvar(:,1) = d
  cvar(:,2) = etot
  cvar(:,3) = du
  cvar(:,4) = A
  cvar(:,5) = d*v
  cvar(:,6) = B
  cvar(:,7) = d*w
  cvar(:,8) = C
  do ivar = 9,nvar
#ifdef KROME
     if (ivar .ne. ichem) then
#endif
        cvar(:,ivar) = d*qvar(:,ivar)
#ifdef KROME
     else  ! Make sure effective adiabatic index is just copied over
        cvar(:,ivar) = qvar(:,ivar)
     endif
#endif
  end do


  ! Current limiter on the velocity flux
  if (cmax2 > 0) then
    lim = cmax2 / (cmax2 + vel*vel)
  else
    lim = one
  endif

  ! Compute fluxes
  ff(:,1) = du
  ff(:,2) = (etot+Ptot)*u-A*(A*u+B*v+C*w)
  ff(:,3) = du*u+P+(emag-A*A)*lim
  ff(:,4) = zero
  ff(:,5) = du*v-A*B*lim
  ff(:,6) = B*u-A*v
  ff(:,7) = du*w-A*C*lim
  ff(:,8) = C*u-A*w
  do ivar = 9,nvar
#ifdef KROME
     if (ivar .ne. ichem) then
#endif
        ff(:,ivar) = du*qvar(:,ivar)
#ifdef KROME
     else  ! Make sure effective adiabatic index has zero flux
        ff(:,ivar) = zero
     endif
#endif
  end do

  call trace_mpi('find_mhd_flux_vec(qvar,cvar,ff,vel)','exit',3)
END SUBROUTINE find_mhd_flux_vec
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_info_vec(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the fastest velocity at which information is exchanged 
  !! at the interface
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  REAL(dp),DIMENSION(1:nvector,1:nvar), target, intent(in) :: qvar  
  REAL(dp),DIMENSION(1:nvector)       :: vel_info
  REAL(dp),DIMENSION(:), pointer:: d,u,v,w,A,B,C,P
  REAL(dp) :: B2,c2,d2,cf,id
  integer::l

  call trace_mpi('find_speed_info_vec(qvar,vel_info)','enter',3)
  
  d=>qvar(:,1); P=>qvar(:,2); u=>qvar(:,3); A=>qvar(:,4)
  v=>qvar(:,5); B=>qvar(:,6); w=>qvar(:,7); C=>qvar(:,8)
  do l=1,nvector
     B2 = A(l)*A(l)+B(l)*B(l)+C(l)*C(l)
     id = one / max(d(l),smallr)
#ifdef KROME
     c2 = merge(qvar(l,ichem),gamma,chemistry)*P(l)*id
#else
     c2 = gamma*P(l)*id
#endif
     d2 = half*(B2*id+c2)
     cf = sqrt( d2 + sqrt(max(d2*d2-c2*A(l)*A(l)*id,zero)) )
     vel_info(l) = cf+abs(u(l))
  enddo
  
  call trace_mpi('find_speed_info_vec(qvar,vel_info)','exit',3)
END SUBROUTINE find_speed_info_vec
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_fast_vec(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the fast magnetosonic velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvector,1:nvar),target:: qvar  
  REAL(dp),DIMENSION(1:nvector)  :: vel_info
  REAL(dp),DIMENSION(:), pointer :: d,P,u,v,w,A,B,C
  REAL(dp),DIMENSION(1:nvector)  :: B2,c2,d2,cf

  call trace_mpi('find_speed_fast_vec(qvar,vel_info)','enter',3)

  d=>qvar(:,1); P=>qvar(:,2); u=>qvar(:,3); A=>qvar(:,4)
  v=>qvar(:,5); B=>qvar(:,6); w=>qvar(:,7); C=>qvar(:,8)
  B2 = A*A+B*B+C*C
#ifdef KROME
  c2 = merge(qvar(:,ichem),gamma,chemistry)*P/d
#else
  c2 = gamma*P/d
#endif
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
  vel_info = cf

  call trace_mpi('find_speed_fast_vec(qvar,vel_info)','exit',3)
END SUBROUTINE find_speed_fast_vec
END MODULE find_for_solvers
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_fast(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the fast magnetosonic velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvar):: qvar  
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf

  call trace_mpi('find_speed_fast(qvar,vel_info)','enter',4)

  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
#ifdef KROME
  c2 = merge(qvar(ichem),gamma,chemistry)*P/d
#else
  c2 = gamma*P/d
#endif
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
  vel_info = cf

  call trace_mpi('find_speed_fast(qvar,vel_info)','exit',4)
END SUBROUTINE find_speed_fast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_sound(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the sound velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvar):: qvar  
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf

  call trace_mpi('find_speed_sound(qvar,vel_info)','enter',4)

  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
#ifdef KROME
  c2 = merge(qvar(ichem),gamma,chemistry)*P/d
#else
  c2 = gamma*P/d
#endif
  vel_info = sqrt(c2)

  call trace_mpi('find_speed_sound(qvar,vel_info)','exit',4)
END SUBROUTINE find_speed_sound
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_mhd_flux(qvar,cvar,ff) 
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! compute the 1D MHD fluxes from the conservative variables
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1, Btransverse1, Vtransverse2, Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i , ivar, ib
  REAL(dp),DIMENSION(1:nvar):: qvar, cvar, ff
  REAL(dp) :: ecin,emag,etot,d,u,v,w,A,B,C,P,Ptot,entho

  call trace_mpi('find_mhd_flux(qvar,cvar,ff)','enter',2)

  ! Local variables
  entho = one/(gamma-one)
  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  ecin = half*(u*u+v*v+w*w)*d
  emag = half*(A*A+B*B+C*C)
#ifdef KROME
  etot = P / (merge(qvar(ichem),gamma,chemistry) - one) + ecin + emag
#else
  etot = P*entho+ecin+emag
#endif
  Ptot = P + emag
  
  ! Compute conservative variables
  cvar(1) = d
  cvar(2) = etot
  cvar(3) = d*u
  cvar(4) = A
  cvar(5) = d*v
  cvar(6) = B
  cvar(7) = d*w
  cvar(8) = C
  do ivar = 9,nvar
#ifdef KROME
     if (ivar .ne. ichem) then
#endif
        cvar(ivar) = d*qvar(ivar)
#ifdef KROME
     else  ! Make sure effective adiabatic index is just copied over
        cvar(ivar) = qvar(ivar)
     endif
#endif
  end do

  ! Compute fluxes
  ff(1) = d*u
  ff(2) = (etot+Ptot)*u-A*(A*u+B*v+C*w)
  ff(3) = d*u*u+Ptot-A*A
  ff(4) = zero
  ff(5) = d*u*v-A*B
  ff(6) = B*u-A*v
  ff(7) = d*u*w-A*C
  ff(8) = C*u-A*w
  do ivar = 9,nvar
#ifdef KROME
     if (ivar .ne. ichem) then
#endif
        ff(ivar) = d*u*qvar(ivar)
#ifdef KROME
     else  ! Make sure effective adiabatic index has zero flux
        ff(ivar) = zero
     endif
#endif
  end do

  call trace_mpi('find_mhd_flux(qvar,cvar,ff)','exit',2)
END SUBROUTINE find_mhd_flux 
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_info(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the fastest velocity at which information is exchanged 
  !! at the interface
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvar):: qvar  
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf

  call trace_mpi('find_speed_info(qvar,vel_info)','enter',4)

  d=qvar(1); P=qvar(2); u=qvar(3); A=qvar(4)
  v=qvar(5); B=qvar(6); w=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
#ifdef KROME
  c2 = merge(qvar(ichem),gamma,chemistry)*P/d
#else
  c2 = gamma*P/d
#endif
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(d2**2-c2*A*A/d) )
  vel_info = cf+abs(u)

  call trace_mpi('find_speed_info(qvar,vel_info)','exit',4)
END SUBROUTINE find_speed_info
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,gg,dx,dt,ncell,kmin,dt_loc,ilevel)
  use amr_commons, only: myid
  use amr_parameters
  use hydro_parameters
  use hydro_commons
  use const
  implicit none
  integer::ncell,kmin,ilevel
  real(kind=8) :: dt_loc
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar+3)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  real(dp),dimension(1:nvector)::a2,B2,rho,ctot
  
  real(dp)::dtcell,smallp,cf,cc,bc,bn
  integer::k,idim
  
  call trace_mpi('cmpdt(uu,gg,dx,dt,ncell,kmin,dt_loc)','enter',3)

  smallp = smallr*smallc**2/gamma

  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
     rho(k)=uu(k,1)
  end do
  do idim = 1,3
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/rho(k)
     end do
  end do

  do k = 1, ncell
     Bc = half*(uu(k,5+1)+uu(k,nvar+1))
     B2(k)=Bc**2
     uu(k,5) = uu(k,5)-half*uu(k,1)*uu(k,1+1)**2-half*Bc**2
  end do
  do idim = 2,3
     do k = 1, ncell
        Bc = half*(uu(k,5+idim)+uu(k,nvar+idim))
        B2(k)=B2(k)+Bc**2
        uu(k,5) = uu(k,5)-half*uu(k,1)*uu(k,idim+1)**2-half*Bc**2
     end do
  end do

  ! Compute thermal sound speed
  do k = 1, ncell
#ifdef KROME
     uu(k,5) = max((merge(uu(k,ichem),gamma,chemistry)-one)*uu(k,5),smallp)
     a2(k)=merge(uu(k,ichem),gamma,chemistry)*uu(k,5)/uu(k,1)
#else
     uu(k,5) = max((gamma-one)*uu(k,5),smallp)
     a2(k)=gamma*uu(k,5)/uu(k,1)
#endif
  end do

  ! Compute maximum wave speed (fast magnetosonic)
  do k = 1, ncell
     ctot(k)=zero
  end do
  if(scheme.eq.'induction')then
     do idim = 1,ndim   ! WARNING: ndim instead of 3  
        do k = 1, ncell
           ctot(k)=ctot(k)+abs(uu(k,idim+1))
        end do
     end do
  else
     do idim = 1,ndim   ! WARNING: ndim instead of 3  
        if (courant_type==0) then                                               ! legacy condition
           do k = 1, ncell
              cc=half*(B2(k)/rho(k)+a2(k))
              BN=half*(uu(k,5+idim)+uu(k,nvar+idim))
              cf=sqrt(cc+sqrt(max(cc**2-a2(k)*BN**2/rho(k),zero)))
              ctot(k)=ctot(k)+cf+abs(uu(k,idim+1))
           end do
        else if (courant_type==1) then                                          ! isotropic condition
           do k = 1, ncell
              cc=half*(B2(k)/rho(k)+a2(k))
              BN=half*(uu(k,5+idim)+uu(k,nvar+idim))
              cf=sqrt(cc+sqrt(max(cc**2-a2(k)*BN**2/rho(k),zero)))
              ctot(k)=ctot(k)+(cf+abs(uu(k,idim+1)))**2
           end do
        end if
     end do
     if (courant_type==2) then                                                  ! velocity + fast mode
        do idim = 1,ndim
           do k = 1, ncell
              ctot(k)=ctot(k)+uu(k,idim+1)**2
           end do
        end do
        do k = 1, ncell
           ctot(k)=sqrt(ctot(k))+sqrt(B2(k)/rho(k)+a2(k))
        end do
     endif
  endif
  if (courant_type==1) ctot(1:ncell)=sqrt(ctot(1:ncell))

  ! Compute gravity strength ratio
  do k = 1, ncell
     rho(k)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell 
        rho(k)=rho(k)+abs(Grho*gg(k,idim))
     end do
  end do
  do k = 1, ncell
     rho(k)=rho(k)*dx/ctot(k)**2
     rho(k)=MAX(rho(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  kmin = -1
  do k = 1,ncell
     dtcell=dx/ctot(k)*(sqrt(one+two*courant_factor*rho(k))-one)/rho(k)
     kmin = merge(kmin,k,dtcell >= dt)
     dt = min(dt,dtcell)
  end do
  if (kmin > -1) then
!$omp critical
    if (dt < dt_loc) then
       courant_number%cg = (sqrt(one+two*courant_factor*rho(kmin))-one)/rho(kmin)
       courant_number%cb = ctot(kmin)
       courant_number%cv = sqrt(sum(uu(kmin,2:4)**2))
       courant_number%cv = sqrt(sum(uu(kmin,2:4)**2))
       courant_number%rho = uu(kmin,1)
       courant_number%umag = sqrt(sum(uu(kmin,2:4)**2))
       courant_number%bmag = sqrt(B2(kmin))
       dt_loc=min(dt_loc,dt)
    endif
!$omp end critical
  endif
  
  call trace_mpi('cmpdt(uu,gg,dx,dt,ncell,kmin,dt_loc)','exit',3)
end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,rok,nn,ilevel)
  use amr_parameters
  use hydro_parameters
  use amr_commons, ONLY: emag_tot, n_refine_d, n_refine_p, n_refine_u, &
                         n_refine_c, n_refine_s, n_refine_b
  use const
  implicit none
  ! dummy arguments
  integer nn,ilevel,n_refine_dp,n_refine_pp
  real(dp), dimension(1:nvector,1:nvar) :: ug, um, ud
  real(dp)::rok(1:nvector)
 
  integer::k,idim
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind,rok0
  real(dp),dimension(1:nvector),save::emagg,emagm,emagd
  real(dp),dimension(1:nvector),save::irhog,irhom,irhod
  !$omp threadprivate(eking,ekinm,ekind,emagg,emagm,emagd,rok0,irhog,irhom,irhod)
  real(dp)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error,emag_loc,ethres,c_err,gam1
  real::n_criteria
  
  call trace_mpi('hydro_refine(ug,um,ud,rok,nn,ilevel)','enter',3)

  gam1 = gamma - one

  ! Convert to primitive variables
  !-----------------------------------
  ! Floor density
  do k = 1,nn
     ug(k,1) = max(ug(k,1),smallr)
     um(k,1) = max(um(k,1),smallr)
     ud(k,1) = max(ud(k,1),smallr)
     irhog(k)= 1.0_dp / ug(k,1)
     irhom(k)= 1.0_dp / um(k,1)
     irhod(k)= 1.0_dp / ud(k,1)
  end do
  ! Velocity
  !IBM* UNROLL(3)
  do idim = 1,3
     do k = 1,nn
        ug(k,idim+1) = ug(k,idim+1)*irhog(k)
        um(k,idim+1) = um(k,idim+1)*irhom(k)
        ud(k,idim+1) = ud(k,idim+1)*irhod(k)
     end do
  end do
  ! Pressure
  do k = 1,nn
     eking(k) = half*ug(k,1)*(ug(k,1+1)**2+ug(k,2+1)**2+ug(k,3+1)**2)
     ekinm(k) = half*um(k,1)*(um(k,1+1)**2+um(k,2+1)**2+um(k,3+1)**2)
     ekind(k) = half*ud(k,1)*(ud(k,1+1)**2+ud(k,2+1)**2+ud(k,3+1)**2)
  end do
  ! The magnetic field contained in 6:8 is already averaged over left and right face
  do k = 1,nn
     emagg(k) = half*( ug(k,1+5)**2 + ug(k,2+5)**2 + ug(k,3+5)**2 )
     emagm(k) = half*( um(k,1+5)**2 + um(k,2+5)**2 + um(k,3+5)**2 )
     emagd(k) = half*( ud(k,1+5)**2 + ud(k,2+5)**2 + ud(k,3+5)**2 )
  end do
  do k = 1,nn
#ifdef KROME
     ug(k,5) = (merge(ug(k,ichem),gamma,chemistry)-one)*(ug(k,5)-eking(k)-emagg(k))
     um(k,5) = (merge(um(k,ichem),gamma,chemistry)-one)*(um(k,5)-ekinm(k)-emagm(k))
     ud(k,5) = (merge(ud(k,ichem),gamma,chemistry)-one)*(ud(k,5)-ekind(k)-emagd(k))
#else
     ug(k,5) = gam1*(ug(k,5)-eking(k)-emagg(k))
     um(k,5) = gam1*(um(k,5)-ekinm(k)-emagm(k))
     ud(k,5) = gam1*(ud(k,5)-ekind(k)-emagd(k))
#endif
  end do  

  ! Compute normalization factor, which makes using one or more 
  ! criteria that are otherwise equivalent (e.g. d and p for isothermal)
  ! have a similar effect.  To avoid changing the effect relative to
  ! the many runs where these were used together, normalize with the
  ! number 4 -- typical of those runs.
  !-----------------------------------
  n_criteria=0.
  if (err_grad_d (ilevel)>0) n_criteria=n_criteria+1.
  if (err_grad_p (ilevel)>0) n_criteria=n_criteria+1.
  if (err_grad_u (ilevel)>0) n_criteria=n_criteria+1.
  if (err_grad_b2(ilevel)>0) n_criteria=n_criteria+1.
  n_criteria=n_criteria/2.

  ! Compute errors
  !-----------------------------------
  if(ilevel < levelmax_d .and. err_grad_d(ilevel) > 0.)then
     c_err = 1./err_grad_d(ilevel)/sqrt(n_criteria)
     n_refine_dp=0
     do k=1,nn
        dg=abs(ug(k,1)); dm=abs(um(k,1)); dd=abs(ud(k,1))
        error=c_err*MAX( &
               ABS(log((dd+floor_d)/(dm+floor_d))) , &
               ABS(log((dg+floor_d)/(dm+floor_d))) )
        rok0(k) = rok(k)
        rok(k) = rok(k)+error**2
        n_refine_dp = n_refine_dp + merge(1,0,rok(k)>1 .and. rok0(k)<1)
     end do
     !$omp atomic
     n_refine_d = n_refine_d + n_refine_dp
  end if

  if(ilevel < levelmax_p .and. err_grad_p(ilevel) > 0.)then
     c_err = 1./err_grad_p(ilevel)/sqrt(n_criteria)
     n_refine_pp=0
     do k=1,nn
        pg=abs(ug(k,5)); pm=abs(um(k,5)); pd=abs(ud(k,5))
        error=c_err*MAX( &
               ABS(log((pd+floor_p)/(pm+floor_p))) , &
               ABS(log((pg+floor_p)/(pm+floor_p))) )
        rok0(k) = rok(k)
        rok(k) = rok(k)+error**2
        n_refine_pp = n_refine_pp + merge(1,0,rok(k)>1 .and. rok0(k)<1)
     end do
     !$omp atomic
     n_refine_p = n_refine_p + n_refine_pp
  end if

  if(ilevel < levelmax_b .and. err_grad_b2(ilevel) > 0.)then
     c_err = 1./err_grad_b2(ilevel)/sqrt(n_criteria)
     n_refine_dp=0
     do k=1,nn
        pg=emagg(k); pm=emagm(k); pd=emagd(k)
        error=c_err*MAX( &
               ABS(log((pd+floor_b2)/(pm+floor_b2))) , &
               ABS(log((pg+floor_b2)/(pm+floor_b2))) )
        rok0(k) = rok(k)
        rok(k) = rok(k)+error**2
        n_refine_dp = n_refine_dp + merge(1,0,rok(k)>1 .and. rok0(k)<1)
     end do
     !$omp atomic
     n_refine_b = n_refine_b + n_refine_dp
  end if

  if(ilevel < levelmax_b .and. err_grad_A(ilevel) > 0.)then
     c_err = 2./err_grad_A(ilevel)/sqrt(n_criteria)
     idim = 1
     do k=1,nn
        vg=ug(k,5+idim)
        vm=um(k,5+idim)
        vd=ud(k,5+idim)
        cg=sqrt(emagg(k))
        cm=sqrt(emagm(k))
        cd=sqrt(emagd(k))
        error=c_err*MAX( &
             & ABS((vd-vm)/(cd+cm+floor_A)) , &
             & ABS((vm-vg)/(cm+cg+floor_A)) )
        rok(k) = rok(k)+error**2
     end do
  end if

  if(ilevel < levelmax_b .and. err_grad_B(ilevel) > 0.)then
     c_err = 2./err_grad_B(ilevel)/sqrt(n_criteria)
     idim = 2
     do k=1,nn
        vg=ug(k,5+idim)
        vm=um(k,5+idim)
        vd=ud(k,5+idim)
        cg=sqrt(emagg(k))
        cm=sqrt(emagm(k))
        cd=sqrt(emagd(k))
        error=c_err*MAX( &
             & ABS((vd-vm)/(cd+cm+floor_B)) , &
             & ABS((vm-vg)/(cm+cg+floor_B)) )
        rok(k) = rok(k)+error**2
     end do
  end if

  if(ilevel < levelmax_b .and. err_grad_C(ilevel) > 0.)then
     c_err = 2./err_grad_C(ilevel)/sqrt(n_criteria)
     idim = 3
     do k=1,nn
        vg=ug(k,5+idim)
        vm=um(k,5+idim)
        vd=ud(k,5+idim)
        cg=sqrt(emagg(k))
        cm=sqrt(emagm(k))
        cd=sqrt(emagd(k))
        error=c_err*MAX( &
             & ABS((vd-vm)/(cd+cm+floor_C)) , &
             & ABS((vm-vg)/(cm+cg+floor_C)) )
        rok(k) = rok(k)+error**2
     end do
  end if

  ! Passive scalars
  if(ilevel < levelmax_d .and. err_grad_s(ilevel) > 0.)then
     c_err = 1./err_grad_s(ilevel)/sqrt(n_criteria)
     n_refine_dp=0
     do k=1,nn
        rok0(k) = rok(k)
     end do
     do idim=9,nvar
        do k=1,nn
           dg=ug(k,idim)*irhog(k)
           dm=um(k,idim)*irhom(k)
           dd=ud(k,idim)*irhod(k)
           error=c_err*MAX( &
               ABS(log((dd+floor_s)/(dm+floor_s))) , &
               ABS(log((dg+floor_s)/(dm+floor_s))) )
           rok(k) = rok(k)+error**2
        end do
     end do
     do k=1,nn
        n_refine_dp = n_refine_dp + merge(1,0,rok(k)>1 .and. rok0(k)<1)
     end do
     !$omp atomic
     n_refine_s = n_refine_s + n_refine_dp
  end if

  if (ilevel < levelmax_u) then
   if(err_grad_u(ilevel) > 0.)then
     n_refine_dp=0
     c_err = 2./err_grad_u(ilevel)/sqrt(n_criteria)
     do k=1,nn
        rok0(k) = rok(k)
     end do
     do idim = 1,3
        do k=1,nn
           vg=ug(k,idim+1); vm=um(k,idim+1); vd=ud(k,idim+1)
#ifdef KROME
           cg=sqrt(max(merge(ug(k,ichem),gamma,chemistry)*ug(k,5)*irhog(k),floor_u**2))
           cm=sqrt(max(merge(um(k,ichem),gamma,chemistry)*um(k,5)*irhom(k),floor_u**2))
           cd=sqrt(max(merge(ud(k,ichem),gamma,chemistry)*ud(k,5)*irhod(k),floor_u**2))
#else
           cg=sqrt(max(gamma*ug(k,5)*irhog(k),floor_u**2))
           cm=sqrt(max(gamma*um(k,5)*irhom(k),floor_u**2))
           cd=sqrt(max(gamma*ud(k,5)*irhod(k),floor_u**2))
#endif
           error=c_err*MAX( &
                & ABS((vd-vm)/(cd+cm+ABS(vd)+ABS(vm)+floor_u)) , &
                & ABS((vm-vg)/(cm+cg+ABS(vm)+ABS(vg)+floor_u)) )
           rok(k) = rok(k)+error**2
        end do
     end do
     do k=1,nn
        n_refine_dp = n_refine_dp + merge(1,0,rok(k)>1 .and. rok0(k)<1)
     end do
     !$omp atomic
     n_refine_u = n_refine_u + n_refine_dp
   else if(err_grad_u(ilevel) < 0.)then
     n_refine_dp=0
     c_err = 2./abs(err_grad_u(ilevel))/sqrt(n_criteria)
     do k=1,nn
        rok0(k) = rok(k)
     end do
     do idim = 1,3
        do k=1,nn
           vg=ug(k,idim+1); vm=um(k,idim+1); vd=ud(k,idim+1)
#ifdef KROME
           cg=sqrt(max(merge(ug(k,ichem),gamma,chemistry)*ug(k,5)*irhog(k),floor_u**2))
           cm=sqrt(max(merge(um(k,ichem),gamma,chemistry)*um(k,5)*irhom(k),floor_u**2))
           cd=sqrt(max(merge(ud(k,ichem),gamma,chemistry)*ud(k,5)*irhod(k),floor_u**2))
#else
           cg=sqrt(max(gamma*ug(k,5)*irhog(k),floor_u**2))
           cm=sqrt(max(gamma*um(k,5)*irhom(k),floor_u**2))
           cd=sqrt(max(gamma*ud(k,5)*irhod(k),floor_u**2))
#endif
           error=c_err*MAX( &
                & ABS((vd-vm)/(cd+cm+floor_u)) , &
                & ABS((vm-vg)/(cm+cg+floor_u)) )
           rok(k) = rok(k)+error**2
        end do
     end do
     do k=1,nn
        n_refine_dp = n_refine_dp + merge(1,0,rok(k)>1 .and. rok0(k)<1)
     end do
     !$omp atomic
     n_refine_u = n_refine_u + n_refine_dp
   end if
  end if

  if(scheme.eq.'induction')then
     if(m_refine(ilevel) > 0.)then
        if(emag_tot==0.0)then
           emag_loc=1.5
        else
           emag_loc=emag_tot
        endif
        ethres=m_refine(ilevel)*emag_loc/boxlen**3
        do k=1,nn
           error = emagm(k)/ethres
           rok(k) = rok(k)+error**2
        end do
     endif
  endif

  call trace_mpi('hydro_refine(ug,um,ud,rok,nn,ilevel)','exit',3)
end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE upwind(qleft,qright,fgdnv,zero_flux)
  USE amr_parameters
  USE const
  USE hydro_parameters
  USE find_for_solvers
  ! 1D Upwind Riemann solver
  IMPLICIT NONE  
  REAL(dp)::zero_flux
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv

  REAL(dp),DIMENSION(1:nvar)::fleft,fright,fmean
  REAL(dp),DIMENSION(1:nvar)::uleft,uright,udiff
  REAL(dp):: vleft,vright,bx_mean
  INTEGER ::l

  call trace_mpi('upwind(qleft,qright,fgdnv,zero_flux)','enter',2)
  
  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright)
  
  ! find the mean flux
  fmean =  half * ( fright + fleft ) * zero_flux
  
  ! find the mean normal velocity
  vleft = half * ( qleft(3) + qright(3) )
  
  ! difference between the 2 states
  udiff = half * ( uright - uleft )
  
  ! the Upwind flux
  fgdnv = fmean - ABS(vleft) * udiff

  call trace_mpi('upwind(qleft,qright,fgdnv,zero_flux)','exit',2)
END SUBROUTINE upwind
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE lax_friedrich_vec(qleft,qright,fgdnv,zero_flux)
  USE amr_parameters
  USE const
  USE hydro_parameters
  USE find_for_solvers
  ! 1D local Lax-Friedrich Riemann solver
  IMPLICIT NONE  
  REAL(dp), intent(in) ::zero_flux
  REAL(dp),DIMENSION(1:nvector,1:nvar),target ::qleft,qright,fgdnv

  REAL(dp),DIMENSION(1:nvector,1:nvar)::fleft,fright
  REAL(dp),DIMENSION(1:nvector,1:nvar)::uleft,uright
  REAL(dp),DIMENSION(1:nvector):: vleft,vright
  REAL(dp):: hz 
  INTEGER ::i,l

  call trace_mpi('lax_friedrich_vec(qleft,qright,fgdnv,zero_flux)','enter',3)
  
  ! Enforce continuity of normal component
  qleft (:,4)=half*(qleft(:,4)+qright(:,4))
  qright(:,4)=qleft(:,4)
  
  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_info_vec(qleft ,vleft )
  CALL find_speed_info_vec(qright,vright)

  ! find the flux vector
  CALL find_mhd_flux_vec(qleft ,uleft ,fleft ,vleft)
  CALL find_mhd_flux_vec(qright,uright,fright,uleft)
    
  if (zero_flux .ne. 0.) then
    vleft = - half * MAX(vleft,vright)
    hz    = half*zero_flux
    do i=1,nvar
    do l=1,nvector
     ! the local Lax-Friedrich flux
     fgdnv(l,i) = hz * ( fright(l,i) + fleft(l,i) ) + vleft(l) * ( uright(l,i) - uleft(l,i) )
    enddo
    enddo
  else
    vleft = - half * MAX(vleft,vright)
    do i=1,nvar
    do l=1,nvector
     ! the local Lax-Friedrich flux
     fgdnv(l,i) = vleft(l) * ( uright(l,i) - uleft(l,i) )
    enddo
    enddo  
  endif
  call trace_mpi('lax_friedrich_vec(qleft,qright,fgdnv,zero_flux)','exit',3)
END SUBROUTINE lax_friedrich_vec
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE lax_friedrich(qleft,qright,fgdnv,zero_flux)
  USE amr_parameters
  USE const
  USE hydro_parameters
  ! 1D local Lax-Friedrich Riemann solver
  IMPLICIT NONE  
  REAL(dp)::zero_flux
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv

  REAL(dp),DIMENSION(1:nvar)::fleft,fright,fmean
  REAL(dp),DIMENSION(1:nvar)::uleft,uright,udiff
  REAL(dp):: vleft,vright,bx_mean
  INTEGER ::l

  call trace_mpi('lax_friedrich(qleft,qright,fgdnv,zero_flux)','enter',3)
  
  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright) 
  
  ! find the mean flux
  fmean =  half * ( fright + fleft ) * zero_flux
  
  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_info(qleft ,vleft )
  CALL find_speed_info(qright,vright)
  
  ! difference between the 2 states
  udiff  = half * ( uright - uleft )
  
  ! the local Lax-Friedrich flux
  fgdnv = fmean - MAX(vleft,vright) * udiff

  call trace_mpi('lax_friedrich(qleft,qright,fgdnv,zero_flux)','exit',3)
END SUBROUTINE lax_friedrich
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE hll(qleft,qright,fgdnv)
  USE amr_parameters
  USE const
  USE hydro_parameters
  USE find_for_solvers
  ! 1D HLL Riemann solver
  IMPLICIT NONE  
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv
  REAL(dp),DIMENSION(1:nvar)::fleft,fright
  REAL(dp),DIMENSION(1:nvar)::uleft,uright,udiff
  REAL(dp):: vleft,vright,bx_mean,cfleft,cfright,SL,SR
  INTEGER ::l

  call trace_mpi('hll(qleft,qright,fgdnv)','enter',3)
  
  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright) 
  
  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_sound(qleft ,cfleft )
  CALL find_speed_sound(qright,cfright)
  vleft =qleft (3)
  vright=qright(3)
  SL=min(min(vleft,vright)-max(cfleft,cfright),zero)  ! the largest negative wave speed
  SR=max(max(vleft,vright)+max(cfleft,cfright),zero)  ! the largest positive wave speed

  ! The HLL flux is constructed below.  If SL = -SR and fleft=fright=flux
  ! then the first two terms give a flux equal to the average of the fluxes
  ! at the left and right interface.  The last term gives a diffusive
  ! flux which attempts to make uright=uleft. This part of the flux scales
  ! with the magnitude of SR and -SL. In cases where SR and -SL are
  ! essentially +- a large Alfven speeds the result is a huge diffusive flux.
  ! This flux attacks, for example, a rapid change in density, when the 
  ! gradient is perpendicular to the magnetic field.  This is in the first
  ! place an extremely unphysical flux, which also does not conserve entropy.
  ! A static or slowly moving contact discontinuity will be very quickly 
  ! washed out by this term, for no good reason.  The term would be useful
  ! and relevant only to avoid oscillations in fast mode waves.

  fgdnv    = (SR*fleft   -SL*fright   +0.25*SR*SL*(uright   -uleft   ))/(SR-SL)
  fgdnv(6) = (SR*fleft(6)-SL*fright(6)+0.50*SR*SL*(uright(6)-uleft(6)))/(SR-SL)
  fgdnv(8) = (SR*fleft(8)-SL*fright(8)+0.50*SR*SL*(uright(8)-uleft(8)))/(SR-SL)

  call trace_mpi('hll(qleft,qright,fgdnv)','exit',3)
END SUBROUTINE hll
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE hlld_vec(qleft,qright,fgdnv)
  USE amr_parameters,   ONLY: chemistry
  USE hydro_parameters, ONLY: nvar, gamma, cmax2, ichem
  USE const
  USE find_for_solvers
  ! HLLD Riemann solver (Miyoshi & Kusano, 2005, JCP, 208, 315)
  IMPLICIT NONE
  REAL(dp),DIMENSION(1:nvector,1:nvar), target ::qleft,qright
  REAL(dp),DIMENSION(1:nvector,1:nvar)         ::fgdnv
  REAL(dp),DIMENSION(1:nvector)::SL,SR,SAL,SAR,A,sgnm
  REAL(dp)::entho
  REAL(dp),DIMENSION(:), pointer ::rl,pl,ul,vl,wl,cl,bl
  REAL(dp),DIMENSION(1:nvector)::ecinl,emagl,etotl,ptotl,vdotbl,el
  REAL(dp),DIMENSION(:), pointer ::rr,pr,ur,vr,wr,cr,br
  REAL(dp),DIMENSION(1:nvector)::ecinr,emagr,etotr,ptotr,vdotbr,er
  REAL(dp),DIMENSION(1:nvector)::cfastl,calfvenl,rcl,rstarl,vstarl,wstarl,bstarl,cstarl,vdotbstarl
  REAL(dp),DIMENSION(1:nvector)::cfastr,calfvenr,rcr,rstarr,vstarr,wstarr,bstarr,cstarr,vdotbstarr
  REAL(dp),DIMENSION(1:nvector)::sqrrstarl,etotstarl,etotstarstarl
  REAL(dp),DIMENSION(1:nvector)::sqrrstarr,etotstarr,etotstarstarr
  REAL(dp),DIMENSION(1:nvector)::ustar,ptotstar,estar,vstarstar,wstarstar,bstarstar,cstarstar,vdotbstarstar
  REAL(dp),DIMENSION(1:nvector)::ro,uo,vo,wo,bo,co,ptoto,etoto,vdotbo
  REAL(dp),DIMENSION(1:nvector)::lim,emago
#ifdef KROME
  REAL(dp),DIMENSION(:), pointer ::gammal,gammar
#endif

  INTEGER ::ivar,i
  
  call trace_mpi('hlld_vec(qleft,qright,fgdnv)','enter',2)
 
  entho = one/(gamma-one)

  ! Enforce continuity of normal component
  A=half*(qleft(:,4)+qright(:,4))
  sgnm=sign(one,A)
  qleft(:,4)=A; qright(:,4)=A

#ifdef KROME
  ! Left and right adiabatic indexes
  gammal=>qleft(:,ichem); gammar=>qright(:,ichem)
#endif

  ! Left variables
  rl=>qleft(:,1); Pl=>qleft(:,2); ul=>qleft(:,3)
  vl=>qleft(:,5); Bl=>qleft(:,6); wl=>qleft(:,7); Cl=>qleft(:,8)
  ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
  emagl = half*(A*A+Bl*Bl+Cl*Cl)
#ifdef KROME
  etotl = Pl/(merge(gammal,gamma,chemistry)-one)+ecinl+emagl
#else
  etotl = Pl*entho+ecinl+emagl
#endif
  Ptotl = Pl + emagl
  vdotBl= ul*A+vl*Bl+wl*cl

  ! Right variables
  rr=>qright(:,1); Pr=>qright(:,2); ur=>qright(:,3)
  vr=>qright(:,5); Br=>qright(:,6); wr=>qright(:,7); Cr=>qright(:,8)
  ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
  emagr = half*(A*A+Br*Br+Cr*Cr)
#ifdef KROME
  etotr = Pr/(merge(gammar,gamma,chemistry)-one)+ecinr+emagr
#else
  etotr = Pr*entho+ecinr+emagr
#endif
  Ptotr = Pr + emagr
  vdotBr= ur*A+vr*Br+wr*Cr

  ! Find the largest eigenvalues in the normal direction to the interface
  CALL find_speed_fast_vec(qleft ,cfastl)
  CALL find_speed_fast_vec(qright,cfastr)

  ! Compute HLL wave speed
  SL=min(ul,ur)-max(cfastl,cfastr)
  SR=max(ul,ur)+max(cfastl,cfastr)

  if (cmax2 > 0) then
    lim = max(abs(SL),abs(SR)) ! use max speed for current limiter
    lim = cmax2 / (cmax2 + lim*lim)
  else
    lim = 1.0_dp
  endif

!  SL=ul-cfastl
!  SR=ur+cfastr
  
  ! Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)

  ! Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

  ! Left star region variables
  rstarl=rl*(SL-ul)/(SL-ustar)
  estar =rl*(SL-ul)*(SL-ustar)-A**2
  el    =rl*(SL-ul)*(SL-ul   )-A**2
  where(abs(estar) < 1d-4*A**2)
     vstarl=vl
     Bstarl=Bl
     wstarl=wl
     Cstarl=Cl
  elsewhere
     vstarl=vl-A*Bl*(ustar-ul)/estar
     Bstarl=Bl*el/estar
     wstarl=wl-A*Cl*(ustar-ul)/estar
     Cstarl=Cl*el/estar
  endwhere
  vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
  etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
  sqrrstarl=sqrt(rstarl)
  calfvenl=abs(A)/sqrrstarl
  SAL=ustar-calfvenl

  ! Right star region variables
  rstarr=rr*(SR-ur)/(SR-ustar)
  estar =rr*(SR-ur)*(SR-ustar)-A**2
  er    =rr*(SR-ur)*(SR-ur   )-A**2
  where(abs(estar) < 1d-4*A**2)
     vstarr=vr
     Bstarr=Br
     wstarr=wr
     Cstarr=Cr
  elsewhere
     vstarr=vr-A*Br*(ustar-ur)/estar
     Bstarr=Br*er/estar
     wstarr=wr-A*Cr*(ustar-ur)/estar
     Cstarr=Cr*er/estar
  endwhere
  vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
  sqrrstarr=sqrt(rstarr)
  calfvenr=abs(A)/sqrrstarr
  SAR=ustar+calfvenr

  ! Double star region variables
  vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
  wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
  Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
  Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
  vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
  etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
  etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)
  
  ! Sample the solution at x/t=0
  do i=1,nvector
    if(SL(i)>zero)then
       ro(i)=rl(i)
       uo(i)=ul(i)
       vo(i)=vl(i)
       wo(i)=wl(i)
       Bo(i)=Bl(i)
       Co(i)=Cl(i)
       Ptoto(i)=Ptotl(i)
       etoto(i)=etotl(i)
       vdotBo(i)=vdotBl(i)
    else if(SAL(i)>zero)then
       ro(i)=rstarl(i)
       uo(i)=ustar(i)
       vo(i)=vstarl(i)
       wo(i)=wstarl(i)
       Bo(i)=Bstarl(i)
       Co(i)=Cstarl(i)
       Ptoto(i)=Ptotstar(i)
       etoto(i)=etotstarl(i)
       vdotBo(i)=vdotBstarl(i)
    else if(ustar(i)>zero)then
       ro(i)=rstarl(i)
       uo(i)=ustar(i)
       vo(i)=vstarstar(i)
       wo(i)=wstarstar(i)
       Bo(i)=Bstarstar(i)
       Co(i)=Cstarstar(i)
       Ptoto(i)=Ptotstar(i)
       etoto(i)=etotstarstarl(i)
       vdotBo(i)=vdotBstarstar(i)
    else if(SAR(i)>zero)then
       ro(i)=rstarr(i)
       uo(i)=ustar(i)
       vo(i)=vstarstar(i)
       wo(i)=wstarstar(i)
       Bo(i)=Bstarstar(i)
       Co(i)=Cstarstar(i)
       Ptoto(i)=Ptotstar(i)
       etoto(i)=etotstarstarr(i)
       vdotBo(i)=vdotBstarstar(i)
    else if (SR(i)>zero)then
       ro(i)=rstarr(i)
       uo(i)=ustar(i)
       vo(i)=vstarr(i)
       wo(i)=wstarr(i)
       Bo(i)=Bstarr(i)
       Co(i)=Cstarr(i)
       Ptoto(i)=Ptotstar(i)
       etoto(i)=etotstarr(i)
       vdotBo(i)=vdotBstarr(i)
    else
       ro(i)=rr(i)
       uo(i)=ur(i)
       vo(i)=vr(i)
       wo(i)=wr(i)
       Bo(i)=Br(i)
       Co(i)=Cr(i)
       Ptoto(i)=Ptotr(i)
       etoto(i)=etotr(i)
       vdotBo(i)=vdotBr(i)
    end if
  end do

  emago=half*(A*A+Bo*Bo+Co*Co)

  ! Compute the Godunov flux
  fgdnv(:,1) = ro*uo
  fgdnv(:,2) = (etoto+Ptoto)*uo-A*vdotBo
  fgdnv(:,3) = ro*uo*uo+Ptoto+emago*(lim-1.0_dp)-A*A*lim
  fgdnv(:,4) = zero
  fgdnv(:,5) = ro*uo*vo-A*Bo*lim
  fgdnv(:,6) = Bo*uo-A*vo
  fgdnv(:,7) = ro*uo*wo-A*Co*lim
  fgdnv(:,8) = Co*uo-A*wo
  do ivar = 9,nvar
     where(fgdnv(:,1)>0)
        fgdnv(:,ivar) = fgdnv(:,1)*qleft (:,ivar)
     elsewhere
        fgdnv(:,ivar) = fgdnv(:,1)*qright(:,ivar)
     endwhere
  end do
#ifdef KROME
  fgdnv(:,ichem) = 0.0_dp
#endif

  call trace_mpi('hlld_vec(qleft,qright,fgdnv)','exit',3)
END SUBROUTINE hlld_vec
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE hlld(qleft,qright,fgdnv,mark_low_order,dbg,xx)
  USE amr_parameters,   ONLY: chemistry
  USE hydro_parameters, ONLY: nvar, gamma, cmax2, cbmax, ichem
  USE const
  ! HLLD Riemann solver (Miyoshi & Kusano, 2005, JCP, 208, 315)
  IMPLICIT NONE
  REAL(dp),DIMENSION(1:nvar)::qleft,qright,fgdnv
  LOGICAL::mark_low_order
  REAL(dp)::SL,SR,SAL,SAR
  REAL(dp)::entho,A,sgnm
  REAL(dp)::rl,pl,ul,vl,wl,cl,ecinl,emagl,etotl,ptotl,vdotbl,bl,el
  REAL(dp)::rr,pr,ur,vr,wr,cr,ecinr,emagr,etotr,ptotr,vdotbr,br,er
  REAL(dp)::cfastl,calfvenl,rcl,rstarl,vstarl,wstarl,bstarl,cstarl,vdotbstarl
  REAL(dp)::cfastr,calfvenr,rcr,rstarr,vstarr,wstarr,bstarr,cstarr,vdotbstarr
  REAL(dp)::sqrrstarl,etotstarl,etotstarstarl
  REAL(dp)::sqrrstarr,etotstarr,etotstarstarr
  REAL(dp)::ustar,ptotstar,estar,vstarstar,wstarstar,bstarstar,cstarstar,vdotbstarstar
  REAL(dp)::ro,uo,vo,wo,bo,co,ptoto,etoto,vdotbo
  REAL(dp)::lim,emago,xx
  LOGICAL dbg

  INTEGER ::ivar,icase
  
  call trace_mpi('hlld','enter',3)
  entho = one/(gamma-one)

  ! Enforce continuity of normal component
  A=half*(qleft(4)+qright(4))
  sgnm=sign(one,A)
  qleft(4)=A; qright(4)=A

  ! Left variables
  rl=qleft(1); Pl=qleft(2); ul=qleft(3)
  vl=qleft(5); Bl=qleft(6); wl=qleft(7); Cl=qleft(8)
  ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
  emagl = half*(A*A+Bl*Bl+Cl*Cl)
#ifdef KROME
  etotl = Pl/(merge(qleft(ichem),gamma,chemistry)-one)+ecinl+emagl
#else
  etotl = Pl*entho+ecinl+emagl
#endif
  Ptotl = Pl + emagl
  vdotBl= ul*A+vl*Bl+wl*cl

  !if (dbg) print 1,'hlld1',rl,vl,ecinl,emagl,etotl,Ptotl,vdotBl
1 format(a,1p,12(g13.6,1h,))

  ! Right variables
  rr=qright(1); Pr=qright(2); ur=qright(3)
  vr=qright(5); Br=qright(6); wr=qright(7); Cr=qright(8)
  ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
  emagr = half*(A*A+Br*Br+Cr*Cr)
#ifdef KROME
  etotr = Pr/(merge(qright(ichem),gamma,chemistry)-one)+ecinr+emagr
#else
  etotr = Pr*entho+ecinr+emagr
#endif
  Ptotr = Pr + emagr
  vdotBr= ur*A+vr*Br+wr*Cr

  !if (dbg) print 1,'hlld2',rr,vr,ecinr,emagr,etotr,Ptotr,vdotBr

  ! Find the largest eigenvalues in the normal direction to the interface
  CALL find_speed_fast(qleft ,cfastl)
  CALL find_speed_fast(qright,cfastr)

  !if (dbg) print 1,'hlld3', cfastl, cfastr

  ! Compute HLL wave speed
  SL=min(ul,ur)-max(cfastl,cfastr)
  SR=max(ul,ur)+max(cfastl,cfastr)

  if (cbmax > 0 .and. mark_low_order) then ! Check if we have to default to llf
    call lax_friedrich(qleft,qright,fgdnv,one)
    call trace_mpi('hlld','exit',3)
    return
  endif

  if (cmax2 > 0) then
    lim = max(abs(SL),abs(SR)) ! use max speed for current limiter
    lim = cmax2 / (cmax2 + lim*lim)
  else
    lim = 1.0_dp
  endif

  ! SL=ul-cfastl
  ! SR=ur+cfastr
  
  ! Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)

  ! Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

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
  !if (dbg) print 1,'hllD',estar,vstarr,vr,A,Br,ustar,ur,estar

  vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
  !if (dbg) print 1,'hlldC',ustar,A,vstarr,Bstarr,wstarr,Cstarr

  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
  !if (dbg) print 1,'hlldB',SR,ur,etotr,Ptotr,Ptotstar,ustar,A,vdotBr,vdotBstarr

  sqrrstarr=sqrt(rstarr)
  calfvenr=abs(A)/sqrrstarr
  SAR=ustar+calfvenr

  ! Double star region variables
  vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
  wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
  Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
  Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
  vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
  etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
  etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)

  !if (dbg) print 1,'hlldA',etotstarstarr,etotstarr,sgnm,sqrrstarr,vdotBstarr,vdotBstarstar


  !if (dbg) print 1,'hlld5', SL,SAL,ustar,SAR,SR
  
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

  ! Compute the Godunov flux
  fgdnv(1) = ro*uo
  fgdnv(2) = (etoto+Ptoto)*uo-A*vdotBo
  if (dbg) print 2,'hlld7', icase,xx,fgdnv(2),etoto,Ptoto,uo,A,vdotBo
2 format(a,i2,f10.6,1p,12(g13.6))
  fgdnv(3) = ro*uo*uo+Ptoto+emago*(lim-1.0_dp)-A*A*lim
  fgdnv(4) = zero
  fgdnv(5) = ro*uo*vo-A*Bo*lim
  fgdnv(6) = Bo*uo-A*vo
  fgdnv(7) = ro*uo*wo-A*Co*lim
  fgdnv(8) = Co*uo-A*wo
  do ivar = 9,nvar
     if(fgdnv(1)>0)then
        fgdnv(ivar) = fgdnv(1)*qleft (ivar)
     else
        fgdnv(ivar) = fgdnv(1)*qright(ivar)
     endif
  end do
#ifdef KROME
  fgdnv(ichem) = 0.0_dp
#endif

  call trace_mpi('hlld','exit',3)
END SUBROUTINE hlld
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_alfven(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the alfven velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal, 
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE
   
  INTEGER :: i ,ib, iv
  REAL(dp),DIMENSION(1:nvar):: qvar  
  REAL(dp) :: vel_info
  REAL(dp) :: d,P,u,v,w,A,B,C,B2,c2,d2,cf

  d=qvar(1); A=qvar(4)
  vel_info = sqrt(A*A/d)

END SUBROUTINE find_speed_alfven
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE athena_roe(qleft,qright,fmean,zero_flux)
  USE amr_parameters,   only : dp, chemistry
  use amr_commons,      only : myid
  USE const,            only : half, zero, two
  USE hydro_parameters, only : nvar
  IMPLICIT NONE

  REAL(dp),DIMENSION(1:nvar)::qleft,qright
  REAL(dp),DIMENSION(1:nvar)::fleft,fright,fmean,qmean,uleft,uright,udiff
  REAL(dp)                  ::vleft,vright

  REAL(dp), DIMENSION(7,7) :: lem, rem
  REAL(dp), DIMENSION(7)        :: lambda, lambdal, lambdar, a

  REAL(dp) :: droe, vxroe, vyroe, vzroe, sqrtdl, sqrtdr
  REAL(dp) :: hroe, byroe, bzroe, Xfactor, Yfactor
  REAL(dp) :: pbr, pbl

  REAL(dp) :: fluxd, fluxmx, fluxmy, fluxmz, fluxe
  REAL(dp) :: fluxby, fluxbz, byl, byr, bzl, bzr, bx
  REAL(dp) :: dl, vxl, vyl, vzl, pl, hl
  REAL(dp) :: dr, vxr, vyr, vzr, pr, hr
  REAL(dp) :: mxl, myl, mzl, el
  REAL(dp) :: mxr, myr, mzr, er
  REAL(dp) :: coef, bx_mean

  REAL(dp) :: zero_flux,dim,pim,eim,mxm,mym,mzm,bym,bzm,etm,l1,l2

  INTEGER :: i, n, m
  LOGICAL :: llf
  
  call trace_mpi('athena_roe(qleft,qright,fmean,zero_flux)','enter',2)

#ifdef KROME
  if (chemistry) then
     if (myid==1) print *, 'athena roe solver is not supported with variable gamma!'
     stop
  endif
#endif

  ! Enforce continuity of normal component
  bx_mean=0.5d0*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  ! compute the fluxes and the conserved variables
  ! remember the convention : rho, E, rhovx, bx, rhovy, by, vz, rhobz
  CALL find_mhd_flux(qleft ,uleft ,fleft )
  CALL find_mhd_flux(qright,uright,fright)

  ! define the primitive quantities explicitly 
  dl  = qleft (1)
  dr  = qright(1)

  pl  = qleft (2)
  pr  = qright(2)

  vxl = qleft (3)
  vxr = qright(3)

  vyl = qleft (5)
  vyr = qright(5)

  byl = qleft (6)
  byr = qright(6)

  vzl = qleft (7)
  vzr = qright(7)

  bzl = qleft (8)
  bzr = qright(8)

  ! finally attribute bx value (left and right are identical)
  bx  = 0.5d0*(qleft(4)+qright(4))

  ! define explicitly the conserved quantities
  el  = uleft (2)
  er  = uright(2)

  mxl = uleft (3)
  mxr = uright(3)

  myl = uleft (5)
  myr = uright(5)

  mzl = uleft (7)
  mzr = uright(7)

  ! the magnetic pressure
  pbl = half * (bx*bx + byl*byl + bzl*bzl)
  pbr = half * (bx*bx + byr*byr + bzr*bzr)
  
  ! the total (specific) enthalpy
  hl = (el+pl+pbl)/dl
  hr = (er+pr+pbr)/dr
  !
  ! Step 2 : Compute Roe-averaged data from left and right states
  !
  sqrtdl  = sqrt(dl)
  sqrtdr  = sqrt(dr)
  droe    = sqrtdl*sqrtdr
  vxroe   = (sqrtdl*vxl + sqrtdr*vxr)/(sqrtdl+sqrtdr)
  vyroe   = (sqrtdl*vyl + sqrtdr*vyr)/(sqrtdl+sqrtdr)
  vzroe   = (sqrtdl*vzl + sqrtdr*vzr)/(sqrtdl+sqrtdr)
  byroe   = (sqrtdr*byl + sqrtdl*byr)/(sqrtdl+sqrtdr)
  bzroe   = (sqrtdr*bzl + sqrtdl*bzr)/(sqrtdl+sqrtdr)
  hroe    = (sqrtdl*hl  + sqrtdr*hr )/(sqrtdl+sqrtdr)
  Xfactor = ((byroe*byroe-byl*byr)+(bzroe*bzroe-bzl*bzr))/(2*droe)
  Yfactor = (dl+dr)/(2*droe)
  !
  ! Step 3 : Compute eigenvalues and eigenmatrices from Roe-averaged values
  !
  call eigen_cons(droe,vxroe,vyroe,vzroe,hroe,bx,byroe,bzroe,Xfactor,Yfactor,lambda,rem,lem)
  !
  ! Step 4: Compute eigenvalues from left and right states
  !  
  call eigenvalues(dl,vxl,vyl,vzl,pl,bx,byl,bzl,lambdal)
  call eigenvalues(dr,vxr,vyr,vzr,pr,bx,byr,bzr,lambdar)
  !
  ! Step 5 : Create intermediate states from eigenmatrices
  !
  do n = 1, 7
     a(n)  = 0.0
     a(n)  = a(n)  + (dr -dl ) * lem(1,n)
     a(n)  = a(n)  + (mxr-mxl) * lem(2,n)
     a(n)  = a(n)  + (myr-myl) * lem(3,n)
     a(n)  = a(n)  + (mzr-mzl) * lem(4,n)
     a(n)  = a(n)  + (er - el) * lem(5,n)
     a(n)  = a(n)  + (byr-byl) * lem(6,n)
     a(n)  = a(n)  + (bzr-bzl) * lem(7,n)
  end do

  llf = .false.
  dim = dl 
  mxm = mxl
  mym = myl
  mzm = mzl
  eim = el 
  bym = byl
  bzm = bzl
  do n = 1, 7
     dim = dim + a(n) * rem(n,1)
     mxm = mxm + a(n) * rem(n,2)
     mym = mym + a(n) * rem(n,3)
     mzm = mzm + a(n) * rem(n,4)
     eim = eim + a(n) * rem(n,5)
     bym = bym + a(n) * rem(n,6)
     bzm = bzm + a(n) * rem(n,7)
     etm = eim-0.5*(mxm*mxm+mym*mym+mzm*mzm)/dim-0.5*(bx*bx+bym*bym+bzm*bzm)
     if(dim.le.zero.or.etm.le.zero)then
        llf=.true.
     endif
  end do

  IF( llf ) THEN
     fmean = half * ( fright + fleft ) * zero_flux
     CALL find_speed_info(qleft ,vleft )
     CALL find_speed_info(qright,vright)
     udiff = half * ( uright - uleft )
     fmean = fmean - MAX(vleft,vright) * udiff
     RETURN
  END IF
  !
  ! Step 6 : Entropy fix for genuinely non linear waves
  !
  do n = 1, 7, 2
     l1 = min(lambdal(n),lambda(n))
     l2 = max(lambdar(n),lambda(n))
     if(l1.lt.zero.and.l2.gt.zero)then
        lambda(n)=(lambda(n)*(l2+l1)-two*l2*l1)/(l2-l1)
     endif
  end do
  !
  ! Step 6 : Compute fluxes at interface using  Roe  solver
  !
  ! add the left and right fluxes
  ! remember the convention : rho, E, rhovx, bx, rhovy, by, vz, rhobz
  fleft  = fleft  * zero_flux
  fright = fright * zero_flux
  fluxd  = fleft(1) + fright(1)
  fluxe  = fleft(2) + fright(2)
  fluxmx = fleft(3) + fright(3)
  fluxmy = fleft(5) + fright(5)
  fluxby = fleft(6) + fright(6)
  fluxmz = fleft(7) + fright(7)
  fluxbz = fleft(8) + fright(8)

  ! now compute the Roe fluxes
  ! remember that the convention of athena's eigenvalues is : 
  ! rho,rhovx,rhovy,rhovz,E,by,bz
  DO n=1,7
      coef = ABS(lambda(n))*a(n)
      fluxd   = fluxd  - coef*rem(n,1)
      fluxe   = fluxe  - coef*rem(n,5)
      fluxmx  = fluxmx - coef*rem(n,2)
      fluxmy  = fluxmy - coef*rem(n,3)
      fluxby  = fluxby - coef*rem(n,6)
      fluxmz  = fluxmz - coef*rem(n,4)
      fluxbz  = fluxbz - coef*rem(n,7)
  ENDDO
  ! take half and put into the fmean variables 
  fmean(1) = half * fluxd 
  fmean(2) = half * fluxe 
  fmean(3) = half * fluxmx
  fmean(4) = zero
  fmean(5) = half * fluxmy
  fmean(6) = half * fluxby
  fmean(7) = half * fluxmz
  fmean(8) = half * fluxbz
  DO n=9,nvar
     if(fmean(1)>0)then
        fmean(n)=qleft (n)*fmean(1)
     else
        fmean(n)=qright(n)*fmean(1)
     endif
  END DO
  call trace_mpi('athena_roe(qleft,qright,fmean,zero_flux)','exit',2)
END SUBROUTINE athena_roe
!###########################################################
!###########################################################
!###########################################################
!########################################################### 
SUBROUTINE hydro_acoustic(qleft,qright,fgdnv)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  real(dp),dimension(1:nvar)::qleft,qright,fgdnv,qgdnv,ugdnv

  ! local variables
  integer::i,n
  real(dp)::smallp, bx_mean
  real(dp)::rl   ,ul   ,pl   ,cl
  real(dp)::rr   ,ur   ,pr   ,cr   
  real(dp)::ro   ,uo   ,po   ,co   
  real(dp)::rstar,ustar,pstar,cstar
  real(dp)::wl   ,wr   ,wo   
  real(dp)::sgnm ,spin ,spout,ushock
  real(dp)::frac

  call trace_mpi('hydro_acoustic(qleft,qright,fgdnv)','enter',2)

  ! constants
  smallp = smallr*smallc**2

  ! Enforce continuity of normal component
  bx_mean=half*(qleft(4)+qright(4))
  qleft (4)=bx_mean
  qright(4)=bx_mean

  ! Initial states pressure, density and velocity
  rl=max(qleft (1),smallr)
  rr=max(qright(1),smallr)

  pl=max(qleft (2),smallp)
  pr=max(qright(2),smallp)

  ul = qleft (3)
  ur = qright(3)

  ! Acoustic star state
#ifdef KROME
  cl = sqrt(merge(qleft (ichem),gamma,chemistry)*pl/rl)
  cr = sqrt(merge(qright(ichem),gamma,chemistry)*pr/rr)
#else
  cl = sqrt(gamma*pl/rl)
  cr = sqrt(gamma*pr/rr)
#endif
  wl = cl*rl
  wr = cr*rr
  pstar = ((wr*pl+wl*pr)+wl*wr*(ul-ur))/(wl+wr)
  ustar = ((wr*ur+wl*ul)+(pl-pr))/(wl+wr)

  ! Left going or right going contact wave
  sgnm = sign(one,ustar)

  ! Left or right unperturbed state
  if(sgnm==one)then
     ro = rl
     uo = ul
     po = pl
     wo = wl
     co = cl
  else
     ro = rr
     uo = ur
     po = pr
     wo = wr
     co = cr
  end if

  ! Star region density and sound speed
  rstar = ro+(pstar-po)/co**2
  rstar = max(rstar,smallr)
#ifdef KROME
  if (chemistry) then
     cstar = sqrt(abs(merge(qleft(ichem),qright(ichem),sgnm==one)*pstar/rstar))
  else
     cstar = sqrt(abs(gamma*pstar/rstar))
  endif
#else
  cstar = sqrt(abs(gamma*pstar/rstar))
#endif
  cstar = max(cstar,smallc)

  ! Head and tail speed of rarefaction
  spout = co   -sgnm*uo   
  spin  = cstar-sgnm*ustar

  ! Shock speed
  ushock = half*(spin+spout)
  ushock = max(ushock,-sgnm*ustar)
  if(pstar>=po)then
     spout=ushock
     spin =spout 
  end if

  ! Sample the solution at x/t=0
  if(spout<zero)then      ! Initial state
     qgdnv(1) = ro
     qgdnv(2) = po
     qgdnv(3) = uo
  else if(spin>=zero)then  ! Star region
     qgdnv(1) = rstar
     qgdnv(2) = pstar
     qgdnv(3) = ustar
  else                        ! Rarefaction
     frac = spout/(spout-spin)
     qgdnv(1) = frac*rstar + (one - frac)*ro
     qgdnv(2) = frac*pstar + (one - frac)*po
     qgdnv(3) = frac*ustar + (one - frac)*uo
  end if

  ! Passive scalars
  do n = 4,nvar
     if(sgnm==one)then
        qgdnv(n) = qleft (n)
     else
        qgdnv(n) = qright(n)
     end if
  end do

  CALL find_mhd_flux(qgdnv,ugdnv,fgdnv)

  call trace_mpi('hydro_acoustic(qleft,qright,fgdnv)','exit',2)
end subroutine hydro_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine eigenvalues(d,vx,vy,vz,p,bx,by,bz,lambda)
!
! MHD adiabatic eigenvalues
!
! Input Arguments:
!   Bx    = magnetic field in sweep direction
!   d     = density
!   vx    = X velocity 
!   vy    = Y velocity 
!   vz    = Z velocity 
!   by    = Y magnetic field
!   bz    = Z magnetic field
!   p     = thermal pressure
!
! Output Arguments:
!
!   lambda  = eigenvalues
! 
  USE amr_parameters
  USE const
  USE hydro_parameters
  IMPLICIT NONE

  real (dp), intent(IN) :: d, vx, vy, vz, p
  real (dp), intent(IN) :: bx, by, bz
  real (dp), dimension(1:7), intent(OUT) :: lambda
  ! local variables
  real (dp) :: btsq, bt, vsq, vax,  vaxsq, hp
  real (dp) :: asq, astarsq, cfsq, cfast, cssq, cslow

  call trace_mpi('eigenvalues(d,vx,vy,vz,p,bx,by,bz,lambda)','enter',3)

  vsq = vx**2+vy**2+vz**2
  btsq = by**2+bz**2
  bt = sqrt(btsq)
  vaxsq = Bx**2/d
  vax = sqrt(vaxsq)
  asq = gamma*p/d
  asq = MAX(asq,smallc**2)
  astarsq = asq+vaxsq+btsq/d
  
  cfsq = .5*(astarsq + sqrt(astarsq**2-4.0*asq*vaxsq))
  cfast = sqrt(cfsq)
  
  cssq = .5*(astarsq - sqrt(astarsq**2-4.0*asq*vaxsq))
  if (cssq .le. 0.) cssq = 0.
  cslow = sqrt(cssq)
  
  lambda(1) = vx - cfast
  lambda(2) = vx - vax
  lambda(3) = vx - cslow
  lambda(4) = vx
  lambda(5) = vx + cslow
  lambda(6) = vx + vax
  lambda(7) = vx + cfast

  call trace_mpi('eigenvalues(d,vx,vy,vz,p,bx,by,bz,lambda)','exit',3)
end subroutine eigenvalues
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE eigen_cons(d,vx,vy,vz,h,Bx,by,bz,Xfac,Yfac,lambda,rem,lem)
!
! Input Arguments:
!   d     = Roe density
!   vx    = Roe X velocity 
!   vy    = Roe Y velocity 
!   vz    = Roe Z velocity 
!   Bx    = magnetic field in sweep direction
!   by    = Roe Y magnetic field
!   bz    = Roe Z magnetic field
!   h     = Roe enthalpy
!   Xfac  = ((by^2-byl*byr)+bz^2-bzl*bzr))/(2*rho)
!   Yfac  = (rho_l+rho_r)/(2*rho)
!
! Output Arguments:
!
!   lambda = eigenvalues
!   rem     = right eigenmatrix
!   lem     = left  eigenmatrix
!
  USE amr_parameters
  USE const
  USE hydro_parameters
  IMPLICIT NONE
  
  REAL (dp) :: d, vx, vy, vz, h
  REAL (dp) :: Bx, by, bz, Xfac, Yfac
  
  REAL(dp), DIMENSION(7,7) :: lem, rem
  REAL(dp), DIMENSION(7)        :: lambda
  
  REAL(dp) :: btsq, bt_starsq, bt, bt_star
  REAL(dp) :: vsq, vax,  vaxsq, hp, twid_asq, q_starsq 
  REAL(dp) :: cfsq, cfast, cssq, cslow
  REAL(dp) :: beta_y, beta_z, beta_ystar, beta_zstar, beta_starsq, vbeta
  REAL(dp) :: alpha_f, alpha_s, droot, s, twid_a
  REAL(dp) :: Qfast, Qslow, af_prime, as_prime, Afpbb, Aspbb, na
  REAL(dp) :: cff, css, af, as, Afpb, Aspb, vqstr, norm
  REAL(dp) :: Q_ystar, Q_zstar
  
  REAL(dp) :: smalle
  
  real(dp),dimension(7)::toto
  integer :: i

  call trace_mpi('eigen_cons(d,vx,vy,vz,h,Bx,by,bz,Xfac,Yfac,lambda,rem,lem)','enter',3)
  
  vsq = vx*vx+vy*vy+vz*vz
  btsq = by*by+bz*bz
  bt_starsq = (gamma-1. - (gamma-2.)*Yfac)*btsq
  bt=sqrt(btsq)
  bt_star = sqrt(bt_starsq)
  
  vaxsq = Bx*Bx/d
  vax = sqrt(vaxsq)
  
  hp = h - (vaxsq + btsq/d)
  twid_asq = ((gamma-1.)*(hp-.5*vsq)-(gamma-2.)*Xfac)
  twid_asq = MAX(twid_asq,smallc*smallc)
  q_starsq = twid_asq+(vaxsq+bt_starsq/d)
  
  cfsq = .5*(q_starsq + sqrt(q_starsq*q_starsq-4.0*twid_asq*vaxsq))
  cfast = sqrt(cfsq)
  
  cssq = .5*(q_starsq - sqrt(q_starsq*q_starsq-4.0*twid_asq*vaxsq))
  if (cssq .le. 0.) cssq = 0.
  cslow = sqrt(cssq)
  
  if (bt .eq. 0) then
     beta_y = .5*sqrt(2.)
     beta_z = .5*sqrt(2.)
     beta_ystar = .5*sqrt(2.)
     beta_zstar = .5*sqrt(2.)
  else
     beta_y = by/bt
     beta_z = bz/bt
     beta_ystar = by/bt_star
     beta_zstar = bz/bt_star
  endif
  beta_starsq = beta_ystar*beta_ystar + beta_zstar*beta_zstar
  vbeta = vy*beta_ystar + vz*beta_zstar
  
  if ( (cfsq - cssq) .eq. 0.) then
     alpha_f = 1.0
     alpha_s = 0.0
  else if ( (twid_asq - cssq) .le. 0.) then
     alpha_f = 0.0
     alpha_s = 1.0
  else if ( (cfsq - twid_asq) .le. 0.) then
     alpha_f = 1.0
     alpha_s = 0.0
  else
     alpha_f = sqrt((twid_asq-cssq)/(cfsq-cssq))
     alpha_s = sqrt((cfsq-twid_asq)/(cfsq-cssq))
  endif
  !
  ! compute Qs and As for eigenmatrices
  !
  droot = sqrt(d)
  s  = SIGN(one,bx)
  twid_a = sqrt(twid_asq)
  Qfast = s*cfast*alpha_f
  Qslow = s*cslow*alpha_s
  af_prime = twid_a*alpha_f/droot
  as_prime = twid_a*alpha_s/droot
  Afpbb = af_prime*bt_star*beta_starsq
  Aspbb = as_prime*bt_star*beta_starsq
  !
  ! eigenvalues
  !
  lambda(1) = vx - cfast
  lambda(2) = vx - vax
  lambda(3) = vx - cslow
  lambda(4) = vx
  lambda(5) = vx + cslow
  lambda(6) = vx + vax
  lambda(7) = vx + cfast
  !
  ! eigenmatrix
  !
  rem(1,1) = alpha_f
  rem(1,2) = alpha_f*(vx-cfast)
  rem(1,3) = alpha_f*vy + Qslow*beta_ystar
  rem(1,4) = alpha_f*vz + Qslow*beta_zstar
  rem(1,5) = alpha_f*(hp-vx*cfast) + Qslow*vbeta + Aspbb
  rem(1,6) = as_prime*beta_ystar
  rem(1,7) = as_prime*beta_zstar

  rem(2,1) = 0.
  rem(2,2) = 0.
  rem(2,3) = -beta_z
  rem(2,4) =  beta_y
  rem(2,5) = -(vy*beta_z - vz*beta_y)
  rem(2,6) = -s*beta_z/droot
  rem(2,7) =  s*beta_y/droot

  rem(3,1) = alpha_s
  rem(3,2) = alpha_s*(vx-cslow)
  rem(3,3) = alpha_s*vy - Qfast*beta_ystar
  rem(3,4) = alpha_s*vz - Qfast*beta_zstar
  rem(3,5) = alpha_s*(hp - vx*cslow) - Qfast*vbeta - Afpbb
  rem(3,6) = -af_prime*beta_ystar
  rem(3,7) = -af_prime*beta_zstar

  rem(4,1) = 1.0
  rem(4,2) = vx
  rem(4,3) = vy
  rem(4,4) = vz
  rem(4,5) = 0.5*vsq + (gamma-2.)*Xfac/(gamma-1.)
  rem(4,6) = 0.
  rem(4,7) = 0.

  rem(5,1) =  alpha_s
  rem(5,2) =  alpha_s*(vx+cslow)
  rem(5,3) =  alpha_s*vy + Qfast*beta_ystar
  rem(5,4) =  alpha_s*vz + Qfast*beta_zstar
  rem(5,5) =  alpha_s*(hp+vx*cslow) + Qfast*vbeta - Afpbb
  rem(5,6) =  rem(3,6)
  rem(5,7) =  rem(3,7)

  rem(6,1) = 0.
  rem(6,2) = 0.
  rem(6,3) =  beta_z
  rem(6,4) = -beta_y
  rem(6,5) = -rem(2,5)
  rem(6,6) =  rem(2,6)
  rem(6,7) =  rem(2,7)

  rem(7,1) =  alpha_f
  rem(7,2) = alpha_f*(vx+cfast)
  rem(7,3) = alpha_f*vy - Qslow*beta_ystar
  rem(7,4) = alpha_f*vz - Qslow*beta_zstar
  rem(7,5) = alpha_f*(hp + vx*cfast) - Qslow*vbeta + Aspbb
  rem(7,6) =  rem(1,6)
  rem(7,7) =  rem(1,7)
  ! 
  ! Left eignematrix
  !
  ! normalize some of the quantities by 1/(2a^2)
  ! some by (gamma-1)/2a^2
  !
  na  = 0.5/twid_asq
  cff = na*alpha_f*cfast
  css = na*alpha_s*cslow
  Qfast = Qfast*na
  Qslow = Qslow*na
  af = na*af_prime*d
  as = na*as_prime*d
  Afpb = na*af_prime*bt_star
  Aspb = na*as_prime*bt_star
  !
  alpha_f = (gamma-1.)*na*alpha_f
  alpha_s = (gamma-1.)*na*alpha_s
  Q_ystar = beta_ystar/beta_starsq
  Q_zstar = beta_zstar/beta_starsq
  vqstr   = (vy*Q_ystar+vz*Q_zstar)
  norm = (gamma-1.)*2.*na

  lem(1,1) = alpha_f*(vsq-hp) + cff*(cfast+vx) - Qslow*vqstr - Aspb
  lem(2,1) = -alpha_f*vx - cff
  lem(3,1) = -alpha_f*vy + Qslow*Q_ystar
  lem(4,1) = -alpha_f*vz + Qslow*Q_zstar
  lem(5,1) =  alpha_f
  lem(6,1) = as*Q_ystar - alpha_f*by
  lem(7,1) = as*Q_zstar - alpha_f*bz

  lem(1,2) =  0.5*(vy*beta_z - vz*beta_y)
  lem(2,2) =  0.
  lem(3,2) = -0.5*beta_z
  lem(4,2) =  0.5*beta_y
  lem(5,2) =  0.
  lem(6,2) = -0.5*droot*beta_z*s
  lem(7,2) =  0.5*droot*beta_y*s

  lem(1,3) =  alpha_s*(vsq-hp) + css*(cslow+vx) + Qfast*vqstr + Afpb
  lem(2,3) = -alpha_s*vx - css
  lem(3,3) = -alpha_s*vy - Qfast*Q_ystar
  lem(4,3) = -alpha_s*vz - Qfast*Q_zstar
  lem(5,3) =  alpha_s
  lem(6,3) = -af*Q_ystar - alpha_s*by
  lem(7,3) = -af*Q_zstar - alpha_s*bz

  ! Attention, il y a une difference de signe avec ramses_mhd...
  lem(1,4) =  1. - norm*(.5*vsq - (gamma-2.)*Xfac/(gamma-1.))
  ! This is the old version...
  !  lem(1,4) =  1. - norm*(.5*vsq + (gamma-2.)*Xfac/(gamma-1.))
  lem(2,4) =  norm*vx
  lem(3,4) =  norm*vy
  lem(4,4) =  norm*vz
  lem(5,4) = -norm
  lem(6,4) =  norm*by
  lem(7,4) =  norm*bz

  lem(1,5) =  alpha_s*(vsq-hp) + css*(cslow-vx) - Qfast*vqstr + Afpb
  lem(2,5) = -alpha_s*vx + css 
  lem(3,5) = -alpha_s*vy + Qfast*Q_ystar 
  lem(4,5) = -alpha_s*vz + Qfast*Q_zstar 
  lem(5,5) =  alpha_s
  lem(6,5) =  lem(6,3)
  lem(7,5) =  lem(7,3)

  lem(1,6) = -lem(1,2)
  lem(2,6) =  0.
  lem(3,6) = -lem(3,2)
  lem(4,6) = -lem(4,2)
  lem(5,6) =  0.
  lem(6,6) =  lem(6,2)
  lem(7,6) =  lem(7,2)

  lem(1,7) =  alpha_f*(vsq-hp) + cff*(cfast-vx) + Qslow*vqstr -Aspb
  lem(2,7) = -alpha_f*vx + cff
  lem(3,7) = -alpha_f*vy - Qslow*Q_ystar
  lem(4,7) = -alpha_f*vz - Qslow*Q_zstar
  lem(5,7) =  alpha_f
  lem(6,7) =  lem(6,1)
  lem(7,7) =  lem(7,1)

  call trace_mpi('eigen_cons(d,vx,vy,vz,h,Bx,by,bz,Xfac,Yfac,lambda,rem,lem)','exit',3)
END SUBROUTINE eigen_cons
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! HLLR3 Flux computation based on mhd_fluxes written
! for FLASH by Knut Waagan
!###########################################################
SUBROUTINE hll3r(qleft,qright,flux,ngrid,speed)
  use amr_parameters,   only : dp, nvector, chemistry
  use amr_commons,      only : myid
  use hydro_parameters, only : nvar, gamma
  implicit none
  integer, intent(in)                          :: ngrid
  REAL(dp),DIMENSION(1:nvector,1:nvar), target :: qleft,qright
  REAL(dp),DIMENSION(1:nvector,1:nvar)         :: flux
  real(dp)                                     :: speed
  REAL(dp),DIMENSION(:), pointer               :: rl,pl,ul,vl,wl,bnl,b1l,b2l, &
                                                  rr,pr,ur,vr,wr,bnr,b1r,b2r
  REAL(dp),DIMENSION(1:nvector)                :: el,er

  integer                             :: i
  real(kind=dp)                       :: alpha 
  real(kind=dp), dimension(nvector)   :: lm, lp, cl, cr, a0l, a0r, cxl2, cxr2, aql, aqr, cfsl2, cfsr2, &
                                         Bx2l, Bx2r, asl2, asr2, Xl, Xr, xsmr, xsml, Bsql, Bsqr, &
                                         ptotl, ptotr, pistar2, rhostarl, rhostarr, u2star, pistar
  real(kind=dp), dimension(nvector,3) :: pil, pir, ustar, pinstar

  call trace_mpi('hll3r(qleft,qright,flux,ngrid,speed)','enter',3)
  
#ifdef KROME
  if (chemistry) then
     if (myid==1) print *, 'HLL3R solver is not supported with variable gamma!'
     stop
  endif
#endif

  ! Left variables
  rl=>qleft(:,1); Pl=>qleft(:,2); ul=>qleft(:,3); Bnl=>qleft(:,4)
  vl=>qleft(:,5); B1l=>qleft(:,6); wl=>qleft(:,7); B2l=>qleft(:,8)

  ! Right variables
  rr=>qright(:,1); Pr=>qright(:,2); ur=>qright(:,3); Bnr=>qright(:,4)
  vr=>qright(:,5); B1r=>qright(:,6); wr=>qright(:,7); B2r=>qright(:,8)

  ! Internal energy :
  if (gamma > 1.0_dp) then
    el = Pl / ((gamma - 1.0_dp)*rl)
    er = Pr / ((gamma - 1.0_dp)*rr)
  else
    el = Pl / rl
    er = Pr / rr
  endif
  
  alpha = 0.5_dp*(gamma+1.0_dp)

  Bx2l=Bnl*Bnl
  Bx2r=Bnr*Bnr

  ! Total pressure and total energy
  Bsql = Bnl*Bnl + B1l*B1l + B2l*B2l
  ptotl=Pl+0.5_dp*Bsql
  
  Bsqr = Bnr*Bnr+B1r*B1r+B2r*B2r
  ptotr=Pr+0.5_dp*Bsqr
  
  pil(:,1) = ptotl - Bnl*Bnl
  pil(:,2) =       - Bnl*B1l
  pil(:,3) =       - Bnl*B2l
  
  pir(:,1) = ptotr - Bnr*Bnr
  pir(:,2) =       - Bnr*B1r
  pir(:,3) =       - Bnr*B2r
  
  ! as[lr] = sound speed
  asl2 = max(gamma*Pl/rl,0.0_dp)
  asr2 = max(gamma*Pr/rr,0.0_dp)
  
  cfsl2=asl2+Bsql/rl 
  cfsr2=asr2+Bsqr/rr
  
  ! aq[lr] = fast magnetosonic speed
  aql  = sqrt(abs(0.5_dp*(cfsl2+sqrt(max(0.0_dp,(cfsl2**2-4.0_dp*asl2*Bx2l/rl))))))
  aqr  = sqrt(abs(0.5_dp*(cfsr2+sqrt(max(0.0_dp,(cfsr2**2-4.0_dp*asr2*Bx2r/rr))))))
  
  Xl = max(ul-ur,0.0_dp)+max((ptotr-Bx2r)-(ptotl-Bx2l),0.0_dp)/(rl*aql+rr*aqr)
  Xr = max(ul-ur,0.0_dp)+max((ptotl-Bx2l)-(ptotr-Bx2r),0.0_dp)/(rl*aql+rr*aqr)
  
  xsml = 1.0_dp-Xl/(aql+alpha*Xl)
  xsmr = 1.0_dp-Xr/(aqr+alpha*Xr)
  cxl2 = asl2*xsml+Bsql/rl
  cxr2 = asr2*xsmr+Bsqr/rr
  
  a0l  = sqrt(abs((0.5_dp/xsml)*(cxl2+sqrt(cxl2**2-4.0_dp*xsml*asl2*Bx2l/rl))))
  a0r  = sqrt(abs((0.5_dp/xsmr)*(cxr2+sqrt(cxr2**2-4.0_dp*xsmr*asr2*Bx2r/rr))))
  
  cl = rl*(a0l + alpha*Xl)
  cr = rr*(a0r + alpha*Xr)
  
  lm = ul-(a0l + alpha*Xl)
  lp = ur+(a0r + alpha*Xr)

  speed = max(speed,max(maxval(abs(lp)),maxval(abs(lm))))

  ustar(:,1) = (cl*ul+cr*ur+pil(:,1)-pir(:,1)) / (cl+cr)
  ustar(:,2) = (cl*vl+cr*vr+pil(:,2)-pir(:,2)) / (cl+cr)
  ustar(:,3) = (cl*wl+cr*wr+pil(:,3)-pir(:,3)) / (cl+cr)
  pinstar(:,1)  = ( cr*pil(:,1) + cl*pir(:,1) - cl*cr*(ur-ul) ) / (cl+cr)
  pinstar(:,2)  = ( cr*pil(:,2) + cl*pir(:,2) - cl*cr*(vr-vl) ) / (cl+cr)
  pinstar(:,3)  = ( cr*pil(:,3) + cl*pir(:,3) - cl*cr*(wr-wl) ) / (cl+cr)
  u2star   = ustar(:,1)*ustar(:,1) + ustar(:,2)*ustar(:,2) + ustar(:,3)*ustar(:,3)
  pistar2  = pinstar(:,1)*pinstar(:,1) + pinstar(:,2)*pinstar(:,2) + pinstar(:,3)*pinstar(:,3)
  rhostarl = 1.0_dp / (1.0_dp / rl + (pil(:,1) - pinstar(:,1)) / (cl*cl))
  rhostarr = 1.0_dp / (1.0_dp / rr + (pir(:,1) - pinstar(:,1)) / (cr*cr))
  
  ! Computation of states and fluxes:
  where (lm > 0.0_dp) ! Super magnetosonic state
    Flux(:,1) = rl*ul
    Flux(:,2) = (rl*(el + 0.5_dp*(ul*ul+vl*vl+wl*wl)) + 0.5_dp*Bsql)*ul &
                  + pil(:,1)*ul + pil(:,2)*vl + pil(:,3)*wl
    
    Flux(:,3) = Flux(:,1)*ul+pil(:,1)
    Flux(:,5) = Flux(:,1)*vl+pil(:,2)
    Flux(:,7) = Flux(:,1)*wl+pil(:,3)
    
    Flux(:,4) = 0.0_dp
    Flux(:,6) = ul*b1l - bnl*vl
    Flux(:,8) = ul*b2l - bnl*wl 
  else where (ustar(:,1) > 0.0_dp) ! inner left state
    Flux(:,1) = rhostarl*ustar(:,1)
    Flux(:,2) = 0.5_dp*rhostarl*(u2star + 2.0_dp*el + Bsql/rl - &
                                 (pil(:,1)*pil(:,1) + pil(:,2)*pil(:,2) + pil(:,3)*pil(:,3) - pistar2) / (cl*cl) ) * &
                  ustar(:,1) + pinstar(:,1)*ustar(:,1) + pinstar(:,2)*ustar(:,2) + pinstar(:,3)*ustar(:,3)
    
    Flux(:,3) = Flux(:,1)*ustar(:,1) + pinstar(:,1)
    Flux(:,5) = Flux(:,1)*ustar(:,2) + pinstar(:,2)
    Flux(:,7) = Flux(:,1)*ustar(:,3) + pinstar(:,3)

    Flux(:,4) = 0.0_dp
    Flux(:,6) = ustar(:,1)* rhostarl*(B1l / rl + Bnl*(pil(:,2)-pinstar(:,2))/(cl*cl)) - Bnl*ustar(:,2)
    Flux(:,8) = ustar(:,1)* rhostarl*(B2l / rl + Bnl*(pil(:,3)-pinstar(:,3))/(cl*cl)) - Bnl*ustar(:,3)
  else where (lp > 0.0_dp) ! inner right state
    Flux(:,1) = rhostarr*ustar(:,1)
    Flux(:,2) = 0.5_dp*rhostarr*(u2star + 2.0_dp*er + Bsqr/rr -  &
                                 (pir(:,1)*pir(:,1) + pir(:,2)*pir(:,2) + pir(:,3)*pir(:,3) - pistar2) / (cr*cr) ) * &
                  ustar(:,1) + pinstar(:,1)*ustar(:,1) + pinstar(:,2)*ustar(:,2) + pinstar(:,3)*ustar(:,3)
    
    Flux(:,3) = Flux(:,1)*ustar(:,1) + pinstar(:,1)
    Flux(:,5) = Flux(:,1)*ustar(:,2) + pinstar(:,2)
    Flux(:,7) = Flux(:,1)*ustar(:,3) + pinstar(:,3)

    Flux(:,4) = 0.0_dp
    Flux(:,6) = ustar(:,1)*rhostarr*(B1r / rr + Bnr*(pir(:,2)-pinstar(:,2))/(cr*cr)) - Bnr*ustar(:,2)
    Flux(:,8) = ustar(:,1)*rhostarr*(B2r / rr + Bnr*(pir(:,3)-pinstar(:,3))/(cr*cr)) - Bnr*ustar(:,3)            
  else where ! Super magnetosonic state
    Flux(:,1) = rr*ur
    Flux(:,2) = (rr*(er + 0.5_dp*(ur**2+vr**2+wr**2)) + 0.5_dp*Bsqr)*ur &
                  + pir(:,1)*ur + pir(:,2)*vr + pir(:,3)*wr
    
    Flux(:,3) = Flux(:,1)*ur+pir(:,1)
    Flux(:,5) = Flux(:,1)*vr+pir(:,2)
    Flux(:,7) = Flux(:,1)*wr+pir(:,3)
    
    Flux(:,4) = 0.0_dp
    Flux(:,6) = ur*b1r - bnr*vr
    Flux(:,8) = ur*b2r - bnr*wr
  endwhere

  call trace_mpi('hll3r(qleft,qright,flux,ngrid,speed)','exit',3)
END SUBROUTINE hll3r
! HLLR3 Wave speeds
!###########################################################
SUBROUTINE hll3r_speeds(qleft,qright,lm,lp)
  use amr_parameters,   only : dp, chemistry
  use amr_commons,      only : myid
  use hydro_parameters, only : nvar, gamma
  implicit none
  real(dp),dimension(1:nvar) :: qleft,qright
  real(dp)                   :: lm, lp
  real(dp)                   :: rl,pl,ul,bnl,b1l,b2l, &
                                rr,pr,ur,bnr,b1r,b2r
  real(kind=dp)               :: cl, cr, a0l, a0r, cxl2, cxr2, aql, aqr
  real(kind=dp)               :: cfsl2, cfsr2, pixl, pixr, Bx2l, Bx2r, asl2, asr2, Xl, Xr
  real(kind=dp)               :: xsmr, xsml, Bsql, Bsqr, alpha, covrhol, covrhor

#ifdef KROME
  if (chemistry) then
     if (myid==1) print *, 'HLL3R solver is not supported with variable gamma!'
     stop
  endif
#endif

  call trace_mpi('hll3r_speeds(qleft,qright,lm,lp)','enter',3)

  ! Left and Right variables
  rl=qleft(1);  Pl=qleft(2);  ul=qleft(3);  Bnl=qleft(4);  B1l=qleft(6);  B2l=qleft(8)
  rr=qright(1); Pr=qright(2); ur=qright(3); Bnr=qright(4); B1r=qright(6); B2r=qright(8)
  
  alpha = 0.5_dp*(gamma+1.0_dp)

  Bx2l=Bnl*Bnl
  Bx2r=Bnr*Bnr

  ! Total pressure and total energy
  Bsql = Bnl*Bnl+B1l*B1l+B2l*B2l
  pixl= Pl+0.5_dp*Bsql-Bx2l
  
  Bsqr = Bnr*Bnr+B1r*B1r+B2r*B2r
  pixr= Pr+0.5_dp*Bsqr-Bx2r
  
  ! as[lr] = sound speed
  asl2 = max(gamma*Pl/rl,0.0_dp)
  asr2 = max(gamma*Pr/rr,0.0_dp)
  
  cfsl2=asl2+Bsql/rl 
  cfsr2=asr2+Bsqr/rr
  
  ! aq[lr] = fast magnetosonic speed
  aql = sqrt(abs(0.5_dp*(cfsl2+sqrt(max(0.0_dp,(cfsl2**2-4.0_dp*asl2*Bx2l/rl))))))
  aqr = sqrt(abs(0.5_dp*(cfsr2+sqrt(max(0.0_dp,(cfsr2**2-4.0_dp*asr2*Bx2r/rr))))))
  
  Xl = max(ul-ur,0.0_dp)+max(pixr-pixl,0.0_dp)/(rl*aql+rr*aqr)
  Xr = max(ul-ur,0.0_dp)+max(pixl-pixr,0.0_dp)/(rl*aql+rr*aqr)
  
  xsml = 1.0_dp-Xl/(aql+alpha*Xl)
  xsmr = 1.0_dp-Xr/(aqr+alpha*Xr)
  cxl2 = asl2*xsml+Bsql/rl
  cxr2 = asr2*xsmr+Bsqr/rr
  
  a0l = sqrt(abs((0.5_dp/xsml)*(cxl2+sqrt(cxl2**2-4.0_dp*xsml*asl2*Bx2l/rl))))
  a0r = sqrt(abs((0.5_dp/xsmr)*(cxr2+sqrt(cxr2**2-4.0_dp*xsmr*asr2*Bx2r/rr))))
  
  covrhol = a0l + alpha*Xl
  covrhor = a0r + alpha*Xr
  
  lm = ul-covrhol
  lp = ur+covrhor
  call trace_mpi('hll3r_speeds(qleft,qright,lm,lp)','exit',3)
END SUBROUTINE hll3r_speeds
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! HLLR5 Flux computation based on mhd_fluxes written
! for FLASH by Knut Waagan
!###########################################################
SUBROUTINE hll5r(qleft,qright,flux,ngrid,speed)
  use amr_parameters,   only : dp, nvector, chemistry
  use amr_commons,      only : myid
  use hydro_parameters, only : nvar, gamma
  implicit none
  integer, intent(in)                          :: ngrid
  REAL(dp),DIMENSION(1:nvector,1:nvar), target :: qleft,qright
  REAL(dp),DIMENSION(1:nvector,1:nvar)         :: flux
  real(dp)                                     :: speed
  REAL(dp),DIMENSION(:), pointer               :: rl,pl,ul,vl,wl,bnl,b1l,b2l, &
                                                  rr,pr,ur,vr,wr,bnr,b1r,b2r
  REAL(dp),DIMENSION(1:nvector)                :: el,er

  integer                     :: i

  real(dp) :: lm,lp,lam,lap,lbm,lbp, &
              cal, car, cbl, cbr, cbovrhol, cbovrhor, caovrhol, caovrhor, aql, aqr, aql2, aqr2, &
              pixl, pixr, Bx2l, Bx2r, asl2, asr2, asl, asr, BxBtl,BxBtr, Xl, Xr, &
              xsmr, xsml, Bsql, Bsqr, alpha, Etotl, Etotr, ptotl, ptotr, &
              pi2l, pi2r,pistar2, rhostarl, rhostarr, eintstarl, eintstarr, &
              B2starl, B2starr, u2star, Etotstarl, Etotstarr
  real(dp), dimension(3) :: pil, pir, ustar,  Bstarl, Bstarr, pinstar, ctabl, ctabr, ustarmid

  call trace_mpi('hll5r(qleft,qright,flux,ngrid,speed)','enter',3)

#ifdef KROME
  if (chemistry) then
     if (myid==1) print *, 'HLL5R solver is not supported with variable gamma!'
     stop
  endif
#endif

  ! Left variables
  rl=>qleft(:,1); Pl=>qleft(:,2); ul=>qleft(:,3); Bnl=>qleft(:,4)
  vl=>qleft(:,5); B1l=>qleft(:,6); wl=>qleft(:,7); B2l=>qleft(:,8)

  ! Right variables
  rr=>qright(:,1); Pr=>qright(:,2); ur=>qright(:,3); Bnr=>qright(:,4)
  vr=>qright(:,5); B1r=>qright(:,6); wr=>qright(:,7); B2r=>qright(:,8)

  ! Internal energy :
  if (gamma > 1.0_dp) then
    el = pl / ((gamma - 1.0_dp)*rl)
    er = pr / ((gamma - 1.0_dp)*rr)
  else
    el = pl / rl
    er = pr / rr
  endif
  
  alpha = 0.5_dp*(gamma+1.0_dp)

  do i = 1, ngrid


    ! Total pressure and total energy
    Bx2l=Bnl(i)*Bnl(i)
    Bsql = Bnl(i)*Bnl(i)+B1l(i)*B1l(i)+B2l(i)*B2l(i)
    ptotl=Pl(i)+0.5_dp*Bsql
    pixl= ptotl-Bx2l
    Etotl = rl(i)*(el(i) + 0.5_dp*(ul(i)**2+vl(i)**2+wl(i)**2)) + 0.5_dp*Bsql
    
    Bx2r=Bnr(i)*Bnr(i)
    Bsqr = Bnr(i)*Bnr(i)+B1r(i)*B1r(i)+B2r(i)*B2r(i)
    ptotr=Pr(i)+0.5_dp*Bsqr
    pixr= ptotr-Bx2r
    Etotr = rr(i)*(er(i) + 0.5_dp*(ur(i)**2+vr(i)**2+wr(i)**2)) + 0.5_dp*Bsqr
    
    pil    =-Bnl(i)*(/Bnl(i),B1l(i),B2l(i)/)
    pil(1) = pil(1) + ptotl
    pi2l   = pil(1)*pil(1) + pil(2)*pil(2) + pil(3)*pil(3)
    
    pir    =-Bnr(i)*(/Bnr(i),B1r(i),B2r(i)/)
    pir(1) = pir(1) + ptotr
    pi2r   = pir(1)*pir(1) + pir(2)*pir(2) + pir(3)*pir(3)

    ! as[lr] = sound speed
    asl2 = max(gamma*Pl(i)/rl(i),0.0_dp)
    asr2 = max(gamma*Pr(i)/rr(i),0.0_dp)
    asl  = sqrt(asl2)
    asr  = sqrt(asr2)

    BxBtl = abs(bnl(i))*sqrt(max(0.0_dp,Bsql-Bx2l))
    BxBtr = abs(bnr(i))*sqrt(max(0.0_dp,Bsqr-Bx2r))

    ! aq[lr] = fast magnetosonic speed
    aql2 = abs(asl2 + (Bsql + BxBtl)/rl(i))
    aqr2 = abs(asr2 + (Bsqr + BxBtr)/rr(i))
    
    aql=sqrt(aql2)
    aqr=sqrt(aqr2)
    
    Xl = max(ul(i)-ur(i),0.0_dp)+max(pixr-pixl,0.0_dp)/(rl(i)*aql+rr(i)*aqr)
    Xr = max(ul(i)-ur(i),0.0_dp)+max(pixl-pixr,0.0_dp)/(rl(i)*aql+rr(i)*aqr)

    xsml = 1.0_dp-Xl/(aql+alpha*Xl)
    xsmr = 1.0_dp-Xr/(aqr+alpha*Xr)

    caovrhol = sqrt((Bx2l + BxBtl)/rl(i)/xsml)
    caovrhor = sqrt((Bx2r + BxBtr)/rr(i)/xsmr)
    cbovrhol = aql/sqrt(xsml) + alpha*Xl 
    cbovrhor = aqr/sqrt(xsmr) + alpha*Xr 

    ! Avoid ill-posedness. 
    caovrhol = max (caovrhol, 1.0e-12_dp*cbovrhol)
    caovrhor = max (caovrhor, 1.0e-12_dp*cbovrhor)

    cal=rl(i)*caovrhol
    car=rr(i)*caovrhor
    cbl=rl(i)*cbovrhol
    cbr=rr(i)*cbovrhor

    ! Velocities
    ustar(1)    = (cbl*ul(i) + cbr*ur(i) + pixl   - pixr)   / (cbl+cbr)
    ustarmid(1) = ustar(1)
    ustarmid(2) = (cal*vl(i) + car*vr(i) + pil(2) - pir(2)) / (cal+car)
    ustarmid(3) = (cal*wl(i) + car*wr(i) + pil(3) - pir(3)) / (cal+car)

    lbm = ul(i) - cbovrhol
    lbp = ur(i) + cbovrhor

    lam = ul(i) - cal/rl(i)
    lap = ur(i) + car/rr(i)

    lm = min(lam,lbm)
    lp = max(lap, lbp)

    speed = max(speed,max(abs(lp),abs(lm)))

    ! Computation of states and fluxes:
    if (lm >= 0.0_dp) then
      Flux(i,1) = rl(i)*ul(i)
      Flux(i,2) = Etotl*ul(i) + pil(1)*ul(i) + pil(2)*vl(i) + pil(3)*wl(i)
      
      Flux(i,3) = Flux(i,1)*ul(i) + pil(1)
      Flux(i,5) = Flux(i,1)*vl(i) + pil(2)
      Flux(i,7) = Flux(i,1)*wl(i) + pil(3)
      
      Flux(i,4) = 0.0_dp
      Flux(i,6) = ul(i)*b1l(i) - bnl(i)*vl(i)
      Flux(i,8) = ul(i)*b2l(i) - bnl(i)*wl(i) 
    else if (ustar(1) > 0.0_dp) then
      if (lbm > 0.0_dp ) then ! lam < 0 
        ustar(1)  = ul(i)
        ustar(2)  = ustarmid(2)
        ustar(3)  = ustarmid(3)
        pinstar(1)= pixl
        pinstar(2)= ( car*pil(2) + cal*pir(2) - cal*car*(vr(i)-vl(i)) ) / (cal+car)
        pinstar(3)= ( car*pil(3) + cal*pir(3) - cal*car*(wr(i)-wl(i)) ) / (cal+car)
        rhostarl  = rl(i)
        Bstarl(1) = Bnl(i)
        Bstarl(2) = B1l(i) + rhostarl*Bnl(i)*(pil(2)-pinstar(2))/(cal*cal)
        Bstarl(3) = B2l(i) + rhostarl*Bnl(i)*(pil(3)-pinstar(3))/(cal*cal)
        u2star    = ustar(1)*ustar(1) + ustar(2)*ustar(2) + ustar(3)*ustar(3)
        B2starl   = Bstarl(1)*Bstarl(1) + Bstarl(2)*Bstarl(2) + Bstarl(3)*Bstarl(3)
        pistar2   = pinstar(1)*pinstar(1) + pinstar(2)*pinstar(2) + pinstar(3)*pinstar(3)
        eintstarl = el(i) + 0.5_dp*((Bsql/rl(i) - B2starl/rhostarl) - (pi2l - pistar2) / (cal*cal))
      else if (lam < 0.0_dp) then
        ctabl    = (/cbl, cal, cal/)
        ctabr    = (/cbr, car, car/)
        ustar    = (ctabl*(/ul(i),vl(i),wl(i)/) + ctabr*(/ur(i),vr(i),wr(i)/) + pil - pir) / (ctabl+ctabr)
        pinstar  = (ctabr*pil + ctabl*pir - ctabl*ctabr*(/ur(i)-ul(i),vr(i)-vl(i),wr(i)-wl(i)/)) / (ctabl+ctabr)
        rhostarl = 1.0_dp / (1.0_dp / rl(i) + (pil(1) - pinstar(1)) / (cbl*cbl))
        Bstarl   = rhostarl * ( (/Bnl(i),B1l(i),B2l(i)/) / rl(i) + Bnl(i)*(pil-pinstar)/(ctabl*ctabl) )
        u2star   = ustar(1)*ustar(1) + ustar(2)*ustar(2) + ustar(3)*ustar(3)
        B2starl  = Bstarl(1)*Bstarl(1) + Bstarl(2)*Bstarl(2) + Bstarl(3)*Bstarl(3)
        pi2l     = pil(1)*pil(1)/(cbl*cbl) + (pil(2)*pil(2) + pil(3)*pil(3))/(cal*cal)
        pistar2  = pinstar(1)*pinstar(1)/(cbl*cbl) + (pinstar(2)*pinstar(2) + pinstar(3)*pinstar(3))/(cal*cal) 
        eintstarl= el(i) + 0.5_dp*((Bsql/rl(i) - B2starl/rhostarl) - (pi2l - pistar2))
      else ! lam > 0, lbm < 0
        ustar(2)  = vl(i)
        ustar(3)  = wl(i)
        pinstar(1)= (cbr*pixl + cbl*pixr - cbl*cbr*(ur(i) - ul(i))) / (cbl + cbr)
        pinstar(2)= pil(2)
        pinstar(3)= pil(3)
        rhostarl  = 1.0_dp / (1.0_dp / rl(i) + (pil(1) - pinstar(1)) / (cbl*cbl))
        Bstarl(1) = Bnl(i)
        Bstarl(2) = rhostarl * B1l(i) / rl(i)
        Bstarl(3) = rhostarl * B2l(i) / rl(i)
        u2star    = ustar(1)*ustar(1) + ustar(2)*ustar(2) + ustar(3)*ustar(3)
        B2starl   = Bstarl(1)*Bstarl(1) + Bstarl(2)*Bstarl(2) + Bstarl(3)*Bstarl(3)
        eintstarl = el(i) + 0.5_dp*((Bsql/rl(i) - B2starl/rhostarl) - (pixl*pixl - pinstar(1)*pinstar(1))/(cbl*cbl))
      endif ! End of left inner states
      ! Calculate flux for left inner states
      Flux(i,1) = rhostarl*ustar(1)
      Etotstarl = rhostarl*(0.5_dp*u2star + eintstarl) + 0.5_dp*B2starl
      Flux(i,2) = Etotstarl*ustar(1) + pinstar(1)*ustar(1) + pinstar(2)*ustar(2) + pinstar(3)*ustar(3)
    
      Flux(i,3) = Flux(i,1)*ustar(1) + pinstar(1)
      Flux(i,5) = Flux(i,1)*ustar(2) + pinstar(2)
      Flux(i,7) = Flux(i,1)*ustar(3) + pinstar(3)

      Flux(i,4) = 0.0_dp
      Flux(i,6) = ustar(1)*Bstarl(2) - Bnl(i)*ustar(2)
      Flux(i,8) = ustar(1)*Bstarl(3) - Bnl(i)*ustar(3)
    else if (lp > 0.0_dp) then
      if (lbp < 0.0_dp ) then ! lap > 0
        ustar(1)  = ur(i)
        ustar(2)  = ustarmid(2)
        ustar(3)  = ustarmid(3)
        pinstar(1)= pixr
        pinstar(2)= (car*pil(2) + cal*pir(2) - cal*car*(vr(i)-vl(i))) / (cal + car)
        pinstar(3)= (car*pil(3) + cal*pir(3) - cal*car*(wr(i)-wl(i))) / (cal + car)  
        rhostarr  = rr(i)
        Bstarr(1) = Bnr(i)
        Bstarr(2) = B1r(i) + rhostarr*Bnr(i)*(pir(2)-pinstar(2))/(car*car)
        Bstarr(3) = B2r(i) + rhostarr*Bnr(i)*(pir(3)-pinstar(3))/(car*car)
        u2star    = ustar(1)*ustar(1) + ustar(2)*ustar(2) + ustar(3)*ustar(3)
        B2starr   = Bstarr(1)*Bstarr(1) + Bstarr(2)*Bstarr(2) + Bstarr(3)*Bstarr(3)
        pistar2   = pinstar(1)*pinstar(1) + pinstar(2)*pinstar(2) + pinstar(3)*pinstar(3)
        eintstarr = er(i) + 0.5_dp*((Bsqr/rr(i) - B2starr/rhostarr) - (pi2r - pistar2) / (car*car))
      else if (lap > 0.0_dp) then ! lbp > 0, lap > 0
        ctabl     = (/cbl, cal, cal/)
        ctabr     = (/cbr, car, car/)
        ustar     = (ctabl*(/ul(i),vl(i),wl(i)/) + ctabr*(/ur(i),vr(i),wr(i)/) + pil - pir) / (ctabl+ctabr)
        pinstar   = (ctabr*pil + ctabl*pir - ctabl*ctabr*(/ur(i)-ul(i),vr(i)-vl(i),wr(i)-wl(i)/)) / (ctabl+ctabr)
        rhostarr  = 1.0_dp / (1.0_dp / rr(i) + (pir(1) - pinstar(1)) / (cbr*cbr))
        Bstarr    = rhostarr * ( (/Bnr(i),B1r(i),B2r(i)/) / rr(i) + Bnr(i)*(pir-pinstar)/(ctabr*ctabr) )
        u2star    = ustar(1)*ustar(1) + ustar(2)*ustar(2) + ustar(3)*ustar(3)
        B2starr   = Bstarr(1)*Bstarr(1) + Bstarr(2)*Bstarr(2) + Bstarr(3)*Bstarr(3)
        pi2r      = pir(1)*pir(1)/(cbr*cbr) + (pir(2)*pir(2) + pir(3)*pir(3))/(car*car)
        pistar2   = pinstar(1)*pinstar(1)/(cbr*cbr) + (pinstar(2)*pinstar(2) + pinstar(3)*pinstar(3))/(car*car) 
        eintstarr = er(i) + 0.5_dp*((Bsqr/rr(i) - B2starr/rhostarr) - (pi2r - pistar2))
      else ! lbp >= 0, lap <= 0
        ustar(2)  = vr(i)
        ustar(3)  = wr(i)
        pinstar(1)= (cbr*pixl + cbl*pixr - cbl*cbr*(ur(i) - ul(i))) / (cbl + cbr)
        pinstar(2)= pir(2)
        pinstar(3)= pir(3)
        rhostarr  = 1.0_dp / (1.0_dp / rr(i) + (pir(1) - pinstar(1)) / (cbr*cbr))
        Bstarr(1) = Bnr(i)
        Bstarr(2) = rhostarr * B1r(i) / rr(i)
        Bstarr(3) = rhostarr * B2r(i) / rr(i)
        u2star    = ustar(1)*ustar(1) + ustar(2)*ustar(2) + ustar(3)*ustar(3)
        B2starr   = Bstarr(1)*Bstarr(1) + Bstarr(2)*Bstarr(2) + Bstarr(3)*Bstarr(3)
        eintstarr = er(i) + 0.5_dp*((Bsqr/rr(i) - B2starr/rhostarr) - (pixr*pixr - pinstar(1)*pinstar(1))/(cbr*cbr))
      endif! End of right inner states 
      ! Calculate flux for right inner states
      Flux(i,1) = rhostarr*ustar(1)
      Etotstarr = rhostarr*(0.5_dp*u2star + eintstarr) + 0.5_dp*B2starr
      Flux(i,2) = Etotstarr*ustar(1) + pinstar(1)*ustar(1) + pinstar(2)*ustar(2) + pinstar(3)*ustar(3)
    
      Flux(i,3) = Flux(i,1)*ustar(1) + pinstar(1)
      Flux(i,5) = Flux(i,1)*ustar(2) + pinstar(2)
      Flux(i,7) = Flux(i,1)*ustar(3) + pinstar(3)

      Flux(i,4) = 0.0_dp
      Flux(i,6) = ustar(1)*Bstarr(2) - Bnr(i)*ustar(2)
      Flux(i,8) = ustar(1)*Bstarr(3) - Bnr(i)*ustar(3)
    else ! Leftwards supermagnetosonic
      Flux(i,1) = rr(i)*ur(i)
      Flux(i,2) = Etotr*ur(i) + pir(1)*ur(i) + pir(2)*vr(i) + pir(3)*wr(i)
      
      Flux(i,3)=Flux(i,1)*ur(i) + pir(1)
      Flux(i,5)=Flux(i,1)*vr(i) + pir(2)
      Flux(i,7)=Flux(i,1)*wr(i) + pir(3)
      
      Flux(i,4)=0.0_dp
      Flux(i,6)=ur(i)*b1r(i) - bnr(i)*vr(i)
      Flux(i,8)=ur(i)*b2r(i) - bnr(i)*wr(i) 
    endif
  end do
  call trace_mpi('hll5r(qleft,qright,flux,ngrid,speed)','exit',3)
END SUBROUTINE hll5r
!###########################################################
!###########################################################
!###########################################################
!###########################################################
