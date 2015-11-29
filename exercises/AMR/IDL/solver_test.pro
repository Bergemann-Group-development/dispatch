PRO update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, dt, gamma
  nu = [0.2,0.6]
  ; downstaggered velocity field
  lnr = alog(rho)
  ldx = xdn(lnr) & ldy = ydn(lnr) & ldz = zdn(lnr)
  rdx = exp(ldx) & rdy = exp(ldy) & rdz = exp(ldz)
  Ux = px / rdx  & Uy = py / rdy  & Uz = pz / rdz
  if (gamma eq 1.) then begin
    Pg = rho
    cs = 1. + rho*0.
  endif else begin
    Pg = rho^gamma*exp(S/rho*(gamma-1.))
    cs = sqrt(gamma*Pg/rho)
  endelse
  ; diffusion and artificial pressure
  ds = min([dx,dy,dz])
  du = nu[0]*cs*ds                                          ; dimension: viscosity
  pa = nu[1]*rho*(ds*(ddxup(Ux)+ddyup(Uy)+ddzup(Uz)<0))^2   ; dimension: pressure
  ; Conservation of mass
  fdifx = rdx*ddxdn(lnr)*xdn(du) ; diffusive mass flux
  fdify = rdy*ddydn(lnr)*ydn(du)
  fdifz = rdz*ddzdn(lnr)*zdn(du)
  ;
  drhodt = - ddxup(px - fdifx) - ddyup(py - fdify) - ddzup(pz - fdifz)
  ; Conservation of entropy
  if (gamma eq 1.0) then begin
    dsdt = 0.*s
  endif else begin
    ss = s / rho
    dsdt = $
      - ddxup1(xdn(ss)*px - rdx*ddxdn(ss)*xdn1(du) - fdifx*xdn(ss)) $
      - ddyup1(ydn(ss)*py - rdy*ddydn(ss)*ydn1(du) - fdify*ydn(ss)) $
      - ddzup1(zdn(ss)*pz - rdz*ddzdn(ss)*zdn1(du) - fdifz*zdn(ss))
  endelse
  ; conservation of momentum. f[xyz][xyz] are the fluxes

  ; diagonal terms - Pm==magnetic pressure
  Pm = 0.5*(xup(Bx)^2 + yup(By)^2 + zup(Bz)^2)
  fxx = rho*xup(Ux)^2 - xup(Bx)^2 + Pg + Pm
  fyy = rho*yup(Uy)^2 - yup(By)^2 + Pg + Pm
  fzz = rho*zup(Uz)^2 - zup(Bz)^2 + Pg + Pm

  ; diagnoal part of viscous stress
  fxx += pa - rho*du*ddxup(Ux) ; rho*U*ds*U/ds=rho*U^2
  fyy += pa - rho*du*ddyup(Uy)
  fzz += pa - rho*du*ddzup(Uz)

  ; off-diagnoal part of Reynolds stress, without the density factor
  fxy = ydn(Ux)*xdn(Uy)
  fyz = zdn(Uy)*ydn(Uz)
  fzx = xdn(Uz)*zdn(Ux)

  ; off-diagonal viscous part
  fxy += - xdn1(ydn1(du))*(ddxdn(Uy)+ddydn(Ux))
  fyz += - ydn1(zdn1(du))*(ddydn(Uz)+ddzdn(Uy))
  fzx += - zdn1(xdn1(du))*(ddzdn(Ux)+ddxdn(Uz))

  ; multiply with edge centred density
  fxy = exp(ydn(ldx))*fxy
  fyz = exp(zdn(ldy))*fyz
  fzx = exp(xdn(ldz))*fzx

  ; off-diagonal part of Maxwell stress
  fxy += - ydn(Bx)*xdn(By)
  fyz += - zdn(By)*ydn(Bz)
  fzx += - xdn(Bz)*zdn(Bx)

  ; consistent addition of mass-flux related momentum-flux
  fxx += - (rho*du)*ddxup(ldx)*xup(Ux)
  fyy += - (rho*du)*ddyup(ldy)*yup(Uy)
  fzz += - (rho*du)*ddzup(ldz)*zup(Uz)
  fxy += - (xdn1(fdify)*ydn(Ux)+ydn1(fdifx)*xdn(Uy))
  fyz += - (ydn1(fdifz)*zdn(Uy)+zdn1(fdify)*ydn(Uz))
  fzx += - (zdn1(fdifx)*xdn(Uz)+xdn1(fdifz)*ydn(Ux))

  dpxdt = - ddxdn1(fxx) - ddyup1(fxy) - ddzup1(fzx)
  dpydt = - ddydn1(fyy) - ddzup1(fyz) - ddxup1(fxy)
  dpzdt = - ddzdn1(fzz) - ddxup1(fzx) - ddyup1(fyz)

  ; electric current and resistive electric field
  Jx = ddydn(Bz) - ddzdn(By)
  Jy = ddzdn(Bx) - ddxdn(Bz)
  Jz = ddxdn(By) - ddydn(Bx)
  Ex = ydn1(zdn1(du))*Jx
  Ey = zdn1(xdn1(du))*Jy
  Ez = xdn1(ydn1(du))*Jz
  ; Joule heating
  QJ = yup1(zup1(Ex*Jx)) $
     + zup1(xup1(Ey*Jy)) $
     + xup1(yup1(Ez*Jz))
  dsdt += QJ*rho
  ; advective E-field
  Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)
  Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)
  Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)

  dBxdt = ddzup(Ey) - ddyup(Ez)
  dBydt = ddxup(Ez) - ddzup(Ex)
  dBzdt = ddyup(Ex) - ddxup(Ey)

  ; Time stepping

  rho += dt * drhodt
  s   += dt * dsdt
  px  += dt * dpxdt
  py  += dt * dpydt
  pz  += dt * dpzdt
  Bx  += dt * dBxdt
  By  += dt * dBydt
  Bz  += dt * dBzdt
END

stagger_2nd ; get stagger operators
;stagger_6th ; get stagger operators
n = 128     ; grid size
s = 1.      ; box size
dx = s/n & dy=1. & dz=1.
init_derivs, dx=dx, dy=dy, dz=dz

; Initial conditions
gamma = 5./3.
; left state
l = lonarr(n/2)
rhol = 1.0d
ul   = [1.0d,0.0d,0.0d]
Bl   = [0.0d,1.0d,0.0d]
Pl   = 1.0d
; right state
r = lonarr(n-n/2)
rhor=1.0d
ur   = [-1.0d,0.0d,0.0d]
Br   = [0.0d,0.0d,1.0d]
Pr   = 1.0d

; set up variables
rho = [l+rhol,r+rhor]

px  = rho * [l+ul[0],r+ur[0]]
py  = rho * [l+ul[1],r+ur[1]]
pz  = rho * [l+ul[2],r+ur[2]]
px = xdn1(px)

Bx  = [l+Bl[0],r+Br[0]]
By  = [l+Bl[1],r+Br[1]]
Bz  = [l+Bl[2],r+Br[2]]

P = [l+Pl,r+Pr]

S = rho / (gamma-1.)*alog(P/rho^gamma)

; Run benchmark
dt = 0.05*dx
t = 0.0
tend = 1.

ws,0
!p.multi = [0,2,2]

time_evol=2
every=10
it=0
while (t lt tend) do begin
  if (time_evol eq 1) then begin ; Euler step
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, 0.5*dt, gamma
  endif else if (time_Evol eq 2) then begin ; mid-point method
    drho=rho & dpx=px & dpy=py & dpz=pz & dBx=Bx & dBy=By & dBz=Bz & dS=S
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, 0.5*dt, gamma
    drho-=rho & dpx-=px & dpy-=py & dpz-=pz & dBx-=Bx & dBy-=By & dBz-=Bz & dS-=S
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, dt, gamma
    rho+=drho & px+=dpx & py+=dpy & pz+=dpz & Bx+=dBx & By+=dBy & Bz+=dBz & S+=dS
  endif else if (time_Evol eq 3) then begin ; Heun's method
    drho=rho & dpx=px & dpy=py & dpz=pz & dBx=Bx & dBy=By & dBz=Bz & dS=S
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, dt, gamma
    drho-=rho & dpx-=px & dpy-=py & dpz-=pz & dBx-=Bx & dBy-=By & dBz-=Bz & dS-=S
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, 0.5*dt, gamma
    rho+=.5*drho & px+=.5*dpx & py+=.5*dpy & pz+=.5*dpz & Bx+=.5*dBx & By+=.5*dBy & Bz+=.5*dBz & S+=.5*dS
  endif

  t = t + dt
  it=it+1

  !x.style=1
  if it mod 10 eq 1 then begin
	  plot, rho, ps=-1,tit='rho time='+str(t,format='(f5.3)'), symsize=0.3,yst=3
	  plot, px, ps=-1, tit='px', symsize=0.3,yst=3
	;  plot, py, ps=-1, tit='py', symsize=0.3
	  ux = px/exp(xdn(alog(rho)))
	  plot, ux, ps=-1, tit='ux', symsize=0.3
	  plot, S, ps=-1, tit='s', symsize=0.3,yst=3
	  wait, 0.1
  end
  ;stat, rho, l='rho'
  ;stat, px , l='Px'
  ;stat, S  , l='S'
  ;stat, P  , l='P'
endwhile
ws,0
!p.multi=[0,2,2]
!x.style=3
!y.style=3
plot, rho, ps=-1,tit='rho time='+str(t,format='(f5.3)'), symsize=0.3
plot, px, ps=-1, tit='px', symsize=0.3
plot, py, ps=-1, tit='py', symsize=0.3
plot, S, ps=-1, tit='s', symsize=0.3,yst=3
ws,1
;!p.multi=[0,2,2]
plot, rho, xr=5+[0,20], psym=-1,yst=3,xst=3
plot, px, xr=5+[0,20], psym=-1,yst=3,xst=3

END
