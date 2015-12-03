@stagger_6th
@stagger_2nd

PRO update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, gamma, vf, drhodt, dpxdt, dpydt, dpzdt, dBxdt, dBydt, dBzdt, dsdt
  nu = [0.5,1.0]
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
  pa = nu[1]*(ds*(ddxup(Ux)+ddyup(Uy)+ddzup(Uz)<0))^2       ; dimension: pressure/density
  ; Conservation of mass
  fdifx = rdx*ddxdn(lnr)*xdn(du) ; diffusive mass flux
  fdify = rdy*ddydn(lnr)*ydn(du)
  fdifz = rdz*ddzdn(lnr)*zdn(du)
  ;
  drhodt = - ddxup(px - fdifx) - ddyup(py - fdify) - ddzup(pz - fdifz)

  ; conservation of momentum. f[xyz][xyz] are the fluxes
  ; diagnoal part of viscous stress
  Sxx = ddxup(Ux)
  Syy = ddyup(Uy)
  Szz = ddzup(Uz)
  fxx = rho*(pa - du*Sxx) ; rho*U*ds*U/ds=rho*U^2
  fyy = rho*(pa - du*Syy)
  fzz = rho*(pa - du*Szz)
  Qv = -(fxx*Sxx+fyy*Syy+fzz*Szz)

  ; diagonal terms - Pm==magnetic pressure
  Pm = 0.5*(xup(Bx)^2 + yup(By)^2 + zup(Bz)^2)
  fxx += rho*xup(Ux)^2 - xup(Bx)^2 + Pg + Pm
  fyy += rho*yup(Uy)^2 - yup(By)^2 + Pg + Pm
  fzz += rho*zup(Uz)^2 - zup(Bz)^2 + Pg + Pm

  ; off-diagonal viscous part
  Sxy = (ddxdn(Uy)+ddydn(Ux))*0.5
  Syz = (ddydn(Uz)+ddzdn(Uy))*0.5
  Szx = (ddzdn(Ux)+ddxdn(Uz))*0.5
  fxy = - xdn1(ydn1(du))*Sxy
  fyz = - ydn1(zdn1(du))*Syz
  fzx = - zdn1(xdn1(du))*Szx
  Qv -= rho*(xup1(yup1(fxy*Sxy))+yup1(zup1(fyz*Syz))+zup1(xup1(fzx*Szx)))

  ; off-diagnoal part of Reynolds stress, without the density factor
  fxy += ydn(Ux)*xdn(Uy)
  fyz += zdn(Uy)*ydn(Uz)
  fzx += xdn(Uz)*zdn(Ux)

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

  ; advective E-field
  Ex = Ex + ydn(Uz)*zdn(By) - zdn(Uy)*ydn(Bz)
  Ey = Ey + zdn(Ux)*xdn(Bz) - xdn(Uz)*zdn(Bx)
  Ez = Ez + xdn(Uy)*ydn(Bx) - ydn(Ux)*xdn(By)

  dBxdt = ddzup(Ey) - ddyup(Ez)
  dBydt = ddxup(Ez) - ddzup(Ex)
  dBzdt = ddyup(Ex) - ddxup(Ey)

  ; Conservation of entropy, from co-moving per unit mass: d(ss) = Q/T = Q*rhp/Pg
  if (gamma eq 1.0) then begin
    dsdt = 0.*s
  endif else begin
    ss = s / rho
    temp = Pg/rho
    dsdt = (Qv+QJ)/temp $
      - ddxup1(xdn(ss)*px - rdx*ddxdn(ss)*xdn1(du) - fdifx*xdn(ss)) $
      - ddyup1(ydn(ss)*py - rdy*ddydn(ss)*ydn1(du) - fdify*ydn(ss)) $
      - ddzup1(zdn(ss)*pz - rdz*ddzdn(ss)*zdn1(du) - fdifz*zdn(ss))
  endelse

  ; fast mode speed relative to moving gas
  vf = sqrt((Pm+gamma*Pg+Pa)/rho)+sqrt(px^2+py^2+pz^2)/rho

END
;
PRO get_components, sol, rho, px, py, pz, Bx, By, Bz, s
  rho=sol.r
  px=sol.px
  py=sol.py
  pz=sol.pz
  Bx=sol.Bx
  By=sol.By
  Bz=sol.Bz
  S=sol.s
END
;
PRO set_components, sol, rho, px, py, pz, Bx, By, Bz, s
  sol = {r:rho, s:s, px:px, py:py, pz:pz, Bx:Bx, By:By, Bz:Bz}
END
; Time stepping
PRO timestep, sol, dx, dy, dz, dt, gamma, vf, time_evol, sol_old=sol_old

  get_components, sol, rho, px, py, pz, Bx, By, Bz, s

  ; Euler timestep
  if (time_evol eq 1) then begin
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, gamma, vf, $
                drhodt, dpxdt, dpydt, dpzdt, dBxdt, dBydt, dBzdt, dsdt
    rho += dt * drhodt
    s   += dt * dsdt
    px  += dt * dpxdt
    py  += dt * dpydt
    pz  += dt * dpzdt
    Bx  += dt * dBxdt
    By  += dt * dBydt
    Bz  += dt * dBzdt
  endif

  ; mid-point method
  if (time_evol eq 2) then begin
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, gamma, vf, $
                drhodt, dpxdt, dpydt, dpzdt, dBxdt, dBydt, dBzdt, dsdt
    rhom += 0.5*dt * drhodt
    sm   += 0.5*dt * dsdt
    pxm  += 0.5*dt * dpxdt
    pym  += 0.5*dt * dpydt
    pzm  += 0.5*dt * dpzdt
    Bxm  += 0.5*dt * dBxdt
    Bym  += 0.5*dt * dBydt
    Bzm  += 0.5*dt * dBzdt
    update_mhd, rhom, pxm, pym, pzm, Bxm, Bym, Bzm, Sm, dx, dy, dz, gamma, vf, $
                drhodt, dpxdt, dpydt, dpzdt, dBxdt, dBydt, dBzdt, dsdt
    rho += dt * drhodt
    s   += dt * dsdt
    px  += dt * dpxdt
    py  += dt * dpydt
    pz  += dt * dpzdt
    Bx  += dt * dBxdt
    By  += dt * dBydt
    Bz  += dt * dBzdt
  endif

  ; 2nd order A-B method
  if (time_evol eq 3) then begin
    update_mhd, rho, px, py, pz, Bx, By, Bz, S, dx, dy, dz, gamma, vf, $
                drhodt, dpxdt, dpydt, dpzdt, dBxdt, dBydt, dBzdt, dsdt
    rho += 1.5*dt * drhodt
    s   += 1.5*dt * dsdt
    px  += 1.5*dt * dpxdt
    py  += 1.5*dt * dpydt
    pz  += 1.5*dt * dpzdt
    Bx  += 1.5*dt * dBxdt
    By  += 1.5*dt * dBydt
    Bz  += 1.5*dt * dBzdt

    if (n_elements(sol_old) eq 0) then sol_old = sol

    get_components, sol_old, rhoo, pxo, pyo, pzo, Bxo, Byo, Bzo, so

    update_mhd, rhoo, pxo, pyo, pzo, Bxo, Byo, Bzo, So, dx, dy, dz, gamma, vf, $
                drhodt, dpxdt, dpydt, dpzdt, dBxdt, dBydt, dBzdt, dsdt

    rho += -0.5*dt * drhodt
    s   += -0.5*dt * dsdt
    px  += -0.5*dt * dpxdt
    py  += -0.5*dt * dpydt
    pz  += -0.5*dt * dpzdt
    Bx  += -0.5*dt * dBxdt
    By  += -0.5*dt * dBydt
    Bz  += -0.5*dt * dBzdt

    sol_old = sol
  endif

  set_components, sol, rho, px, py, pz, Bx, By, Bz, s
END

stagger_2nd ; get stagger operators
;stagger_6th ; get stagger operators
n = 128     ; grid size
s = 1.      ; box size
dx = s/n & dy = s/n & dz=1.
init_derivs, dx=dx, dy=dy, dz=dz

; Initial conditions
gamma = 1.000001d0
gamma = 5./3.

; launch shock from a central gaussian blob in entropy
w_gauss = 0.1
p_gauss = 4.
amplitude = 1.0

xx = (dindgen(n)-n/2+0.5)/n*s
yy = xx

Bx = dblarr(n,n) & By=Bx & Bz=Bx
px = Bx          & py = Bx & pz = Bx
rho = Bx + 1.
s = Bx
for i=0l,n-1l do $
  s[*,i] = amplitude*exp(-(xx^2+yy[i]^2)/w_gauss^2)^p_gauss

; Run benchmark
dt = 0.0007
t = 0.0
tend = 0.1

window, 0, xsize=512, ysize=512
imsz, 256, 256

time_evol=3
every=1
it=0
vfmax = 0.

set_components, sol, rho, px, py, pz, Bx, By, Bz, s

while (t lt tend) do begin
  timestep, sol, dx, dy, dz, dt, gamma, vf, time_evol, sol_old=sol_old
  vfmax = max(vf) > vfmax
  t = t + dt

  if it mod every eq 0 then begin
    image, sol.r, 0, title='rho t='+str(t), /minmax
    image, sol.px, 1,  title='px', /minmax
    image, sol.py, 2,  title='py', /minmax
    image, sol.s, 3,   title='s', /minmax
    wait, 0.01
  end
  it=it+1

  ;stat, rho, l='rho'
  ;stat, px , l='Px'
  ;stat, S  , l='S'
  ;stat, P  , l='P'
endwhile

sol_old=0

print,'max Courant number:',vfmax*dt/dx

if (0b) then begin
ws,0
!p.multi=[0,2,2]
!x.style=3
!y.style=3
plot, sol.rho, ps=-1,tit='rho time='+str(t,format='(f5.3)'), symsize=0.3
plot, sol.px, ps=-1, tit='px', symsize=0.3
;ux = px/exp(xdn(alog(rho)))
;plot, ux, ps=-1, tit='ux', symsize=0.3
plot, sol.by, ps=-1, tit='by', symsize=0.3
plot, sol.S/sol.rho, ps=-1, tit='s', symsize=0.3,yst=3

ws,1
;!p.multi=[0,2,2]
plot, sol.rho, xr=5+[0,20], psym=-1,yst=3,xst=3
plot, sol.px, xr=5+[0,20], psym=-1,yst=3,xst=3
endif

delvar, sol_old
END
