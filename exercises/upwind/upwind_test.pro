; $Id$
; ------------------------------------------------------------------------------
; The initial condition is a density which is a Gaussian of width w
; ------------------------------------------------------------------------------
FUNCTION initial_condition,n,w,center
  x = findgen(n)-center[0]
  y = findgen(n)-center[1]
  x = ((x + 10*n) mod n) - n/2
  y = ((y + 10*n) mod n) - n/2
  x = rebin(reform(x,n,1),n,n)
  y = rebin(reform(y,1,n),n,n)
  rho = exp(-(x^2+y^2)/w^2)
  ux = 1.0
  uy = 0.0
  return, {rho:rho, ux:ux, uy:uy}
END
; ------------------------------------------------------------------------------
; Derivatives evaluated with central differences
; ------------------------------------------------------------------------------
FUNCTION ddxup, f
  return, (shift(f,-1,0)-f)
END
FUNCTION ddxdn, f
  return, (f-shift(f,1,0))
END
FUNCTION ddx, f
  return, (shift(f,-1,0)-shift(f,1,0))*0.5
END
FUNCTION ddy, f
  return, (shift(f,0,-1)-shift(f,0,1))*0.5
END
; ------------------------------------------------------------------------------
; Simplest possible time step -- first order in time
; ------------------------------------------------------------------------------
FUNCTION time_step,f
  umax = sqrt(max(f.ux^2+f.uy^2))
  dt = f.courant*f.dx/umax
  f.t = f.t+dt
  delta = dt*f.ux*1.5
  a = (1.0+delta)*0.5
  b = (1.0-delta)*0.5
  drho0 = - (dt/f.dx)*((a*ddxdn(f.rho*f.ux)+b*ddxup(f.rho*f.ux)) + ddy(f.rho*f.uy))
  ;drho0 = - (dt/f.dx)*(ddx(f.rho*f.ux) + ddy(f.rho*f.uy))
  rho0 = f.rho+drho0
  drho1 = - (dt/f.dx)*(ddx(rho0*f.ux) + ddy(rho0*f.uy))
;  f.rho = f.rho + 0.5*(drho0+drho1)
  f.rho = f.rho + drho0
  return, f
END
; ------------------------------------------------------------------------------
; This test should just move the Gaussian, without changing the shape, and w/o
; adding any erroneous new features.  Run it first for 100 steps, then 300, then
; 1000 steps.  Vary the Courant number (the fraction of a mesh the profiles is
; moved per step), and try to move the shape as far as possible (disregarding
; the number of time steps needed).
; ------------------------------------------------------------------------------
PRO upwind_test, courant_number=courant_number, n_step=nstep, width=width
  n=64
  if n_elements(width) eq 0 then width=10.
  if n_elements(n_step) eq 0 then n_step=10000
  if n_elements(courant_number) eq 0 then courant_number=0.2
  ic = initial_condition(n,width,[0.,0.])
  f = {rho:ic.rho, ux:ic.ux, uy:ic.uy, dx:1.0, courant:courant_number,t:0.0}
  wsz=256
  device,decompose=0
  window,5,xs=wsz,ys=wsz & imsz,wsz
  umax=sqrt(f.ux^2+f.uy^2)
  window,0,xs=600,ys=400
  plot,[0,1]*umax*width^2*3/(1+f.courant/0.12),xst=1,[0,0.11],yst=1,/nodata
  oplot,[0,n_step],[1,1]*0.1,line=1
  loadct,39
  for i=0,n_step do begin
    f = time_step(f)
    fc = initial_condition(n,width,[f.ux,f.uy]*f.t)
    ws,5 & image,/gr,f.rho,title=str(i) & ws,0
    e=rms(f.rho-fc.rho)/rms(fc.rho)
    plots,f.t,e,psy=1,syms=.1
    if e gt 0.1 then break
    wait,0.01
    ;print,i,e
  end
;  print,(i+1)*f.courant,e,f.courant,(i+1)
  print,(i+1)*f.courant,f.courant,(i+1)
  window,1,xs=600,ys=400
  plot,fc.rho[*,n/2],xst=3,line=2,yst=3 & zeroline
  oplot,f.rho[*,n/2],psy=-1
  ;window,xsize=wsz,ysiz=wsz
  ;imsz,wsz
  ;image,f.rho,title=str(i)
  ;wdupl,1
  ;image,fc.rho
END
