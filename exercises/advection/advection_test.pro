; $Id$
; ------------------------------------------------------------------------------
; The initial condition is a density which is a Gaussian of width w
; ------------------------------------------------------------------------------
FUNCTION initial_condition,n,w,center
  x = rebin(reform(findgen(n)-n/2,n,1),n,n)-center[0]
  y = rebin(reform(findgen(n)-n/2,1,n),n,n)-center[1]
  rho = exp(-(x^2+y^2)/w^2)
  ux = 1.0
  uy = 0.5
  return, {rho:rho, ux:ux, uy:uy}
END
; ------------------------------------------------------------------------------
; Derivatives evaluated with central differences
; ------------------------------------------------------------------------------
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
  drho0 = - (dt/f.dx)*(ddx(f.rho*f.ux) + ddy(f.rho*f.uy))
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
PRO advection_test, n_step, courant_number
  ic = initial_condition(64,10.0,[0.,0.])
  if n_elements(n_step) eq 0 then n_step=300
  if n_elements(courant_number) eq 0 then courant_number=0.10
  f = {rho:ic.rho, ux:ic.ux, uy:ic.uy, dx:1.0, courant:courant_number,t:0.0}
  wsz=512
  device,decompose=0
  ws,4
  umax=sqrt(f.ux^2+f.uy^2)
  plot,[0,n_step]*f.courant*umax,xst=1,[0,0.11],yst=1,/nodata
  oplot,[0,n_step],[1,1]*0.1,line=1
  loadct,39
  for i=0,n_step do begin
    f = time_step(f)
    fc = initial_condition(64,10.0,[f.ux,f.uy]*f.t)
    ;image,f.rho,title=str(i)
    e=rms(f.rho-fc.rho)/rms(fc.rho)
    plots,f.t,e,psy=1
    if e gt 0.1 then break
    wait,0.01
  end
  print,(i+1)*f.courant,e,f.courant,(i+1)
  window,xsize=wsz,ysiz=wsz
  imsz,wsz
  image,f.rho,title=str(i)
  wdupl,1
  image,fc.rho
END
