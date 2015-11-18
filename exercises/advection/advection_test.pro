; $Id$
; ------------------------------------------------------------------------------
; The initial condition is a density which is a Gaussian of width w
; ------------------------------------------------------------------------------
FUNCTION initial_condition,n,w
  x = rebin(reform(findgen(n)-n/2,n,1),n,n)
  y = rebin(reform(findgen(n)-n/2,1,n),n,n)
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
  f.rho = f.rho - (dt/f.dx)*(ddx(f.rho*f.ux) + ddy(f.rho*f.uy))
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
  ic = initial_condition(64,10.0)
  if n_elements(n_step) eq 0 then n_step=100
  if n_elements(courant_number) eq 0 then courant_number=0.2
  f = {rho:ic.rho, ux:ic.ux, uy:ic.uy, dx:1.0, courant:courant_number}
  wsz=512
  device,decompose=0
  loadct,39
  window,xsize=wsz,ysiz=wsz
  imsz,wsz
  for i=0,n_step do begin
    f = time_step(f)
    image,f.rho,title=str(i)
    wait,0.01
  end
END
