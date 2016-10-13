; $Id: tmake.pro,v 1.11 2012/08/03 15:22:31 aake Exp $
;+
 PRO tmake,tt,lnr,ee,debug=debug,plot=plot,itmax=itmax,eps=eps,safe=safe
;  Compute, by iterative correction, the energy 'ee' consistent
;  with the input temperature 'tt' and log density 'lnr'.  The
;  table lookup must be activated first, using 'IDL> eos_table,tablefile'.
;-
common ctable,ttmean,rhm,drhm,eemin,eemax,deetab,mtab,itab,tab
common ctable,ttmean,rhm,drhm,eemin,eemax,deetab,mtab,itab,tab
   default,itmax,100
   default,eps,1e-4
   default,safe,100
   if keyword_set(plot) then plottablim
   nr=n_elements(rhm)
   ntt=n_elements(tt)
   if n_elements(ee) eq 0 then ee=make_array(size=size(tt))
   w=where(ee eq 0,nw)
   if nw gt 0 then begin
     ir=0>long((lnr(w)-alog(rhm(0)))/drhm)<(nr-1)
     ee(w)=0.5*(eemin(ir)+eemax(ir))
   endif
   for it=0,itmax do begin
     tt1=lookup(lnr,ee,dttdee,dttdlnr,iv=2)
     dif=tt1-tt
     dee=-dif/dttdee
     ir=((long((lnr-alog(rhm(0)))/drhm)+1) > 0) < (nr-1)
     dee=dee/(1.+abs(safe*dee/(eemax(ir)-eemin(ir))))
     ee=ee+dee
     if keyword_set(plot) then oplot,exp(lnr),ee
     if keyword_set(debug) then print, max(abs(dif))
     err=max(abs(dif)/tt)
     if (err lt eps) then return
   endfor
   print,'warning, not converged:', err, eps
END
