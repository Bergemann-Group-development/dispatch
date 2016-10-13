PRO plottablim
;
   common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab
   
   nr=n_elements(rhm)
   plot,[rhm(0),rhm(nr-1)],[min(eemin),max(eemax)],/nodata,/xlog,yst=0
   for ir=0,nr-1 do begin
     oplot,rhm(ir)*[1,1],[eemin(ir),eemax(ir)]
   end
   
   tt1=12000.
   rho1=10.
   
   ir=alog(rho1/rhm(0))/drhm+.5
   ee=0.5*(eemin(ir)+eemax(ir))
   tmake,tt1,alog(rho1),ee
   ;print,ir,eemin(ir),ee
   
   xyouts,/norm,0.25,0.75,'E(!7q!X=10,T=12000) ='+string(ee,format='(f6.3)')
END