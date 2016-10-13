;+
 PRO sstab
;
;  Replace the 2nd slot (index 1) in the lookup table with S; the entropy
;
;-
    common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab

    nrho=n_elements(rhm)
    ;plot_oi,rhm,/nodata,[-.001,.0005],ytit='S',xtit='rho',yst=1

    dlnr=alog(rhm(nrho-1)/rhm(0))/(nrho-1)					; delta-lnr
    for nr=0,nrho-1 do begin
      ntab=mtab(nr)
      ltab=3*ntab
      kpg=itab(nr)-1+indgen(ntab)						; index to pg-data
      kab=kpg+ltab								; index to abs data
      ktt=kab+ltab								; index to temperature data
      tt=tab(ktt)								; temperatures
      dttdee=tab(ktt+ntab)							; (dT/dE)_rho
      pg=exp(tab(kpg))								; Pg
      dssdee=1./tt								; (dS/dE)_rho
      d2ssdee2=-dttdee/tt^2							; 2nd deriv
      dssdlnr=-pg/rhm(nr)/tt							; (dS/dlnr)_E
      ss=fltarr(ntab)
      dss=dee(nr)*(0.5*(dssdee+dssdee(1:*))+(d2ssdee2-d2ssdee2(1:*))/12.)

      if nr gt 0 then begin
        for k=1,ntab-1 do ss(k)=ss(k-1)+dss(k-1)				; variation with E
        if eemin(nr) gt eemin(nr-1) then begin
          k1=long((eemin(nr)-eemin(nr-1))/dee(nr)+.5)
          ss=ss+aver(ssp(k1:*)-ss(*)+0.5*dlnr*(dssdlnrp(k1:*)+dssdlnr(*)))
          ;if nr eq nrho-1 then stop
	  ; plot,ss,yr=[min(ss<ssp),max(ss>ssp)],tit=nr & oplot,ssp[k1:*] & wait,.1
        end else begin
          k1=long((eemin(nr-1)-eemin(nr))/dee(nr)+.5)
          ss=ss+aver(ssp(*)-ss(k1:*)+0.5*dlnr*(dssdlnrp(*)+dssdlnr(k1:*)))
	  ; plot,ss[k1:*],yr=[min(ss<ssp),max(ss>ssp)],tit=nr & oplot,ssp & wait,.1
        end
      end
      ssp=ss
      dssdlnrp=dssdlnr

      tab(kab)=ss
      dssdee=(ss(2:*)-ss)/2.
      dssdee=[2.*dssdee(0)-dssdee(1),dssdee,2.*dssdee(ntab-3)-dssdee(ntab-4)]
      tab(kab+ntab)=dssdee
      tab(kab+2*ntab)=dssdlnr*dlnr

      ;oplot,replicate(rhm(nr),ntab),ss,psym=3
    end
END
