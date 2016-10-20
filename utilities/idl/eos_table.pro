pro chtab,nz1
    common cdat,x,y,z,nx,ny,nz,mp,ntmax,date0,time0
    common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab

    nz=nz1
    rhm=rhm(0:n)
    eemin=eemin(0:n)
    eemax=eemax(0:n)
    dee=dee(0:n)
    mtab=mtab(0:n)
    itab=itab(0:n)
end

;+
PRO eos_table,file,sstab=sstab,quiet=quiet,plot=plot, $
    structure=str,_extra=_extra
;
;  Load EOS table from file compatible with the STAGGER code format
;
;-
    common cdat,x,y,z,nx,ny,nz,mp,ntmax,date0,time0
    common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab
    common mycom,cray,one,two
    common ctable2,mrho,iupdte,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
    common ctable3,tmean,tamp,xcorr,thmin,thmax,dth

;    default,file,'table.dat'
    close,2
    openr,2,file,_extra=_extra,/f77
    mrho=0L
    iupdte=0L
    nvar=0L
    mbox=0L
    xmin=0.
    dbox=0.
    ul=0.
    ur=0.
    ut=0.
    eps=0.
    tff=0.
    grv=0.
    abnd=0.
    readu,2,mrho,iupdte,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
    if not keyword_set(quiet) then print,'mrho,mbox,nvar=',mrho,mbox,nvar
    a=fltarr(mrho)
    tmean=a
    tamp=a
    rhm=a
    xcorr=a
    thmin=a
    thmax=a
    dth=a
    eemin=a
    eemax=a
    dee=a
    itab=lonarr(mrho)
    mtab=itab
    readu,2,tmean,tamp,rhm,xcorr,thmin,thmax,dth,eemin,eemax,dee,itab,mtab
    ntab=0L
    readu,2,ntab
    tab=fltarr(ntab)
    readu,2,tab
    close,2

    njon=nvar-2*mbox
    if njon eq 0 then njon=nvar-mbox

    str={nrho:mrho,nbox:mbox,nvar:nvar,njon:njon}

    rhm1=alog(rhm(0))
    rhm2=alog(rhm(mrho-1))
    drhm=(rhm2-rhm1)/(mrho-1)

    if keyword_set(sstab) then sstab
    ;if keyword_set(plot) then plottablim
    if keyword_set(plot) then begin
      i=itab[plot]-1
      n=mtab[plot]
      window,0,xsize=1200,ysize=600
      !p.multi=[0,4,2]
      !p.charsize=1.6
      ee=eemin[plot]+dee[0]*findgen(n)
      print,'ee =',eemin[plot],ee[0],dee[0],ee[n-1],eemax[plot]
      print,'rho =',rhm[plot]
      for iv=0,nvar-1 do begin
        plot,ee,tab[i:i+n-1],psym=-1,xstyle=3
        i=i+3*n
      end
      return
    end
END
