PRO test_plot1, ilr, iv
  common ctab2,nlr,lr0,dlr_new
  common ctab3,nee,ee0,dee_new
  common ctab4,table
  default,ilr,100
  default,iv,2
  plot,table[*,ilr,iv]
END

;+
 PRO eos_square_table,file,plot=plot,_extra=_extra
;
;  Load EOS table from file compatible with the STAGGER code format
;
;-
  common ctab0,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
  common ctab1,mlr,lnr0,dlr,eemin,eemax,dee
  common ctab2,nlr,lr0,dlr_new
  common ctab3,nee,ee0,dee_new
  common ctab4,table
  ;-----------------------------------------------------------------------------
  ; Read the table  file
  ;-----------------------------------------------------------------------------
  default,file,'table_square.dat'
  openr,u,/get,/f77,file
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
  readu,u,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd  ; basic params
  mlr=0L
  lnr0=0.0
  dlr=0.0
  readu,u,mlr,lnr0,dlr                                   ; table limits
  eemin=fltarr(mlr)
  eemax=eemin
  dee=0.0
  readu,u,eemin,eemax,dee                                ; table limits
  nlr=0L
  lr0=0.0
  dlr_new=0.0
  readu,u,nlr,lr0,dlr_new                                ; ln(rho) scale
  nee=0L
  ee0=0.0
  dee_new=0.0
  readu,u,nee,ee0,dee_new                                ; ee scale
  table=fltarr(nee,nlr,nvar)
  readu,u,table                                          ; table[nee,nlr,nvar]
  free_lun,u
  ;
  if keyword_set(plot) then begin
    window,12,xsize=1200,ysize=600
    !p.multi=[0,4,2]
    lnr=lnr0+plot*dlr
    ilr=long((lnr-lr0)/dlr_new)
    i0=long((eemin[plot]-ee0)/dee_new+1.0)
    i1=long((eemax[plot]-ee0)/dee_new)
    ee=ee0+dee_new*(i0+findgen(i1-i0+1))
    print,eemin[plot],ee[0],dee[0],ee[i1-i0],eemax[plot]
    for iv=0,nvar-1 do begin
      f=table[i0:i1,ilr,iv]
      if iv ne 2 then f=alog(f)
      plot,ee,f,psym=-1
    end
  end
END
