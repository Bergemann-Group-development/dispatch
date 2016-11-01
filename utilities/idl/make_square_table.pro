;+ -----------------------------------------------------------------------------
 PRO make_square_table, infile, outfile, dlr=dlr_new, dee=dee_new, plot=plot
; make a table with linear instead of log variables
;- -----------------------------------------------------------------------------
  common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab
  common ctable2,mrho,iupdte,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
  common ctable3,tmean,tamp,xcorr,thmin,thmax,dth
  ;-----------------------------------------------------------------------------
  ; For those variables in the table that are logarithmic, change to linear
  ;-----------------------------------------------------------------------------
  infile = 'table_linear.dat'
  outfile = 'table_square.dat'
  dlr_new = 0.1
  dee_new = 0.01
  eos_table, infile
  print,'dee:',dee[0]
  mlr=mrho
  lr0=alog(rhm[0])
  lr1=alog(rhm[mlr-1])
  dlr=(lr1-lr0)/(mlr-1)
  print,'dlr:',dlr
  print,'dee_new:',dee_new
  print,'dlr_new:',dlr_new
  ;-----------------------------------------------------------------------------
  ; Find the smallest eemin and the largest eemax
  ;-----------------------------------------------------------------------------
  ee0 = min(eemin)
  ee1 = max(eemax)
  nee = long((ee1-ee0)/dee_new+1.5)
  nlr = long((lr1-lr0)/dlr_new+1.5)
  table = fltarr(nee,nlr,nvar)
  print,'table dims:',nee,nlr,nvar
  print,'table size:',nee*float(nlr)*nvar*4d0/1024.^2,' MB'
  ;-----------------------------------------------------------------------------
  ; For those variables in the table that are logarithmic, change to linear
  ;-----------------------------------------------------------------------------
  for ilr=0,nlr-1 do begin
    lr = lr0 + ilr*dlr_new
    i1 = long((lr-lr0)/drhm)
    ee2 = (eemin[i1] > eemin[i1+1])
    ee3 = (eemax[i1] < eemax[i1+1])
    iee0 = long((ee2-ee0)/dee_new+1.0)
    iee1 = long((ee3-ee0)/dee_new)
    n = iee1-iee0+1
    iee = iee0+lindgen(n)
    ee = ee0 + dee_new*iee
    lr = replicate(lr0 + ilr*dlr_new,n)
    for iv=0,nvar-1 do begin
      v = lookup(lr,ee,iv=iv)
      table[iee,ilr,iv] = v
    end
  end
  if keyword_set(plot) then begin
    window,8,xsize=1200,ysize=600
    !p.multi=[0,4,2]
    ilr=long((alog(rhm[plot])-lr0)/dlr_new)
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
  ;-----------------------------------------------------------------------------
  ; Write the output file
  ;-----------------------------------------------------------------------------
  openw,u,/get,/f77,outfile
  writeu,u,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd  ; basic params
  writeu,u,mlr,lr0,dlr                                    ; table limits
  writeu,u,eemin,eemax,dee[0]                             ; table limits
  writeu,u,nlr,lr0,dlr_new                                ; ln(rho) scale
  writeu,u,nee,ee0,dee_new                                ; ee scale
  writeu,u,table                                          ; table[nee,nlr,nvar]
  free_lun,u
END
