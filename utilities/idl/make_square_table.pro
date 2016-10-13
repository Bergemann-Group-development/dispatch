;+ -----------------------------------------------------------------------------
 PRO make_square_table, infile, outfile, dlr=dlr_new, dee=dee_new
; make a table with linear instead of log variables
;- -----------------------------------------------------------------------------
  common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab
  common ctable2,mrho,iupdte,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
  common ctable3,tmean,tamp,xcorr,thmin,thmax,dth
  ;-----------------------------------------------------------------------------
  ; For those variables in the table that are logarithmic, change to linear
  ;-----------------------------------------------------------------------------
  default,infile ,'table_linear.dat'
  default,outfile,'table_square.dat'
  default, dlr_new, 0.1
  default, dee_new, 0.01
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
    i1 = long((lr-lr0)/drhm+0.5)
    ee2 = eemin[i1]
    ee3 = eemax[i1]
    iee0 = long((ee2-ee0)/dee_new+0.5)
    iee1 = long((ee3-ee0)/dee_new+0.5)
    nee = iee1-iee0+1
    iee = iee0+lindgen(nee)
    ee = ee0 + dee_new*iee
    lr = replicate(lr0 + ilr*dlr_new,nee)
    for iv=0,nvar-1 do begin
      v = lookup(lr,ee,iv=iv)
      table[iee,ilr,iv] = v
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
