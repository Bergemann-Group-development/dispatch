;+ -----------------------------------------------------------------------------
 PRO make_linear_table, infile, outfile, plot=plot
; make a table with linear instead of log variables
;- -----------------------------------------------------------------------------
  common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab
  common ctable2,mrho,iupdte,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
  common ctable3,tmean,tamp,xcorr,thmin,thmax,dth
  ;-----------------------------------------------------------------------------
  ; For those variables in the table that are logarithmic, change to linear
  ;-----------------------------------------------------------------------------
  default,infile,'table.dat'
  default,outfile,'table_linear.dat'
  eos_table, infile
  ;-----------------------------------------------------------------------------
  ; For those variables in the table that are logarithmic, change to linear
  ;-----------------------------------------------------------------------------
  for irho=0,mrho-1 do begin
    ntab=mtab[irho]
    i=itab[irho]
    for iv=0,nvar-1 do begin
      if iv eq 2 then begin
        i=i+ntab*3
      end else begin
        tab[i:i+ntab-1]=exp(tab[i:i+ntab-1])
        j=i+ntab
        tab[j:j+ntab-1]=tab[j:j+ntab-1]*tab[i:i+ntab-1]
        j=j+ntab
        tab[j:j+ntab-1]=tab[j:j+ntab-1]*tab[i:i+ntab-1]
        i=j+ntab
      end
    end
  end
  if keyword_set(plot) then begin
    i=itab[plot]-1
    n=mtab[plot]
    window,4,xsize=1200,ysize=600
    !p.multi=[0,4,2]
    !p.charsize=1.6
    ee=eemin[plot]+dee[0]*findgen(n)
    print,'rho =',rhm[plot]
    print,'ee =',eemin[plot],ee[0],dee[0],ee[n-1],eemax[plot]
    for iv=0,nvar-1 do begin
      if iv eq 2 then f=tab[i:i+n-1] else f=alog(tab[i:i+n-1])
      plot,ee,f,psym=-1
      i=i+3*n
    end
  end
  ;-----------------------------------------------------------------------------
  ; Write the output file
  ;-----------------------------------------------------------------------------
  openw,u,/get,/f77,outfile
  writeu,u,mrho,iupdte,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
  writeu,u,tmean,tamp,rhm,xcorr,thmin,thmax,dth,eemin,eemax,dee,itab,mtab
  writeu,u,n_elements(tab)
  writeu,u,tab
  free_lun,u
END
