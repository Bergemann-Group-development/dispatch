;+
; $Id: lookup.pro,v 1.13 2015/01/08 21:23:28 aake Exp $
;-------------------------------------------------------------------------------
 FUNCTION lookup_square,lnrho,ee,dvdecr,dvdrce,iv=iv,mz=mz,status=status, $
    w1=w1, w2=w2, mu_corona=mu_corona, plot=plot,vcorners=corners, px=px, $
    py=py, bilinear=bilinear, simple=simple, verbose=verbose
;-
  common ctab0,nvar,mbox,xmin,dbox,ul,ut,ur,eps,tff,grv,abnd
  common ctab1,mlr,lnr0,dlr,eemin,eemax,dee
  common ctab2,nlr,lr0,dlr_new
  common ctab3,nee,ee0,dee_new
  common ctab4,table

  if n_params() eq 0 then begin
    print,'lookup_square,lnrho,ee[,dvdecr,dvdrce],mz=mz,iv=iv'
    print,'  dvdecr is the partial der of the answ wrt ee at const lnr'
    print,'  dvdrce is the partial der of the answ wrt lnr at const ee'
    return,0
  endif

  s=size(lnrho)
  if s[s[0]+2] gt 15e6 then begin
    if s[0] eq 3 then begin
      print,s[3],' slices'
      f      = make_array(size=s,/nozero)
      dvdecr = make_array(size=s,/nozero)
      dvdrce = make_array(size=s,/nozero)
      for i=0,s[3]-1 do begin
        f[*,*,i]=lookup(lnrho[*,*,i],ee[*,*,i],iv=iv,mz=mz,d1,d2)
        dvdecr[*,*,i] = d1
        dvdrce[*,*,i] = d2
      end
      return, f
    end else if s[0] eq 2 then begin
      f =      make_array(size=s,/nozero)
      dvdecr = make_array(size=s,/nozero)
      dvdrce = make_array(size=s,/nozero)
      for i=0,s[2]-1 do begin
        f[*,i]=lookup(lnrho[*,i],ee[*,i],iv=iv,mz=mz ,d1,d2)
        dvdecr[*,i] = d1
        dvdrce[*,i] = d2
      end
      return, f
    end
  end

  if n_elements(iv) eq 0 then iv=0
  ;
  ; Density index
  ;
  if s[0] gt 0 then nout=product(s[1:s[0]]) else nout=1
  lr=reform(lnrho,nout)
  plr=(lr-lr0)/dlr_new
  ilr=((long(plr) > 0) < (nlr-2))
  plr=ilr>plr<(ilr+1)
  plr=plr-ilr
  qlr=1.0-plr
  ;
  ; Energy index
  ;
  ee1=reform(ee,nout)
  pee=(ee1-ee0)/dee_new
  iee=((long(pee) > 0) < (nee-2))
  pee=iee>pee<(iee+1)
  pee=pee-iee
  qee=1.0-pee
  ;
  ; Bi-linear interpolation
  ;
  out=qlr*(qee*table[iee,ilr  ,iv]+pee*table[iee+1,ilr  ,iv]) $
     +plr*(qee*table[iee,ilr+1,iv]+pee*table[iee+1,ilr+1,iv])
  if          s[0] eq 0 then begin
    out=out[0]
  end else if s[0] eq 1 then begin
    out=reform(out,s[1])
  end else if s[0] eq 2 then begin
    out=reform(out,s[1],s[2])
  end else if s[0] eq 3 then begin
    out=reform(out,s[1],s[2],s[3])
  end else begin
    out=0.0
    print,'ERROR: cannot handle '+strtrim(string(s[0]),2)+' dimensions'
  end

  status=0
  if status eq 0 then return,out
  ;
  ;  Ideal gas extrapolations
  ;
  ;  The EOS/UMA code counts energy relative to ionized hydrogen.  To ensure
  ;  positive values of the internal energy (for technical reasons), the tabcmp.f
  ;  code adds 1.5 times an estimate of the ionization energy of hydrogen, thus
  ;  compensating for recombination and also for molecule formation (mainly H2).
  ;  However, things are complicated by other contributions to the internal
  ;  energy, so the simplest approach is to use the last table value (which is
  ;  automatically available) and extraplate from that.
  ;
  ;  To compute temperatures for ionized gas we need to subtract off that offset
  ;  again, and then assume that the gas has a molecular weight corresponding to
  ;  fully ionized solar composition gas, which means 1.2 electrons and 1.1 atoms
  ;  carry a weight of 1.43 atomic units, so mu = 1.43/2.3 = 0.62. When helium is
  ;  neutral the ratio is 1.4/2.1 = 0.68.  The 0.03 extra weight per proton is
  ;  from elements heavier than helium.
  ;
  default, T0_ideal, 3e4                                ; critical choice!
  default, T_margin, 100.                               ; inside margin
  default, mu_corona, 0.593                             ; see above
  ul=1e8                                                ; length unit
  ut=1e2                                                ; time unit
  ur=1e-7                                               ; density unit
  up=ur*(ul/ut)^2                                       ; pressure unit
  kb=1.38e-16                                           ; boltzmann
  mp=1.67e-24                                           ; proton mass
  dttdee=2./3.*(ul/ut)^2*mu_corona*1.67e-24/1.38e-16    ; ideal gas outside
  eeoffset=eemax[np]                                    ; offset to start from
  T_margin=100.                                         ; 100 K overlap (inside!)
  ee_margin=T_margin/dttdee                             ; overlap in ee

  arg=(-75.) > ((ee-eeoffset)/ee_margin+4.) < 75.       ; is 4 at eeoffset
  wf = 0.5*(1.0-tanh(arg))                              ; goes to ~0 at eeoffset

  if iv eq 2 then begin                                 ; temperature
      tabvar0 = tabvar                                  ; for plotting
      dvdecr0 = dvdecr

      T_eemax = tabvar                                  ; T at end of table
      T_ideal = T_eemax + (ee-eeoffset)*dttdee          ; ideal gas extrapolation
      tabvar  = tabvar*wf + (1.0-wf)*T_ideal            ; combine

      dTdee_ideal = dttdee                              ; ee derivative
  ;          dvdecr  = dvdecr*wf + (1.0-wf)*dTdee_ideal       ; values outside are OK

      dTdlnr_ideal = 0.                                 ; lnr derivative
      dvdrce = dvdrce*wf + (1.0-wf)*dTdlnr_ideal

      if keyword_set(plot) and n_elements(ee) gt 1 then begin
        multi=!p.multi & !p.multi=[0,2,1]
        plot,ee,xst=3,wf,yr=[-2,2]
        oplot,ee,0.5*(1.0-tanh(wf1/w1)),line=1,psy=-1
        oplot,ee,0.5*(1.0-tanh(wf2/w2)),line=2,psy=-2
        plot,ee,tabvar & oplot,ee,tabvar0,line=2
        n=n_elements(ee) & ee0=ee[n/2]
        plots,ee0,tabvar[n/2],psym=6
        oplot,ee,tabvar[n/2]+(ee-ee0)*dvdecr0[n/2],line=1
        oplot,ee,tabvar[n/2]+(ee-ee0)*dvdecr [n/2],line=2
        !p.multi=multi
      end
  end else if iv eq 0 then begin                        ; pressure
      ee0 = 17.5                                        ; empirical
      temp_fact = (ee-ee0)/(eemax[np]-ee0) > 0.1        ; temperature factor
      lnp_ideal = tabvar + alog(temp_fact > 0.1)        ; table end + extrapolation
      tabvar    = tabvar*wf + (1.0-wf)*lnp_ideal        ; combine

  end else if iv eq 3 then begin
      nel_ideal = lnrho+alog(ur*0.5/(mu_corona*mp))     ; 50% electrons assumed
      tabvar    = tabvar*wf + (1.0-wf)*nel_ideal        ; combine
  end
  ;
  return,tabvar
end
