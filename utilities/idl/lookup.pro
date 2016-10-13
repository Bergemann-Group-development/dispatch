; $Id: lookup.pro,v 1.13 2015/01/08 21:23:28 aake Exp $
;-------------------------------------------------------------------------------
FUNCTION nint,x
  ix=floor(x+.5)
  return,ix
END

;-------------------------------------------------------------------------------
PRO constant_entropy, lnr, ee, deedlnr_s, dlnpdlnr_s
;  T dS = dee - P dV = dee + P/rho * dlnr
;  e = rho*ee
;  dedr_s = rho*deedr_s + ee = deedlnr_s + ee

  pg = exp(lookup(lnr,ee,dlnpdlnr_ee,dlnpdee_lnr))
  deedlnr_s = pg/exp(lnr)
  dlnpdlnr_s = dlnpdlnr_ee + dlnpdee_lnr*deedlnr_s
END

;-------------------------------------------------------------------------------
PRO constant_pressure, lnr, ee, deedlnr_pg,  dssdlnr_pg, tt=tt
;  T dS = dee - P dV = dee + P/rho * dlnr
;  e  = rho*ee
;  de = rho*dee + drho*ee = rho*dee + rho*ee*dlnrho

  pg = exp(lookup(lnr,ee,dlnpdlnr_ee,dlnpdee_lnr))
  deedlnr_pg = -dlnpdlnr_ee/dlnpdee_lnr
  if n_elements(tt) gt 0 then dssdlnr_pg = (pg/exp(lnr) + deedlnr_pg)/tt
END

;-------------------------------------------------------------------------------
FUNCTION lookup,lnrho,ee,dvdecr,dvdrce,iv=iv,mz=mz,status=status, $
                w1=w1, w2=w2, mu_corona=mu_corona, plot=plot, $
                corners=corners, px=px, py=py, bilinear=bilinear, $
                simple=simple,verbose=verbose

      common cdat,x,y,z,nx,ny,nz,mp,ntmax,date0,time0
      common ctable,ttmean,rhm,drhm,eemin,eemax,dee,mtab,itab,tab

      if n_params() eq 0 then begin
      print,'lookup,lnrho,ee[,dvdecr,dvdrce],mz=mz,iv=iv'
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
      if n_elements(mz) eq 0 then mz=n_elements(rhm)

      rhm1=alog(rhm(0))
      rhm2=alog(rhm(mz-1))
      drhm=(rhm2-rhm1)/(mz-1)

      rhm11=lnrho
      drhm1=rhm11
      if s(0) eq 3 then begin
        rhm11(*,*,*)=rhm1
        drhm1(*,*,*)=drhm
      endif else if s(0) eq 2 then begin
        rhm11(*,*)=rhm1
        drhm1(*,*)=drhm
      endif else begin
        rhm11=rhm1
        drhm1=drhm
      endelse
;
;  density index
;
      algrk=(lnrho-rhm11)/drhm1
      np=((long(algrk) > 0) < (mz-2))
      if iv eq 0 or iv eq 2 then algrk=np>algrk<(np+1)                  ; bracket
      px=algrk-np
      qx=1.-px
      pxqx=px*qx

      pxmax = max(px)
      qxmax = max(qx)
      status=0

      if pxmax gt 1.01 or qxmax gt 1.01 then begin
        status=status+1
        print,'warning: extrapolation in lnrho, pxmax =',pxmax,' qxmax =',qxmax
      endif
;
;  Interpolate the vars and dvar/dlnrho at smaller lnrho
;
      ntab=mtab(np)
      ltab=3*ntab
      eek=(ee-eemin(np))/dee(np)
      kee=((long(eek) > 0) < (ntab-2))
      if iv eq 0 or iv eq 2 then eek=kee>eek<(kee+1)                    ; bracket
      py=eek-kee
      qy=1.-py
      pyqy=py*qy
      if keyword_set(verbose) then print,'py ,kee  =',py,kee

      pymax = max(py)
      qymax = max(qy)
      if pymax gt 1.01 or qymax gt 1.01 then begin
        status=status+1
        ;print,'warning: extrapolation in ee, pymax =',pymax,' qymax =',qymax
      endif

      ik=itab(np)+kee-1
      ntab1=mtab(np+1)

      eek1=(ee-eemin(np+1))/dee(np+1)
      kee1=kee+floor((eemin(np)-eemin(np+1))/dee(np+1)+0.5)
      if iv eq 0 or iv eq 2 then eek1=kee1>eek1<(kee1+1)                ; bracket
      py1=eek1-kee1
      if keyword_set(verbose) then print,'py1,kee1 =',py1,kee1

      if max(abs(((py1-py + 1.5) mod 1.0)-0.5)) gt 0.1 then begin
        print,'address problem in lookup'
        ;stop
      end
      kee1=(kee1 > 0) < (ntab1-2)
      ik1=itab(np+1)+kee1-1
;
;  Get ik for variable desired
;
      if keyword_set(simple) then begin
        ik=ik+iv*ntab
        ik1=ik1+iv*ntab1
      end else begin
        ik=ik+iv*3*ntab
        ik1=ik1+iv*3*ntab1
      end
;
      dfx1=tab(ik1)-tab(ik)
      dfx2=tab(ik1+1)-tab(ik+1)
      dfy1=tab(ik+1)-tab(ik)
      dfy2=tab(ik1+1)-tab(ik1)
;
      if keyword_set(bilinear) or keyword_set(simple) then begin
        tabvar = $
          qy * (qx * tab(ik   ) + $
                px * tab(ik1  ))+ $
          py * (px * tab(ik1+1) + $
                qx * tab(ik +1))
      end else begin
        tabvar = $
          qy * (qx * (tab(ik           ) + $
               pxqx* (tab(ik +2*ntab   ) - dfx1) + $
               pyqy* (tab(ik +  ntab   ) - dfy1) ) + $
                px * (tab(ik1          ) - $
               pxqx* (tab(ik1+2*ntab1  ) - dfx1) + $
               pyqy* (tab(ik1+  ntab1  ) - dfy2) ) ) + $
          py * (px * (tab(ik1        +1) - $
               pxqx* (tab(ik1+2*ntab1+1) - dfx2) - $
               pyqy* (tab(ik1+  ntab1+1) - dfy2) ) + $
                qx * (tab(ik         +1) + $
               pxqx* (tab(ik +2*ntab +1) - dfx2) - $
               pyqy* (tab(ik +  ntab +1) - dfy1) ) )
      end
      if n_elements(ee) eq 1 then begin
        corners = [tab(ik), tab(ik+1), tab(ik1), tab(ik1+1)]
      end
;
;  Energy derivative
;
      w=where(ee gt eemax[np],status)
      if n_params() lt 3 and status eq 0 then return,tabvar

      dvdecr = $
               qx * (   - (tab(ik           )       ) + $
     (1.-4.*py+3.*py^2) * (tab(ik +  ntab   ) - dfy1) + $
                (-pxqx) * (tab(ik +2*ntab   ) - dfx1) + $
                          (tab(ik         +1)       ) - $
        (2.*py-3.*py^2) * (tab(ik +  ntab +1) - dfy1) + $
                   pxqx * (tab(ik +2*ntab +1) - dfx2) ) + $
               px * (     (tab(ik1        +1)       ) - $
        (2.*py-3.*py^2) * (tab(ik1+  ntab1+1) - dfy2) - $
                   pxqx * (tab(ik1+2*ntab1+1) - dfx2) - $
                          (tab(ik1          )       ) + $
     (1.-4.*py+3.*py^2) * (tab(ik1+  ntab1  ) - dfy2) - $
                (-pxqx) * (tab(ik1+2*ntab1  ) - dfx1) )
;
      dvdecr = dvdecr/dee(np)
      if n_params() lt 4 and status eq 0 then return,tabvar
;
;  Density derivative
;
      dvdrce = $
               qy * (   - (tab(ik           )       ) + $
     (1.-4.*px+3.*px^2) * (tab(ik +2*ntab   ) - dfx1) + $
                (-pyqy) * (tab(ik +  ntab   ) - dfy1) + $
                          (tab(ik1          )       ) - $
        (2.*px-3.*px^2) * (tab(ik1+2*ntab1  ) - dfx1) + $
                   pyqy * (tab(ik1+  ntab1  ) - dfy2) ) + $
               py * (     (tab(ik1        +1)       ) - $
        (2.*px-3.*px^2) * (tab(ik1+2*ntab1+1) - dfx2) - $
                   pyqy * (tab(ik1+  ntab1+1) - dfy2) - $
                          (tab(ik         +1)       ) + $
     (1.-4.*px+3.*px^2) * (tab(ik +2*ntab +1) - dfx2) - $
                (-pyqy) * (tab(ik +  ntab +1) - dfy1) )
;
;  Dimensional factors, because table values for density derivatives are
;  in table grid units.  The output values are d(f)/d(ee) and d(f)/d(lnr).
;
      dlnrho = drhm
      dvdrce = dvdrce/dlnrho
      if status eq 0 then return,tabvar
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
