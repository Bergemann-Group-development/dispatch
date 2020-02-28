; $Id$
;-----------------------------------------------------------------------
; These procedures may be used to illustrate the growth of random
; velocity dispersion with time, due to particle-grid interaction.
; With a good code and good parameter values the growth should be
; very slow, even when the initial Debye length is much smaller than
; the mesh size.
;-----------------------------------------------------------------------

PRO plot_saved_growth, files, psyms=psyms, max=max, _extra=_extra
  default, max, 1e-3
  restore,files[0]
  n=n_elements(noise)
  plot, [0,n],[1e-6,1]*max, /ylog, /nodata, _extra=_extra
  for i=0,n_elements(files)-1 do begin
    restore,files[i]
    oplot, noise, psym=psyms[i]
    oplot, (noise-noise[0]) > 1e-30, psym=psyms[i]
  end
END

pro plot_growth, n, max, oplot=oplot, save=save, _extra=_extra
default, n, 100
default, max, 1e-1
!p.charsize=1.5
if not keyword_set(oplot) then begin
  p=rp(0)
  ie=(where(trim(p.s.pname) eq "electron"))[0]  ; electron index
  temp=p.s.temperature[ie]                      ; electron temperature
  mass=p.s.mass[ie]                             ; electron mass
  vrms=sqrt(temp/mass)                          ; rms velocity
  help,temp,vrms
  default,max,3*vrms
  plot, [0,n],[1e-6,1]*max, /ylog, /nodata, _extra=_extra
endif
p=rp(0)
noise=fltarr(n+1)
noise[0]=(rms0=rms(p.p[*,2,0]))
for i=1,n do begin
  p=rp(i)
  noise[i]=rms(p.p[*,2,0])
  plots,i,(noise[i]-rms0) > 1e-30,psy=1
  plots,i,(noise[i]) > 1e-30,psy=1
end
if n_elements(save) gt 0 then save,noise,file=save
END
