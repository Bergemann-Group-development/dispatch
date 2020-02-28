; $Id$
;-----------------------------------------------------------------------
; This procedure illustrates the degree to which particles that have
; an initial random velocity on top of the bulk motion continue in their
; path (or not).
;-----------------------------------------------------------------------
PRO plot_traces, ii
  default,ii,[0,100]                            ; dump range
  pgam=0.5                                      ; v0 in input.nml
  gam=sqrt(1.+pgam^2)                           ; rel. gamma
  v=pgam/gam                                    ; actual speed
  toutp=2.                                      ; particle dump interval
  dzdi=pgam/gam*toutp                           ; shift needed per dump
  p1=rp(ii[0])                                  ; particle data
  w=where(trim(p1.s.pname) eq "electron")       ; electron index
  ie=w[0]                                       ; scalar index
  temp=p1.s.temperature[ie]                     ; electron temperature
  mass=p1.s.mass[ie]                            ; electron mass
  vrms=sqrt(temp/mass)                          ; rms velocity
  range=5.*vrms*toutp*(ii[1]-ii[0])             ; plot range
  ind=indgen(100)                               ; number to plot
  plot,[-1,1]*range,[-1,1]*range,/nodata        ; frame
  ds=p1.s.ds                                    ; mesh size
  s1=sort(p1.i[*,0])                            ; sorted indices
  ;f='(10i5)'
  ;print,p1.i[s1[ind],0],format=f

; Plot the path of particles, relative to their starting point, and
; compensating for the bulk motion
  for i=ii[0],ii[1] do begin
    p2=rp(i)                                    ; particle data
    s2=sort(p2.i[*,0])                          ; sorted indices
    ;print,p2.i[s2[ind],0],format=f
    plots, (p2.r[s2[ind],2,0]-p1.r[s1[ind],2,0]-dzdi*(i-ii[0]))/ds[2], $
           (p2.r[s2[ind],0,0]-p1.r[s1[ind],0,0])/ds[0],psym=3
  end
END