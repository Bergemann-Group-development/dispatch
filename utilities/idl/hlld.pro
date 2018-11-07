;-------------------------------------------------------------------------------
PRO find_speed_fast,qvar,vel_info,gamma
  half=0.5d0
  zero=0.0d0
  d=qvar[0] & P=qvar[1] & u=qvar[2] & A=qvar[3]
  v=qvar[4] & B=qvar[5] & w=qvar[6] & C=qvar[7]
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(max(d2*d2-c2*A*A/d,zero)) )
  vel_info = cf
END

;-------------------------------------------------------------------------------
PRO hlld,qleft,qright,fgdnv,gamma=gamma,dbg=dbg
  one=1d0
  default,gamma,1.4d0
  entho = one/(gamma-one)
  half=0.5d0
  zero=0.0d0

  ; Enforce continuity of normal component
  A=half*(qleft[3]+qright[3])
  sgnm=A/(abs(A)+1d-60)
  qleft[3]=A & qright[3]=A

  ; Left variables
  rl=qleft[0] & Pl=qleft[1] & ul=qleft[2]
  vl=qleft[4] & Bl=qleft[5] & wl=qleft[6] & Cl=qleft[7]
  ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
  emagl = half*(A*A+Bl*Bl+Cl*Cl)
  etotl = Pl*entho+ecinl+emagl
  Ptotl = Pl + emagl
  vdotBl= ul*A+vl*Bl+wl*cl

  fmt='(a,12g14.6)'
  if keyword_set(dbg) then print,form=fmt,'hlld1',rl,vl,ecinl,emagl,etotl,Ptotl,vdotBl

  ; Right variables
  rr=qright[0] & Pr=qright[1] & ur=qright[2]
  vr=qright[4] & Br=qright[5] & wr=qright[6] & Cr=qright[7]
  ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
  emagr = half*(A*A+Br*Br+Cr*Cr)
  etotr = Pr*entho+ecinr+emagr
  Ptotr = Pr + emagr
  vdotBr= ur*A+vr*Br+wr*Cr

  if keyword_set(dbg) then print,form=fmt,'hlld2',rr,vr,ecinr,emagr,etotr,Ptotr,vdotBr

  ; Find the largest eigenvalues in the normal direction to the interface
  find_speed_fast,qleft ,cfastl,gamma
  find_speed_fast,qright,cfastr,gamma

  if keyword_set(dbg) then print,form=fmt,'hlld3', cfastl, cfastr

  ; Compute HLL wave speed
  SL=min([ul,ur])-max([cfastl,cfastr])
  SR=max([ul,ur])+max([cfastl,cfastr])

  lim = 1.0d0

  ; Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)

  ; Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

  ; Left star region variables
  rstarl=rl*(SL-ul)/(SL-ustar)
  estar =rl*(SL-ul)*(SL-ustar)-A^2
  el    =rl*(SL-ul)*(SL-ul   )-A^2
  if keyword_set(dbg) then print,form=fmt,'hlld4', estar,rl,SL,UL,ustar,A

  if (estar eq 0) then begin
     vstarl=vl
     Bstarl=Bl
     wstarl=wl
     Cstarl=Cl
  end else begin
     vstarl=vl-A*Bl*(ustar-ul)/estar
     Bstarl=Bl*el/estar
     wstarl=wl-A*Cl*(ustar-ul)/estar
     Cstarl=Cl*el/estar
  end
  vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
  etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
  sqrrstarl=sqrt(rstarl)
  calfvenl=abs(A)/sqrrstarl
  SAL=ustar-calfvenl

  ; Right star region variables
  rstarr=rr*(SR-ur)/(SR-ustar)
  estar =rr*(SR-ur)*(SR-ustar)-A^2
  er    =rr*(SR-ur)*(SR-ur   )-A^2
  if (estar le zero) then begin
     vstarr=vr
     Bstarr=Br
     wstarr=wr
     Cstarr=Cr
  end else begin
     vstarr=vr-A*Br*(ustar-ur)/estar
     Bstarr=Br*er/estar
     wstarr=wr-A*Cr*(ustar-ur)/estar
     Cstarr=Cr*er/estar
  end
  if keyword_set(dbg) then print,form=fmt,'hlldD',estar,vstarr,vr,A,Br,ustar,ur,estar

  vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
  if keyword_set(dbg) then print,form=fmt,'hlldC',ustar,A,vstarr,Bstarr,wstarr,Cstarr

  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
  if keyword_set(dbg) then print,form=fmt,'hlldB',SR,ur,etotr,Ptotr,Ptotstar,ustar,A,vdotBr,vdotBstarr

  sqrrstarr=sqrt(rstarr)
  calfvenr=abs(A)/sqrrstarr
  SAR=ustar+calfvenr

  ; Double star region variables
  vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
  wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
  Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
  Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
  vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
  etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
  etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)

  if keyword_set(dbg) then print,form=fmt,'hlldA',etotstarstarr,etotstarr,sgnm,sqrrstarr,vdotBstarr,vdotBstarstar

  if keyword_set(dbg) then print,form=fmt,'hlld5', SL,SAL,ustar,SAR,SR

  ; Sample the solution at x/t=0
  if (SL gt 0d0) then begin
     ro=rl
     uo=ul
     vo=vl
     wo=wl
     Bo=Bl
     Co=Cl
     Ptoto=Ptotl
     etoto=etotl
     vdotBo=vdotBl
     emago=emagl
  end else if (SAL gt 0d0) then begin
     ro=rstarl
     uo=ustar
     vo=vstarl
     wo=wstarl
     Bo=Bstarl
     Co=Cstarl
     Ptoto=Ptotstar
     etoto=etotstarl
     vdotBo=vdotBstarl
     emago=half*(A*A+Bo*Bo+Co*Co)
  end else if (ustar gt 0d0) then begin
     ro=rstarl
     uo=ustar
     vo=vstarstar
     wo=wstarstar
     Bo=Bstarstar
     Co=Cstarstar
     Ptoto=Ptotstar
     etoto=etotstarstarl
     vdotBo=vdotBstarstar
     emago=half*(A*A+Bo*Bo+Co*Co)
  end else if (SAR gt 0d0) then begin
     ro=rstarr
     uo=ustar
     vo=vstarstar
     wo=wstarstar
     Bo=Bstarstar
     Co=Cstarstar
     Ptoto=Ptotstar
     etoto=etotstarstarr
     help,etotstarstarr
     vdotBo=vdotBstarstar
     emago=half*(A*A+Bo*Bo+Co*Co)
  end else if (SR gt 0d0) then begin
     ro=rstarr
     uo=ustar
     vo=vstarr
     wo=wstarr
     Bo=Bstarr
     Co=Cstarr
     Ptoto=Ptotstar
     etoto=etotstarr
     vdotBo=vdotBstarr
     emago=half*(A*A+Bo*Bo+Co*Co)
  end else begin
     ro=rr
     uo=ur
     vo=vr
     wo=wr
     Bo=Br
     Co=Cr
     Ptoto=Ptotr
     etoto=etotr
     vdotBo=vdotBr
     emago=emagr
  end

  ; Compute the Godunov flux
  fgdnv = make_array(size=size(qleft))
  fgdnv[0] = ro*uo
  fgdnv[1] = (etoto+Ptoto)*uo-A*vdotBo
  if keyword_set(dbg) then print,form=fmt,'hlld6', fgdnv[1],etoto,Ptoto,uo,A,vdotBo
  fgdnv[2] = ro*uo*uo+Ptoto+emago*(lim-1.0d0)-A*A*lim
  fgdnv[3] = zero
  fgdnv[4] = ro*uo*vo-A*Bo*lim
  fgdnv[5] = Bo*uo-A*vo
  fgdnv[6] = ro*uo*wo-A*Co*lim
  fgdnv[7] = Co*uo-A*wo
  nvar = n_elements(qleft)
  for ivar=9,nvar do begin
     if (fgdnv[0] gt 0) then begin
        fgdnv[ivar-1] = fgdnv[0]*qleft [ivar-1]
     end else begin
        fgdnv[ivar-1] = fgdnv[0]*qright[ivar-1]
     end
  end

END
