from numpy import sqrt, zeros

#-------------------------------------------------------------------------------
def find_speed_fast (qvar, gamma):
  half=0.5
  zero=0.0
  d=qvar[0]; P=qvar[1]; A=qvar[3]; B=qvar[5]; C=qvar[7]
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
  d2 = half*(B2/d+c2)
  cf = sqrt(d2 + sqrt(max(d2*d2-c2*A*A/d,zero)))
  vel_info = cf
  return vel_info

#-------------------------------------------------------------------------------
def hlld (qleft, qright, gamma=1.4, dbg=0):
  one=1.0
  entho=one/(gamma-one)
  half=0.5
  zero=0.0

  # Enforce continuity of normal component
  A=half*(qleft[3]+qright[3])
  sgnm=A/(abs(A)+1e-60)
  qleft[3]=A; qright[3]=A

  # Left variables
  rl=qleft[0]; Pl=qleft[1]; ul=qleft[2]
  vl=qleft[4]; Bl=qleft[5]; wl=qleft[6]
  Cl=qleft[7]
  ecinl = half*(ul*ul+vl*vl+wl*wl)*rl
  emagl = half*(A*A+Bl*Bl+Cl*Cl)
  etotl = Pl*entho+ecinl+emagl
  Ptotl = Pl + emagl
  vdotBl= ul*A+vl*Bl+wl*Cl

  # Right variables
  rr=qright[0]; Pr=qright[1]; ur=qright[2]
  vr=qright[4]; Br=qright[5]; wr=qright[6]; Cr=qright[7]
  ecinr = half*(ur*ur+vr*vr+wr*wr)*rr
  emagr = half*(A*A+Br*Br+Cr*Cr)
  etotr = Pr*entho+ecinr+emagr
  Ptotr = Pr + emagr
  vdotBr= ur*A+vr*Br+wr*Cr

  # Find the largest eigenvalues in the normal direction to the interface
  cfastl = find_speed_fast (qleft, gamma)
  cfastr = find_speed_fast (qright, gamma)

  # Compute HLL wave speed
  SL=min([ul,ur])-max([cfastl,cfastr])
  SR=max([ul,ur])+max([cfastl,cfastr])

  lim = 1.0

  # Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)

  # Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

  # Left star region variables
  rstarl=rl*(SL-ul)/(SL-ustar)
  estar =rl*(SL-ul)*(SL-ustar)-A**2
  el    =rl*(SL-ul)*(SL-ul   )-A**2

  if (estar == 0):
     vstarl=vl
     Bstarl=Bl
     wstarl=wl
     Cstarl=Cl
  else:
     vstarl=vl-A*Bl*(ustar-ul)/estar
     Bstarl=Bl*el/estar
     wstarl=wl-A*Cl*(ustar-ul)/estar
     Cstarl=Cl*el/estar

  vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
  etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
  sqrrstarl=sqrt(rstarl)
  calfvenl=abs(A)/sqrrstarl
  SAL=ustar-calfvenl

  # Right star region variables
  rstarr=rr*(SR-ur)/(SR-ustar)
  estar =rr*(SR-ur)*(SR-ustar)-A**2
  er    =rr*(SR-ur)*(SR-ur   )-A**2
  if (estar <= zero):
     vstarr=vr
     Bstarr=Br
     wstarr=wr
     Cstarr=Cr
  else:
     vstarr=vr-A*Br*(ustar-ur)/estar
     Bstarr=Br*er/estar
     wstarr=wr-A*Cr*(ustar-ur)/estar
     Cstarr=Cr*er/estar


  vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr

  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)

  sqrrstarr=sqrt(rstarr)
  calfvenr=abs(A)/sqrrstarr
  SAR=ustar+calfvenr

  # Double star region variables
  vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
  wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
  Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
  Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
  vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
  etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
  etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)

  # Sample the solution at x/t=0
  if (SL > 0.0):
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
  elif (SAL > 0.0):
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
  elif (ustar > 0.0):
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
  elif (SAR > 0.0):
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
  elif (SR > 0.0):
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
  else:
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


  # Compute the Godunov flux
  fgdnv = zeros(len(qleft))
  fgdnv[0] = ro*uo
  fgdnv[1] = (etoto+Ptoto)*uo-A*vdotBo
  fgdnv[2] = ro*uo*uo+Ptoto+emago*(lim-1.0)-A*A*lim
  fgdnv[3] = zero
  fgdnv[4] = ro*uo*vo-A*Bo*lim
  fgdnv[5] = Bo*uo-A*vo
  fgdnv[6] = ro*uo*wo-A*Co*lim
  fgdnv[7] = Co*uo-A*wo
  nvar = len(qleft)
  if nvar > 8:
      for ivar in (9,nvar):
         if fgdnv[0] > 0:
            fgdnv[ivar] = fgdnv[0]*qleft [ivar]
         else:
            fgdnv[ivar] = fgdnv[0]*qright[ivar]

  return fgdnv
