
;@stagger_2nd
dir=['','']
run=['','']
dir[0]='C:/cygwin64/home/codes/dispatch2_clean/experiments/turbulence/data/
dir[1]='E:/Paolo/'
run[0]='t_turn=0_2'
nt0=1
nt1=200
dir[0]='C:/cygwin64/home/codes/dispatch2_clean/experiments/turbulence/data/
dir[1]='C:/cygwin64/home/codes/dispatch2/experiments/turbulence/data/
run[1]='datahub'
run[0]='download'
run[0]='t_turn=0_1' ;& dir[1]='E:/Paolo/'
n=96
nt0=1
nt1=93
;
dt=0.005
rms=fltarr(2)
u0=fltarr(3,2)
u1=fltarr(2)
for io=nt0,nt1 do begin
   for i=0,1 do begin
      if i eq 1 then begin
	  	 file=dir[i]+run[i]+'/0000_'+string(io,format='(i4.4)')+'.dat'
      end else begin
	 	 file=dir[i]+run[i]+'/'+string(io,format='(i5.5)')+'/00000.dat'
      end
	  open,a,/quiet,file,dim=[n,n,n]
	  d=a[0]
	  px=a[1]
 	  py=a[2]
	  pz=a[3]
	  lnd=alog(d)
	  ad=aver(d)
	  u0[0,i]=aver(px)/ad
	  u0[1,i]=aver(py)/ad
	  u0[2,i]=aver(pz)/ad
	  ux=px/exp(xdn(lnd))
	  uy=py/exp(ydn(lnd))
	  uz=pz/exp(zdn(lnd))
	  t=dt*io
	  ux=ux-u0[0,i]
	  uy=uy-u0[1,i]
	  uz=uz-u0[2,i]
	  u2=ux[0:31,0:31,0:31]^2+uy[0:31,0:31,0:31]^2+uz[0:31,0:31,0:31]^2
	  u1[i]=sqrt(aver(u2))
	  u1[i]=u1[i]/t
	  rms[i]=sqrt(aver(ux^2+uy^2+uz^2))
   end
   print,io,t,rms,u1,reform(u0,6),format='(i3,f7.3,2x,2f7.3,2x,2f7.1,2(2x,3f6.2))'
end

END