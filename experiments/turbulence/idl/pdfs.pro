
dir=['','']
run=['','']
dir[0]='C:/cygwin64/home/codes/dispatch2/experiments/turbulence/data/
dir[1]='C:/cygwin64/home/codes/dispatch2/experiments/turbulence/data/
;dir[1]='E:/Paolo/'
run[0]='datahub'
run[1]='download'
dir[0]='C:/cygwin64/home/codes/dispatch2_clean/experiments/turbulence/data/
dir[1]='E:/Paolo/'
run[0]='t_turn=0_2'
run[1]='t_turn=0_1'
nt0=1
nt1=90
nt2=200
dir[0]='C:/cygwin64/home/codes/dispatch2_clean/experiments/turbulence/data/
dir[1]='C:/cygwin64/home/codes/dispatch2_clean/experiments/turbulence/data/
run[0]='datahub'
run[1]='t_turn=0_1'
run[1]='mach4'
nt0=1
nt1=55
nt2=120
dt=0.005
n=96
min=-15.
max=10.
opl=0
pdf=0
yr=[11,0.3*float(n)^3]
!y.style=1
if pdf then yr=yr/yr[1]
for i=0,1 do begin
  ws,i
  bw
  hh=0
  nh=0
  for io=nt0,nt2 do begin
    if io gt nt1 then opl=1 else opl=0
    if i eq 3 then begin
      file=dir[i]+run[i]+'/0000_'+string(io,format='(i4.4)')+'.dat'
    end else begin
      file=dir[i]+run[i]+'/'+string(io,format='(i5.5)')+'/00000.dat'
    end
    open,a,/quiet,file,dim=[n,n,n]
    d=a[0]
    t=io*dt
    !p.title=str(t,format='(f6.3)')
    !p.title=str(io)
    if io ge nt1 then !p.title=run[i]
    plhist,alog(d),/ylog,pdf=pdf,h=h,min=min,max=max,opl=opl,scale=lnd, $
      xr=[min,max],yr=yr
    if io ge nt1 then begin
      hh=hh+h
      nh=nh+1
    end
    wait,0.1
  end
  color
  if i eq 0 then begin
   oplot,lnd,hh/nh,thick=3,color=thecolor('blue')
  end else begin
   oplot,lnd,hh/nh,thick=3,color=thecolor('red')
  end
  ws,2
  bw
  if i eq 0 then begin
   plot,/ylog,yr=yr,lnd,hh/nh,/nodata
   color
   oplot,lnd,hh/nh,thick=3,color=thecolor('blue')
  end else begin
   color
   oplot,lnd,hh/nh,thick=3,color=thecolor('red')
  end
end

END