
ws,xsize=512,ysize=512
erase
imsz,256
device,decompose=0
loadct,39
n=16
nout=10
nid=16
!order=0

for iout=0,nout-1 do begin
	for id=1,nid do begin
	  name='patch_'+string(id,format='(i3.3)')+'_'+string(iout,format='(i4.4)')
	  file='data/'+name+'.dat'
	  close,1
	  openr,1,file
	  a=assoc(1,fltarr(n,n))
	  ;image,a[0],title=file,id-1
	  y=(id-1) mod 4
	  x=(id-1)/4
	  tv,bytscl(rebin(a[0],128,128),top=254,max=2.0,min=1.0)+1,x+(3-y)*4
	  ;print,id,min(a[0],max=max),max
	  x=0.125+0.25*x
	  y=0.125+0.25*y
    xyouts,x,y,name,align=0.5
	end
  wait,1
end
close,1
END