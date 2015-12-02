
ws,xsize=512,ysize=512
erase
imsz,256
device,decompose=0
loadct,39
n=32
nout=1
nid=4
nx=long(sqrt(nid)+0.1)
!order=0
iv=1

for iout=0,nout-1 do begin
	for id=1,nid do begin
	  name='patch_'+string(id,format='(i3.3)')+'_'+string(iout,format='(i4.4)')
	  file='data/'+name+'.dat'
	  close,1
	  openr,1,file
	  a=assoc(1,fltarr(n,n))
	  ;image,a[0],title=file,id-1
	  iy=(id-1) mod nx
	  ix=(id-1)/nx
	  tv,bytscl(rebin(a[iv],512/nx,512/nx),top=254)+1,ix+(nx-1-iy)*nx
	  ;print,id,min(a[iv],max=max),max
	  sz=1./nx
	  x=sz/2+sz*ix
	  y=sz/2+sz*iy
	  print,file,ix,iy,x,y
    xyouts,x,y,name,align=0.5
	end
  wait,1
end
close,1
END