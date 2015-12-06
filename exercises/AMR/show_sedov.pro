
sz=512
ws,xsize=sz,ysize=sz
device,decompose=0
loadct,39
iv=0
im=0
k=1

mmax=0
for iout=0,nout-1 do begin
	for id=1,1000 do begin
	  name='data/sedov_'+string(id,format='(i3.3)')+'_'+string(iout,format='(i4.4)')
	  print,name
	  close,1 & openr,1,name+'.txt',error=error
	  if error ne 0 then break
	  time=0.0 & readf,1,time
	  m=lonarr(3) & n=m & li=m & ui=m & readf,1,m,n,li,ui
	  mmax=(mmax > m)
  end
end
nout=iout+1
print,'max resolution:', mmax, iout

for iout=0,nout-1 do begin
	for id=1,1000 do begin
	  name='data/sedov_'+string(id,format='(i3.3)')+'_'+string(iout,format='(i4.4)')
	  close,1 & openr,1,name+'.txt',error=error
	  if error ne 0 then break
	  ;print,name
	  time=0.0 & readf,1,time
	  m=lonarr(3) & n=m & li=m & ui=m & readf,1,m,n,li,ui
	  k=mmax[0]/m[0]
	  imsz=size(im)
	  if n_elements(im) le 1 then begin
	    im=fltarr(k*m[0],k*m[1])
	    im1=im
	  end
	  close,1 & openr,1,name+'.dat'
	  a=assoc(1,fltarr(n[0],n[1]))
	  new=rebin(a[iv],k*n[0],k*n[1])
	  im[k*li[0]:k*(ui[0]+1)-1,k*li[1]:k*(ui[1]+1)-1]=new
	  b=4.
	  new[*,0]=b & new[0,*]=b & new[*,k*n[1]-1]=b & new[k*n[0]-1,*]=b
	  im1[k*li[0]:k*(ui[0]+1)-1,k*li[1]:k*(ui[1]+1)-1]=new
	  tvscl,rebin(im1,sz,sz),top=254
	  wait,.01
  end
	xyouts,0.5,0.95,align=0.5,/normal,'time ='+string(time,format='(f7.4)')
  print,name,time
  wait,0.1
end
;close,1
END