
sz=512
window,xsize=sz,ysize=sz
device,decompose=0
loadct,39
nout=10
iv=1
im=0

for iout=0,nout-1 do begin
	for id=1,1000 do begin
	  name='data/sedov_'+string(id,format='(i3.3)')+'_'+string(iout,format='(i4.4)')
	  print,name
	  close,1 & openr,1,name+'.txt',error=error
	  if error ne 0 then break
	  time=0.0 & readf,1,time
	  m=lonarr(3) & n=m & li=m & ui=m & readf,1,m,n,li,ui
	  if n_elements(im) le 1 then im=fltarr(m[0],m[1])
	  close,1 & openr,1,name+'.dat'
	  a=assoc(1,fltarr(n[0],n[1]))
	  im[li[0]:ui[0],li[1]:ui[1]]=a[iv]
	end
	tvscl,rebin(im,sz,sz),top=254
  xyouts,0.5,0.95,align=0.5,/normal,'time ='+string(time,format='(f5.2)')
  wait,0.1
end
close,1
END