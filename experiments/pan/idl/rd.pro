FUNCTION rd,run,iv=iv,iout=iout,n=n,gn=gn,ix=ix,iy=iy,iz=iz,nv=nv,verbose=verbose,order=order,time=time
  default,iout,0L
  default,n,32L
  default,gn,1024L
  default,iv,0
  default,nv,5L
  default,ix,[0,gn/n-1]
  default,iy,[0,gn/n-1]
  default,iz,[0,gn/n-1]
  default,verbose,0
  default,order,1
  if keyword_set(time) then wc=systime(1)
  if n_elements(ix) eq 1 then ix=[ix,ix]
  if n_elements(iy) eq 1 then iy=[iy,iy]
  if n_elements(iz) eq 1 then iz=[iz,iz]
  file='data/'+run+'/'+string(iout,format='(i5.5)')+'/snapshot.dat'
  file='data/'+run+'/00000/snapshot.dat'
  openr,u,/get,file
  a=assoc(u,fltarr(n,n,n))
  m=long(gn)/n
  np=m^3
  f=fltarr(n*(ix[1]-ix[0]+1),n*(iy[1]-iy[0]+1),n*(iz[1]-iz[0]+1),/nozer)
  help,f
  if verbose gt 0 then begin
    print,np
    print,ix
    print,iy
    print,iz
  end 
  for i=0L,np-1 do begin
    ii=i
    jz=i/m^2
    if jz lt iz[0] then continue
    if jz gt iz[1] then continue
    ii=i-jz*m^2
    jy=ii/m
    if jy lt iy[0] then continue
    if jy gt iy[1] then continue
    jx=ii-jy*m
    if jx lt ix[0] then continue
    if jx gt ix[1] then continue
    kx=n*(jx-ix[0])
    ky=n*(jy-iy[0])
    kz=n*(jz-iz[0])
    nx=n
    ny=n
    nz=n
    if verbose eq 2 then print,i,kx,ky,kz,jx,jy,jz
    if verbose gt 2 then begin
      if order eq 0 then ii=iv+i*nv+iout*nv*np else ii=i+(iv+iout*nv)*np
      ff=a[ii]
      print,i,kx,ky,kz,nx,ny,nz,min(ff,max=max),max
      f[kx:kx+nx-1,ky:ky+ny-1,kz:kz+nx-1]=ff
    end else begin
      if order eq 0 then ii=iv+i*nv+iout*nv*np else ii=i+(iv+iout*nv)*np
      f[kx:kx+nx-1,ky:ky+ny-1,kz:kz+nx-1]=a[ii]
    end
  end
  close,u
  free_lun,u
  if keyword_set(time) then print,systime(1)-wc,' sec'
  return,f
END
