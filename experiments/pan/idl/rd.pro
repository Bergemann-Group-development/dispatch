; ==============================================================================
; Convert a zero-based integer address to a 3-D address
; ==============================================================================
FUNCTION from_1D_to_3D, i, n
  ii=i
  ix=ii mod n[0]
  ii=(ii-ix)/n[0]
  iy=ii mod n[1]
  ii=(ii-iy)/n[1]
  iz=ii
  return, [ix,iy,iz]
END

; ==============================================================================
; Convert a zero-based patch number to a 3-D space address, with MPI-numbering
; ==============================================================================
FUNCTION id_to_cartesian, id, n_mpi, n_patch
  n_rank = product(n_mpi,/integer)
  rank = id mod n_rank
  mpi_3d = from_1d_to_3d (rank, n_mpi)
  id_1d = (id-rank)/n_rank
  id_3d = from_1d_to_3d(id_1d, n_patch)
  cart_3d = id_3d + mpi_3d*n_patch
  print,'   rank:',rank
  print,' mpi_3d:',mpi_3d
  print,'  id_1d:',id_1d
  print,'  id_3d:',id_3d
  print,'cart_3d:',cart_3d
  return, cart_3d
END

FUNCTION rd,run,mpi=mpi,iv=iv,iout=iout,nn=n,gn=gn,ix=ix,iy=iy,iz=iz,nv=nv,verbose=verbose,order=order,time=time
  default,mpi,[1,1,1]
  default,iout,0L
  default,n,32L
  default,gn,1024L
  default,iv,0
  default,nv,5L
  default,ix,[0,gn/n-1]
  default,iy,[0,gn/n-1]
  default,iz,[0,gn/n-1]
  default,verbose,0
  default,order,0
  if keyword_set(time) then wc=systime(1)
  if n_elements(ix) eq 1 then ix=[ix,ix]
  if n_elements(iy) eq 1 then iy=[iy,iy]
  if n_elements(iz) eq 1 then iz=[iz,iz]
  ;
  ; Open file for reading patch variables
  file='data/'+run+'/snapshots.dat'
  openr,u,/get,file
  a=assoc(u,fltarr(n,n,n))
  ;
  ; Compute 3D dimensions
  m_cell=[gn,gn,gn]                     ; 3D total cells
  n_cell=[n,n,n]                        ; 3D patch cells
  m_patch=m_cell/n_cell                 ; 3D patch count
  n_mpi=long(mpi)                       ; 3D MPI rank dims
  n_patch=m_patch/n_mpi                 ; 3D patches per MPI rank
  ;
  np=product(m_patch,/integer)
  if verbose gt 0 then begin
    print,' m_cell:',m_cell
    print,' n_cell:',n_cell
    print,'m_patch:',m_patch
    print,'n_patch:',n_patch
    print,'  n_mpi:',n_mpi
  end
  f=fltarr(gn,gn,gn,/nozero)
  ;
  ; Loop over all patches (IDL_id = Fortan_id - 1)
  for id1=0L,np-1 do begin
    cart=id_to_cartesian(id1,n_mpi,n_patch)
    jx=cart[0]
    jy=cart[1]
    jz=cart[2]
    ;
    ; At this point [jx,jy,jz] are the coordinates in patch 3D space
    kx=n*jx
    ky=n*jy
    kz=n*jz
    nx=n
    ny=n
    nz=n
    ;
    ; Read patch and insert values
    nmpi=product(n_mpi,/integer)
    if id1 mod nmpi eq 0 then begin
      id=id1
    end else begin
      id=id1+nmpi
    end
    if order eq 0 then ii=iv+id*nv+iout*nv*np else ii=id+(iv+iout*nv)*np
    ff=a[ii]
    f[kx:kx+nx-1,ky:ky+ny-1,kz:kz+nx-1]=ff
    if verbose eq 2 then print,id,jx,jy,jz,kx,ky,kz
    if verbose gt 2 then print,id,kx,ky,kz,min(ff,max=max),max
  end
  close,u
  free_lun,u
  if keyword_set(time) then print,systime(1)-wc,' sec'
  return,f
END
