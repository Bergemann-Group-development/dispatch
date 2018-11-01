;===============================================================================
; Read the RHS value of a namelist variable in a file with NAME = value
;===============================================================================
FUNCTION read_str,file,name
  openr,u,file,/get
  str=''
  while not eof(u) do begin
    readf,u,str
    if stregex(str,strupcase(name)+' *=') ne -1 then begin
      a=strsplit(str,'=',/extract)
      b=strsplit(a[1],',',/extract)
      free_lun,u
      return,strtrim(b[0],2)
    end
  end
  free_lun,u
  return,''
END

; ==============================================================================
; Convert a zero-based integer address to a 3-D address -- optional /mpi order
; ==============================================================================
FUNCTION from_1D_to_3D, i, n, mpi=mpi
  if keyword_set(mpi) then begin
    ii=i
    iz=ii mod n[2]
    ii=(ii-iz)/n[2]
    iy=ii mod n[1]
    ii=(ii-iy)/n[1]
    ix=ii
  end else begin
    ii=i
    ix=ii mod n[0]
    ii=(ii-ix)/n[0]
    iy=ii mod n[1]
    ii=(ii-iy)/n[1]
    iz=ii
  end
  return, [ix,iy,iz]
END

; ==============================================================================
; Convert a zero-based patch number to a 3-D space address, with MPI-numbering
; ==============================================================================
FUNCTION id_to_cartesian, id, n_mpi, n_patch, verbose=verbose
  n_rank = product(n_mpi,/integer)
  rank = id mod n_rank
  mpi_3d = from_1d_to_3d (rank, n_mpi, /mpi)
  id_1d = (id-rank)/n_rank
  id_3d = from_1d_to_3d(id_1d, n_patch)
  cart_3d = id_3d + mpi_3d*n_patch
  if verbose gt 1 then begin
    print,'   rank:',rank
    print,' mpi_3d:',mpi_3d
    print,'  id_1d:',id_1d
    print,'  id_3d:',id_3d
    print,'cart_3d:',cart_3d
  end
  return, cart_3d
END

; ==============================================================================
; Read a DISPATCH unigrid snapshot in compact format.  Assumes that the patches
; in each rank are collected and written out in sequentail order.
; ==============================================================================
FUNCTION rd_single,u,mpi,iout,iv,nv,gn,n,verbose=verbose
  v=fltarr(gn,gn,gn,/nozero)
  n_patch=[1,1,1]*gn/n
  n_per_rank=n_patch/mpi
  nrank=product(mpi,/integer)
  ntask=product(n_per_rank,/integer)
  a=assoc(u,fltarr(n,n,n,ntask,nrank))
  b=a[iv+nv*iout]
  free_lun,u
  for ir=0,nrank-1 do begin
    i_mpi=from_1d_to_3d(ir,mpi,/mpi)
    for ip=0,ntask-1 do begin
      i_pos=from_1d_to_3d(ip,n_per_rank)
      i_patch=i_mpi*n_per_rank+i_pos
      ix=i_patch[0]*n & jx=ix+n-1
      iy=i_patch[1]*n & jy=iy+n-1
      iz=i_patch[2]*n & jz=iz+n-1
      if verbose gt 1 then print,ir,ip
      v[ix:jx,iy:jy,iz:jz]=b[*,*,*,ip,ir]
    end
  end
  return,v
END

; ==============================================================================
; Read a DISPATCH unigrid snapshot
; ==============================================================================
FUNCTION rd,run,mpi=mpi,iv=iv,iout=iout,nn=n,gn=gn,nv=nv,order=order,time=time, $
            hack=hack,verbose=verbose
  default,mpi,[1,1,1]
  default,iout,0L
  default,n,32L
  default,gn,1024L
  default,iv,0
  default,nv,5L
  default,verbose,0
  default,order,0
  default,hack,0
  if keyword_set(time) then wc=systime(1)
  ;
  ; Open file for reading patch variables
  file='data/'+run+'/snapshots.dat'
  meta='data/'+run+'/'+string(iout,format='(i5.5)')+'/snapshot.nml'
  ioformat=long(read_str(meta,'format'))
  if verbose gt 0 then print,'ioformat:', ioformat
  openr,u,/get,file
  if ioformat gt 9 and ioformat lt 14 then begin
    start=systime(1)
    v=rd_single(u,mpi,iout,iv,nv,gn,n,verbose=verbose)
    if keyword_set(time) then print,systime(1)-start,' sec'
    return,v
  end
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
    cart=id_to_cartesian(id1,n_mpi,n_patch,verbose=verbose)
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
    if keyword_set(hack) and id1 mod nmpi ne 0 then begin
      id=id1+nmpi
    end else begin
      id=id1
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
