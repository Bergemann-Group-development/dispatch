!===============================================================================
!> IC module for reading in legacy Stagger snapshots into existing patches.
!>
!> The purpose with this separate module is to avoid having to load the entire
!> snapshot on each MPI rank --  there isn't enough memory to do that on the
!> largest experiments from the archive.  In this module one can loop over either
!> variable index iv, or even iv and depth index, reading only one xy-slice at a
!> time.
!>
!> To save startup timme, separate OMP tasks are spawned for each patch, to work
!> in parallel on the interpolation of the ICs to the new patch mesh.
!>
!> The separate reader also saves memory because it operates before the RT solver
!> has started to run, allocating additional scratch memory.
!===============================================================================

MODULE read_snapshot_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE mesh_mod
  USE link_mod
  USE task_list_mod
  USE lagrange_mod
  USE experiment_mod
  USE omp_mod
  USE scaling_mod
  implicit none
  private
  type, public:: rs_t
  contains
    procedure, nopass:: init
  end type
  character(len=64), save:: meshfile='mhd63.msh'
  character(len=64), save:: snapfile='mhd63.dat'
  character(len=64), save:: model='bifrost'
  integer, save:: order=2
  integer, save:: verbose=1
  logical, save:: flip_z = .false.
  real, save:: b0=0.0
  type(rs_t), public:: read_snapshot
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (task_list)
  type(task_list_t), pointer:: task_list
  integer:: iostat
  namelist /snapshot_params/ verbose, meshfile, order, verbose, snapfile, model, b0, flip_z
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, snapshot_params, iostat=iostat)
  if (io%master) then
    write (*, snapshot_params)
  else
    verbose = 0
  end if
  if (io%restart < 0) then 
    select case(trim(model))
      case('bifrost')
        call read_bifrost_snap(task_list)
      case('stagger')
        call read_stagger_snap(task_list)
      case default
        call io%abort('read_snapshot_mod::init: unknown input model')
    end select
  end if
END SUBROUTINE init

!===============================================================================
!> Read a legacy Bifrost raw data file, with auxiliary grid file
!===============================================================================
SUBROUTINE read_bifrost_snap (task_list)
  USE kinds_mod
  type(task_list_t), pointer:: task_list
  !.............................................................................
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  class(mesh_t), pointer, dimension(:):: m
  integer:: n(3), iv, iiv, recl, ix, iy, iz
  integer:: zi
  integer:: mx, my, mz, nx, ny, nz, nv, ii(8)
  real:: sc_u, sc_b
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn
  real, dimension(:,:,:), allocatable:: rbuf, buffer
  integer, save:: itimer=0
  integer:: iostat
  !-----------------------------------------------------------------------------
  call trace%begin('read_bifrost_snap', itimer=itimer)
  call io%header('begin ic_mod%read_snapshot: Reading Bifrost legacy data')
  !-----------------------------------------------------------------------------
  ! Pick up patch info from task_list
  !-----------------------------------------------------------------------------
  patch => task2patch(task_list%head%task)
  m => patch%mesh
  nv = 8
  n = nint(m%b/m%d)
  !---------------------------------------------------------------------------
  ! Read the mesh data from the mesh file
  !---------------------------------------------------------------------------
  open (io_unit%datain, file=trim(meshfile),status='old', form='formatted')
  ! X direction
  read(io_unit%datain, *, iostat=iostat) mx
  print *, 'mx',mx
  allocate( dxm(mx), dxmdn(mx), xm(mx+1), xmdn(mx+1))
  read(io_unit%datain, *, iostat=iostat) (xm(ix),ix=1,mx), (xmdn(ix),ix=1,mx), dxm, dxmdn
  
  ! Y direction
  read(io_unit%datain, *, iostat=iostat) my
  print *, 'my',my
  allocate( dym(my), dymdn(my), ym(my+1), ymdn(my+1))
  read(io_unit%datain, *, iostat=iostat) (ym(iy),iy=1,my), (ymdn(iy),iy=1,my), dym, dymdn
  
  ! Z direction
  read(io_unit%datain, *, iostat=iostat) mz
  print *, 'mz',mz
  allocate( dzm(mz), dzmdn(mz), zm(mz  ), zmdn(mz  ))
  read(io_unit%datain, *, iostat=iostat) (zm(iz),iz=1,mz), (zmdn(iz),iz=1,mz), dzm, dzmdn
  
  if (io%master) then
    print *,'reading mesh data from ' //trim(meshfile)
    print *, 'setup dimensions =', n, nv
    print *, ' file dims =', mx, my, mz
    print *, ' mesh dims =', m%n, m%gn
    print *, 'setup dims =', nint(m%b/m%d), m%b
  end if

  ! periodic extensions
  xm(mx+1) = xm(mx) + (xm(mx)-xm(mx-1))
  ym(my+1) = ym(my) + (ym(my)-ym(my-1))
  xmdn(mx+1) = xmdn(mx) + (xmdn(mx)-xmdn(mx-1))
  ymdn(my+1) = ymdn(my) + (ymdn(my)-ymdn(my-1))
  
  if (verbose>1) then
    do iz=1,size(dzm)
      print *, 'initial_t%stagger: iz, min(dx),min(dy), dz', &
        iy, minval(dxm), minval(dym), dzm(iz)
    end do
  end if
  if (io%master) then
    print *, 'min(dx),min(dy),min(dz)', minval(dxm), minval(dym), minval(dzm)
    print *, mz, zm(1), size(dzm)
    print *, 'nz for max res:', int((zm(mz)-zm(1))/minval(dzm))+1
  end if
  close (io_unit%datain)
  
  !---------------------------------------------------------------------------
  ! Allocate a buffer for reading Bifrost legacy snapshots
  !---------------------------------------------------------------------------
  allocate (rbuf(mx,my,mz))
  if (allocated(rbuf))call io%bits_mem (storage_size (rbuf), mx*my*mz)
  allocate (buffer(mx+1,my+1,mz))
  call io%bits_mem (storage_size (buffer), product(shape(buffer)), 'buffer')
  recl = kind(1.0)*mx*my*mz
  if (io%master) &
    print *,'reading MHD data from ',trim(snapfile)
  open (io_unit%datain, file=trim(snapfile), form='unformatted', access='direct', recl=recl)
  ! Indices
  ii=[patch%idx%d,patch%idx%px,patch%idx%py,patch%idx%pz, &
      patch%idx%e,patch%idx%bx,patch%idx%by,patch%idx%bz]
  !-----------------------------------------------------------------------------
  ! Read data and optionally swap z
  !-----------------------------------------------------------------------------
  do iv=1,min(nv,8)
    iiv = ii(iv)
    if (io%master) print *, 'reading variable index iv =', iv
    read (io_unit%datain, rec=iv) rbuf(1:mx,1:my,1:mz)
    
      print *, 'iv minmax in snapshot', minval(rbuf(1:mx,1:my,1:mz)), maxval(rbuf(1:mx,1:my,1:mz))
      print *, 'iv minmax location in snapshot', minloc(rbuf(1:mx,1:my,1:mz)), maxloc(rbuf(1:mx,1:my,1:mz))
        
    if (flip_z) then
      do iz=1,mz; do iy=1,my; do ix=1,mx
      zi = mz -iz + 1
      buffer(ix,iy,iz) = rbuf(ix,iy,zi)
      end do; end do; end do
    else
      do iz=1,mz; do iy=1,my; do ix=1,mx
      buffer(ix,iy,iz) = rbuf(ix,iy,iz)
      end do; end do; end do
    end if
    if (verbose == 1) &
      print *, iv, minval(buffer(1:mx,1:my,1:mz)), & 
                   minloc(buffer(1:mx,1:my,1:mz)),'||', &
                   maxval(buffer(1:mx,1:my,1:mz)), & 
                   maxloc(buffer(1:mx,1:my,1:mz))
    
    !---------------------------------------------------------------------------
    ! Add periodic extension, if needed
    !---------------------------------------------------------------------------
    ! for X
    if (patch%periodic(1)) then
      buffer(mx+1,:,:) = buffer(1,:,:)
    else
      buffer(mx+1,:,:) = buffer(mx,:,:)
    end if
    ! for Y
    if (patch%periodic(2)) then
      buffer(:,my+1,:) = buffer(:,1,:)
    else
      buffer(:,my+1,:) = buffer(:,my,:)
    end if
    ! Z is not expected to be periodic, for now.
    if (patch%periodic(3)) then
      call io%abort("read_snapshot_mod::read_bifrost: setup is periodic in Z. That not implemented or tested.")
    end if
    !---------------------------------------------------------------------------
    ! Pick out the piece that belongs to each patch, based on its position
    !---------------------------------------------------------------------------
    !$omp parallel default (shared) private(link, patch)
    !$omp single
    link => task_list%head
    do while (associated(link))
      patch => task2patch(link%task)
      if (verbose==4) print *, 'read_bifrost_snap spawning task for patch%id =', patch%id
      !-------------------------------------------------------------------------
      ! Spawn OMP tasks to interpolate the patch data
      !-------------------------------------------------------------------------
      !$omp task default(shared) firstprivate(patch,iiv)
        call interpolate (iiv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
       !$omp end task
      link => link%next
    end do
    !$omp end single
    !$omp end parallel
  end do
  close (io_unit%datain)
  call io%bits_mem (-storage_size (buffer), product(shape(buffer)), 'buffer')
  deallocate (buffer)
  if (allocated(rbuf))call io%bits_mem (-storage_size (rbuf), mx*my*mz)
  if (allocated(rbuf)) deallocate(rbuf)
  deallocate (xm,ym,zm,xmdn,ymdn,zmdn)
  !-----------------------------------------------------------------------------
  ! Reset the vertical magnetic field
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    patch => task2patch(link%task)
    if (b0 >= 0.0) &
      patch%mem(:,:,:,patch%idx%bz,1,1) = b0
    link => link%next
  end do
  if (io%master) &
    print *,'read_bifrost_snap finished'
  call trace%end(itimer)
END SUBROUTINE read_bifrost_snap

!===============================================================================
!> Read a legacy Stagger raw data file, with auxiliary grid file
!===============================================================================
SUBROUTINE read_stagger_snap (task_list)
  type(task_list_t), pointer:: task_list
  !.............................................................................
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  class(mesh_t), pointer, dimension(:):: m
  integer:: n(3), iv, iiv, recl, ix, iy, iz
  integer:: zi
  integer:: mx, my, mz, nx, ny, nz, nv, ii(9), biv
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn
  real, dimension(:,:,:), allocatable:: rbuf, buffer
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_snap', itimer=itimer)
  call io%header('begin ic_mod%read_stagger: Reading legacy Stagger Code data')
  !-----------------------------------------------------------------------------
  ! Pick up patch info from task_list
  !-----------------------------------------------------------------------------
  patch => task2patch(task_list%head%task)
  m => patch%mesh
  nv = patch%nv
  n = nint(m%b/m%d)
  n(3) = n(3)+4
  !---------------------------------------------------------------------------
  ! Read the mesh meta-data from the .msh file
  !---------------------------------------------------------------------------
  open (io_unit%datain, file=trim(meshfile), form='unformatted')
  read(io_unit%datain) mx, mz, my
  if (io%master) then
    print *,'reading mesh data from ',trim(meshfile)
    print *, 'setup dimensions =', n, nv
    print *, ' file dims =', mx, mz, my
    print *, ' mesh dims =', m%n, m%gn
    print *, 'setup dims =', nint(m%b/m%d), m%b
  end if
  ! allocate mesh arrays
  allocate( dxm(mx), dxmdn(mx), xm(mx+1), xmdn(mx+1))
  allocate( dym(my), dymdn(my), ym(my+1), ymdn(my+1))
  allocate( dzm(mz), dzmdn(mz), zm(mz  ), zmdn(mz  ))
  ! read meash properties, in order y,z,x -- to make z vertical
  read(io_unit%datain) dym, dymdn, (ym(iy),iy=1,my), (ymdn(iy),iy=1,my)
  read(io_unit%datain) dzm, dzmdn, (zm(iz),iz=1,mz), (zmdn(iz),iz=1,mz)
  read(io_unit%datain) dxm, dxmdn, (xm(ix),ix=1,mx), (xmdn(ix),ix=1,mx)
  ! periodic extensions
  xm(mx+1) = xm(mx) + (xm(mx)-xm(mx-1))
  ym(my+1) = ym(my) + (ym(my)-ym(my-1))
  xmdn(mx+1) = xmdn(mx) + (xmdn(mx)-xmdn(mx-1))
  ymdn(my+1) = ymdn(my) + (ymdn(my)-ymdn(my-1))
  if (verbose>1) then
    do iz=1,size(dzm)
      print *, 'initial_t%stagger: iz, min(dx),min(dy), dz', &
        iy, minval(dxm), minval(dym), dzm(iz)
    end do
  end if
  if (io%master) then
    print *, 'min(dx),min(dy),min(dz)', minval(dxm), minval(dym), minval(dzm)
    print *, 'nz for max res:', int((zm(mz)-zm(1))/minval(dzm))+1
  end if
  close (io_unit%datain)
  !---------------------------------------------------------------------------
  ! Allocate buffers for reading and for transposing to z vertical
  !---------------------------------------------------------------------------
  allocate (rbuf(my,mz,mx))
  call io%bits_mem (storage_size (rbuf), mx*my*mz)
  allocate (buffer(mx+1,my+1,mz))
  call io%bits_mem (storage_size (buffer), product(shape(buffer)), 'buffer')
  recl = 4*product(shape(rbuf))     ! FIXME; this could be compiler dependent
  if (io%master) &
    print *,'reading MHD data from ',trim(snapfile)
  open (io_unit%datain, file=trim(snapfile), form='unformatted', access='direct', recl=recl)
  !-----------------------------------------------------------------------------
  ! Read data and swap indices, so z becomes vertical
  !-----------------------------------------------------------------------------
  ii = [patch%idx%d, patch%idx%py, patch%idx%pz, patch%idx%px, &
        patch%idx%e, patch%idx%tt, patch%idx%by, patch%idx%bz, patch%idx%bx]
  biv=1
  do iv=1,min(nv,9)
    if ((iv==6).and.(patch%idx%tt<1)) then
      biv = biv+1
      print *, 'skipping temperature', iv
      !cycle
    end if
    iiv = ii(biv)
    if (io%master) &
      print *, 'reading variable index iv =', biv,'stored in', ii(biv)
    read (io_unit%datain, rec=biv) rbuf
    if (flip_z) then
      do iz=1,mz; do iy=1,my; do ix=1,mx
        zi = mz - iz + 1
        buffer(ix,iy,iz) = rbuf(iy,zi,ix)
      end do; end do; end do
    else
      do iz=1,mz; do iy=1,my; do ix=1,mx
        buffer(ix,iy,iz) = rbuf(iy,iz,ix)
      end do; end do; end do
    end if
    if (verbose == 1) &
      print *, biv, minval(rbuf), maxval(rbuf)
    !---------------------------------------------------------------------------
    ! Add periodic extension
    !---------------------------------------------------------------------------
    buffer(mx+1,:,:) = buffer(1,:,:)
    buffer(:,my+1,:) = buffer(:,1,:)
    !---------------------------------------------------------------------------
    ! Pick out the piece that belongs to each patch, based on its position
    !---------------------------------------------------------------------------
    !$omp parallel default (shared) private(link, patch)
    !$omp single
    link => task_list%head
    do while (associated(link))
      patch => task2patch(link%task)
      if (verbose==4) print *, 'read_stagger_snap spawning task for patch%id =', patch%id
      !-------------------------------------------------------------------------
      ! Spawn OMP tasks to interpolate the patch data
      !-------------------------------------------------------------------------
      !$omp task default(shared) firstprivate(patch,iiv)
        call interpolate (iiv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
       !$omp end task
      link => link%next
    end do
    !$omp end single
    !$omp end parallel
    biv = biv+1
  end do
  close (io_unit%datain)
  deallocate (rbuf, buffer)
  deallocate (xm,ym,zm,xmdn,ymdn,zmdn)
  !-----------------------------------------------------------------------------
  ! Reset the vertical magnetic field
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    patch => task2patch(link%task)
    if (b0 >= 0.0) &
      patch%mem(:,:,:,patch%idx%bz,1,1) = b0
    link => link%next
  end do
  if (io%master) &
    print *,'read_stagger_snap finished'
  call trace%end(itimer)
END SUBROUTINE read_stagger_snap

!===============================================================================
!> Thread safe interpolation work, where only iv, patch, ff changes in the calls
!===============================================================================
SUBROUTINE interpolate (iv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
  class(patch_t), pointer:: patch
  real, dimension(:,:,:), allocatable:: buffer
  integer:: iv, mx, my, mz
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn
  !.............................................................................
  class(mesh_t), pointer, dimension(:):: m
  integer:: nx, ny, nz, o(3)
  real, dimension(:), pointer:: xi, yi, zi, xo, yo, zo
  !-----------------------------------------------------------------------------
  m => patch%mesh
  nx = m(1)%gn
  ny = m(2)%gn
  nz = m(3)%gn
  allocate (xo(nx), yo(ny), zo(nz))
  associate (m1=>m(1), m2=>m(2), m3=>m(3))
  xo = m1%p + m1%r + m1%h(iv)*m1%d
  yo = m2%p + m2%r + m2%h(iv)*m2%d
  zo = m3%p + m3%r + m3%h(iv)*m3%d
  end associate
  xi => xm
  yi => ym
  zi => zm
  if (iv==patch%idx%px.or.iv==patch%idx%bx) xi => xmdn
  if (iv==patch%idx%py.or.iv==patch%idx%by) yi => ymdn
  if (iv==patch%idx%pz.or.iv==patch%idx%bz) zi => zmdn
  xo = modulo(xo-xi(1),xm(mx+1)-xm(1))
  yo = modulo(yo-yi(1),ym(my+1)-ym(1))
  associate (ff => patch%mem(:,:,:,iv,patch%it,1))
  if (order==1) then
    call lagrange%trilinear3d (xi, yi, zi, buffer, &
                               xo, yo, zo, ff, order)
  else
    call lagrange%interpolate3d (xi, yi, zi, buffer, &
                                 xo, yo, zo, ff, order)
  end if
  !-----------------------------------------------------------------------
  ! Diagnostic output
  !-----------------------------------------------------------------------
  if (verbose==2) then
    print '(i5,i3,2x,2(2x,a,3(2x,2f7.2)))', patch%id, iv, &
    'xyz(in):', &
    minval(xi), maxval(xi), minval(yi), maxval(yi), minval(zi), maxval(zi), &
    'xyz(out):', &
    minval(xo), maxval(xo), minval(yo), maxval(yo), minval(zo), maxval(zo)
  else if (verbose==3) then
      o = maxloc(ff)
      print '(1x,a,2i3,3i4,2x,3f8.3,1p,2(2x,2e12.4))', &
        'omp,iv,ix,iy,iz,pos,bufmin,bufmax,fmin,fmax =', omp%thread, iv, o, m%p, &
        minval(buffer), maxval(buffer), &
        minval(    ff), maxval(    ff)
  end if
  end associate
  deallocate (xo, yo, zo)
END SUBROUTINE

END MODULE read_snapshot_mod
