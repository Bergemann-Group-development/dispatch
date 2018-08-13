!> =============================================================================
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
!> =============================================================================

MODULE IC_mod
  USE io_mod
  USE trace_mod
  USE patch_mod
  USE mesh_mod
  USE link_mod
  USE task_list_mod
  USE lagrange_mod
  USE experiment_mod
  USE omp_mod
  implicit none
  private
  type, public:: ic_t
  contains
    procedure, nopass:: init
  end type
  integer, save:: verbose=1
  integer, save:: order=2
  type(ic_t), public:: ic
CONTAINS

!> =============================================================================
!> Initialize from a legacy raw data file, with auxiliary grid file
!> =============================================================================
SUBROUTINE init (task_list)
  type(task_list_t), pointer:: task_list
  !.............................................................................
  class(link_t), pointer:: link
  class(patch_t), pointer:: patch
  class(mesh_t), pointer, dimension(:):: m
  integer:: n(3), iv, iiv, recl, ix, iy, iz, iostat
  integer:: mx, my, mz, nx, ny, nz, nv, ii(9)
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn
  real, dimension(:,:,:), allocatable:: rbuf, buffer
  character(len=64):: file='mhd63'
  namelist /ic_params/ file, order, verbose
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace_begin('ic_mod%init', itimer=itimer)
  rewind (io%input)
  read (io%input, ic_params, iostat=iostat)
  if (io%master) then
    write (*, ic_params)
  else
    verbose = 0
  end if
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
  open (io_unit%data, file=trim(file)//'.msh', form='unformatted')
  read(io_unit%data) mx, mz, my
  if (io%master) then
    print *,'reading mesh data from ' //trim(file)//'.msh'
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
  read(io_unit%data) dym, dymdn, (ym(iy),iy=1,my), (ymdn(iy),iy=1,my)
  read(io_unit%data) dzm, dzmdn, (zm(iz),iz=1,mz), (zmdn(iz),iz=1,mz)
  read(io_unit%data) dxm, dxmdn, (xm(ix),ix=1,mx), (xmdn(ix),ix=1,mx)
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
    print *, 'nz for max res:', int((zm(my)-zm(1))/minval(dzm))+1
  end if
  close (io_unit%data)
  !---------------------------------------------------------------------------
  ! Allocate buffers for reading and for transposing to z vertical
  !---------------------------------------------------------------------------
  allocate (rbuf(my,mz,mx))
  call io%bits_mem (storage_size (rbuf), mx*my*mz)
  allocate (buffer(mx+1,my+1,mz))
  call io%bits_mem (storage_size (buffer), product(shape(buffer)), 'buffer')
  recl = 4*product(shape(rbuf))     ! FIXME; this could be compiler dependent
  if (io%master) &
    print *,'reading MHD data from ' //trim(file)//'.dat'
  open (io_unit%data, file=trim(file)//'.dat', form='unformatted', access='direct', recl=recl)
  ii=[patch%idx%d,patch%idx%py,patch%idx%pz,patch%idx%px, &
      patch%idx%e,patch%idx%tt,patch%idx%by,patch%idx%bz,patch%idx%bx]
  !-----------------------------------------------------------------------------
  ! Read data and swap indices, so z becomes vertical
  !-----------------------------------------------------------------------------
  do iv=1,min(nv,9)
    iiv = ii(iv)
    if (io%master) &
      print *, 'reading variable index iv =', iv
    read (io_unit%data, rec=iv) rbuf
    do iz=1,mz; do iy=1,my; do ix=1,mx
      buffer(ix,iy,iz) = rbuf(iy,iz,ix)
    end do; end do; end do
    if (verbose > 0) &
      print *, iv, minval(rbuf), maxval(rbuf)
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
      if (verbose==3) print *, 'ic_t%init: spawning task for patch%id =', patch%id
      !-------------------------------------------------------------------------
      ! Spawn OMP tasks to interpolate the patch data
      !-------------------------------------------------------------------------
      !$omp task default(shared) firstprivate(patch,iiv)
        call interpolate (iiv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
        call patch%init_bdries
        select type (patch)
        class is (experiment_t)
          call patch%boundary_conditions
        end select
      !$omp end task
      link => link%next
    end do
    !$omp end single
    !$omp end parallel
  end do
  close (io_unit%data)
  deallocate (rbuf, buffer)
  deallocate (xm,ym,zm,xmdn,ymdn,zmdn)
  if (io%master) &
    print *,'ic_t%init finished'
  call trace_end
END SUBROUTINE init

!> =============================================================================
!> Thread safe interpolation work, where only iv, patch, ff changes in the calls
!> =============================================================================
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
  if (verbose==1) then
    print '(i5,i3,2x,2(2x,a,3(2x,2f7.2)))', patch%id, iv, &
    'xyz(in):', &
    minval(xi), maxval(xi), minval(yi), maxval(yi), minval(zi), maxval(zi), &
    'xyz(out):', &
    minval(xo), maxval(xo), minval(yo), maxval(yo), minval(zo), maxval(zo)
  else if (verbose==2) then
      o = maxloc(ff)
      print '(1x,a,2i3,3i4,2x,3f8.3,1p,2(2x,2e12.4))', &
        'omp,iv,ix,iy,iz,pos,bufmin,bufmax,fmin,fmax =', omp%thread, iv, o, m%p, &
        minval(buffer), maxval(buffer), &
        minval(    ff), maxval(    ff)
  end if
  end associate
  deallocate (xo, yo, zo)
END SUBROUTINE

END MODULE ic_mod
