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
  USE task_mod
  USE patch_mod
  USE mesh_mod
  USE link_mod
  USE task_list_mod
  USE lagrange_mod
  USE experiment_mod
  USE omp_mod
  USE bifrost_scaling_mod
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
  real, dimension(:), allocatable, save:: d_mean, e_mean
  integer, save:: order=2
  integer, save:: verbose=1
  integer, save:: ext=1
  logical, save:: flip_z = .false.
  logical, save:: do_means = .true.
  logical, save:: force_mean = .false.
  real, save:: b0(3)=-999
  type(rs_t), public:: read_snapshot
CONTAINS

!===============================================================================
!===============================================================================
SUBROUTINE init (task_list)
  type(task_list_t), pointer:: task_list
  integer:: iostat
  integer, save:: extension=1
  namelist /snapshot_params/ verbose, meshfile, order, verbose, snapfile, &
            model, b0, flip_z, do_means, force_mean, extension
  !-----------------------------------------------------------------------------
  rewind (io%input)
  read (io%input, snapshot_params, iostat=iostat)
  if (io%master) then
    call io%namelist_warning ('snapshot_params', iostat)
    write (io%namelist, snapshot_params)
    ext=extension
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
  class(task_t), pointer:: task
  class(experiment_t), pointer:: patch
  class(mesh_t), pointer, dimension(:):: m
  integer:: n(3), iv, iiv, recl, ix, iy, iz
  integer:: zi, iext
  integer:: mx, my, mz, nx, ny, nz, nv, ii(8)
  real:: sb(8), sd(8), sc_u, sc_b, minv, maxv
  real(8):: zpos
  real:: ee_top
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn, fl, fldn
  real, dimension(:,:,:), allocatable:: rbuf, buffer
  integer, save:: itimer=0
  integer:: iostat
  logical:: searching
  !-----------------------------------------------------------------------------
  call trace%begin('read_bifrost_snap', itimer=itimer)
  call io%header('begin ic_mod%read_snapshot: Reading Bifrost legacy data')
  call bifrost_scaling%init
  !-----------------------------------------------------------------------------
  ! Pick up patch info from task_list
  !-----------------------------------------------------------------------------
  searching=.true.
  link => task_list%head
  do while (associated(link).and.searching)
    task => link%task
    select type(task)
    class is(experiment_t)
      patch => task
      searching=.false.
    end select
  end do
  m => patch%mesh
  nv = 8
  n = nint(m%b/m%d)
  !---------------------------------------------------------------------------
  ! Read the mesh data from the mesh file
  !---------------------------------------------------------------------------
  open (io_unit%datain, file=trim(meshfile),status='old', form='formatted')
  ! X direction
  read(io_unit%datain, *, iostat=iostat) mx
   if (io%master) print *, 'mx',mx
  allocate( dxm(mx), dxmdn(mx), xm(mx+ext), xmdn(mx+ext))
  read(io_unit%datain, *, iostat=iostat) (xm(ix),ix=1,mx), (xmdn(ix),ix=1,mx), dxm, dxmdn
  
  ! Y direction
  read(io_unit%datain, *, iostat=iostat) my
   if (io%master) print *, 'my',my
  allocate( dym(my), dymdn(my), ym(my+ext), ymdn(my+ext))
  read(io_unit%datain, *, iostat=iostat) (ym(iy),iy=1,my), (ymdn(iy),iy=1,my), dym, dymdn
  
  ! Z direction
  read(io_unit%datain, *, iostat=iostat) mz
   if (io%master) print *, 'mz',mz
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
  if (ext>0) then
    do iext=1,ext
      xm(mx+iext) = xm(mx+iext-1) + (xm(mx+iext-1)-xm(mx+iext-2))
      ym(my+iext) = ym(my+iext-1) + (ym(my+iext-1)-ym(my+iext-2))
      xmdn(mx+iext) = xmdn(mx+iext-1) + (xmdn(mx+iext-1)-xmdn(mx+iext-2))
      ymdn(my+iext) = ymdn(my+iext-1) + (ymdn(my+iext-1)-ymdn(my+iext-2))
    end do
  end if
  if (flip_z) then
    allocate(fl(mz), fldn(mz))
    fl = zm
    fldn = zmdn
    do iz=1,mz
      zi = mz - iz + 1
      zm(iz) = -fl(zi)
      zmdn(iz) = -fldn(zi)
    end do
    deallocate(fl, fldn)
  end if
  
  if (verbose>1) then
    do iz=1,size(dzm)
       if (io%master) print *, 'initial_t%stagger: iz, min(dx),min(dy), dz', &
        iy, minval(dxm), minval(dym), dzm(iz)
    end do
  end if
  if (io%master) then
    print *, 'min(dx),min(dy),min(dz)', minval(dxm), minval(dym), minval(dzm)
    print *, mz, zm(1), size(dzm)
    print *, 'nz for max res:', int((zm(mz)-zm(1))/minval(dzm))+1
  end if
  close (io_unit%datain)

  if (do_means.and..not.allocated(d_mean)) then
    allocate (d_mean(mz))
    allocate (e_mean(mz))
  end if
  
  !---------------------------------------------------------------------------
  ! Allocate a buffer for reading Bifrost legacy snapshots
  !---------------------------------------------------------------------------
  allocate (rbuf(mx,my,mz))
  if (allocated(rbuf))call io%bits_mem (storage_size (rbuf), mx*my*mz)
  allocate (buffer(mx+ext,my+ext,mz))
  call io%bits_mem (storage_size (buffer), product(shape(buffer)), 'buffer')
  recl = kind(1.0)*mx*my*mz
  if (io%master) &
    print *,'reading MHD data from ',trim(snapfile)
  open (io_unit%datain, file=trim(snapfile), form='unformatted', access='direct', recl=recl)
  ! Indices
  ii=[patch%idx%d,patch%idx%px,patch%idx%py,patch%idx%pz, &
      patch%idx%e,patch%idx%bx,patch%idx%by,patch%idx%bz]
  ! Bifrost scaling
  sb=[bifrost_scaling%d, bifrost_scaling%u*bifrost_scaling%d, &
      bifrost_scaling%u*bifrost_scaling%d, bifrost_scaling%u*bifrost_scaling%d, &
      bifrost_scaling%e, bifrost_scaling%b, bifrost_scaling%b, bifrost_scaling%b]
  ! Dispatch scaling
  sd=[scaling%d, scaling%u*scaling%d, scaling%u*scaling%d, scaling%u*scaling%d, &
      scaling%e, scaling%b, scaling%b, scaling%b]
  !-----------------------------------------------------------------------------
  ! Read data and optionally swap z
  !-----------------------------------------------------------------------------
  do iv=1,min(nv,8)
    iiv = ii(iv)
    if (io%master) print *, 'reading variable index iv =', iv
    read (io_unit%datain, rec=iv) rbuf(1:mx,1:my,1:mz)
              
    if (flip_z) then
      do iz=1,mz
        zi = mz -iz + 1
        do iy=1,my
          do ix=1,mx
            buffer(ix,iy,iz) = rbuf(ix,iy,zi) * sb(iv) / sd(iv)
          end do
        end do
      end do
      if ((iiv==patch%idx%pz).or.(iiv==patch%idx%py).or. &
        (iiv==patch%idx%bz).or.(iiv==patch%idx%by)) then
        if (io%master) print  *, 'flipping buffer..'
        do iz=1,mz
          do iy=1,my
            do ix=1,mx
              buffer(ix,iy,iz) = -buffer(ix,iy,iz)
            end do
          end do
        end do
      end if
    else
      do iz=1,mz; do iy=1,my; do ix=1,mx
      buffer(ix,iy,iz) = rbuf(ix,iy,iz) * sb(iv) / sd(iv)
      end do; end do; end do
    end if
      minv = minval(buffer(1:mx,1:my,1:mz))
      maxv = maxval(buffer(1:mx,1:my,1:mz))
       if (io%master) then
      if ((minv>0.0).and.(maxv>0.0)) then
        print *, 'iv minmax log minmax in snapshot', minv, maxv, log(minv), log(maxv)
      else
        print *, 'iv minmax log minmax in snapshot', minv, maxv
      end if
      print *, 'iv minmax location in snapshot', minloc(buffer(1:mx,1:my,1:mz)), maxloc(buffer(1:mx,1:my,1:mz))
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
      do iext=1,ext
        buffer(mx+iext,:,:) = buffer(iext,:,:)
      end do
    else
      do iext=1,ext
        buffer(mx+iext,:,:) = buffer(mx+iext-1,:,:)
      end do
    end if
    ! for Y
    if (patch%periodic(2)) then
      do iext=1,ext
        buffer(:,my+iext,:) = buffer(:,iext,:)
      end do
    else
      do iext=1,ext
        buffer(:,my+iext,:) = buffer(:,my+iext-1,:)
      end do
    end if
    ! Z is not expected to be periodic, for now.
    if (patch%periodic(3)) then
      call io%abort("read_snapshot_mod::read_bifrost: setup is periodic in Z. That not implemented or tested.")
    end if
    ! --------------------------------------------------------------------------
    ! Estimate layer mean values, if needed
    ! --------------------------------------------------------------------------
    if (do_means) then
      if (iiv==patch%idx%d) then
        do iz=1,mz
          d_mean(iz) = sum(buffer(:mx,:my,iz))/(mx*my)
        end do
      end if
      if (iiv==patch%idx%e) then
        do iz=1,mz
          e_mean(iz) = sum(buffer(:mx,:my,iz))/(mx*my)
        end do
      end if
    end if
    !---------------------------------------------------------------------------
    ! Pick out the piece that belongs to each patch, based on its position
    !---------------------------------------------------------------------------
    !$omp parallel default (shared) private(link, patch)
    !$omp single
    link => task_list%head
    do while (associated(link))
      task => link%task
      select type(task)
      class is(experiment_t)
        patch => task
        if (verbose==4) print *, &
          'read_bifrost_snap spawning task for patch%id =', patch%id
        !-------------------------------------------------------------------------
        ! Spawn OMP tasks to interpolate the patch data
        !-------------------------------------------------------------------------
        !$omp task default(shared) firstprivate(patch,iiv)
        if (force_mean) then
          call interpolate_mean &
            (iiv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
        else
          call interpolate &
            (iiv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
        end if
         !$omp end task
      end select
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
  !-----------------------------------------------------------------------------
  ! Finalize snapshots. If needed, do the following:
  ! 1) Store mean values of a given layer
  ! 2) Reset the vertical magnetic field
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type(task)
    class is(experiment_t)
      patch => task
      if (do_means) then
        m => patch%mesh
        if (m(3)%lower_boundary) then
          zpos = m(3)%r(m(3)%lo)+m(3)%p
          patch%boundaries%d_bot = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(d_mean,8))
          patch%boundaries%e_bot = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(e_mean,8))
        end if
        if (m(3)%upper_boundary) then
          zpos = m(3)%r(m(3)%ui)+m(3)%p 
          ee_top = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(e_mean,8)) &
                 / lagrange%interpolate1d8(zpos,mz,real(zm,8),real(d_mean,8))
          patch%newton%ee_top = ee_top
          patch%boundaries%ee_top = ee_top
        else
          ee_top = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(e_mean,8)) &
                 / lagrange%interpolate1d8(zpos,mz,real(zm,8),real(d_mean,8))
          patch%newton%ee_top = ee_top
          patch%boundaries%ee_top = ee_top
        end if
      end if
      if (b0(1) >= -99) patch%mem(:,:,:,patch%idx%bx,1,1) = b0(1)
      if (b0(2) >= -99) patch%mem(:,:,:,patch%idx%by,1,1) = b0(2)
      if (b0(3) >= -99) patch%mem(:,:,:,patch%idx%bz,1,1) = b0(3)
      if (patch%kind(1:6)=='ramses') call patch%to_total_e(patch%it)
    end select
    link => link%next
  end do
  deallocate (xm,ym,zm,xmdn,ymdn,zmdn)
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
  class(task_t), pointer:: task
  class(experiment_t), pointer:: patch
  class(mesh_t), pointer, dimension(:):: m
  integer:: n(3), iv, iiv, recl, ix, iy, iz
  integer:: zi, iext
  integer:: mx, my, mz, nx, ny, nz, nv, ii(9), biv
  real(8):: zpos
  real   :: ee_top
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn, fl, fldn
  real, dimension(:,:,:), allocatable:: rbuf, buffer
  logical:: searching
  integer, save:: itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_snap', itimer=itimer)
  call io%header('begin ic_mod%read_stagger: Reading legacy Stagger Code data')
  !-----------------------------------------------------------------------------
  ! Pick up patch info from task_list
  !-----------------------------------------------------------------------------
  searching=.true.
  link => task_list%head
  do while (associated(link).and.searching)
    task => link%task
    select type(task)
    class is(experiment_t)
      patch => task
      searching=.false.
    end select
  end do
  m => patch%mesh
  nv = patch%nv
  n = nint(m%b/m%d)
  !---------------------------------------------------------------------------
  ! Read the mesh meta-data from the .msh file
  !---------------------------------------------------------------------------
  open (io_unit%datain, file=trim(meshfile), form='unformatted')
  read(io_unit%datain) my, mz, mx
  if (io%master) then
    print *,'reading mesh data from ',trim(meshfile)
    print *, 'setup dimensions =', n, nv
    print *, ' file dims =', mx, my, mz
    print *, ' mesh dims =', m%n, m%gn
    print *, 'setup dims =', nint(m%b/m%d), m%b
    print *, 'periodic extension =', ext
  end if
  ! allocate mesh arrays
  allocate( dxm(mx), dxmdn(mx), xm(mx+ext), xmdn(mx+ext))
  allocate( dym(my), dymdn(my), ym(my+ext), ymdn(my+ext))
  allocate( dzm(mz), dzmdn(mz), zm(mz  ), zmdn(mz  ))
  allocate( fl(mz  ), fldn(mz  ))
  ! read meash properties, in order y,z,x -- to make z vertical
  read(io_unit%datain) dym, dymdn, (ym(iy),iy=1,my), (ymdn(iy),iy=1,my)
  read(io_unit%datain) dzm, dzmdn, (zm(iz),iz=1,mz), (zmdn(iz),iz=1,mz)
  read(io_unit%datain) dxm, dxmdn, (xm(ix),ix=1,mx), (xmdn(ix),ix=1,mx)
  ! periodic extensions
  if (ext>0) then
    do iext=1,ext
      xm(mx+iext) = xm(mx+iext-1) + (xm(mx+iext-1)-xm(mx+iext-2))
      ym(my+iext) = ym(my+iext-1) + (ym(my+iext-1)-ym(my+iext-2))
      xmdn(mx+iext) = xmdn(mx+iext-1) + (xmdn(mx+iext-1)-xmdn(mx+iext-2))
      ymdn(my+iext) = ymdn(my+iext-1) + (ymdn(my+iext-1)-ymdn(my+iext-2))
    end do
  end if
  if (flip_z) then
    fl = zm
    fldn = zmdn
    do iz=1,mz
      zi = mz - iz + 1
      zm(iz) = -fl(zi)
      zmdn(iz) = -fldn(zi)
    end do
  end if
  if (do_means.and..not.allocated(d_mean)) then
    allocate (d_mean(mz))
    allocate (e_mean(mz))
  end if
  if (verbose>1) then
    do iz=1,size(dzm)
      print *, 'initial_t%stagger: iz, zm, zmdn, min(dx),min(dy), dz', &
        iz, zm(iz), zmdn(iz),minval(dxm), minval(dym), dzm(iz)
    end do
  end if
  if (io%master) then
    print *, 'min(dx),min(dy),min(dz)', minval(dxm), minval(dym), minval(dzm)
    print *, 'nz for max res:', int((zm(1)-zm(mz))/minval(dzm))+1
  end if
  close (io_unit%datain)
  !---------------------------------------------------------------------------
  ! Allocate buffers for reading and for transposing to z vertical
  !---------------------------------------------------------------------------
  allocate (rbuf(my,mz,mx))
  call io%bits_mem (storage_size (rbuf), mx*my*mz)
  allocate (buffer(mx+ext,my+ext,mz))
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
    ! --------------------------------------------------------------------------
    ! Flip the Z and Y axis
    ! --------------------------------------------------------------------------
    if (flip_z) then
      do iz=1,mz; do iy=1,my; do ix=1,mx
        zi = mz - iz + 1
        buffer(ix,iy,iz) = rbuf(iy,zi,ix)
      end do; end do; end do
      if ((iiv==patch%idx%pz).or.(iiv==patch%idx%py).or. &
        (iiv==patch%idx%bz).or.(iiv==patch%idx%by)) then
        !do iz=1,mz; do iy=1,my; do ix=1,mx
        if (io%master) print  *, 'flipping buffer..'
        buffer = -buffer
        !end do; end do; end do
      end if
    else
      do iz=1,mz; do iy=1,my; do ix=1,mx
        buffer(ix,iy,iz) = rbuf(iy,iz,ix)
      end do; end do; end do
    end if
    !
    if (verbose > 0) &
      print *, biv, minval(rbuf), maxval(rbuf)
    !---------------------------------------------------------------------------
    ! Add periodic extension
    !---------------------------------------------------------------------------
    ! for X
    if (patch%periodic(1)) then
      do iext=1,ext
        buffer(mx+iext,:,:) = buffer(iext,:,:)
      end do
    else
      do iext=1,ext
        buffer(mx+iext,:,:) = buffer(mx+iext-1,:,:)
      end do
    end if
    ! for Y
    if (patch%periodic(2)) then
      do iext=1,ext
        buffer(:,my+iext,:) = buffer(:,iext,:)
      end do
    else
      do iext=1,ext
        buffer(:,my+iext,:) = buffer(:,my+iext-1,:)
      end do
    end if
    if (verbose > 0) &
      print *, biv, minval(buffer), maxval(buffer)
    ! --------------------------------------------------------------------------
    ! Estimate layer mean values, if needed
    ! --------------------------------------------------------------------------
    if (do_means) then
      if (iiv==patch%idx%d) then
        do iz=1,mz
          d_mean(iz) = sum(buffer(:mx,:my,iz))/(mx*my)
        end do
      end if
      if (iiv==patch%idx%e) then
        do iz=1,mz
          e_mean(iz) = sum(buffer(:mx,:my,iz))/(mx*my)
        end do
      end if
    end if
    !---------------------------------------------------------------------------
    ! Pick out the piece that belongs to each patch, based on its position
    !---------------------------------------------------------------------------
    !$omp parallel default (shared) private(link, patch)
    !$omp single
    link => task_list%head
    do while (associated(link))
      task => link%task
      select type(task)
      class is(experiment_t)
        patch => task
        if (verbose==4) print *, 'read_stagger_snap spawning task for patch%id =', patch%id
        !-------------------------------------------------------------------------
        ! Spawn OMP tasks to interpolate the patch data
        !-------------------------------------------------------------------------
        !$omp task default(shared) firstprivate(patch,iiv)
         if (force_mean) then
          call interpolate_mean (iiv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
        else
          call interpolate (iiv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
        end if
        !$omp end task
      end select
      link => link%next
    end do
    !$omp end single
    !$omp end parallel
    biv = biv+1
  end do
  close (io_unit%datain)
  if (allocated(buffer)) then
    call io%bits_mem (-storage_size (buffer), product(shape(buffer)), 'buffer')
    deallocate (buffer)
  end if
  if (allocated(rbuf)) then
    call io%bits_mem (-storage_size (rbuf), mx*my*mz)
    deallocate(rbuf)
  end if
  !-----------------------------------------------------------------------------
  ! Finalize snapshots. If needed, do the following:
  ! 1) Store mean values of a given layer
  ! 2) Reset the vertical magnetic field
  !-----------------------------------------------------------------------------
  link => task_list%head
  do while (associated(link))
    task => link%task
    select type(task)
    class is(experiment_t)
      patch => task
      if (do_means) then
        m => patch%mesh
        !if (m(3)%lower_boundary) then
          zpos = m(3)%r(m(3)%lo)+m(3)%p
          !zpos = patch%origin(3)+ 0.5*m(3)%d
          patch%boundaries%d_bot = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(d_mean,8),order=1)
          patch%boundaries%e_bot = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(e_mean,8),order=1)
        !end if
        if (m(3)%upper_boundary) then
          !zpos = m(3)%r(m(3)%uo)+m(3)%p
          zpos = m(3)%r(m(3)%ui)+m(3)%p 
          ee_top = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(e_mean,8),order=1) &
                 / lagrange%interpolate1d8(zpos,mz,real(zm,8),real(d_mean,8),order=1)
          patch%newton%ee_top = ee_top
          patch%boundaries%ee_top = ee_top
        else
          ee_top = lagrange%interpolate1d8(zpos,mz,real(zm,8),real(e_mean,8),order=1) &
                 / lagrange%interpolate1d8(zpos,mz,real(zm,8),real(d_mean,8),order=1)
          patch%newton%ee_top = ee_top
          patch%boundaries%ee_top = ee_top
        end if
      end if
      !
      if (b0(1) >= -99) patch%mem(:,:,:,patch%idx%bx,1,1) = b0(1)
      if (b0(2) >= -99) patch%mem(:,:,:,patch%idx%by,1,1) = b0(2)
      if (b0(3) >= -99) patch%mem(:,:,:,patch%idx%bz,1,1) = b0(3)
      if (patch%kind(1:6)=='ramses') call patch%to_total_e(patch%it)
    end select
    link => link%next
  end do
  deallocate (xm,ym,zm,xmdn,ymdn,zmdn, fl, fldn)
  if (io%master) &
    print *,'read_stagger_snap finished'
  call trace%end(itimer)
END SUBROUTINE read_stagger_snap

!===============================================================================
!> Thread safe interpolation work, where only iv, patch, ff changes in the calls
!===============================================================================
SUBROUTINE interpolate (iv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
  class(experiment_t), pointer:: patch
  real, dimension(:,:,:), allocatable:: buffer
  integer:: iv, mx, my, mz
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn
  !.............................................................................
  class(mesh_t), pointer:: m(:), m1, m2, m3
  integer:: nx, ny, nz, o(3)
  real, dimension(:), pointer:: xi, yi, zi, xo, yo, zo
  !-----------------------------------------------------------------------------
  m => patch%mesh
  nx = m(1)%gn
  ny = m(2)%gn
  nz = m(3)%gn
  allocate (xo(nx), yo(ny), zo(nz))
  m1=>m(1)
  m2=>m(2)
  m3=>m(3)
  xo = m1%p + m1%r + m1%h(iv)*m1%d
  yo = m2%p + m2%r + m2%h(iv)*m2%d
  zo = m3%p + m3%r + m3%h(iv)*m3%d
  xi => xm
  yi => ym
  zi => zm
  if (iv==patch%idx%px.or.iv==patch%idx%bx) xi => xmdn
  if (iv==patch%idx%py.or.iv==patch%idx%by) yi => ymdn
  if (iv==patch%idx%pz.or.iv==patch%idx%bz) zi => zmdn
  xo = modulo(xo-xi(1),xm(mx+ext)-xm(1))
  yo = modulo(yo-yi(1),ym(my+ext)-ym(1))
  associate (ff => patch%mem(:,:,:,iv,patch%it,1))
  if (order==-2) then
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

!===============================================================================
!> Thread safe interpolation work, where only iv, patch, ff changes in the calls
!===============================================================================
SUBROUTINE interpolate_mean (iv, patch, buffer, mx, my, mz, xm, ym, zm, xmdn, ymdn, zmdn)
  class(experiment_t), pointer:: patch
  real, dimension(:,:,:), allocatable:: buffer
  integer:: iv, mx, my, mz, iz
  real(8):: zpos
  real, dimension(:), allocatable, target:: dxm, dxmdn, xm, xmdn
  real, dimension(:), allocatable, target:: dym, dymdn, ym, ymdn
  real, dimension(:), allocatable, target:: dzm, dzmdn, zm, zmdn
  !.............................................................................
  class(mesh_t), pointer:: m(:), m1, m2, m3
   real, dimension(:,:,:), pointer :: ff
  integer:: nx, ny, nz, o(3)
  real, dimension(:), pointer:: xi, yi, zi, xo, yo, zo, mean_buff
  !-----------------------------------------------------------------------------
  m => patch%mesh
  nx = m(1)%gn
  ny = m(2)%gn
  nz = m(3)%gn
  allocate (xo(nx), yo(ny), zo(nz), mean_buff(mz))
  m1=>m(1)
  m2=>m(2)
  m3=>m(3)
  xo = m1%p + m1%r + m1%h(iv)*m1%d
  yo = m2%p + m2%r + m2%h(iv)*m2%d
  zo = m3%p + m3%r + m3%h(iv)*m3%d
  xi => xm
  yi => ym
  zi => zm
  if (iv==patch%idx%px.or.iv==patch%idx%bx) xi => xmdn
  if (iv==patch%idx%py.or.iv==patch%idx%by) yi => ymdn
  if (iv==patch%idx%pz.or.iv==patch%idx%bz) zi => zmdn
  xo = modulo(xo-xi(1),xm(mx+ext)-xm(1))
  yo = modulo(yo-yi(1),ym(my+ext)-ym(1))
  ff => patch%mem(:,:,:,iv,patch%it,1)
  if (iv==patch%idx%pz.or.iv==patch%idx%py.or.iv==patch%idx%px) then
    ff(:,:,:) = 0.0
  else
    do iz = 1, mz
      mean_buff(iz) = sum(buffer(:,:,iz))/(mx*my)
    end do
    do iz = 1, nz
      zpos = m3%r(iz)+m3%p + m3%h(iv)*m3%d
      ff(:,:,iz) = lagrange%interpolate1d8(zpos,mz,real(zi,8),real(mean_buff,8))
    end do
  end if
  deallocate (xo, yo, zo, mean_buff)
END SUBROUTINE interpolate_mean

END MODULE read_snapshot_mod
