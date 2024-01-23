MODULE read_stagger_snapshot_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE handyfuncs_mod
  USE op_math_mod, only: linspace, linterpolate_3d
  USE m3d_patch_mod, only: pivot_z, ray_offset
  implicit none
  private

  ! stagger ghost zones
  integer, parameter :: gz=5

  ! stagger coefficients
  real(4), parameter :: uul = 1.0e+8      ! length
  real(4), parameter :: uur = 1.0e-7      ! density
  real(4), parameter :: luul = log(uul)
  real(4), parameter :: luur = log(uur)
  real(4), parameter :: uut = 1.0e2       ! time
  real(4), parameter :: uuv = uul / uut   ! velocity
  real(4), parameter :: uum = uuv * uur   ! momentum
  real(4), parameter :: uue = uuv**2      ! internal energy per unit mass
  real(4), parameter :: uup = uur * uue   ! pressure
  real(4), parameter :: uua = 1./uul      ! extinction coefficient (alpha)
  !character(len=128), public :: stagger_file='/lustre/astro/rhoppe/multi-3d/input/atmos/t5777g44m00/t5777g44m0005'
  character(len=256), save :: mesh_file, snap_file, data_file, aux_file
  integer, save :: dims(3)!=[5,5,230], snap=0
  real(4), save :: maxTemp=huge(1.0)
  real(4), save :: ee_offset=3.2122
  character(len=128), save :: refinement=''
  logical, public, save :: legacy_stagger_interp=.false.!, irregular=.false.
  !integer, dimension(:), allocatable, save, public  :: stagger_levels
  real(8), dimension(:), allocatable, save:: level_ratio
  integer, save :: maxlvl, minlvl
  real(8), save :: norm, mindz
  real(4), dimension(:,:,:), allocatable, save :: iie, iltemp, ilne, ilrho
  real(4), dimension(:,:,:), allocatable, save :: ipx, ipy, ipz

  public :: read_stagger_params
  public :: read_stagger_size
  public :: read_stagger_snapshot
  public :: load_stagger_subregion
  public :: get_stagger_subregion
  public :: deallocate_stagger
    
contains

!===============================================================================
SUBROUTINE read_stagger_params(stagger_file, snap, nx, ny, nz)
  character(len=*), intent(in) :: stagger_file
  integer, intent(in) :: snap, nx, ny, nz
  integer :: iostat, nsub
  logical, save :: first_time=.true.
  namelist /stagger_params/ refinement
  !-----------------------------------------------------------------------------
  if (first_time) then
    first_time = .false.
    rewind (io%input)
    read (io%input, stagger_params, iostat=iostat)
    if (iostat > 0) then
      rewind (io%input)
      read (io%input, stagger_params)
    end if
    call io%namelist_warning ('stagger_params', iostat)
    write (io%namelist, stagger_params)

    dims = [nx, ny, nz]

    mesh_file = trim(stagger_file) // '.msh'
    write(snap_file, '(A,A1,I5.5)') trim(stagger_file), '_', snap
    write(data_file, '(A,A)') trim(snap_file), '.dat'
    write(aux_file,  '(A,A)') trim(snap_file), '.aux'
    irregular = refinement /= ''

    call count_words(refinement, nsub)
    allocate(stagger_levels(nsub))
    allocate(level_ratio(nsub))
    read(refinement, *) stagger_levels
    maxlvl = maxval(stagger_levels)
    minlvl = minval(stagger_levels)

    !level_ratio = 2**(maxlvl-stagger_levels)
    !norm = sum(level_ratio)
    !level_ratio = level_ratio / norm
    call read_stagger_data(ilrho, ipx, ipy, ipz, iie, iltemp, ilne)
  end if
END SUBROUTINE read_stagger_params

!===============================================================================
SUBROUTINE read_stagger_snapshot(ds, d, px, py, pz, ee, tt, ne)
  real(8), intent(out) :: ds(3)
  real(4), pointer, contiguous, dimension(:,:,:), intent(out) :: d, px, py, pz
  real(4), pointer, contiguous, dimension(:,:,:), intent(out) :: ee, tt, ne
  real(4), pointer, contiguous, dimension(:) :: xx, yy, zz
  real(4), pointer, contiguous, dimension(:) :: newx, newy, newz
  real(8) :: box(3)
  integer :: i, nx, ny, nz, m(3)
  integer, save :: itimer
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_snapshot_mod%read_stagger', itimer=itimer)
  call read_stagger_mesh(nx, ny, nz, xx, yy, zz)
  call read_stagger_size(box)

  ds = box / (uul*dims)

  call os%allocate(newx, dims(1))
  call os%allocate(newy, dims(2))
  call os%allocate(newz, dims(3))

  if (legacy_stagger_interp) then
    newx = linspace(xx(1), xx(nx), dims(1))
    newy = linspace(yy(1), yy(ny), dims(2))
    newz = linspace(zz(1), zz(nz), dims(3))

    ds(1) = (newx(dims(1)) - newx(1)) / (dims(1) - 1)
    ds(2) = (newy(dims(2)) - newy(1)) / (dims(2) - 1)
    ds(3) = (newz(dims(3)) - newz(1)) / (dims(3) - 1)
  else
    newx = [ (ds(1)*i, i=0, dims(1)-1) ]
    newy = [ (ds(2)*i, i=0, dims(2)-1) ]
    newz = [ (ds(3)*i, i=0, dims(3)-1) ]
  endif

  ds = ds * uul

  call os%allocate(d,  dims)
  call os%allocate(px, dims)
  call os%allocate(py, dims)
  call os%allocate(pz, dims)
  call os%allocate(ee, dims)
  call os%allocate(tt, dims)
  call os%allocate(ne, dims)

  call interpolate_stagger(xx, yy, zz, newx, newy, newz, d, px, py, pz, ee=ee, tt=tt, ne=ne)

  call os%deallocate(newx)
  call os%deallocate(newy)
  call os%deallocate(newz)

  call os%deallocate(xx)
  call os%deallocate(yy)
  call os%deallocate(zz)

  call deallocate_stagger()

  call trace%end(itimer)
END SUBROUTINE read_stagger_snapshot

!===============================================================================
SUBROUTINE load_stagger_subregion(n, sub, pos, d, px, py, pz, tt, ne, bot2pivot, vec)
  integer, intent(in) :: n(3)
  real(8), intent(in) :: sub(3), pos(3)
  real(8), intent(in), optional :: bot2pivot, vec(3)
  real(4), intent(out), dimension(:,:,:) :: d, px, py, pz
  real(4), intent(out), dimension(:,:,:), optional :: tt, ne
  real(4), pointer, contiguous, dimension(:) :: xx, yy, zz
  real(4), dimension(:), allocatable :: newx, newy, newz
  real(4) :: lower(3), ds(3)
  integer :: i, nx, ny, nz
  integer, save :: itimer
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_snapshot_mod%load_stagger_subregion', itimer=itimer)
  call read_stagger_mesh(nx, ny, nz, xx, yy, zz)

  allocate(newx(n(1)))
  allocate(newy(n(2)))
  allocate(newz(n(3)))

  lower = (pos - 0.5d0*sub) / uul
  ds = sub / (uul*n)

  newx = [ (ds(1)*i + lower(1), i=0, n(1)-1) ]
  newy = [ (ds(2)*i + lower(2), i=0, n(2)-1) ]
  newz = [ (ds(3)*i + lower(3), i=0, n(3)-1) ]

  call interpolate_stagger(xx, yy, zz, newx, newy, newz, d, px, py, pz, tt=tt, ne=ne, bot2pivot=bot2pivot, vec=vec)

  call os%deallocate(xx)
  call os%deallocate(yy)
  call os%deallocate(zz)

  call trace%end(itimer)
END SUBROUTINE load_stagger_subregion

!===============================================================================
SUBROUTINE deallocate_stagger()
  !-----------------------------------------------------------------------------
  deallocate(ilrho)
  deallocate(iltemp)
  deallocate(ilne)
  deallocate(ipx)
  deallocate(ipy)
  deallocate(ipz)
  deallocate(iie)
END SUBROUTINE deallocate_stagger

!===============================================================================
SUBROUTINE read_stagger_dim(nx, ny, nz, lunit_out)
  integer, intent(out) :: nx, ny, nz
  integer :: lunit, iostat
  integer, optional :: lunit_out
  !-----------------------------------------------------------------------------
  call os%check_file(mesh_file)
  open(newunit=lunit, file=mesh_file, iostat=iostat, convert='big_endian', form='unformatted')

  ! read first entry in mesh file
  read(unit=lunit, iostat=iostat) nx, nz, ny ! y is z in stagger

  ! Index the ghost zones to be ommited 
  !nz = nnz - (2 * gz)

  if (present(lunit_out)) then
    lunit_out = lunit
  else
    close(lunit)
  endif
END SUBROUTINE read_stagger_dim

!===============================================================================
!> Reading mesh file
!===============================================================================
SUBROUTINE read_stagger_mesh(nx, ny, nz, xx, yy, zz)
  integer, intent(out) :: nx, ny, nz
  real(4), intent(out), dimension(:), pointer, contiguous :: xx, yy, zz
  real(4), allocatable, dimension(:) :: xdn, dxidxup, dxidxdn, dx, dxdn
  real(4), allocatable, dimension(:) :: ydn, dyidyup, dyidydn, dy, dydn
  real(4), allocatable, dimension(:) :: zdn, dzidzup, dzidzdn, dz, dzdn
  real(4), allocatable, dimension(:) :: zm
  integer :: lunit, iostat
  !-----------------------------------------------------------------------------
  call read_stagger_dim(nx, ny, nz, lunit_out=lunit)

  call os%allocate(xx, nx)
  call os%allocate(yy, ny)
  call os%allocate(zz, nz)

  allocate(xdn, dxidxup, dxidxdn, dx, dxdn, mold=xx)
  allocate(ydn, dyidyup, dyidydn, dy, dydn, mold=yy)
  allocate(zdn, dzidzup, dzidzdn, dz, dzdn, mold=zz)

  read (unit=lunit, iostat=iostat) dx, dxdn, xx, xdn, dxidxup, dxidxdn
  read (unit=lunit, iostat=iostat) dz, dzdn, zz, zdn, dzidzup, dzidzdn ! y is z in stagger
  read (unit=lunit, iostat=iostat) dy, dydn, yy, ydn, dyidyup, dyidydn
  ! close mesh file
  close(unit=lunit)

  zz = -zz(nz:1:-1) ! flipping z-axis

  xx = xx - xx(1) ! let the origin be zero
  yy = yy - yy(1) ! let the origin be zero
  zz = zz - zz(gz+1) ! let the origin be the layer above ghost zone
  mindz = minval(dz)
END SUBROUTINE read_stagger_mesh

!===============================================================================
!> Reading data and aux file
!===============================================================================
SUBROUTINE read_stagger_data(lrho, px, py, pz, ie, ltemp, lne)
  real(4), intent(out), dimension(:,:,:), allocatable :: ie, ltemp, lne, lrho
  real(4), intent(out), dimension(:,:,:), allocatable :: px, py, pz
  real(4), dimension(:,:,:), allocatable :: ppx, ppy, ppz
  integer :: iostat, lunit, recl
  integer :: nx, ny, nz, m(3), i, j, zs, ze
  integer, save :: itimer
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_data', itimer=itimer)
  call read_stagger_dim(nx, ny, nz)

  m = [nx, ny, nz]
  allocate(lrho (nx,ny,nz))
  allocate(ie   (nx,ny,nz))
  allocate(ltemp(nx,ny,nz))
  allocate(px   (nx,ny,nz))
  allocate(py   (nx,ny,nz))
  allocate(pz   (nx,ny,nz))
  allocate(lne  (nx,ny,nz))
  allocate(ppx, mold=px)
  allocate(ppy, mold=py)
  allocate(ppz, mold=pz)

  recl = kind(1.0)*product(m)

  !=============================================================================
  ! Reading data file
  !-----------------------------------------------------------------------------
  call os%check_file(data_file)
  open(newunit=lunit, file=data_file, iostat=iostat, convert='big_endian', &
       form='unformatted', status='old', access='direct', recl=recl)

    ! 1: density
    ! 2: momentum per unit volume in direction x
    ! 3: momentum per unit volume in direction y
    ! 4: momentum per unit volume in direction z
    ! 5: internal energy per unit volume
    ! 6: temperature

  read (unit=lunit, iostat=iostat, rec=1) lrho
  read (unit=lunit, iostat=iostat, rec=2) ppx
  read (unit=lunit, iostat=iostat, rec=3) ppz ! y is z in stagger
  read (unit=lunit, iostat=iostat, rec=4) ppy ! z is y in stagger
  read (unit=lunit, iostat=iostat, rec=5) ie
  !read (unit=lunit, iostat=iostat, rec=6) temp
  ! close data file
  close(unit=lunit)

  !=============================================================================
  ! Reading aux file
  !-----------------------------------------------------------------------------
  call os%check_file(aux_file)
  open(newunit=lunit, file=aux_file, iostat=iostat, convert='big_endian', &
        form='unformatted', status='old', access='direct', recl=recl)

    ! 1: ln total pressure (gas + rad)
    ! 2: ln Rosseland extinction coeff.
    ! 3: ln Temperature
    ! 4: ln electron number density
    ! 5: ln Planck-averaged extinction coeff.
    ! 6: ln extinction coeff at 500 nm

  read (unit=lunit, iostat=iostat, rec=3) ltemp
  read (unit=lunit, iostat=iostat, rec=4) lne
  ! close aux file
  close(unit=lunit)
  !=============================================================================

  ! "un"-staggering vectors relative to scalars
  call xup_set(ppx, px)

  ! transposing data (y to z)
  call transpose_stagger(lrho)
  call transpose_stagger(ie)
  call transpose_stagger(ltemp)
  call transpose_stagger(lne)

  call transpose_stagger(px)
  call transpose_stagger(ppy)
  call transpose_stagger(ppz)

  call yup_set(ppy, py)
  call zup_set(ppz, pz)

  do j=1, ny
  do i=1, nx
    lrho (i,j,:) = log(lrho(i,j, nz:1:-1))
    ie   (i,j,:) =       ie(i,j, nz:1:-1)
    ltemp(i,j,:) =    ltemp(i,j, nz:1:-1)
    lne  (i,j,:) =      lne(i,j, nz:1:-1)
    px   (i,j,:) =       px(i,j, nz:1:-1)
    py   (i,j,:) =       py(i,j, nz:1:-1)
    pz   (i,j,:) =     - pz(i,j, nz:1:-1)
  end do
  end do

  call trace%end(itimer)
END SUBROUTINE read_stagger_data

!===============================================================================
SUBROUTINE read_stagger_size(box)
  real(8), intent(out) :: box(3)
  real(4), pointer, contiguous, dimension(:) :: xx, yy, zz
  real(8) :: ds(3)
  integer :: n(3), nx, ny, nz
  !-----------------------------------------------------------------------------
  call read_stagger_mesh(nx, ny, nz, xx, yy, zz)
  n = [nx, ny, nz-2*gz]
  box(1) = dble(xx(nx) - xx(1))
  box(2) = dble(yy(ny) - yy(1))
  box(3) = dble(zz(nz-gz) - zz(1+gz))
  ds = box / (n-1)
  box = ds * n
  box = box * uul
END SUBROUTINE read_stagger_size

!===============================================================================
SUBROUTINE get_stagger_subregion(ipos, nbox, nsub, ds, pos, zs)
  integer, intent(in) :: ipos
  integer, intent(in) :: nbox(3)
  integer, intent(inout) :: nsub(3)
  integer, intent(out) :: zs
  real(8), intent(out) :: ds(3), pos(3)
  real(8) :: box(3), sub(3)
  real(4), pointer, contiguous, dimension(:) :: xx, yy, zz
  integer :: i, nx, ny, nz, lvl, n(3)
  !-----------------------------------------------------------------------------
  call read_stagger_mesh(nx, ny, nz, xx, yy, zz)
  n = [nx, ny, nz]
  box(1) = xx(nx) - xx(1)
  box(2) = yy(ny) - yy(1)
  box(3) = zz(nz) - zz(1)
  box = box * uul
  ds = box / (n-1)
  box = ds * n
  sub = box

  !lvl = stagger_levels(ipos+1)
  lvl = 1
  nsub(1:2) = nbox(1:2) * 2**(lvl-1)
  sub(3) = level_ratio(ipos+1) * box(3)
  zs = nsub(3)*ipos

  ds = sub / nsub

  pos = 0.5d0 * sub
  pos(3) = pos(3) + sum(level_ratio(1:ipos)) * box(3)
END SUBROUTINE get_stagger_subregion

!===============================================================================
SUBROUTINE interpolate_stagger(xx, yy, zz, newx, newy, newz, d, px, py, pz, ee, tt, ne, vec, bot2pivot)
  real(4), intent(in), dimension(:) :: xx, yy, zz
  real(4), intent(in), dimension(:) :: newx, newy, newz
  real(4), intent(out), dimension(:,:,:) :: d, px, py, pz
  real(4), intent(out), dimension(:,:,:), optional :: tt, ne, ee
  real(8), intent(in), optional :: vec(3), bot2pivot
  real(4) :: dx, dy, dz
  integer, save :: itimer
  !-----------------------------------------------------------------------------
  call trace%begin('interpolate_stagger', itimer=itimer)

  if (present(vec)) then
    dx = xx(2)
    dy = yy(2)
    call linterpolate_3d_rotated(ilrho, d,  dx, dy, zz, newx, newy, newz, vec, bot2pivot, ray_offset)
    call linterpolate_3d_rotated( ipx, px, dx, dy, zz, newx, newy, newz, vec, bot2pivot, ray_offset)
    call linterpolate_3d_rotated( ipy, py, dx, dy, zz, newx, newy, newz, vec, bot2pivot, ray_offset)
    call linterpolate_3d_rotated( ipz, pz, dx, dy, zz, newx, newy, newz, vec, bot2pivot, ray_offset)

    if (present(ee)) then
      call linterpolate_3d_rotated( iie, ee, dx, dy, zz, newx, newy, newz, vec, bot2pivot, ray_offset)
      ee = ee / d + ee_offset
    endif

    if (present(tt)) then
      call linterpolate_3d_rotated(iltemp, tt, dx, dy, zz, newx, newy, newz, vec, bot2pivot, ray_offset)
    endif

    if (present(ne)) then
      call linterpolate_3d_rotated(ilne, ne, dx, dy, zz, newx, newy, newz, vec, bot2pivot, ray_offset)
    endif
  else
    call linterpolate_3d(ilrho, d,  xx, yy, zz, newx, newy, newz)
    call linterpolate_3d( ipx, px, xx, yy, zz, newx, newy, newz)
    call linterpolate_3d( ipy, py, xx, yy, zz, newx, newy, newz)
    call linterpolate_3d( ipz, pz, xx, yy, zz, newx, newy, newz)

    if (present(ee)) then
      call linterpolate_3d( iie, ee, xx, yy, zz, newx, newy, newz)
      ee = ee / d + ee_offset
    endif

    if (present(tt)) then
      call linterpolate_3d(iltemp, tt, xx, yy, zz, newx, newy, newz)
    endif

    if (present(ne)) then
      call linterpolate_3d(ilne, ne, xx, yy, zz, newx, newy, newz)
    endif
  endif

  d  = exp(d + luur)
  tt = exp(tt)
  ne = exp(ne - 3*luul)
  px = px * uum
  py = py * uum
  pz = pz * uum
  call trace%end(itimer)
END SUBROUTINE interpolate_stagger

!===============================================================================
SUBROUTINE xup_set (f, fp)
  ! f is centered on (i-.5,j,k), fifth order stagger
  real(4), dimension(:,:,:), intent(in) :: f
  real(4), dimension(:,:,:), intent(out) :: fp
  integer :: j, k
  !--------------------------------------------------------------------------
  do k=1, size(f, 3)
    do j=1, size(f, 2)
      call destagger_vector(f(:,j,k), fp(:,j,k))
    end do
  end do
END SUBROUTINE xup_set
    
!===============================================================================
SUBROUTINE yup_set (f, fp)
  ! f is centered on (i-.5,j,k), fifth order stagger
  real(4), dimension(:,:,:), intent(in) :: f
  real(4), dimension(:,:,:), intent(out) :: fp
  integer :: i, k
  !--------------------------------------------------------------------------
  do k=1, size(f, 3)
    do i=1, size(f, 1)
      call destagger_vector(f(i,:,k), fp(i,:,k))
    end do
  end do
END SUBROUTINE yup_set

!===============================================================================
SUBROUTINE zup_set (f, fp)
  ! f is centered on (i-.5,j,k), fifth order stagger
  real(4), dimension(:,:,:), intent(in) :: f
  real(4), dimension(:,:,:), intent(out) :: fp
  integer :: i, j
  !--------------------------------------------------------------------------
  do j=1, size(f, 2)
    do i=1, size(f, 1)
      call destagger_vector(f(i,j,:), fp(i,j,:))
    end do
  end do
END SUBROUTINE zup_set

!===============================================================================
SUBROUTINE destagger_vector(f, fout)
  real, dimension(:), intent(in)  :: f
  real, dimension(:), intent(out) :: fout
  integer :: i, n, j, idx
  real(4), parameter :: c = 3./256.
  real(4), parameter :: b = -25./256.
  real(4), parameter :: a = 0.5-b-c
  real(4) :: wts(-2:3)=[c,b,a,a,b,c]
  ! ----------------------------------------------------------------------------
  n = size(f,1)
  do i=1, n
    do j=-2,3
      idx = i + j
      if (idx<1) then
        idx = n + idx
      elseif (idx>n) then
        idx = idx-n
      endif  
      fout(i) = fout(i) + wts(j) * f(idx)
    end do
  end do
END SUBROUTINE destagger_vector

!===============================================================================
PURE SUBROUTINE transpose_stagger(f)
  real, dimension(:,:,:), intent(inout)  :: f
  integer :: i
  ! ----------------------------------------------------------------------------
  do i=1, size(f, 1)
    f(i, :, :) = transpose(f(i, :, :))
  end do
END SUBROUTINE transpose_stagger

END MODULE read_stagger_snapshot_mod