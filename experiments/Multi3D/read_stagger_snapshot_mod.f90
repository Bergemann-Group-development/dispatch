MODULE read_stagger_snapshot_mod
  USE io_mod
  USE os_mod
  USE trace_mod
  USE handyfuncs_mod
  USE const_mod
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
  real(4), save :: maxTemp=huge(1.0)
  real(4), save :: ee_offset=3.2122
  logical, public, save :: legacy_stagger_interp=.false.
  real(4), dimension(:,:,:), allocatable, save :: iie, iltemp, ilpe, ilrho
  real(4), dimension(:,:,:), allocatable, save :: ipx, ipy, ipz
  real(4), dimension(:), allocatable, save :: ixx, iyy, izz
  real(4), public, dimension(:), allocatable, save :: amr_ref
  real(8), public, dimension(:), allocatable, save :: sub_pos, sub_size

  public :: read_stagger
  public :: read_co5bold
  public :: deallocate_stagger
  public :: interpolate_atmos
  public :: setup_amr_mesh
    
contains

!===============================================================================
SUBROUTINE read_stagger(stagger_file, snap, box)
  character(len=*), intent(in) :: stagger_file
  integer, intent(in) :: snap
  real(8), intent(out) :: box(3)
  real(8) :: ds(3)
  integer :: n(3), nx, ny, nz
  !-----------------------------------------------------------------------------
  mesh_file = trim(stagger_file) // '.msh'
  write(snap_file, '(A,A1,I5.5)') trim(stagger_file), '_', snap
  write(data_file, '(A,A)') trim(snap_file), '.dat'
  write(aux_file,  '(A,A)') trim(snap_file), '.aux'
  call read_stagger_mesh(nx, ny, nz)
  n = [nx, ny, nz-2*gz]
  box(1) = dble(ixx(nx) - ixx(1))
  box(2) = dble(iyy(ny) - iyy(1))
  box(3) = dble(izz(nz-gz) - izz(1+gz))
  ds = box / (n-1)
  box = ds * n
  call read_stagger_data()
END SUBROUTINE read_stagger

!===============================================================================
SUBROUTINE read_co5bold(atmos_file, box)
  character(len=*), intent(in) :: atmos_file
  real(8), intent(out) :: box(3)
  integer :: lunit, iostat, nsub
  real(4) :: model_time
  real(8) :: ds(3)
  integer :: n(3), nx, ny, nz
  !-----------------------------------------------------------------------------
  open(newunit=lunit, file=atmos_file, form='unformatted', access='stream')
  read(lunit) model_time
  read(lunit) nx, ny, nz

  allocate(ixx(nx))
  allocate(iyy(ny))
  allocate(izz(nz))
  allocate(iltemp(nx,ny,nz))
  allocate(ilrho (nx,ny,nz))
  allocate(ilpe  (nx,ny,nz))
  allocate(ipx   (nx,ny,nz))
  allocate(ipy   (nx,ny,nz))
  allocate(ipz   (nx,ny,nz))

  read(lunit) ixx, iyy, izz
  read(lunit) iltemp
  read(lunit) ilrho 
  read(lunit) ilpe  
  read(lunit) ipx
  read(lunit) ipy
  read(lunit) ipz
  close(lunit)

  ixx = ixx - ixx(1) ! let the origin be zero
  iyy = iyy - iyy(1) ! let the origin be zero
  izz = izz - izz(1) ! let the origin be zero

  n = [nx, ny, nz]
  box(1) = dble(ixx(nx))
  box(2) = dble(iyy(ny))
  box(3) = dble(izz(nz))
  ds = box / (n-1)
  box = ds * n
END SUBROUTINE read_co5bold

!===============================================================================
SUBROUTINE deallocate_stagger()
  !-----------------------------------------------------------------------------
  if (allocated(ilrho))  deallocate(ilrho)
  if (allocated(iltemp)) deallocate(iltemp)
  if (allocated(ilpe))   deallocate(ilpe)
  if (allocated(ipx))    deallocate(ipx)
  if (allocated(ipy))    deallocate(ipy)
  if (allocated(ipz))    deallocate(ipz)
  if (allocated(iie))    deallocate(iie)
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
SUBROUTINE read_stagger_mesh(nx, ny, nz)
  integer, intent(out) :: nx, ny, nz
  real(4), allocatable, dimension(:) :: xdn, dxidxup, dxidxdn, dx, dxdn
  real(4), allocatable, dimension(:) :: ydn, dyidyup, dyidydn, dy, dydn
  real(4), allocatable, dimension(:) :: zdn, dzidzup, dzidzdn, dz, dzdn
  real(4), allocatable, dimension(:) :: zm
  integer :: lunit, iostat
  !-----------------------------------------------------------------------------
  call read_stagger_dim(nx, ny, nz, lunit_out=lunit)

  allocate(ixx(nx))
  allocate(iyy(ny))
  allocate(izz(nz))

  allocate(xdn, dxidxup, dxidxdn, dx, dxdn, mold=ixx)
  allocate(ydn, dyidyup, dyidydn, dy, dydn, mold=iyy)
  allocate(zdn, dzidzup, dzidzdn, dz, dzdn, mold=izz)

  read (unit=lunit, iostat=iostat) dx, dxdn, ixx, xdn, dxidxup, dxidxdn
  read (unit=lunit, iostat=iostat) dz, dzdn, izz, zdn, dzidzup, dzidzdn ! y is z in stagger
  read (unit=lunit, iostat=iostat) dy, dydn, iyy, ydn, dyidyup, dyidydn
  ! close mesh file
  close(unit=lunit)

  izz = -izz(nz:1:-1) ! flipping z-axis

  ixx = uul * (ixx - ixx(1))    ! let the origin be zero
  iyy = uul * (iyy - iyy(1))    ! let the origin be zero
  izz = uul * (izz - izz(gz+1)) ! let the origin be the layer above ghost zone
END SUBROUTINE read_stagger_mesh

!===============================================================================
!> Reading data and aux file
!===============================================================================
SUBROUTINE read_stagger_data()
  real(4), dimension(:,:,:), allocatable :: ppx, ppy, ppz
  integer :: iostat, lunit, recl
  integer :: nx, ny, nz, m(3), i, j, zs, ze
  integer, save :: itimer
  !-----------------------------------------------------------------------------
  call trace%begin('read_stagger_data', itimer=itimer)
  call read_stagger_dim(nx, ny, nz)

  m = [nx, ny, nz]
  allocate(ilrho (nx,ny,nz))
  allocate(iie   (nx,ny,nz))
  allocate(iltemp(nx,ny,nz))
  allocate(ipx   (nx,ny,nz))
  allocate(ipy   (nx,ny,nz))
  allocate(ipz   (nx,ny,nz))
  allocate(ilpe  (nx,ny,nz))
  allocate(ppx, mold=ipx)
  allocate(ppy, mold=ipy)
  allocate(ppz, mold=ipz)

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

  read (unit=lunit, iostat=iostat, rec=1) ilrho
  read (unit=lunit, iostat=iostat, rec=2) ppx
  read (unit=lunit, iostat=iostat, rec=3) ppz ! y is z in stagger
  read (unit=lunit, iostat=iostat, rec=4) ppy ! z is y in stagger
  read (unit=lunit, iostat=iostat, rec=5) iie
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

  read (unit=lunit, iostat=iostat, rec=3) iltemp
  read (unit=lunit, iostat=iostat, rec=4) ilpe
  ! close aux file
  close(unit=lunit)
  !=============================================================================

  ! "un"-staggering vectors relative to scalars
  call xup_set(ppx, ipx)

  ! transposing data (y to z)
  call transpose_stagger(ilrho)
  call transpose_stagger(iie)
  call transpose_stagger(iltemp)
  call transpose_stagger(ilpe)

  call transpose_stagger(ipx)
  call transpose_stagger(ppy)
  call transpose_stagger(ppz)

  call yup_set(ppy, ipy)
  call zup_set(ppz, ipz)

  do j=1, ny
  do i=1, nx
    ilrho (i,j,:) = log(ilrho(i,j, nz:1:-1)) + luur
    iie   (i,j,:) =       iie(i,j, nz:1:-1) 
    iltemp(i,j,:) =    iltemp(i,j, nz:1:-1) 
    ilpe  (i,j,:) =      ilpe(i,j, nz:1:-1)  - 3*luul + log(BK) + iltemp(i,j,:)
    ipx   (i,j,:) =       ipx(i,j, nz:1:-1)  * uum
    ipy   (i,j,:) =       ipy(i,j, nz:1:-1)  * uum
    ipz   (i,j,:) =     - ipz(i,j, nz:1:-1)  * uum
  end do
  end do

  call trace%end(itimer)
END SUBROUTINE read_stagger_data

!===============================================================================
SUBROUTINE interpolate_atmos(n, sub, pos, d, px, py, pz, ee, tt, pe, vec, bot2pivot)
  integer, intent(in) :: n(3)
  real(8), intent(in) :: sub(3), pos(3)
  real(4), intent(out), dimension(:,:,:) :: d, px, py, pz
  real(4), intent(out), dimension(:,:,:), optional :: tt, pe, ee
  real(8), intent(in), optional :: vec(3), bot2pivot
  real(8), dimension(:), allocatable :: newx, newy, newz
  real(8) :: lower(3), ds(3)
  real(8) :: dx, dy, dz
  integer :: i, itimer=0
  !-----------------------------------------------------------------------------
  call trace%begin('interpolate_atmos', itimer=itimer)

  allocate(newx(n(1)))
  allocate(newy(n(2)))
  allocate(newz(n(3)))

  lower = pos - 0.5d0*sub
  ds = sub / n

  newx = [ (ds(1)*i + lower(1), i=0, n(1)-1) ]
  newy = [ (ds(2)*i + lower(2), i=0, n(2)-1) ]

  !newz = [ (ds(3)*i + lower(3), i=0, n(3)-1) ]
  newz = [ (ds(3)*i + lower(3), i=1, n(3)) ]
  ! this way the distance to the patch below is ds of the current patch, otherwise it is ds of the patch below
  ! this is needed for the amr setup in order to avoid additional interpolations for the final spectrum.
  
  dx = ixx(2)
  dy = iyy(2)

  if (present(vec)) then
    call linterpolate_3d_rotated(ilrho, d,  dx, dy, izz, newx, newy, newz, vec, bot2pivot, ray_offset)
    call linterpolate_3d_rotated(  ipx, px, dx, dy, izz, newx, newy, newz, vec, bot2pivot, ray_offset)
    call linterpolate_3d_rotated(  ipy, py, dx, dy, izz, newx, newy, newz, vec, bot2pivot, ray_offset)
    call linterpolate_3d_rotated(  ipz, pz, dx, dy, izz, newx, newy, newz, vec, bot2pivot, ray_offset)

    if (present(ee)) then
      call linterpolate_3d_rotated( iie, ee, dx, dy, izz, newx, newy, newz, vec, bot2pivot, ray_offset)
      ee = ee / d + ee_offset
    endif

    if (present(tt)) then
      call linterpolate_3d_rotated(iltemp, tt, dx, dy, izz, newx, newy, newz, vec, bot2pivot, ray_offset)
    endif

    if (present(pe)) then
      call linterpolate_3d_rotated(ilpe, pe, dx, dy, izz, newx, newy, newz, vec, bot2pivot, ray_offset)
    endif
  else
    call linterpolate_3d(ilrho, d,  dx, dy, izz, newx, newy, newz)
    call linterpolate_3d(  ipx, px, dx, dy, izz, newx, newy, newz)
    call linterpolate_3d(  ipy, py, dx, dy, izz, newx, newy, newz)
    call linterpolate_3d(  ipz, pz, dx, dy, izz, newx, newy, newz)

    if (present(ee)) then
      call linterpolate_3d( iie, ee, dx, dy, izz, newx, newy, newz)
      ee = ee / d + ee_offset
    endif

    if (present(tt)) then
      call linterpolate_3d(iltemp, tt, dx, dy, izz, newx, newy, newz)
    endif

    if (present(pe)) then
      call linterpolate_3d(ilpe, pe, dx, dy, izz, newx, newy, newz)
    endif
  endif

  d  = exp(d)
  tt = exp(tt)
  pe = exp(pe)
  call trace%end(itimer)
END SUBROUTINE interpolate_atmos

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
SUBROUTINE transpose_stagger(f)
  real, dimension(:,:,:), intent(inout)  :: f
  integer :: i, itimer=0
  ! ----------------------------------------------------------------------------
  call trace%begin('transpose_stagger', itimer=itimer)
  do i=1, size(f, 1)
    f(i, :, :) = transpose(f(i, :, :))
  end do
  call trace%end(itimer)
END SUBROUTINE transpose_stagger

!===============================================================================
SUBROUTINE setup_amr_mesh(nsub)
  integer, intent(in) :: nsub
  integer :: i, j, k, nx, ny, nz
  real(8), dimension(nsub+1) :: var_edge
  real(4), dimension(nsub+1) :: sub_edge
  real(4) :: vmin, vmax
  !-----------------------------------------------------------------------------
  nx = size(ixx)
  ny = size(iyy)
  nz = size(izz)-2*gz
  allocate(amr_ref(nz), source=0.0)
  allocate(sub_pos(nsub), sub_size(nsub))

  do k=1, nz
  do j=1, ny
  do i=1, nx
    ! we use minus here, because interp1d_sorted expects increasing values
    amr_ref(k) = amr_ref(k) - ilpe(i,j,k+gz)
  end do
  end do
  end do
  amr_ref = amr_ref / (nx*ny)

  vmin = amr_ref(1)
  vmax = amr_ref(nz)
  var_edge = linspace(vmin, vmax, nsub+1)
  call interp1d_sorted(amr_ref, izz(gz+1:nz+gz), var_edge, sub_edge)
  sub_size = sub_edge(2:nsub+1) - sub_edge(1:nsub)
  sub_pos = sub_edge(1:nsub) + 0.5d0*sub_size
END SUBROUTINE setup_amr_mesh

END MODULE read_stagger_snapshot_mod